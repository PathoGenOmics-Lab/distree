use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::process::ExitCode;

use distree::lca::build_lca_structure;
use distree::midpoint::midpoint_root;
use distree::parser::{flatten_raw, parse_newick};
use distree::tree::Node;
use distree::{compute_distance, DistMode};

/// Largest accepted `--precision`.
///
/// An f64 carries about 17 significant decimal digits, so anything past this
/// is zero padding. The cap also keeps the formatting machinery inside the
/// range it accepts: `{:.prec$}` panics outright on very large values.
const MAX_PRECISION: usize = 30;

/// Target number of cells per parallel batch, which caps the batch buffer at
/// 8 MB and gives each fork/join enough work to be worth its cost.
const BATCH_CELLS: usize = 1 << 20;

/// Output buffer size. The default 8 KB turns a multi-gigabyte matrix into
/// hundreds of thousands of write syscalls.
const WRITE_BUFFER: usize = 1 << 20;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // A closed downstream pipe ("distree tree.nwk | head") is how a
            // reader says it has seen enough, not a failure to report.
            if let Some(io_err) = e.downcast_ref::<io::Error>() {
                if io_err.kind() == io::ErrorKind::BrokenPipe {
                    return ExitCode::SUCCESS;
                }
            }
            eprintln!("Error: {}", e);
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("distree")
        .version(env!("CARGO_PKG_VERSION"))
        .author("Paula Ruiz-Rodriguez")
        .about("Extracts a distance matrix from a phylogeny (parallel, low-memory)")
        .arg(
            Arg::new("phylogeny")
                .help("Path to the tree file in Newick format (use '-' for stdin)")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("midpoint")
                .long("midpoint")
                .help("Midpoint-root the tree before computing distances")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("lmm")
                .long("lmm")
                .help("Produce the var-covar matrix C (depth of the MRCA in branch lengths)")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("topology")
                .long("topology")
                .help("Ignore branch lengths; use purely topological distances")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .help("Path to write the TSV output file (defaults to stdout)")
                .value_name("FILE"),
        )
        .arg(
            Arg::new("precision")
                .long("precision")
                .short('p')
                .help("Number of decimal places for output values")
                .default_value("10")
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .help("Number of threads for parallel computation (default: all cores)")
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(
            Arg::new("lower")
                .long("lower")
                .help("Output PHYLIP lower-triangle format (taxa count header, row labels, no diagonal)")
                .action(ArgAction::SetTrue),
        )
        .get_matches();

    let tree_path = matches
        .get_one::<String>("phylogeny")
        .ok_or("Tree file path is required")?
        .to_string();
    let do_midpoint = *matches.get_one::<bool>("midpoint").unwrap_or(&false);
    let do_lmm = *matches.get_one::<bool>("lmm").unwrap_or(&false);
    let do_topology = *matches.get_one::<bool>("topology").unwrap_or(&false);
    let do_lower = *matches.get_one::<bool>("lower").unwrap_or(&false);
    let output_path = matches.get_one::<String>("output");
    let precision = *matches.get_one::<usize>("precision").unwrap_or(&10);

    if precision > MAX_PRECISION {
        return Err(format!(
            "--precision must be between 0 and {}, got {}. A 64-bit float holds about \
             17 significant digits, so more decimals only add zeros.",
            MAX_PRECISION, precision
        )
        .into());
    }

    // Determine distance mode
    let mode = if do_lmm {
        if do_topology {
            eprintln!("Warning: --lmm and --topology are mutually exclusive. Using --lmm.");
        }
        DistMode::Lmm
    } else if do_topology {
        DistMode::Topology
    } else {
        DistMode::Patristic
    };

    // Configure thread pool
    if let Some(&num_threads) = matches.get_one::<usize>("threads") {
        if num_threads == 0 {
            // Rayon reads 0 as "pick the default", which silently ignores what
            // the user asked for.
            return Err("--threads must be at least 1. Omit the flag to use all cores.".into());
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .map_err(|e| format!("Failed to initialize thread pool: {}", e))?;
    }

    // Read input from file or stdin
    let source = if tree_path == "-" { "standard input" } else { tree_path.as_str() };
    let mut raw_input: Vec<u8> = Vec::new();
    if tree_path == "-" {
        io::stdin().read_to_end(&mut raw_input)?;
    } else {
        File::open(&tree_path)
            .map_err(|e| format!("Cannot open '{}': {}", tree_path, e))?
            .read_to_end(&mut raw_input)?;
    }

    // Large trees are usually shipped compressed, and a gzip file reaching the
    // UTF-8 check produced "stream did not contain valid UTF-8", which says
    // nothing about what to do next.
    if raw_input.starts_with(&[0x1f, 0x8b]) {
        return Err(format!(
            "{} is gzip-compressed. distree reads plain text; decompress it on the way in:\n\
             \x20   gunzip -c {} | distree -",
            source,
            if tree_path == "-" { "FILE.nwk.gz" } else { tree_path.as_str() }
        )
        .into());
    }

    let newick_str = String::from_utf8(raw_input).map_err(|_| {
        format!(
            "{} is not valid UTF-8 text, so it cannot be a Newick tree.",
            source
        )
    })?;

    // Parse the Newick string
    let raw_root = parse_newick(&newick_str)
        .map_err(|e| format!("Failed to parse Newick tree: {}", e))?;

    // Flatten into Vec<Node>
    let mut nodes: Vec<Node> = Vec::new();
    let mut root_idx = flatten_raw(&raw_root, None, &mut nodes);

    // Nothing below reads the input text or the recursive parse tree, and on a
    // large tree they are a second copy of it. Holding them to the end of the
    // run would carry that alongside the LCA table, which is the peak.
    drop(raw_root);
    drop(newick_str);

    // Warn about negative branch lengths
    if nodes.iter().any(|n| n.length < 0.0) {
        eprintln!(
            "Warning: negative branch lengths detected in the tree. Some distances may come \
             out negative, and --midpoint cannot locate the diameter reliably."
        );
    }

    // Midpoint-root if requested.
    //
    // Topological distance counts edges and does not depend on where the tree
    // is rooted, so rooting can only distort it: the node inserted at the
    // midpoint splits one edge in two and adds 1 to every pair whose path
    // crosses it. Skip the rooting rather than report those inflated counts.
    if do_midpoint {
        if mode == DistMode::Topology {
            eprintln!(
                "Warning: --midpoint is ignored in --topology mode. Edge counts do not depend \
                 on the root, and the node inserted at the midpoint would add 1 to every \
                 distance crossing that edge."
            );
        } else {
            root_idx = midpoint_root(root_idx, &mut nodes);
        }
    }

    // Build LCA
    let lca_data = build_lca_structure(root_idx, &nodes);

    // A distance is d_i + d_j - 2*d_m, so a root-to-tip depth past half of
    // f64's range overflows that sum to infinity, and inf - inf is NaN. The
    // matrix comes out as inf with NaN down the diagonal, which is not a
    // failure any downstream tool would notice.
    let deepest = lca_data.depth_len.iter().copied().fold(0.0_f64, f64::max);
    if !(deepest * 2.0).is_finite() {
        return Err(format!(
            "Branch lengths are too large to compute distances with: the deepest leaf is {:e} \
             from the root, and summing two such depths overflows a 64-bit float. Rescale the \
             tree.",
            deepest
        )
        .into());
    }

    // Collect leaves
    let leaf_indices: Vec<usize> = nodes
        .iter()
        .enumerate()
        .filter(|(_, nd)| nd.children.is_empty() && nd.name.is_some())
        .map(|(i, _)| i)
        .collect();

    if leaf_indices.is_empty() {
        return Err("No labeled leaves found in the tree.".into());
    }

    // Leaves without a label cannot be named in the matrix, so they are left
    // out entirely. Say so rather than quietly returning a smaller matrix.
    let unlabeled_leaves = nodes
        .iter()
        .filter(|nd| nd.children.is_empty() && nd.name.is_none())
        .count();
    if unlabeled_leaves > 0 {
        eprintln!(
            "Warning: {} leaf/leaves have no label and were excluded from the matrix.",
            unlabeled_leaves
        );
    }

    // Check for duplicate leaf names and for characters that would break the
    // row/column structure of the output
    {
        let mut seen = HashSet::with_capacity(leaf_indices.len());
        for &i in &leaf_indices {
            let name = nodes[i].name.as_ref().expect("leaf has name");
            if !seen.insert(name.as_str()) {
                return Err(format!(
                    "Duplicate leaf name '{}' found. Leaf names must be unique.",
                    name
                )
                .into());
            }
            if let Some(bad) = name.chars().find(|c| matches!(c, '\t' | '\n' | '\r')) {
                let what = match bad {
                    '\t' => "a tab character",
                    '\n' => "a newline",
                    _ => "a carriage return",
                };
                return Err(format!(
                    "Leaf name '{}' contains {}, which would corrupt the output by splitting \
                     the row. Use an underscore or rename the leaf.",
                    name.escape_debug(),
                    what
                )
                .into());
            }
            // A PHYLIP row is "name<whitespace>values", so a reader that splits
            // on whitespace reads "Taxon A<TAB>4.5" as three fields and every
            // row after it lands one column out. The square TSV is delimited by
            // tabs alone and does not have the problem.
            if do_lower && name.chars().any(char::is_whitespace) {
                return Err(format!(
                    "Leaf name '{}' contains whitespace, which PHYLIP readers treat as the \
                     end of the name, so --lower would produce a file they misread. Use an \
                     underscore, or drop --lower for the square TSV.",
                    name
                )
                .into());
            }
        }
    }

    // Sort by label
    let mut leaf_label_pairs: Vec<(String, usize)> = leaf_indices
        .iter()
        .map(|&i| (nodes[i].name.clone().expect("leaf has name"), i))
        .collect();
    leaf_label_pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    let sorted_labels: Vec<&str> = leaf_label_pairs
        .iter()
        .map(|(lab, _)| lab.as_str())
        .collect();
    let sorted_leaf_indices: Vec<usize> =
        leaf_label_pairs.iter().map(|(_, idx)| *idx).collect();
    let n_leaves = sorted_leaf_indices.len();

    // Warn if no branch lengths and not topology mode
    if mode == DistMode::Patristic && nodes.iter().all(|n| n.length == 0.0) {
        eprintln!(
            "Warning: no branch lengths detected, all patristic distances will be zero. Consider --topology."
        );
    }

    // Open the output only once the tree is known to be usable: creating it
    // earlier truncated the previous results before a parse error could be
    // reported, leaving the user with an empty file and nothing to fall back on.
    let mut writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::with_capacity(
            WRITE_BUFFER,
            File::create(path).map_err(|e| format!("Cannot write to '{}': {}", path, e))?,
        ))
    } else {
        Box::new(BufWriter::with_capacity(WRITE_BUFFER, io::stdout()))
    };

    // The lower triangle has no diagonal. For a distance matrix that loses
    // nothing, since it is all zeros, but an LMM diagonal is each tip's own
    // root-to-tip length: real data, and what the off-diagonal covariances have
    // to be read against.
    if do_lower && mode == DistMode::Lmm {
        eprintln!(
            "Warning: --lower omits the diagonal, which in --lmm mode holds each leaf's \
             root-to-tip length rather than zeros. Drop --lower to keep it."
        );
    }

    // Print header
    if do_lower {
        // PHYLIP format: first line is the number of taxa
        writeln!(writer, "{}", n_leaves)?;
    } else {
        writer.write_all(b"\t")?;
        for (i, lab) in sorted_labels.iter().enumerate() {
            writer.write_all(lab.as_bytes())?;
            if i + 1 < n_leaves {
                writer.write_all(b"\t")?;
            }
        }
        writer.write_all(b"\n")?;
    }

    // Compute and print the distance matrix, a batch of rows at a time.
    //
    // The batch is sized so each one carries roughly BATCH_CELLS cells of work,
    // which bounds the buffer at 8 MB and, more to the point, keeps the
    // fork/join out of the inner loop. A cell is an MRCA query and three array
    // reads, so a single row of a small matrix is nowhere near enough work to
    // pay for synchronising a thread pool: on an 8,000-leaf tree, one job per
    // row ran *slower* on 14 cores than on one.
    // Each worker formats its rows as it computes them, so the writer only ever
    // hands finished bytes to the operating system. Formatting a float to a
    // fixed number of decimals costs more than the distance it prints, so
    // leaving it in the serial write loop capped the whole run: at six decimals
    // it was about two thirds of the work no number of cores could touch.
    let batch_rows = (BATCH_CELLS / n_leaves).clamp(1, n_leaves);
    let mut batch: Vec<Vec<u8>> = vec![Vec::new(); batch_rows];
    let row_width = |row_i: usize| if do_lower { row_i } else { n_leaves };

    let mut first_row = 0;
    while first_row < n_leaves {
        let rows = (n_leaves - first_row).min(batch_rows);

        batch[..rows].par_iter_mut().enumerate().for_each(|(k, out)| {
            let row_i = first_row + k;
            let leaf_i = sorted_leaf_indices[row_i];
            out.clear();
            out.extend_from_slice(sorted_labels[row_i].as_bytes());
            for &leaf_j in &sorted_leaf_indices[..row_width(row_i)] {
                let dist = compute_distance(leaf_i, leaf_j, mode, &lca_data);
                out.push(b'\t');
                format_distance(out, dist, mode, precision);
            }
            out.push(b'\n');
        });

        for row in &batch[..rows] {
            writer.write_all(row)?;
        }

        first_row += rows;
    }

    // BufWriter flushes on drop but discards whatever error it hits, so a full
    // disk or a failing filesystem produced a truncated matrix and exit code 0.
    if let Err(e) = writer.flush() {
        if e.kind() == io::ErrorKind::BrokenPipe {
            return Ok(());
        }
        return Err(match output_path {
            Some(path) => format!("Failed to write '{}': {}", path, e).into(),
            None => format!("Failed to write to stdout: {}", e).into(),
        });
    }

    Ok(())
}

/// Append a single distance value to a row buffer.
///
/// Topology mode outputs integers; patristic and LMM use `precision` decimal
/// places. Writing into a `Vec<u8>` cannot fail, which is what lets this run
/// inside the parallel row builder.
#[inline]
fn format_distance(out: &mut Vec<u8>, dist: f64, mode: DistMode, precision: usize) {
    let result = if mode == DistMode::Topology {
        write!(out, "{}", dist as i64)
    } else {
        write!(out, "{:.prec$}", dist, prec = precision)
    };
    result.expect("writing to a Vec<u8> cannot fail");
}
