use clap::{Arg, ArgAction, Command};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::process::ExitCode;
use std::time::Instant;

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
    let started = Instant::now();
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
        .arg(
            Arg::new("taxa")
                .long("taxa")
                .help("Restrict the matrix to the leaf labels listed in FILE, one per line")
                .value_name("FILE"),
        )
        .arg(
            Arg::new("stats")
                .long("stats")
                .help("Print a summary of the run to stderr")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("npy")
                .long("npy")
                .help("Write a NumPy .npy array instead of text (requires -o; labels go to <FILE>.labels.txt)")
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

    let do_npy = *matches.get_one::<bool>("npy").unwrap_or(&false);

    if do_npy && output_path.is_none() {
        return Err(
            "--npy writes binary, so it needs a file: add -o FILE. The labels go to \
             FILE.labels.txt beside it."
                .into(),
        );
    }
    if do_npy && do_topology && !do_lmm {
        // Edge counts are integers; an f64 array of them is honest but odd, so
        // say what is happening rather than let it surprise anyone.
        eprintln!("Warning: --npy writes 64-bit floats, so --topology edge counts are stored as floats.");
    }

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

    // Large trees are usually shipped compressed, so read gzip directly rather
    // than making the caller pipe it through gunzip. Detected by magic bytes
    // rather than by extension, so it works from stdin too.
    if raw_input.starts_with(&[0x1f, 0x8b]) {
        let mut decoded = Vec::new();
        GzDecoder::new(&raw_input[..])
            .read_to_end(&mut decoded)
            .map_err(|e| format!("{} is gzip-compressed but could not be read: {}", source, e))?;
        raw_input = decoded;
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

    // Restrict to the requested subset, if any. The distances themselves are
    // unaffected: the path between two leaves does not depend on which other
    // leaves are in the matrix, so this filters the output rather than pruning
    // the tree.
    if let Some(path) = matches.get_one::<String>("taxa") {
        let wanted = read_taxa_list(path)?;
        let present: HashSet<&str> = leaf_label_pairs.iter().map(|(l, _)| l.as_str()).collect();

        let missing: Vec<&str> = wanted
            .iter()
            .filter(|w| !present.contains(w.as_str()))
            .map(String::as_str)
            .collect();
        if !missing.is_empty() {
            // Silently dropping these would give a matrix smaller than asked
            // for, which is the kind of thing nobody notices until much later.
            let shown = missing.iter().take(5).cloned().collect::<Vec<_>>().join(", ");
            return Err(format!(
                "{} of the {} labels in '{}' are not leaves of the tree: {}{}",
                missing.len(),
                wanted.len(),
                path,
                shown,
                if missing.len() > 5 { ", ..." } else { "" }
            )
            .into());
        }

        let keep: HashSet<&str> = wanted.iter().map(String::as_str).collect();
        leaf_label_pairs.retain(|(label, _)| keep.contains(label.as_str()));
        if leaf_label_pairs.is_empty() {
            return Err(format!("'{}' selected no leaves.", path).into());
        }
    }

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
    if do_npy {
        // NumPy needs the shape up front, which is known, so the array still
        // streams out a row at a time behind this.
        let shape = if do_lower {
            // The condensed form scipy.spatial.distance.squareform expects
            format!("({},)", n_leaves * (n_leaves - 1) / 2)
        } else {
            format!("({}, {})", n_leaves, n_leaves)
        };
        writer.write_all(&npy_header(&shape))?;

        // .npy has nowhere to put the labels, so they go beside it in the
        // order of the rows.
        let labels_path = format!("{}.labels.txt", output_path.expect("--npy requires -o"));
        let mut labels = String::with_capacity(n_leaves * 12);
        for lab in &sorted_labels {
            labels.push_str(lab);
            labels.push('\n');
        }
        std::fs::write(&labels_path, labels)
            .map_err(|e| format!("Cannot write the labels to '{}': {}", labels_path, e))?;
    } else if do_lower {
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
    let mut batch: Vec<(Vec<u8>, Stats)> = vec![(Vec::new(), Stats::default()); batch_rows];

    // Which columns a row carries. The two triangular forms are not the same
    // triangle: PHYLIP wants row i to hold columns 0..i, while the condensed
    // vector SciPy reads wants columns i+1..n, which concatenated over the rows
    // is its (0,1), (0,2), ... (1,2), ... ordering. Emitting the PHYLIP order
    // into a .npy would hand squareform the right values in the wrong places.
    let row_range = |row_i: usize| -> std::ops::Range<usize> {
        match (do_lower, do_npy) {
            (false, _) => 0..n_leaves,
            (true, false) => 0..row_i,
            (true, true) => (row_i + 1).min(n_leaves)..n_leaves,
        }
    };
    let mut totals = Stats::default();

    let mut first_row = 0;
    while first_row < n_leaves {
        let rows = (n_leaves - first_row).min(batch_rows);

        batch[..rows]
            .par_iter_mut()
            .enumerate()
            .for_each(|(k, (out, stats))| {
                let row_i = first_row + k;
                let leaf_i = sorted_leaf_indices[row_i];
                out.clear();
                if !do_npy {
                    out.extend_from_slice(sorted_labels[row_i].as_bytes());
                }
                *stats = Stats::default();
                let cols = row_range(row_i);
                let first_col = cols.start;
                for (offset, &leaf_j) in sorted_leaf_indices[cols].iter().enumerate() {
                    let col = first_col + offset;
                    let dist = compute_distance(leaf_i, leaf_j, mode, &lca_data);
                    if do_npy {
                        out.extend_from_slice(&dist.to_le_bytes());
                    } else {
                        out.push(b'\t');
                        format_distance(out, dist, mode, precision);
                    }
                    // The diagonal is zero by construction in the distance
                    // modes and the leaf's own depth under --lmm, so it would
                    // only drag the summary around.
                    if col != row_i {
                        stats.add(dist);
                    }
                }
                if !do_npy {
                    out.push(b'\n');
                }
            });

        for (row, stats) in &batch[..rows] {
            writer.write_all(row)?;
            totals.merge(stats);
        }

        first_row += rows;
    }

    if *matches.get_one::<bool>("stats").unwrap_or(&false) {
        let unit = match mode {
            DistMode::Patristic => "patristic distance",
            DistMode::Topology => "edge count",
            DistMode::Lmm => "root-to-MRCA depth",
        };
        eprintln!("--- Statistics ---");
        eprintln!("Leaves in matrix:  {}", n_leaves);
        if n_leaves != leaf_indices.len() {
            eprintln!("Leaves in tree:    {}", leaf_indices.len());
        }
        eprintln!("Nodes in tree:     {}", nodes.len());
        eprintln!("Mode:              {}", unit);
        if do_midpoint && mode != DistMode::Topology {
            eprintln!("Rooting:           midpoint");
        }
        eprintln!("Cells written:     {}", totals.count + n_leaves as u64 * u64::from(!do_lower));
        if totals.count > 0 {
            // Off the diagonal, which is zero by construction in the distance
            // modes and would only pull the summary towards it.
            if mode == DistMode::Topology {
                eprintln!("Minimum:           {}", totals.min as i64);
                eprintln!("Maximum:           {}", totals.max as i64);
                eprintln!("Mean:              {:.3}", totals.mean());
            } else {
                eprintln!("Minimum:           {:.*}", precision, totals.min);
                eprintln!("Maximum:           {:.*}", precision, totals.max);
                eprintln!("Mean:              {:.*}", precision, totals.mean());
            }
        }
        eprintln!("Time:              {:.3}s", started.elapsed().as_secs_f64());
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

/// Build the header of a NumPy `.npy` file (format 1.0) for a little-endian
/// f64 array of the given shape.
///
/// Magic, version, then a two-byte header length and an ASCII dict, the whole
/// thing padded to a multiple of 64 bytes so the array data lands aligned.
fn npy_header(shape: &str) -> Vec<u8> {
    let dict = format!(
        "{{'descr': '<f8', 'fortran_order': False, 'shape': {}, }}",
        shape
    );

    const PREAMBLE: usize = 6 + 2 + 2; // magic + version + header length
    let unpadded = PREAMBLE + dict.len() + 1; // + the closing newline
    let padding = (64 - unpadded % 64) % 64;
    let header_len = dict.len() + padding + 1;

    let mut out = Vec::with_capacity(PREAMBLE + header_len);
    out.extend_from_slice(b"\x93NUMPY");
    out.extend_from_slice(&[0x01, 0x00]);
    out.extend_from_slice(&(header_len as u16).to_le_bytes());
    out.extend_from_slice(dict.as_bytes());
    out.extend(std::iter::repeat_n(b' ', padding));
    out.push(b'\n');
    out
}

/// Running summary of the off-diagonal cells, for `--stats`.
///
/// Accumulated per row by the worker that computed it, then merged in the
/// writing loop, so the numbers cost one pass and no second traversal.
#[derive(Clone, Copy)]
struct Stats {
    min: f64,
    max: f64,
    sum: f64,
    count: u64,
}

impl Default for Stats {
    fn default() -> Self {
        Stats { min: f64::INFINITY, max: f64::NEG_INFINITY, sum: 0.0, count: 0 }
    }
}

impl Stats {
    #[inline]
    fn add(&mut self, value: f64) {
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
        self.sum += value;
        self.count += 1;
    }

    fn merge(&mut self, other: &Stats) {
        if other.count == 0 {
            return;
        }
        if other.min < self.min {
            self.min = other.min;
        }
        if other.max > self.max {
            self.max = other.max;
        }
        self.sum += other.sum;
        self.count += other.count;
    }

    fn mean(&self) -> f64 {
        if self.count == 0 {
            f64::NAN
        } else {
            self.sum / self.count as f64
        }
    }
}

/// Read the leaf labels for `--taxa`: one per line, blanks and `#` comments
/// skipped, duplicates collapsed.
fn read_taxa_list(path: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read the taxa list '{}': {}", path, e))?;

    let mut seen = HashSet::new();
    let mut wanted = Vec::new();
    for line in text.lines() {
        // Trailing \r survives a file written on Windows and would stop every
        // label matching.
        let label = line.trim_end_matches(['\r', '\n']).trim();
        if label.is_empty() || label.starts_with('#') {
            continue;
        }
        if seen.insert(label.to_string()) {
            wanted.push(label.to_string());
        }
    }

    if wanted.is_empty() {
        return Err(format!("The taxa list '{}' is empty.", path).into());
    }
    Ok(wanted)
}

/// Powers of ten that are exact as f64, which is all of them up to 1e22.
const POW10: [f64; 16] = [
    1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
];

/// Write `value` with `precision` decimals, or return false and write nothing.
///
/// `write!("{:.p$}")` goes through `core::fmt` and an exact decimal expansion
/// of the float, which is correct and slow: it was around half the cost of a
/// text run. Scaling by a power of ten and emitting the digits by hand is
/// roughly seven times faster.
///
/// It is only safe where the answer is not in doubt. Multiplying by 10^p
/// rounds once, by at most `|scaled| * 2^-53`, and if that is enough to carry
/// the product across the nearest `.5` boundary then which way the exact value
/// rounds is unknowable from here. Those cases return false and go to the
/// formatter, which is why this can be fast without being approximate.
#[inline]
fn write_fixed_fast(out: &mut Vec<u8>, value: f64, precision: usize) -> bool {
    if precision >= POW10.len() || !value.is_finite() {
        return false;
    }
    let scaled = value * POW10[precision];
    // Past 2^53 the integers are no longer consecutive, so the digits below
    // would be invented.
    if scaled.abs() >= 9.007_199_254_740_992e15 {
        return false;
    }

    let distance_from_tie = ((scaled - scaled.floor()) - 0.5).abs();
    let rounding_error = scaled.abs() * f64::EPSILON + f64::MIN_POSITIVE;
    if distance_from_tie <= rounding_error {
        return false;
    }

    let rounded = scaled.round() as i64;
    // The sign comes from the input rather than the rounded integer, so a small
    // negative value that rounds to zero still prints as "-0.00", matching the
    // standard formatter.
    let negative = value.is_sign_negative();

    let mut buf = [0u8; 32];
    let mut at = buf.len();
    let mut digits = rounded.unsigned_abs();
    for _ in 0..precision {
        at -= 1;
        buf[at] = b'0' + (digits % 10) as u8;
        digits /= 10;
    }
    if precision > 0 {
        at -= 1;
        buf[at] = b'.';
    }
    if digits == 0 {
        at -= 1;
        buf[at] = b'0';
    }
    while digits > 0 {
        at -= 1;
        buf[at] = b'0' + (digits % 10) as u8;
        digits /= 10;
    }
    if negative {
        at -= 1;
        buf[at] = b'-';
    }
    out.extend_from_slice(&buf[at..]);
    true
}

/// Append a single distance value to a row buffer.
///
/// Topology mode outputs integers; patristic and LMM use `precision` decimal
/// places. Writing into a `Vec<u8>` cannot fail, which is what lets this run
/// inside the parallel row builder.
#[inline]
fn format_distance(out: &mut Vec<u8>, dist: f64, mode: DistMode, precision: usize) {
    if mode == DistMode::Topology {
        write!(out, "{}", dist as i64).expect("writing to a Vec<u8> cannot fail");
        return;
    }
    if write_fixed_fast(out, dist, precision) {
        return;
    }
    write!(out, "{:.prec$}", dist, prec = precision).expect("writing to a Vec<u8> cannot fail");
}

#[cfg(test)]
mod tests {
    use super::*;

    /// xorshift64, so a failure is reproducible.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// The fast path either produces exactly what the standard formatter would,
    /// or declines. Anything else is a silently wrong number in the matrix.
    fn assert_agrees(value: f64, precision: usize) {
        let mut fast = Vec::new();
        if !write_fixed_fast(&mut fast, value, precision) {
            return; // declined, so the formatter handles it
        }
        let slow = format!("{:.prec$}", value, prec = precision);
        assert_eq!(
            String::from_utf8(fast).unwrap(),
            slow,
            "fast path disagrees for {:e} at {} decimals",
            value,
            precision
        );
    }

    #[test]
    fn test_fast_formatter_matches_std_on_random_values() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for _ in 0..200_000 {
            let bits = rng.next_u64();
            let magnitude = (bits % 24) as i32 - 12;
            let sign = if bits & 1 == 0 { 1.0 } else { -1.0 };
            let value = rng.unit() * 10f64.powi(magnitude) * sign;
            assert_agrees(value, (rng.next_u64() % 16) as usize);
        }
    }

    #[test]
    fn test_fast_formatter_matches_std_near_rounding_boundaries() {
        // Exact halves are where a scaled-integer shortcut goes wrong, so walk
        // straight down the middle of them at every precision.
        for precision in 0..16usize {
            let step = 10f64.powi(-(precision as i32));
            for i in 0..4_000u64 {
                for delta in [0.0, 0.5, -0.5, 0.499_999_999, 0.500_000_001, 1.0 / 3.0] {
                    assert_agrees((i as f64 + delta) * step, precision);
                    assert_agrees(-(i as f64 + delta) * step, precision);
                }
            }
        }
    }

    #[test]
    fn test_fast_formatter_matches_std_on_distance_shaped_values() {
        // What distree actually emits: small positive distances, at the
        // precisions people use.
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        for _ in 0..100_000 {
            let value = rng.unit() * 1e-4;
            for precision in [0usize, 3, 6, 9, 10, 12, 15] {
                assert_agrees(value, precision);
                assert_agrees(value * 4_411_532.0, precision);
            }
        }
    }

    #[test]
    fn test_fast_formatter_edge_cases() {
        let awkward = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            0.5,
            -0.5,
            0.05,
            0.005,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            5e-324, // the smallest subnormal
            1e-300,
            1e300,
            9.007_199_254_740_992e15, // 2^53
            9.007_199_254_740_991e15,
            0.1,
            0.2,
            0.3,
            1.0 / 3.0,
            2.0 / 3.0,
            1e15,
            1e16,
            123_456.789_012_345,
        ];
        for &value in &awkward {
            for precision in 0..=MAX_PRECISION {
                assert_agrees(value, precision);
            }
        }
    }

    #[test]
    fn test_fast_formatter_declines_rather_than_guesses() {
        // Anything it cannot do exactly must be refused, not approximated
        let mut out = Vec::new();
        assert!(!write_fixed_fast(&mut out, 1.0, 16), "precision past the table");
        assert!(!write_fixed_fast(&mut out, f64::NAN, 6), "NaN");
        assert!(!write_fixed_fast(&mut out, f64::INFINITY, 6), "infinity");
        assert!(!write_fixed_fast(&mut out, 1e300, 6), "past 2^53 once scaled");
        assert!(out.is_empty(), "a declined value must write nothing");
    }
}
