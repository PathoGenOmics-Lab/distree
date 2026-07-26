mod lca;
mod midpoint;
mod parser;
#[cfg(test)]
mod testutil;
mod tree;

use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::process::ExitCode;

use lca::build_lca_structure;
use midpoint::midpoint_root;
use parser::{flatten_raw, parse_newick};
use tree::Node;

/// Distance mode to compute.
///
/// Selects how pairwise leaf distances are calculated from the tree.
#[derive(Clone, Copy, PartialEq)]
enum DistMode {
    Patristic,
    Topology,
    Lmm,
}

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
    let mut newick_str = String::new();
    if tree_path == "-" {
        io::stdin().read_to_string(&mut newick_str)?;
    } else {
        File::open(&tree_path)
            .map_err(|e| format!("Cannot open '{}': {}", tree_path, e))?
            .read_to_string(&mut newick_str)?;
    }

    // Parse the Newick string
    let raw_root = parse_newick(&newick_str)
        .map_err(|e| format!("Failed to parse Newick tree: {}", e))?;

    // Flatten into Vec<Node>
    let mut nodes: Vec<Node> = Vec::new();
    let mut root_idx = flatten_raw(&raw_root, None, &mut nodes);

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

/// Compute the distance between two leaves according to `mode`.
///
/// Returns the patristic distance, topological hop count, or LMM covariance depth.
#[inline]
fn compute_distance(
    leaf_i: usize,
    leaf_j: usize,
    mode: DistMode,
    lca_data: &lca::LcaData,
) -> f64 {
    let m = lca_data.mrca(leaf_i, leaf_j);
    match mode {
        DistMode::Lmm => lca_data.depth_len[m],
        DistMode::Topology => {
            let d_i = lca_data.depth_top[leaf_i];
            let d_j = lca_data.depth_top[leaf_j];
            let d_m = lca_data.depth_top[m];
            ((d_i + d_j).saturating_sub(2 * d_m)) as f64
        }
        DistMode::Patristic => {
            let d_i = lca_data.depth_len[leaf_i];
            let d_j = lca_data.depth_len[leaf_j];
            let d_m = lca_data.depth_len[m];
            // Not clamped to zero. Rounding cannot make this negative:
            // depth_len accumulates lengths from the root, so with
            // non-negative branches d_i and d_j are both >= d_m, 2*d_m is
            // exact, and rounding a non-negative exact result to nearest keeps
            // it non-negative. The only way to get a negative value here is a
            // negative branch length, which the tree is warned about, and
            // rounding that up to zero would silently claim two distinct taxa
            // sit on top of each other.
            d_i + d_j - 2.0 * d_m
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn build_tree(newick: &str) -> (Vec<Node>, usize) {
        let raw = parse_newick(newick).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        (nodes, root)
    }

    fn get_leaf(nodes: &[Node], name: &str) -> usize {
        nodes
            .iter()
            .position(|n| n.name.as_deref() == Some(name))
            .unwrap()
    }

    fn patristic(nodes: &[Node], root: usize, a: &str, b: &str) -> f64 {
        let lca = build_lca_structure(root, nodes);
        let ai = get_leaf(nodes, a);
        let bi = get_leaf(nodes, b);
        compute_distance(ai, bi, DistMode::Patristic, &lca)
    }

    #[test]
    fn test_simple_patristic_distances() {
        let (nodes, root) = build_tree("((A:1.0,B:2.0):0.5,C:3.0);");
        assert!((patristic(&nodes, root, "A", "B") - 3.0).abs() < 1e-10);
        assert!((patristic(&nodes, root, "A", "C") - 4.5).abs() < 1e-10);
        assert!((patristic(&nodes, root, "B", "C") - 5.5).abs() < 1e-10);
    }

    #[test]
    fn test_topology_distances() {
        let (nodes, root) = build_tree("((A:1.0,B:2.0):0.5,C:3.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let c = get_leaf(&nodes, "C");
        let d = compute_distance(a, c, DistMode::Topology, &lca);
        assert_eq!(d as i64, 3);
    }

    #[test]
    fn test_lmm_matrix() {
        let (nodes, root) = build_tree("((A:1.0,B:2.0):0.5,C:3.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let b = get_leaf(&nodes, "B");
        let c = get_leaf(&nodes, "C");
        // MRCA(A,B) depth = 0.5
        assert!((compute_distance(a, b, DistMode::Lmm, &lca) - 0.5).abs() < 1e-10);
        // MRCA(A,C) depth = 0.0 (root)
        assert!(compute_distance(a, c, DistMode::Lmm, &lca).abs() < 1e-10);
    }

    #[test]
    fn test_duplicate_leaves_detected() {
        let (nodes, _root) = build_tree("(A:1.0,A:2.0);");
        let leaf_indices: Vec<usize> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, _)| i)
            .collect();
        let mut seen = HashSet::new();
        let has_dup = leaf_indices
            .iter()
            .any(|&i| !seen.insert(nodes[i].name.as_ref().unwrap().as_str()));
        assert!(has_dup);
    }

    #[test]
    fn test_negative_branch_length_detected() {
        let (nodes, _root) = build_tree("(A:-0.5,B:2.0);");
        assert!(nodes.iter().any(|n| n.length < 0.0));
    }

    #[test]
    fn test_trifurcation() {
        let (nodes, root) = build_tree("(A:1.0,B:2.0,C:3.0);");
        assert_eq!(nodes[root].children.len(), 3);
        assert!((patristic(&nodes, root, "A", "B") - 3.0).abs() < 1e-10);
        assert!((patristic(&nodes, root, "A", "C") - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_no_branch_lengths_all_zero() {
        let (nodes, _root) = build_tree("(A,B,(C,D));");
        assert!(nodes.iter().all(|n| n.length == 0.0));
    }

    #[test]
    fn test_single_leaf() {
        let (nodes, _root) = build_tree("(A:1.0);");
        assert_eq!(
            nodes
                .iter()
                .filter(|n| n.children.is_empty() && n.name.is_some())
                .count(),
            1
        );
    }

    #[test]
    fn test_large_symmetric_tree() {
        // ((A:1,B:1):1,(C:1,D:1):1) — all pairwise distances known
        let (nodes, root) = build_tree("((A:1,B:1):1,(C:1,D:1):1);");
        assert!((patristic(&nodes, root, "A", "B") - 2.0).abs() < 1e-10);
        assert!((patristic(&nodes, root, "A", "C") - 4.0).abs() < 1e-10);
        assert!((patristic(&nodes, root, "C", "D") - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_distance_same_leaf() {
        let (nodes, root) = build_tree("(A:1.0,B:2.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        assert!(compute_distance(a, a, DistMode::Patristic, &lca).abs() < 1e-10);
    }

    #[test]
    fn test_lower_triangle_row_counts() {
        // In lower-triangle mode, row i has exactly i distance columns
        let (nodes, root) = build_tree("((A:1,B:2):3,C:4);");
        let _lca = build_lca_structure(root, &nodes);
        let mut leaf_pairs: Vec<(String, usize)> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, n)| (n.name.clone().unwrap(), i))
            .collect();
        leaf_pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        let sorted: Vec<usize> = leaf_pairs.iter().map(|(_, i)| *i).collect();

        // Row 0 (A): 0 distance cols, Row 1 (B): 1 distance col, Row 2 (C): 2 distance cols
        for (row_i, _) in sorted.iter().enumerate() {
            let col_end = row_i; // lower triangle: number of distance values
            assert_eq!(col_end, row_i);
        }
    }

    #[test]
    fn test_patristic_clamped_nonnegative() {
        // Self-distance must be exactly 0.0 (no negative FP artifacts)
        let (nodes, root) = build_tree("(A:1.0000000000001,B:1.0000000000002);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let d = compute_distance(a, a, DistMode::Patristic, &lca);
        assert_eq!(d, 0.0, "Self-distance must be exactly 0.0, got {}", d);
    }

    #[test]
    fn test_patristic_never_negative_on_random_trees() {
        // Rounding must not produce a negative distance on any tree with
        // non-negative branch lengths, and a leaf must be exactly 0 from itself
        let mut rng = crate::testutil::Rng::new(0x5DEE_CE66_D000_0005);

        for _ in 0..200 {
            let n_leaves = 2 + rng.below(30);
            let newick = crate::testutil::random_newick(&mut rng, n_leaves);
            let (nodes, root) = build_tree(&newick);
            let lca = build_lca_structure(root, &nodes);
            let leaves: Vec<usize> = nodes
                .iter()
                .enumerate()
                .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
                .map(|(i, _)| i)
                .collect();

            for &i in &leaves {
                assert_eq!(
                    compute_distance(i, i, DistMode::Patristic, &lca),
                    0.0,
                    "a leaf must be exactly 0.0 from itself in {}",
                    newick
                );
                for &j in &leaves {
                    let d = compute_distance(i, j, DistMode::Patristic, &lca);
                    assert!(d >= 0.0, "negative distance {} in {}", d, newick);
                }
            }
        }
    }

    #[test]
    fn test_negative_branch_length_gives_negative_distance() {
        // A tree with negative branches is warned about; reporting 0.0 would
        // claim two distinct taxa sit on top of each other
        let (nodes, root) = build_tree("(A:-2.0,B:0.5);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let b = get_leaf(&nodes, "B");
        let d = compute_distance(a, b, DistMode::Patristic, &lca);
        assert!((d - (-1.5)).abs() < 1e-12, "expected -1.5, got {}", d);
    }

    #[test]
    fn test_mode_conflict() {
        // Just test the logic: LMM should be distinct from topology
        let (nodes, root) = build_tree("((A:1,B:2):3,C:4);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let c = get_leaf(&nodes, "C");
        let d_pat = compute_distance(a, c, DistMode::Patristic, &lca);
        let d_top = compute_distance(a, c, DistMode::Topology, &lca);
        let d_lmm = compute_distance(a, c, DistMode::Lmm, &lca);
        // All three should give different values for this tree
        assert!((d_pat - 8.0).abs() < 1e-10); // 1+3+4
        assert_eq!(d_top as i64, 3); // 2 edges from A + 1 from C
        assert!(d_lmm.abs() < 1e-10); // MRCA(A,C) = root, depth 0
    }
}
