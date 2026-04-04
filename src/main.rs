mod lca;
mod midpoint;
mod parser;
mod tree;

use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};

use lca::build_lca_structure;
use midpoint::midpoint_root;
use parser::{flatten_raw, parse_subtree};
use tree::Node;

fn main() -> io::Result<()> {
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
        .get_matches();

    let tree_path = matches
        .get_one::<String>("phylogeny")
        .expect("Tree file path is required")
        .to_string();
    let do_midpoint = *matches.get_one::<bool>("midpoint").unwrap();
    let do_lmm = *matches.get_one::<bool>("lmm").unwrap();
    let do_topology = *matches.get_one::<bool>("topology").unwrap();
    let output_path = matches.get_one::<String>("output");
    let precision = *matches.get_one::<usize>("precision").unwrap();

    // Configure thread pool
    if let Some(&num_threads) = matches.get_one::<usize>("threads") {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    let mut writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(io::stdout())
    };

    // Read input from file or stdin
    let mut newick_str = String::new();
    if tree_path == "-" {
        io::stdin().read_to_string(&mut newick_str)?;
    } else {
        let mut f = File::open(&tree_path)?;
        f.read_to_string(&mut newick_str)?;
    }

    // Parse the Newick string
    let mut chars = newick_str.trim().chars().peekable();
    let raw_root = parse_subtree(&mut chars).expect("Failed to parse Newick tree");
    if let Some(&c) = chars.peek() {
        if c == ';' {
            chars.next();
        }
    }

    // Flatten into Vec<Node>
    let mut nodes: Vec<Node> = Vec::new();
    let mut root_idx = flatten_raw(&raw_root, None, &mut nodes);

    // Warn about negative branch lengths
    let has_negative = nodes.iter().any(|n| n.length < 0.0);
    if has_negative {
        eprintln!("Warning: negative branch lengths detected in the tree.");
    }

    // Midpoint-root if requested
    if do_midpoint {
        root_idx = midpoint_root(root_idx, &mut nodes);
    }

    // Build LCA
    let lca_data = build_lca_structure(root_idx, &nodes);

    // Collect leaves
    let mut leaf_indices: Vec<usize> = Vec::new();
    for (i, nd) in nodes.iter().enumerate() {
        if nd.children.is_empty() && nd.name.is_some() {
            leaf_indices.push(i);
        }
    }
    if leaf_indices.is_empty() {
        eprintln!("No labeled leaves found in the tree. Exiting.");
        std::process::exit(1);
    }

    // Check for duplicate leaf names
    {
        let mut seen = HashSet::new();
        for &i in &leaf_indices {
            let name = nodes[i].name.as_ref().unwrap();
            if !seen.insert(name.clone()) {
                eprintln!("Error: duplicate leaf name '{}' found. Leaf names must be unique.", name);
                std::process::exit(1);
            }
        }
    }

    // Sort by label
    let mut leaf_label_pairs: Vec<(String, usize)> = leaf_indices
        .iter()
        .map(|&i| (nodes[i].name.clone().unwrap(), i))
        .collect();
    leaf_label_pairs.sort_by(|a, b| a.0.cmp(&b.0));

    let sorted_labels: Vec<String> = leaf_label_pairs.iter().map(|(lab, _)| lab.clone()).collect();
    let sorted_leaf_indices: Vec<usize> =
        leaf_label_pairs.iter().map(|(_, idx)| *idx).collect();
    let n_leaves = sorted_leaf_indices.len();

    // Warn if no branch lengths and not topology mode
    if !do_topology {
        let all_zero = nodes.iter().all(|n| n.length == 0.0);
        if all_zero {
            eprintln!(
                "Warning: no branch lengths detected, all patristic distances will be zero. Consider --topology."
            );
        }
    }

    // Print TSV header (leading tab for R/Python compatibility)
    writer.write_all(b"\t")?;
    for (i, lab) in sorted_labels.iter().enumerate() {
        writer.write_all(lab.as_bytes())?;
        if i + 1 < n_leaves {
            writer.write_all(b"\t")?;
        }
    }
    writer.write_all(b"\n")?;

    // Compute and print distance matrix
    let prec = precision;
    for (row_i, &leaf_i) in sorted_leaf_indices.iter().enumerate() {
        let this_row: Vec<f64> = sorted_leaf_indices
            .par_iter()
            .map(|&leaf_j| {
                if do_lmm {
                    let m = lca_data.mrca(leaf_i, leaf_j);
                    lca_data.depth_len[m]
                } else if do_topology {
                    let m = lca_data.mrca(leaf_i, leaf_j);
                    let d_i = lca_data.depth_top[leaf_i];
                    let d_j = lca_data.depth_top[leaf_j];
                    let d_m = lca_data.depth_top[m];
                    ((d_i + d_j).saturating_sub(2 * d_m)) as f64
                } else {
                    let m = lca_data.mrca(leaf_i, leaf_j);
                    let d_i = lca_data.depth_len[leaf_i];
                    let d_j = lca_data.depth_len[leaf_j];
                    let d_m = lca_data.depth_len[m];
                    d_i + d_j - 2.0 * d_m
                }
            })
            .collect();

        writer.write_all(sorted_labels[row_i].as_bytes())?;
        for dist in this_row.iter() {
            writer.write_all(b"\t")?;
            write!(writer, "{:.prec$}", dist, prec = prec)?;
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_tree(newick: &str) -> (Vec<Node>, usize) {
        let mut chars = newick.trim().chars().peekable();
        let raw = parse_subtree(&mut chars).unwrap();
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
        let m = lca.mrca(ai, bi);
        lca.depth_len[ai] + lca.depth_len[bi] - 2.0 * lca.depth_len[m]
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
        let m = lca.mrca(a, c);
        let topo = lca.depth_top[a] + lca.depth_top[c] - 2 * lca.depth_top[m];
        // A is depth 2, C is depth 1, root is depth 0 => 2+1-0 = 3
        assert_eq!(topo, 3);
    }

    #[test]
    fn test_lmm_matrix() {
        let (nodes, root) = build_tree("((A:1.0,B:2.0):0.5,C:3.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = get_leaf(&nodes, "A");
        let b = get_leaf(&nodes, "B");
        let c = get_leaf(&nodes, "C");

        // MRCA(A,B) is the inner node at depth 0.5
        let m_ab = lca.mrca(a, b);
        assert!((lca.depth_len[m_ab] - 0.5).abs() < 1e-10);

        // MRCA(A,C) is root at depth 0
        let m_ac = lca.mrca(a, c);
        assert!((lca.depth_len[m_ac]).abs() < 1e-10);
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
        let mut has_dup = false;
        for &i in &leaf_indices {
            let name = nodes[i].name.as_ref().unwrap();
            if !seen.insert(name.clone()) {
                has_dup = true;
                break;
            }
        }
        assert!(has_dup);
    }

    #[test]
    fn test_negative_branch_length_detected() {
        let (nodes, _root) = build_tree("(A:-0.5,B:2.0);");
        let has_negative = nodes.iter().any(|n| n.length < 0.0);
        assert!(has_negative);
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
        let all_zero = nodes.iter().all(|n| n.length == 0.0);
        assert!(all_zero);
    }
}
