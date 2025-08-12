use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::str::Chars;
use std::iter::Peekable;
use std::num::ParseFloatError;

/// Our internal representation of a tree node:
/// - `name`: if this is a leaf, its label; otherwise `None`.
/// - `length`: branch length from its parent (0.0 if unspecified).
/// - `parent`: index of the parent node in `nodes`, or `None` if this is the root.
/// - `children`: list of indices of child nodes.
struct Node {
    name: Option<String>,
    length: f64,
    parent: Option<usize>,
    children: Vec<usize>,
}

/// Precomputed data for LCA (binary lifting):
/// - `up[k][u]`: the 2^k-th ancestor of node `u`, or `None` if above the root.
/// - `depth_len[u]`: distance (sum of branch lengths) from node `u` up to the root.
/// - `depth_top[u]`: number of edges from node `u` up to the root.
struct LcaData {
    up: Vec<Vec<Option<usize>>>,
    depth_len: Vec<f64>,
    depth_top: Vec<usize>,
}

/// A temporary structure used during Newick parsing,
/// before flattening into `Vec<Node>`.
struct RawNode {
    name: Option<String>,
    length: f64,
    children: Vec<RawNode>,
}

fn main() -> io::Result<()> {
    // ========== 1. Parse command-line arguments ==========
    let matches = Command::new("distree")
        .version("1.0.0")
        .author("Paula Ruiz-Rodriguez")
        .about("Extracts a distance matrix from a phylogeny (parallel, low-memory)")
        .arg(
            Arg::new("phylogeny")
                .help("Path to the tree file in Newick format")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .help("Tree file format (only 'newick' is supported)")
                .default_value("newick")
                .value_parser(["newick"]),
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
        .get_matches();

    let tree_path = matches
        .get_one::<String>("phylogeny")
        .expect("Tree file path is required")
        .to_string();
    let _format = matches.get_one::<String>("format").unwrap().as_str(); // only "newick"
    let do_midpoint = *matches.get_one::<bool>("midpoint").unwrap();
    let do_lmm = *matches.get_one::<bool>("lmm").unwrap();
    let do_topology = *matches.get_one::<bool>("topology").unwrap();
    let output_path = matches.get_one::<String>("output");

    // Prepare writer: either a file or stdout
    let mut writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(io::stdout())
    };

    // ========== 2. Read the entire Newick file into a String ==========
    let mut newick_str = String::new();
    {
        let mut f = File::open(&tree_path)?;
        f.read_to_string(&mut newick_str)?;
    }

    // ========== 3. Parse the Newick string with our minimal parser ==========
    let mut chars = newick_str.trim().chars().peekable();
    let raw_root = parse_subtree(&mut chars)
        .expect("Failed to parse Newick tree");
    // After the tree, there may be a trailing semicolon
    if let Some(&c) = chars.peek() {
        if c == ';' {
            chars.next();
        }
    }

    // ========== 4. Flatten RawNode into Vec<Node> ==========
    let mut nodes: Vec<Node> = Vec::new();
    let root_idx = flatten_raw(&raw_root, None, &mut nodes);

    // ========== 5. Midpoint-root optional ==========
    if do_midpoint {
        midpoint_root(root_idx, &mut nodes);
    }

    // ========== 6. Build the LCA data structure ==========
    let lca_data = build_lca_structure(root_idx, &nodes);

    // ========== 7. Collect all leaves (taxa) ==========
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

    // Build Vec<(label, index)> and sort by label
    let mut leaf_label_pairs: Vec<(String, usize)> = leaf_indices
        .iter()
        .map(|&i| (nodes[i].name.clone().unwrap(), i))
        .collect();
    leaf_label_pairs.sort_by(|a, b| a.0.cmp(&b.0));

    let sorted_labels: Vec<String> = leaf_label_pairs.iter().map(|(lab, _)| lab.clone()).collect();
    let sorted_leaf_indices: Vec<usize> =
        leaf_label_pairs.iter().map(|(_, idx)| *idx).collect();
    let n_leaves = sorted_leaf_indices.len();

    // ========== 8. Print the TSV header row ==========
    writer.write_all(b"\t")?;
    for (i, lab) in sorted_labels.iter().enumerate() {
        writer.write_all(lab.as_bytes())?;
        if i + 1 < n_leaves {
            writer.write_all(b"\t")?;
        }
    }
    writer.write_all(b"\n")?;

    // ========== 9. For each leaf (row), compute distances to all leaves in parallel ==========
    for (row_i, &leaf_i) in sorted_leaf_indices.iter().enumerate() {
        let this_row: Vec<f64> = sorted_leaf_indices
            .par_iter()
            .map(|&leaf_j| {
                if do_lmm {
                    // LMM: depth of MRCA in branch-length units
                    let m = lca_data.mrca(leaf_i, leaf_j);
                    lca_data.depth_len[m]
                } else if do_topology {
                    // Topological distance = edge count
                    let m = lca_data.mrca(leaf_i, leaf_j);
                    let d_i = lca_data.depth_top[leaf_i];
                    let d_j = lca_data.depth_top[leaf_j];
                    let d_m = lca_data.depth_top[m];
                    ((d_i + d_j).saturating_sub(2 * d_m)) as f64
                } else {
                    // Patristic distance = sum of branch lengths
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
            writer.write_all(format!("{}", dist).as_bytes())?;
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

/// Recursively parse a Newick subtree and return a `RawNode`.
fn parse_subtree(chars: &mut Peekable<Chars>) -> Result<RawNode, String> {
    let mut node = RawNode {
        name: None,
        length: 0.0,
        children: Vec::new(),
    };

    // If the next character is '(', this is an internal node with children
    if let Some(&c) = chars.peek() {
        if c == '(' {
            // Consume '('
            chars.next();

            // Parse each child until we see a ')'
            loop {
                let child = parse_subtree(chars)?;
                node.children.push(child);
                match chars.peek() {
                    Some(',') => {
                        chars.next();
                        continue;
                    }
                    Some(')') => {
                        chars.next();
                        break;
                    }
                    Some(other) => {
                        return Err(format!("Expected ',' or ')', found '{}'", other));
                    }
                    None => return Err("Unexpected end of input after '('".to_string()),
                }
            }

            // After ')', we may have an internal node label
            if let Some(&c2) = chars.peek() {
                if c2 != ':' && c2 != ',' && c2 != ')' && c2 != ';' {
                    let name = parse_label(chars);
                    if !name.is_empty() {
                        node.name = Some(name);
                    }
                }
            }

            // Then optionally a colon followed by a branch length
            if let Some(&':') = chars.peek() {
                chars.next();
                let length = parse_length(chars)?;
                node.length = length;
            }

            return Ok(node);
        }
    }

    // Otherwise, it's a leaf: parse the label
    let name = parse_label(chars);
    if !name.is_empty() {
        node.name = Some(name);
    }

    // Optionally, a colon followed by a branch length
    if let Some(&':') = chars.peek() {
        chars.next();
        let length = parse_length(chars)?;
        node.length = length;
    }

    Ok(node)
}

/// Parse a node label (until we hit ':', ',', ')', ';', or whitespace).
fn parse_label(chars: &mut Peekable<Chars>) -> String {
    let mut label = String::new();
    while let Some(&c) = chars.peek() {
        if c == ':' || c == ',' || c == ')' || c == ';' || c.is_whitespace() {
            break;
        }
        label.push(c);
        chars.next();
    }
    label
}

/// Parse a floating-point branch length.
fn parse_length(chars: &mut Peekable<Chars>) -> Result<f64, String> {
    let mut numstr = String::new();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-' {
            numstr.push(c);
            chars.next();
        } else {
            break;
        }
    }
    numstr
        .parse::<f64>()
        .map_err(|e: ParseFloatError| format!("Failed to parse branch length '{}': {}", numstr, e))
}

/// Recursively flatten a `RawNode` into a flat `Vec<Node>`, returning the index of the newly-added node.
fn flatten_raw(raw: &RawNode, parent: Option<usize>, nodes: &mut Vec<Node>) -> usize {
    let idx = nodes.len();
    nodes.push(Node {
        name: raw.name.clone(),
        length: raw.length,
        parent,
        children: Vec::new(),
    });
    if let Some(p) = parent {
        nodes[p].children.push(idx);
    }
    for child in &raw.children {
        flatten_raw(child, Some(idx), nodes);
    }
    idx
}

/// Midpoint-rooting: find two leaves at the ends of the diameter (maximum branch-length distance),
/// reconstruct that path, locate its midpoint, and insert a new root node there.
fn midpoint_root(root_idx: usize, nodes: &mut Vec<Node>) {
    // Helper: from `start`, do a DFS to find the farthest leaf and its distance.
    fn farthest_from(start: usize, nodes: &Vec<Node>) -> (usize, f64, Vec<Option<usize>>) {
        let mut best_leaf = start;
        let mut best_dist = 0.0;
        let mut parent_trace = vec![None; nodes.len()];
        let mut visited = vec![false; nodes.len()];
        let mut stack = vec![(start, 0.0)];
        visited[start] = true;

        while let Some((u, dist_u)) = stack.pop() {
            if nodes[u].children.is_empty() {
                // If it's a leaf, check distance
                if dist_u > best_dist {
                    best_dist = dist_u;
                    best_leaf = u;
                }
            }
            // Visit the parent, if any
            if let Some(p) = nodes[u].parent {
                if !visited[p] {
                    visited[p] = true;
                    parent_trace[p] = Some(u);
                    stack.push((p, dist_u + nodes[u].length));
                }
            }
            // Visit children
            for &v in &nodes[u].children {
                if !visited[v] {
                    visited[v] = true;
                    parent_trace[v] = Some(u);
                    stack.push((v, dist_u + nodes[v].length));
                }
            }
        }
        (best_leaf, best_dist, parent_trace)
    }

    // 1) Find any leaf to start (descend to a leaf if root isn't one)
    let mut any_leaf = root_idx;
    if !nodes[any_leaf].children.is_empty() {
        let mut cur = any_leaf;
        while !nodes[cur].children.is_empty() {
            cur = nodes[cur].children[0];
        }
        any_leaf = cur;
    }

    // 2) From any_leaf, find leaf_a (one endpoint of the diameter)
    let (leaf_a, _, _) = farthest_from(any_leaf, nodes);
    // 3) From leaf_a, find leaf_b and the total distance dist_ab
    let (leaf_b, dist_ab, parent_trace_b) = farthest_from(leaf_a, nodes);

    // Reconstruct the path from leaf_b back to leaf_a
    let mut path = Vec::new();
    {
        let mut cur = leaf_b;
        loop {
            path.push(cur);
            if cur == leaf_a {
                break;
            }
            cur = parent_trace_b[cur].unwrap();
        }
    }
    // The midpoint lies at dist_ab/2 along that path (starting from leaf_b)
    let half = dist_ab / 2.0;
    let mut accum = 0.0;
    let mut midpoint_node = path[0];

    for i in 0..path.len() - 1 {
        let u = path[i];
        let v = path[i + 1];
        // Edge length between u ↔ v
        let edge_len = if nodes[v].children.contains(&u) {
            nodes[u].length
        } else {
            nodes[v].length
        };
        if accum + edge_len >= half {
            // The midpoint falls on this edge (u–v)
            let dist_into_edge = half - accum;
            // Determine which is parent and which is child in the current tree
            let (parent_node, child_node, parent_to_child_len) =
                if nodes[v].parent == Some(u) {
                    (u, v, nodes[v].length)
                } else {
                    (v, u, nodes[u].length)
                };
            let new_root_idx = nodes.len();
            let dist_to_old_parent = parent_to_child_len - dist_into_edge;
            let dist_to_old_child = dist_into_edge;

            // Create the new root node R
            nodes.push(Node {
                name: None,
                length: 0.0,
                parent: None,
                children: Vec::new(),
            });

            // 1) Detach child_node from parent_node, attach child_node → R
            let old_parent_of_child = nodes[child_node].parent.take();
            assert_eq!(old_parent_of_child, Some(parent_node));
            nodes[child_node].parent = Some(new_root_idx);

            // 2) In parent_node.children, replace child_node with new_root_idx
            if let Some(pos) = nodes[parent_node]
                .children
                .iter()
                .position(|&x| x == child_node)
            {
                nodes[parent_node].children[pos] = new_root_idx;
            }

            // 3) Now set R.children = [child_node], and R.parent = Some(parent_node)
            nodes[new_root_idx].parent = Some(parent_node);
            nodes[new_root_idx].children.push(child_node);

            // 4) Adjust branch lengths:
            //    - child_node.length = dist_to_old_child
            //    - R.length = dist_to_old_parent
            nodes[child_node].length = dist_to_old_child;
            nodes[new_root_idx].length = dist_to_old_parent;

            // 5) Insert parent_node above R
            let old_parent_of_parent = nodes[parent_node].parent.take();
            nodes[parent_node].parent = Some(new_root_idx);
            nodes[new_root_idx].children.push(parent_node);
            if let Some(grand) = old_parent_of_parent {
                // If parent_node wasn't originally the root, replace in its old parent's children
                if let Some(pos2) = nodes[grand].children.iter().position(|&x| x == parent_node) {
                    nodes[grand].children[pos2] = new_root_idx;
                }
            }
            // Finally, R must become the new root
            nodes[new_root_idx].parent = None;

            midpoint_node = new_root_idx;
            break;
        }
        accum += edge_len;
        midpoint_node = v;
    }

    // If for some reason we never split (degenerate case), force midpoint_node to be root
    if let Some(p) = nodes[midpoint_node].parent {
        if let Some(pos) = nodes[p].children.iter().position(|&x| x == midpoint_node) {
            nodes[p].children.remove(pos);
        }
        nodes[midpoint_node].parent = None;
    }
}

/// Build the `LcaData` for binary lifting:
/// 1) Do a single DFS from `root_idx` to fill `depth_len`, `depth_top`, and `up[0][u]`.
/// 2) Fill `up[k][u] = up[k-1][ up[k-1][u] ]` for k = 1..⌈log₂(n)⌉.
fn build_lca_structure(root_idx: usize, nodes: &Vec<Node>) -> LcaData {
    let n = nodes.len();
    let mut depth_len = vec![0.0; n];
    let mut depth_top = vec![0; n];
    let max_log = ((n as f64).log2().ceil() as usize) + 1;
    let mut up: Vec<Vec<Option<usize>>> = vec![vec![None; n]; max_log];

    // 1) DFS to set depth and immediate parent
    {
        let mut stack = vec![root_idx];
        depth_len[root_idx] = 0.0;
        depth_top[root_idx] = 0;
        up[0][root_idx] = None;

        while let Some(u) = stack.pop() {
            for &v in &nodes[u].children {
                up[0][v] = Some(u);
                depth_len[v] = depth_len[u] + nodes[v].length;
                depth_top[v] = depth_top[u] + 1;
                stack.push(v);
            }
        }
    }

    // 2) Build all 2^k ancestors
    for k in 1..max_log {
        for u in 0..n {
            up[k][u] = up[k - 1][u].and_then(|mid| up[k - 1][mid]);
        }
    }

    LcaData {
        up,
        depth_len,
        depth_top,
    }
}

impl LcaData {
    /// Return the index of the lowest common ancestor (LCA) of nodes `u` and `v` in O(log n).
    fn mrca(&self, mut u: usize, mut v: usize) -> usize {
        if u == v {
            return u;
        }
        // 1) Lift the deeper node up until both are at the same depth_top
        if self.depth_top[u] < self.depth_top[v] {
            std::mem::swap(&mut u, &mut v);
        }
        let diff = self.depth_top[u] - self.depth_top[v];
        let mut x = diff;
        let mut k = 0;
        while x > 0 {
            if (x & 1) == 1 {
                u = self.up[k][u].unwrap();
            }
            x >>= 1;
            k += 1;
        }
        if u == v {
            return u;
        }
        // 2) Lift both in powers of two until their parents differ
        for k in (0..self.up.len()).rev() {
            if let (Some(au), Some(av)) = (self.up[k][u], self.up[k][v]) {
                if au != av {
                    u = au;
                    v = av;
                }
            }
        }
        // Now u and v have the same parent, which is the LCA
        self.up[0][u].unwrap()
    }
}
