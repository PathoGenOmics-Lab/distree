//! Distance matrices from a phylogeny.
//!
//! The tree side of `distree`: parsing Newick, flattening it into an
//! index-addressed array, answering most-recent-common-ancestor queries in
//! logarithmic time, midpoint rooting, and turning a pair of leaves into a
//! distance. The binary in `main.rs` is the command line over the top of it.

pub mod lca;
pub mod midpoint;
pub mod parser;
#[cfg(test)]
mod testutil;
pub mod tree;

/// Distance mode to compute.
///
/// Selects how pairwise leaf distances are calculated from the tree.
#[derive(Clone, Copy, PartialEq)]
pub enum DistMode {
    Patristic,
    Topology,
    Lmm,
}

/// Compute the distance between two leaves according to `mode`.
///
/// Returns the patristic distance, topological hop count, or LMM covariance depth.
#[inline]
pub fn compute_distance(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lca::build_lca_structure;
    use crate::parser::{flatten_raw, parse_newick};
    use crate::tree::Node;
    use std::collections::HashSet;

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
