use crate::tree::Node;

/// Marks "no such ancestor" in the binary-lifting table.
///
/// `Option<usize>` would be the obvious type, but it has no spare niche and so
/// costs 16 bytes per entry against 8 for a bare `usize`. The table holds
/// n·log n entries, which on a million-leaf tree is the difference between
/// roughly 700 MB and 350 MB.
const NO_ANCESTOR: usize = usize::MAX;

/// Precomputed data for O(log n) LCA queries using binary lifting.
pub struct LcaData {
    /// Binary lifting table, stored flat: `up[k * n + u]` is the 2^k-th
    /// ancestor of node `u`, or [`NO_ANCESTOR`] if `u` has none.
    up: Vec<usize>,
    /// Node count, and therefore the stride of a row of `up`.
    n: usize,
    /// Number of levels held in `up`.
    levels: usize,
    /// Cumulative branch-length distance from the root to each node.
    pub depth_len: Vec<f64>,
    /// Topological depth (number of edges) from the root to each node.
    pub depth_top: Vec<usize>,
}

/// Build the [`LcaData`] binary-lifting structure rooted at `root_idx`.
///
/// Runs in O(n log n) time and O(n log n) space.
pub fn build_lca_structure(root_idx: usize, nodes: &[Node]) -> LcaData {
    let n = nodes.len();
    let levels = if n <= 1 { 1 } else { ((n as f64).log2().ceil() as usize) + 1 };
    let mut depth_len = vec![0.0; n];
    let mut depth_top = vec![0; n];
    let mut up: Vec<usize> = vec![NO_ANCESTOR; levels * n];

    {
        let mut stack = vec![root_idx];
        depth_len[root_idx] = 0.0;
        depth_top[root_idx] = 0;

        while let Some(u) = stack.pop() {
            for &v in &nodes[u].children {
                up[v] = u;
                depth_len[v] = depth_len[u] + nodes[v].length;
                depth_top[v] = depth_top[u] + 1;
                stack.push(v);
            }
        }
    }

    for k in 1..levels {
        for u in 0..n {
            let mid = up[(k - 1) * n + u];
            up[k * n + u] = if mid == NO_ANCESTOR {
                NO_ANCESTOR
            } else {
                up[(k - 1) * n + mid]
            };
        }
    }

    LcaData {
        up,
        n,
        levels,
        depth_len,
        depth_top,
    }
}

impl LcaData {
    /// The 2^k-th ancestor of `u`, or [`NO_ANCESTOR`].
    #[inline]
    fn ancestor(&self, k: usize, u: usize) -> usize {
        self.up[k * self.n + u]
    }

    /// Return the index of the Most Recent Common Ancestor (MRCA) of nodes `u` and `v`.
    ///
    /// Uses binary lifting for O(log n) queries.
    pub fn mrca(&self, mut u: usize, mut v: usize) -> usize {
        if u == v {
            return u;
        }
        if self.depth_top[u] < self.depth_top[v] {
            std::mem::swap(&mut u, &mut v);
        }
        // Lift the deeper node until both sit at the same depth
        let mut diff = self.depth_top[u] - self.depth_top[v];
        let mut k = 0;
        while diff > 0 {
            if (diff & 1) == 1 {
                u = self.ancestor(k, u);
                assert_ne!(
                    u, NO_ANCESTOR,
                    "LCA: depth table disagrees with the tree; tree may be malformed"
                );
            }
            diff >>= 1;
            k += 1;
        }
        if u == v {
            return u;
        }
        // Climb in step as long as the ancestors differ; they meet one above
        for k in (0..self.levels).rev() {
            let au = self.ancestor(k, u);
            let av = self.ancestor(k, v);
            if au != av && au != NO_ANCESTOR && av != NO_ANCESTOR {
                u = au;
                v = av;
            }
        }
        let m = self.ancestor(0, u);
        assert_ne!(
            m, NO_ANCESTOR,
            "LCA: could not find common ancestor; tree may be malformed"
        );
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{flatten_raw, parse_newick};
    use crate::testutil::{random_newick, Rng};
    use std::collections::HashSet;

    fn make(newick: &str) -> (Vec<crate::tree::Node>, usize) {
        let raw = parse_newick(newick).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        (nodes, root)
    }

    #[test]
    fn test_mrca_siblings() {
        // (A:1,B:2); — MRCA(A,B) = root
        let (nodes, root) = make("(A:1.0,B:2.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
        let b = nodes.iter().position(|n| n.name.as_deref() == Some("B")).unwrap();
        assert_eq!(lca.mrca(a, b), root);
    }

    #[test]
    fn test_mrca_deeper() {
        // ((A:1,B:2):3,C:4); — MRCA(A,B) = inner node, MRCA(A,C) = root
        let (nodes, root) = make("((A:1.0,B:2.0):3.0,C:4.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
        let b = nodes.iter().position(|n| n.name.as_deref() == Some("B")).unwrap();
        let c = nodes.iter().position(|n| n.name.as_deref() == Some("C")).unwrap();
        let ab = lca.mrca(a, b);
        assert_ne!(ab, root, "MRCA(A,B) should be the inner node, not root");
        assert_eq!(lca.mrca(a, c), root);
        // depth of inner node = 3.0
        assert!((lca.depth_len[ab] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_mrca_self() {
        // MRCA(x, x) must be x itself
        let (nodes, root) = make("(A:1.0,B:2.0);");
        let lca = build_lca_structure(root, &nodes);
        let a = nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
        assert_eq!(lca.mrca(a, a), a);
        assert_eq!(lca.mrca(root, root), root);
    }

    /// Walk to the root from `u`, then from `v` until the paths meet.
    fn mrca_by_walking(nodes: &[Node], u: usize, v: usize) -> usize {
        let mut seen = HashSet::new();
        let mut cur = Some(u);
        while let Some(x) = cur {
            seen.insert(x);
            cur = nodes[x].parent;
        }
        let mut cur = Some(v);
        while let Some(x) = cur {
            if seen.contains(&x) {
                return x;
            }
            cur = nodes[x].parent;
        }
        panic!("no common ancestor between {} and {}", u, v);
    }

    #[test]
    fn test_mrca_matches_walking_up_on_random_trees() {
        let mut rng = Rng::new(0x2545_F491_4F6C_DD1D);

        for case in 0..200 {
            let n_leaves = 2 + rng.below(30);
            let newick = random_newick(&mut rng, n_leaves);
            let raw = parse_newick(&newick).unwrap();
            let mut nodes = Vec::new();
            let root = flatten_raw(&raw, None, &mut nodes);
            let lca = build_lca_structure(root, &nodes);

            for u in 0..nodes.len() {
                for v in 0..nodes.len() {
                    let expected = mrca_by_walking(&nodes, u, v);
                    assert_eq!(
                        lca.mrca(u, v),
                        expected,
                        "case {}: mrca({}, {}) in {}",
                        case,
                        u,
                        v,
                        newick
                    );
                }
            }
        }
    }

    #[test]
    fn test_depth_top() {
        // ((A,B),C); — A and B are at depth 2, C at depth 1
        let (nodes, root) = make("((A:1,B:1):1,C:1);");
        let lca = build_lca_structure(root, &nodes);
        let a = nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
        let b = nodes.iter().position(|n| n.name.as_deref() == Some("B")).unwrap();
        let c = nodes.iter().position(|n| n.name.as_deref() == Some("C")).unwrap();
        assert_eq!(lca.depth_top[root], 0);
        assert_eq!(lca.depth_top[c], 1);
        assert_eq!(lca.depth_top[a], 2);
        assert_eq!(lca.depth_top[b], 2);
    }
}
