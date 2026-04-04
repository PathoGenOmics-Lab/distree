use crate::tree::Node;

/// Precomputed data for O(log n) LCA queries using binary lifting.
pub struct LcaData {
    /// Binary lifting table: `up[k][u]` is the 2^k-th ancestor of node `u`.
    pub up: Vec<Vec<Option<usize>>>,
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
    let max_log = if n <= 1 { 1 } else { ((n as f64).log2().ceil() as usize) + 1 };
    let mut depth_len = vec![0.0; n];
    let mut depth_top = vec![0; n];
    let mut up: Vec<Vec<Option<usize>>> = vec![vec![None; n]; max_log];

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
        for k in (0..self.up.len()).rev() {
            if let (Some(au), Some(av)) = (self.up[k][u], self.up[k][v]) {
                if au != av {
                    u = au;
                    v = av;
                }
            }
        }
        self.up[0][u].expect("LCA: could not find common ancestor; tree may be malformed")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{flatten_raw, parse_newick};

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
