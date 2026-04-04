use crate::tree::Node;

/// Precomputed data for LCA (binary lifting).
pub struct LcaData {
    pub up: Vec<Vec<Option<usize>>>,
    pub depth_len: Vec<f64>,
    pub depth_top: Vec<usize>,
}

/// Build the `LcaData` for binary lifting from `root_idx`.
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
    /// Return the index of the MRCA of nodes `u` and `v` in O(log n).
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
        self.up[0][u].unwrap()
    }
}
