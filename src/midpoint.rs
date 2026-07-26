use crate::tree::Node;

/// Midpoint-root the tree: find the diameter, split the midpoint edge, re-root there.
/// Returns the index of the new root.
pub fn midpoint_root(root_idx: usize, nodes: &mut Vec<Node>) -> usize {
    // BFS/DFS on the unrooted tree to find the farthest leaf from `start`.
    fn farthest_from(start: usize, nodes: &[Node]) -> (usize, f64, Vec<Option<usize>>) {
        let mut best_leaf = start;
        let mut best_dist = 0.0;
        let mut parent_trace: Vec<Option<usize>> = vec![None; nodes.len()];
        let mut visited = vec![false; nodes.len()];
        let mut stack = vec![(start, 0.0)];
        visited[start] = true;

        while let Some((u, dist_u)) = stack.pop() {
            if nodes[u].children.is_empty() && u != start && dist_u > best_dist {
                best_dist = dist_u;
                best_leaf = u;
            }
            if let Some(p) = nodes[u].parent {
                if !visited[p] {
                    visited[p] = true;
                    parent_trace[p] = Some(u);
                    stack.push((p, dist_u + nodes[u].length));
                }
            }
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

    // 1) Find diameter endpoints
    let mut any_leaf = root_idx;
    while !nodes[any_leaf].children.is_empty() {
        any_leaf = nodes[any_leaf].children[0];
    }

    let (leaf_a, _, _) = farthest_from(any_leaf, nodes);
    let (leaf_b, dist_ab, parent_trace) = farthest_from(leaf_a, nodes);

    if dist_ab == 0.0 {
        return root_idx;
    }

    // 2) Reconstruct path from leaf_b back to leaf_a
    let mut path = Vec::new();
    {
        let mut cur = leaf_b;
        loop {
            path.push(cur);
            if cur == leaf_a {
                break;
            }
            cur = parent_trace[cur]
                .expect("midpoint: path from leaf_b to leaf_a is disconnected; tree may be malformed");
        }
    }

    // 3) Walk the path and find the edge where the midpoint falls
    let half = dist_ab / 2.0;
    let mut accum = 0.0;

    for i in 0..path.len() - 1 {
        let u = path[i]; // closer to leaf_b
        let v = path[i + 1]; // closer to leaf_a

        let edge_len = if nodes[u].parent == Some(v) {
            nodes[u].length
        } else {
            nodes[v].length
        };

        if accum + edge_len >= half {
            let dist_to_u = half - accum;
            let dist_to_v = edge_len - dist_to_u;

            // Determine parent-child in the current tree
            let (parent_side, child_side, len_child, len_parent) =
                if nodes[u].parent == Some(v) {
                    // v is parent of u
                    (v, u, dist_to_u, dist_to_v)
                } else {
                    // u is parent of v
                    (u, v, dist_to_v, dist_to_u)
                };

            // Insert new root node on this edge
            let new_root = nodes.len();
            nodes.push(Node {
                name: None,
                length: 0.0,
                parent: None,
                children: Vec::new(),
            });

            // Detach child_side from parent_side
            nodes[parent_side].children.retain(|&x| x != child_side);

            // Attach child_side to new_root
            nodes[child_side].parent = Some(new_root);
            nodes[child_side].length = len_child;
            nodes[new_root].children.push(child_side);

            // Collect path from parent_side to old root
            let mut path_to_old_root = vec![];
            {
                let mut c = parent_side;
                loop {
                    path_to_old_root.push(c);
                    if let Some(p) = nodes[c].parent {
                        c = p;
                    } else {
                        break;
                    }
                }
            }

            // Attach parent_side to new_root
            let old_len = nodes[parent_side].length;
            nodes[parent_side].parent = Some(new_root);
            nodes[parent_side].length = len_parent;
            nodes[new_root].children.push(parent_side);

            // Reverse parent-child relationships from parent_side up to old root
            let mut prev_old_length = old_len;
            for j in 1..path_to_old_root.len() {
                let child_now = path_to_old_root[j];
                let parent_now = path_to_old_root[j - 1];

                // Remove parent_now from child_now's children
                nodes[child_now].children.retain(|&x| x != parent_now);
                // Add child_now as child of parent_now
                nodes[parent_now].children.push(child_now);
                // Set child_now's parent
                nodes[child_now].parent = Some(parent_now);

                // Swap branch lengths
                std::mem::swap(&mut nodes[child_now].length, &mut prev_old_length);
            }

            return new_root;
        }
        accum += edge_len;
    }

    root_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{flatten_raw, parse_newick};
    use crate::lca::build_lca_structure;

    fn patristic_distance(
        leaf_i: usize,
        leaf_j: usize,
        lca: &crate::lca::LcaData,
    ) -> f64 {
        let m = lca.mrca(leaf_i, leaf_j);
        lca.depth_len[leaf_i] + lca.depth_len[leaf_j] - 2.0 * lca.depth_len[m]
    }

    /// xorshift64, so the randomised checks stay reproducible without pulling
    /// in a dependency.
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

        fn below(&mut self, n: usize) -> usize {
            (self.next_u64() % n as u64) as usize
        }

        /// Zero one branch in seven: collapsing weakly supported branches
        /// leaves plenty of them in real trees.
        fn length(&mut self) -> f64 {
            if self.below(7) == 0 {
                return 0.0;
            }
            0.01 + (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 * 2.0
        }
    }

    /// Join random groups of 2 or 3 until one tree remains, so the sample
    /// covers polytomies as well as strictly bifurcating trees.
    fn random_newick(rng: &mut Rng, n_leaves: usize) -> String {
        let mut pool: Vec<String> = (0..n_leaves)
            .map(|i| format!("L{}:{:.6}", i, rng.length()))
            .collect();
        while pool.len() > 1 {
            let arity = if pool.len() > 2 && rng.below(4) == 0 { 3 } else { 2 };
            let children: Vec<String> = (0..arity)
                .map(|_| pool.swap_remove(rng.below(pool.len())))
                .collect();
            pool.push(format!("({}):{:.6}", children.join(","), rng.length()));
        }
        format!("{};", pool.pop().unwrap())
    }

    /// Every node reachable exactly once, every child pointing back at its parent.
    fn assert_tree_is_consistent(nodes: &[Node], root: usize) {
        assert!(nodes[root].parent.is_none(), "the root must have no parent");
        let mut visited = vec![false; nodes.len()];
        let mut stack = vec![root];
        visited[root] = true;
        while let Some(u) = stack.pop() {
            for &v in &nodes[u].children {
                assert_eq!(
                    nodes[v].parent,
                    Some(u),
                    "child {} does not point back at parent {}",
                    v,
                    u
                );
                assert!(!visited[v], "node {} reached twice: the graph has a cycle", v);
                visited[v] = true;
                stack.push(v);
            }
        }
        for (i, seen) in visited.iter().enumerate() {
            assert!(seen, "node {} is unreachable from the root", i);
        }
    }

    fn leaf_indices(nodes: &[Node]) -> Vec<usize> {
        nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, _)| i)
            .collect()
    }

    #[test]
    fn test_midpoint_random_trees() {
        let mut rng = Rng(0x9E3779B97F4A7C15);

        for case in 0..300 {
            let n_leaves = 2 + rng.below(24);
            let newick = random_newick(&mut rng, n_leaves);

            let raw = parse_newick(&newick).unwrap();
            let mut before = Vec::new();
            let root_before = flatten_raw(&raw, None, &mut before);
            let lca_before = build_lca_structure(root_before, &before);
            let leaves_before = leaf_indices(&before);

            // Names in flatten order, so the two trees can be lined up by label
            let names: Vec<String> = leaves_before
                .iter()
                .map(|&i| before[i].name.clone().unwrap())
                .collect();

            let mut diameter: f64 = 0.0;
            let mut dist_before = vec![vec![0.0; leaves_before.len()]; leaves_before.len()];
            for (a, &i) in leaves_before.iter().enumerate() {
                for (b, &j) in leaves_before.iter().enumerate() {
                    let d = patristic_distance(i, j, &lca_before);
                    dist_before[a][b] = d;
                    diameter = diameter.max(d);
                }
            }

            let raw2 = parse_newick(&newick).unwrap();
            let mut after = Vec::new();
            let root_after = flatten_raw(&raw2, None, &mut after);
            let new_root = midpoint_root(root_after, &mut after);

            assert_tree_is_consistent(&after, new_root);

            let lca_after = build_lca_structure(new_root, &after);
            let leaves_after: Vec<usize> = names
                .iter()
                .map(|name| {
                    after
                        .iter()
                        .position(|n| n.name.as_deref() == Some(name.as_str()))
                        .unwrap_or_else(|| panic!("case {}: leaf {} vanished", case, name))
                })
                .collect();

            for a in 0..leaves_after.len() {
                for b in 0..leaves_after.len() {
                    let d = patristic_distance(leaves_after[a], leaves_after[b], &lca_after);
                    assert!(
                        (d - dist_before[a][b]).abs() < 1e-9,
                        "case {}: rooting changed d({},{}) from {} to {}\n{}",
                        case,
                        names[a],
                        names[b],
                        dist_before[a][b],
                        d,
                        newick
                    );
                }
            }

            // The defining property: the two longest root-to-tip paths are
            // equal, so the deepest tip sits exactly half a diameter away.
            let deepest = leaves_after
                .iter()
                .map(|&i| lca_after.depth_len[i])
                .fold(0.0_f64, f64::max);
            assert!(
                (deepest - diameter / 2.0).abs() < 1e-9,
                "case {}: deepest tip is {} from the root, expected {} (diameter {})\n{}",
                case,
                deepest,
                diameter / 2.0,
                diameter,
                newick
            );
        }
    }

    #[test]
    fn test_midpoint_simple() {
        // ((A:1,B:1):1,(C:1,D:1):1);
        // Diameter = 4 (e.g. A→root→D), midpoint at distance 2
        let input = "((A:1,B:1):1,(C:1,D:1):1);";
        let raw = parse_newick(input).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        let new_root = midpoint_root(root, &mut nodes);

        // Verify the tree is valid: new root has no parent
        assert!(nodes[new_root].parent.is_none());

        // Build LCA and check distances are preserved
        let lca = build_lca_structure(new_root, &nodes);
        let leaves: Vec<usize> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, _)| i)
            .collect();

        // A-B should be 2.0, A-C should be 4.0, etc.
        for &i in &leaves {
            for &j in &leaves {
                let d = patristic_distance(i, j, &lca);
                assert!(d >= 0.0, "Negative distance between leaves");
            }
        }

        // Check specific: A to B = 2.0
        let a = leaves.iter().find(|&&i| nodes[i].name.as_deref() == Some("A")).unwrap();
        let b = leaves.iter().find(|&&i| nodes[i].name.as_deref() == Some("B")).unwrap();
        let c = leaves.iter().find(|&&i| nodes[i].name.as_deref() == Some("C")).unwrap();
        assert!((patristic_distance(*a, *b, &lca) - 2.0).abs() < 1e-10);
        assert!((patristic_distance(*a, *c, &lca) - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_midpoint_asymmetric() {
        // (A:1.0,B:3.0);
        // Diameter = 4, midpoint at distance 2 from each
        let input = "(A:1.0,B:3.0);";
        let raw = parse_newick(input).unwrap();
        let mut nodes = Vec::new();
        let root = flatten_raw(&raw, None, &mut nodes);
        let new_root = midpoint_root(root, &mut nodes);

        assert!(nodes[new_root].parent.is_none());
        let lca = build_lca_structure(new_root, &nodes);

        let a = nodes.iter().position(|n| n.name.as_deref() == Some("A")).unwrap();
        let b = nodes.iter().position(|n| n.name.as_deref() == Some("B")).unwrap();
        assert!((patristic_distance(a, b, &lca) - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_midpoint_preserves_distances() {
        let input = "(((A:0.5,B:0.3):0.4,C:0.9):0.1,D:1.2);";
        let raw = parse_newick(input).unwrap();

        // Compute distances before midpoint rooting
        let mut nodes_before = Vec::new();
        let root_before = flatten_raw(&raw, None, &mut nodes_before);
        let lca_before = build_lca_structure(root_before, &nodes_before);

        let leaf_names = ["A", "B", "C", "D"];
        let leaves_before: Vec<usize> = leaf_names
            .iter()
            .map(|name| {
                nodes_before
                    .iter()
                    .position(|n| n.name.as_deref() == Some(name))
                    .unwrap()
            })
            .collect();

        let mut dists_before = vec![vec![0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                dists_before[i][j] = patristic_distance(leaves_before[i], leaves_before[j], &lca_before);
            }
        }

        // Now midpoint root
        let raw2 = parse_newick(input).unwrap();
        let mut nodes_after = Vec::new();
        let root_after = flatten_raw(&raw2, None, &mut nodes_after);
        let new_root = midpoint_root(root_after, &mut nodes_after);
        let lca_after = build_lca_structure(new_root, &nodes_after);

        let leaves_after: Vec<usize> = leaf_names
            .iter()
            .map(|name| {
                nodes_after
                    .iter()
                    .position(|n| n.name.as_deref() == Some(name))
                    .unwrap()
            })
            .collect();

        for i in 0..4 {
            for j in 0..4 {
                let d = patristic_distance(leaves_after[i], leaves_after[j], &lca_after);
                assert!(
                    (d - dists_before[i][j]).abs() < 1e-10,
                    "Distance mismatch for {}-{}: before={}, after={}",
                    leaf_names[i],
                    leaf_names[j],
                    dists_before[i][j],
                    d
                );
            }
        }
    }
}
