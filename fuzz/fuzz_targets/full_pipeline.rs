#![no_main]

//! Run whatever parses all the way through to a distance.
//!
//! The parser rejecting bad input is only half of it. A tree that parses can
//! still be a shape nothing downstream expects, and flattening, midpoint
//! rooting and the binary-lifting queries all index into flat arrays. This
//! target takes every input the parser accepts and pushes it through the rest
//! of the pipeline, in every mode, with and without rooting.

use libfuzzer_sys::fuzz_target;

use distree::lca::build_lca_structure;
use distree::midpoint::midpoint_root;
use distree::parser::{flatten_raw, parse_newick};
use distree::tree::Node;
use distree::{compute_distance, DistMode};

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    if text.len() > 16 * 1024 {
        return;
    }
    let Ok(raw) = parse_newick(text) else {
        return;
    };

    for root_it in [false, true] {
        let mut nodes: Vec<Node> = Vec::new();
        let mut root = flatten_raw(&raw, None, &mut nodes);

        // Midpoint rooting is undefined with negative branch lengths, and the
        // binary warns rather than refusing, so keep those out of this arm.
        if root_it {
            if nodes.iter().any(|n| n.length < 0.0 || !n.length.is_finite()) {
                continue;
            }
            root = midpoint_root(root, &mut nodes);
        }

        let lca = build_lca_structure(root, &nodes);

        let leaves: Vec<usize> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, _)| i)
            .collect();

        // The matrix is quadratic; a few hundred leaves is plenty to exercise
        // the query paths without starving the fuzzer of executions.
        let sample = leaves.len().min(24);
        for &i in leaves.iter().take(sample) {
            for &j in leaves.iter().take(sample) {
                for mode in [DistMode::Patristic, DistMode::Topology, DistMode::Lmm] {
                    let d = compute_distance(i, j, mode, &lca);
                    if i == j && mode != DistMode::Lmm {
                        assert_eq!(d, 0.0, "a leaf must be exactly 0 from itself");
                    }
                }
            }
        }
    }
});
