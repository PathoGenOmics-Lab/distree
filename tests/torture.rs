//! Push malformed and mutated input through the parser and the pipeline.
//!
//! This is the fuzzing that runs on every push. It is not coverage-guided, so
//! it will not find what `cargo fuzz` finds given hours, but it is
//! deterministic, needs no nightly toolchain, and covers the mutations that
//! actually happen to a Newick file: truncation partway through a download, a
//! byte flipped in transfer, a chunk duplicated, a delimiter dropped.
//!
//! The contract under test is narrow and absolute: **any** input at all comes
//! back as `Ok` or `Err`. Never a panic, never an out-of-bounds index, never a
//! loop that does not advance. Whatever parses then has to survive flattening,
//! rooting and every distance mode.
//!
//! For the deeper, coverage-guided version:
//!
//! ```text
//! cargo +nightly fuzz run parse_newick
//! cargo +nightly fuzz run full_pipeline
//! ```

use distree::lca::build_lca_structure;
use distree::midpoint::midpoint_root;
use distree::parser::{flatten_raw, parse_newick};
use distree::tree::Node;
use distree::{compute_distance, DistMode};

/// xorshift64, so a failure is reproducible from the seed in the message.
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
        if n == 0 {
            return 0;
        }
        (self.next_u64() % n as u64) as usize
    }
}

/// The bytes that mean something in Newick, plus a few that do not, so the
/// mutations land on structure rather than on label text.
const ALPHABET: &[u8] = b"(),:;[]'\"0123456789.eE+- \t\nABxy";

fn seed_trees() -> Vec<String> {
    vec![
        "(A:1.0,B:2.0);".to_string(),
        "((A:1.0,B:2.0):0.5,C:3.0);".to_string(),
        "(A:1,B:2,C:3);".to_string(),
        "[&R] ((A:0.1[&&NHX:S=human],B:0.2):0.3,C:0.4);".to_string(),
        "('Taxon A':1.0,\"Taxon B\":2.0)'my clade':0.5;".to_string(),
        "(A:1.5e-3,B:2.0E+1,(C:0.0,D:-0.5):1e10);".to_string(),
        "((((A:1):1):1):1);".to_string(),
        "(it''s:1,B:2);".to_string(),
        "(A,B,(C,D));".to_string(),
        ";".to_string(),
    ]
}

/// Apply one random mutation to `bytes`.
fn mutate(rng: &mut Rng, bytes: &mut Vec<u8>) {
    if bytes.is_empty() {
        bytes.push(ALPHABET[rng.below(ALPHABET.len())]);
        return;
    }
    match rng.below(6) {
        // Truncate: an interrupted write or a partial download
        0 => {
            let at = rng.below(bytes.len());
            bytes.truncate(at);
        }
        // Flip one byte to another meaningful one
        1 => {
            let at = rng.below(bytes.len());
            bytes[at] = ALPHABET[rng.below(ALPHABET.len())];
        }
        // Insert
        2 => {
            let at = rng.below(bytes.len() + 1);
            bytes.insert(at, ALPHABET[rng.below(ALPHABET.len())]);
        }
        // Delete
        3 => {
            let at = rng.below(bytes.len());
            bytes.remove(at);
        }
        // Duplicate a run
        4 => {
            let start = rng.below(bytes.len());
            let end = (start + 1 + rng.below(bytes.len() - start)).min(bytes.len());
            let chunk = bytes[start..end].to_vec();
            let at = rng.below(bytes.len() + 1);
            for (k, b) in chunk.into_iter().enumerate() {
                if bytes.len() < 4096 {
                    bytes.insert(at + k, b);
                }
            }
        }
        // Swap two bytes
        _ => {
            let a = rng.below(bytes.len());
            let b = rng.below(bytes.len());
            bytes.swap(a, b);
        }
    }
}

/// Everything the binary does after a successful parse.
fn exercise_pipeline(text: &str) {
    let Ok(raw) = parse_newick(text) else {
        return;
    };

    for do_midpoint in [false, true] {
        let mut nodes: Vec<Node> = Vec::new();
        let mut root = flatten_raw(&raw, None, &mut nodes);

        if do_midpoint {
            // Midpoint rooting assumes non-negative, finite lengths; the binary
            // warns about the rest rather than refusing.
            if nodes.iter().any(|n| n.length < 0.0 || !n.length.is_finite()) {
                continue;
            }
            root = midpoint_root(root, &mut nodes);
        }

        let lca = build_lca_structure(root, &nodes);

        // The binary refuses a tree whose depths are too large to subtract,
        // because d_i + d_j overflows and inf - inf is NaN. Hold the pipeline
        // to the same precondition rather than asserting past it.
        let deepest = lca.depth_len.iter().copied().fold(0.0_f64, f64::max);
        if !(deepest * 2.0).is_finite() {
            continue;
        }

        let leaves: Vec<usize> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.children.is_empty() && n.name.is_some())
            .map(|(i, _)| i)
            .collect();

        for &i in leaves.iter().take(16) {
            for &j in leaves.iter().take(16) {
                for mode in [DistMode::Patristic, DistMode::Topology, DistMode::Lmm] {
                    let d = compute_distance(i, j, mode, &lca);
                    if i == j && mode != DistMode::Lmm {
                        assert_eq!(d, 0.0, "a leaf must be exactly 0.0 from itself");
                    }
                }
            }
        }
    }
}

#[test]
fn test_mutated_trees_never_panic() {
    let seeds = seed_trees();

    for (s, seed_tree) in seeds.iter().enumerate() {
        let mut rng = Rng(0x243F_6A88_85A3_08D3 ^ (s as u64) << 32);

        for round in 0..4_000 {
            let mut bytes = seed_tree.as_bytes().to_vec();
            // One to four mutations: one finds the shallow cases, four gets far
            // enough from a valid tree to be interesting.
            for _ in 0..=rng.below(4) {
                mutate(&mut rng, &mut bytes);
            }
            // Every byte in ALPHABET is ASCII, so the result is still UTF-8
            let text = match std::str::from_utf8(&bytes) {
                Ok(t) => t,
                Err(_) => continue,
            };

            // A panic here fails the test with the seed and round that caused it
            let result = std::panic::catch_unwind(|| {
                let _ = parse_newick(text);
                exercise_pipeline(text);
            });
            assert!(
                result.is_ok(),
                "panicked on seed {} round {} with input {:?}",
                s,
                round,
                text
            );
        }
    }
}

#[test]
fn test_random_bytes_never_panic() {
    let mut rng = Rng(0x13198A2E_03707344);

    for round in 0..4_000 {
        let len = rng.below(64);
        let bytes: Vec<u8> = (0..len).map(|_| ALPHABET[rng.below(ALPHABET.len())]).collect();
        let text = std::str::from_utf8(&bytes).expect("ALPHABET is ASCII");

        let result = std::panic::catch_unwind(|| {
            let _ = parse_newick(text);
            exercise_pipeline(text);
        });
        assert!(result.is_ok(), "panicked on round {} with input {:?}", round, text);
    }
}

#[test]
fn test_pathological_shapes_never_panic() {
    // Shapes that are valid enough to reach the pipeline but are not what any
    // of it was written with in mind.
    let cases: Vec<String> = vec![
        // A ladder, which makes the LCA depth equal to the node count
        format!("{}A:1{};", "(".repeat(2_000), ")".repeat(2_000)),
        // One node with thousands of children
        format!("({});", (0..2_000).map(|i| format!("L{}:0.1", i)).collect::<Vec<_>>().join(",")),
        // Every branch zero, so the diameter is zero and midpoint has nothing to find
        "((A:0,B:0):0,(C:0,D:0):0);".to_string(),
        // Enormous and tiny lengths in the same tree
        "((A:1e308,B:1e-308):1e300,C:0.0);".to_string(),
        // A single leaf, and a single leaf with no length
        "(A:1.0);".to_string(),
        "A;".to_string(),
        // Internal nodes carrying the labels instead of the tips
        "((:1,:2)inner:3,:4)root;".to_string(),
        // Deeply nested comments
        format!("(A:1[{}],B:2);", "[".repeat(500) + &"]".repeat(500)),
    ];

    for (i, case) in cases.iter().enumerate() {
        let result = std::panic::catch_unwind(|| exercise_pipeline(case));
        assert!(result.is_ok(), "panicked on pathological case {}", i);
    }
}
