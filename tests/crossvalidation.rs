//! Check distree's output against reference matrices computed by ape.
//!
//! The rest of the suite is self-consistent: it asserts that midpoint rooting
//! preserves distances, that the binary-lifting MRCA agrees with walking up
//! from both nodes, that a patristic distance is never negative. All of that is
//! written against the same understanding of the problem as the code, so none
//! of it would catch a systematic error in what a distance is defined to be.
//!
//! The fixtures here come from somewhere else. `tests/fixtures/*.nwk` are
//! random trees, and the matrices beside them were computed by R's ape package
//! (`cophenetic.phylo` for patristic, `vcv.phylo` for the variance-covariance
//! matrix, and `cophenetic.phylo` over unit branch lengths for edge counts).
//!
//! Regenerate them, and cross-validate over hundreds more trees than are
//! committed here, with:
//!
//! ```text
//! Rscript scripts/crossvalidate.R 250 target/release/distree
//! Rscript scripts/crossvalidate.R 6 target/release/distree --write-fixtures
//! ```

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_distree")
}

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// Parse a square TSV into labels and a row-major matrix.
fn parse_matrix(text: &str) -> (Vec<String>, HashMap<(String, String), f64>) {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .expect("matrix has a header")
        .split('\t')
        .skip(1) // leading empty cell
        .map(str::to_string)
        .collect();

    let mut cells = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let row = fields.next().expect("row label").to_string();
        for (col, value) in header.iter().zip(fields) {
            let v: f64 = value
                .parse()
                .unwrap_or_else(|_| panic!("cell '{}' is not a number", value));
            cells.insert((row.clone(), col.clone()), v);
        }
    }
    (header, cells)
}

fn run_distree(tree: &PathBuf, flags: &[&str]) -> String {
    let out = Command::new(bin())
        .arg(tree)
        .args(flags)
        .args(["-p", "12"])
        .output()
        .expect("failed to run distree");
    assert!(
        out.status.success(),
        "distree failed on {:?} with {:?}: {}",
        tree,
        flags,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("distree writes UTF-8")
}

/// Compare distree's matrix for `flags` against the fixture for `mode`.
fn check_mode(stem: &str, mode: &str, flags: &[&str], tolerance: f64) {
    let dir = fixture_dir();
    let tree = dir.join(format!("{}.nwk", stem));
    let expected_path = dir.join(format!("{}.{}.tsv", stem, mode));

    let expected_text = std::fs::read_to_string(&expected_path)
        .unwrap_or_else(|e| panic!("cannot read {:?}: {}", expected_path, e));

    let (want_labels, want) = parse_matrix(&expected_text);
    let (got_labels, got) = parse_matrix(&run_distree(&tree, flags));

    assert_eq!(
        got_labels, want_labels,
        "{} {}: label ordering differs from ape's",
        stem, mode
    );

    let mut worst = 0.0_f64;
    let mut worst_at = (String::new(), String::new());
    for a in &want_labels {
        for b in &want_labels {
            let key = (a.clone(), b.clone());
            let g = got[&key];
            let w = want[&key];
            let diff = (g - w).abs();
            if diff > worst {
                worst = diff;
                worst_at = key;
            }
        }
    }

    assert!(
        worst <= tolerance,
        "{} {}: distree and ape disagree by {:.3e} at ({}, {}); tolerance is {:.0e}",
        stem,
        mode,
        worst,
        worst_at.0,
        worst_at.1,
        tolerance
    );
}

fn stems() -> Vec<String> {
    let mut found: Vec<String> = std::fs::read_dir(fixture_dir())
        .expect("tests/fixtures exists")
        .filter_map(|e| {
            let name = e.ok()?.file_name().into_string().ok()?;
            name.strip_suffix(".nwk").map(str::to_string)
        })
        .collect();
    found.sort();
    assert!(!found.is_empty(), "no fixture trees found");
    found
}

/// Both tools write 12 decimals, so a value in the hundreds carries about
/// 1e-10 of round-trip error. Anything larger is a real disagreement.
const TOLERANCE: f64 = 1e-9;

#[test]
fn test_patristic_matches_ape_cophenetic() {
    for stem in stems() {
        check_mode(&stem, "patristic", &[], TOLERANCE);
    }
}

#[test]
fn test_lmm_matches_ape_vcv() {
    for stem in stems() {
        check_mode(&stem, "lmm", &["--lmm"], TOLERANCE);
    }
}

#[test]
fn test_topology_matches_ape_unit_branch_cophenetic() {
    for stem in stems() {
        // Edge counts are integers on both sides, so they must agree exactly
        check_mode(&stem, "topology", &["--topology"], 0.0);
    }
}

#[test]
fn test_midpoint_leaves_patristic_alone() {
    // ape has no opinion to check here, but the fixtures are a broader set of
    // shapes than the hand-written trees, and rooting must not move a
    // patristic distance on any of them.
    for stem in stems() {
        let tree = fixture_dir().join(format!("{}.nwk", stem));
        assert_eq!(
            run_distree(&tree, &[]),
            run_distree(&tree, &["--midpoint"]),
            "{}: --midpoint changed the patristic matrix",
            stem
        );
    }
}
