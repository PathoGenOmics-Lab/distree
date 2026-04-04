/// Integration tests: run the `distree` binary and check stdout/exit codes.
use std::io::Write;
use std::process::{Command, Stdio};

/// Build the binary once and return its path.
fn bin() -> std::path::PathBuf {
    // `cargo test` runs with the project root as cwd; the binary ends up here.
    let mut p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("target/debug/distree");
    p
}

fn run(args: &[&str], stdin: Option<&str>) -> (i32, String, String) {
    let mut cmd = Command::new(bin());
    cmd.args(args);
    if stdin.is_some() {
        cmd.stdin(Stdio::piped());
    }
    cmd.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = cmd.spawn().expect("failed to spawn distree binary");
    if let Some(input) = stdin {
        child
            .stdin
            .take()
            .unwrap()
            .write_all(input.as_bytes())
            .unwrap();
    }
    let out = child.wait_with_output().unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn test_binary_basic_patristic() {
    // Write a temp tree file
    let dir = tempfile::tempdir().expect("tempdir");
    let tree = dir.path().join("tree.nwk");
    std::fs::write(&tree, "((A:1.0,B:2.0):0.5,C:3.0);").unwrap();

    let (code, stdout, _stderr) = run(&[tree.to_str().unwrap()], None);
    assert_eq!(code, 0, "exit code should be 0");

    // Header row: first field is empty, then A B C (sorted)
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 4, "header + 3 rows");
    assert!(lines[0].starts_with('\t'), "header starts with tab");
    assert!(lines[1].starts_with("A\t"), "first data row starts with A");
}

#[test]
fn test_binary_stdin() {
    let newick = "(A:1.0,B:3.0);";
    let (code, stdout, _) = run(&["-"], Some(newick));
    assert_eq!(code, 0);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3); // header + 2 rows
}

#[test]
fn test_binary_topology_flag() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1,B:2):3,C:4);").unwrap();

    let (code, stdout, _) = run(&["--topology", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    // Topology values should be integers (no decimal point)
    for line in stdout.lines().skip(1) {
        for cell in line.split('\t').skip(1) {
            assert!(
                !cell.contains('.'),
                "Topology output should be integers, got: {}",
                cell
            );
        }
    }
}

#[test]
fn test_binary_lower_triangle() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1,B:2):3,C:4);").unwrap();

    let (code, stdout, _) = run(&["--lower", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    let lines: Vec<&str> = stdout.lines().collect();
    // Lower triangle: 3 leaves → 3 rows, no header
    assert_eq!(lines.len(), 3, "should have 3 rows");
    // Row 0 is empty (no columns for first leaf)
    assert!(lines[0].is_empty(), "first row should be empty, got: {:?}", lines[0]);
    // Row 1 has exactly 1 value
    assert_eq!(lines[1].split('\t').count(), 1, "row 1 should have 1 column");
    // Row 2 has exactly 2 values separated by a tab
    assert_eq!(lines[2].split('\t').count(), 2, "row 2 should have 2 columns");
}

#[test]
fn test_binary_duplicate_leaf_error() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "(A:1.0,A:2.0);").unwrap();

    let (code, _stdout, stderr) = run(&[tree.to_str().unwrap()], None);
    assert_ne!(code, 0, "should exit with error on duplicate leaves");
    assert!(stderr.contains("Duplicate"), "stderr should mention 'Duplicate': {}", stderr);
}

#[test]
fn test_binary_empty_input_error() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("empty.nwk");
    std::fs::write(&tree, "").unwrap();

    let (code, _stdout, stderr) = run(&[tree.to_str().unwrap()], None);
    assert_ne!(code, 0);
    assert!(stderr.contains("Empty") || stderr.contains("empty"), "stderr: {}", stderr);
}

#[test]
fn test_binary_precision_flag() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "(A:1.0,B:3.0);").unwrap();

    let (code, stdout, _) = run(&["-p", "3", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    // With precision=3, values like 4.000 should appear (exactly 3 decimal places)
    for line in stdout.lines().skip(1) {
        for cell in line.split('\t').skip(1) {
            let dot_pos = cell.find('.');
            if let Some(pos) = dot_pos {
                let decimals = cell.len() - pos - 1;
                assert_eq!(decimals, 3, "Expected 3 decimal places, got '{}' in line '{}'", cell, line);
            }
        }
    }
}
