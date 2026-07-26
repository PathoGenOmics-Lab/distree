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
fn test_binary_midpoint_does_not_inflate_topology() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1,B:1):1,(C:1,D:1):1);").unwrap();

    let (code, plain, _) = run(&["--topology", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    let (code, rooted, stderr) = run(
        &["--topology", "--midpoint", tree.to_str().unwrap()],
        None,
    );
    assert_eq!(code, 0);
    assert_eq!(
        plain, rooted,
        "rooting must not change the number of edges between two leaves"
    );
    assert!(stderr.contains("--midpoint"), "stderr: {}", stderr);
    // A-C crosses the root: 2 edges up + 2 edges down
    assert!(plain.contains("\t4\t"), "unexpected matrix:\n{}", plain);
}

#[test]
fn test_binary_midpoint_preserves_patristic() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "(((A:0.5,B:0.3):0.4,C:0.9):0.1,D:1.2);").unwrap();

    let (code, plain, _) = run(&[tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    let (code, rooted, _) = run(&["--midpoint", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    assert_eq!(
        plain, rooted,
        "patristic distances do not depend on the root"
    );
}

#[test]
fn test_binary_lower_triangle() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1,B:2):3,C:4);").unwrap();

    let (code, stdout, _) = run(&["--lower", tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    let lines: Vec<&str> = stdout.lines().collect();
    // PHYLIP lower triangle: header line (taxa count) + 3 leaf rows
    assert_eq!(lines.len(), 4, "should have 1 header + 3 rows, got: {:?}", lines);
    // First line: taxa count
    assert_eq!(lines[0].trim(), "3", "first line should be taxa count");
    // Row 0 (A): label only, 0 distance columns
    assert_eq!(lines[1].trim(), "A", "first data row is just label A");
    // Row 1 (B): label + 1 distance
    let cols_b: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(cols_b.len(), 2, "row B: label + 1 distance");
    assert_eq!(cols_b[0], "B");
    // Row 2 (C): label + 2 distances
    let cols_c: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(cols_c.len(), 3, "row C: label + 2 distances");
    assert_eq!(cols_c[0], "C");
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
fn test_binary_rejects_newline_in_label() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1.0,'bad\nlabel':2.0):0.5,C:3.0);").unwrap();

    let (code, _stdout, stderr) = run(&[tree.to_str().unwrap()], None);
    assert_ne!(code, 0, "a newline in a label should be rejected");
    assert!(stderr.contains("newline"), "stderr: {}", stderr);
}

#[test]
fn test_binary_warns_about_unlabeled_leaves() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "((A:1.0,:2.0):0.5,C:3.0);").unwrap();

    let (code, stdout, stderr) = run(&[tree.to_str().unwrap()], None);
    assert_eq!(code, 0);
    assert!(stderr.contains("no label"), "stderr: {}", stderr);
    // The unlabeled leaf is absent from the matrix, so say so
    assert_eq!(stdout.lines().count(), 3, "header + A + C");
}

#[test]
fn test_binary_rejects_out_of_range_precision() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "(A:1,B:2);").unwrap();

    let (code, _stdout, stderr) = run(&["-p", "50000000", tree.to_str().unwrap()], None);
    assert_ne!(code, 0, "an out-of-range precision should be rejected");
    assert!(
        stderr.contains("--precision"),
        "stderr should name the flag: {}",
        stderr
    );
    assert!(
        !stderr.contains("panicked"),
        "should fail cleanly, not panic: {}",
        stderr
    );
}

#[test]
fn test_binary_rejects_zero_threads() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("t.nwk");
    std::fs::write(&tree, "(A:1,B:2);").unwrap();

    let (code, _stdout, stderr) = run(&["-t", "0", tree.to_str().unwrap()], None);
    assert_ne!(code, 0, "--threads 0 should be rejected");
    assert!(stderr.contains("--threads"), "stderr: {}", stderr);
}

#[test]
fn test_binary_keeps_output_file_on_parse_error() {
    let dir = tempfile::tempdir().unwrap();
    let tree = dir.path().join("bad.nwk");
    std::fs::write(&tree, "('unclosed:1.0,B:2.0);").unwrap();
    let out = dir.path().join("results.tsv");
    std::fs::write(&out, "previous results\n").unwrap();

    let (code, _stdout, _stderr) = run(
        &[tree.to_str().unwrap(), "-o", out.to_str().unwrap()],
        None,
    );
    assert_ne!(code, 0, "should exit with error on unclosed quote");
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "previous results\n",
        "a failed run must not clobber an existing output file"
    );
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
