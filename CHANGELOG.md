# Changelog

## [1.0.1] - 2026-07-29

### Fixed
- Branch lengths that are not finite are rejected. `A:1e910` parses as infinity rather than failing, and an infinite depth made every distance `inf` with `NaN` down the diagonal, reported with a success exit code. The run also stops when finite depths are large enough that `d_i + d_j` overflows
- A leaf label containing a space is rejected in `--lower` mode: PHYLIP readers treat whitespace as the end of the name, so every row after it landed a column out
- `--lower` with `--lmm` now warns that the omitted diagonal holds each leaf's root-to-tip length rather than zeros
- Parallel computation now speeds the run up instead of slowing it down. One parallel job per row could not pay for synchronising the thread pool, and float formatting sat in the serial write loop; rows are now batched and each worker formats its own. A 20,000-tip matrix went from 10.6 s to 1.6 s, and an 8,000-tip one now scales from 1.26 s on one core to 0.21 s on eight
- Trees prefixed with a `[&R]` / `[&U]` rooting marker, or carrying a comment before a label, no longer fail to parse
- Non-ASCII leaf labels (accents, Greek, CJK) are preserved instead of being mangled into mojibake
- Truncated trees, trailing content, unclosed comments and free-form text are rejected instead of yielding a plausible but wrong matrix
- Files holding more than one tree are rejected instead of silently using the first
- Apostrophes in unquoted labels (`O'Brien`) no longer break parsing, and a doubled quote keeps the whitespace around it
- Deeply nested trees no longer overflow the stack when the parse tree is freed
- `-o FILE` no longer truncates an existing file before the tree has parsed
- Write failures are reported instead of being discarded along with a success exit code; a closed pipe exits quietly
- `--precision` out of range is rejected instead of panicking inside the formatter
- `--threads 0` is rejected instead of quietly starting one thread per core
- Newlines and carriage returns in labels are rejected like tabs, and unlabeled leaves are reported rather than dropped in silence
- `--midpoint` no longer adds a hop to every topological distance crossing the midpoint edge
- Negative patristic distances are reported as computed instead of being rounded up to zero, which claimed distinct taxa were identical
- Midpoint rooting (`--midpoint`) rewritten to fix infinite loop caused by cyclic graph construction
- NHX/bracket comment annotations (e.g., `[&&NHX:S=human]`) no longer crash the parser
- Single-quoted labels (e.g., `'Taxon A'`) now parsed correctly
- Double-quoted labels (e.g., `"Taxon A"`) now parsed correctly
- Whitespace and newlines in Newick strings no longer crash the parser
- Stack overflow on deeply nested trees (>5,000 levels) — parser and flattener rewritten iteratively
- Duplicate leaf names are now detected with a clear error message
- `--lmm --topology` conflict now warns instead of silently choosing LMM
- Topology mode now outputs integers instead of float decimals
- Floating point output uses configurable precision instead of raw representation
- Header row always prints leading tab for R/Python distance matrix compatibility

### Added
- `--precision` / `-p` flag to control decimal places in output (default: 10)
- `-t` / `--threads` flag to control number of parallel threads
- Stdin support: use `-` as the phylogeny argument to read from stdin
- Warning when no branch lengths are detected in patristic mode
- Warning when negative branch lengths are found in the tree
- `--npy` writes the matrix as a NumPy array of 64-bit floats instead of text: exact, half the size, and 2.6x faster. With `--lower` it writes the condensed vector SciPy reads. Labels go to `<FILE>.labels.txt`
- `--taxa FILE` restricts the matrix to a list of leaf labels, keeping the distances the full tree gives rather than pruning it
- `--stats` prints a summary to stderr: leaves, nodes, mode, cells, and the minimum, maximum and mean off the diagonal
- Gzipped input is read directly, detected by content rather than by file name, so it works from stdin too
- CITATION.cff with DOI
- MkDocs Material documentation site with a tutorial, published to GitHub Pages
- Comprehensive test suite (96 tests), including randomised checks of midpoint rooting and of MRCA queries against a brute-force walk

### Changed
- Codebase split into modules: `parser.rs`, `tree.rs`, `lca.rs`, `midpoint.rs`
- Fixed-precision float formatting is done by hand where the rounding is unambiguous and handed to the standard formatter where it is not, which cut a text run by a third to a half with byte-identical output
- Cross-validated against R's ape (`cophenetic.phylo`, `vcv.phylo`) and `phangorn::midpoint` over 250 generated trees; reference matrices for six of them are committed so `cargo test` checks the same thing without R
- The parser and the whole pipeline are fuzzed, in the test suite on every push and with `cargo fuzz` for longer runs
- The tree logic moved into a library target, so it can be fuzzed, documented and used from other Rust code
- LCA binary-lifting table stores plain `usize` rather than `Option<usize>`, halving the memory it needs
- Output buffer raised from 8 KB to 1 MB, so a multi-gigabyte matrix is not hundreds of thousands of write syscalls
- The input text and the recursive parse tree are released once the flat node array is built, rather than held to the end of the run
- Version string now derived from Cargo.toml via `env!("CARGO_PKG_VERSION")`
- Removed unused `--format` flag
- Fixed clippy warnings (`&Vec<Node>` → `&[Node]`)
- Fixed contributors link in README
- README output examples replaced with the real output of a worked example; the previous ones lost their header row to the code fence and the topological matrix was not a realisable tree
- CI runs on every branch, lints the test code, and runs the test suite in release as well as debug

## [1.0.0] - 2025-06-01

### Added
- Initial release
- Patristic distance matrix extraction
- Topological distance computation
- LMM (var-covar) matrix output
- Midpoint rooting option
- Parallel computation with rayon
