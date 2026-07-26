# Changelog

## [1.0.1] - 2026-04-04

### Fixed
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
- CITATION.cff with DOI
- Comprehensive test suite (73 tests), including randomised checks of midpoint rooting and of MRCA queries against a brute-force walk

### Changed
- Codebase split into modules: `parser.rs`, `tree.rs`, `lca.rs`, `midpoint.rs`
- LCA binary-lifting table stores plain `usize` rather than `Option<usize>`, halving the memory it needs
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
