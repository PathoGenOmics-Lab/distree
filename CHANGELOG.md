# Changelog

## [1.0.1] - 2026-04-04

### Fixed
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
- Comprehensive test suite (24 tests)

### Changed
- Codebase split into modules: `parser.rs`, `tree.rs`, `lca.rs`, `midpoint.rs`
- Version string now derived from Cargo.toml via `env!("CARGO_PKG_VERSION")`
- Removed unused `--format` flag
- Fixed clippy warnings (`&Vec<Node>` → `&[Node]`)
- Fixed contributors link in README

## [1.0.0] - 2025-06-01

### Added
- Initial release
- Patristic distance matrix extraction
- Topological distance computation
- LMM (var-covar) matrix output
- Midpoint rooting option
- Parallel computation with rayon
