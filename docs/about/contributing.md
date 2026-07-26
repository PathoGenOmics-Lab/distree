# Contributing

Issues and pull requests are welcome at
[PathoGenOmics-Lab/distree](https://github.com/PathoGenOmics-Lab/distree).

## Reporting a problem

The most useful bug report is one somebody else can run. Include:

- the exact command, with every flag;
- the version (`distree --version`) and how it was installed;
- the tree, or a smaller tree that still shows the problem;
- what you expected and what you got.

The dangerous failures here are the quiet ones. A matrix that is wrong looks
exactly like a matrix that is right, so a report saying "these distances
disagree with `ape::cophenetic` and here are both" is worth more than a stack
trace.

If the tree carries anything identifiable, do not attach it. Tip labels can be
replaced with `S1`, `S2` and so on without changing the behaviour, and a
generated tree of the same shape usually reproduces the problem.

## Building and checking

```bash
cargo build --release
cargo test
cargo test --release
cargo clippy --all-targets -- -D warnings
```

CI runs all four on every push and every pull request. A change that fails any
of them will not go in.

## Tests

The suite is in three places:

- **Unit tests** next to the code they cover, in `src/parser.rs`, `src/lca.rs`,
  `src/midpoint.rs` and `src/main.rs`.
- **Integration tests** in `tests/integration.rs`, which run the built binary
  and check its stdout, its stderr and its exit code.
- **Randomised property checks**, which are where the confidence comes from.
  `src/testutil.rs` generates trees with polytomies and zero-length branches
  from a seeded xorshift, and the checks assert properties rather than
  particular numbers: midpoint rooting preserves every pairwise distance and
  leaves a valid tree with the deepest tip half a diameter from the root; the
  binary-lifting MRCA agrees with walking up from both nodes, over every pair of
  nodes; and a patristic distance is never negative on a tree with non-negative
  branch lengths.
- **Cross-validation against ape** in `tests/crossvalidation.rs`, checked
  against reference matrices committed under `tests/fixtures/`.
- **Mutation testing** in `tests/torture.rs`, which throws truncated, flipped
  and duplicated Newick at the parser and the pipeline.

### Cross-validating against ape

Everything above is written against the same understanding of the problem as
the code, so none of it would catch a distance being *defined* wrongly. ape is
the independent check. With R, `ape` and `phangorn` installed:

```bash
cargo build --release
Rscript scripts/crossvalidate.R 250 target/release/distree
```

It compares all four modes against `cophenetic.phylo`, `vcv.phylo`,
`cophenetic.phylo` over unit branch lengths, and `phangorn::midpoint`. The
expected result is agreement to under 1e-9, which is the 12-decimal text
round-trip, and exactly zero for edge counts.

To refresh the committed fixtures after a deliberate change:

```bash
Rscript scripts/crossvalidate.R 6 target/release/distree --write-fixtures
```

The same comparison runs as a weekly GitHub Actions job, and on demand from the
Actions tab.

### Fuzzing

`tests/torture.rs` runs on every push and needs nothing extra. The
coverage-guided version needs a nightly toolchain:

```bash
cargo +nightly install cargo-fuzz
cargo +nightly fuzz run parse_newick fuzz/seeds
cargo +nightly fuzz run full_pipeline fuzz/seeds
```

`parse_newick` throws arbitrary text at the parser; `full_pipeline` takes
everything that parses and pushes it through flattening, rooting and every
distance mode. Both assert only that nothing panics and nothing hangs, which is
a low bar and has already caught a real bug: a branch length of `1e910` parsed
as infinity and turned the matrix into `inf` with `NaN` down the diagonal.

A parser change wants a case in `src/parser.rs` for what should now parse **and**
one for what should still be rejected. The failures worth guarding are the ones
where a malformed tree produced a plausible matrix instead of an error.

## Writing a fix

A few conventions the codebase follows:

- **Refuse rather than guess.** If the input is ambiguous or truncated, error
  with the position. A matrix that looks right and is not costs more than a
  failed run.
- **Warnings go to stderr**, never into the matrix, so a piped run stays clean.
- **Keep the tree traversals iterative.** The parser, the flattener, the LCA
  build and even `RawNode`'s `Drop` are loops with explicit stacks, because a
  ladder-shaped tree hundreds of thousands of levels deep is a real input.
- **Do not put work in the serial write loop.** Whatever a worker can do to its
  own row, it should, including formatting the numbers.

## Documentation

This site is MkDocs Material, and the pages live in `docs/`:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

```bash
mkdocs build --strict
```

`--strict` turns a broken internal link or a dead anchor into a failed build,
which is what CI runs. Numbers quoted in the docs, in the benchmark tables above
all, should be measured rather than estimated, and the page should say which
when it is not obvious.
