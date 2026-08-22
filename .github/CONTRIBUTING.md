# Contributing

Issues and pull requests are welcome at
[PathoGenOmics-Lab/distree](https://github.com/PathoGenOmics-Lab/distree).

The fullest version of this is on the documentation site, under
[Contributing](https://pathogenomics-lab.github.io/distree/about/contributing/).
What follows is the short form.

## Reporting a problem

The most useful bug report is one somebody else can run: the exact command, the
version, a tree that still shows it, and what you expected against what you got.

The dangerous failures here are the quiet ones. A matrix that is wrong looks
exactly like a matrix that is right, so a report saying "these distances
disagree with `ape::cophenetic` and here are both" is worth more than a stack
trace.

If the tree carries anything identifiable, do not attach it. Tip labels can be
replaced with `S1`, `S2` and so on without changing the behaviour.

## Building and checking

```bash
cargo build --release
cargo test
cargo test --release
cargo clippy --all-targets -- -D warnings
```

CI runs all four on every push and every pull request. A change that fails any
of them will not go in.

Two deeper checks are run on purpose rather than on every push:

```bash
# Against ape, over 250 generated trees. Needs R with ape and phangorn.
Rscript scripts/crossvalidate.R 250 target/release/distree

# Coverage-guided fuzzing. Needs a nightly toolchain.
cargo +nightly fuzz run parse_newick fuzz/seeds
```

## Conventions the code follows

- **Refuse rather than guess.** If the input is ambiguous or truncated, error
  with the position. A matrix that looks right and is not costs more than a
  failed run.
- **Warnings go to stderr**, never into the matrix, so a piped run stays clean.
- **Keep the tree traversals iterative.** The parser, the flattener, the LCA
  build and `RawNode`'s `Drop` are all loops with explicit stacks, because a
  ladder-shaped tree hundreds of thousands of levels deep is a real input.
- **Do not put work in the serial write loop.** Whatever a worker can do to its
  own row, it should, including formatting the numbers.
- **Output is byte-identical unless the change is about the output.** If it
  changes, say so in the pull request: it is a change to everyone's files.

## Code of conduct

By taking part you agree to the [Code of Conduct](CODE_OF_CONDUCT.md).
