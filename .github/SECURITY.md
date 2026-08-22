# Security policy

## Supported versions

The latest release. distree is a single binary with no network access and no
server component, so "supported" here means fixes go into the next release
rather than into a maintenance branch.

## Reporting a vulnerability

Report it privately, through
[GitHub's advisory form](https://github.com/PathoGenOmics-Lab/distree/security/advisories/new),
rather than in a public issue.

What is in scope is narrower than it looks for a tool of this shape, and worth
naming. distree reads a Newick file and writes a matrix. It opens no sockets,
runs no subprocesses, and executes nothing from its input. The realistic
concerns are:

- A crafted tree that makes it read or write out of bounds. The parser is
  hand-written and scans bytes, so this is the one worth looking hardest at. It
  is fuzzed for exactly this, in `tests/torture.rs` on every push and with
  `cargo fuzz` for longer runs.
- A crafted tree that makes it allocate without bound, or loop forever.
- A path handling bug in `-o` or `--taxa`.

A tree that produces a *wrong matrix* is a correctness bug rather than a
security one, and belongs in a normal issue where it can be discussed openly.
It is still the most valuable kind of report this project gets.

## What to expect

An acknowledgement within a week. If it is confirmed, a fix and a release, with
credit in the advisory unless you would rather not have it.
