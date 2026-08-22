<!--
Fill in what applies and delete what does not. A short pull request needs a
short description; nobody is asking for an essay to fix a typo.
-->

## What this changes

<!-- One or two sentences. What is different after this is merged? -->

## Why

<!-- The problem it solves. Link an issue with "Closes #123" if there is one. -->

## How it was verified

<!--
The important part, and the one a reviewer cannot reconstruct.

Not "it should work", but what you ran and what it printed. For example: the
test you added and the fact that it fails without the fix, the tree you ran
through it and the matrix you got before and after, or the ape cross-validation
if you touched anything the distances depend on.

If it is a change CI already covers, say which check covers it.
-->

## Checklist

- [ ] `cargo test` passes, in debug and release.
- [ ] `cargo clippy --all-targets -- -D warnings` is clean.
- [ ] New behaviour has a test, or I have said above why it does not.
- [ ] Documentation under `docs/` is updated if the change is user-visible.
- [ ] If the output changed, I have said so: byte-identical output is the
      default expectation, and a change to it is a change to everyone's files.
