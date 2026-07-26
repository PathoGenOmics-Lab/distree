# Input

distree reads one Newick tree, from a file or from stdin:

```bash
distree tree.nwk
gunzip -c tree.nwk.gz | distree -
```

## What is accepted

**Branch lengths** in any form Rust's float parser takes, including scientific
notation and negatives: `A:0.1`, `B:1.5e-3`, `C:2.0E+1`, `D:-0.5`. A missing
length is 0.

**Polytomies.** `(A:1,B:2,C:3);` is a node with three children, and nothing in
distree assumes a bifurcating tree.

**Internal node labels**, quoted or not, with or without a length:
`((A:1,B:2)clade:0.5,C:3);`. Bootstrap values live in this position too and are
read as labels, which is to say ignored.

**Quoted tip labels**, single or double, so a label may hold spaces:
`('Taxon A':1.0,"Taxon B":2.0);`. Per the Newick convention a doubled quote
inside a quoted label is one literal quote, so `'it''s'` is the label `it's`.

**Comments** in square brackets, before or after a label, nested, and at the
start of the file:

```
[&R] ((A:0.1[&&NHX:S=human],B:0.2[&rate=0.5]):0.3,C:0.4);
```

The `[&R]` / `[&U]` rooting marker that IQ-TREE, MrBayes and BEAST write is a
comment like any other. So is NHX and BEAST per-branch metadata, wherever it
sits relative to the label. All of it is skipped; the labels and the branch
lengths survive.

**Non-ASCII labels.** Accents, Greek letters and CJK text pass through
unchanged, so a label keeps matching the sample name it came from.

**Whitespace and newlines** anywhere outside a quoted label, which is how most
tools wrap a long tree over several lines.

**A missing trailing semicolon**, which some tools omit.

## What is rejected

The guiding rule is that distree would rather stop than hand back a matrix that
looks right and is not. Each of these fails with a message naming the offending
position, and nothing is written to `-o`.

| Input | Message |
|:--|:--|
| Empty file | `Empty input: no Newick tree found.` |
| `((A:1,B:2),C:3` | `Unexpected end of input: 1 unclosed '(' remain. The tree is truncated.` |
| `(A:1,B:2);(C:1,D:2);` | `the file appears to hold more than one tree` |
| `(A:1,B:2);garbage` | `Unexpected content after the end of the tree` |
| `('unclosed:1,B:2);` | `Unclosed quote starting at position 1.` |
| `(A:1,B:2)[oops;` | `Unclosed comment starting at position 9: no matching ']'.` |
| `(A:,B:2);` | `Expected a numeric branch length` |
| `(A:1,A:2);` | `Duplicate leaf name 'A' found. Leaf names must be unique.` |
| A label holding a tab, newline or carriage return | `contains a tab character, which would corrupt the output by splitting the row` |
| A tree with no labelled tips | `No labeled leaves found in the tree.` |

!!! warning "One tree per file"

    A file of bootstrap replicates or a posterior sample holds many trees.
    distree computes a matrix for one tree, and earlier versions silently used
    the first one. It now refuses, because a matrix labelled with your dataset's
    name that describes replicate 1 of 1,000 is worse than an error. Split the
    file first:

    ```bash
    split -l 1 posterior.trees tree_
    for t in tree_*; do distree "$t" -o "${t}.tsv"; done
    ```

!!! question "Truncated trees"

    A tree cut short by an interrupted download or a full disk used to parse:
    the open parentheses were closed at end of input, and every internal node
    left dangling silently lost its branch length. The matrix came out looking
    entirely reasonable. That is now an error naming how many parentheses were
    still open.

## Warnings

These do not stop the run. They go to stderr, so they stay out of a piped
matrix.

**Negative branch lengths.**

```
Warning: negative branch lengths detected in the tree. Some distances may come
out negative, and --midpoint cannot locate the diameter reliably.
```

Neighbour-joining and BioNJ produce these routinely. distree reports the
distances as computed rather than clamping them to zero, since a clamp would
claim two distinct taxa sit on top of each other. Midpoint rooting, on the other
hand, finds the longest path by a two-pass sweep that assumes non-negative
lengths, so `--midpoint` on such a tree is not reliable.

**No branch lengths at all**, in patristic mode:

```
Warning: no branch lengths detected, all patristic distances will be zero.
Consider --topology.
```

A cladogram has topology and nothing else. `--topology` is the mode that reads
it.

**Unlabelled tips.**

```
Warning: 1 leaf/leaves have no label and were excluded from the matrix.
```

A tip with no label cannot be named in a row or a column, so it is left out. The
matrix is correct for the tips that remain, but it is smaller than the tree, and
the warning is there so that is not a surprise.

**Conflicting modes.** `--lmm` and `--topology` ask for different things;
passing both uses `--lmm` and says so. `--midpoint` with `--topology` is ignored,
for the reason in [Midpoint rooting](rooting.md).

## Label rules

Tip labels become row and column headers in a TSV, which puts two requirements
on them.

**They must be unique.** Two tips called `A` would give two identical rows with
no way to tell which is which, so a duplicate is an error rather than a warning.

**They must not hold a tab, a newline or a carriage return.** Any of the three
splits a row and leaves a file whose row count no longer matches its header.
A quoted label can hold them, so this is checked rather than assumed. Replace
them with underscores before running distree.

Spaces are fine, as long as the label is quoted in the Newick. So is anything
else UTF-8 can express.

## Next

- [Distance modes](distances.md), for what distree does with the tree.
- [Output](output.md), for the shape of what comes back.
