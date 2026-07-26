# CLI reference

```
Usage: distree [OPTIONS] <phylogeny>
```

## Argument

`<phylogeny>`

:   Path to the tree file in Newick format, plain or gzipped. Use `-` to read
    from stdin. The file must hold exactly one tree, and its tip labels must be
    unique. See [Input](input.md).

## Options

| Flag | Default | Effect |
|:--|:--|:--|
| `--midpoint` | off | Midpoint-root the tree before computing distances. Ignored with `--topology`. See [Midpoint rooting](rooting.md) |
| `--lmm` | off | Write the variance-covariance matrix: each cell is the root-to-MRCA distance. Takes precedence over `--topology` |
| `--topology` | off | Ignore branch lengths and count edges. Values are written as integers |
| `--lower` | off | Write a PHYLIP lower triangle: a taxa count, then one row per taxon with no diagonal. With `--npy`, the condensed vector instead |
| `--npy` | off | Write a NumPy `.npy` array of 64-bit floats instead of text. Requires `-o`; labels go to `<FILE>.labels.txt` |
| `--taxa FILE` | all tips | Restrict the matrix to the labels listed in `FILE`, one per line |
| `--stats` | off | Print a summary of the run to stderr |
| `-o, --output FILE` | stdout | Write the matrix to a file. Not created until the tree has parsed |
| `-p, --precision N` | `10` | Decimal places, from 0 to 30. Ignored by `--topology` and by `--npy` |
| `-t, --threads N` | all cores | Threads for the parallel row computation. Must be at least 1 |
| `-h, --help` | | Print help |
| `-V, --version` | | Print version |

Without `--lmm` or `--topology`, distree computes patristic distances.

## Selecting a subset

```bash
distree tree.nwk --taxa keep.txt -o subset.tsv
```

`keep.txt` holds one leaf label per line. Blank lines and lines starting with
`#` are skipped, and duplicates are collapsed. Order does not matter: the matrix
comes out alphabetical as always.

A label that is not a leaf of the tree is an error naming it, rather than a
matrix quietly smaller than the list.

The distances themselves are the ones from the full tree. A subset filters the
output, it does not prune the tree, so the path between two tips is the same
whether or not the tips in between are in the matrix. That is the opposite of
what pruning the tree first would give you.

## The run summary

```bash
distree tree.nwk --stats -o distances.tsv
```

```
--- Statistics ---
Leaves in matrix:  7
Nodes in tree:     13
Mode:              patristic distance
Cells written:     49
Minimum:           0.0000004534
Maximum:           0.0000612030
Mean:              0.0000358688
Time:              0.001s
```

It goes to stderr, so a piped run is unaffected. The minimum, maximum and mean
are taken off the diagonal, which is zero by construction in the distance modes
and would only pull the summary towards it. When `--taxa` is in use, the summary
reports both counts.

The numbers are accumulated by the workers as they compute, so the summary costs
one pass and nothing measurable.

## How the modes interact

| Combination | Result |
|:--|:--|
| `--lmm --topology` | `--lmm` wins, with a warning. They ask for different things |
| `--midpoint --topology` | The rooting is skipped, with a warning. Edge counts do not depend on the root, and the inserted node would add a hop |
| `--midpoint --lmm` | Both apply. This is the combination `--midpoint` exists for |
| `--midpoint` alone | Applies, and changes nothing: a patristic matrix is the same under any rooting |
| `--lower` with any mode | Applies to all three |
| `-p` with `--topology` | Ignored; edge counts are integers |
| `--npy` without `-o` | An error. Binary output needs a file |
| `--npy` with `--lower` | Writes the condensed vector SciPy reads, not PHYLIP's triangle |
| `--npy` with `-p` | The precision is ignored; a 64-bit float is exact |
| `--npy` with `--topology` | Works, with a warning that integers are stored as floats |
| `--taxa` with any mode | Applies to all three, and to both output formats |

## Exit codes

| Code | Meaning |
|:--|:--|
| 0 | Success, including a downstream pipe closing early |
| 1 | Failure. The message is on stderr and begins with `Error:` |

## Messages

Everything below goes to stderr, never into the matrix.

### Errors

| Message | Cause |
|:--|:--|
| `Cannot open '<path>': ...` | The tree file is missing or unreadable |
| `... is gzip-compressed but could not be read: ...` | The gzip stream is truncated or corrupt |
| `... is not valid UTF-8 text` | Not a text file, and not gzip either |
| `Empty input: no Newick tree found.` | The file is empty or only whitespace |
| `Branch length '<x>' at position N is not a finite number.` | An exponent past the range of a 64-bit float, which would make every distance infinite |
| `Branch lengths are too large to compute distances with` | The depths are finite but summing two of them overflows |
| `<N> of the <M> labels in '<path>' are not leaves of the tree` | `--taxa` names something the tree does not have |
| `Cannot read the taxa list '<path>': ...` | `--taxa` file missing or unreadable |
| `--npy writes binary, so it needs a file` | `--npy` without `-o` |
| `Failed to parse Newick tree: ...` | Malformed tree. The rest of the message names the position |
| `No labeled leaves found in the tree.` | The tree parsed but no tip carries a label |
| `Duplicate leaf name '<name>' found.` | Two tips share a label; the matrix could not be indexed |
| `Leaf name '<name>' contains a tab character / a newline / a carriage return` | The label would split a row |
| `--precision must be between 0 and 30, got N.` | `-p` out of range |
| `--threads must be at least 1.` | `-t 0` |
| `Cannot write to '<path>': ...` | The output path is not writable |
| `Failed to write '<path>': ...` | The write failed partway, typically a full disk |
| `Failed to initialize thread pool: ...` | Rayon could not start the requested threads |

### Warnings

| Message | Meaning |
|:--|:--|
| `negative branch lengths detected` | Some distances may come out negative, and `--midpoint` cannot locate the diameter reliably |
| `no branch lengths detected, all patristic distances will be zero` | A cladogram. Use `--topology` |
| `N leaf/leaves have no label and were excluded from the matrix` | Unlabelled tips cannot be named in a row, so they are left out |
| `--lmm and --topology are mutually exclusive. Using --lmm.` | Both were passed |
| `--midpoint is ignored in --topology mode.` | Both were passed |
| `--lower omits the diagonal, which in --lmm mode holds each leaf's root-to-tip length` | The lower triangle drops data that is not zeros under `--lmm` |
| `--npy writes 64-bit floats, so --topology edge counts are stored as floats` | Both were passed |

## Examples

```bash
# Patristic distances to a file, six decimals
distree tree.nwk -p 6 -o distances.tsv

# PHYLIP lower triangle for neighbor
distree tree.nwk --lower -p 6 -o infile

# Variance-covariance matrix for PGLS, midpoint-rooted
distree tree.nwk --midpoint --lmm -p 8 -o varcovar.tsv

# Edge counts from a cladogram
distree cladogram.nwk --topology -o topo.tsv

# A gzipped tree, capped at 4 threads
distree tree.nwk.gz -t 4 -o distances.tsv

# Just the 200 samples in a list, as a NumPy array
distree tree.nwk --taxa cohort.txt --npy -o cohort.npy

# The condensed vector scipy.cluster.hierarchy takes directly
distree tree.nwk --lower --npy -o condensed.npy

# Just the header, to check the tip ordering
distree tree.nwk | head -1 | tr '\t' '\n' | tail -n +2
```
