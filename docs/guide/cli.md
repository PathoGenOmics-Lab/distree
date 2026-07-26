# CLI reference

```
Usage: distree [OPTIONS] <phylogeny>
```

## Argument

`<phylogeny>`

:   Path to the tree file in Newick format. Use `-` to read from stdin. The file
    must hold exactly one tree, and its tip labels must be unique. See
    [Input](input.md).

## Options

| Flag | Default | Effect |
|:--|:--|:--|
| `--midpoint` | off | Midpoint-root the tree before computing distances. Ignored with `--topology`. See [Midpoint rooting](rooting.md) |
| `--lmm` | off | Write the variance-covariance matrix: each cell is the root-to-MRCA distance. Takes precedence over `--topology` |
| `--topology` | off | Ignore branch lengths and count edges. Values are written as integers |
| `--lower` | off | Write a PHYLIP lower triangle: a taxa count, then one row per taxon with no diagonal |
| `-o, --output FILE` | stdout | Write the matrix to a file. Not created until the tree has parsed |
| `-p, --precision N` | `10` | Decimal places, from 0 to 30. Ignored by `--topology` |
| `-t, --threads N` | all cores | Threads for the parallel row computation. Must be at least 1 |
| `-h, --help` | | Print help |
| `-V, --version` | | Print version |

Without `--lmm` or `--topology`, distree computes patristic distances.

## How the modes interact

| Combination | Result |
|:--|:--|
| `--lmm --topology` | `--lmm` wins, with a warning. They ask for different things |
| `--midpoint --topology` | The rooting is skipped, with a warning. Edge counts do not depend on the root, and the inserted node would add a hop |
| `--midpoint --lmm` | Both apply. This is the combination `--midpoint` exists for |
| `--midpoint` alone | Applies, and changes nothing: a patristic matrix is the same under any rooting |
| `--lower` with any mode | Applies to all three |
| `-p` with `--topology` | Ignored; edge counts are integers |

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
| `Empty input: no Newick tree found.` | The file is empty or only whitespace |
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

# From stdin, capped at 4 threads
gunzip -c tree.nwk.gz | distree - -t 4 -o distances.tsv

# Just the header, to check the tip ordering
distree tree.nwk | head -1 | tr '\t' '\n' | tail -n +2
```
