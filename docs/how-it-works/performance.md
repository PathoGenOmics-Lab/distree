# Performance

## Complexity

For a tree with `N` labelled tips and `M` nodes:

| Stage | Time | Memory |
|:--|:--|:--|
| Parse and flatten | `O(M)` | `O(M)` |
| Midpoint rooting | `O(M)` | `O(M)` |
| LCA structure | `O(M log M)` | `O(M log M)` |
| The matrix | `O(N² log M)` | `O(1)` beyond the buffer |

Everything before the matrix is linear or near-linear and, in practice, free.
On a 100,000-tip tree, parsing and building the LCA structure takes 0.05 s. The
run is the matrix, and the matrix is quadratic in the number of tips no matter
what: `N²` cells have to be computed and written.

The `log M` in the cell cost is the MRCA query. It is a handful of array reads
against a table that fits in cache for anything under a few thousand nodes, so
the constant is small; the format of the number costs more than finding it.

## Measured

Balanced binary trees, `-p 6 --lower`, output discarded, release build on an
Apple M4 Pro with 14 cores:

| Tips | Cells | Time | Peak memory |
|--:|--:|--:|--:|
| 1,000 | 0.5 M | 0.00 s | 11 MB |
| 2,000 | 2.0 M | 0.01 s | 24 MB |
| 4,000 | 8.0 M | 0.06 s | 35 MB |
| 8,000 | 32 M | 0.23 s | 32 MB |
| 20,000 | 200 M | 1.8 s | 35 MB |

Time is quadratic in the tips, as it must be. Memory is not: it flattens out
around 35 MB, because past a few thousand tips it is the fixed batch buffer plus
the LCA table rather than anything that grows with the matrix.

## Threads

`-t` caps the thread count; the default is every core. On the 8,000-tip tree:

| `-t` | Time | Speedup |
|--:|--:|--:|
| 1 | 1.52 s | 1.0x |
| 2 | 0.66 s | 2.3x |
| 4 | 0.36 s | 4.2x |
| 8 | 0.24 s | 6.3x |
| 14 | 0.23 s | 6.6x |

Each worker computes and formats whole rows, so the parallel section covers both
halves of the cost and the curve holds up until it runs into memory bandwidth
and the single writer.

!!! note "This was not always true"

    Before version 1.0.1, each row was its own parallel job. A cell is one MRCA
    query and three array reads, so a row of a small matrix could not pay for
    synchronising the thread pool, and the same 8,000-tip run took 2.52 s on 14
    cores against 1.52 s on one. Formatting the floats sat in the serial write
    loop on top of that. Both are fixed; rows are now batched, and each worker
    formats its own.

Set `-t` when you are sharing a machine, or when the run is part of a pipeline
that is already parallel. Leaving it unset is right for a dedicated node.

## Memory

Nothing in distree holds the matrix. The peak is:

```
  LCA table       M × (⌈log₂ M⌉ + 1) × 8 bytes
  Node depths     M × 16 bytes
  The tree        roughly M × 72 bytes, plus the labels
  Batch buffer    about 15 MB, fixed
  Output buffer   1 MB
```

For a bifurcating tree, `M ≈ 2N`. Worked through:

| Tips | Nodes | LCA table | Total, roughly |
|--:|--:|--:|--:|
| 10,000 | 20,000 | 2.6 MB | 35 MB |
| 100,000 | 200,000 | 30 MB | 100 MB |
| 1,000,000 | 2,000,000 | 350 MB | 600 MB |

The 100,000-tip figure is measured; the million-tip one is the arithmetic.

The LCA table is the term that grows, and it is stored as plain `usize` with a
sentinel rather than `Option<usize>`, which halves it. The batch buffer holds
about a million cells' worth of formatted rows regardless of tree size, which is
what keeps the total flat as the matrix grows.

## Output size

The output is usually larger than anything in memory, and precision is the lever
on it. At `-p 10` a cell is 12 to 14 bytes; at `-p 6` it is 8 to 10.

| Tips | Square, `-p 10` | Square, `-p 6` | `--lower`, `-p 6` |
|--:|--:|--:|--:|
| 1,000 | 13 MB | 9 MB | 4.5 MB |
| 10,000 | 1.3 GB | 900 MB | 450 MB |
| 50,000 | 33 GB | 22 GB | 11 GB |

Two things follow. Use `--lower` when the reader accepts it, which halves the
file. And do not ask for more decimals than the branch lengths carry: `-p 6`
against the default `-p 10` is a third off the file for no loss on any realistic
tree.

Rows are written as they are computed, so a pipeline that consumes them as they
arrive never holds the whole matrix either:

```bash
distree big.nwk -p 6 | gzip > distances.tsv.gz
```

## Where the ceiling is

The quadratic term decides. Around 50,000 tips a square matrix is tens of
gigabytes and the run is minutes; at 100,000 it is hundreds of gigabytes and
the file, not distree, is the problem.

If the tips are more than that, the question is usually not really "what is the
whole matrix". Reduce the tree first, with something like
[axetree](https://github.com/PathoGenOmics-Lab/axetree), and compute the matrix
over the subset.

## Reproducibility

The row order is the tip labels sorted, and the cell order within a row is the
same sort, neither of which depends on the thread count or on how the Newick was
written. Two runs of the same tree produce byte-identical output, on any number
of cores.

## Next

- [How it works](algorithm.md), for what these stages actually do.
- [Output](../guide/output.md), for the formats and the precision flag.
