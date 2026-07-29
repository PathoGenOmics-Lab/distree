# Output

The matrix goes to stdout unless `-o` names a file. Warnings go to stderr, so a
piped run stays clean:

```bash
distree tree.nwk | head -1
distree tree.nwk -o distances.tsv
```

## The square TSV

The default layout is a full tab-separated matrix:

```
	LeafA	LeafB	LeafC	LeafD
LeafA	0.000	5.200	3.600	4.000
LeafB	5.200	0.000	4.900	5.300
LeafC	3.600	4.900	0.000	2.000
LeafD	4.000	5.300	2.000	0.000
```

- The first line is a header of tip labels, preceded by an **empty cell** so the
  labels sit over the data columns rather than one to the left. This is what
  `read.table(header=TRUE, row.names=1)` and `pandas.read_csv(index_col=0)`
  expect.
- Every following line begins with a tip label and holds one value per tip.
- **Tips are sorted alphabetically**, by byte order, in both the rows and the
  columns. The order is the same in every run and does not depend on how the
  tree was written, so two matrices from two trees over the same tips line up
  cell for cell.
- The matrix is symmetric and the diagonal is exactly zero, except under
  [`--lmm`](distances.md#variance-covariance), where the diagonal is each tip's
  root-to-tip length.

Reading it back:

=== "R"

    ```r
    m <- as.matrix(read.table("distances.tsv", header = TRUE,
                              row.names = 1, sep = "\t", check.names = FALSE))
    ```

    `check.names = FALSE` keeps labels like `Sample-1` and `2024_isolate` intact;
    without it R rewrites them.

=== "Python"

    ```python
    import pandas as pd
    m = pd.read_csv("distances.tsv", sep="\t", index_col=0)
    ```

=== "Shell"

    ```bash
    # the distance between two named tips
    awk -F'\t' -v a=LeafA -v b=LeafC '
      NR==1 { for (i=2; i<=NF; i++) if ($i==b) col=i; next }
      $1==a { print $col }' distances.tsv
    ```

## PHYLIP lower triangle

```bash
distree tree.nwk --lower -p 4
```

```
4
LeafA
LeafB	5.2000
LeafC	3.6000	4.9000
LeafD	4.0000	5.3000	2.0000
```

A taxa count on the first line, then one row per taxon holding its label and the
distances to every taxon **above** it in the ordering. No diagonal, and none of
the mirrored half.

This is the format PHYLIP's `neighbor` and `fitch` read, and what Mash and a
number of clustering tools expect. It also halves the output: a 50,000-tip tree
is 2.5 billion cells as a square matrix and 1.25 billion as a lower triangle.

The first data row is the label alone, with no values after it, which is correct
and is what those readers expect.

!!! note "Relaxed PHYLIP"

    Fields are tab-separated and labels are written in full. Strict PHYLIP pads
    every name to exactly ten characters and truncates anything longer, which
    would collide two tips whose names share a prefix. distree writes the
    relaxed form that whitespace-splitting readers, including modern PHYLIP
    builds, accept.

## NumPy arrays

```bash
distree tree.nwk --npy -o distances.npy
```

Text is an expensive way to move a large matrix. At twelve decimals a cell is
fourteen bytes and producing them is most of the run; as a 64-bit float it is
eight bytes, exact, and costs nothing to write. On an 8,000-tip tree:

| | Time | Size |
|:--|--:|--:|
| `--lower -p 12` | 0.26 s | 458 MB |
| `--lower --npy` | 0.10 s | 244 MB |
| `-p 12` | 0.56 s | 916 MB |
| `--npy` | 0.19 s | 488 MB |

`.npy` has nowhere to put labels, so they go to `<FILE>.labels.txt` in row
order. That is why `--npy` needs `-o` rather than writing to stdout.

```python
import numpy as np
m = np.load("distances.npy")
labels = open("distances.npy.labels.txt").read().split()
```

With `--lower`, the array is the **condensed vector** rather than PHYLIP's
lower triangle, because those are two different triangles: SciPy reads the
upper one row by row. It goes straight into `scipy`:

```python
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

v = np.load("condensed.npy")          # distree tree.nwk --lower --npy -o condensed.npy
z = linkage(v, method="single")       # takes the condensed form directly
full = squareform(v)                  # or expand it to the square matrix
```

`-p` has no effect here; a 64-bit float carries the value exactly.

## Precision

```bash
distree tree.nwk -p 6
```

`-p` sets the number of decimal places, defaulting to 10 and capped at 30. It
applies to patristic and `--lmm`; `--topology` writes integers and ignores it.

A 64-bit float carries about 17 significant decimal digits, so anything past
that is zero padding, and the cap is there to keep an accidental `-p 50000000`
from being an error inside the formatter rather than an error about the flag.

Precision is also the main lever on output size. At `p=10` each cell is 12 to 14
bytes; at `p=6` it is 8 to 10. Across a 20,000-tip square matrix that is the
difference between roughly 5 GB and 3.5 GB.

## Writing to a file

```bash
distree tree.nwk -o distances.tsv
```

The file is not created until the tree has parsed, midpoint rooting has run and
the tip labels have been checked, so a run that fails on a malformed tree leaves
whatever was at that path untouched.

Write errors are reported. A full disk fails with the path in the message and a
non-zero exit code, rather than leaving a truncated matrix behind and reporting
success.

## Piping

`distree tree.nwk | head` closes the pipe early, which is a normal way to look
at the first few rows. distree treats that as the reader having seen enough and
exits cleanly rather than reporting a broken pipe.

Because rows are written as they are computed, the first row of a very large
matrix appears long before the last one. A pipeline that consumes rows as they
arrive does not wait for the whole matrix:

```bash
distree big.nwk -p 6 | awk -F'\t' 'NR>1 { for (i=2; i<=NF; i++) if ($i+0 < 1e-5 && i-1 != NR-1) print $1, i }'
```

## Exit codes

| Code | Meaning |
|:--|:--|
| 0 | The matrix was written and flushed. Also the code when a downstream pipe closed early |
| 1 | Anything went wrong: an unreadable file, a malformed tree, a bad flag, a failed write. The message is on stderr and begins with `Error:` |

## Next

- [Recipes](../recipes.md), for what to do with the matrix once you have it.
- [CLI reference](cli.md), for every flag in one table.
