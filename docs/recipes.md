# Recipes

Short, complete examples of what people actually do with the matrix.

## Ordination: PCoA and MDS

=== "R"

    ```r
    m <- as.matrix(read.table("distances.tsv", header = TRUE,
                              row.names = 1, sep = "\t", check.names = FALSE))
    fit <- cmdscale(as.dist(m), k = 2, eig = TRUE)
    plot(fit$points, pch = 19, xlab = "PCo1", ylab = "PCo2")
    text(fit$points, labels = rownames(m), pos = 3, cex = 0.6)

    # variance explained by the first two axes
    round(100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0]), 1)
    ```

=== "Python"

    ```python
    import pandas as pd
    from sklearn.manifold import MDS

    m = pd.read_csv("distances.tsv", sep="\t", index_col=0)
    coords = MDS(n_components=2, dissimilarity="precomputed",
                 random_state=0, normalized_stress="auto").fit_transform(m.values)
    ```

    For UMAP, pass the same matrix with `metric="precomputed"`.

Generate it with plenty of precision; ordination is sensitive to the small
distances:

```bash
distree tree.nwk -p 10 -o distances.tsv
```

## Transmission clusters by a SNP threshold

A patristic distance from a tree built on an alignment is in substitutions per
site. Multiply by the alignment length to get SNP-scale distances, then cut the
matrix at your threshold:

=== "R"

    ```r
    m <- as.matrix(read.table("distances.tsv", header = TRUE,
                              row.names = 1, sep = "\t", check.names = FALSE))
    snps <- m * 4411532          # H37Rv genome length
    clusters <- cutree(hclust(as.dist(snps), method = "single"), h = 5)
    split(names(clusters), clusters)
    ```

=== "Python"

    ```python
    import pandas as pd
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    m = pd.read_csv("distances.tsv", sep="\t", index_col=0) * 4411532
    z = linkage(squareform(m.values, checks=False), method="single")
    labels = fcluster(z, t=5, criterion="distance")
    ```

Single linkage is the right choice for a transmission threshold: it clusters
samples that are within the threshold **of someone** in the cluster, which is
what a chain of transmission looks like.

## PGLS and phylogenetic mixed models

`--lmm` writes the `C` matrix these models want. Root the tree deliberately
first, because [rooting is part of the answer](guide/rooting.md):

```bash
distree tree.nwk --midpoint --lmm -p 10 -o varcovar.tsv
```

```r
C <- as.matrix(read.table("varcovar.tsv", header = TRUE,
                          row.names = 1, sep = "\t", check.names = FALSE))
traits <- read.csv("traits.csv", row.names = 1)
C <- C[rownames(traits), rownames(traits)]      # align the orderings

library(nlme)
fit <- gls(y ~ x, data = traits,
           correlation = corSymm(C[lower.tri(C)] / max(C), fixed = TRUE))
```

Two things to check before fitting. The row and column order is alphabetical by
tip label, which is almost certainly not the order of your trait table, so
reindex rather than assume. And a non-constant diagonal means the tree is not
ultrametric, which several comparative methods assume.

## Neighbour-joining from an existing tree

PHYLIP's `neighbor` reads the lower triangle directly:

```bash
distree tree.nwk --lower -p 6 -o infile
printf 'L\nY\n' | neighbor
```

The `L` selects lower-triangular input; `neighbor` writes `outtree` and
`outfile`.

## The closest relative of every sample

```bash
distree tree.nwk -p 8 | awk -F'\t' '
  NR == 1 { for (i = 2; i <= NF; i++) name[i] = $i; next }
  {
    best = ""; bestd = 1e308
    for (i = 2; i <= NF; i++) if (i - 1 != NR - 1 && $i + 0 < bestd) { bestd = $i + 0; best = name[i] }
    printf "%s\t%s\t%s\n", $1, best, bestd
  }'
```

Skipping `i - 1 == NR - 1` skips the diagonal, which is otherwise always the
closest.

## Comparing two trees over the same tips

Because the tip ordering is alphabetical in every run, two matrices over the
same tip set line up cell for cell:

```bash
distree iqtree.nwk  -p 10 --lower -o a.phy
distree raxml.nwk   -p 10 --lower -o b.phy
```

```python
import numpy as np, pandas as pd
from scipy.stats import pearsonr

def lower(path):
    vals = []
    with open(path) as fh:
        next(fh)                                  # taxa count
        for line in fh:
            vals += [float(v) for v in line.split("\t")[1:]]
    return np.array(vals)

a, b = lower("a.phy"), lower("b.phy")
print(pearsonr(a, b))
```

A high correlation with a slope away from 1 means the two trees agree on the
shape and disagree on the rate.

## Cladograms and unreliable branch lengths

A tree with no branch lengths gives an all-zero patristic matrix, and says so.
Count edges instead:

```bash
distree cladogram.nwk --topology -o topo.tsv
```

The same applies when the lengths exist but are not comparable across the tree,
for instance from concatenated loci with different rates.

## Very large trees

Rows are written as they are computed, so nothing needs the whole matrix in one
piece. Compress on the way out:

```bash
distree big.nwk -p 6 --lower | gzip > distances.phy.gz
```

Or filter as it streams, keeping only the pairs that matter:

```bash
distree big.nwk -p 8 --lower | awk -F'\t' '
  NR == 1 { next }
  { for (i = 2; i <= NF; i++) if ($i + 0 < 1e-5) print $1, i - 1, $i }'
```

Past 50,000 tips or so the file becomes the constraint rather than the
computation. See [Performance](how-it-works/performance.md#where-the-ceiling-is).

## Batch runs

One matrix per tree, in parallel, one core each:

```bash
ls trees/*.nwk | xargs -P 8 -I{} sh -c \
  'distree "$1" -t 1 -p 6 --lower -o "${1%.nwk}.phy"' _ {}
```

`-t 1` matters here: without it every distree would try to use every core and
they would fight each other.

## Splitting a multi-tree file

distree takes one tree per file, and refuses a file holding several rather than
silently using the first:

```bash
split -l 1 -d posterior.trees rep_ --additional-suffix=.nwk
ls rep_*.nwk | xargs -P 8 -I{} sh -c \
  'distree "$1" -t 1 -p 6 --lower -o "${1%.nwk}.phy"' _ {}
```
