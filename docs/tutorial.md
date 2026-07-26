# Tutorial

A worked outbreak investigation, start to finish. Seven *Mycobacterium
tuberculosis* isolates, one tree, and the question every such dataset opens
with: **which of these are transmission links, and which are unrelated cases
that happen to be in the same clinic?**

Everything here runs in a few seconds and needs nothing but distree, so you can
paste along. R and Python appear in tabs; pick whichever you use.

## 1. The tree

Save this as `outbreak.nwk`. It is what a maximum-likelihood tool would hand
you: branch lengths in **substitutions per site**, in scientific notation.

```
((((TB_001:2.266786e-07,TB_002:2.266786e-07):2.266786e-07,TB_003:1.586750e-06):8.613788e-06,(TB_004:3.400179e-07,TB_005:3.400179e-07):9.860520e-06):1.360072e-05,(TB_006:1.813429e-05,TB_007:2.153447e-05):1.586750e-05);
```

```bash
distree --version
```

## 2. The first matrix

```bash
distree outbreak.nwk -p 9
```

```
	TB_001	TB_002	TB_003	TB_004	TB_005	TB_006	TB_007
TB_001	0.000000000	0.000000453	0.000002040	0.000019268	0.000019268	0.000056670	0.000060070
TB_002	0.000000453	0.000000000	0.000002040	0.000019268	0.000019268	0.000056670	0.000060070
TB_003	0.000002040	0.000002040	0.000000000	0.000020401	0.000020401	0.000057803	0.000061203
TB_004	0.000019268	0.000019268	0.000020401	0.000000000	0.000000680	0.000057803	0.000061203
TB_005	0.000019268	0.000019268	0.000020401	0.000000680	0.000000000	0.000057803	0.000061203
TB_006	0.000056670	0.000056670	0.000057803	0.000057803	0.000057803	0.000000000	0.000039669
TB_007	0.000060070	0.000060070	0.000061203	0.000061203	0.000061203	0.000039669	0.000000000
```

That is the whole tool, really. Each cell is the sum of the branch lengths on
the path between two tips: a **patristic distance**. The rest of this page is
about reading it.

Some structure is already visible. TB_001 and TB_002 are an order of magnitude
closer to each other than to anything else, TB_004 and TB_005 likewise, and
TB_006 and TB_007 sit far from everyone.

!!! note "Why the empty first cell"

    The header line starts with a tab, so the labels sit over the data columns
    rather than one to the left. That is what `read.table(row.names = 1)` and
    `pandas.read_csv(index_col=0)` expect, and it is why the matrix loads
    without any fiddling in [step 4](#4-load-it).

## 3. Substitutions to SNPs

Nobody thinks in substitutions per site. Multiply by the length of the reference
the tree was built against, and the numbers become countable differences.
H37Rv is 4,411,532 bases:

```bash
distree outbreak.nwk -p 12 | awk 'NR == 1 { print; next }
  { printf "%s", $1; for (i = 2; i <= NF; i++) printf "\t%.0f", $i * 4411532; print "" }'
```

```
	TB_001	TB_002	TB_003	TB_004	TB_005	TB_006	TB_007
TB_001	0	2	9	85	85	250	265
TB_002	2	0	9	85	85	250	265
TB_003	9	9	0	90	90	255	270
TB_004	85	85	90	0	3	255	270
TB_005	85	85	90	3	0	255	270
TB_006	250	250	255	255	255	0	175
TB_007	265	265	270	270	270	175	0
```

Now it reads. TB_001 and TB_002 differ by **2 SNPs**. TB_003 is **9** from both.
TB_004 and TB_005 are **3** apart but **85** from the first group. TB_006 and
TB_007 are hundreds of SNPs from everything, including each other.

!!! warning "Use enough precision before you scale"

    `-p 12` above, not the default `-p 10`. Multiplying by four million
    multiplies the rounding too: at `-p 6` a distance of 2 SNPs would round to 0
    before you ever saw it. Scale from the most precise output you can, and
    round at the end.

## 4. Load it

```bash
distree outbreak.nwk -p 12 -o distances.tsv
```

=== "R"

    ```r
    m <- as.matrix(read.table("distances.tsv", header = TRUE, row.names = 1,
                              sep = "\t", check.names = FALSE))
    snps <- m * 4411532
    round(snps)
    ```

    `check.names = FALSE` matters. Without it R rewrites labels like `TB-001`
    or `2024_isolate` and they stop matching your metadata.

=== "Python"

    ```python
    import pandas as pd

    m = pd.read_csv("distances.tsv", sep="\t", index_col=0)
    snps = m * 4411532
    print(snps.round())
    ```

## 5. Find the clusters

Public health uses SNP thresholds. Two are standard for *M. tuberculosis*:
**5 SNPs** for recent transmission, and **12 SNPs** for a plausible
epidemiological link.[^walker] Single linkage is the right rule, because it
groups isolates that are within the threshold of *someone* in the group, which
is what a chain of transmission looks like.

=== "R"

    ```r
    clusters <- function(h) split(rownames(snps),
                                  cutree(hclust(as.dist(snps), method = "single"), h = h))
    clusters(5)
    clusters(12)
    ```

=== "Python"

    ```python
    import collections
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    z = linkage(squareform(snps.values, checks=False), method="single")

    def clusters(t):
        out = collections.defaultdict(list)
        for name, c in zip(snps.index, fcluster(z, t=t, criterion="distance")):
            out[c].append(name)
        return list(out.values())

    print(clusters(5))
    print(clusters(12))
    ```

```
threshold  5 SNPs:  TB_001,TB_002 | TB_003 | TB_004,TB_005 | TB_006 | TB_007
threshold 12 SNPs:  TB_001,TB_002,TB_003 | TB_004,TB_005 | TB_006 | TB_007
```

**The threshold is the finding.** TB_003 is 9 SNPs away: outside recent
transmission, inside a plausible link. Whether it belongs in the cluster is an
epidemiological question, not a computational one, and the matrix will not
answer it for you. Report which threshold you used.

Note also what did *not* change. TB_004 and TB_005 are a pair at either
threshold, and TB_006 and TB_007 are singletons at either. Findings that survive
both thresholds are the ones to lead with.

## 6. Who is closest to whom

Before any clustering, the cheapest useful question is each isolate's nearest
neighbour:

```bash
distree outbreak.nwk -p 12 | awk -F'\t' '
  NR == 1 { for (i = 2; i <= NF; i++) name[i] = $i; next }
  { best = ""; bestd = 1e308
    for (i = 2; i <= NF; i++)
      if (i - 1 != NR - 1 && $i + 0 < bestd) { bestd = $i + 0; best = name[i] }
    printf "%s\t%s\t%.0f\n", $1, best, bestd * 4411532 }'
```

```
TB_001	TB_002	2
TB_002	TB_001	2
TB_003	TB_001	9
TB_004	TB_005	3
TB_005	TB_004	3
TB_006	TB_007	175
TB_007	TB_006	175
```

The `i - 1 != NR - 1` skips the diagonal, which is otherwise always the closest.

The last two rows are the interesting ones. TB_006 and TB_007 are each other's
nearest neighbours and still 175 SNPs apart, which is not a transmission link;
it just means nobody closer was sampled.

## 7. When branch lengths are not the point

Sometimes the tree's branch lengths are not comparable, or there are none at
all. `--topology` counts edges instead:

```bash
distree outbreak.nwk --topology
```

```
	TB_001	TB_002	TB_003	TB_004	TB_005	TB_006	TB_007
TB_001	0	2	3	5	5	6	6
TB_002	2	0	3	5	5	6	6
TB_003	3	3	0	4	4	5	5
TB_004	5	5	4	0	2	5	5
TB_005	5	5	4	2	0	5	5
TB_006	6	6	5	5	5	0	2
TB_007	6	6	5	5	5	2	0
```

Compare it with the SNP matrix and the difference in what they measure is
plain. Topologically TB_006 and TB_007 are 2 apart, the same as TB_001 and
TB_002; in SNPs they are 175 apart against 2. Edge counts describe the shape of
the tree, not the divergence across it, so use them for tree-shape questions
and never as a proxy for relatedness when you have usable branch lengths.

## 8. The comparative-analysis matrix

A different question: you have a trait for each isolate and want to test whether
it associates with something, without the shared ancestry inflating your
significance. That needs the variance-covariance matrix, and here rooting is
part of the answer, so root deliberately:

```bash
distree outbreak.nwk --midpoint --lmm -p 8 -o varcovar.tsv
```

```
	TB_001	TB_002	TB_003	TB_004	TB_005	TB_006	TB_007
TB_001	0.00002947	0.00002924	0.00002901	0.00002040	0.00002040	0.00000000	0.00000000
TB_002	0.00002924	0.00002947	0.00002901	0.00002040	0.00002040	0.00000000	0.00000000
...
```

Entry `(i, j)` is how much evolutionary history the two isolates share: the
distance from the root down to their common ancestor. The diagonal is each
isolate's own root-to-tip length. Pairs whose ancestor is the root share
nothing and score 0.

Three things to check before you fit anything with it:

- **It is not a distance matrix.** The diagonal is not zero, and bigger means
  more similar. Do not feed it to something expecting distances.
- **Reindex it.** Rows and columns are in alphabetical label order, which is
  almost certainly not the order of your trait table.
- **Root on purpose.** `--midpoint` is one defensible choice; an outgroup is
  another. Unlike the patristic matrix, this one changes if you change your
  mind. See [Midpoint rooting](guide/rooting.md).

## 9. Making it bigger

Seven isolates is a tutorial. Seven thousand is a Tuesday, and two flags matter
then.

```bash
distree big.nwk --lower -p 9 -t 8 -o distances.phy
```

`--lower` writes the PHYLIP lower triangle, halving the file:

```
7
TB_001
TB_002	0.000000453
TB_003	0.000002040	0.000002040
TB_004	0.000019268	0.000019268	0.000020401
TB_005	0.000019268	0.000019268	0.000020401	0.000000680
TB_006	0.000056670	0.000056670	0.000057803	0.000057803	0.000057803
TB_007	0.000060070	0.000060070	0.000061203	0.000061203	0.000061203	0.000039669
```

`-t` caps the threads; leave it off on a machine of your own and set it when
you are sharing one.

Rows are written as they are computed, so a large run streams rather than
building the whole matrix first. Compress on the way out, or filter as it goes:

```bash
distree big.nwk -p 12 --lower | gzip > distances.phy.gz
```

For what a run costs at various sizes, see
[Performance](how-it-works/performance.md).

## 10. When it refuses

distree would rather stop than hand back a matrix that looks right and is not.
The three you are most likely to meet:

```
Error: /path/tree.nwk is gzip-compressed. distree reads plain text; decompress it on the way in:
    gunzip -c /path/tree.nwk | distree -
```

```
Error: Failed to parse Newick tree: Unexpected content at position 84: the file appears to
hold more than one tree. distree processes a single Newick tree; split the file first.
```

```
Error: Duplicate leaf name 'TB_001' found. Leaf names must be unique.
```

A duplicate label is worth pausing over. It usually means the same sample went
into the tree twice under two accessions, and no matrix can tell the two rows
apart afterwards. The full list is in [Input](guide/input.md#what-is-rejected).

## What to take away

- The **patristic** matrix is what almost every downstream question wants, and
  it does not care how the tree is rooted.
- **Scale it once**, from high precision, and round at the end.
- **The threshold is a decision**, not an output. State it.
- **`--lmm` is a different object**: a covariance, not a distance, and rooting
  changes it.
- **`--topology` measures shape**, not divergence.

## Next

- [Recipes](recipes.md) for PCoA, PGLS and neighbour-joining, at the same length
  as this section but without the story.
- [Distance modes](guide/distances.md) for the three modes in full.
- [CLI reference](guide/cli.md) for every flag.

[^walker]:
    Walker TM, Ip CLC, Harrell RH, et al. *Whole-genome sequencing to delineate
    Mycobacterium tuberculosis outbreaks: a retrospective observational study.*
    The Lancet Infectious Diseases. 2013;13(2):137-146.
