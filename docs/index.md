---
hide:
  - navigation
---

<div class="hero" markdown>

![distree](assets/distree_wordmark.svg){ .hero-logo }

# distree

Pairwise distance matrices from a phylogeny, in Rust.

</div>

Once a tree is built, most of the questions that follow are pairwise. Which
isolates sit within five substitutions of each other. What the ordination looks
like. What covariance structure a comparative model should assume. All of them
want the same thing first: the distance between every pair of tips.

Computing that naively means walking the tree once per pair, and holding an
N × N matrix while you do it. distree instead spends `O(M log M)` up front on a
binary-lifting structure that answers *"where do these two tips meet?"* in
`O(log M)` time, then streams the matrix out one row at a time. A row is
computed in parallel across cores and written before the next one starts, so
memory tracks the number of tips, not the square of it.

<div class="grid cards" markdown>

-   :material-download: **[Install](getting-started/installation.md)**

    Bioconda, a prebuilt binary, or `cargo build --release`. No runtime
    dependencies.

-   :material-rocket-launch: **[Quick start](getting-started/quickstart.md)**

    The handful of commands that cover most of what people want a matrix for.

-   :material-file-tree: **[Input](guide/input.md)**

    What Newick distree accepts, what it rejects, and why it rejects rather
    than guesses.

-   :material-ruler: **[Distance modes](guide/distances.md)**

    Patristic, topological and variance-covariance, on one worked example.

-   :material-format-vertical-align-center: **[Midpoint rooting](guide/rooting.md)**

    When the root changes the answer, and when it cannot.

-   :material-table: **[Output](guide/output.md)**

    The TSV layout, PHYLIP lower triangle, precision, and reading it back in.

-   :material-school: **[Tutorial](tutorial.md)**

    A worked outbreak investigation, from the tree to transmission clusters.

-   :material-cog: **[How it works](how-it-works/algorithm.md)**

    The parser, the LCA structure, and the streaming loop.

-   :material-speedometer: **[Performance](how-it-works/performance.md)**

    What the run costs in time and memory, and where the ceiling is.

</div>

## In one command

```bash
distree tree.nwk --lower -p 6 -t 8 -o distances.phy
```

Patristic distances for every pair of tips, as a PHYLIP lower triangle, six
decimal places, across eight threads. Swap `--topology` for edge counts,
`--lmm` for the variance-covariance matrix a comparative model wants, and drop
`--lower` for a full square TSV that `read.table` and `pandas` open directly.

## What it does

| | |
|:--|:--|
| **Patristic** | The sum of branch lengths on the path between two tips |
| **Topological** | The number of edges on that path, ignoring branch lengths |
| **Variance-covariance** | The root-to-MRCA distance for each pair, the `C` matrix of a PGLS or a phylogenetic mixed model |
| **Rooting** | Optional midpoint rooting, which matters for `--lmm` and cannot matter for the other two |
| **Formats** | Square TSV with labels, a PHYLIP lower triangle, or a NumPy array |
| **Input** | Newick from a file or stdin, plain or gzipped, with quoted labels, polytomies, NHX and BEAST comments, and UTF-8 tip names |
| **Subsets** | `--taxa` restricts the matrix to a list of tips, with the distances the full tree gives |
| **Scale** | Memory grows with the number of tips, not with the matrix; rows stream out as they are computed |
| **Parallelism** | Each row is computed across all cores, or as many as `-t` allows |

## What it does not do

distree works on the tree you give it. It does not build one, it does not read
an alignment, and it does not cluster the matrix for you.

It also does not guess. A truncated tree, a file holding several trees, an
unclosed quote or a label with a tab in it are all rejected with a message
naming the position, rather than turned into a matrix that looks right and is
not. The one case that is reported and carried through rather than refused is a
negative branch length, because neighbour-joining trees legitimately have them
and the resulting negative distances are information, not noise.

## Citation

If you use distree, please cite the archived release. See
[Citation](about/citation.md) for the DOI and the version to record.
