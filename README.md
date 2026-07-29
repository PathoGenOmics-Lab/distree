<p align="center">
  <img src=".github/logo/distree.png" title="pdistree" style="width:750px; height: auto;">
</p>

# distree
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/PathoGenOmics-Lab/distree/blob/main/LICENSE)
[![distree](https://img.shields.io/badge/distree-rust-%23ff8000)](https://github.com/PathoGenOmics-Lab/distree)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/distree.svg)](https://anaconda.org/bioconda/distree)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/distree/badges/version.svg)](https://anaconda.org/bioconda/distree)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?)](https://github.com/PathoGenOmics-Lab)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.16811766-%23ff0077)](https://doi.org/10.5281/zenodo.16811766)

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>

`distree` reads a phylogenetic tree in Newick format and writes the pairwise
distance between every pair of tips: **patristic** (summed branch lengths),
**topological** (edge counts), or the **variance-covariance** matrix a
comparative model wants. Rows stream out as they are computed, so memory grows
with the number of tips rather than with the square of it.

📖 **Full documentation: <https://pathogenomics-lab.github.io/distree/>**

| Page | Description |
|:-----|:------------|
| **[Quick start](https://pathogenomics-lab.github.io/distree/getting-started/quickstart/)** | The handful of commands most runs actually use |
| **[Tutorial](https://pathogenomics-lab.github.io/distree/tutorial/)** | A worked outbreak investigation, tree to transmission clusters |
| **[Input](https://pathogenomics-lab.github.io/distree/guide/input/)** | What Newick is accepted, what is rejected, and why |
| **[Distance modes](https://pathogenomics-lab.github.io/distree/guide/distances/)** | Patristic, topological and var-covar, on one worked example |
| **[Midpoint rooting](https://pathogenomics-lab.github.io/distree/guide/rooting/)** | When the root changes the answer, and when it cannot |
| **[Output](https://pathogenomics-lab.github.io/distree/guide/output/)** | TSV, PHYLIP lower triangle, NumPy arrays, precision |
| **[CLI reference](https://pathogenomics-lab.github.io/distree/guide/cli/)** | Every flag, every default, every message |
| **[How it works](https://pathogenomics-lab.github.io/distree/how-it-works/algorithm/)** | The parser, the LCA structure, the streaming loop |
| **[Performance](https://pathogenomics-lab.github.io/distree/how-it-works/performance/)** | Measured speed, memory and thread scaling |
| **[Recipes](https://pathogenomics-lab.github.io/distree/recipes/)** | PCoA, transmission clusters, PGLS, neighbour-joining |

---

## Installation

```bash
# Bioconda
conda install -c bioconda distree

# From source
cargo build --release
```

Prebuilt binaries for Linux and macOS, x86-64 and arm64, are on the
[releases page](https://github.com/PathoGenOmics-Lab/distree/releases).

## Quick start

```bash
# Patristic distances, as a tab-separated matrix
distree tree.nwk -o distances.tsv

# Six decimals across eight threads, PHYLIP lower triangle
distree tree.nwk --lower -p 6 -t 8 -o distances.phy

# Edge counts instead of branch lengths
distree tree.nwk --topology -o topology.tsv

# The variance-covariance matrix for PGLS, midpoint-rooted
distree tree.nwk --midpoint --lmm -o varcovar.tsv

# A NumPy array: exact, half the size of text, and 2.6x faster
distree tree.nwk --npy -o distances.npy

# The condensed vector scipy.cluster.hierarchy reads directly
distree tree.nwk --lower --npy -o condensed.npy

# Just the samples in a list, from a gzipped tree
distree tree.nwk.gz --taxa cohort.txt --stats -o cohort.tsv

# From stdin
gunzip -c tree.nwk.gz | distree - -o distances.tsv
```

## What it does

| | |
|:--|:--|
| **Patristic** | The sum of branch lengths on the path between two tips |
| **Topological** | The number of edges on that path, ignoring branch lengths |
| **Variance-covariance** | The root-to-MRCA distance for each pair, the `C` matrix of a PGLS or phylogenetic mixed model |
| **Rooting** | Optional midpoint rooting, which matters for `--lmm` and cannot matter for the other two |
| **Formats** | Square TSV, PHYLIP lower triangle, or a NumPy `.npy` array |
| **Input** | Newick from a file or stdin, plain or gzipped, with quoted labels, polytomies, NHX and BEAST comments, and UTF-8 tip names |
| **Subsets** | `--taxa` restricts the matrix to a list of tips, with the distances the full tree gives |
| **Scale** | Rows stream out as they are computed; memory tracks the tips, not the matrix |
| **Parallelism** | Distances and their formatting both spread across cores |

It works on the tree you give it: it does not build one, does not read an
alignment, and does not cluster the matrix for you. It also does not guess. A
truncated tree, a file holding several trees, an unclosed quote or a label with
a tab in it are rejected with the position, rather than turned into a matrix
that looks right and is not.

## Performance

Balanced trees, `-p 6 --lower`, release build on an Apple M4 Pro with 14 cores.

| Tips | Cells | Time | Peak memory |
|--:|--:|--:|--:|
| 1,000 | 0.5 M | 0.00 s | 11 MB |
| 8,000 | 32 M | 0.13 s | 32 MB |
| 20,000 | 200 M | 1.02 s | 35 MB |

Memory flattens out because nothing holds the matrix: past a few thousand tips
it is the LCA table plus one batch of rows. `--npy` is another 2.6x on top, at
full precision and half the file size. Full numbers, thread scaling and the
memory arithmetic are
[here](https://pathogenomics-lab.github.io/distree/how-it-works/performance/).

## Correctness

Distances are cross-validated against R's
[ape](https://cran.r-project.org/package=ape) over generated trees:
`cophenetic.phylo` for patristic, `vcv.phylo` for variance-covariance,
`cophenetic.phylo` over unit branch lengths for edge counts, and
`phangorn::midpoint` for the rooting. All four agree to under 1e-9, which is the
text round-trip, and exactly for edge counts. The parser and the pipeline are
fuzzed. See
[Contributing](https://pathogenomics-lab.github.io/distree/about/contributing/).

## Citation

> Ruiz-Rodriguez P, Coscolla M. *distree: distance matrices from a phylogeny.*
> PathoGenOmics Lab. [doi:10.5281/zenodo.16811766](https://doi.org/10.5281/zenodo.16811766)

Please record the version and the mode you used: a patristic matrix, a
topological one and a variance-covariance matrix are three different objects.

## License

[GPL-3.0](https://github.com/PathoGenOmics-Lab/distree/blob/main/LICENSE).

---
<h2 id="contributors" align="center">

✨ [Contributors](https://github.com/PathoGenOmics-Lab/distree/graphs/contributors)
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
distree is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->
