# Citation

distree is archived on Zenodo. Cite the release you used:

> Ruiz-Rodriguez P, Coscolla M. *distree: distance matrices from a phylogeny.*
> PathoGenOmics Lab. [doi:10.5281/zenodo.16811766](https://doi.org/10.5281/zenodo.16811766)

The repository carries a [`CITATION.cff`](https://github.com/PathoGenOmics-Lab/distree/blob/main/CITATION.cff),
so GitHub's "Cite this repository" button and most reference managers will
produce the entry for you.

Please record the version (`distree --version`) and the mode you used. A
patristic matrix, a topological one and a variance-covariance matrix are three
different objects, and "distances from the tree" does not say which.

## Authors

**Paula Ruiz-Rodriguez** and **Mireia Coscolla**, I²SysBio, University of
Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health,
Valencia, Spain.

## Background

The quantities distree computes are standard, and the terms are worth pinning
down.

**Patristic distance** is the sum of branch lengths along the path between two
tips. The name is Farris's.

> Farris JS. *A Method for Computing Wagner Trees.* Systematic Zoology.
> 1970;19(1):83-92.

**The variance-covariance matrix** that `--lmm` writes is the covariance
structure implied by Brownian-motion trait evolution on a tree, which is what
makes it the `C` of a phylogenetic generalised least squares fit.

> Felsenstein J. *Phylogenies and the Comparative Method.* The American
> Naturalist. 1985;125(1):1-15.

> Grafen A. *The phylogenetic regression.* Philosophical Transactions of the
> Royal Society B. 1989;326(1233):119-157.

**Midpoint rooting** places the root at the middle of the longest tip-to-tip
path.

> Farris JS. *Estimating Phylogenetic Trees from Distance Matrices.* The
> American Naturalist. 1972;106(951):645-668.

**Binary lifting** is the technique behind the `O(log n)` ancestor queries. The
sparse-table formulation of the least-common-ancestor problem goes back to:

> Bender MA, Farach-Colton M. *The LCA Problem Revisited.* LATIN 2000, Lecture
> Notes in Computer Science 1776, pp. 88-94.

## Formats and downstream tools

The `--lower` output is PHYLIP's lower-triangular distance format:

> Felsenstein J. *PHYLIP: Phylogeny Inference Package (Version 3.2).* Cladistics.
> 1989;5:164-166.

## License

[GPL-3.0](https://github.com/PathoGenOmics-Lab/distree/blob/main/LICENSE).
