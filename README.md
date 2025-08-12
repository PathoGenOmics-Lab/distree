<p align="center">
  <img src=".github/logo/distree.png" title="pdistree" style="width:750px; height: auto;">
</p>

# distree
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/PathoGenOmics-Lab/distree/blob/main/LICENSE)
[![distree](https://img.shields.io/badge/distree-rust-%23ff8000)](https://github.com/PathoGenOmics-Lab/distree)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/distree.svg)](https://anaconda.org/bioconda/distree)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/distree/badges/version.svg)](https://anaconda.org/bioconda/distree)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?)](https://github.com/PathoGenOmics-Lab)

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  

`distree` is a command-line tool written in Rust that extracts a distance matrix from a phylogenetic tree in Newick format. It is designed to handle large trees with thousands of sequences by using a low-memory, parallelized approach.

## Features

* **Patristic distances**: Computes the sum of branch lengths between every pair of leaves (taxa).
* **Topological distances**: Ignores branch lengths and calculates the number of edges (nodes) between leaves.
* **LMM (var-covar) matrix**: Outputs the depth (in branch-length units) of the lowest common ancestor (MRCA) for each pair of leaves.
* **Midpoint rooting**: Optionally re-roots the tree at its midpoint before computing distances.
* **Low memory footprint**: Streams each row to output immediately, without holding the full matrix in RAM.
* **Parallel computation**: Uses multiple CPU cores to compute each row in parallel, reducing runtime for large datasets.

## Typical Applications

1. **Comparative Genomics & Evolutionary Studies**

   * **Patristic Distances**: When analyzing evolutionary divergence, the sum of branch lengths between two sequences reflects their accumulated genetic change. This is useful for constructing distance-based phylogenetic trees (e.g., Neighbor-Joining), clustering sequences into operational taxonomic units (OTUs), or computing pairwise divergence metrics for downstream analyses.
   * **Topological Distances**: In cases where branch-length estimates are unreliable or not meaningful (e.g., when raw sequence counts are compared), using topological distances (number of nodes between leaves) can provide a rough measure of relatedness. This can be helpful for broad clustering or when using simpler distance-based algorithms that only require tree shape.

2. **Epidemiology & Public Health**

   * **Rapid Outbreak Tracking**: For pathogens such as bacteria or viruses, building a phylogenetic tree from whole-genome data can be computationally expensive to revisit. Extracting a distance matrix allows quick pairwise comparisons to identify clusters of closely related strains (e.g., potential transmission clusters) without recomputing distances from raw alignments.
   * **Contact Tracing & Transmission Networks**: If you have a large outbreak dataset, computing patristic distances between samples quickly helps identify subclusters or “clades” of interest (e.g., to infer likely transmission chains). Similarly, topological distances might approximate epidemiological closeness if branch-length estimates vary widely.

3. **Microbiome & Environmental Sequencing**

   * **OTU/ASV Clustering**: In 16S rRNA amplicon studies, one often constructs a phylogenetic tree of all amplicon sequence variants (ASVs). A distance matrix (patristic or topological) can feed into beta-diversity metrics (e.g., UniFrac requires branch lengths). Here, `distree` can quickly compute pairwise distances after the tree is built, facilitating ordination (PCoA) or clustering of samples by their phylogenetic composition.
   * **Phylogenetic Diversity**: Calculating the sum of branch lengths between taxa supports metrics like Faith’s Phylogenetic Diversity or UniFrac. `distree` can produce the matrix needed for those algorithms without re-traversing the tree multiple times.

4. **Machine Learning & Dimensionality Reduction**

   * **Input for MDS/t-SNE/UMAP**: Many machine learning or visualization methods (Multidimensional Scaling, t-SNE, UMAP) accept a distance matrix as input. Converting a large phylogeny into a pairwise distance matrix enables embedding taxa into low-dimensional space, highlighting evolutionary relationships or clusters.
   * **Distance-Based Feature Engineering**: In trait-based prediction models (e.g., predicting phenotype from genotype), patristic distances can be used as kernel features. `distree` can produce the kernel matrix needed for Gaussian Process Regression or support vector machines (SVM) with a custom kernel.

5. **Benchmarking & Simulation Studies**

   * **Comparing Tree-Building Methods**: When evaluating different phylogenetic inference methods, it is often necessary to compute distances from each resulting tree. `distree` can generate distance matrices for multiple trees rapidly, enabling head-to-head comparisons.
   * **Simulation Validation**: In simulation frameworks (e.g., Seq-Gen), one generates a tree, simulates sequence data, reconstructs a tree, and compares distance matrices. Having a fast, consistent tool to extract matrices accelerates simulation benchmarks.

## Installation

### 🐍 Using conda
```
conda install -c bioconda distree
```
### 🐍 Using mamba
```
mamba install -c bioconda distree
```
### Compilation
1. Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed.

2. Clone or download the `distree` repository.

3. In the project root, run:

   ```bash
   cargo build --release
   ```

4. The optimized binary will be available at:

   ```bash
   target/release/distree
   ```

## Usage

```bash
Usage: distree [OPTIONS] <phylogeny>

Arguments:
  <phylogeny>            Path to the tree file in Newick format

Options:
  --format <FORMAT>      Tree file format (only 'newick' is supported) [default: newick]
  --midpoint             Midpoint-root the tree before computing distances
  --lmm                  Produce the var-covar matrix C (depth of the MRCA)
  --topology             Ignore branch lengths; use purely topological distances
  -o, --output <FILE>     Path to write the TSV output file (defaults to stdout)
  -h, --help             Print help information
  -V, --version          Print version information
```

### Argument & Option Details

* `<phylogeny>`: Path to the input tree in Newick format. Leaf labels must be unique and not contain tabs or newline characters.

* `--format <FORMAT>`

  * Currently only `newick` is supported. Future versions may support other formats (e.g., Nexus, PhyloXML).
  * This option exists to maintain CLI consistency; it does not change parsing for now.

* `--midpoint`

  * Before computing any distances, the tree will be re-rooted at its midpoint (the point halfway along the longest path between any two leaves). Useful when no outgroup is known or when you want a balanced root for downstream analyses.
  * Use this if your downstream distance metric expects an unrooted, centrically-rooted tree.

* `--lmm`

  * “LMM” stands for var-covar matrix (matrix C) in Phylogenetic Comparative Methods. Each entry (i, j) equals the depth (distance from root) of the lowest common ancestor of leaf i and leaf j. This matrix is often used in linear mixed models (LAMM, PGLS) to account for phylogenetic covariance.
  * When specified, LMM distances override `--topology`. The output is a matrix of MRCA depths, not pairwise distances.

* `--topology`

  * Ignores branch lengths. Each entry (i, j) equals the number of edges between leaf i and leaf j:

  * Use this mode if branch-lengths are not meaningful or if you only care about tree shape.

* `-o, --output <FILE>`

  * Write the TSV distance matrix to the specified path. If omitted, the matrix is printed to standard output.
  * Example: `-o distances.tsv`.

## Output Format

The tool outputs a tab-separated values (TSV) matrix.

1. The first line is a header with leaf labels sorted alphabetically, preceded by an empty cell (for row labels).

2. Each subsequent line begins with a leaf label (sorted alphabetically), followed by N columns of distances (depending on the chosen mode) to every other leaf in the same sorted order.

Example (patristic distances):

```LeafA LeafB	LeafC
LeafA	0.000	5.200	3.100
LeafB	5.200	0.000	4.400
LeafC	3.100	4.400	0.000
```

Example (topological distances):

```LeafA LeafB	LeafC	LeafD
LeafA	0	2	3	1
LeafB	2	0	1	3
LeafC	3	1	0	4
LeafD	1	3	4	0
```

Example (LMM depths):

```LeafA LeafB	LeafC
LeafA	7.000	3.000	5.000
LeafB	3.000	7.000	2.500
LeafC	5.000	2.500	7.000
```

## Detailed Use Cases

### 1. Producing a Patristic Distance Matrix

**Scenario**: You have a large set of bacterial genomes, build a phylogenetic tree with reliable branch lengths (e.g., using RAxML or IQ-TREE), and now need the pairwise patristic distances to feed into clustering, PCoA, or hierarchical analyses.

**Command**:

```bash
./distree tree.nwk -o patristic.tsv
```

**Why**: Downstream tools like SciKit-Learn (for MDS) or R’s `ape::cmdscale()` expect a distance matrix. Patristic distances reflect evolutionary time or change.

### 2. Computing Topological Distances Only

**Scenario**: You have a phylogenetic tree where branch lengths come from different sources or are not directly comparable (e.g., concatenated multi-locus data). You prefer to measure closeness by number of shared nodes instead.

**Command**:

```bash
./distree --topology tree.nwk -o topo_distances.tsv
```

**Why**: Topological distances emphasize tree structure without scaling by substitution rate or time. Useful when comparing tree shapes or for rapid, coarse clustering when exact branch lengths are noisy.

### 3. Generating an LMM (Var-Covar) Matrix for Comparative Methods

**Scenario**: You are performing phylogenetic generalized least squares (PGLS) in R (`caper::pgls` or `nlme::gls`) and need the phylogenetic variance-covariance matrix. Each entry C\[i,j] is the distance from root to MRCA(i, j).

**Command**:

```bash
./distree --lmm tree.nwk -o varcovar.tsv
```

**Why**: In comparative models, traits shared due to common ancestry introduce covariance. This LMM matrix directly encodes that covariance structure for all taxa.

### 4. Using Midpoint Rooting Before Distance Extraction

**Scenario**: Your input tree is unrooted or rooted arbitrarily (e.g., by outgroup choice). You want a symmetric distance matrix that does not depend on outgroup, so you midpoint-root the tree first.

**Command**:

```bash
./distree --midpoint tree.nwk -o midrooted_distances.tsv
```

**Why**: Midpoint rooting places the root at a balanced position. Downstream patristic or LMM calculations become more interpretable if the root is centrally placed.

### 5. Downstream Dimensionality Reduction & Visualization

**Scenario**: You intend to visualize relationships among taxa via MDS or t-SNE. You need a distance matrix as input.

**Command**:

```bash
./distree tree.nwk | Rscript -e "dist <- as.matrix(read.table('file:///dev/stdin', header=TRUE, row.names=1)); mds <- cmdscale(dist); plot(mds)"
```

**Why**: Feeding the distance matrix directly into R for MDS (multidimensional scaling) or UMAP allows you to view clustering of taxa in two or three dimensions.

## Performance and Resource Considerations

* **Memory**: `distree` streams one row at a time. At any given moment, only a single vector of length N (number of leaves) resides in memory, plus O(M log M) for LCA structures, where M is the total number of nodes. For trees with tens of thousands of taxa, memory usage remains low.

* **Parallelism**: Each row’s distance computations are parallelized across available CPU cores via Rayon. For N taxa, computing N rows (each of size N) takes O(N^2 / #cores) time.

* **Disk I/O**: If writing to a file via `--output`, a buffered writer (`BufWriter`) minimizes I/O calls. Streaming directly to stdout also remains efficient.

## Troubleshooting and Tips

* **Invalid Newick**: Ensure your Newick tree is syntactically correct (matching parentheses, semicolon at end). `distree` will error if parsing fails.
* **Whitespace in Labels**: Leaf labels should not contain spaces, tabs, or newline characters, as these break the TSV format. Replace spaces with underscores if needed.
* **Choosing Distance Type**:

  * Use `--lmm` if performing phylogenetic comparative analyses (e.g., trait evolution, PGLS).
  * Use default patristic distances for standard evolutionary distance-based methods (e.g., Neighbor-Joining, clustering).
  * Use `--topology` for quick, coarse tree-shape comparisons when branch lengths are unreliable.
* **Large Trees**: For very large trees (e.g., >50,000 leaves), ensure you have enough CPU cores. You may run `distree` on a compute node or multi-core server to exploit parallel speedups.

---

*distree* is a versatile tool for extracting various phylogenetic distance matrices from large trees. By combining midpoint rooting, multiple distance modes, and parallel streaming, it caters to evolutionary, comparative, and clustering analyses without overwhelming memory.

---
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
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
