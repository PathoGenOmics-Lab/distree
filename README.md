# distree

`distree` is a command-line tool written in Rust that extracts a distance matrix from a phylogenetic tree in Newick format. It is designed to handle large trees with thousands of sequences by using a low-memory, parallelized approach.

## Features

* **Patristic distances**: Computes the sum of branch lengths between every pair of leaves.
* **Topological distances**: Ignores branch lengths and calculates the number of edges between leaves.
* **LMM (var-covar) matrix**: Outputs the depth of the lowest common ancestor (in branch-length units) for each pair of leaves.
* **Midpoint rooting**: Optionally re-roots the tree at its midpoint before computing distances.
* **Low memory footprint**: Does not store the entire matrix in memory; it outputs each row as soon as it is computed.
* **Parallel computation**: Uses multiple CPU cores to compute each row in parallel, speeding up large matrices.

## Installation

1. Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed.

2. Clone or download the `distree` repository.

3. In the project root, run:

   ```bash
   cargo build --release
   ```

4. The optimized binary will be available at:

   ```
   target/release/distree
   ```

## Usage

```bash
distree [OPTIONS] <phylogeny>
```

### Required Argument

* `<phylogeny>`: Path to the input tree file in Newick format.

### Optional Flags and Options

* `--format <FORMAT>`

  * Description: Specify the tree file format.
  * Default: `newick`
  * Example: `--format newick`

* `--midpoint`

  * Description: Midpoint-root the tree before computing distances.
  * Usage: Include this flag if you want to root the tree at its midpoint.
  * Example: `--midpoint`

* `--lmm`

  * Description: Produce the var-covar matrix C. For each pair of leaves, outputs the depth (in branch-length units) of their lowest common ancestor.
  * Usage: Include this flag to generate LMM distances instead of patristic distances.
  * Example: `--lmm`

* `--topology`

  * Description: Ignore branch lengths and compute purely topological distances (number of edges between leaves).
  * Usage: Include this flag if you want edge-count distances instead of branch-length sums.
  * Example: `--topology`

* `--output <FILE>` or `-o <FILE>`

  * Description: Write the output matrix to the specified file instead of standard output.
  * Usage: Provide a path to the desired output file.
  * Example: `-o distances.tsv`

## Output

The tool outputs a tab-separated values (TSV) matrix where both rows and columns represent leaf labels sorted alphabetically. The first line is a header with leaf labels. Each subsequent line begins with a leaf label, followed by its distance to every other leaf in the same sorted order.

Example (patristic distances):

```LeafA LeafB	LeafC
LeafA	0.0	5.2	3.1
LeafB	5.2	0.0	4.4
LeafC	3.1	4.4	0.0
```

## Examples

1. **Compute patristic distances (default)**

   ```bash
   distree tree.nwk > distances.tsv
   ```

2. **Compute topological distances**

   ```bash
   distree --topology tree.nwk > topo_distances.tsv
   ```

3. **Compute LMM (var-covar) distances**

   ```bash
   distree --lmm --midpoint tree.nwk -o lmm_matrix.tsv
   ```

4. **Use midpoint rooting before computing distances**

   ```bash
   distree --midpoint tree.nwk > mid_distances.tsv
   ```

## Notes

* If both `--lmm` and `--topology` are specified, `--lmm` takes precedence (the tool will compute LMM distances).
* If neither `--lmm` nor `--topology` is specified, the default is patristic distances.
* Ensure the input Newick file is valid and that leaf labels do not contain whitespace or special characters that could disrupt TSV formatting.

---

*distree* is optimized for large phylogenies. By computing each row in parallel and streaming results, it minimizes RAM usage while leveraging multi-core CPUs for faster performance.
