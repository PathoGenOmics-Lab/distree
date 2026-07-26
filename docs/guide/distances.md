# Distance modes

distree computes three different things, and they answer three different
questions. All of them come from the same primitive: for a pair of tips `i` and
`j`, find their most recent common ancestor `m`, then read off depths that were
precomputed in a single pass.

| Mode | Flag | Cell `(i, j)` |
|:--|:--|:--|
| Patristic | *(default)* | `depth(i) + depth(j) - 2 · depth(m)` |
| Topological | `--topology` | `hops(i) + hops(j) - 2 · hops(m)` |
| Variance-covariance | `--lmm` | `depth(m)` |

`depth` is the summed branch length from the root; `hops` is the number of edges
from the root. Every mode is one MRCA query plus three array reads, which is why
the cost of a cell does not grow with how far apart the tips are.

The examples below all use:

```
((LeafA:1.95,LeafB:3.25):0.35,(LeafC:0.80,LeafD:1.20):0.50);
```

## Patristic

```bash
distree tree.nwk -p 3
```

<div class="matrix" markdown>

| | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 0.000 | 5.200 | 3.600 | 4.000 |
| **LeafB** | 5.200 | 0.000 | 4.900 | 5.300 |
| **LeafC** | 3.600 | 4.900 | 0.000 | 2.000 |
| **LeafD** | 4.000 | 5.300 | 2.000 | 0.000 |

</div>

The sum of the branch lengths along the path between the two tips: how much
change has accumulated between them along the tree. This is the mode nearly
every downstream use wants, from transmission clustering to ordination.

If the tree was built with a substitution model, a patristic distance is in
substitutions per site. Multiply by the alignment length for a SNP-scale
distance:

```bash
distree tree.nwk -p 10 | awk 'NR==1 {print; next} {printf "%s", $1;
  for (i=2; i<=NF; i++) printf "\t%.1f", $i * 4411532; print ""}'
```

The diagonal is exactly `0.0` and the matrix is exactly symmetric, both by
construction rather than by rounding.

!!! note "Negative distances"

    With non-negative branch lengths, the subtraction cannot come out negative:
    `depth` accumulates from the root, so `depth(i)` and `depth(j)` are both at
    least `depth(m)`, and the arithmetic preserves that. A negative cell
    therefore means the tree has a negative branch length, which
    neighbour-joining produces routinely. distree reports it rather than
    rounding it up to zero, since a zero would claim two distinct taxa are the
    same sample.

## Topological

```bash
distree tree.nwk --topology
```

<div class="matrix" markdown>

| | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 0 | 2 | 4 | 4 |
| **LeafB** | 2 | 0 | 4 | 4 |
| **LeafC** | 4 | 4 | 0 | 2 |
| **LeafD** | 4 | 4 | 2 | 0 |

</div>

The number of edges on the path, written as integers with no decimal point and
no `-p`. LeafA and LeafB are siblings, so 2. LeafA to LeafC goes up two and back
down two, so 4.

Use it when the branch lengths are not comparable across the tree, when they
come from concatenated loci with different rates, or when only the shape of the
tree is meant to matter. A cladogram with no lengths at all is exactly this
case: patristic mode would return an all-zero matrix and warns as much.

Two things follow from counting edges rather than lengths. A polytomy shortens
distances, because the tips under it are one hop apart rather than several. And
the count does not depend on where the tree is rooted, which is why
[`--midpoint` is ignored here](rooting.md#why-midpoint-is-ignored-in-topology-mode).

## Variance-covariance

```bash
distree tree.nwk --lmm -p 3
```

<div class="matrix" markdown>

| | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 2.300 | 0.350 | 0.000 | 0.000 |
| **LeafB** | 0.350 | 3.600 | 0.000 | 0.000 |
| **LeafC** | 0.000 | 0.000 | 1.300 | 0.500 |
| **LeafD** | 0.000 | 0.000 | 0.500 | 1.700 |

</div>

Entry `(i, j)` is the distance from the root down to the MRCA of the two tips:
the length of the evolutionary history they share before they diverge. Under a
Brownian-motion model of trait evolution, that is exactly the covariance between
their trait values, which is what makes this the `C` matrix of a phylogenetic
generalised least squares fit or a phylogenetic mixed model.

Reading the matrix above:

- `C[LeafA][LeafA] = 2.300` is LeafA's own root-to-tip length, `0.35 + 1.95`.
  The diagonal is always the tip's total path from the root, because a tip's
  MRCA with itself is itself.
- `C[LeafA][LeafB] = 0.350` is the depth of the node the two share.
- `C[LeafA][LeafC] = 0.000`, because their MRCA is the root, and they share no
  history at all.

Three consequences are worth keeping in mind:

- **The matrix is not a distance.** The diagonal is not zero and larger values
  mean more similarity, not less. Do not feed it to something expecting
  distances.
- **Rooting changes it.** `depth` is measured from the root, so where the root
  sits is part of the answer. This is the only mode where `--midpoint` does
  anything.
- **An ultrametric tree gives a constant diagonal.** If every tip is the same
  distance from the root, every `C[i][i]` is that distance, which is what most
  comparative methods assume.

Patristic and variance-covariance are two views of the same tree and convert
into each other: `d(i,j) = C[i][i] + C[j][j] - 2·C[i][j]`. Check it on the
matrices above: `2.300 + 3.600 - 2 × 0.350 = 5.200`, which is what the patristic
matrix says for LeafA and LeafB.

## Choosing

| If you are | Use |
|:--|:--|
| Clustering isolates by genetic distance | Patristic |
| Running PCoA, MDS, t-SNE or UMAP on the tree | Patristic |
| Building a distance-based tree from an existing one | Patristic, with [`--lower`](output.md#phylip-lower-triangle) |
| Fitting PGLS, `caper::pgls`, `nlme::gls` or a phylogenetic mixed model | `--lmm` |
| Comparing tree shapes, or working from a cladogram | `--topology` |
| Working with branch lengths from incomparable sources | `--topology` |

## Next

- [Midpoint rooting](rooting.md), for when the root is part of the answer.
- [Output](output.md), for the formats these matrices come out in.
