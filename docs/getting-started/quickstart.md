# Quick start

Every example on this page uses the same four-tip tree, so the matrices can be
compared directly:

```
((LeafA:1.95,LeafB:3.25):0.35,(LeafC:0.80,LeafD:1.20):0.50);
```

## The default: patristic distances

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

Each cell is the sum of the branch lengths on the path between the two tips.
LeafA to LeafB is `1.95 + 3.25 = 5.20`; LeafA to LeafC crosses the root,
`1.95 + 0.35 + 0.50 + 0.80 = 3.60`.

`-p` sets the decimal places. The default is 10, which is more than most branch
lengths carry any information at; 6 is a reasonable working value and 3 is
enough to read by eye.

## Reading from stdin

`-` reads the tree from standard input, so distree drops into a pipeline:

```bash
gunzip -c tree.nwk.gz | distree - -o distances.tsv
```

## Ignoring branch lengths

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

The number of edges between the two tips, as integers. Useful when the branch
lengths come from different sources, or when only the shape of the tree is
meant to matter.

## The variance-covariance matrix

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

Entry `(i, j)` is the distance from the root down to the most recent common
ancestor of `i` and `j`: how much evolutionary history the two tips share. The
diagonal is each tip's own root-to-tip length, and any pair whose ancestor is
the root scores 0. This is the `C` matrix that PGLS and phylogenetic mixed
models expect, and it is the one mode where [where you root the
tree](../guide/rooting.md) changes the answer.

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

A taxa count, then one row per taxon holding its label and the distances to
every taxon above it. Half the cells and no diagonal, which is what PHYLIP,
Mash and most distance-matrix readers want.

## Midpoint rooting

```bash
distree tree.nwk --midpoint --lmm -p 3
```

<div class="matrix" markdown>

| | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 2.550 | 0.000 | 0.600 | 0.600 |
| **LeafB** | 0.000 | 2.650 | 0.000 | 0.000 |
| **LeafC** | 0.600 | 0.000 | 2.250 | 1.450 |
| **LeafD** | 0.600 | 0.000 | 1.450 | 2.650 |

</div>

The longest path in this tree runs from LeafB to LeafD and is 5.30 long, so the
new root goes 2.65 from each: partway along LeafB's own branch. Every tip is now
at most 2.65 from the root, and the shared history of the other three tips has
been reorganised around it.

!!! tip "It cannot change a patristic matrix"

    Rooting moves the root, not the tips. The path between two tips is the same
    path whatever you call the top of the tree, so `--midpoint` leaves a
    patristic or topological matrix exactly as it was. Reach for it when you are
    producing `--lmm`. See [Midpoint rooting](../guide/rooting.md).

## Controlling threads and output

```bash
distree tree.nwk -t 8 -o distances.tsv
```

`-t` caps the thread count; without it distree uses every core. `-o` writes to a
file instead of stdout, and is not created until the tree has parsed, so a
failed run leaves whatever was there before untouched.

## Next

- [Input](../guide/input.md), for what distree accepts and what it turns down.
- [Distance modes](../guide/distances.md), for the three modes in full.
- [Recipes](../recipes.md), for feeding the matrix to R, Python or PHYLIP.
