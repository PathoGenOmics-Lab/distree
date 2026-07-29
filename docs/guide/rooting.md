# Midpoint rooting

```bash
distree tree.nwk --midpoint --lmm -p 3
```

`--midpoint` re-roots the tree at the middle of its longest tip-to-tip path
before any distance is computed. It exists for one mode.

## When rooting matters

A [patristic](distances.md#patristic) distance is the length of the path between
two tips. A [topological](distances.md#topological) distance is the number of
edges on that path. Neither depends on where you call the top of the tree: it is
the same path either way. Re-rooting a tree cannot change either matrix, and
distree's own tests assert exactly that.

A [variance-covariance](distances.md#variance-covariance) entry is the distance
from **the root** down to the MRCA of two tips. Move the root and every entry
moves. That makes `--lmm` the mode where rooting is a modelling decision rather
than a formality, and `--midpoint` the flag that makes it explicitly.

## What it does

The longest path between any two tips is the tree's diameter. Midpoint rooting
puts the new root exactly halfway along it, so the two deepest tips end up
equidistant from the root, and no tip is further from the root than half the
diameter.

On the tree used throughout this guide:

```
((LeafA:1.95,LeafB:3.25):0.35,(LeafC:0.80,LeafD:1.20):0.50);
```

the diameter runs LeafB to LeafD and measures `3.25 + 0.35 + 0.50 + 1.20 = 5.30`.
Half of that is 2.65, which falls **inside LeafB's own branch**, 2.65 up from
LeafB and 0.60 below the node it hangs from. distree inserts a new root there,
splitting the 3.25 branch into 2.65 and 0.60, and reverses the parent-child
links from that point up to the old root so the tree hangs correctly from its
new top.

The `--lmm` matrix before and after:

<div class="matrix" markdown>

| Original root | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 2.300 | 0.350 | 0.000 | 0.000 |
| **LeafB** | 0.350 | 3.600 | 0.000 | 0.000 |
| **LeafC** | 0.000 | 0.000 | 1.300 | 0.500 |
| **LeafD** | 0.000 | 0.000 | 0.500 | 1.700 |

</div>

<div class="matrix" markdown>

| Midpoint root | LeafA | LeafB | LeafC | LeafD |
|:--|--:|--:|--:|--:|
| **LeafA** | 2.550 | 0.000 | 0.600 | 0.600 |
| **LeafB** | 0.000 | 2.650 | 0.000 | 0.000 |
| **LeafC** | 0.600 | 0.000 | 2.250 | 1.450 |
| **LeafD** | 0.600 | 0.000 | 1.450 | 2.650 |

</div>

The whole structure has changed. LeafB now splits off at the root and shares
nothing with anyone, while LeafA, LeafC and LeafD share 0.600 of history that
the original rooting attributed differently. The deepest tips, LeafB and LeafD,
sit at 2.650, which is half the diameter, as they should.

The patristic matrix, meanwhile, is byte-identical with and without the flag.

## Why midpoint is ignored in topology mode

Adding `--midpoint` to `--topology` used to inflate the answer. The new root is
a **new node**, inserted partway along an existing edge, and it splits that edge
into two. Every pair of tips whose path crosses that edge therefore gained one
hop:

```
((A:1,B:1):1,(C:1,D:1):1);
```

| | Without `--midpoint` | With `--midpoint` |
|:--|--:|--:|
| A to C | 4 | 5 |
| A to B | 2 | 2 |

The 5 was an artefact of the extra node, not a property of the tree. Since edge
counts cannot depend on the root in the first place, there is nothing
`--midpoint` could correctly contribute here, so distree skips the rooting and
says so:

```
Warning: --midpoint is ignored in --topology mode. Edge counts do not depend on
the root, and the node inserted at the midpoint would add 1 to every distance
crossing that edge.
```

The same insertion happens in patristic mode, but there it is harmless: the two
pieces of the split edge sum to the original length, so no path changes.

## Caveats

!!! warning "Negative branch lengths"

    distree finds the diameter with the standard two-pass sweep: walk to the
    furthest tip from an arbitrary tip, then to the furthest tip from that one.
    The argument that this finds the true diameter assumes non-negative edge
    weights, and a neighbour-joining tree does not satisfy it. Midpoint rooting
    such a tree may pick the wrong path, and the run warns about the negative
    lengths for that reason among others.

**A zero-length diameter leaves the tree alone.** If every branch is zero there
is no midpoint to find, and distree returns the tree as it was.

**Only tips count.** The diameter is measured tip to tip, matching what `ape`
and `phytools` do. A root with a single child is a degree-one node in the
unrooted tree, but it is not a tip and does not extend the diameter.

**The new root is unlabelled** and adds one node to the tree. That is invisible
in the output, which only ever has one row per labelled tip.

## Next

- [Distance modes](distances.md), for what rooting does and does not affect.
- [How it works](../how-it-works/algorithm.md#midpoint-rooting), for the
  re-rooting procedure itself.
