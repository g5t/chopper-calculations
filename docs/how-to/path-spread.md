# Allow for a guide's path spread

A neutron in a guide travels further than the straight line between two points, and how
much further depends on where it bounced. That smears arrival times, and so widens what a
train passes.

## Apply it

`transmit` takes one extra flight path per chopper, in metres:

```python
--8<-- "path_spread.py:spread"
```

```
straight: [(1.1458362867034736, 3.0515432206219723)]
spread:   [(1.1446919..., 3.0515432206219723)]
area grew by 1.185%
```

A tenth of a percent of extra path is a generous guide; `1e-4` is a realistic figure. The
deviation accumulates along the flight path, so scaling each disc's own `path` is the
natural way to express it.

## It is not an even growth

The deviation is in **path**, so what it does to an arrival time is
`deviation × inverse_velocity` — larger for a slow neutron, and nothing at all to the
inverse velocity itself.

So it does not dilate the region evenly in every direction. It tilts one of the two lines
bounding each slab, opening it into a wedge that widens as the inverse velocity grows.

Only the slow edge moves:

```python
--8<-- "path_spread.py:one_sided"
```

```
slow edge moved 1.1443 mAA
fast edge moved 0.0000 mAA
```

That is right, and a good sign the model is doing what it should: a guide can only **add**
path, so a longer route rescues a neutron that would have arrived too early and can do
nothing for one arriving too late.

If you would rather have a symmetric `[L − δ, L + δ]`, that is a chopper-lib call away —
the wedge takes both paths independently — but `transmit` exposes only the one-sided form,
because that is what a guide does.

## What it gives you, and what it does not

The wedge is the **support** of the transmitted set under path uncertainty, not a weighting
over it: a neutron is passed if *some* path in the range would have got it through, with no
account of how likely that path is.

So it re-introduces an over-estimate — but a principled one, bounded by a number you chose
from guide physics, rather than an artefact of an algorithm. It is the same bargain
[`Chopper.aperture`](describe-a-chopper.md#the-aperture) makes for the width of the beam.

A real path-length distribution is not a hard interval, and polygons cannot carry a weight
across their interior. If you need the distribution rather than its support, use the region
as an importance-sampling envelope and apply the weight per ray afterwards.

## The limit

A spread wide enough to reach from one turn of a disc into the next would stop the
transmitted pieces being disjoint, and their areas would count twice. That is **refused**
rather than computed:

```python
source.transmit(train, [0.5 * chopper.path for chopper in train])
# RuntimeError: chopper-lib refused to transmit this region ...
```

The condition is `spread × inverse_velocity_maximum < period − opening`, checked per disc
against the slowest neutron in the region. For BIFROST against a 45 Å ceiling the tightest
discs break at `dL/L ≈ 0.036`; a real guide is two orders of magnitude inside that.
