# Three answers to one question

chopcal has three ways to answer "what does this chopper train pass?", and they do not
agree. This is why, and which to believe.

## The question

A neutron leaves the source at time `t` with inverse velocity `a`. It reaches a disc `L`
metres downstream at `t + L·a`, and gets through if that lands in an opening. The set of
`(a, t)` that survive every disc is the **transmitted region**, and every question anyone
asks of a chopper train is a question about it — its extent in `a`, its area, a point's
membership, a draw from it.

## Answer 1: the window functions

`wavelength_windows`, `wavelength_limits`, and their inverse-velocity counterparts.

For each disc separately, work out which inverse velocities could get through *for some*
emission time in the pulse. Intersect those sets across discs.

That is fast, and it is an **over**-estimate. Each disc is allowed its own emission time,
so the answer can contain an inverse velocity that no *single* emission time delivers.
Formally it is an intersection of projections, which contains the projection of the
intersection. [The next page](over-reporting.md) works through the BIFROST case, where it
reports a 0.03 Å band near 38 Å that nothing passes.

It is tight wherever the discs leave wide, overlapping emission windows — which is the
usual case, and includes every band an instrument is deliberately set up to pass.

## Answer 2: the mask

`inverse_velocity_time_mask` keeps both coordinates, so it cannot make that error. It
samples the region onto a grid instead, and a grid is wrong in two directions at once:

* a channel thinner than a bin is **lost**;
* a bin the region only partly covers is kept **whole**.

For counting purposes the second dominates, so the mask's area — and `MaskSampler`'s
`acceptance`, which multiplies every drawn ray's weight — comes out too large. The error
falls as the cell count, so halving it costs four times the memory, and it is never zero.

## Answer 3: the region

`Region` builds the transmitted set itself.

A disc open on `[lower, upper]` accepts `lower ≤ t + L·a ≤ upper`: a slab between two
parallel lines. A train's acceptance is an intersection of unions of slabs, one union per
disc. Intersection distributes over union, so

> the transmitted region **is** a union of convex pieces — one per choice of which opening
> and which turn of each disc a neutron goes through.

Nothing is approximated. Two consequences make it cheap: every piece is an intersection of
half-planes and therefore convex, so no general polygon intersection is needed; and one
half-plane clip adds at most one vertex, so the vertex count is bounded in advance. The six
BIFROST discs leave **one polygon of five vertices**.

## Which to believe

Where the mask and the windows disagree, the mask is right. Where the region and the mask
disagree, the region is right.

| you want | use |
|---|---|
| a quick check that a train passes roughly the band you meant | `wavelength_windows` |
| an array to plot, or a mask from a McStas run | `inverse_velocity_time_mask` |
| the area, the acceptance, a membership test, or an unbiased sampler | `Region` |

## Why keep the other two

The window functions are the cheapest thing to call and need no source rectangle — just a
train. The mask is the only one that gives you pixels, and it is the format McStas
components write.

Neither is being deprecated. They are approximations whose direction and size are known,
which is a different thing from being wrong.
