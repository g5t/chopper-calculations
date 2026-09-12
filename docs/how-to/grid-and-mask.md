# Work on a grid instead

Sometimes you want an array to plot rather than a region to reason about. The mask is that
array — and it is the one part of chopcal that needs **numpy**.

## Build a mask

```python
--8<-- "grid_and_mask.py:mask"
```

```
(1500, 2000) time by inverse velocity
9814 of 3000000 bins allowed
```

!!! warning "The shape is (time, inverse velocity)"

    Time is the slow axis — the reverse of the argument order. That is how the C stores it,
    and it is easy to get backwards.

`inverse_velocities` and `times` are bin **edges**, so the mask is one smaller in each
direction. Pass lists if that is what you have; they are converted.

`grow=n` expands each allowed region by `n` bins in every direction. It is a blunt
instrument — it admits inverse velocities no disc ever passes — so prefer
[`Chopper.aperture`](describe-a-chopper.md#the-aperture), which opens the windows in time
only, where the beam's width actually acts.

## Transmission through it

```python
--8<-- "grid_and_mask.py:probability"
```

`unmasked_probability` is the allowed fraction of a **weighted** signal — the transmission
an instrument sees. It is *not* the weight correction a sampler needs; see
[Sample without wasting rays](sample-the-region.md).

## What the grid costs

```python
--8<-- "grid_and_mask.py:versus"
```

```
grid    4.905000e-03 from 9814 cells
region  3.929119e-03 exactly
grid is 24.84% high
```

The error is one-signed — a partly covered cell is counted whole, so a grid can only
**over**-accept — and it falls only as fast as the cell count, so halving it costs four
times the memory. Three million cells is still a quarter high here, on a deliberately wide
source; on a tight one it does much better, but it is never exact.

That matters most for `MaskSampler.acceptance`, which multiplies every drawn ray's weight.
A grid figure that is 25% high makes every weight 25% wrong.

## When the grid is still the right tool

* You want a picture. `Region` has no pixels.
* You are reading a mask written by a McStas component.
* The over-estimate is harmless and you would rather not carry a polygon set around.

Otherwise reach for [`Region`](transmitted-region.md).
