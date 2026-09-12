# Sample without wasting rays

A chopper train that passes a percent of the plane makes a source spend ninety-nine percent
of its rays on neutrons the discs stop. Drawing from the transmitted region instead keeps
all of them — and is the *same measurement*, provided every weight is corrected.

## Draw

```python
--8<-- "sample_the_region.py:sampler"
```

```
3 triangles, acceptance 3.929119e-03
(200000,) (200000,)
```

Each polygon is fanned into triangles from its first vertex, so `count` is triangles rather
than polygons. A draw costs one binary search and three uniform deviates, and **never
rejects**.

`sample(n, generator)` needs numpy. The compiled `draw` and `draw_many` take their deviates
as arguments, so that chopper-lib needs no generator of its own and a McStas TRACE can hand
over `rand01()`:

```python
inverse_velocity, time = sampler.draw(0.5, 0.25, 0.25)
```

## Correct the weight

**Multiply every drawn ray's weight by `sampler.acceptance`.** Not only the ones that
needed redrawing — the ones that would have landed inside anyway are a draw from the same
restricted distribution and carry the same correction.

The argument: a ray drawn from proposal `q` carrying weight `w` estimates a downstream
tally as `N·E_q[w·f]`, and `f` is zero outside the allowed set `A` because the discs stop
those rays. Restricting `q` to `A` multiplies that by `1/Q`, where `Q` is the probability a
draw from `q` lands in `A` — so multiplying by `Q` puts it back.

Here `Q` is a ratio of two areas, both known in closed form, so it is **exact**:

```python
sampler.acceptance == region.area / source.area
```

## Check it yourself

Two estimators of the same quantity — uniform-with-rejection, and draw-from-the-region with
the correction:

```python
--8<-- "sample_the_region.py:unbiased"
```

```
reject: 1.462563e-03   direct: 1.461009e-03
```

They agree, and the direct one has a far smaller error for the same number of draws.

!!! danger "Two things `acceptance` is not"

    **Not the transmission.** That is the allowed fraction of a *weighted* signal —
    `unmasked_probability` — which is what an instrument sees and the wrong normalisation
    for this. `acceptance` counts draws, not intensity.

    **Not recoverable from a rejection loop.** For a geometric number of trials `k`,
    `E[1/k]` is not `1/E[k]`, so counting attempts per ray biases the answer.

## When the correction does not apply

`acceptance` is the right factor only while **both coordinates are drawn uniformly and
independently** — of each other, and of everything else the caller samples.

A source that picks its emission time from a window centred on the neutron's own velocity
breaks that. McStas' `ESS_butterfly` under time focusing does exactly that, and no single
number corrects it.

That is not a loss: time focusing already *is* this trick, done exactly. Its window traces
a band of slope `−tfocus_dist` in the plane — the shape of a chopper band — and
`w_tfocus` is already the compensating weight.

## Against the grid

`MaskSampler` does the same job over a grid, and its `acceptance` can only be too **large**,
because a partly covered cell is weighted whole. See [Work on a grid
instead](grid-and-mask.md) for the size of that.
