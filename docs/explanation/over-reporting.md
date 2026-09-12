# Why the window functions over-report

`wavelength_windows` on the BIFROST train reports two bands:

```
1.1454 to 3.0515 AA
37.9630 to 37.9950 AA
```

The second is not a band the instrument delivers. This is why it appears.

## The relaxation

`chopper_inverse_velocity_windows` builds, for each disc, the set of inverse velocities
that could pass it — and it lets the emission time range over the whole pulse while doing
so. In the C:

```c
double wiv_min = (t_open[opening] + n_tau - latest_emission) / path;
double wiv_max = (t_close[opening] + n_tau) / path;
```

`wiv_min` uses emission at the *end* of the pulse and `wiv_max` emission at the *start*.
So each disc's set is its slice of the transmitted region **projected onto the inverse
velocity axis**, with the emission time integrated away.

Then those sets are intersected across discs.

## Why that is not the same question

Write `A_i` for the region disc `i` passes, and `π` for projection onto the inverse
velocity axis. What the code computes is

```
⋂ᵢ π(Aᵢ)
```

What was asked for is

```
π( ⋂ᵢ Aᵢ )
```

and these are not equal. An inverse velocity survives the first if *each disc* admits it
at *some* emission time — the discs being allowed different times. It survives the second
only if *one* emission time gets it through *all* of them.

Intersection of projections always contains the projection of the intersection, so the
computed answer is a superset: it can name bands that nothing passes, and never miss one
that something does. The error is one-signed.

## Where it bites

The relaxation is tight when the discs leave wide, overlapping emission windows: then "some
time for each" and "one time for all" coincide, and the two expressions agree.

It fails when the windows are narrow and disjoint. At 38 Å on BIFROST, with the beam at its
real width:

| disc | emission time it permits |
|---|---|
| ps1 | 387 µs |
| ps2 | 2376 µs |
| **fo1** | **32 µs** |
| **fo2** | **64 µs** |
| bw1, bw2 | the whole 3 ms |

Every disc permits *something*, so every disc's projection contains 38 Å and the
intersection keeps it. But the frame-overlap pair's slivers **do not overlap each other**,
so no single emission time gets through both. The intersection of the regions is empty.

## Three checks that it is empty

Independent of each other and of the polygon code:

1. **A brute-force scan** over 200 001 emission times across the pulse: every one blocked.
2. **An exact interval intersection** — each disc's permitted emission times as intervals,
   intersected: empty.
3. **`inverse_velocity_time_mask`**, which keeps both coordinates, over a fine grid
   spanning 37.90–38.05 Å: **0 allowed cells** out of six million.

And `Region`, which projects the intersection rather than intersecting projections, reports
one band.

## It is not only the spurious band

The same relaxation widens the *real* band slightly. The window function's low edge is
1.145358 Å; the region's is 1.145836 Å — 0.48 mÅ of inverse velocities that no emission
time delivers.

## What to do about it

Nothing, usually. If you are asking "roughly what band does this train pass", the answer is
right to a fraction of a milliångström, and the spurious bands appear at wavelengths no
spectrometer is looking at. Take `wavelength_windows(...)[0]` and get on with it.

If a number is going into a calculation — an acceptance, an area, a source's sampling
bounds — use [`Region`](../how-to/transmitted-region.md), which never separates the two
coordinates and so cannot make this error.
