# Weight the region by a source

`Region` says which neutrons a train passes, not how many: whether a neutron gets
through depends only on where it is in (inverse velocity, emission time), not on how
intense the source is there. So the transmitted intensity is the source intensity
integrated over the region, and the region only has to be built once.

## A gridded source

Give the source as a histogram: the amount in each (time, inverse velocity) cell, spread
evenly over it — what a Monte Carlo source or a tabulated moderator spectrum gives you.
Shape it as the masks are, **(time, inverse velocity)**, with bin edges in s and s/m.

```python
--8<-- "weight_by_source.py:grid"
```

```
(60, 400) 1.000
transmitted 28.1285 of 7600.9
```

`coverage` is the fraction of each cell the region covers, found by clipping each polygon
to each cell it reaches. `integrate` is the signal times the coverage, summed. Both are
exact for a histogram, where `inverse_velocity_time_mask` and `unmasked_probability` count
a partly covered cell whole.

## At a monitor or the sample

A neutron emitted at inverse velocity `a` and time `t` reaches distance `L` at
`t + L*a`. That shear keeps areas, so `arrivals` carries each cell's share of the region
into arrival-time bins exactly: the density times the area between two lines of constant
arrival time.

```python
--8<-- "weight_by_source.py:arrivals"
```

```
(714, 400) peak at 51.8 ms
```

The result is **(arrival time, inverse velocity)**, because the shear leaves inverse
velocity alone: sum the last axis for the time profile, or keep it for the wavelength as
well. `period` folds arrival times into one frame, as a detector's `event_time_offset`
does; without it, give arrival bins that reach the slowest neutron. Distance 0 gives the
emission times back.

Every neutron is given one flight path. A guide's spread of paths smears arrival times,
and is not included here.

## A finite beam

With disc `aperture`s, which `chopcal.bifrost` sets by default, a `Region`'s boundary is
where transmission falls to zero: everything inside counts as passing, including neutrons
that only the edge of the beam carries through. Inside, an edge sweeping a uniformly lit
beam ramps the transmission linearly, and `transmit_weighted` carries that ramp:

```python
--8<-- "weight_by_source.py:weighted"
```

```
support 1.3289e-07, point beam 9.2192e-08, weighted 8.9711e-08 s^2/m
transmitted 17.1654
```

The result is a `WeightedRegion`, with the same `coverage`, `integrate` and `arrivals`;
each piece counts by the integral of its transmission rather than its area. Here the
source passes 17.2 rather than the 28.1 the support counts.

The choppers handed to `transmit_weighted` carry `beam_aperture(..., at='centre')`. The
region's edge is where the last of the beam is cleared, at the window's inner corner,
which is what `chopcal.bifrost` and the default `at='inner'` give. The ramp in between is
set by the beam centre: while an edge spans the window's height, the fraction it has
passed grows linearly at the rate the centre radius sets. Against McStas with the real
window geometry, the centre angle fits a BIFROST train's arrival times to within counting
error; the inner one leaves the tails of the profile 6-18% low.

`groups` says which discs share a housing. A neutron crosses those at the same place
across the beam, so a group passes the part of the beam open in all of them. For a pair
turning the same way that is the smaller of the two transmissions. For a counter-rotating
pair, like BIFROST's bandwidth discs, equal openings sweeping from opposite sides stay shut
until the edges meet at the half-way time and then open twice as fast. Different groups
multiply, as reflections in a guide scramble where across the beam a neutron is. Without
`groups` every disc is its own group.

The point beam, `apertures=False`, has its boundary where each ramp is half open. Its area
is the first-order answer, within 2.8 % of the weighted one here; the difference is where
ramps of different discs overlap.
