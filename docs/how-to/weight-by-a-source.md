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

The region's boundary is where transmission reaches zero for the beam width the discs'
`aperture` describes. With `apertures=False`, the point beam, the boundary is instead
where a linear ramp across the beam is half open. If the source changes little over one
ramp, integrating over that point-beam region gives the transmitted total correct to first
order. The ramps mostly soften the edges of the arrival-time profile rather than change
how much gets through.
