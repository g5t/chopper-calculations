# Find the band a train passes

Given any set of `Chopper`s, which wavelengths get through?

## The quick answer

```python
--8<-- "find_the_band.py:windows"
```

```
1.1454 to 3.0515 AA
37.9630 to 37.9950 AA
```

`wavelength_windows` gives the bands; `wavelength_limits` gives the envelope around them
as `(count, (low, high))`:

```python
--8<-- "find_the_band.py:limits"
```

```
2 band(s), envelope 1.1454 to 37.9950 AA
```

Reach for the window list whenever there is more than one band. An envelope spanning two
bands also spans the **gap** between them, which the train does not pass — here it claims
everything from 1.1 Å to 38 Å.

Both take an optional search range and emission window:

```python
wavelength_windows(train, wavelength_min=0.5, wavelength_max=10.0, latest_emission=3e-3)
```

The defaults search 1e-4 Å to 100 Å with a 3 ms pulse.

## The catch

That second band is not real. **Nothing gets through at 38 Å** — three independent checks
agree, including chopper-lib's own mask.

`chopper_inverse_velocity_windows` works out each disc's admissible inverse velocities
letting the emission time range over the whole pulse *independently per disc*, then
intersects those ranges. An intersection of projections is a superset of the projection of
the intersection, so it can report a band that no single emission time delivers.

It is tight wherever the discs leave wide, overlapping emission windows — which is the
usual case, and includes every band an instrument is set up to pass. BIFROST with its real
beam is the exception. [Why the window functions over-report](../explanation/over-reporting.md)
has the derivation.

## The exact answer

```python
--8<-- "find_the_band.py:region"
```

```
[(1.1458362867034736, 3.0515432206219723)]
```

One band, and the low edge 0.5 mÅ above what the quick answer claimed. `Region` builds the
region itself rather than projecting it, and cannot make that error.

## Which to use

| you want | use |
|---|---|
| a quick check that a train passes roughly the band you meant | `wavelength_windows` |
| a number going into a calculation | [`Region`](transmitted-region.md) |
| an array to plot | [`inverse_velocity_time_mask`](grid-and-mask.md) |
