# Place the BIFROST band

`chopcal.bifrost` sets all six discs from one number.

```python
--8<-- "quickstart.py:settings"
```

BIFROST passes a fixed band about **1.91 Å** wide, so one number decides where it sits.
Give it either way round:

```python
chopcal.bifrost(wavelength_max=3.0)   # longest wavelength to pass, angstrom
chopcal.bifrost(energy_min=7.0)       # lowest energy to pass, meV
```

A positive `wavelength_max` wins if you give both. Giving neither raises, rather than
returning settings with infinite delays.

!!! warning "The number you give is the *slow* edge"

    `energy_min=7` does not cap the incident energy at 7 meV. It puts 7 meV at the
    **bottom** of the band and lets everything faster through, up to about 31 meV. Both
    names this argument has carried said the opposite, which is why the tests pin the
    direction rather than trusting the name.

## What comes back

A `ChopperSet`: a `dict` keyed by disc name, in beam order, that prints as a table.

| name | what it does |
|---|---|
| `ps1`, `ps2` | the pulse shaping pair, co-rotating at a harmonic of 14 Hz |
| `fo1`, `fo2` | frame overlap suppression |
| `bw1`, `bw2` | the bandwidth pair, counter-rotating |

Everything a `dict` does, it does. `settings['bw1']`, `settings.values()`, `**settings`.

## Change the burst length

`shaping_time` is how long the pulse shaping pair bursts for, in seconds — 0.2 ms by
default. The two discs co-rotate with one lagging the other, so the burst is one slit
crossing less the lag.

```python
chopcal.bifrost(wavelength_max=3.0, shaping_time=1e-3)
```

Ask for longer than a slit crossing at 196 Hz and the pair drops to a lower harmonic, and
says so on stdout.

## Turn the beam width off

By default the discs carry their real `aperture`, from BIFROST's own disc geometry. Pass
`apertures=False` for the pencil beam every release before 0.6 described:

```python
chopcal.bifrost(wavelength_max=3.0, apertures=False)
```

That narrows the band from `BIFROST_BANDWIDTH` (1.9061 Å) to
`BIFROST_BANDWIDTH_POINT_BEAM` (1.7715 Å). Only the `aperture` field changes — the speeds,
delays, edges and paths are identical, because a beam width does not move a chopper.

## Read the band back

```python
--8<-- "quickstart.py:band"
```

The band's top edge sits about 0.05 Å past what you asked for. That is not an error: the
settings place the band's *centre*, and a real beam widens it about that centre. With
`apertures=False` it falls about 0.015 Å short instead.

## Next

* [Describe a chopper](describe-a-chopper.md) — the fields on each disc.
* [Find the band a train passes](find-the-band.md) — for any discs, not only BIFROST's.
