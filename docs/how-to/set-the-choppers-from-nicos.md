# Set the choppers from NICOS

An ESS chopper controller is not set with chopcal's `delay`. It locks its
top-dead-centre sensor to a total that the EPICS module `evr_seq_calc` forms from
three parts, in ns:

    TotDly = MechDly / 360 / |Spd_S| * 1e9 + BeamPosDly + ChopDly

| PV | NICOS | what it is |
|---|---|---|
| `MechDly-S` | — | measured for each disc, so that with the other two at zero its first opening edge, NeXus's `slit_edges[0]`, is at the beam at the timing reference |
| `BeamPosDly-S` | — | the pulse delay: moves every disc together; fixed by the instrument team |
| `ChopDly-S` | `<disc>_chopper_delay` | the experiment delay: what an experiment sets |

So the first opening edge is at the beam `pulse_delay + experiment_delay` after the
reference, and how the real disc is mounted never enters what an experiment sets.
`chopcal.nicos` turns chopcal's settings into those numbers:

```pycon
>>> from chopcal.nicos import bifrost
>>> bifrost(pulse_delay=10_000, energy_min=7.0)
{'psc1_chopper_speed': 196.0, 'psc1_chopper_delay': 641654.8...,
 'psc2_chopper_speed': 196.0, 'psc2_chopper_delay': 2850951.9...,
 'foc1_chopper_speed': 14.0, 'foc1_chopper_delay': 10828029.4...,
 'foc2_chopper_speed': 14.0, 'foc2_chopper_delay': 16281415.9...,
 'bwc1_chopper_speed': 14.0, 'bwc1_chopper_delay': 67096439.9...,
 'bwc2_chopper_speed': -14.0, 'bwc2_chopper_delay': 35151995.4...}
```

Each key is the NICOS device that takes the value: a speed in Hz, signed (`bwc2`
counter-rotates), and an experiment delay in ns, folded into one turn of its disc.
`pulse_delay` is the `BeamPosDly-S` the controllers hold, in ns. It has no default,
because no value is right for every installation.

The pulse shaping pair turns at the fastest harmonic of 14 Hz that both the requested
`shaping_time` and the installed disks allow. The disks BIFROST was designed for are
rated to 196 Hz; a substitute set rated lower needs its own limit, which NICOS knows as
the largest speed its `psc1_chopper_speed` device maps:

```pycon
>>> bifrost(pulse_delay=10_000, energy_min=7.0, pulse_shaping_speed_max=126)['psc1_chopper_speed']
126.0
```

`chopcal.nicos.experiment_delay(chopper, pulse_delay)` does the same for one
`Chopper`, in seconds. chopcal's `delay` is when the centre of a disc's one opening is
on the beam, so the edge `slit_edges[0]` reaches the beam `width / 2` of a turn later
for a disc turning anticlockwise (it closes the slit) and earlier for one turning
clockwise (it opens it). A disc of several openings is refused: which one comes first
depends on where its top-dead-centre mark is, and chopcal does not know.

## What this assumes

- chopcal's times are measured from the proton pulse, and that is the reference the
  choppers and the detectors are timed from (NiSync rather than event14).
- `MechDly-S` is measured to `slit_edges[0]`, the geometric first opening edge, rather
  than to whichever edge reaches the beam first.
- Each real disc has the one opening chopcal describes, of the same width. NICOS's
  BIFROST setups say so: 170, 170, 38.26, 52.01, 161 and 161 degrees.
