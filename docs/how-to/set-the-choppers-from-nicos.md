# Set the choppers from NICOS

An ESS chopper controller is not set with chopcal's `delay`. It locks its
top-dead-centre sensor to a total that the EPICS module `evr_seq_calc` forms from
three parts, in ns:

    TotDly = MechDly / 360 / |Spd_S| * 1e9 + BeamPosDly + ChopDly

| PV | NICOS | what it is |
|---|---|---|
| `MechDly-S` | — | measured for each disc, from its TDC to the centre of its reference opening, so that with the other two at zero that centre is on the beam at the timing reference |
| `BeamPosDly-S` | — | the pulse delay: moves every disc together; fixed by the instrument team |
| `ChopDly-S` | `<disc>_chopper_delay` | the experiment delay: what an experiment sets |

So the reference opening is centred on the beam `pulse_delay + experiment_delay` after
the reference, as ECDC's Choppers page defines the delays, and how the real disc is
mounted never enters what an experiment sets.
`chopcal.nicos` turns chopcal's settings into those numbers:

```pycon
>>> from chopcal.nicos import bifrost
>>> bifrost(pulse_delay=10_000, energy_min=7.0)
{'psc1_chopper_speed': -196.0, 'psc1_chopper_delay': 4539047.1...,
 'psc2_chopper_speed': -196.0, 'psc2_chopper_delay': 1646303.3...,
 'foc1_chopper_speed': -14.0, 'foc1_chopper_delay': 7032394.5...,
 'foc2_chopper_speed': -14.0, 'foc2_chopper_delay': 11121693.7...,
 'bwc1_chopper_speed': -14.0, 'bwc1_chopper_delay': 51124217.6...,
 'bwc2_chopper_speed': 14.0, 'bwc2_chopper_delay': 51124217.6...}
```

Each key is the NICOS device that takes the value: a speed in Hz, and an experiment
delay in ns, folded into one turn of its disc. The speed is signed as the controller's
`Spd_S` is, positive clockwise facing away from the source. That is the opposite of
chopcal's and NeXus's sign, so a disc chopcal turns anticlockwise gets a negative speed
(`bwc2` counter-rotates).
`pulse_delay` is the `BeamPosDly-S` the controllers hold, in ns. It has no default,
because no value is right for every installation.

The pulse shaping pair turns at the fastest harmonic of 14 Hz that both the requested
`shaping_time` and the installed disks allow. The disks BIFROST was designed for are
rated to 196 Hz; a substitute set rated lower needs its own limit, which NICOS knows as
the largest speed its `psc1_chopper_speed` device maps:

```pycon
>>> bifrost(pulse_delay=10_000, energy_min=7.0, pulse_shaping_speed_max=126)['psc1_chopper_speed']
-126.0
```

`chopcal.nicos.experiment_delay(chopper, pulse_delay)` does the same for one
`Chopper`, in seconds. chopcal's `delay` is when the centre of a disc's one opening is
on the beam, which is the time the controller wants, so the experiment delay is that
less the pulse delay, whichever way the disc turns. A disc of several openings is
refused: which one is the reference opening depends on where its top-dead-centre mark
is, and chopcal does not know.

## What this assumes

- chopcal's times are measured from the proton pulse, and that is the reference the
  choppers and the detectors are timed from (NiSync rather than event14).
- `MechDly-S` is measured to the centre of the reference opening, as ECDC's Choppers
  page defines it.
- chopcal's direction for each disc is the real one. With BIFROST's NICOS setups mapping
  only positive `Spd_S` for the pulse-shaping and frame-overlap discs, those discs would
  turn clockwise facing away from the source, the opposite of chopcal's.
- Each real disc has the one opening chopcal describes, of the same width. NICOS's
  BIFROST setups say so: 170, 170, 38.26, 52.01, 161 and 161 degrees.
