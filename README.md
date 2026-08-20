# Chopcal

## Description
Exposes chopper calculations to Python which were otherwise hidden in McStas instruments.


## Supported components

| source          | name                                                            | component                              | type      |
|-----------------|-----------------------------------------------------------------|----------------------------------------|-----------|
| Instr           | BIFROST                                                         | `chopcal.bifrost`                      | function  |
|                 | a set of choppers, printed as a table                           | `chopcal.lib.ChopperSet`               | class     |
| runtime library | [chopper-lib](https://github.com/mcdotstar/mcstas-chopper-lib/) | `chopcal.lib`                          | submodule |
|                 | chopper parameters                                              | `chopcal.lib.Chopper`                  | class     |
|                 | transmitted inverse-velocity phase space                        | `chopcal.lib.inverse_velocity_windows` | function  |
|                 | transmitted inverse-velocity extrema                            | `chopcal.lib.inverse_velocity_limits`  | function  |
|                 | transmitted wavelength extrema                                  | `chopcal.lib.wavelength_limits`        | function  |


## Placing the band

BIFROST passes a fixed bandwidth of about 1.77 Å, so one number decides where that band
sits. Give it as either the **longest wavelength** or the **lowest energy** you want to
reach the sample — the band runs from there to about 1.77 Å shorter.

```pycon
>>> import chopcal
>>> chopcal.bifrost(wavelength_max=3.0)
name  speed [Hz]  delay [ms]  opening [deg]  open [ms]  path [m]
 ps1         196     3.87771            170      2.409     6.342
 ps2         196       6.087            170      2.409     6.362
 fo1          14     6.13892          38.26      7.591      8.53
 fo2          14     9.54467          52.01      10.32    14.973
 bw1          14     42.8605            161      31.94        78
 bw2         -14     42.8605            161      31.94     78.02
>>> chopcal.lib.wavelength_limits(list(chopcal.bifrost(wavelength_max=3.0).values()))
(1, (1.2115739130694068, 2.9831347854552717))
```

Mind the sense of it: the number you give is the *slow* end of the band. Asking for
`energy_min=7` does not cap the incident energy at 7 meV — it puts 7 meV at the bottom
and lets everything faster through, up to about 31 meV.

The choppers come back by name in beam order, each with a `speed` in Hz, a `delay` in
seconds, an opening `angle` in degrees and a flight `path` in metres. `delay` is when an
opening is centred on the beam, so a chopper passes neutrons for `angle/360/|speed|`
seconds around it and again every `1/speed` after that.



## Installation

```bash
pip install chopcal
```

Or, from the source repository to get the latest development version
```bash
pip install git+https://github.com/g5t/chopper-calculations.git
```

## Describing a chopper

A `Chopper` is `speed`, `delay`, `angle` and `path`: how fast the disk turns in Hz, when
one of its openings is on the beam in seconds, how wide that opening is in degrees, and
how far the disk sits from the source in metres.

`delay` is a time, so it means the same thing at any speed and in either direction of
rotation — which is how a real chopper is set, and what McStas' `DiskChopper` acts on.
It replaces the `phase` in degrees that `chopper-lib` took before its version 2.0.0, and
which it converted to a delay at every point of use anyway.

## Developing

The build fetches `chopper-lib` from GitHub. To build against a local checkout instead —
which you need when changing both together — point `FetchContent` at it:

```bash
pip install --no-build-isolation -ve . \
  -Ccmake.define.FETCHCONTENT_SOURCE_DIR_CHOPPER_LIB=/path/to/mcstas-chopper-lib
pytest test
```

Nothing else changes: the same sources are compiled into the extension modules, so a
modification to `chopper-lib.c` shows up in `chopcal` on the next build. `chopper-lib`
carries a `CHOPPER_LIB_VERSION` macro and `src/choppers.cpp` asserts on it, so building
against a version that means something else by a chopper's fields fails rather than
returning different numbers.
