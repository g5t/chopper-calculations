# Chopcal

## Description
Exposes chopper calculations to Python which were otherwise hidden in McStas instruments.


## Supported components

| source          | name                                                            | component                                 | type      |
|-----------------|-----------------------------------------------------------------|-------------------------------------------|-----------|
| Instr           | BIFROST                                                         | `chopcal.bifrost`                         | function  |
| runtime library | [chopper-lib](https://github.com/mcdotstar/mcstas-chopper-lib/) | `chopcal.mcstas`                          | submodule |
|                 | chopper parameters                                              | `chopper.mcstas.Chopper`                  | class     |
|                 | transmitted inverse-velocity phase space                        | `chopper.mcstas.inverse_velocity_windows` | function  |
|                 | transmitted inverse-velocity extremea                           | `chopper.mcstas.inverse_velocity_limits`  | function  |
|                 | transmitted wavelength extremea                                 | `chopper.mcstas.wavelength_limits`        | function  |



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
