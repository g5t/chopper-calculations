# Constants

`chopcal.constants` holds every number the BIFROST calculation depends on, with its
provenance. They are derived from the SI/CODATA primitives rather than tabulated, so they
cannot disagree with each other — which they used to: the neutron `h²/2m` appeared as
`81.82` in one line and as `0.1106` (its inverse square root) in the next, and
`0.1106 × √81.82` is `1.000426`, not 1. Asking for a 3 Å band produced settings for
2.99872 Å.

```pycon
>>> from chopcal import constants
>>> constants.H_OVER_M
3956.034006119441
```

## Physical

| name | unit | source |
|---|---|---|
| `PLANCK` | J s | exact by definition of the SI |
| `ELEMENTARY_CHARGE` | C | exact by definition of the SI |
| `NEUTRON_MASS` | kg | CODATA 2022 |
| `PI` | — | spelled out, not from `<numbers>` |

## Derived, in the units the calculations work in

| name | unit | meaning |
|---|---|---|
| `H_OVER_M` | Å m/s | wavelength times velocity |
| `H2_OVER_2M` | meV Å² | energy times wavelength squared |
| `V2K` | — | velocity (m/s) to wavevector (1/Å) |
| `K2V` | — | the inverse, derived as `1/V2K` |

## ESS source

| name | unit |
|---|---|
| `SOURCE_DURATION` | s |
| `SOURCE_FREQUENCY` | Hz |
| `PULSE_HIGH_FLUX_OFFSET` | s |

## BIFROST geometry

Distances from the moderator, opening angles, and the disc geometry the beam aperture is
computed from.

| name | unit |
|---|---|
| `INSTRUMENT_LENGTH`, `PULSE_SHAPING_DISTANCE` | m |
| `FRAME_OVERLAP_1_DISTANCE`, `FRAME_OVERLAP_2_DISTANCE` | m |
| `BANDWIDTH_DISTANCE`, `PAIR_SEPARATION` | m |
| `PULSE_SHAPING_ANGLE`, `FRAME_OVERLAP_1_ANGLE`, `FRAME_OVERLAP_2_ANGLE`, `BANDWIDTH_ANGLE` | degrees |
| `DISK_RADIUS` | m |
| `NARROW_WINDOW_WIDTH`, `NARROW_WINDOW_HEIGHT`, `NARROW_SLIT_HEIGHT` | m |
| `BANDWIDTH_WINDOW_WIDTH`, `BANDWIDTH_WINDOW_HEIGHT`, `BANDWIDTH_SLIT_HEIGHT` | m |

`*_SLIT_HEIGHT` is the radial extent of the opening cut in the disc — McStas
`DiskChopper`'s `yheight` — and `*_WINDOW_*` is the beam inside it, which is smaller.

!!! warning "These are not chopper-lib's conversion constants"

    `constants.V2K` and `constants.K2V` are derived from SI/CODATA and are exact
    reciprocals. chopper-lib compiles against the **McStas runtime's** values, which are
    rounded as separate literals whose product is `0.99999999891`.

    Converting by hand with `constants.V2K` therefore disagrees with `wavelength_windows`
    and `Region.wavelength_ranges` in the eighth digit. Use
    [`chopcal.lib.wavelength_to_inverse_velocity`](api/lib.md) and its inverse, which are
    the library's own — and which, for the same reason, round-trip 1.1 parts in a billion
    short.
