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
|                 | transmitted wavelength bands                                    | `chopcal.lib.wavelength_windows`       | function  |
|                 | transmitted wavelength extrema                                  | `chopcal.lib.wavelength_limits`        | function  |
|                 | beam width on a disk, in degrees                                | `chopcal.lib.beam_aperture`            | function  |
|                 | whether a parked disk stands open                               | `chopcal.lib.Chopper.parked_is_open`   | method    |
|                 | allowed (inverse velocity, time) bins                           | `chopcal.lib.inverse_velocity_time_mask` | function |
|                 | the allowed fraction of a signal                                | `chopcal.lib.unmasked_probability`     | function  |
|                 | direct sampling over a finished mask                            | `chopcal.lib.MaskSampler`              | class     |


## Placing the band

BIFROST passes a fixed bandwidth of about 1.91 Å, so one number decides where that band
sits. Give it as either the **longest wavelength** or the **lowest energy** you want to
reach the sample — the band runs from there to about 1.77 Å shorter.

```pycon
>>> import chopcal
>>> chopcal.bifrost(wavelength_max=3.0)
name  speed [Hz]  delay [ms]  beam [deg]  opening [deg]  open [ms]  aperture [deg]  path [m]
 ps1         196      3.8781           0            170      2.409           5.651     6.342
 ps2         196      6.0874           0            170      2.409           5.651     6.362
 fo1          14     6.13997           0          38.26      7.591           5.651      8.53
 fo2          14     9.54764           0          52.01      10.32           5.651    14.973
 bw1          14     42.8823           0            161      31.94           13.38        78
 bw2         -14     42.8823           0            161      31.94           13.38     78.02
>>> chopcal.lib.wavelength_limits(list(chopcal.bifrost(wavelength_max=3.0).values()))
(1, (1.2126790814161053, 2.984239670498114))
```

Mind the sense of it: the number you give is the *slow* end of the band. Asking for
`energy_min=7` does not cap the incident energy at 7 meV — it puts 7 meV at the bottom
and lets everything faster through, up to about 31 meV.

The choppers come back by name in beam order. Each of BIFROST's is a single opening
centred on the beam, so `beam` is 0 and `edges` is `[-width/2, +width/2]`; the table's
`open` column is how long the beam spends inside the openings each turn, and `delay` is
when the disk point at `beam` is on the beam path. A chopper therefore passes neutrons
for `open` seconds around `delay`, and again every `1/speed` after that — a little
longer than that, in fact, because `aperture` is how wide the beam is on the disk and
every window opens half of it early and closes half of it late.

Note the second window in the output above. It is not a band BIFROST delivers; see
[The band list, and a caveat](#the-band-list-and-a-caveat).



## Installation

```bash
pip install chopcal
```

Or, from the source repository to get the latest development version
```bash
pip install git+https://github.com/g5t/chopper-calculations.git
```

## Describing a chopper

A `Chopper` is `speed`, `delay`, `beam`, `edges`, `path` and `aperture` — the disk
described the way the NeXus `NXdisk_chopper` standard and McStas' `CollectorDiskChopper`
component describe one, so the disk's own numbers go in unchanged.

`edges` are the slit edges in degrees from the disk's zero mark, two per opening,
strictly increasing, opening edge first; `beam` is the angle from that mark to where the
beam crosses the disk. An edge at angle `a` is on the beam at

```
t(a) = delay + (beam - a) / (360 * speed)
```

and every `1/|speed|` seconds after that. Only `beam - a` appears, so there is no
conversion to do. A slit straddling the mark is written past 360 — `[350, 370]` rather
than `[350, 10]` — so the pairs stay ordered and each width is a difference. Negative
angles are fine: one opening astride the mark reads better as `[-85, 85]`.

`Chopper.centred(speed, delay, width, path)` writes a single opening centred on the beam,
which is what the `angle` field of `chopper-lib` 2.x meant and all `chopcal.bifrost`
produces.

`aperture` is the angular width of the beam where it crosses the disk, which widens every
window by half of it at each end; 0, the default, is the point beam every other field
describes. It is an angle rather than a width in metres because that is what the disk
sees — see `chopper-lib`'s header for the conversion, which is not a plain division.

`delay` is a time, so it means the same thing at any speed and in either direction of
rotation — which is how a real chopper is set, and what McStas' `DiskChopper` acts on.
It replaces the `phase` in degrees that `chopper-lib` took before its version 2.0.0, and
which it converted to a delay at every point of use anyway.

A disk that is not turning is open or shut for good, and the two mean opposite things: the
window and mask functions leave one parked *open* out of the calculation, having no period
to constrain anything, while one parked *shut* is a beam stop and empties the answer.
`Chopper.parked_is_open` tells them apart, which an empty result cannot.

## Masking and sampling

`inverse_velocity_time_mask` returns which (inverse velocity, time) bins the train passes,
as an array shaped **(time, inverse velocity)** — time is the slow axis, the opposite of
the argument order. It takes bin *edges*, so the mask is one smaller in each direction.

For BIFROST that mask allows well under a percent of the plane, so a source that draws
uniformly and discards what the mask excludes spends almost its whole ray budget. A
`MaskSampler` draws from the allowed cells instead, never rejects, and carries the
`acceptance` every drawn ray's weight must be multiplied by to keep the estimate
unbiased. Over the same 200 × 300 grid that is a ~300× smaller error for the same number
of rays.

`acceptance` is a ratio of areas, so it is exact — but only while both coordinates are
drawn uniformly and independently of each other and of everything else the caller samples.
A source that picks its emission time from a window centred on the neutron's own velocity,
as McStas' `ESS_butterfly` does under time focusing, breaks that, and no single number
corrects it. It is also not `unmasked_probability`, which is the allowed fraction of a
*weighted* signal — the transmission an instrument sees, and the wrong normalisation for
this.

These three need `numpy`, and only when called; nothing else in `chopcal` does.

## The band list, and a caveat

`wavelength_windows` gives the bands the train passes; `wavelength_limits` gives the
envelope around them. Reach for the first when there is more than one, because an envelope
spanning two bands also spans the gap between them.

Both are an **over**-approximation. `chopper-lib` works out each disk's admissible inverse
velocities letting the emission time range over the whole pulse independently, then
intersects those ranges across disks — and an intersection of projections is a superset of
the projection of the intersection. It can therefore report a band that no single emission
time actually delivers.

It is tight wherever the disks leave wide, overlapping emission windows, which is the
usual case and includes every band an instrument is set up to pass. BIFROST with its real
beam is the exception: a ~0.03 Å sliver near 38 Å, where the frame overlap disks each
leave only tens of microseconds and those windows do not overlap. Three independent checks
agree nothing gets through — a brute-force scan over emission time, an exact interval
intersection, and `inverse_velocity_time_mask`, which keeps the two coordinates coupled and
cannot make the error.

So: where the mask and the windows disagree, the mask is right. `test_mask.py` records
this one.

## Describing a disk's beam

`Chopper.aperture` is the angular width of the beam where it crosses the disk, and
`beam_aperture` computes it from the geometry:

```pycon
>>> from chopcal.lib import beam_aperture
>>> beam_aperture(radius=0.35, slit_height=0.09846, window_width=0.060, window_height=0.090)
13.379...
```

`slit_height` is the radial extent of the opening — McStas `DiskChopper`'s `yheight`,
which puts the beam centre at `radius - yheight/2` — and the window is the beam inside it.
The largest angle belongs to the window's *inner* corners, nearest the spindle, where a
given width subtends the most angle. Taking the width over the beam-crossing radius
instead misses the height entirely: 11.4° where the real figure is 13.4°.

`chopcal.bifrost` fills this in from the instrument's own disk geometry. Pass
`apertures=False` for the pencil beam every release before this one described, and the
narrower band that went with it — `BIFROST_BANDWIDTH_POINT_BEAM` rather than
`BIFROST_BANDWIDTH`.

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
carries a `CHOPPER_LIB_VERSION` macro and `src/chopper.h` asserts on it, so building
against a version that means something else by a chopper's fields fails rather than
returning different numbers. `chopcal` needs 4.1.0 or newer.

Only `chopper-lib.c` is compiled with `CHOPPER_LIB_DEFINITIONS`, not the whole target:
those are `V2K`, `K2V` and `PI` as macros, and `chopcal::constants` declares C++ variables
by the same names.
