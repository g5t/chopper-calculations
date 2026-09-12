# Install and first band

## Install

```shell
pip install chopcal
```

Wheels are published for Linux, macOS and Windows on CPython 3.11 and newer. There is
nothing to compile and nothing else to install: chopcal has **no runtime dependencies**.

Two optional packages unlock extras, and only when you call something that needs them:

| package | what it adds |
|---|---|
| `numpy` | the grid mask, and drawing samples in bulk |
| `scipp` | `Chopper.quantities`, the fields with units attached |

## Set a train and read its band

BIFROST passes a fixed band about 1.91 Å wide, so a single number decides where it sits.
Give it as the **longest wavelength** you want to reach the sample:

```python
--8<-- "quickstart.py:settings"
```

```
name  speed [Hz]  delay [ms]  beam [deg]  opening [deg]  open [ms]  aperture [deg]  path [m]
 ps1         196      3.8781           0            170      2.409           5.651     6.342
 ps2         196      6.0874           0            170      2.409           5.651     6.362
 fo1          14     6.13997           0          38.26      7.591           5.651      8.53
 fo2          14     9.54764           0          52.01      10.32           5.651    14.973
 bw1          14     42.8823           0            161      31.94           13.38        78
 bw2         -14     42.8823           0            161      31.94           13.38     78.02
```

It is a `dict` in every respect — `settings['ps1']`, `.values()`, `**settings` all work —
that prints as a table rather than as six memory addresses.

Now ask what those discs let through:

```python
--8<-- "quickstart.py:band"
```

```
1.1454 to 3.0515 AA, 1.9062 AA wide
```

Note the sense of the argument: the number you give is the **slow** end. Asking for
`energy_min=7` does not cap the incident energy at 7 meV — it puts 7 meV at the bottom and
lets everything faster through, up to about 31 meV.

## Read one disc

```python
--8<-- "quickstart.py:one_chopper"
```

```
14 Hz anticlockwise, 161 deg opening centred on the beam on the beam at 42.88 ms, 78 m from the source
0.04288227154711843
[-80.5, 80.5]
```

`delay` is in seconds, not the milliseconds the table prints; `edges` are the disc's own
slit angles in degrees. [Describe a chopper](how-to/describe-a-chopper.md) goes through
the rest of the fields.

## Where next

* [Core concepts](concepts.md) — the coordinates, the delay, the aperture.
* [Find the band a train passes](how-to/find-the-band.md) — for any discs, not just BIFROST's.
* [Build the transmitted region](how-to/transmitted-region.md) — the exact answer, not an
  approximation of it.
