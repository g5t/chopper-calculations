# Describe a chopper

A `Chopper` is `speed`, `delay`, `beam`, `edges`, `path` and `aperture` — the disc
described the way the NeXus `NXdisk_chopper` standard and McStas' `CollectorDiskChopper`
describe one, so the disc's own numbers go in unconverted.

## One opening, centred on the beam

The common case, and what every BIFROST disc is:

```python
--8<-- "describe_a_chopper.py:centred"
```

`Chopper.centred` writes `beam = 0` and `edges = [-width/2, +width/2]`. That pair is
exactly what the single `angle` field meant before chopper-lib 4.0.0.

## Several openings

Give the disc's own slit angles: two per opening, strictly increasing, opening edge first.

```python
--8<-- "describe_a_chopper.py:edges"
```

`beam` is the angle from the disc's mark to where the beam crosses it, so the openings and
the beam crossing are described separately — there is no conversion to do.

An opening straddling the mark is written **past 360**, never wrapped:

```python
--8<-- "describe_a_chopper.py:across_the_mark"
```

`[350, 370]`, not `[350, 10]`, so the pairs stay ordered and each width is a difference.
Negative angles are fine too: a single opening astride the mark reads better as
`[-85, 85]`.

!!! tip "Two mistakes that are refused rather than computed"

    An **odd** number of edges, and edges that do not strictly increase. chopper-lib checks
    neither: it divides the count by two, so an odd list silently drops its last edge, and
    it reads unordered edges as negative widths. Both come back as a `ValueError` here.

    After assigning to `edges` directly, call `chopper.validate()` — a plain attribute
    cannot check itself.

## The aperture

`aperture` is how wide the beam is on the disc, in degrees about its spindle.
`beam_aperture` works it out from the geometry:

```python
--8<-- "describe_a_chopper.py:aperture"
```

```
13.3796 degrees
```

`slit_height` is the radial extent of the opening cut in the disc — McStas `DiskChopper`'s
`yheight`, which puts the beam centre at `radius - yheight/2`. The window is the beam
inside it, and is smaller.

The largest angle belongs to the window's **inner** corners, nearest the spindle, where a
given width subtends the most angle. Taking the width over the beam-crossing radius instead
misses the height entirely: 11.4° where the real figure is 13.4°.

Leave `aperture` at 0 for a pencil beam, which is what every field but this one describes.

## A parked disc

A disc that is not turning is open or shut for good, and the two mean opposite things: one
parked **open** has no period and so constrains nothing, while one parked **shut** is a
beam stop that empties the answer.

```python
--8<-- "describe_a_chopper.py:parked"
```

An empty result does not say which happened. `parked_is_open` does. It ignores `speed`,
because a turning disc stands open once a period whatever its angles are.

## Units, at a glance

| field | unit |
|---|---|
| `speed` | Hz, signed — negative turns the disc the other way |
| `delay` | s |
| `beam` | degrees from the disc's mark |
| `edges` | degrees from the disc's mark |
| `path` | m from the source |
| `aperture` | degrees |

`chopper.quantities` gives the same values as `scipp` variables with those units attached,
if scipp is installed. It is not a dependency and does not become one.
