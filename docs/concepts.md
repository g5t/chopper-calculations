# Core concepts

Five ideas carry the whole library. None is difficult; all of them are places where a
plausible-looking guess is wrong.

## The plane everything is drawn in

Every answer chopcal gives is about one two-dimensional space: **inverse velocity** against
**emission time**.

Inverse velocity `a` (s/m) rather than wavelength, because it makes the arithmetic linear.
A neutron emitted at time `t` reaches a disc `L` metres downstream at

```
t + L * a
```

so *a disc's acceptance is a straight band in this plane* — bounded by two parallel lines
of slope `−L`. Two discs at different distances cross their bands at an angle, and what a
train passes is the patch they share.

In wavelength the same statement needs a curve. That is the only reason for the choice; the
conversion is one multiplication, and
[`wavelength_to_inverse_velocity`](reference/api/lib.md) does it the way the library does.

!!! note "Emission time, not arrival time"

    `t` is when the neutron left the source, not when it reached anything. A monitor
    downstream sees arrival times; those are this plane sheared by the flight time, and the
    two cannot be compared directly.

## A delay is a time, not a phase

`Chopper.delay` is **when** the disc's mark is on the beam, in seconds. It used to be a
phase in degrees, which the library divided by `360 × |speed|` at every point of use.

A delay says the same thing without needing to know the speed — and, unlike a phase, it
does not wrap at one revolution, so it can exceed a period. It is also what a real chopper
is set with, and what McStas' `DiskChopper` takes.

Openings recur at `delay + n / speed` for integer `n`. The sign of `speed` sets the
direction of rotation and does **not** move `delay`: a time means the same thing whichever
way the disc turns.

## A disc is its slit edges

`Chopper.edges` are the disc's own slit angles in degrees from its top-dead-centre mark,
two per opening, strictly increasing, opening edge first. `Chopper.beam` is the angle from
that mark to where the beam crosses the disc.

This is the NeXus `NXdisk_chopper` description, and McStas'
`CollectorDiskChopper` takes the same numbers — so a disc goes in unconverted. An edge at
angle `a` is on the beam at

```
t(a) = delay + (beam - a) / (360 * speed)
```

Only `beam − a` appears, which is what makes the conversion unnecessary.

A single opening `w` degrees wide centred on the beam is `beam = 0` and
`edges = [-w/2, +w/2]`, which is what `Chopper.centred` writes and what the old single
`angle` field meant.

## The beam has width, and the width acts on time only

`Chopper.aperture` is how wide the beam is **on the disc**, in degrees about its spindle.
It widens every window by half of it at each end.

It is an angle rather than a width in metres because that is what the disc sees. A neutron
crossing `w` metres to one side of the beam centre reaches an opening edge `w / d` radians
early or late, `d` being the distance from the spindle. The largest such angle belongs to
the window's *inner* corners — nearest the spindle, where a given width subtends the most
angle — and [`beam_aperture`](how-to/describe-a-chopper.md#the-aperture) works it out from
the geometry.

The important half: **a wide beam changes when a neutron gets through, never how fast it
is going.** Widening the answer in both directions instead — growing a grid mask by whole
bins, say — admits inverse velocities no disc ever passes.

## Three answers, and only one is exact

Asking "what does this train pass?" has three implementations, and they disagree:

| | what it gives | how it is wrong |
|---|---|---|
| `wavelength_windows`, `wavelength_limits` | the bands, quickly | **over**-reports: can name a band no single emission time delivers |
| `inverse_velocity_time_mask` | a picture on a grid | quantises: loses thin channels, keeps partly covered bins whole |
| `Region` | the region itself, as polygons | it is not |

For BIFROST with a real beam the first reports a band near 38 Å that nothing passes.
[Three answers to one question](explanation/three-answers.md) is the whole story, and
[Why the window functions over-report](explanation/over-reporting.md) is the proof.

Reach for `Region` for anything feeding a number into a calculation; the other two remain
the cheapest way to sanity-check a band and the only one that gives you an array to plot.
