# chopcal

Chopper calculations that used to be locked inside McStas instruments, in Python.

A chopper train decides two things an instrument scientist cares about: **where the
wavelength band sits**, and **which neutrons a source need not bother emitting**. chopcal
answers both, using the same C the instrument uses — [chopper-lib][cl] — so the numbers
here and the numbers in a simulation cannot drift apart.

```pycon
>>> import chopcal
>>> settings = chopcal.bifrost(wavelength_max=3.0)
>>> chopcal.lib.wavelength_windows(list(settings.values()))[0]
(1.1453582740163402, 3.0515432206219723)
```

[cl]: https://github.com/mcdotstar/mcstas-chopper-lib

## Where to start

<div class="grid cards" markdown>

-   :material-rocket-launch: **[Getting started](getting-started.md)**

    Install it and get a band out of it.

-   :material-lightbulb: **[Core concepts](concepts.md)**

    The two coordinates everything is drawn in, why a delay is not a phase, and what a
    disc's `aperture` is for.

-   :material-wrench: **[How-to guides](how-to/place-the-band.md)**

    Task by task: place a band, describe a disc, build the transmitted region, sample it
    without wasting rays.

-   :material-book-open-variant: **[Explanation](explanation/three-answers.md)**

    Three functions answer "what does this train pass?" and they do not agree. Which to
    trust, and why.

</div>

## What it is for

**Setting a chopper train.** `chopcal.bifrost` gives the six BIFROST discs their speeds
and delays from one number — the slow edge of the band you want.

**Asking what a train passes.** Given any set of discs, chopcal will tell you the
wavelength bands that get through, on a grid, or exactly.

**Not simulating what cannot arrive.** A train that passes a percent of the plane makes a
source spend ninety-nine percent of its rays on neutrons the discs stop. chopcal computes
the region that does get through, and a sampler that draws from it — the same machinery
McStas' own `Polygon_ESS_butterfly` uses.

**Asking how much arrives.** Give it a gridded source, such as a moderator's tabulated
spectrum, and chopcal integrates it over the transmitted region exactly: the total, what
reaches a monitor or the sample by arrival time and wavelength, or the mean and spread of
wavelength in each arrival-time bin, the numbers a wavelength lookup table stores. With a finite beam it also
weights each neutron by how likely it is to get through while a disc's edge is still
sweeping the beam. See [Weight the region by a source](how-to/weight-by-a-source.md).

**Setting the real choppers.** `chopcal.nicos` gives BIFROST's speeds and experiment
delays keyed by the NICOS devices that take them. See
[Set the choppers from NICOS](how-to/set-the-choppers-from-nicos.md).

## What it is not

It is not a disc. chopcal asks *when a chopper is open on the beam axis*, widened by the
beam's own width. Its weighted transmission models an edge crossing a uniformly lit beam
as a linear ramp; it does not trace a neutron through a slit. That agrees with McStas for
BIFROST's wedge-shaped openings to within counting error. For openings of another shape,
or a beam that is not evenly lit, simulate the disc.

It is not a beamline either. Every neutron flies one path from source to sample; a
guide's spread of path lengths smears arrival times, and only `Region.transmit`'s
`path_spreads` accounts for it, as a support rather than a distribution.
