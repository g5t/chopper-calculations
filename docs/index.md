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

## What it is not

It is not a disc. chopcal asks *when a chopper is open on the beam axis*, widened by the
beam's own width; it does not trace a neutron through a slit. For what a real disc passes,
simulate one.

It is not a transmission calculation either. The answers here are geometric — which
`(wavelength, emission time)` pairs can get through — with no spectrum in them. Fold in
the source's own brilliance afterwards.
