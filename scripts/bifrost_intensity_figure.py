"""Pre-compute the BIFROST source-intensity figure for the documentation.

Writes, into ``docs/assets/bifrost-intensity/``:

- ``source.png``: the ESS source BIFROST sees, as an image filling the source rectangle,
  the background of the interactive figure's left panel;
- ``frames.json``: for each setting on a grid of lowest energy, pulse-shaping time and
  pulse delay, the transmitted region's outline and what reaches the sample, by arrival
  time and by wavelength -- weighted by transmission, and counting the whole region open;
- ``static.png``: the whole figure for the default setting, for readers without JavaScript.

`docs/javascripts/bifrost-intensity.js` draws the interactive version from the first two.
This needs ``tof`` (for the source library, which it downloads) and ``matplotlib``, so it
is run by hand rather than by the test suite:

    python scripts/bifrost_intensity_figure.py
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / 'docs' / 'assets' / 'bifrost-intensity'

#: The source rectangle: wide enough for every setting on the grid below.
WAVELENGTH = (0.5, 7.0)          # Å
LATEST_EMISSION = 3.5e-3         # s
WAVELENGTH_BIN = 0.025           # Å
TIME_BIN = 20e-6                 # s
ARRIVAL_BINS = 200

ENERGIES = (2.0, 3.0, 4.5, 7.0, 10.0)          # meV, the slow edge of the band
SHAPING_TIMES = (0.1, 0.2, 0.5, 1.0)          # ms
PULSE_DELAYS = (-300.0, 0.0, 300.0)            # µs
#: The pulse-shaping disks BIFROST was designed for, and the substitutes installed now.
PSC_SPEEDS = (196.0, 126.0)                    # Hz
DEFAULT = (4.5, 0.2, 0.0, 196.0)

#: The ESS source distributions in `tof` are sampled this far from the moderator, and
#: their emission times counted there; chopcal measures paths from the moderator.
SOURCE_PLANE = 0.05  # m
GROUPS = [0, 0, 1, 2, 3, 3]   # the pulse-shaping and bandwidth pairs share housings
NAMES = ('ps1', 'ps2', 'fo1', 'fo2', 'bw1', 'bw2')
ENERGY_TICKS = (1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 80, 120, 200, 300)


def source_grid(wavelength_edges, time_edges):
    """The tabulated source summed into these bins, (time, wavelength)."""
    import scipp as sc
    from tof.facilities import get_source_path
    data = sc.io.load_hdf5(get_source_path('ess-bifrost'))
    for dim in ('wavelength', 'birth_time'):
        centres = data.coords[dim]
        step = centres[1] - centres[0]
        data.coords[dim] = sc.concat([centres - 0.5 * step, centres[-1:] + 0.5 * step], dim)
    data.coords['birth_time'] = data.coords['birth_time'].to(unit='s')
    out = data.rebin(
        wavelength=sc.array(dims=['wavelength'], values=wavelength_edges, unit='angstrom'),
        birth_time=sc.array(dims=['birth_time'], values=time_edges, unit='s'))
    return np.ascontiguousarray(out.transpose(['birth_time', 'wavelength']).values)


def trains(energy_min, shaping_time, pulse_delay, psc_speed_max):
    """BIFROST's discs from the source plane: inner-corner apertures for the region's
    edge, centre apertures for the ramps, and every delay moved by `pulse_delay`."""
    import chopcal
    from chopcal import constants as c
    from chopcal.lib import Chopper, beam_aperture
    settings = chopcal.bifrost(energy_min=energy_min, shaping_time=shaping_time,
                               pulse_shaping_speed_max=psc_speed_max)
    centre = {
        False: beam_aperture(c.DISK_RADIUS, c.NARROW_SLIT_HEIGHT, c.NARROW_WINDOW_WIDTH,
                             c.NARROW_WINDOW_HEIGHT, at='centre'),
        True: beam_aperture(c.DISK_RADIUS, c.BANDWIDTH_SLIT_HEIGHT,
                            c.BANDWIDTH_WINDOW_WIDTH, c.BANDWIDTH_WINDOW_HEIGHT, at='centre'),
    }

    def moved(disc, aperture):
        return Chopper(disc.speed, disc.delay + pulse_delay, disc.beam, list(disc.edges),
                       disc.path - SOURCE_PLANE, aperture)

    return ([moved(settings[n], settings[n].aperture) for n in NAMES],
            [moved(settings[n], centre[n.startswith('bw')]) for n in NAMES])


def frame(amounts, wavelength_edges, time_edges, energy_min, shaping_ms, pulse_us, psc):
    from chopcal import constants as c
    from chopcal.lib import Region, wavelength_to_inverse_velocity
    per_angstrom = wavelength_to_inverse_velocity(1.0)
    u = wavelength_edges * per_angstrom
    rectangle = Region.rectangle(u[0], u[-1] - u[0], time_edges[0], time_edges[-1] - time_edges[0])
    support_train, ramped_train = trains(energy_min, shaping_ms * 1e-3, pulse_us * 1e-6, psc)
    support = rectangle.transmit(support_train)
    weighted = rectangle.transmit_weighted(ramped_train, GROUPS)
    distance = c.INSTRUMENT_LENGTH - SOURCE_PLANE
    if not support:
        return {'outline': [], 'fraction': 0.0}
    reached = [t + distance * a for p in support.polygons for a, t in p.vertices]
    first, last = min(reached), max(reached)
    margin = 0.03 * (last - first)
    arrivals = np.linspace(first - margin, last + margin, ARRIVAL_BINS + 1)
    step_ms = (arrivals[1] - arrivals[0]) * 1e3

    def at_sample(region):
        return region.arrivals(amounts, u, time_edges, distance, arrivals).sum(axis=1) / step_ms

    def by_wavelength(region):
        return (amounts * region.coverage(u, time_edges)).sum(axis=0) / WAVELENGTH_BIN

    wl_support, wl_weighted = by_wavelength(support), by_wavelength(weighted)
    nonzero = np.nonzero(wl_support)[0]
    lo, hi = max(nonzero[0] - 2, 0), min(nonzero[-1] + 3, len(wl_support))

    def rounded(values):
        return [float(f'{v:.4g}') for v in values]

    return {
        'outline': [[[round(t * 1e3, 5), round(a / per_angstrom, 5)] for a, t in p.vertices]
                    for p in support.polygons],
        'arrival': {'start': round(arrivals[0] * 1e3, 4), 'step': round(step_ms, 6),
                    'weighted': rounded(at_sample(weighted)),
                    'support': rounded(at_sample(support))},
        'wavelength': {'start': round(float(wavelength_edges[lo]), 5), 'step': WAVELENGTH_BIN,
                       'weighted': rounded(wl_weighted[lo:hi]),
                       'support': rounded(wl_support[lo:hi])},
        'fraction': float(f'{weighted.integrate(amounts, u, time_edges) / amounts.sum():.4g}'),
        'psc_speed': abs(support_train[0].speed),
    }


def energy(wavelength):
    """meV from Å, safe at the zero a matplotlib secondary axis probes with."""
    from chopcal import constants as c
    x = np.asarray(wavelength, dtype=float)
    return np.where(x > 0, c.H2_OVER_2M / np.where(x > 0, x, 1.0) ** 2, np.inf)


def wavelength_of(energy_mev):
    """Å from meV, the inverse of `energy`."""
    from chopcal import constants as c
    x = np.asarray(energy_mev, dtype=float)
    return np.where(x > 0, np.sqrt(c.H2_OVER_2M / np.where(x > 0, x, 1.0)), np.inf)


def write_source_image(density):
    """The source as a bare image, top row the longest wavelength, filling the rectangle."""
    import matplotlib.pyplot as plt
    rgba = plt.get_cmap('viridis')(density.T[::-1] / density.max())
    plt.imsave(OUT / 'source.png', rgba)


def write_static(amounts, wavelength_edges, time_edges, density, data):
    """The whole figure for the default setting, drawn by matplotlib."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure(figsize=(12, 6), layout='constrained')
    grid = fig.add_gridspec(2, 2, width_ratios=(1.15, 1))
    left, top, bottom = (fig.add_subplot(grid[:, 0]), fig.add_subplot(grid[0, 1]),
                         fig.add_subplot(grid[1, 1]))
    t, w = time_edges * 1e3, wavelength_edges
    mesh = left.pcolormesh(t, w, density.T, cmap='viridis', rasterized=True)
    for polygon in data['outline']:
        left.add_patch(Polygon(polygon, fill=False, edgecolor='white', linewidth=1.6))
        left.add_patch(Polygon(polygon, fill=False, edgecolor='crimson', linewidth=0.9))
    left.set(xlabel='emission time at the source [ms]', ylabel='wavelength [Å]',
             title='source, and the region the train transmits')
    secondary = left.secondary_yaxis('right', functions=(energy, wavelength_of))
    ticks = [e for e in ENERGY_TICKS if energy(w[-1]) <= e <= energy(w[0])]
    secondary.set_yticks(ticks, labels=[f'{e:g}' for e in ticks])
    secondary.set_ylabel('energy [meV]')
    fig.colorbar(mesh, ax=left, location='bottom', shrink=0.8,
                 label='source intensity [per ms per Å, arb.]')
    for ax, key, colour, label in ((top, 'arrival', 'C0', 'arrival time at the sample [ms]'),
                                    (bottom, 'wavelength', 'C1', 'wavelength [Å]')):
        d = data[key]
        edges = d['start'] + d['step'] * np.arange(len(d['weighted']) + 1)
        ax.stairs(d['support'], edges, color='0.6', linestyle='--', label='whole region open')
        ax.stairs(d['weighted'], edges, color=colour, linewidth=1.5,
                  label='weighted by transmission')
        ax.set_xlim(edges[0], edges[-1])
        ax.set_xlabel(label)
        ax.legend(loc='upper right', fontsize='small')
    top.set_ylabel('intensity [per ms, arb.]')
    bottom.set_ylabel('intensity [per Å, arb.]')
    top.set_title(f'lowest energy {DEFAULT[0]} meV, shaping {DEFAULT[1]} ms, pulse delay '
                  f'{DEFAULT[2]:g} µs, PSC {DEFAULT[3]:g} Hz:\n{data["fraction"]:.2%} of the '
                  f'source rectangle reaches the sample', fontsize='medium')
    fig.savefig(OUT / 'static.png', dpi=90)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wavelength_edges = np.arange(WAVELENGTH[0], WAVELENGTH[1] + 0.5 * WAVELENGTH_BIN,
                                 WAVELENGTH_BIN)
    time_edges = np.arange(0.0, LATEST_EMISSION + 0.5 * TIME_BIN, TIME_BIN)
    amounts = source_grid(wavelength_edges, time_edges)
    density = amounts / (TIME_BIN * 1e3 * WAVELENGTH_BIN)
    write_source_image(density)

    frames = {}
    for e, s, p, m in product(ENERGIES, SHAPING_TIMES, PULSE_DELAYS, PSC_SPEEDS):
        frames[f'{e:g}|{s:g}|{p:g}|{m:g}'] = frame(amounts, wavelength_edges, time_edges,
                                                  e, s, p, m)
    lo, hi = WAVELENGTH
    data = {
        'time': [0.0, LATEST_EMISSION * 1e3], 'wavelength': list(WAVELENGTH),
        'energy_ticks': [[float(wavelength_of(e)), f'{e:g}'] for e in ENERGY_TICKS
                         if energy(hi) <= e <= energy(lo)],
        'energies': ENERGIES, 'shaping_times': SHAPING_TIMES, 'pulse_delays': PULSE_DELAYS,
        'psc_speeds': PSC_SPEEDS,
        'default': DEFAULT, 'frames': frames,
    }
    (OUT / 'frames.json').write_text(json.dumps(data, separators=(',', ':')))
    write_static(amounts, wavelength_edges, time_edges, density,
                 frames['|'.join(f'{v:g}' for v in DEFAULT)])
    for path in sorted(OUT.iterdir()):
        print(f'{path.name:12s} {path.stat().st_size / 1024:7.1f} kB')


if __name__ == '__main__':
    main()
