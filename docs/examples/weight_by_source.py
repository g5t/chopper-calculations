"""A gridded source intensity through the transmitted region, and its arrival times."""
import numpy as np
import chopcal
from chopcal.lib import Region

train = list(chopcal.bifrost(wavelength_max=3.0).values())
source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)
region = source.transmit(train)

# --8<-- [start:grid]
# bin edges over the source rectangle: inverse velocity in s/m, emission time in s
corners = [u for u, _ in source.polygons[0].vertices]
inverse_velocities = np.linspace(min(corners), max(corners), 401)
times = np.linspace(0.0, 3e-3, 61)

# the amount in each (time, inverse velocity) cell -- here a stand-in that peaks early
centres = 0.5 * (times[1:] + times[:-1])
signal = np.exp(-centres / 1e-3)[:, None] * np.ones(len(inverse_velocities) - 1)

coverage = region.coverage(inverse_velocities, times)
print(coverage.shape, f'{coverage.max():.3f}')
print(f'transmitted {region.integrate(signal, inverse_velocities, times):.4f} '
      f'of {signal.sum():.1f}')
# --8<-- [end:grid]

# --8<-- [start:arrivals]
frame = 1 / 14
arrival_times = np.linspace(0.0, frame, 715)        # 0.1 ms bins over one frame
at_sample = region.arrivals(signal, inverse_velocities, times, 162.0, arrival_times,
                            period=frame)
profile = at_sample.sum(axis=1)                     # over inverse velocity
peak = arrival_times[np.argmax(profile)]
print(at_sample.shape, f'peak at {peak * 1e3:.1f} ms')
# --8<-- [end:arrivals]

assert np.isclose(at_sample.sum(), region.integrate(signal, inverse_velocities, times))

# --8<-- [start:moments]
unfolded = np.linspace(0.10, 0.18, 801)               # 0.1 ms bins, no folding
moments = region.arrival_moments(signal, inverse_velocities, times, 162.0, unfolded)
busiest = np.nanargmax(moments.amount)
print(f'busiest bin {unfolded[busiest] * 1e3:.1f} ms: '
      f'{moments.wavelength[busiest]:.4f} ± {moments.sigma[busiest]:.4f} Å')
# --8<-- [end:moments]

assert np.allclose(moments.amount,
                   region.arrivals(signal, inverse_velocities, times, 162.0,
                                   unfolded).sum(axis=1))

# --8<-- [start:weighted]
from chopcal import constants as c
from chopcal.lib import Chopper, beam_aperture

# the ramp is set by the angle the beam subtends at its centre, not at its inner corner
narrow = beam_aperture(c.DISK_RADIUS, c.NARROW_SLIT_HEIGHT, c.NARROW_WINDOW_WIDTH,
                       c.NARROW_WINDOW_HEIGHT, at='centre')
wide = beam_aperture(c.DISK_RADIUS, c.BANDWIDTH_SLIT_HEIGHT, c.BANDWIDTH_WINDOW_WIDTH,
                     c.BANDWIDTH_WINDOW_HEIGHT, at='centre')
ramped = [Chopper(d.speed, d.delay, d.beam, list(d.edges), d.path,
                  wide if name.startswith('bw') else narrow)
          for name, d in chopcal.bifrost(wavelength_max=3.0).items()]
# the pulse-shaping pair and the bandwidth pair each share a housing
weighted = source.transmit_weighted(ramped, groups=[0, 0, 1, 2, 3, 3])
point = source.transmit(list(chopcal.bifrost(wavelength_max=3.0, apertures=False).values()))
print(f'support {region.area:.4e}, point beam {point.area:.4e}, '
      f'weighted {weighted.integral:.4e} s^2/m')
print(f'transmitted {weighted.integrate(signal, inverse_velocities, times):.4f}')
softer = weighted.arrivals(signal, inverse_velocities, times, 162.0, arrival_times,
                           period=frame)
# --8<-- [end:weighted]

assert np.all(softer <= at_sample + 1e-12)
