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
