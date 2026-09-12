"""Drawing from the region instead of throwing rays away, without biasing the answer."""
import numpy as np

import chopcal
from chopcal.lib import Region

train = list(chopcal.bifrost(wavelength_max=3.0).values())
source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)
region = source.transmit(train)

# --8<-- [start:sampler]
sampler = region.sampler(source.area)
print(f'{sampler.count} triangles, acceptance {sampler.acceptance:.6e}')

inverse_velocity, time = sampler.sample(200_000, np.random.default_rng(1))
print(inverse_velocity.shape, time.shape)
# --8<-- [end:sampler]

# --8<-- [start:unbiased]
def tally(inverse_velocity, time):
    """Some smooth downstream quantity; nothing depends on which one."""
    return 1.0 + 10.0 * inverse_velocity + 3.0 * time


# uniform over the whole rectangle, keeping only what the train passes
rng = np.random.default_rng(0)
low = np.array([v[0] for p in [source.polygons[0]] for v in p.vertices])
times = np.array([v[1] for p in [source.polygons[0]] for v in p.vertices])
draws = 200_000
a = low.min() + (low.max() - low.min()) * rng.random(draws)
t = times.min() + (times.max() - times.min()) * rng.random(draws)
kept = np.array([region.contains(x, y) for x, y in zip(a, t)])
rejection = (tally(a, t) * kept).mean()

# drawn from the region, every weight corrected by the acceptance
direct = sampler.acceptance * tally(inverse_velocity, time).mean()

print(f'reject: {rejection:.6e}   direct: {direct:.6e}')
# --8<-- [end:unbiased]


def main():
    assert sampler.count == sum(len(p) - 2 for p in region.polygons)
    assert abs(sampler.acceptance - region.area / source.area) < 1e-15
    assert inverse_velocity.size == 200_000

    # every draw is inside; spot-check rather than test all 200k
    assert all(region.contains(float(x), float(y))
               for x, y in zip(inverse_velocity[:2000], time[:2000]))

    # the two estimators agree, and the direct one is far quieter
    error = (tally(a, t) * kept).std() / np.sqrt(draws)
    assert abs(rejection - direct) < 4 * error
    direct_error = sampler.acceptance * tally(inverse_velocity, time).std() / np.sqrt(draws)
    assert direct_error < error / 10


main()
