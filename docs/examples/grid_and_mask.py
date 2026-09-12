"""The grid answer: a mask, its transmission, and how it compares with the region."""
import numpy as np

import chopcal
from chopcal.lib import Region, inverse_velocity_time_mask, wavelength_to_inverse_velocity

train = list(chopcal.bifrost(wavelength_max=3.0).values())

# --8<-- [start:mask]
inverse_velocities = np.linspace(wavelength_to_inverse_velocity(0.4),
                                 wavelength_to_inverse_velocity(45.0), 2001)
times = np.linspace(0.0, 3e-3, 1501)

mask, allowed = inverse_velocity_time_mask(train, inverse_velocities, times)
print(mask.shape, 'time by inverse velocity')     # note the order
print(f'{allowed} of {mask.size} bins allowed')
# --8<-- [end:mask]

# --8<-- [start:probability]
from chopcal.lib import unmasked_probability

signal = np.ones(mask.shape)
print(f'flat signal transmits {unmasked_probability(signal, mask):.6e}')
# --8<-- [end:probability]

# --8<-- [start:versus]
from chopcal.lib import MaskSampler

grid = MaskSampler(mask, inverse_velocities, times,
                   inverse_velocities[0],
                   inverse_velocities[-1] - inverse_velocities[0],
                   times[0], times[-1] - times[0])

source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)
exact = source.transmit(train).area / source.area

print(f'grid    {grid.acceptance:.6e} from {grid.count} cells')
print(f'region  {exact:.6e} exactly')
print(f'grid is {grid.acceptance / exact - 1:.2%} high')
# --8<-- [end:versus]


def main():
    assert mask.shape == (times.size - 1, inverse_velocities.size - 1)
    assert allowed > 0
    # a flat signal transmits the allowed fraction of cells
    assert abs(unmasked_probability(signal, mask) - allowed / mask.size) < 1e-12
    # a grid can only over-accept: a partly covered cell is kept whole
    assert grid.acceptance > exact


main()
