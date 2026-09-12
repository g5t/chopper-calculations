"""A guide is not a straight line, and what that does to the region."""
import chopcal
from chopcal.lib import Region

train = list(chopcal.bifrost(wavelength_max=3.0).values())
source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)

# --8<-- [start:spread]
straight = source.transmit(train)

# a tenth of a percent of extra path, accumulating along the flight path
spreads = [1e-3 * chopper.path for chopper in train]
spread = source.transmit(train, spreads)

print('straight:', straight.wavelength_ranges())
print('spread:  ', spread.wavelength_ranges())
print(f'area grew by {spread.area / straight.area - 1:.3%}')
# --8<-- [end:spread]

# --8<-- [start:one_sided]
(tight_low, tight_high), = straight.wavelength_ranges()
(loose_low, loose_high), = spread.wavelength_ranges()
print(f'slow edge moved {1e3 * (tight_low - loose_low):.4f} mAA')
print(f'fast edge moved {1e3 * (loose_high - tight_high):.4f} mAA')
# --8<-- [end:one_sided]


def main():
    assert spread.area > straight.area
    # only the slow edge moves: a longer route rescues a neutron arriving early and can
    # do nothing for one arriving late
    assert loose_low < tight_low
    assert abs(loose_high - tight_high) < 1e-12

    # a spread wide enough to overlap one turn of a disk with the next is refused
    try:
        source.transmit(train, [0.5 * chopper.path for chopper in train])
    except RuntimeError:
        pass
    else:                        # pragma: no cover
        raise AssertionError('an absurd spread should have been refused')


main()
