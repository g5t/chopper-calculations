"""What a train passes, and where the quick answer is loose."""
import chopcal
from chopcal.lib import wavelength_limits, wavelength_windows

train = list(chopcal.bifrost(wavelength_max=3.0).values())

# --8<-- [start:windows]
bands = wavelength_windows(train)
for low, high in bands:
    print(f'{low:.4f} to {high:.4f} AA')
# --8<-- [end:windows]

# --8<-- [start:limits]
count, (low, high) = wavelength_limits(train)
print(f'{count} band(s), envelope {low:.4f} to {high:.4f} AA')
# --8<-- [end:limits]

# --8<-- [start:region]
from chopcal.lib import Region

source = Region.from_wavelengths(0.4, 45.0, 0.0, 3e-3)
region = source.transmit(train)
print(region.wavelength_ranges())
# --8<-- [end:region]


def main():
    # The quick answer reports a band near 38 AA that nothing passes, so the envelope
    # spans a gap the train does not deliver.
    assert len(bands) == 2
    assert 37.9 < bands[1][0] < 38.1
    assert count == 2
    assert high > 37.0

    # The region does not have it, and its band is inside the reported one.
    exact = region.wavelength_ranges()
    assert len(exact) == 1
    assert bands[0][0] <= exact[0][0]
    assert bands[0][1] >= exact[0][1] - 1e-12


main()
