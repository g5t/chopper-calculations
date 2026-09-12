"""Chopper settings for BIFROST, and the band they pass."""
# --8<-- [start:settings]
import chopcal

settings = chopcal.bifrost(wavelength_max=3.0)
print(settings)
# --8<-- [end:settings]

# --8<-- [start:band]
from chopcal.lib import wavelength_windows

bands = wavelength_windows(list(settings.values()))
low, high = bands[0]
print(f'{low:.4f} to {high:.4f} AA, {high - low:.4f} AA wide')
# --8<-- [end:band]

# --8<-- [start:one_chopper]
bandwidth = settings['bw1']
print(bandwidth)                 # said out loud
print(bandwidth.delay)           # seconds, not milliseconds
print(list(bandwidth.edges))     # degrees from the disk's zero mark
# --8<-- [end:one_chopper]


def main():
    assert set(settings) == {'ps1', 'ps2', 'fo1', 'fo2', 'bw1', 'bw2'}
    # the band is placed where it was asked for, to within the beam's own width
    assert abs(high - 3.0) < 0.06
    assert abs((high - low) - chopcal.BIFROST_BANDWIDTH) < 0.01
    # the settings are a dict, whatever else they print like
    assert settings['bw1'].speed == -settings['bw2'].speed


main()
