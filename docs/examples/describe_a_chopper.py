"""Describing a disk: slit edges, where the beam crosses, and how wide it is."""
# --8<-- [start:centred]
from chopcal.lib import Chopper

# a single opening 170 degrees across, centred on the beam
disk = Chopper.centred(speed=196.0, delay=3.878e-3, width=170.0, path=6.342)
print(list(disk.edges))          # [-85.0, 85.0]
print(disk.beam)                 # 0.0: the zero mark is on the beam at `delay`
# --8<-- [end:centred]

# --8<-- [start:edges]
# the disk's own numbers, unconverted: two edges per opening, increasing
multi = Chopper(speed=14.0, delay=0.0, beam=12.5,
                edges=[0.0, 10.0, 90.0, 105.0, 180.0, 195.0], path=8.53)
print(multi.openings)            # [(0.0, 10.0), (90.0, 105.0), (180.0, 195.0)]
print(multi.opening)             # 40.0 degrees of opening in total
# --8<-- [end:edges]

# --8<-- [start:across_the_mark]
# an opening straddling the zero mark is written past 360, never wrapped
straddling = Chopper(speed=14.0, edges=[350.0, 370.0], path=1.0)
print(straddling.opening)        # 20.0
# --8<-- [end:across_the_mark]

# --8<-- [start:aperture]
from chopcal.lib import beam_aperture

# a 60 x 90 mm beam through a 98.46 mm slit on a 350 mm disk
aperture = beam_aperture(radius=0.350, slit_height=0.09846,
                         window_width=0.060, window_height=0.090)
print(f'{aperture:.4f} degrees')
# --8<-- [end:aperture]

# --8<-- [start:parked]
shut = Chopper(speed=0.0, beam=180.0, edges=[-10.0, 10.0], path=20.0)
print(shut.parked_is_open())     # False: a beam stop
# --8<-- [end:parked]


def main():
    assert list(disk.edges) == [-85.0, 85.0]
    assert disk.opening == 170.0
    assert len(multi.openings) == 3
    assert abs(multi.opening - 40.0) < 1e-12
    assert straddling.opening == 20.0
    assert abs(aperture - 13.3796) < 1e-4
    assert not shut.parked_is_open()
    assert Chopper(speed=0.0, beam=0.0, edges=[-10.0, 10.0]).parked_is_open()

    # the two checks that catch a disk described wrongly
    for bad in ([0.0, 10.0, 20.0], [10.0, 0.0]):
        try:
            Chopper(speed=14.0, edges=bad, path=1.0)
        except ValueError:
            pass
        else:                    # pragma: no cover
            raise AssertionError(f'{bad} should have been refused')


main()
