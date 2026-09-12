import pytest
from platform import system


def get_factors():
    from importlib.util import find_spec
    if find_spec('scipp') and find_spec('scipy'):
        from scipp.constants import Planck, neutron_mass
        h_over_m = (Planck / neutron_mass).to(unit='angstrom m/s').value
        h2_over_2m = (Planck * Planck / 2 / neutron_mass).to(unit='millielectronvolt angstrom * angstrom').value
    else:
        # No scipp here -- which is the case when cibuildwheel builds a wheel, since
        # test-requires deliberately leaves out scipp and scipy. These are the same
        # values chopcal was compiled with; test_constants holds them to that.
        from chopcal import constants
        h_over_m, h2_over_2m = constants.H_OVER_M, constants.H2_OVER_2M
    return h_over_m, h2_over_2m


def steps(start, stop, step):
    """The points numpy.arange would give, without needing numpy to give them."""
    from math import ceil
    return [start + i * step for i in range(max(0, ceil((stop - start) / step)))]


def do_test_choplib(apertures=True):
    """Where the band lands, over the whole range of settings.

    Two invariants, both tighter than the tolerances this used to carry:

    * once the requested wavelength is past the bandwidth, the band saturates at exactly
      that width -- 1.91 AA for the real beam, 1.77 AA for a pencil one;
    * the top of the band tracks the request with a fixed offset, which is the only place
      the beam width shows up in where the band *sits* rather than how wide it is.
    """
    import chopcal as cc
    lambdas = (steps(0.5, 1.75, 0.05) + steps(1.76, 1.9, 0.001) + steps(2., 5, 0.1))
    settings = [cc.bifrost(0., x, apertures=apertures) for x in lambdas]
    # The first window, not the envelope. With the beam at its real width the train also
    # leaks a ~0.03 AA sliver near 38 AA -- an artefact of how chopper-lib intersects, not
    # something that gets through: see MaskContradictsTheWindowsTestCase in test_mask.py.
    calc = [cc.lib.wavelength_windows(list(setting.values())) for setting in settings]
    assert all(c for c in calc), 'every setting should pass a band'

    minimum = [c[0][0] for c in calc]
    maximum = [c[0][1] for c in calc]
    width = [b - a for a, b in zip(minimum, maximum)]

    expected_bandwidth = (cc.BIFROST_BANDWIDTH if apertures
                          else cc.BIFROST_BANDWIDTH_POINT_BEAM)

    # the nominal width from the geometry alone, which the 161 degree bandwidth opening
    # cuts into: ~1.815, between the two figures above rather than above both
    from chopcal import constants
    h_m, _ = get_factors()
    active_length = constants.INSTRUMENT_LENGTH - constants.PULSE_SHAPING_DISTANCE
    nominal = 1 / active_length / constants.SOURCE_FREQUENCY * h_m
    assert 1.7 < nominal < 1.9

    # saturated: far enough past the bandwidth that the band is not clipped at the bottom
    saturated = [w for w, l in zip(width, lambdas) if l > expected_bandwidth + 0.1]
    assert saturated, 'no saturated settings to check'
    assert all(abs(w - expected_bandwidth) < 0.01 for w in saturated)

    # the top edge tracks the request at a fixed offset
    offsets = [b - l for b, l in zip(maximum, lambdas)]
    assert max(offsets) - min(offsets) < 0.005, 'the offset should not depend on the band'
    if apertures:
        assert all(0.04 < o < 0.07 for o in offsets), 'a real beam spills past the request'
    else:
        assert all(-0.03 < o < 0.0 for o in offsets), 'a pencil beam falls just short'

    return lambdas, minimum, maximum, expected_bandwidth


@pytest.mark.skipif(system().lower().startswith('win'), reason="Test tolerance too tight for windows")
def test_choplib():
    # hide the return values from pytest
    do_test_choplib()


@pytest.mark.skipif(system().lower().startswith('win'), reason="Test tolerance too tight for windows")
def test_choplib_point_beam():
    """The same, for the pencil beam every release before 0.6 described."""
    do_test_choplib(apertures=False)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    lambdas, minimum, maximum, expected_bandwidth = do_test_choplib()
    plt.plot(lambdas, minimum, label='minimum')
    plt.plot(lambdas, maximum, label='maximum')
    plt.plot(lambdas, lambdas, label='target')
    plt.plot(lambdas, [x - y for x, y in zip(maximum, minimum)], label='bandwidth')
    plt.plot(lambdas, [expected_bandwidth] * len(lambdas), label='expected bandwidth')
    plt.setp(plt.gca(), 'xlabel', r'set $\lambda/\mathrm{\AA}$', 'ylabel', r'$\lambda/\mathrm{\AA}$')
    plt.legend()
    plt.show()
