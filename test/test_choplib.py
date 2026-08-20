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


def do_test_choplib():
    import chopcal as cc
    lambdas = (steps(0.5, 1.75, 0.05) + steps(1.76, 1.9, 0.001) + steps(2., 5, 0.1))
    settings = [cc.bifrost(0., x) for x in lambdas]
    calc = [cc.lib.wavelength_limits(list(setting.values())) for setting in settings]

    assert all(c[0] == 1 for c in calc)
    minimum = [c[1][0] for c in calc]
    maximum = [c[1][1] for c in calc]

    assert all(x - y < 0.025 for x, y in zip(lambdas, maximum))

    h_m, h2_2m = get_factors()
    # the same geometry the calculation itself uses, rather than a third copy of it
    from chopcal import constants
    active_length = constants.INSTRUMENT_LENGTH - constants.PULSE_SHAPING_DISTANCE
    maximum_bandwidth = 1 / active_length / constants.SOURCE_FREQUENCY * h_m  # ~= 1.815, not 1.77:
    # the band is cut by the 161 degree bandwidth opening, not by this nominal width
    expected_bandwidth = 1.77  # why isn't this 1.815??

    difference = [x - y for x, y in zip(maximum, minimum)]
    over = [d for d, l in zip(difference, lambdas) if l > maximum_bandwidth]
    assert all(abs(x - expected_bandwidth) < 0.01 for x in over)
    under = [(d, l) for d, l in zip(difference, lambdas) if l < maximum_bandwidth]
    assert all(abs(x - y + maximum_bandwidth - expected_bandwidth) < 0.05 for x, y in under)

    return lambdas, minimum, maximum, expected_bandwidth


@pytest.mark.skipif(system().lower().startswith('win'), reason="Test tolerance too tight for windows")
def test_choplib():
    # hide the return values from pytest
    do_test_choplib()


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
