"""Every constant the chopper calculations use, checked against an independent source.

These used to be magic numbers, and the same quantity was written more than one way. The
neutron h^2/2m appeared as 81.82 on one line and as 0.1106 -- its inverse square root --
on the next, and 0.1106 * sqrt(81.82) is 1.000426, not 1. So asking for a 3 angstrom band
produced settings for 2.99872 angstrom, and no test could see it because the expected
values were derived from the same wrong numbers.

They are derived from the SI primitives now, in one header, so a disagreement of that kind
cannot be written. What is left is checking that the primitives themselves are right, and
scipp is the independent source for that.

scipp is deliberately absent from `test-requires` in pyproject.toml -- cibuildwheel would
have to build it and scipy for every wheel -- so this module is skipped when building
wheels and runs for a developer. The literals in test_choplib.get_factors are the branch
CI exercises instead, and one test here holds them to the same standard.
"""
import pytest

import chopcal

scipp_constants = pytest.importorskip(
    'scipp.constants', reason='scipp is not installed; it is not in test-requires')

TO_LAST_BIT = 1e-15


def test_the_primitives_are_the_defined_si_values():
    """Two of the three are exact by definition, so these are equalities, not tolerances."""
    from scipp.constants import Planck, elementary_charge, neutron_mass

    assert chopcal.constants.PLANCK == Planck.to(unit='J s').value
    assert chopcal.constants.ELEMENTARY_CHARGE == elementary_charge.to(unit='C').value
    assert chopcal.constants.NEUTRON_MASS == pytest.approx(
        neutron_mass.to(unit='kg').value, rel=TO_LAST_BIT)


def test_h_over_m_converts_wavelength_to_speed():
    """A neutron of wavelength L angstrom travels H_OVER_M / L metres per second."""
    from scipp.constants import Planck, neutron_mass

    expected = (Planck / neutron_mass).to(unit='angstrom m/s').value
    assert chopcal.constants.H_OVER_M == pytest.approx(expected, rel=TO_LAST_BIT)


def test_h2_over_2m_converts_wavelength_to_energy():
    """A neutron of wavelength L angstrom has energy H2_OVER_2M / (L * L) meV."""
    from scipp.constants import Planck, neutron_mass

    expected = (Planck * Planck / 2 / neutron_mass).to(
        unit='millielectronvolt angstrom * angstrom').value
    assert chopcal.constants.H2_OVER_2M == pytest.approx(expected, rel=TO_LAST_BIT)


def test_v2k_and_k2v_convert_speed_to_wavevector():
    from scipp.constants import hbar, neutron_mass

    expected = (neutron_mass / hbar).to(unit='1/angstrom/(m/s)').value
    assert chopcal.constants.V2K == pytest.approx(expected, rel=TO_LAST_BIT)
    assert chopcal.constants.K2V == pytest.approx(1 / expected, rel=TO_LAST_BIT)


def test_the_two_wavelength_conversions_agree():
    """The bug that motivated all of this.

    H2_OVER_2M turns a wavelength into an energy and back, and H_OVER_M turns it into a
    speed. Both are h over a neutron mass, so they have to be consistent -- 0.1106 and
    81.82 were not.
    """
    from math import sqrt

    c = chopcal.constants
    # E = H2_OVER_2M / L^2 and v = H_OVER_M / L, so E = m v^2 / 2 must come back out
    for wavelength in (0.5, 1.0, 3.0, 10.0):
        energy = c.H2_OVER_2M / (wavelength * wavelength)
        assert sqrt(c.H2_OVER_2M / energy) == pytest.approx(wavelength, rel=TO_LAST_BIT)

    # and K2V * 2 pi is the same h/m that H_OVER_M is
    assert c.K2V * 2 * c.PI == pytest.approx(c.H_OVER_M, rel=TO_LAST_BIT)


def test_pi_is_pi():
    """Spelled out rather than taken from <numbers>.

    That header is C++20 *library*, and the musllinux images cibuildwheel targets ship a
    libstdc++ without it even though the compiler takes C++20 language features -- so
    including it broke the cp311 and cp312 musllinux wheels. The literal carries more
    digits than a double holds, so it is bit-identical to std::numbers::pi_v<double>.
    """
    from math import pi

    assert chopcal.constants.PI == pi
    # and the conversions built on it stay consistent with each other
    assert chopcal.constants.V2K * chopcal.constants.K2V == pytest.approx(1.0, rel=TO_LAST_BIT)


def test_the_ess_values_are_the_ones_mcstas_simulates_with():
    """From mcstas-comps/share/ESS_butterfly-lib.h.

    The calculations used to round the pulse duration to 2.86e-3, which disagreed with
    the source actually being simulated.
    """
    assert chopcal.constants.SOURCE_DURATION == 2.857e-3
    assert chopcal.constants.SOURCE_FREQUENCY == 14.0


def test_the_no_scipp_fallbacks_are_not_stale():
    """test_choplib falls back to literals when scipp is missing, which is what CI does.

    They were frozen against an older CODATA revision and had drifted by 1.5e-9. Nothing
    would have noticed, because the branch that uses them is the branch that runs where
    scipp is not installed.
    """
    from test_choplib import get_factors

    h_over_m, h2_over_2m = get_factors()
    assert h_over_m == pytest.approx(chopcal.constants.H_OVER_M, rel=TO_LAST_BIT)
    assert h2_over_2m == pytest.approx(chopcal.constants.H2_OVER_2M, rel=TO_LAST_BIT)


def test_the_geometry_is_named_but_unchanged():
    """This change names BIFROST's geometry; it does not restate it.

    Other tests assert some of these by exact equality, so a typo here would surface as a
    puzzling failure elsewhere.
    """
    c = chopcal.constants
    assert c.INSTRUMENT_LENGTH == 162.0
    assert c.PULSE_SHAPING_DISTANCE == 4.41 + 0.032 + 2.0 - 0.1
    assert c.FRAME_OVERLAP_1_DISTANCE == 8.530
    assert c.FRAME_OVERLAP_2_DISTANCE == 14.973
    assert c.BANDWIDTH_DISTANCE == 78.0
    assert (c.PULSE_SHAPING_ANGLE, c.FRAME_OVERLAP_1_ANGLE,
            c.FRAME_OVERLAP_2_ANGLE, c.BANDWIDTH_ANGLE) == (170.0, 38.26, 52.01, 161.0)
