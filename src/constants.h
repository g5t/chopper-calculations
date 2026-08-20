//
// Every number these calculations depend on, in one place, with its provenance.
//
// The chopper calculations used to write the same physical constant several ways -- the
// neutron h^2/2m appeared as 81.82 in one line and as 0.1106 (its inverse square root) in
// the next, and 0.1106 * sqrt(81.82) is 1.000426, not 1. Asking for a 3 angstrom band
// therefore produced settings for 2.99872 angstrom. Deriving everything from the SI
// primitives makes that class of disagreement impossible to write.
//
#ifndef CHOPCAL_CONSTANTS_H
#define CHOPCAL_CONSTANTS_H

#include <numbers>

namespace chopcal::constants {

// -- SI 2019 and CODATA 2022 -------------------------------------------------
// The first two are exact by definition of the SI; the neutron mass is measured.

/// Planck constant, J s. Exact: the SI defines the kilogram through it.
inline constexpr double PLANCK = 6.62607015e-34;
/// Elementary charge, C. Exact: the SI defines the ampere through it.
inline constexpr double ELEMENTARY_CHARGE = 1.602176634e-19;
/// Neutron rest mass, kg. CODATA 2022.
inline constexpr double NEUTRON_MASS = 1.67492750056e-27;

inline constexpr double PI = std::numbers::pi_v<double>;

// -- derived, in the units the chopper calculations work in ------------------
// Each is the SI expression scaled once into angstrom/meV, so the value and the way it
// was reached are the same line of code. They agree with scipp.constants to the last
// bit, which test/test_constants.py asserts.

/// h/m, angstrom m/s. A neutron of wavelength L travels H_OVER_M / L metres per second.
inline constexpr double H_OVER_M = PLANCK / NEUTRON_MASS * 1e10;

/// h^2/2m, meV angstrom^2. A neutron of wavelength L has energy H2_OVER_2M / (L * L),
/// and a neutron of energy E has wavelength sqrt(H2_OVER_2M / E).
inline constexpr double H2_OVER_2M =
    PLANCK * PLANCK / (2 * NEUTRON_MASS) / (ELEMENTARY_CHARGE * 1e-3) * 1e20;

/// Speed to wavevector, (1/angstrom) per (m/s). McStas spells this V2K.
inline constexpr double V2K = 2 * PI * NEUTRON_MASS / PLANCK * 1e-10;
/// Wavevector to speed, (m/s) per (1/angstrom). McStas spells this K2V.
inline constexpr double K2V = 1.0 / V2K;

// -- the ESS source ----------------------------------------------------------
// From mcstas-comps/share/ESS_butterfly-lib.h, which is what the simulated source
// itself uses. The calculations previously rounded the duration to 2.86e-3.

/// Length of the high-flux plateau, s. McStas: ESS_SOURCE_DURATION.
inline constexpr double SOURCE_DURATION = 2.857e-3;
/// Source repetition rate, Hz. McStas: ESS_SOURCE_FREQUENCY.
inline constexpr double SOURCE_FREQUENCY = 14.0;

/// Time from the proton pulse to the high-flux plateau, s. A BIFROST working figure
/// rather than a published one; it has no McStas counterpart.
inline constexpr double PULSE_HIGH_FLUX_OFFSET = 2.0e-4;

// -- BIFROST geometry --------------------------------------------------------
// Values unchanged: tests assert some of these by exact equality, and this section is
// naming rather than correcting.

/// Moderator to sample, m.
inline constexpr double INSTRUMENT_LENGTH = 162.0;
/// Moderator to the first pulse-shaping disk, m. Left as its four surveyed parts.
inline constexpr double PULSE_SHAPING_DISTANCE = 4.41 + 0.032 + 2.0 - 0.1;
/// Moderator to the first and second frame-overlap disks, m.
inline constexpr double FRAME_OVERLAP_1_DISTANCE = 8.530;
inline constexpr double FRAME_OVERLAP_2_DISTANCE = 14.973;
/// Moderator to the first bandwidth disk, m.
inline constexpr double BANDWIDTH_DISTANCE = 78.0;
/// Axial gap between the two disks of a counter-rotating pair, m.
inline constexpr double PAIR_SEPARATION = 0.02;

/// Opening angles, degrees.
inline constexpr double PULSE_SHAPING_ANGLE = 170.0;
inline constexpr double FRAME_OVERLAP_1_ANGLE = 38.26;
inline constexpr double FRAME_OVERLAP_2_ANGLE = 52.01;
inline constexpr double BANDWIDTH_ANGLE = 161.0;

/// Degrees in a turn -- a unit conversion, not a measurement.
inline constexpr double DEGREES_PER_TURN = 360.0;

}  // namespace chopcal::constants

#endif  // CHOPCAL_CONSTANTS_H
