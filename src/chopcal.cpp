//
// Created by gst on 24/10/23.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>

#include "choppers.h"
#include "constants.h"

namespace nb = nanobind;
using namespace nb::literals;

/// Expose the constants so they can be checked against an independent source.
/// test/test_constants.py compares every one of these to scipp.constants.
static void bind_constants(nb::module_ &m) {
    namespace c = chopcal::constants;
    auto constants = m.def_submodule(
        "constants",
        "Physical and facility constants, and where each of them comes from.\n\n"
        "The neutron quantities are derived from the SI/CODATA primitives rather than\n"
        "tabulated, so they cannot disagree with each other -- which they used to.");
    constants.attr("PLANCK") = c::PLANCK;
    constants.attr("ELEMENTARY_CHARGE") = c::ELEMENTARY_CHARGE;
    constants.attr("NEUTRON_MASS") = c::NEUTRON_MASS;
    constants.attr("PI") = c::PI;
    constants.attr("H_OVER_M") = c::H_OVER_M;
    constants.attr("H2_OVER_2M") = c::H2_OVER_2M;
    constants.attr("V2K") = c::V2K;
    constants.attr("K2V") = c::K2V;
    constants.attr("SOURCE_DURATION") = c::SOURCE_DURATION;
    constants.attr("SOURCE_FREQUENCY") = c::SOURCE_FREQUENCY;
    constants.attr("PULSE_HIGH_FLUX_OFFSET") = c::PULSE_HIGH_FLUX_OFFSET;
    constants.attr("INSTRUMENT_LENGTH") = c::INSTRUMENT_LENGTH;
    constants.attr("PULSE_SHAPING_DISTANCE") = c::PULSE_SHAPING_DISTANCE;
    constants.attr("FRAME_OVERLAP_1_DISTANCE") = c::FRAME_OVERLAP_1_DISTANCE;
    constants.attr("FRAME_OVERLAP_2_DISTANCE") = c::FRAME_OVERLAP_2_DISTANCE;
    constants.attr("BANDWIDTH_DISTANCE") = c::BANDWIDTH_DISTANCE;
    constants.attr("PAIR_SEPARATION") = c::PAIR_SEPARATION;
    constants.attr("PULSE_SHAPING_ANGLE") = c::PULSE_SHAPING_ANGLE;
    constants.attr("FRAME_OVERLAP_1_ANGLE") = c::FRAME_OVERLAP_1_ANGLE;
    constants.attr("FRAME_OVERLAP_2_ANGLE") = c::FRAME_OVERLAP_2_ANGLE;
    constants.attr("BANDWIDTH_ANGLE") = c::BANDWIDTH_ANGLE;
    constants.attr("DEGREES_PER_TURN") = c::DEGREES_PER_TURN;
}

NB_MODULE(_chopcal_impl, m) {
bind_constants(m);
m.def("bifrost", &bifrost, "energy_min"_a=0, "wavelength_max"_a=0, "shaping_time"_a=0.0002,
      "Chopper settings for BIFROST, as a speed, delay, opening angle and flight path\n"
      "for each of its six choppers.\n\n"
      "The band is a fixed ~1.77 angstrom wide, so one number places it: either the\n"
      "longest wavelength (angstrom) or the lowest energy (meV) to pass, the band\n"
      "running from there to ~1.77 angstrom shorter. A positive wavelength_max wins\n"
      "over energy_min; giving neither leaves every delay infinite.\n\n"
      "Prefer chopcal.bifrost, which rejects that case and returns a ChopperSet.");
}