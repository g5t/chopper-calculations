//
// Created by gst on 24/10/23.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>

#include "choppers.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_chopcal_impl, m) {
m.def("bifrost", &bifrost, "energy_min"_a=0, "wavelength_max"_a=0, "shaping_time"_a=0.0002,
      "Chopper settings for BIFROST, as a speed, delay, opening angle and flight path\n"
      "for each of its six choppers.\n\n"
      "The band is a fixed ~1.77 angstrom wide, so one number places it: either the\n"
      "longest wavelength (angstrom) or the lowest energy (meV) to pass, the band\n"
      "running from there to ~1.77 angstrom shorter. A positive wavelength_max wins\n"
      "over energy_min; giving neither leaves every delay infinite.\n\n"
      "Prefer chopcal.bifrost, which rejects that case and returns a ChopperSet.");
}