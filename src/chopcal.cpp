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
m.def("bifrost", &bifrost, "e_max"_a, "lambda_min"_a=0, "shaping_time"_a=0.0002,
      "Calculate chopper parameters for the BIFROST spectrometer. Each of (phase, delay, phase) for the six choppers.");
}