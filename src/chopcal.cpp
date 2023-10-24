//
// Created by gst on 24/10/23.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>

#include "choppers.h"

NB_MODULE(_chopcal_impl, m) {
m.def("bifrost", &bifrost);
}