#include <sstream>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>

extern "C" {
#include <chopper-lib.h>
}

// Use the McStas defines if possible, or define them ourselves
#ifndef V2K
#define V2K 1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#endif
#ifndef K2V
#define K2V 629.622368        /* Convert k[1/AA] to v[m/s] */
#endif
#ifndef PI
#define PI 3.14159265358979323846
#endif

#define LAMBDA_MIN 1e-4
#define LAMBDA_MAX 1e2
#define L2INVV(L) (L / 2 / PI / K2V)
#define INVERSE_V_MIN L2INVV(LAMBDA_MIN)
#define INVERSE_V_MAX L2INVV(LAMBDA_MAX)
#define LATEST_EMISSION 0.003

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_chopper_lib_impl, m) {
nb::class_<chopper_parameters>(m, "Chopper",
      "One disk chopper: how fast it turns, when it is open, how wide, and how far away.\n\n"
      "An opening is centred on the beam at `delay` and again every 1/`speed` after\n"
      "that, staying there for `angle`/360/|`speed`| seconds each time. A negative\n"
      "`speed` turns the disk the other way, which changes none of that -- a delay is a\n"
      "time, so it means the same thing in either direction.")
      .def(nb::init<double, double, double, double>(), "speed"_a=0, "delay"_a=0, "angle"_a=0, "path"_a=0)
      .def_rw("speed", &chopper_parameters::speed,
              "Rotation frequency in Hz; negative turns the disk the other way")
      .def_rw("delay", &chopper_parameters::delay,
              "When an opening is centred on the beam, in seconds")
      .def_rw("angle", &chopper_parameters::angle, "Angular width of one opening, in degrees")
      .def_rw("path", &chopper_parameters::path, "Flight path from the source to the disk, in metres")
      ;

m.def("inverse_velocity_windows",
      [](const std::vector<chopper_parameters> & choppers, const double inv_v_min, const double inv_v_max, const double latest_emission) {
            const auto size = static_cast<unsigned>(choppers.size());
            const auto * data = choppers.data();
            const auto rs = chopper_inverse_velocity_windows(size, data, inv_v_min, inv_v_max, latest_emission);
            std::vector<std::pair<double, double>> out;
            out.reserve(rs.count);
            for (int i=0; i<rs.count; ++i) out.emplace_back(rs.ranges[i].minimum, rs.ranges[i].maximum);
            return out;
      }, "choppers"_a, "inv_v_min"_a=INVERSE_V_MIN, "inv_v_max"_a=INVERSE_V_MAX, "latest_emission"_a=LATEST_EMISSION
);

m.def("inverse_velocity_limits",
      [](const std::vector<chopper_parameters> & choppers, const double inv_v_min, const double inv_v_max, const double latest_emission) {
            std::pair<double, double> out;
            auto no = chopper_inverse_velocity_limits(&out.first, &out.second, choppers.size(), choppers.data(), inv_v_min, inv_v_max, latest_emission);
            return std::make_tuple(no, out);
      }, "choppers"_a, "inv_v_min"_a=INVERSE_V_MIN, "inv_v_max"_a=INVERSE_V_MAX, "latest_emission"_a=LATEST_EMISSION
);

m.def("wavelength_limits",
      [](const std::vector<chopper_parameters> & choppers, const double lambda_min, const double lambda_max, const double latest_emission) {
      std::pair<double, double> out;
      auto no = chopper_wavelength_limits(&out.first, &out.second, choppers.size(), choppers.data(), lambda_min, lambda_max, latest_emission);
return std::make_tuple(no, out);
      }, "choppers"_a, "wavelength_min"_a=LAMBDA_MIN, "wavelength_max"_a=LAMBDA_MAX, "latest_emission"_a=LATEST_EMISSION
);

}