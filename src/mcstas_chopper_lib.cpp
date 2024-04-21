#include <sstream>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>

extern "C" {
#include <chopper-lib.h>
}

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_chopper_lib_impl, m) {
nb::class_<chopper_parameters>(m, "Chopper")
      .def(nb::init<double, double, double, double>(), "speed"_a=0, "phase"_a=0, "angle"_a=0, "distance"_a=0)
      .def_rw("speed", &chopper_parameters::speed, "Disk rotation speed in Hz")
      .def_rw("phase", &chopper_parameters::phase, "Disk rotation phase in degrees")
      .def_rw("angle", &chopper_parameters::angle, "Disk opening angle in degrees")
      .def_rw("distance", &chopper_parameters::path, "Source to disk path length distance in meters")
      .def("__str__", [](const chopper_parameters * c) {
            std::stringstream s;
            s << "Chopper(" << c->speed << " Hz, " << c->phase << " deg, " << c->angle << " deg, " << c->path << " m)";
            return s.str();
      })
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
      }, "choppers"_a, "inv_v_min"_a=1e-9, "inv_v_max"_a=1e3, "latest_emission"_a=0.003
);

m.def("inverse_velocity_limits",
      [](const std::vector<chopper_parameters> & choppers, const double inv_v_min, const double inv_v_max, const double latest_emission) {
            std::pair<double, double> out;
            auto no = chopper_inverse_velocity_limits(&out.first, &out.second, choppers.size(), choppers.data(), inv_v_min, inv_v_max, latest_emission);
            return std::make_tuple(no, out);
      }, "choppers"_a, "inv_v_min"_a=1e-9, "inv_v_max"_a=1e3, "latest_emission"_a=0.003
);

m.def("wavelength_limits",
      [](const std::vector<chopper_parameters> & choppers, const double lambda_min, const double lambda_max, const double latest_emission) {
      std::pair<double, double> out;
      auto no = chopper_wavelength_limits(&out.first, &out.second, choppers.size(), choppers.data(), lambda_min, lambda_max, latest_emission);
return std::make_tuple(no, out);
      }, "choppers"_a, "wavelength_min"_a=1e-4, "wavelength_max"_a=1e2, "latest_emission"_a=0.003
);

}