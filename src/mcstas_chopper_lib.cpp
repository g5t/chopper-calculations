#include <algorithm>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>

#include "chopper.h"
#include "constants.h"

// These fell back to McStas' own macros when the file was built inside an instrument,
// which never happens: this is a nanobind extension, so the literals were always what
// compiled. Since chopper-lib 4.0.0 the C file defines nothing of its own and refuses to
// compile without V2K, K2V and PI -- CMake passes the McStas runtime's values in, which
// is what an embedded build would have used. What follows only sets the default search
// bounds below, and is deliberately the SI/CODATA-derived set instead.
using chopcal::constants::K2V;
using chopcal::constants::PI;

/// Wavelength bounds searched when a caller does not give their own, angstrom.
inline constexpr double LAMBDA_MIN = 1e-4;
inline constexpr double LAMBDA_MAX = 1e2;
inline constexpr double INVERSE_V_MIN = LAMBDA_MIN / 2 / PI / K2V;
inline constexpr double INVERSE_V_MAX = LAMBDA_MAX / 2 / PI / K2V;
/// How long after t=0 a neutron may still be emitted, s.
inline constexpr double LATEST_EMISSION = 0.003;

// The library's own conversion numbers, passed in by CMake from the same
// CHOPPER_LIB_DEFINITIONS the C file is compiled with, under prefixed names so they do
// not collide with the chopcal::constants variables called PI, V2K and K2V. Converting
// with anything else would put `wavelength_windows` and `wavelength_limits` -- which the
// C converts for itself -- into disagreement.
#if !defined(CHOPPER_LIB_V2K) || !defined(CHOPPER_LIB_K2V) || !defined(CHOPPER_LIB_PI)
#error "Build with CHOPCAL_RUNTIME_DEFINITIONS; see CMakeLists.txt"
#endif
/// Wavelength in angstrom to inverse velocity in s/m, as chopper_wavelength_limits does it.
inline constexpr double wavelength_to_inverse_velocity(const double lambda) {
  return lambda * CHOPPER_LIB_V2K / 2 / CHOPPER_LIB_PI;
}
inline constexpr double inverse_velocity_to_wavelength(const double inverse_velocity) {
  return inverse_velocity * CHOPPER_LIB_K2V * 2 * CHOPPER_LIB_PI;
}

namespace nb = nanobind;
using namespace nb::literals;
using chopcal::Chopper;

namespace {

using f64_1d = nb::ndarray<const double, nb::ndim<1>, nb::c_contig>;
using f64_2d = nb::ndarray<const double, nb::ndim<2>, nb::c_contig>;
using i32_2d = nb::ndarray<const int, nb::ndim<2>, nb::c_contig>;
using mask_out = nb::ndarray<nb::numpy, int, nb::ndim<2>>;

/// Bin edges, checked. The C exits the process on a count mismatch, so it never sees one.
unsigned bin_count(const f64_1d & edges, const char * name) {
  if (edges.size() < 2) {
    throw std::invalid_argument(
        std::string(name) + " needs at least two bin edges to make one bin; got "
        + std::to_string(edges.size()));
  }
  const auto * v = edges.data();
  for (size_t i = 1; i < edges.size(); ++i) {
    if (!(v[i] > v[i - 1])) {
      throw std::invalid_argument(std::string(name) + " bin edges must strictly increase");
    }
  }
  return static_cast<unsigned>(edges.size() - 1);
}

/// A freshly allocated 2-D int array Python owns.
mask_out make_mask_array(int * data, const size_t rows, const size_t columns) {
  nb::capsule owner(data, [](void * p) noexcept { delete[] static_cast<int *>(p); });
  const size_t shape[2] = {rows, columns};
  return mask_out(data, 2, shape, owner);
}

/** chopper-lib's mask sampler, owning what `chopper_mask_sampler_make` allocated.
 *
 * The C hands back six raw pointers and a `_free`; this ties them to a Python object's
 * lifetime, which is the whole reason the sampler cannot be exposed as a plain struct.
 */
class MaskSampler {
  chopper_mask_sampler sampler_{};

 public:
  MaskSampler(const i32_2d & mask, const f64_1d & inverse_velocities, const f64_1d & times,
              const double inverse_velocity_minimum, const double inverse_velocity_range,
              const double time_minimum, const double time_range) {
    chopper_mask_sampler_empty(&sampler_);
    const auto time_bins = bin_count(times, "times");
    const auto inverse_velocity_bins = bin_count(inverse_velocities, "inverse_velocities");
    if (mask.shape(0) != time_bins || mask.shape(1) != inverse_velocity_bins) {
      throw std::invalid_argument(
          "A mask over these edges is (" + std::to_string(time_bins) + ", "
          + std::to_string(inverse_velocity_bins) + "); got ("
          + std::to_string(mask.shape(0)) + ", " + std::to_string(mask.shape(1)) + ")");
    }
    sampler_ = chopper_mask_sampler_make(
        mask.data(), inverse_velocity_bins, time_bins,
        inverse_velocities.data(), times.data(),
        inverse_velocity_minimum, inverse_velocity_range, time_minimum, time_range);
  }

  MaskSampler(const MaskSampler &) = delete;
  MaskSampler & operator=(const MaskSampler &) = delete;

  ~MaskSampler() { chopper_mask_sampler_free(&sampler_); }

  [[nodiscard]] unsigned count() const { return sampler_.count; }
  [[nodiscard]] double acceptance() const { return sampler_.acceptance; }

  void require_cells() const {
    if (sampler_.count == 0) {
      throw std::runtime_error(
          "This sampler has no cells to draw from: the mask allows nothing inside the "
          "sampled region. Check `count` before drawing.");
    }
  }

  [[nodiscard]] std::pair<double, double> draw(const double cell_deviate,
                                               const double inverse_velocity_deviate,
                                               const double time_deviate) const {
    require_cells();
    std::pair<double, double> out;
    chopper_mask_sampler_draw(&sampler_, cell_deviate, inverse_velocity_deviate,
                              time_deviate, &out.first, &out.second);
    return out;
  }

  /// The same draw over three arrays of deviates, so the caller keeps its own generator.
  [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> draw_many(
      const f64_1d & cell_deviates, const f64_1d & inverse_velocity_deviates,
      const f64_1d & time_deviates) const {
    require_cells();
    const auto n = cell_deviates.size();
    if (inverse_velocity_deviates.size() != n || time_deviates.size() != n) {
      throw std::invalid_argument("The three deviate arrays must be the same length");
    }
    std::vector<double> inverse_velocity(n), time(n);
    for (size_t i = 0; i < n; ++i) {
      chopper_mask_sampler_draw(&sampler_, cell_deviates.data()[i],
                                inverse_velocity_deviates.data()[i],
                                time_deviates.data()[i], &inverse_velocity[i], &time[i]);
    }
    return {std::move(inverse_velocity), std::move(time)};
  }
};

/** chopper-lib's transmitted region, owning the polygons it holds.
 *
 * `chopper_polygon_set` is a count and a pointer with a `_free`; this ties them to a
 * Python object's lifetime. A `chopper_polygon` needs no such treatment -- it is a
 * fixed-size value -- so it is bound as it stands.
 */
class Region {
  chopper_polygon_set set_;

  void append(const chopper_polygon_set & other) {
    for (unsigned i = 0; i < other.count; ++i) {
      if (!chopper_polygon_set_add(&set_, &other.polygon[i])) {
        throw std::runtime_error("Out of memory copying a transmitted region");
      }
    }
  }

 public:
  Region() : set_(chopper_polygon_set_empty()) {}
  Region(const Region & other) : set_(chopper_polygon_set_empty()) { append(other.set_); }
  Region & operator=(const Region & other) {
    if (this != &other) {
      chopper_polygon_set_free(&set_);
      append(other.set_);
    }
    return *this;
  }
  ~Region() { chopper_polygon_set_free(&set_); }

  /// The rectangle a source draws from, which is where every region starts.
  static Region rectangle(const double inverse_velocity_minimum,
                          const double inverse_velocity_range,
                          const double time_minimum, const double time_range) {
    if (!(inverse_velocity_range > 0) || !(time_range > 0)) {
      throw std::invalid_argument(
          "A source rectangle needs a positive width in both coordinates; got "
          + std::to_string(inverse_velocity_range) + " s/m by "
          + std::to_string(time_range) + " s");
    }
    Region region;
    const chopper_polygon polygon = chopper_polygon_rectangle(
        inverse_velocity_minimum, inverse_velocity_range, time_minimum, time_range);
    if (!chopper_polygon_set_add(&region.set_, &polygon)) {
      throw std::runtime_error("Out of memory building a source rectangle");
    }
    return region;
  }

  [[nodiscard]] const chopper_polygon_set & raw() const { return set_; }
  [[nodiscard]] unsigned count() const { return set_.count; }
  [[nodiscard]] double area() const { return chopper_polygon_set_area(&set_); }

  [[nodiscard]] std::vector<chopper_polygon> polygons() const {
    return {set_.polygon, set_.polygon + set_.count};
  }

  [[nodiscard]] bool contains(const double inverse_velocity, const double time) const {
    return chopper_polygon_set_contains(&set_, inverse_velocity, time) != 0;
  }

  /// What is left of this region after a train. A new region; this one is unchanged.
  [[nodiscard]] Region transmit(const std::vector<Chopper> & choppers,
                                const std::vector<double> & path_spreads) const {
    if (!path_spreads.empty() && path_spreads.size() != choppers.size()) {
      throw std::invalid_argument(
          "One path spread per chopper, or none at all; got "
          + std::to_string(path_spreads.size()) + " for "
          + std::to_string(choppers.size()) + " choppers");
    }
    Region out(*this);
    const auto pars = chopcal::c_structs(choppers);
    if (!chopper_polygon_set_transmit_train(
            &out.set_, static_cast<unsigned>(pars.size()), pars.data(),
            path_spreads.empty() ? nullptr : path_spreads.data())) {
      throw std::runtime_error(
          "chopper-lib refused to transmit this region and printed the reason: either a "
          "path spread wide enough to make one turn of a disk overlap the next, which "
          "would make the transmitted areas count twice, or a polygon needing more "
          "vertices than it has room for.");
    }
    return out;
  }

  [[nodiscard]] std::vector<std::pair<double, double>> inverse_velocity_ranges() const {
    const range_set ranges = chopper_polygon_set_inverse_velocity_ranges(&set_);
    std::vector<std::pair<double, double>> out;
    out.reserve(ranges.count);
    for (unsigned i = 0; i < ranges.count; ++i) {
      out.emplace_back(ranges.ranges[i].minimum, ranges.ranges[i].maximum);
    }
    if (ranges.ranges) free(ranges.ranges);
    return out;
  }

  [[nodiscard]] std::vector<std::pair<double, double>> wavelength_ranges() const {
    auto out = inverse_velocity_ranges();
    for (auto & band : out) {
      band.first = inverse_velocity_to_wavelength(band.first);
      band.second = inverse_velocity_to_wavelength(band.second);
    }
    return out;
  }

  void write_json(const std::string & path,
                  const std::optional<chopper_polygon> & sampled) const {
    /* The C splits a path into directory, name and extension and joins them back; give it
     * the whole thing as the name and nothing else, so what lands on disk is exactly the
     * path asked for. */
    if (!chopper_write_polygons_to_file(nullptr, path.c_str(), nullptr, "/", &set_,
                                        sampled ? &*sampled : nullptr)) {
      throw std::runtime_error("Could not open " + path + " for writing");
    }
  }
};

/** chopper-lib's sampler over a transmitted region, owning what it allocated. */
class RegionSampler {
  chopper_polygon_sampler sampler_{};

 public:
  RegionSampler(const Region & region, const double sampled_area) {
    chopper_polygon_sampler_empty(&sampler_);
    if (!(sampled_area > 0.0)) {
      throw std::invalid_argument(
          "The sampled area sets the acceptance, so it has to be positive; got "
          + std::to_string(sampled_area));
    }
    sampler_ = chopper_polygon_sampler_make(&region.raw(), sampled_area);
  }

  RegionSampler(const RegionSampler &) = delete;
  RegionSampler & operator=(const RegionSampler &) = delete;
  ~RegionSampler() { chopper_polygon_sampler_free(&sampler_); }

  [[nodiscard]] unsigned count() const { return sampler_.count; }
  [[nodiscard]] double acceptance() const { return sampler_.acceptance; }

  void require_triangles() const {
    if (sampler_.count == 0) {
      throw std::runtime_error(
          "This sampler has no triangles to draw from: the region is empty. Check "
          "`count` before drawing.");
    }
  }

  [[nodiscard]] std::pair<double, double> draw(const double triangle_deviate,
                                               const double first_deviate,
                                               const double second_deviate) const {
    require_triangles();
    std::pair<double, double> out;
    chopper_polygon_sampler_draw(&sampler_, triangle_deviate, first_deviate,
                                 second_deviate, &out.first, &out.second);
    return out;
  }

  [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> draw_many(
      const f64_1d & triangle_deviates, const f64_1d & first_deviates,
      const f64_1d & second_deviates) const {
    require_triangles();
    const auto n = triangle_deviates.size();
    if (first_deviates.size() != n || second_deviates.size() != n) {
      throw std::invalid_argument("The three deviate arrays must be the same length");
    }
    std::vector<double> inverse_velocity(n), time(n);
    for (size_t i = 0; i < n; ++i) {
      chopper_polygon_sampler_draw(&sampler_, triangle_deviates.data()[i],
                                   first_deviates.data()[i], second_deviates.data()[i],
                                   &inverse_velocity[i], &time[i]);
    }
    return {std::move(inverse_velocity), std::move(time)};
  }
};

}  // namespace

NB_MODULE(_chopper_lib_impl, m) {

nb::class_<Chopper>(m, "Chopper",
      "One disk chopper: how fast it turns, where its openings are, and how far away.\n\n"
      "The disk is described the way the NeXus `NXdisk_chopper` standard and McStas'\n"
      "`CollectorDiskChopper` describe one. `edges` are the slit edges in degrees from\n"
      "the disk's own zero mark, two per opening, strictly increasing, opening edge\n"
      "first; `beam` is the angle from that mark to where the beam crosses the disk.\n"
      "An edge at angle `a` is on the beam at\n\n"
      "    t(a) = delay + (beam - a) / (360 * speed)\n\n"
      "and every 1/|speed| seconds after that. Only `beam - a` appears, so a caller can\n"
      "hand over the disk's own numbers with no conversion to do. A negative `speed`\n"
      "turns the disk the other way, which does not move `delay` -- a delay is a time,\n"
      "so it means the same thing in either direction -- but does swap which side of it\n"
      "an opening arrives on.\n\n"
      "This replaces the `angle` of chopper-lib 2.x, which could only describe one\n"
      "opening centred on the beam. `Chopper.centred` still writes one of those.")
      .def(nb::init<double, double, double, std::vector<double>, double, double>(),
           "speed"_a = 0, "delay"_a = 0, "beam"_a = 0,
           "edges"_a = std::vector<double>{}, "path"_a = 0, "aperture"_a = 0)
      .def_static("centred", &Chopper::centred,
                  "speed"_a = 0, "delay"_a = 0, "width"_a = 0, "path"_a = 0,
                  "aperture"_a = 0,
                  "A disk of one opening `width` degrees across, centred on the beam.\n\n"
                  "Equivalent to beam=0 and edges=[-width/2, +width/2], which is what\n"
                  "chopper-lib built internally from the `angle` of its 2.x structure.")
      .def_rw("speed", &Chopper::speed,
              "Rotation frequency in Hz; negative turns the disk the other way")
      .def_rw("delay", &Chopper::delay,
              "When the disk point at angle `beam` is on the beam path, in seconds")
      .def_rw("beam", &Chopper::beam,
              "From the disk's zero mark to where the beam crosses it, in degrees")
      .def_rw("edges", &Chopper::edges,
              "Slit edges in degrees from the zero mark, two per opening, increasing.\n"
              "A slit across the mark is written past 360 -- [350, 370], not [350, 10].")
      .def_rw("path", &Chopper::path,
              "Flight path from the source to the disk, in metres")
      .def_rw("aperture", &Chopper::aperture,
              "Angular width of the beam where it crosses the disk, in degrees.\n"
              "Widens every window by half of it at each end; 0 is a point beam.")
      .def_prop_ro("opening", &Chopper::opening,
                   "Total opening angle in degrees, summed over the openings")
      .def("validate", &Chopper::validate,
           "Raise if the edges are not an increasing, even-length list, or the\n"
           "aperture is negative. Called on construction; call it again after\n"
           "assigning to `edges`, which cannot check itself.")
      .def("parked_is_open",
           [](const Chopper & chopper) {
             return chopper_parked_is_open(chopper.c_struct()) != 0;
           },
           "Whether this disk stands open on the beam when it is not turning.\n\n"
           "`speed` is not read: a turning disk stands open once a period whatever its\n"
           "angles are. It matters for a parked one, because the window and mask\n"
           "functions leave a disk parked *open* out of the calculation -- no period,\n"
           "so it constrains nothing -- while one parked *shut* is a beam stop and\n"
           "empties the answer. An empty result does not say which happened; this does.")
      ;

m.def("inverse_velocity_windows",
      [](const std::vector<Chopper> & choppers, const double inv_v_min,
         const double inv_v_max, const double latest_emission) {
            const auto pars = chopcal::c_structs(choppers);
            const auto rs = chopper_inverse_velocity_windows(
                static_cast<unsigned>(pars.size()), pars.data(), inv_v_min, inv_v_max,
                latest_emission);
            std::vector<std::pair<double, double>> out;
            out.reserve(rs.count);
            for (unsigned i = 0; i < rs.count; ++i)
              out.emplace_back(rs.ranges[i].minimum, rs.ranges[i].maximum);
            if (rs.ranges) free(rs.ranges);
            return out;
      }, "choppers"_a, "inv_v_min"_a=INVERSE_V_MIN, "inv_v_max"_a=INVERSE_V_MAX,
         "latest_emission"_a=LATEST_EMISSION,
      "Every inverse velocity band the train passes, in s/m.\n\n"
      "A disk parked open is left out, having no period to constrain anything; one\n"
      "parked shut empties the list -- see `Chopper.parked_is_open`."
);

m.def("inverse_velocity_limits",
      [](const std::vector<Chopper> & choppers, const double inv_v_min,
         const double inv_v_max, const double latest_emission) {
            const auto pars = chopcal::c_structs(choppers);
            std::pair<double, double> out;
            auto no = chopper_inverse_velocity_limits(
                &out.first, &out.second, static_cast<unsigned>(pars.size()), pars.data(),
                inv_v_min, inv_v_max, latest_emission);
            return std::make_tuple(no, out);
      }, "choppers"_a, "inv_v_min"_a=INVERSE_V_MIN, "inv_v_max"_a=INVERSE_V_MAX,
         "latest_emission"_a=LATEST_EMISSION,
      "The envelope of those bands, as (count, (low, high)) in s/m.\n\n"
      "A count above one means the envelope also contains inverse velocities the train\n"
      "does not pass; use `inverse_velocity_windows` for the bands themselves."
);

m.def("wavelength_windows",
      [](const std::vector<Chopper> & choppers, const double lambda_min,
         const double lambda_max, const double latest_emission) {
            const auto pars = chopcal::c_structs(choppers);
            const auto rs = chopper_inverse_velocity_windows(
                static_cast<unsigned>(pars.size()), pars.data(),
                wavelength_to_inverse_velocity(lambda_min),
                wavelength_to_inverse_velocity(lambda_max), latest_emission);
            std::vector<std::pair<double, double>> out;
            out.reserve(rs.count);
            for (unsigned i = 0; i < rs.count; ++i)
              out.emplace_back(inverse_velocity_to_wavelength(rs.ranges[i].minimum),
                               inverse_velocity_to_wavelength(rs.ranges[i].maximum));
            if (rs.ranges) free(rs.ranges);
            return out;
      }, "choppers"_a, "wavelength_min"_a=LAMBDA_MIN, "wavelength_max"_a=LAMBDA_MAX,
         "latest_emission"_a=LATEST_EMISSION,
      "Every wavelength band the train passes, in angstrom.\n\n"
      "The band list rather than the envelope `wavelength_limits` returns, and the one to\n"
      "reach for when there is more than one: an envelope spanning two bands also spans\n"
      "the gap between them, which the train does not pass.\n\n"
      "See `wavelength_limits` on why a band here may be narrower than it looks."
);

m.def("wavelength_limits",
      [](const std::vector<Chopper> & choppers, const double lambda_min,
         const double lambda_max, const double latest_emission) {
            const auto pars = chopcal::c_structs(choppers);
            std::pair<double, double> out;
            auto no = chopper_wavelength_limits(
                &out.first, &out.second, static_cast<unsigned>(pars.size()), pars.data(),
                lambda_min, lambda_max, latest_emission);
            return std::make_tuple(no, out);
      }, "choppers"_a, "wavelength_min"_a=LAMBDA_MIN, "wavelength_max"_a=LAMBDA_MAX,
         "latest_emission"_a=LATEST_EMISSION,
      "The same envelope in angstrom, as (count, (low, high)).\n\n"
      "A count above one means the envelope also spans wavelengths the train does not\n"
      "pass; `wavelength_windows` gives the bands themselves.\n\n"
      "This and `wavelength_windows` are an *over*-approximation. chopper-lib works out\n"
      "each disk's admissible inverse velocities letting the emission time range over the\n"
      "whole pulse independently, then intersects those ranges, so it can report a band\n"
      "that no single emission time actually delivers -- an intersection of projections\n"
      "is a superset of the projection of the intersection. It is tight where the disks\n"
      "leave wide, overlapping emission windows, which is the usual case and includes\n"
      "every band an instrument is set up to pass. `inverse_velocity_time_mask` keeps the\n"
      "two coordinates coupled and does not make this error; where the two disagree, the\n"
      "mask is right."
);

m.def("beam_aperture", &chopcal::beam_aperture,
      "radius"_a, "slit_height"_a, "window_width"_a, "window_height"_a,
      "The angular width of a beam where it crosses a disk, in degrees.\n\n"
      "This is what `Chopper.aperture` wants. An opening is angular and a beam is not, so\n"
      "a neutron crossing to one side of the beam centre meets an edge early or late; the\n"
      "largest such angle belongs to the inner corners of the beam window, nearest the\n"
      "spindle, where a given width subtends the most angle.\n\n"
      "`slit_height` is the radial extent of the opening cut in the disk -- McStas\n"
      "`DiskChopper`'s `yheight`, which puts the beam centre at `radius - yheight/2`.\n"
      "`window_width` and `window_height` describe the beam inside it, which is smaller.\n\n"
      "All four lengths in metres; the answer is in degrees.")
      ;

m.attr("MASK_EXCLUDED") = static_cast<int>(CHOPPER_MASK_EXCLUDED);
m.attr("MASK_INCLUDED") = static_cast<int>(CHOPPER_MASK_INCLUDED);
m.attr("MASK_GROWN") = static_cast<int>(CHOPPER_MASK_GROWN);

m.def("inverse_velocity_time_mask",
      [](const std::vector<Chopper> & choppers, const f64_1d & inverse_velocities,
         const f64_1d & times, const int grow) {
            const auto inverse_velocity_bins = bin_count(inverse_velocities,
                                                         "inverse_velocities");
            const auto time_bins = bin_count(times, "times");
            const auto pars = chopcal::c_structs(choppers);
            auto data = std::make_unique<int[]>(
                static_cast<size_t>(inverse_velocity_bins) * time_bins);
            const auto allowed = chopper_inverse_velocity_time_mask(
                data.get(), inverse_velocity_bins, time_bins,
                inverse_velocities.data(),
                static_cast<unsigned>(inverse_velocities.size()),
                times.data(), static_cast<unsigned>(times.size()),
                pars.data(), static_cast<unsigned>(pars.size()), grow);
            return std::make_pair(
                make_mask_array(data.release(), time_bins, inverse_velocity_bins),
                allowed);
      }, "choppers"_a, "inverse_velocities"_a, "times"_a, "grow"_a = 0,
      "Which (inverse velocity, time) bins the train passes.\n\n"
      "`inverse_velocities` (s/m) and `times` (s, at the source) are bin *edges*, so the\n"
      "mask is one smaller in each direction. It comes back shaped\n"
      "**(time, inverse velocity)** -- time is the slow axis, which is the order the C\n"
      "stores it in and the opposite of the argument order.\n\n"
      "`grow` expands each allowed region by that many bins in every direction, which\n"
      "is a blunt instrument: it admits inverse velocities no disk ever passes. Prefer\n"
      "`Chopper.aperture` for a beam of finite width, which opens the windows in time\n"
      "only, where the width actually acts.\n\n"
      "Returns (mask, allowed_bin_count). Requires numpy."
);

m.def("unmasked_probability",
      [](const f64_2d & signal, const i32_2d & mask) {
            if (signal.shape(0) != mask.shape(0) || signal.shape(1) != mask.shape(1)) {
              throw std::invalid_argument("signal and mask must have the same shape");
            }
            return chopper_unmasked_probability(
                signal.data(), mask.data(), static_cast<unsigned>(mask.shape(1)),
                static_cast<unsigned>(mask.shape(0)));
      }, "signal"_a, "mask"_a,
      "The fraction of `signal` lying in the allowed bins of `mask`.\n\n"
      "Both are shaped (time, inverse velocity), as `inverse_velocity_time_mask`\n"
      "returns. This is the transmission an instrument sees -- a *weighted* fraction --\n"
      "and is not the weight correction a sampler needs; that is\n"
      "`MaskSampler.acceptance`, which counts draws rather than intensity."
);

m.attr("chopper_lib_version") = std::to_string(CHOPPER_LIB_VERSION_MAJOR) + "."
                              + std::to_string(CHOPPER_LIB_VERSION_MINOR) + "."
                              + std::to_string(CHOPPER_LIB_VERSION_PATCH);

/* The library's own conversion, exposed so a caller converting by hand agrees with
 * `wavelength_limits` and `Region.wavelength_ranges` exactly rather than in the eighth
 * digit. These are the McStas runtime's numbers; chopcal.constants.V2K is derived from
 * SI/CODATA and is deliberately a hair different.
 *
 * They are also not exact reciprocals of each other: the runtime rounds each literal on
 * its own, and V2K * K2V is 0.99999999891, so converting a wavelength to an inverse
 * velocity and back is 1.1e-9 short. chopcal.constants derives K2V as 1/V2K and does not
 * have that; these are here to agree with the library, not to be self-consistent. */
m.def("wavelength_to_inverse_velocity", &wavelength_to_inverse_velocity, "wavelength"_a,
      "Angstrom to s/m, the way chopper-lib does it.");
m.def("inverse_velocity_to_wavelength", &inverse_velocity_to_wavelength,
      "inverse_velocity"_a, "s/m to angstrom, the way chopper-lib does it.");

nb::class_<chopper_polygon>(m, "Polygon",
      "One convex piece of a transmitted region, in (inverse velocity, time).\n\n"
      "A value, not a handle: it owns nothing and copies freely. Its vertices are the\n"
      "disc edges that bound it, at the precision of a double -- nothing here is sampled\n"
      "onto a grid.")
      .def_prop_ro("vertices",
                   [](const chopper_polygon & polygon) {
                     std::vector<std::pair<double, double>> out;
                     out.reserve(polygon.count);
                     for (unsigned i = 0; i < polygon.count; ++i) {
                       out.emplace_back(polygon.vertex[i].inverse_velocity,
                                        polygon.vertex[i].time);
                     }
                     return out;
                   },
                   "(inverse velocity, time) pairs in order, s/m and s")
      .def_prop_ro("area", [](const chopper_polygon & p) { return chopper_polygon_area(&p); },
                   "Area in s^2/m -- phase space, not a count")
      .def("contains",
           [](const chopper_polygon & p, double inverse_velocity, double time) {
             return chopper_polygon_contains(&p, inverse_velocity, time) != 0;
           },
           "inverse_velocity"_a, "time"_a,
           "Whether a point lies inside, its boundary counting as inside")
      .def("__len__", [](const chopper_polygon & p) { return p.count; })
      ;

nb::class_<Region>(m, "Region",
      "The region of (inverse velocity, emission time) a chopper train transmits.\n\n"
      "A union of convex polygons, and exact: a neutron emitted at inverse velocity `a`\n"
      "and time `t` reaches path `L` at `t + L*a`, so a disc open on [lower, upper]\n"
      "accepts `lower <= t + L*a <= upper` -- a slab between two parallel lines. A\n"
      "train's acceptance is an intersection of unions of such slabs, and intersection\n"
      "distributes over union, so the region *is* a union of convex pieces: one per\n"
      "choice of which opening and which turn of each disc a neutron goes through.\n\n"
      "This is what `inverse_velocity_windows` and `inverse_velocity_time_mask`\n"
      "approximate, each in its own direction. The window functions project each disc\n"
      "separately and intersect the projections, which can report a band no single\n"
      "emission time delivers. The mask quantises, losing channels thinner than a bin\n"
      "and keeping partly covered bins whole. This does neither.\n\n"
      "Start from `Region.rectangle`, the region a source draws from, and `transmit` it\n"
      "through the train.")
      .def(nb::init<>(), "An empty region, transmitting nothing.")
      .def_static("rectangle", &Region::rectangle,
                  "inverse_velocity_minimum"_a, "inverse_velocity_range"_a,
                  "time_minimum"_a, "time_range"_a,
                  "The rectangle a source draws from: inverse velocities in s/m by\n"
                  "emission times in s. Both ranges must be positive.")
      .def_prop_ro("polygons", &Region::polygons, "The convex pieces, as a list")
      .def_prop_ro("area", &Region::area,
                   "Total area in s^2/m. A plain sum: the pieces are disjoint, because\n"
                   "two turns of one disc cannot both pass the same neutron.")
      .def("contains", &Region::contains, "inverse_velocity"_a, "time"_a,
           "Whether a neutron emitted at this (inverse velocity, time) gets through")
      .def("transmit", &Region::transmit, "choppers"_a,
           "path_spreads"_a = std::vector<double>{},
           "What a train leaves of this region, as a new Region; this one is unchanged.\n\n"
           "`path_spreads` is one extra flight path per chopper, in metres, for a guide\n"
           "that is not a straight line. The deviation is in *path*, so its effect on an\n"
           "arrival time is `deviation * inverse_velocity` -- larger for a slow neutron,\n"
           "and nothing at all to the inverse velocity -- so it does not grow the region\n"
           "evenly, it opens each slab into a wedge. It gives the *support* under that\n"
           "uncertainty rather than a distribution over it, the same bargain\n"
           "`Chopper.aperture` makes for the width of the beam.")
      .def("inverse_velocity_ranges", &Region::inverse_velocity_ranges,
           "The bands the region covers in inverse velocity, sorted and merged, s/m.\n\n"
           "The projection of the intersection, where `inverse_velocity_windows` gives\n"
           "the intersection of the projections -- a superset of this one.")
      .def("wavelength_ranges", &Region::wavelength_ranges,
           "The same bands in angstrom, converted as chopper-lib converts them.")
      .def("sampler",
           [](const Region & region, double sampled_area) {
             return std::make_unique<RegionSampler>(region, sampled_area);
           },
           "sampled_area"_a,
           "A sampler over this region, drawing uniformly and never rejecting.\n\n"
           "`sampled_area` is the area of the region the caller draws from -- normally\n"
           "the `area` of the source rectangle this was transmitted from -- and sets the\n"
           "`acceptance`.")
      .def("write_json", &Region::write_json, "path"_a,
           "sampled"_a = nb::none(),
           "Write the region to `path` as JSON, the same file the McStas component\n"
           "`Polygon_ESS_butterfly` writes: the polygons' vertices, the area, the bands\n"
           "and the acceptance, each number with enough digits to read back bit-exact.\n\n"
           "`sampled` is the source rectangle's single polygon, which sets `acceptance`;\n"
           "without it both that and `sampled` are written as null. Use `to_dict` for a\n"
           "region in memory -- this exists so a Python caller can produce byte-identical\n"
           "output to the component.")
      .def("__len__", &Region::count)
      .def("__bool__", [](const Region & r) { return r.count() != 0; })
      ;

nb::class_<RegionSampler>(m, "RegionSampler",
      "Draws (inverse velocity, time) pairs uniformly from a transmitted region.\n\n"
      "The same job as `MaskSampler` and the same shape -- three uniform deviates, one\n"
      "binary search, never rejects -- over the exact region rather than a grid\n"
      "approximation of it. Each polygon is fanned into triangles from its first vertex,\n"
      "so `count` here is triangles, not polygons.")
      .def(nb::init<const Region &, double>(), "region"_a, "sampled_area"_a)
      .def_prop_ro("count", &RegionSampler::count,
                   "Triangles this draws from; 0 means nothing to draw")
      .def_prop_ro("acceptance", &RegionSampler::acceptance,
              "The factor every drawn ray's weight must be multiplied by.\n\n"
              "Exact: the transmitted area over the sampled area, both in closed form\n"
              "rather than counted in cells. `MaskSampler.acceptance` estimates the same\n"
              "number on a grid and can only over-estimate it, since a partly covered\n"
              "cell is weighted whole.\n\n"
              "Everything the note on `MaskSampler.acceptance` says about independence\n"
              "still applies: the correction is right only while both coordinates are\n"
              "drawn uniformly and independently of each other and of everything else.")
      .def("draw", &RegionSampler::draw,
           "triangle_deviate"_a, "first_deviate"_a, "second_deviate"_a,
           "One (inverse_velocity, time) pair from three uniform deviates on [0, 1).")
      .def("draw_many", &RegionSampler::draw_many,
           "triangle_deviates"_a, "first_deviates"_a, "second_deviates"_a,
           "The same draw over three equal-length arrays, returning two lists.")
      ;

nb::class_<MaskSampler>(m, "MaskSampler",
      "Draws (inverse velocity, time) pairs from the allowed cells of a finished mask.\n\n"
      "Sampling the whole region and discarding what the mask excludes spends the ray\n"
      "budget to keep the fraction the choppers pass. Drawing from the allowed cells\n"
      "keeps all of it and is the same distribution, provided every drawn weight is\n"
      "multiplied by `acceptance`.\n\n"
      "Cells are clipped to the sampled region before they are weighted, so a grid\n"
      "sized with `ceil` -- whose last row and column run past the region -- does not\n"
      "inflate `acceptance`.")
      .def(nb::init<const i32_2d &, const f64_1d &, const f64_1d &, double, double,
                    double, double>(),
           "mask"_a, "inverse_velocities"_a, "times"_a,
           "inverse_velocity_minimum"_a, "inverse_velocity_range"_a,
           "time_minimum"_a, "time_range"_a,
           "`mask` is shaped (time, inverse velocity) and `inverse_velocities`/`times`\n"
           "are its bin edges. The four bounds describe the region the caller samples\n"
           "uniformly, which is not in general the region the grid covers.")
      .def_prop_ro("count", &MaskSampler::count,
                   "Allowed cells this draws from; 0 means nothing to draw")
      .def_prop_ro("acceptance", &MaskSampler::acceptance,
              "The factor every drawn ray's weight must be multiplied by.\n\n"
              "It is the allowed area over the sampled area -- a ratio of areas, so it\n"
              "is exact, but exact only while both coordinates are drawn uniformly and\n"
              "independently of each other and of everything else the caller samples. A\n"
              "source that picks its emission time from a window centred on the\n"
              "neutron's own velocity (McStas' `ESS_butterfly` under time focusing)\n"
              "breaks that, and no single number corrects it.\n\n"
              "It is neither `unmasked_probability`, which weights by intensity, nor\n"
              "anything recoverable from a rejection loop's trial count: for a geometric\n"
              "number of trials k, E[1/k] is not this.")
      .def("draw", &MaskSampler::draw,
           "cell_deviate"_a, "inverse_velocity_deviate"_a, "time_deviate"_a,
           "One (inverse_velocity, time) pair from three uniform deviates on [0, 1).\n\n"
           "The deviates are arguments rather than drawn here so the library needs no\n"
           "generator of its own and a caller keeps whatever it was already using.")
      .def("draw_many", &MaskSampler::draw_many,
           "cell_deviates"_a, "inverse_velocity_deviates"_a, "time_deviates"_a,
           "The same draw over three equal-length arrays, returning two lists.")
      ;

}
