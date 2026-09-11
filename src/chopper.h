//
// An owning counterpart to chopper-lib's `chopper_parameters`.
//
#ifndef CHOPCAL_CHOPPER_H
#define CHOPCAL_CHOPPER_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include <chopper-lib.h>
}

// The structure gained a pointer and a count in 4.0.0, and `edges` mean something
// different from the `windows` they replaced. A build against an older library would
// fail on the field names anyway, but say why.
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40100
#error "chopcal describes disks by their slit edges and samples finished masks; chopper-lib 4.1.0 or newer is required"
#endif

namespace chopcal {

/** A disk chopper, owning its slit edges.
 *
 * `chopper_parameters` holds a bare `double *`, which is the right thing for McStas --
 * the array is a local in the instrument's `init()` and outlives everything that reads
 * it -- and the wrong thing to hand to Python, where the lifetime of whatever allocated
 * it is nobody's business. This owns the edges and builds the C structure on demand.
 */
struct Chopper {
  double speed{0};
  double delay{0};
  double beam{0};
  std::vector<double> edges;
  double path{0};
  double aperture{0};

  Chopper() = default;

  Chopper(const double speed, const double delay, const double beam,
          std::vector<double> edges, const double path, const double aperture)
      : speed{speed}, delay{delay}, beam{beam}, edges{std::move(edges)},
        path{path}, aperture{aperture} {
    validate();
  }

  /** A disk of one opening `width` degrees across, centred on the beam.
   *
   * Every chopper `chopcal.bifrost` describes is one of these, and it is what the
   * structure's `angle` field meant before 4.0.0 -- `single_to_multi_chopper` built
   * exactly this pair internally.
   */
  static Chopper centred(const double speed, const double delay, const double width,
                         const double path, const double aperture) {
    return Chopper{speed, delay, 0.0, {-width / 2.0, width / 2.0}, path, aperture};
  }

  /** What the library checks for itself, and what it does not.
   *
   * It guards an edge count below two and a null pointer, and stops there: an odd count
   * is divided by two to get the opening count, so the last edge is dropped and a
   * plausible wrong answer comes back. Unordered edges are read as widths and go
   * negative. Neither is recoverable downstream, so refuse both here.
   */
  void validate() const {
    if (edges.size() % 2 != 0) {
      throw std::invalid_argument(
          "A chopper needs two edges per opening, so an even number of them; got "
          + std::to_string(edges.size()) + ". chopper-lib divides the count by two and "
          "would silently ignore the last edge.");
    }
    for (size_t i = 1; i < edges.size(); ++i) {
      if (!(edges[i] > edges[i - 1])) {
        throw std::invalid_argument(
            "Slit edges must strictly increase, opening edge of each slit first; edge "
            + std::to_string(i) + " is " + std::to_string(edges[i]) + " after "
            + std::to_string(edges[i - 1]) + ". A slit across the zero mark is written "
            "past 360 -- {350, 370} rather than {350, 10}.");
      }
    }
    if (edges.size() >= 2 && (edges.back() - edges.front()) > 360.0) {
      throw std::invalid_argument(
          "Slit edges span " + std::to_string(edges.back() - edges.front())
          + " degrees, more than one turn of the disk.");
    }
    if (aperture < 0) {
      throw std::invalid_argument(
          "An aperture is a width, so it has no sign; got " + std::to_string(aperture)
          + " degrees. Zero is the point beam the other fields describe.");
    }
    if (!std::isfinite(speed) || !std::isfinite(delay) || !std::isfinite(beam)
        || !std::isfinite(path) || !std::isfinite(aperture)) {
      throw std::invalid_argument("Chopper fields must all be finite");
    }
  }

  /** The library's view of this disk.
   *
   * `chopper_parameters::edges` is `double *` rather than `const double *`; every use of
   * it in chopper-lib 4.1.0 is a read, so the cast is safe, and it lives here rather than
   * at each call site so there is one place to check that claim against a new release.
   */
  [[nodiscard]] chopper_parameters c_struct() const {
    return chopper_parameters{speed, delay, beam,
                              static_cast<unsigned>(edges.size()),
                              edges.empty() ? nullptr
                                            : const_cast<double *>(edges.data()),
                              path, aperture};
  }

  /** Total opening angle, degrees -- the sum over the openings, not the largest. */
  [[nodiscard]] double opening() const {
    double total = 0.0;
    for (size_t i = 0; i + 1 < edges.size(); i += 2) total += edges[i + 1] - edges[i];
    return total;
  }
};

/** The C structures for a train, valid only while `choppers` outlives them. */
inline std::vector<chopper_parameters> c_structs(const std::vector<Chopper> & choppers) {
  std::vector<chopper_parameters> out;
  out.reserve(choppers.size());
  for (const auto & chopper : choppers) out.push_back(chopper.c_struct());
  return out;
}

}  // namespace chopcal

#endif  // CHOPCAL_CHOPPER_H
