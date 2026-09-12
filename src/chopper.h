//
// An owning counterpart to chopper-lib's `chopper_parameters`.
//
#ifndef CHOPCAL_CHOPPER_H
#define CHOPCAL_CHOPPER_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants.h"

extern "C" {
#include <chopper-lib.h>
}

// The structure gained a pointer and a count in 4.0.0, and `edges` mean something
// different from the `windows` they replaced. A build against an older library would
// fail on the field names anyway, but say why.
// 4.2.1 rather than 4.2.0, which has the polygons but merges overlapping ranges wrongly:
// it loses the extent of a range containing the one after it, and leaves the answer to
// whichever order `qsort` puts tied lower edges in, which is not fixed across platforms.
// `wavelength_windows` on a real train therefore returned a different band on Windows
// than on Linux. Nothing here calls the broken function directly, so this guard is the
// only thing standing between a caller and a platform-dependent answer.
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40201
#error "chopcal builds transmitted regions as polygons, and needs range merging that does not depend on the platform's qsort; chopper-lib 4.2.1 or newer is required"
#endif

namespace chopcal {

/** The angular width of a beam where it crosses a disk, in degrees.
 *
 * An opening is angular and a beam is not, so a neutron crossing the disk to one side of
 * the beam centre meets an edge before or after one crossing at the centre does. The
 * angle between them is what `chopper_parameters::aperture` widens every window by, half
 * at each end.
 *
 * The largest such angle belongs to the *inner* corners of the beam window -- nearest the
 * spindle, where a given width subtends the most angle -- so that is the radius the
 * aperture is taken at:
 *
 *     beam_centre = radius - slit_height / 2      (McStas DiskChopper's delta_y)
 *     inner       = beam_centre - window_height / 2
 *     aperture    = 2 * atan2(window_width / 2, inner)
 *
 * `slit_height` is the radial extent of the opening and `window_height` that of the beam
 * inside it; the beam is the smaller of the two, and it is the beam that bounds where a
 * neutron can be.
 *
 * This follows McStas `DiskChopper`, whose slit is bounded on the inside by a circle of
 * radius `radius - yheight` -- it absorbs below that -- so a depth is measured from the
 * rim at the beam centre. chopper-lib's header instead measures it from where the rim has
 * dropped to at the *edge* of the window, `sqrt(radius^2 - (width/2)^2)`, which is
 * `radius - sqrt(radius^2 - (width/2)^2)` further in: 0.3 mm for the narrow BIFROST disks
 * and 1.3 mm for the bandwidth pair, worth 0.006 and 0.067 degrees of aperture. The
 * header's is the more conservative of the two by that margin. They are the same
 * quantity, computed against the two different inner boundaries those conventions
 * describe; this one matches the disks as the instrument defines them.
 *
 * Taking the width over the beam-crossing radius instead, `window_width / beam_centre`,
 * misses the height entirely and comes out low -- 11.4 degrees against 13.4 for the
 * BIFROST bandwidth disks.
 *
 * @param radius Disk radius, m
 * @param slit_height Radial extent of the opening cut in the disk, m
 * @param window_width Width of the beam where it crosses the disk, m
 * @param window_height Radial extent of that beam, m
 * @return The angular width of the beam on the disk, in degrees
 */
inline double beam_aperture(const double radius, const double slit_height,
                            const double window_width, const double window_height) {
  if (window_width <= 0) return 0.0;
  const double beam_centre = radius - slit_height / 2.0;
  const double inner = beam_centre - window_height / 2.0;
  if (inner <= 0) {
    throw std::invalid_argument(
        "A beam window " + std::to_string(window_height) + " m deep inside a "
        + std::to_string(slit_height) + " m slit on a " + std::to_string(radius)
        + " m disk reaches the spindle; there is no radius to take an aperture at.");
  }
  return 2.0 * constants::DEGREES_PER_TURN / 2.0 / constants::PI
         * std::atan2(window_width / 2.0, inner);
}

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
   * it in chopper-lib 4.2.1 is a read, so the cast is safe, and it lives here rather than
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
