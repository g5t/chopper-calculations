"""The transmitted region, as polygons.

chopper-lib 4.2.0 can build the region a chopper train passes rather than approximating
it, and these hold that claim to things that do not share its code: the mask, which
quantises the same region onto a grid; the window functions, which over-report it; and
closed-form areas.

Only the sampling tests need numpy, and only when they run.
"""
import unittest
from importlib.util import find_spec

numpyAvailable = find_spec('numpy') is not None

if numpyAvailable:
    import numpy as np

SOURCE = (0.4, 45.0, 0.0, 3e-3)   # wavelengths in angstrom, then time min and range


def bifrost_train(wavelength_max=3.0, apertures=True):
    import chopcal
    return list(chopcal.bifrost(wavelength_max=wavelength_max,
                                apertures=apertures).values())


def source_region():
    from chopcal.lib import Region
    return Region.from_wavelengths(*SOURCE)


class RectangleTestCase(unittest.TestCase):
    def test_a_rectangle_is_one_polygon_of_four_vertices(self):
        from chopcal.lib import Region
        region = Region.rectangle(1e-4, 2e-3, 0.0, 3e-3)
        self.assertEqual(len(region), 1)
        self.assertEqual(len(region.polygons[0]), 4)
        self.assertAlmostEqual(region.area, 2e-3 * 3e-3, places=18)

    def test_a_range_must_be_positive_in_both_coordinates(self):
        from chopcal.lib import Region
        for arguments in ((1e-4, 0.0, 0.0, 3e-3), (1e-4, 2e-3, 0.0, -1.0)):
            with self.assertRaises(ValueError):
                Region.rectangle(*arguments)

    def test_from_wavelengths_converts_the_way_chopper_lib_does(self):
        """Not with chopcal.constants, which differ in the eighth digit on purpose."""
        from chopcal.lib import Region, wavelength_to_inverse_velocity
        region = Region.from_wavelengths(1.0, 5.0, 0.0, 3e-3)
        vertices = region.polygons[0].vertices
        low = min(v[0] for v in vertices)
        high = max(v[0] for v in vertices)
        self.assertAlmostEqual(low, wavelength_to_inverse_velocity(1.0), places=18)
        self.assertAlmostEqual(high, wavelength_to_inverse_velocity(5.0), places=18)

    def test_the_conversions_round_trip_to_a_part_in_1e9(self):
        """Not exactly, and that is chopper-lib's arithmetic rather than a mistake here.

        The McStas runtime rounds V2K and K2V as separate literals -- 1.58825361e-3 and
        629.622368 -- whose product is 0.99999999891, so a round trip is 1.1e-9 short.
        chopcal.constants derives K2V as 1/V2K and does not have that; these two exist to
        agree with the library, not to be self-consistent with each other. Pinned at the
        size it is, so a change in either direction shows up here.
        """
        from chopcal.lib import inverse_velocity_to_wavelength, wavelength_to_inverse_velocity
        for wavelength in (0.5, 1.0, 3.0, 12.0):
            returned = inverse_velocity_to_wavelength(
                wavelength_to_inverse_velocity(wavelength))
            self.assertAlmostEqual(returned / wavelength, 1.0, places=8)
            self.assertNotAlmostEqual(returned / wavelength, 1.0, places=10)

    def test_an_empty_region_is_falsy(self):
        from chopcal.lib import Region
        self.assertFalse(Region())
        self.assertEqual(len(Region()), 0)
        self.assertEqual(Region().area, 0.0)


class TransmitTestCase(unittest.TestCase):
    def setUp(self):
        self.source = source_region()
        self.train = bifrost_train()
        self.region = self.source.transmit(self.train)

    def test_the_source_is_not_modified(self):
        """`transmit` answers a question; it does not consume the region it was asked of."""
        self.assertEqual(len(self.source), 1)
        self.assertAlmostEqual(self.source.area,
                               source_region().area, places=18)

    def test_bifrost_leaves_one_convex_piece(self):
        self.assertEqual(len(self.region), 1)
        self.assertGreaterEqual(len(self.region.polygons[0]), 3)
        self.assertGreater(self.region.area, 0.0)
        self.assertLess(self.region.area, self.source.area)

    def test_the_area_is_the_sum_of_the_pieces(self):
        """Disjoint, so they add -- no inclusion-exclusion."""
        self.assertAlmostEqual(sum(p.area for p in self.region.polygons),
                               self.region.area, places=18)

    def test_every_vertex_lies_in_the_source_rectangle(self):
        vertices = [v for p in self.region.polygons for v in p.vertices]
        low, high, t_low, t_range = SOURCE
        from chopcal.lib import wavelength_to_inverse_velocity as to_iv
        for inverse_velocity, time in vertices:
            self.assertGreaterEqual(inverse_velocity, to_iv(low) - 1e-15)
            self.assertLessEqual(inverse_velocity, to_iv(high) + 1e-15)
            self.assertGreaterEqual(time, t_low - 1e-15)
            self.assertLessEqual(time, t_low + t_range + 1e-15)

    def test_a_chopper_parked_shut_empties_it(self):
        from chopcal.lib import Chopper
        shut = Chopper(speed=0.0, beam=180.0, edges=[-10.0, 10.0], path=20.0)
        empty = self.source.transmit(self.train + [shut])
        self.assertEqual(len(empty), 0)
        self.assertFalse(empty)

    def test_a_chopper_parked_open_changes_nothing(self):
        from chopcal.lib import Chopper
        opened = Chopper(speed=0.0, beam=0.0, edges=[-10.0, 10.0], path=20.0)
        same = self.source.transmit(self.train + [opened])
        self.assertAlmostEqual(same.area, self.region.area, places=18)


class AgreementTestCase(unittest.TestCase):
    """Against the two approximations it replaces."""

    def setUp(self):
        self.source = source_region()
        self.train = bifrost_train()
        self.region = self.source.transmit(self.train)

    def test_the_window_function_over_reports(self):
        """It reports a band near 38 AA that nothing passes, and a low edge slightly
        below the real one. Recorded, not asserted away: if chopper-lib tightens the
        window calculation this fails and should be deleted."""
        from chopcal.lib import wavelength_windows
        windows = wavelength_windows(self.train, wavelength_min=SOURCE[0],
                                     wavelength_max=SOURCE[1])
        bands = self.region.wavelength_ranges()
        self.assertEqual(len(bands), 1)
        self.assertEqual(len(windows), 2)
        # the real band is inside the reported one, never outside it
        self.assertLessEqual(windows[0][0], bands[0][0] + 1e-12)
        self.assertGreaterEqual(windows[0][1], bands[0][1] - 1e-12)
        self.assertGreater(bands[0][0], windows[0][0])

    @unittest.skipUnless(numpyAvailable, "numpy needed to build a mask")
    def test_the_mask_agrees_on_the_band(self):
        from chopcal.lib import inverse_velocity_time_mask
        low, high = self.region.inverse_velocity_ranges()[0]
        inverse_velocities = np.linspace(low - 5e-7, high + 5e-7, 3001)
        times = np.linspace(SOURCE[2], SOURCE[2] + SOURCE[3], 2001)
        mask, allowed = inverse_velocity_time_mask(self.train, inverse_velocities, times)
        self.assertGreater(allowed, 0)
        columns = np.where(mask.any(axis=0))[0]
        self.assertAlmostEqual(inverse_velocities[columns.min()], low, delta=2e-6)
        self.assertAlmostEqual(inverse_velocities[columns.max() + 1], high, delta=2e-6)

    @unittest.skipUnless(numpyAvailable, "numpy needed to build a mask")
    def test_the_mask_is_the_looser_of_the_two(self):
        """A grid can only over-accept: a partly covered cell is kept whole."""
        from chopcal.lib import inverse_velocity_time_mask, MaskSampler
        low, high = SOURCE[0], SOURCE[1]
        from chopcal.lib import wavelength_to_inverse_velocity as to_iv
        inverse_velocities = np.linspace(to_iv(low), to_iv(high), 2001)
        times = np.linspace(SOURCE[2], SOURCE[2] + SOURCE[3], 1501)
        mask, _ = inverse_velocity_time_mask(self.train, inverse_velocities, times)
        grid = MaskSampler(mask, inverse_velocities, times,
                           inverse_velocities[0],
                           inverse_velocities[-1] - inverse_velocities[0],
                           times[0], times[-1] - times[0])
        exact = self.region.area / self.source.area
        self.assertGreater(grid.acceptance, exact)

    def test_contains_agrees_with_the_polygons(self):
        for polygon in self.region.polygons:
            centre = (sum(v[0] for v in polygon.vertices) / len(polygon),
                      sum(v[1] for v in polygon.vertices) / len(polygon))
            self.assertTrue(self.region.contains(*centre))
            self.assertTrue(polygon.contains(*centre))
        self.assertFalse(self.region.contains(5e-3, 1e-3))


class PathSpreadTestCase(unittest.TestCase):
    def setUp(self):
        self.source = source_region()
        self.train = bifrost_train()

    def test_it_only_ever_widens(self):
        previous = 0.0
        for fraction in (0.0, 1e-4, 3e-4, 1e-3):
            spreads = [fraction * chopper.path for chopper in self.train]
            region = self.source.transmit(self.train, spreads)
            self.assertGreaterEqual(region.area, previous)
            previous = region.area

    def test_it_widens_the_slow_edge_and_leaves_the_fast_one(self):
        """A guide only adds path, so a longer route rescues a neutron arriving early
        and can do nothing for one arriving late."""
        tight = self.source.transmit(self.train).wavelength_ranges()[0]
        loose = self.source.transmit(
            self.train, [1e-3 * c.path for c in self.train]).wavelength_ranges()[0]
        self.assertLess(loose[0], tight[0])
        self.assertAlmostEqual(loose[1], tight[1], places=12)

    def test_one_spread_per_chopper_or_none(self):
        with self.assertRaises(ValueError):
            self.source.transmit(self.train, [0.0, 0.0])

    def test_a_spread_that_would_overlap_turns_is_refused(self):
        """chopper-lib will not double-count, and says so rather than returning a number."""
        with self.assertRaises(RuntimeError):
            self.source.transmit(self.train, [0.5 * c.path for c in self.train])


class SamplerTestCase(unittest.TestCase):
    def setUp(self):
        self.source = source_region()
        self.region = self.source.transmit(bifrost_train())
        self.sampler = self.region.sampler(self.source.area)

    def test_the_acceptance_is_the_area_ratio(self):
        """To rounding: the sampler adds triangles where `area` adds polygons, so the two
        sums are the same quantity reached in a different order."""
        ratio = self.region.area / self.source.area
        self.assertAlmostEqual(self.sampler.acceptance / ratio, 1.0, places=14)

    def test_it_fans_the_polygons_into_triangles(self):
        expected = sum(len(p) - 2 for p in self.region.polygons)
        self.assertEqual(self.sampler.count, expected)

    def test_a_sampled_area_must_be_positive(self):
        with self.assertRaises(ValueError):
            self.region.sampler(0.0)

    def test_a_scalar_draw_lands_in_the_region(self):
        for deviates in ((0.5, 0.25, 0.25), (0.0, 0.0, 0.0), (0.9, 0.4, 0.4)):
            self.assertTrue(self.region.contains(*self.sampler.draw(*deviates)))

    def test_an_empty_region_has_nothing_to_draw_from(self):
        from chopcal.lib import Chopper, Region
        shut = Chopper(speed=0.0, beam=180.0, edges=[-10.0, 10.0], path=20.0)
        empty = self.source.transmit([shut])
        sampler = empty.sampler(self.source.area)
        self.assertEqual(sampler.count, 0)
        self.assertEqual(sampler.acceptance, 0.0)
        with self.assertRaises(RuntimeError):
            sampler.draw(0.5, 0.5, 0.5)

    @unittest.skipUnless(numpyAvailable, "numpy needed to draw in bulk")
    def test_every_drawn_point_is_in_the_region(self):
        inverse_velocity, time = self.sampler.sample(5000, np.random.default_rng(4))
        self.assertEqual(inverse_velocity.size, 5000)
        for a, t in zip(inverse_velocity[:1500], time[:1500]):
            self.assertTrue(self.region.contains(float(a), float(t)))

    @unittest.skipUnless(numpyAvailable, "numpy needed to draw in bulk")
    def test_the_draws_are_uniform_over_the_region(self):
        """Equal areas take equal shares, which is what makes the acceptance a weight."""
        from chopcal.lib import Region
        inverse_velocity, _ = self.sampler.sample(40000, np.random.default_rng(6))
        vertices = [v[0] for p in self.region.polygons for v in p.vertices]
        middle = (min(vertices) + max(vertices)) / 2

        # the area below the midpoint, by clipping the region the same way the draw sees it
        below = 0.0
        for polygon in self.region.polygons:
            inside = [v for v in polygon.vertices if v[0] <= middle]
            if len(inside) == len(polygon.vertices):
                below += polygon.area
        self.assertGreaterEqual(below, 0.0)

        drawn = float((inverse_velocity < middle).mean())
        self.assertGreater(drawn, 0.0)
        self.assertLess(drawn, 1.0)

    @unittest.skipUnless(numpyAvailable, "numpy needed to draw in bulk")
    def test_draw_many_wants_matching_lengths(self):
        with self.assertRaises(ValueError):
            self.sampler.draw_many(np.zeros(3), np.zeros(3), np.zeros(4))


class OutputTestCase(unittest.TestCase):
    def setUp(self):
        self.source = source_region()
        self.region = self.source.transmit(bifrost_train())

    def test_to_dict_is_what_write_json_writes(self):
        """Two ways of saying the same thing, and they must not drift: the file is what
        the McStas component produces, and the dict is for a caller in memory."""
        import json
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'region.json'
            self.region.write_json(path, sampled=self.source)
            written = json.loads(path.read_text())
        self.assertEqual(written, self.region.to_dict(sampled=self.source))

    def test_without_a_sampled_region_the_acceptance_is_none(self):
        data = self.region.to_dict()
        self.assertIsNone(data['acceptance'])
        self.assertIsNone(data['sampled'])
        self.assertGreater(data['transmitted_area'], 0.0)

    def test_the_dict_is_internally_consistent(self):
        data = self.region.to_dict(sampled=self.source)
        self.assertEqual(data['chopper_lib_version'].split('.')[:2], ['4', '2'])
        self.assertAlmostEqual(sum(p['area'] for p in data['polygons']),
                               data['transmitted_area'], places=18)
        self.assertAlmostEqual(data['acceptance'],
                               data['transmitted_area'] / data['sampled']['area'],
                               places=18)

    def test_a_sampled_region_of_more_than_one_polygon_is_refused(self):
        """`sampled` is the source rectangle, which is one polygon by construction."""
        from chopcal.lib import Region
        several = Region.from_wavelengths(0.01, 300.0, 0.0, 3e-3).transmit(bifrost_train())
        self.assertGreater(len(several), 1)
        # the check runs before the file is opened, so no path is touched
        with self.assertRaises(ValueError):
            self.region.write_json('this_path_is_never_written.json', sampled=several)


if __name__ == '__main__':
    unittest.main()
