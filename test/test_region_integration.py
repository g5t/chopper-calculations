"""A gridded source through the transmitted region: coverage, totals and arrival times.

The region is exact, so these are held to things that do not share its code: closed-form
areas and convolutions, the mask, and points tested one at a time with `contains`.
"""
import unittest
from importlib.util import find_spec

numpyAvailable = find_spec('numpy') is not None
if numpyAvailable:
    import numpy as np

SOURCE = (0.4, 45.0, 0.0, 3e-3)   # wavelengths in angstrom, then time min and range


def bifrost_region():
    import chopcal
    from chopcal.lib import Region
    train = list(chopcal.bifrost(wavelength_max=3.0).values())
    return Region.from_wavelengths(*SOURCE).transmit(train)


def enclosing_grid(region, columns=80, rows=60):
    """Bin edges a little wider than the region, so every piece lands in some cell."""
    vertices = np.array([v for p in region.polygons for v in p.vertices])
    (u0, t0), (u1, t1) = vertices.min(axis=0), vertices.max(axis=0)
    return (np.linspace(u0 - 0.1 * (u1 - u0), u1 + 0.1 * (u1 - u0), columns + 1),
            np.linspace(t0 - 0.1 * (t1 - t0), t1 + 0.1 * (t1 - t0), rows + 1))


def cell_areas(inverse_velocities, times):
    return np.outer(np.diff(times), np.diff(inverse_velocities))


@unittest.skipUnless(numpyAvailable, 'needs numpy')
class CoverageTestCase(unittest.TestCase):
    def test_a_grid_matching_a_rectangle_is_covered_everywhere(self):
        from chopcal.lib import Region
        region = Region.rectangle(1e-4, 2e-3, 0.0, 3e-3)
        coverage = region.coverage(np.linspace(1e-4, 2.1e-3, 5), np.linspace(0, 3e-3, 4))
        self.assertEqual(coverage.shape, (3, 4))     # (time, inverse velocity)
        np.testing.assert_allclose(coverage, 1.0, rtol=1e-12)

    def test_a_cell_half_outside_is_half_covered(self):
        from chopcal.lib import Region
        region = Region.rectangle(0.0, 1.0, 0.0, 1.0)
        coverage = region.coverage([0.5, 1.5], [0.0, 1.0])
        self.assertAlmostEqual(coverage[0, 0], 0.5, places=12)

    def test_the_covered_area_is_the_region_area(self):
        region = bifrost_region()
        u, t = enclosing_grid(region)
        covered = (region.coverage(u, t) * cell_areas(u, t)).sum()
        self.assertAlmostEqual(covered / region.area, 1.0, places=10)

    def test_every_cell_is_covered_between_nothing_and_all_of_it(self):
        region = bifrost_region()
        coverage = region.coverage(*enclosing_grid(region))
        self.assertGreaterEqual(coverage.min(), 0.0)
        self.assertLessEqual(coverage.max(), 1.0 + 1e-12)
        self.assertTrue(((coverage > 0) & (coverage < 1)).any())

    def test_the_mask_counts_partly_covered_cells_whole(self):
        """The mask is the looser of the two: wherever coverage is non-zero it allows the
        cell, so its count is at least the covered cells' and its area at least ours."""
        from chopcal.lib import inverse_velocity_time_mask
        region = bifrost_region()
        u, t = enclosing_grid(region)
        import chopcal
        mask, _ = inverse_velocity_time_mask(
            list(chopcal.bifrost(wavelength_max=3.0).values()), u, t)
        coverage = region.coverage(u, t)
        self.assertTrue(np.all(mask[coverage > 1e-9] != 0))


@unittest.skipUnless(numpyAvailable, 'needs numpy')
class IntegrateTestCase(unittest.TestCase):
    def test_a_uniform_density_gives_the_area_times_it(self):
        region = bifrost_region()
        u, t = enclosing_grid(region)
        density = 7.0
        total = region.integrate(density * cell_areas(u, t), u, t)
        self.assertAlmostEqual(total / (density * region.area), 1.0, places=10)

    def test_the_signal_must_match_the_grid(self):
        region = bifrost_region()
        u, t = enclosing_grid(region)
        with self.assertRaises(ValueError):
            region.arrivals(np.ones((3, 3)), u, t, 0.0, t)


@unittest.skipUnless(numpyAvailable, 'needs numpy')
class ArrivalsTestCase(unittest.TestCase):
    def test_at_the_source_the_arrivals_are_the_covered_signal(self):
        region = bifrost_region()
        u, t = enclosing_grid(region)
        signal = np.random.default_rng(3).uniform(size=(len(t) - 1, len(u) - 1))
        arrived = region.arrivals(signal, u, t, 0.0, t)
        np.testing.assert_allclose(arrived, signal * region.coverage(u, t),
                                   rtol=1e-9, atol=1e-12)

    def test_a_square_shears_into_a_triangle(self):
        """Uniform in u and t over a 1 ms by 1 ms/m square, 1 m on: t + u is the sum of two
        uniforms, a triangle over 2 ms. Quarters of it hold 1/8, 3/8, 3/8 and 1/8."""
        from chopcal.lib import Region
        region = Region.rectangle(0.0, 1e-3, 0.0, 1e-3)
        u, t = [0.0, 1e-3], [0.0, 1e-3]
        arrived = region.arrivals([[1.0]], u, t, 1.0, np.linspace(0, 2e-3, 5))
        np.testing.assert_allclose(arrived[:, 0], [1 / 8, 3 / 8, 3 / 8, 1 / 8], rtol=1e-12)

    def test_a_frame_folds_late_arrivals_onto_early_ones(self):
        from chopcal.lib import Region
        region = Region.rectangle(0.0, 1e-3, 0.0, 1e-3)
        arrived = region.arrivals([[1.0]], [0.0, 1e-3], [0.0, 1e-3], 1.0,
                                  [0.0, 5e-4, 1e-3], period=1e-3)
        np.testing.assert_allclose(arrived[:, 0], [0.5, 0.5], rtol=1e-12)

    def test_nothing_is_lost_or_made_on_the_way(self):
        region = bifrost_region()
        u, t = enclosing_grid(region)
        signal = np.random.default_rng(4).uniform(size=(len(t) - 1, len(u) - 1))
        frame = 1 / 14
        arrived = region.arrivals(signal, u, t, 162.0, np.linspace(0, frame, 200),
                                  period=frame)
        self.assertAlmostEqual(arrived.sum() / region.integrate(signal, u, t), 1.0,
                               places=10)

    def test_the_sample_sees_what_neutrons_tested_one_by_one_see(self):
        """Points drawn uniformly over the source, kept when `contains` says they pass,
        and binned by arrival time: the same profile within counting error."""
        region = bifrost_region()
        u, t = enclosing_grid(region, columns=40, rows=30)
        rng = np.random.default_rng(5)
        n = 400_000
        a = rng.uniform(u[0], u[-1], n)
        e = rng.uniform(t[0], t[-1], n)
        kept = np.fromiter((region.contains(x, y) for x, y in zip(a, e)), bool, n)
        distance = 162.0
        edges = np.linspace((e + distance * a)[kept].min(), (e + distance * a)[kept].max(), 25)
        counted, _ = np.histogram((e + distance * a)[kept], bins=edges)
        # a uniform density over the grid, n points over the whole grid area
        density = n / ((u[-1] - u[0]) * (t[-1] - t[0]))
        arrived = region.arrivals(density * cell_areas(u, t), u, t, distance, edges).sum(axis=1)
        sigma = np.sqrt(np.maximum(counted, 1))
        self.assertLess(np.max(np.abs(arrived - counted) / sigma), 5.0)


if __name__ == '__main__':
    unittest.main()
