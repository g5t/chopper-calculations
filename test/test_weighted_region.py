"""The transmitted region weighted by how likely each neutron is to get through.

chopper-lib's own tests hold the weights to their definition; these check what chopcal
adds on top -- the binding, the grid integration and the arrival times -- and what the
weighting should and should not change about BIFROST.
"""
import unittest
from importlib.util import find_spec

numpyAvailable = find_spec('numpy') is not None
if numpyAvailable:
    import numpy as np

SOURCE = (0.4, 45.0, 0.0, 3e-3)   # wavelengths in angstrom, then time min and range
#: the pulse-shaping pair and the bandwidth pair each share a housing
BIFROST_GROUPS = [0, 0, 1, 2, 3, 3]


def train(apertures=True):
    import chopcal
    return list(chopcal.bifrost(wavelength_max=3.0, apertures=apertures).values())


def source():
    from chopcal.lib import Region
    return Region.from_wavelengths(*SOURCE)


class WeightedRegionTestCase(unittest.TestCase):
    def test_a_point_beam_has_nothing_to_ramp(self):
        """Without apertures every edge is sharp, so the weights are all one."""
        point = source().transmit(train(apertures=False))
        weighted = source().transmit_weighted(train(apertures=False), BIFROST_GROUPS)
        self.assertAlmostEqual(weighted.integral / point.area, 1.0, places=12)
        self.assertTrue(all(not factors for _, factors in weighted.pieces))

    def test_the_weights_pass_less_than_the_support_and_about_the_point_beam(self):
        """The support is where anything gets through, so the weighted integral is less.
        It is close to the point beam's area, which is the first-order answer, and not
        equal to it, because ramps overlap at BIFROST's narrow pulse-shaping openings."""
        support = source().transmit(train())
        point = source().transmit(train(apertures=False))
        weighted = source().transmit_weighted(train(), BIFROST_GROUPS)
        self.assertLess(weighted.integral, support.area)
        self.assertLess(abs(weighted.integral / point.area - 1.0), 0.2)

    def test_the_transmission_is_a_probability(self):
        weighted = source().transmit_weighted(train(), BIFROST_GROUPS)
        support = source().transmit(train())
        polygon = support.polygons[0]
        centre = [sum(v[k] for v in polygon.vertices) / len(polygon) for k in (0, 1)]
        self.assertGreater(weighted.transmission(*centre), 0.0)
        self.assertLessEqual(weighted.transmission(*centre), 1.0)
        self.assertEqual(weighted.transmission(5e-3, 1e-3), 0.0)

    def test_one_group_per_chopper_or_none(self):
        with self.assertRaises(ValueError):
            source().transmit_weighted(train(), [0, 0, 1])

    def test_a_group_must_be_listed_together(self):
        with self.assertRaises(RuntimeError):
            source().transmit_weighted(train(), [0, 1, 0, 2, 3, 3])

    def test_without_groups_every_disc_is_independent(self):
        """Then the counter-rotating bandwidth pair multiplies rather than intersecting,
        which passes a different amount."""
        grouped = source().transmit_weighted(train(), BIFROST_GROUPS)
        independent = source().transmit_weighted(train())
        self.assertNotAlmostEqual(grouped.integral, independent.integral, places=12)


@unittest.skipUnless(numpyAvailable, 'needs numpy')
class WeightedGridTestCase(unittest.TestCase):
    def setUp(self):
        self.weighted = source().transmit_weighted(train(), BIFROST_GROUPS)
        vertices = np.array([v for p, _ in self.weighted.pieces for v in p.vertices])
        (u0, t0), (u1, t1) = vertices.min(axis=0), vertices.max(axis=0)
        self.u = np.linspace(u0 - 0.1 * (u1 - u0), u1 + 0.1 * (u1 - u0), 61)
        self.t = np.linspace(t0 - 0.1 * (t1 - t0), t1 + 0.1 * (t1 - t0), 41)
        self.areas = np.outer(np.diff(self.t), np.diff(self.u))

    def test_the_weighted_coverage_adds_up_to_the_integral(self):
        coverage = self.weighted.coverage(self.u, self.t)
        self.assertAlmostEqual((coverage * self.areas).sum() / self.weighted.integral, 1.0,
                               places=10)
        self.assertLessEqual(coverage.max(), 1.0 + 1e-12)

    def test_a_uniform_density_gives_the_integral_times_it(self):
        total = self.weighted.integrate(3.0 * self.areas, self.u, self.t)
        self.assertAlmostEqual(total / (3.0 * self.weighted.integral), 1.0, places=10)

    def test_the_arrivals_lose_nothing_and_never_exceed_the_support(self):
        signal = np.random.default_rng(6).uniform(size=self.areas.shape)
        frame = 1 / 14
        edges = np.linspace(0, frame, 300)
        weighted = self.weighted.arrivals(signal, self.u, self.t, 162.0, edges, period=frame)
        support = source().transmit(train()).arrivals(signal, self.u, self.t, 162.0, edges,
                                                      period=frame)
        self.assertAlmostEqual(weighted.sum() / self.weighted.integrate(signal, self.u, self.t),
                               1.0, places=10)
        self.assertTrue(np.all(weighted <= support + 1e-12))

    def test_the_weighted_arrivals_are_softer_at_the_edges(self):
        """Fewer neutrons near the edges of the support's arrival window, the same in the
        middle: the ramps reshape the profile at its edges."""
        signal = np.ones(self.areas.shape)
        edges = np.linspace(0.04, 0.075, 400)
        weighted = self.weighted.arrivals(signal, self.u, self.t, 162.0, edges).sum(axis=1)
        support = source().transmit(train()).arrivals(signal, self.u, self.t, 162.0,
                                                      edges).sum(axis=1)
        reached = np.nonzero(support)[0]
        first, last = reached[0], reached[-1]
        self.assertLess(weighted[first], support[first])
        self.assertLess(weighted[last], support[last])


if __name__ == '__main__':
    unittest.main()
