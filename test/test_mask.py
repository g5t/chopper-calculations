"""The (inverse velocity, time) mask, and the sampler that draws from it.

These need numpy, and cibuildwheel deliberately installs only pytest when it tests a
wheel, so they skip rather than fail there. Nothing else in chopcal needs it -- the
window and limit functions return plain tuples.
"""
import unittest
from importlib.util import find_spec

numpyAvailable = find_spec('numpy') is not None

if numpyAvailable:
    import numpy as np


def bifrost_train(wavelength_max=3.0):
    from chopcal import bifrost
    return list(bifrost(wavelength_max=wavelength_max).values())


def grid(inverse_velocity_bins=200, time_bins=300):
    """Bin *edges* covering one source period and a generous inverse velocity range."""
    return (np.linspace(1e-4, 6e-3, inverse_velocity_bins + 1),
            np.linspace(0.0, 1 / 14, time_bins + 1))


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class MaskShapeTestCase(unittest.TestCase):
    """The mask is one bin smaller than its edges, and time is the slow axis.

    The C stores it as `mask[time_index * inverse_velocity_count + velocity_index]`, so
    the shape is the reverse of the argument order. Easy to get backwards, hence a test.
    """

    def setUp(self):
        self.train = bifrost_train()
        self.inverse_velocities, self.times = grid()

    def test_the_shape_is_time_by_inverse_velocity(self):
        from chopcal.lib import inverse_velocity_time_mask
        mask, _ = inverse_velocity_time_mask(self.train, self.inverse_velocities, self.times)
        self.assertEqual(mask.shape, (len(self.times) - 1, len(self.inverse_velocities) - 1))

    def test_the_count_is_the_allowed_bins(self):
        from chopcal.lib import inverse_velocity_time_mask, MaskValue
        mask, allowed = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                                   self.times)
        self.assertEqual(allowed, int((mask == MaskValue.INCLUDED).sum()))
        self.assertGreater(allowed, 0)

    def test_a_finished_mask_holds_only_excluded_and_included(self):
        """GROWN is the growing pass's own marker and is folded back before returning."""
        from chopcal.lib import inverse_velocity_time_mask, MaskValue
        mask, _ = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                             self.times, grow=3)
        self.assertEqual(set(np.unique(mask)) - {MaskValue.EXCLUDED, MaskValue.INCLUDED},
                         set())

    def test_growing_admits_more(self):
        from chopcal.lib import inverse_velocity_time_mask
        _, plain = inverse_velocity_time_mask(self.train, self.inverse_velocities, self.times)
        _, grown = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                              self.times, grow=2)
        self.assertGreater(grown, plain)

    def test_bad_edges_raise_rather_than_exit(self):
        """The C calls exit(-1) on a count mismatch, which from Python would take the
        interpreter with it. The binding derives the counts, so it can never see one --
        but the edges themselves still have to be usable."""
        from chopcal.lib import inverse_velocity_time_mask
        with self.assertRaises(ValueError):
            inverse_velocity_time_mask(self.train, np.array([1e-3]), self.times)
        with self.assertRaises(ValueError):
            inverse_velocity_time_mask(self.train, self.inverse_velocities[::-1], self.times)


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class MaskAgreesWithTheWindowsTestCase(unittest.TestCase):
    """The mask and the window functions answer the same question two ways."""

    def test_the_allowed_rows_lie_inside_the_inverse_velocity_limits(self):
        from chopcal.lib import inverse_velocity_time_mask, inverse_velocity_limits, MaskValue
        train = bifrost_train()
        inverse_velocities, times = grid()
        mask, _ = inverse_velocity_time_mask(train, inverse_velocities, times)
        count, (low, high) = inverse_velocity_limits(train)
        self.assertEqual(count, 1)

        allowed_columns = np.where((mask == MaskValue.INCLUDED).any(axis=0))[0]
        self.assertGreater(allowed_columns.size, 0)
        # a bin is allowed if any part of it is, so compare the bin's edges
        self.assertGreaterEqual(inverse_velocities[allowed_columns.max() + 1], low)
        self.assertLessEqual(inverse_velocities[allowed_columns.min()], high)


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class UnmaskedProbabilityTestCase(unittest.TestCase):
    def setUp(self):
        self.train = bifrost_train()
        self.inverse_velocities, self.times = grid()
        from chopcal.lib import inverse_velocity_time_mask
        self.mask, self.allowed = inverse_velocity_time_mask(
            self.train, self.inverse_velocities, self.times)

    def test_a_flat_signal_gives_the_allowed_fraction(self):
        from chopcal.lib import unmasked_probability
        flat = np.ones(self.mask.shape)
        self.assertAlmostEqual(unmasked_probability(flat, self.mask),
                               self.allowed / self.mask.size)

    def test_it_weights_by_intensity(self):
        """This is the transmission an instrument sees, which is why it is not the
        sampler's weight correction: put all the signal in one allowed bin and it is 1."""
        from chopcal.lib import unmasked_probability, MaskValue
        signal = np.zeros(self.mask.shape)
        first = np.argwhere(self.mask == MaskValue.INCLUDED)[0]
        signal[first[0], first[1]] = 1.0
        self.assertAlmostEqual(unmasked_probability(signal, self.mask), 1.0)

    def test_mismatched_shapes_raise(self):
        from chopcal.lib import unmasked_probability
        with self.assertRaises(ValueError):
            unmasked_probability(np.ones((3, 4)), self.mask)


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class MaskSamplerTestCase(unittest.TestCase):
    def setUp(self):
        from chopcal.lib import inverse_velocity_time_mask, MaskSampler
        self.train = bifrost_train()
        self.inverse_velocities, self.times = grid()
        self.mask, self.allowed = inverse_velocity_time_mask(
            self.train, self.inverse_velocities, self.times)
        self.sampler = MaskSampler(
            self.mask, self.inverse_velocities, self.times,
            self.inverse_velocities[0],
            self.inverse_velocities[-1] - self.inverse_velocities[0],
            self.times[0], self.times[-1] - self.times[0])

    def test_it_draws_from_every_allowed_cell(self):
        self.assertEqual(self.sampler.count, self.allowed)

    def test_acceptance_is_the_allowed_area_over_the_sampled_area(self):
        """Uniform bins and a sampled region equal to the grid, so it is the allowed
        fraction of cells outright."""
        self.assertAlmostEqual(self.sampler.acceptance,
                               self.allowed / self.mask.size, places=12)

    def test_every_draw_lands_in_an_allowed_cell(self):
        from chopcal.lib import MaskValue
        inverse_velocity, time = self.sampler.sample(20_000, np.random.default_rng(7))
        column = np.clip(np.searchsorted(self.inverse_velocities, inverse_velocity) - 1,
                         0, self.mask.shape[1] - 1)
        row = np.clip(np.searchsorted(self.times, time) - 1, 0, self.mask.shape[0] - 1)
        self.assertTrue(bool((self.mask[row, column] == MaskValue.INCLUDED).all()))

    def test_it_never_rejects(self):
        """Every draw returns a point, which is the whole reason to use it."""
        inverse_velocity, time = self.sampler.sample(5_000, np.random.default_rng(1))
        self.assertEqual(inverse_velocity.size, 5_000)
        self.assertEqual(time.size, 5_000)

    def test_a_scalar_draw_takes_its_own_deviates(self):
        """The deviates are arguments so that chopper-lib needs no generator of its own
        and a McStas TRACE can hand over rand01()."""
        low = self.sampler.draw(0.0, 0.0, 0.0)
        high = self.sampler.draw(1.0 - 1e-12, 1.0 - 1e-12, 1.0 - 1e-12)
        self.assertEqual(len(low), 2)
        self.assertNotEqual(low, high)

    def test_draw_many_wants_matching_lengths(self):
        with self.assertRaises(ValueError):
            self.sampler.draw_many(np.zeros(3), np.zeros(3), np.zeros(4))

    def test_a_sampler_over_nothing_has_no_cells_and_refuses_to_draw(self):
        from chopcal.lib import MaskSampler
        empty = np.zeros(self.mask.shape, dtype=np.int32)
        sampler = MaskSampler(empty, self.inverse_velocities, self.times,
                              self.inverse_velocities[0], 1e-3, self.times[0], 1e-3)
        self.assertEqual(sampler.count, 0)
        self.assertEqual(sampler.acceptance, 0.0)
        with self.assertRaises(RuntimeError):
            sampler.draw(0.5, 0.5, 0.5)

    def test_a_mask_of_the_wrong_shape_is_refused(self):
        from chopcal.lib import MaskSampler
        with self.assertRaises(ValueError):
            MaskSampler(np.zeros((3, 4), dtype=np.int32), self.inverse_velocities,
                        self.times, 0.0, 1e-3, 0.0, 1e-3)

    def test_clipping_keeps_acceptance_a_fraction_of_the_sampled_region(self):
        """Sample half the time range and the acceptance is measured against that half,
        not against the whole grid."""
        from chopcal.lib import MaskSampler
        half = (self.times[-1] - self.times[0]) / 2
        narrow = MaskSampler(
            self.mask, self.inverse_velocities, self.times,
            self.inverse_velocities[0],
            self.inverse_velocities[-1] - self.inverse_velocities[0],
            self.times[0], half)
        self.assertLessEqual(narrow.count, self.sampler.count)
        self.assertGreaterEqual(narrow.acceptance, 0.0)
        self.assertLessEqual(narrow.acceptance, 1.0)


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class CoercionTestCase(unittest.TestCase):
    """The bindings take typed arrays; the Python layer takes whatever you have.

    The compiled functions read the caller's buffer rather than a copy, so they insist on
    contiguous C doubles and 32-bit ints. A caller should not have to know that, and a
    mask read back from a file is unlikely to arrive in the right width.
    """

    def setUp(self):
        self.train = bifrost_train()
        self.inverse_velocities = [1e-4 + i * 5.9e-5 for i in range(101)]
        self.times = [i * (1 / 14) / 150 for i in range(151)]

    def test_plain_lists_work_as_bin_edges(self):
        from chopcal.lib import inverse_velocity_time_mask
        mask, allowed = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                                   self.times)
        self.assertEqual(mask.shape, (150, 100))
        self.assertGreater(allowed, 0)

    def test_a_narrower_float_type_is_converted(self):
        from chopcal.lib import inverse_velocity_time_mask
        _, from_list = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                                  self.times)
        _, from_float32 = inverse_velocity_time_mask(
            self.train, np.asarray(self.inverse_velocities, dtype=np.float32),
            np.asarray(self.times))
        self.assertEqual(from_list, from_float32)

    def test_a_mask_of_another_int_width_is_converted(self):
        from chopcal.lib import inverse_velocity_time_mask, unmasked_probability, MaskSampler
        mask, allowed = inverse_velocity_time_mask(self.train, self.inverse_velocities,
                                                   self.times)
        wide = mask.astype(np.int64)
        self.assertAlmostEqual(unmasked_probability(np.ones(mask.shape), wide),
                               allowed / mask.size)
        sampler = MaskSampler(wide, self.inverse_velocities, self.times,
                              self.inverse_velocities[0],
                              self.inverse_velocities[-1] - self.inverse_velocities[0],
                              self.times[0], self.times[-1] - self.times[0])
        self.assertEqual(sampler.count, allowed)


@unittest.skipUnless(numpyAvailable, "numpy needed for mask tests")
class AcceptanceIsUnbiasedTestCase(unittest.TestCase):
    """The claim `acceptance` exists to make.

    A ray drawn from proposal q carrying weight w estimates a downstream tally as
    N E_q[w f], and f is zero outside the allowed set because the disks stop those rays.
    Restricting q to the allowed set multiplies that by 1/Q, so multiplying every
    accepted ray's weight by Q puts it back. If that is right, sampling the whole region
    and rejecting must give the same answer as drawing from the allowed cells and
    scaling -- and the second must be far quieter, since it keeps every ray.

    Fixed seeds, so this is a regression test rather than a coin toss.
    """

    SAMPLES = 200_000

    def setUp(self):
        from chopcal.lib import inverse_velocity_time_mask, MaskSampler
        self.train = bifrost_train()
        self.inverse_velocities, self.times = grid()
        self.mask, _ = inverse_velocity_time_mask(
            self.train, self.inverse_velocities, self.times)
        self.low = (self.inverse_velocities[0], self.times[0])
        self.span = (self.inverse_velocities[-1] - self.inverse_velocities[0],
                     self.times[-1] - self.times[0])
        self.sampler = MaskSampler(self.mask, self.inverse_velocities, self.times,
                                   self.low[0], self.span[0], self.low[1], self.span[1])

    @staticmethod
    def tally(inverse_velocity, time):
        """Some smooth downstream quantity; nothing depends on which one."""
        return 1.0 + 10.0 * inverse_velocity + 3.0 * time

    def uniform_estimate(self, seed):
        from chopcal.lib import MaskValue
        rng = np.random.default_rng(seed)
        inverse_velocity = self.low[0] + self.span[0] * rng.random(self.SAMPLES)
        time = self.low[1] + self.span[1] * rng.random(self.SAMPLES)
        column = np.clip(np.searchsorted(self.inverse_velocities, inverse_velocity) - 1,
                         0, self.mask.shape[1] - 1)
        row = np.clip(np.searchsorted(self.times, time) - 1, 0, self.mask.shape[0] - 1)
        kept = self.tally(inverse_velocity, time) * (self.mask[row, column]
                                                     == MaskValue.INCLUDED)
        return kept.mean(), kept.std() / np.sqrt(self.SAMPLES)

    def sampled_estimate(self, seed):
        inverse_velocity, time = self.sampler.sample(self.SAMPLES,
                                                     np.random.default_rng(seed))
        weighted = self.sampler.acceptance * self.tally(inverse_velocity, time)
        return weighted.mean(), weighted.std() / np.sqrt(self.SAMPLES)

    def test_the_two_estimators_agree(self):
        uniform, uniform_error = self.uniform_estimate(0)
        sampled, sampled_error = self.sampled_estimate(1)
        combined = np.hypot(uniform_error, sampled_error)
        self.assertLess(abs(uniform - sampled), 4 * combined,
                        f'{uniform:.6g} vs {sampled:.6g}, combined error {combined:.3g}')

    def test_the_sampler_is_far_quieter_for_the_same_ray_count(self):
        """Rejection spends the budget to keep the fraction the choppers pass, which for
        this train is well under a percent."""
        _, uniform_error = self.uniform_estimate(0)
        _, sampled_error = self.sampled_estimate(1)
        self.assertGreater(uniform_error / sampled_error, 10.0)


if __name__ == '__main__':
    unittest.main()
