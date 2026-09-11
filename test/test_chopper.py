"""The disk itself: how chopcal describes one, and what it refuses to describe.

chopper-lib 4.0.0 replaced the single `angle` of a chopper with the disk's own slit
`edges` and the `beam` angle they are measured against. These hold that description to
what the library means by it, and hold the translation of BIFROST's single-opening disks
to the pair `angle` always stood for.
"""
import unittest

from chopcal.lib import Chopper


class EdgeValidationTestCase(unittest.TestCase):
    """What the library does not check for itself.

    It guards an edge count below two and a null pointer and stops there. An odd count is
    divided by two to get the opening count, so the last edge is dropped and a plausible
    wrong answer comes back; unordered edges are read as widths and go negative. Neither
    is recoverable downstream, so both are refused here.
    """

    def test_an_odd_number_of_edges_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            Chopper(speed=14.0, edges=[0.0, 10.0, 20.0], path=1.0)
        self.assertIn('even', str(caught.exception))

    def test_edges_must_strictly_increase(self):
        for edges in ([10.0, 0.0], [-5.0, 5.0, 3.0, 8.0], [0.0, 0.0]):
            with self.assertRaises(ValueError, msg=f'{edges}'):
                Chopper(speed=14.0, edges=edges, path=1.0)

    def test_edges_may_not_span_more_than_a_turn(self):
        with self.assertRaises(ValueError):
            Chopper(speed=14.0, edges=[-10.0, 400.0], path=1.0)

    def test_a_slit_across_the_mark_is_written_past_360(self):
        """{350, 370}, not {350, 10} -- so the pairs stay ordered and a width subtracts."""
        chopper = Chopper(speed=14.0, edges=[350.0, 370.0], path=1.0)
        self.assertEqual(chopper.opening, 20.0)

    def test_negative_edges_are_allowed(self):
        """NXdisk_chopper wants them positive; the library and the component do not.

        A single opening astride the mark reads better as {-85, 85} than as {275, 445}.
        """
        self.assertEqual(Chopper(speed=14.0, edges=[-85.0, 85.0], path=1.0).opening, 170.0)

    def test_an_aperture_has_no_sign(self):
        with self.assertRaises(ValueError):
            Chopper(speed=14.0, edges=[-10.0, 10.0], path=1.0, aperture=-1.0)

    def test_a_disk_with_no_openings_is_describable(self):
        """A blank disk is a beam stop, not a malformed one."""
        self.assertEqual(Chopper(speed=14.0, edges=[], path=1.0).opening, 0.0)


class CentredTestCase(unittest.TestCase):
    """`Chopper.centred` is what the pre-4.0.0 `angle` field meant."""

    def test_it_is_the_pair_the_library_used_to_build_itself(self):
        """chopper-lib's own `single_to_multi_chopper` built {-angle/2, +angle/2}."""
        chopper = Chopper.centred(speed=14.0, delay=0.01, width=170.0, path=6.342)
        self.assertEqual(chopper.beam, 0.0)
        self.assertEqual(list(chopper.edges), [-85.0, 85.0])

    def test_it_leaves_a_point_beam_behind(self):
        """Zero is what a caller that does not set an aperture gets, and what every
        field but that one describes."""
        self.assertEqual(Chopper.centred(speed=14.0, width=10.0).aperture, 0.0)

    def test_opening_is_the_total_over_the_openings(self):
        """It replaces `angle`, which could only name one. For one opening they agree."""
        self.assertEqual(Chopper.centred(speed=14.0, width=38.26).opening, 38.26)
        several = Chopper(speed=14.0, edges=[0.0, 10.0, 90.0, 105.0, 180.0, 195.0])
        self.assertAlmostEqual(several.opening, 40.0)


class ParkedTestCase(unittest.TestCase):
    """A disk that is not turning is open or shut for good.

    Every disk with a speed of zero used to drop out of the calculation, which is right
    for one parked open -- no period, so it constrains nothing -- and wrong for one
    parked shut, which passes nothing at any time and was reported as a band the
    instrument would not deliver. 4.0.0 split the two cases; this is chopcal's first
    sight of the difference.
    """

    def test_the_beam_inside_an_opening_is_open(self):
        self.assertTrue(Chopper(speed=0.0, beam=0.0, edges=[-10.0, 10.0]).parked_is_open())

    def test_the_beam_on_the_disk_body_is_shut(self):
        self.assertFalse(Chopper(speed=0.0, beam=180.0, edges=[-10.0, 10.0]).parked_is_open())

    def test_angles_fold(self):
        """An opening written across the mark, and a negative beam, both work."""
        self.assertTrue(Chopper(speed=0.0, beam=-5.0, edges=[350.0, 370.0]).parked_is_open())
        self.assertTrue(Chopper(speed=0.0, beam=355.0, edges=[-10.0, 10.0]).parked_is_open())

    def test_a_disk_with_no_openings_is_shut(self):
        self.assertFalse(Chopper(speed=0.0, beam=0.0, edges=[]).parked_is_open())

    def test_speed_is_not_read(self):
        """A turning disk stands open once a period whatever its angles are, so this
        answers about the parked disk regardless of what `speed` happens to hold."""
        turning = Chopper(speed=14.0, beam=180.0, edges=[-10.0, 10.0])
        self.assertFalse(turning.parked_is_open())


class ParkedShutEmptiesTheBandTestCase(unittest.TestCase):
    """A disk parked shut is a beam stop, and the train admits nothing."""

    def setUp(self):
        from chopcal import bifrost
        self.train = list(bifrost(wavelength_max=3.0).values())

    def test_the_train_passes_a_band_to_begin_with(self):
        from chopcal.lib import wavelength_limits
        count, _ = wavelength_limits(self.train)
        self.assertEqual(count, 1)

    def test_a_disk_parked_shut_empties_it(self):
        from chopcal.lib import inverse_velocity_windows
        shut = Chopper(speed=0.0, beam=180.0, edges=[-10.0, 10.0], path=20.0)
        self.assertFalse(shut.parked_is_open())
        self.assertEqual(inverse_velocity_windows(self.train + [shut]), [])

    def test_a_disk_parked_open_changes_nothing(self):
        """No period, so it constrains no inverse velocity and drops out."""
        from chopcal.lib import inverse_velocity_windows
        opened = Chopper(speed=0.0, beam=0.0, edges=[-10.0, 10.0], path=20.0)
        self.assertTrue(opened.parked_is_open())
        self.assertEqual(inverse_velocity_windows(self.train + [opened]),
                         inverse_velocity_windows(self.train))


class DisplayTestCase(unittest.TestCase):
    def setUp(self):
        from chopcal import bifrost
        self.settings = bifrost(wavelength_max=3.0)

    def test_repr_names_the_new_fields(self):
        text = repr(self.settings['bw1'])
        for field in ('speed', 'delay', 'beam', 'edges', 'path', 'aperture'):
            self.assertIn(field, text)
        self.assertNotIn('angle', text)

    def test_str_describes_a_centred_opening_as_centred(self):
        self.assertIn('centred on the beam', str(self.settings['bw1']))

    def test_str_places_an_off_centre_opening(self):
        chopper = Chopper(speed=14.0, edges=[10.0, 30.0], path=1.0)
        self.assertIn('20 deg opening at 20 deg', str(chopper))

    def test_str_counts_several_openings(self):
        chopper = Chopper(speed=14.0, edges=[0.0, 10.0, 90.0, 100.0], path=1.0)
        self.assertIn('2 openings', str(chopper))

    def test_openings_pairs_the_edges(self):
        self.assertEqual(self.settings['bw1'].openings, [(-80.5, 80.5)])

    def test_the_table_shows_beam_and_opening(self):
        table = repr(self.settings)
        self.assertIn('beam [deg]', table)
        self.assertIn('opening [deg]', table)


if __name__ == '__main__':
    unittest.main()
