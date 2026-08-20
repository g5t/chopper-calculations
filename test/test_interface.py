"""What `chopcal.bifrost` asks for and what it hands back."""
import unittest


def band(settings):
    """The wavelength range the train passes, in angstrom."""
    from chopcal.lib import wavelength_limits
    count, (low, high) = wavelength_limits(list(settings.values()))
    assert count == 1, f'expected one band, got {count}'
    return low, high


def to_energy(wavelength):
    return 81.8042103582802 / (wavelength * wavelength)


class ArgumentSenseTestCase(unittest.TestCase):
    """The one number given is the *slow* edge of the band, not a cap on it.

    Both names this used to carry said the opposite -- `e_max` for what is really the
    lowest energy passed, `lambda_min` for the longest wavelength -- so these pin the
    direction rather than trusting a name.
    """

    def test_wavelength_max_is_the_top_of_the_band(self):
        from chopcal import bifrost
        for requested in (2.0, 3.0, 4.0):
            low, high = band(bifrost(wavelength_max=requested))
            self.assertAlmostEqual(high, requested, delta=0.03)
            self.assertLess(low, requested)

    def test_energy_min_is_the_bottom_of_the_band(self):
        from chopcal import bifrost
        for requested in (2.0, 5.0, 7.0, 14.0):
            low, high = band(bifrost(energy_min=requested))
            # long wavelength is low energy, so the top of the band in angstrom is the
            # bottom of it in meV
            self.assertAlmostEqual(to_energy(high), requested,
                                   delta=0.02 * requested)
            self.assertGreater(to_energy(low), requested)

    def test_the_band_is_a_fixed_width(self):
        from chopcal import bifrost, BIFROST_BANDWIDTH
        for requested in (2.5, 3.0, 4.0):
            low, high = band(bifrost(wavelength_max=requested))
            self.assertAlmostEqual(high - low, BIFROST_BANDWIDTH, delta=0.01)

    def test_either_unit_places_the_same_band(self):
        from chopcal import bifrost
        wavelength = 3.0
        by_wavelength = band(bifrost(wavelength_max=wavelength))
        by_energy = band(bifrost(energy_min=to_energy(wavelength)))
        for a, b in zip(by_wavelength, by_energy):
            self.assertAlmostEqual(a, b, delta=0.01)

    def test_a_wavelength_overrides_an_energy(self):
        from chopcal import bifrost
        both = band(bifrost(energy_min=14.0, wavelength_max=3.0))
        self.assertEqual(both, band(bifrost(wavelength_max=3.0)))

    def test_asking_for_neither_is_an_error(self):
        from chopcal import bifrost
        # the compiled function answers this with infinite delays
        with self.assertRaises(ValueError):
            bifrost()
        with self.assertRaises(ValueError):
            bifrost(0.0, 0.0)

    def test_the_positional_call_still_means_what_it_did(self):
        from chopcal import bifrost
        self.assertEqual(band(bifrost(7.0, 0, 0.0002)),
                         band(bifrost(energy_min=7.0, shaping_time=0.0002)))


class ChopperSetTestCase(unittest.TestCase):
    def setUp(self):
        from chopcal import bifrost
        self.settings = bifrost(energy_min=14.0)

    def test_it_is_still_a_dict(self):
        self.assertIsInstance(self.settings, dict)
        self.assertEqual(len(self.settings), 6)
        self.assertEqual(sorted({**self.settings}),
                         ['bw1', 'bw2', 'fo1', 'fo2', 'ps1', 'ps2'])
        self.assertEqual(self.settings['fo1'].path, 8.530)

    def test_it_comes_out_in_beam_order(self):
        paths = [chopper.path for chopper in self.settings.values()]
        self.assertEqual(paths, sorted(paths))
        self.assertEqual(list(self.settings),
                         ['ps1', 'ps2', 'fo1', 'fo2', 'bw1', 'bw2'])

    def test_it_prints_a_table_rather_than_addresses(self):
        text = repr(self.settings)
        self.assertNotIn('0x', text)
        self.assertNotIn('object at', text)
        for name in ('ps1', 'fo2', 'bw2'):
            self.assertIn(name, text)
        # one header line and one line per chopper
        self.assertEqual(len(text.splitlines()), 7)

    def test_ipython_uses_that_table_too(self):
        """A dict subclass is pretty-printed as a dict unless it says otherwise."""
        pretty = __import__('IPython.lib.pretty', fromlist=['pretty']).pretty
        self.assertEqual(pretty(self.settings), repr(self.settings))

    def test_it_offers_a_table_to_notebooks(self):
        html = self.settings._repr_html_()
        self.assertIn('<table>', html)
        self.assertIn('ps1', html)


class ChopperTestCase(unittest.TestCase):
    def setUp(self):
        from chopcal import bifrost
        self.settings = bifrost(energy_min=14.0)

    def test_repr_names_the_fields_and_their_stored_values(self):
        chopper = self.settings['bw2']
        text = repr(chopper)
        self.assertNotIn('0x', text)
        for field in ('speed', 'delay', 'angle', 'path'):
            self.assertIn(field, text)
        # the values are exact, not rounded for display
        self.assertIn(repr(chopper.delay), text)

    def test_repr_is_what_a_container_shows(self):
        self.assertEqual(repr([self.settings['ps1']]),
                         f'[{repr(self.settings["ps1"])}]')

    def test_str_says_which_way_a_disk_turns(self):
        self.assertIn('clockwise', str(self.settings['bw2']))       # speed is negative
        self.assertIn('anticlockwise', str(self.settings['bw1']))


if __name__ == '__main__':
    unittest.main()
