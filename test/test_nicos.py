import unittest


class NicosTestCase(unittest.TestCase):
    def test_the_opening_is_centred_on_the_beam_after_the_two_delays(self):
        """ESS centres the reference opening on the beam pulse_delay + experiment_delay
        after the reference. chopper-lib puts a disc point at angle a on the beam at
        delay + (beam - a) / (360 speed), every turn."""
        from chopcal import bifrost
        from chopcal.nicos import experiment_delay
        pulse = 1e-5
        for name, chopper in bifrost(7.0, 0, 0.0002).items():
            period = 1 / abs(chopper.speed)
            given = experiment_delay(chopper, pulse)
            self.assertGreaterEqual(given, 0.0, name)
            self.assertLess(given, period, name)
            centre = (chopper.edges[0] + chopper.edges[1]) / 2
            at_beam = chopper.delay + (chopper.beam - centre) / (360 * chopper.speed)
            turns = (pulse + given - at_beam) / period
            self.assertAlmostEqual(turns, round(turns), places=9, msg=name)

    def test_the_counter_rotating_disc_is_centred_at_the_same_time(self):
        """bw1 and bw2 share a delay, the middle of their opening, and turn opposite
        ways. Timing the centre rather than an edge, they share an experiment delay too."""
        from chopcal import bifrost
        from chopcal.nicos import experiment_delay
        d = bifrost(7.0, 0, 0.0002)
        self.assertEqual(d['bw1'].speed, -d['bw2'].speed)
        self.assertAlmostEqual(experiment_delay(d['bw1']), experiment_delay(d['bw2']),
                               places=12)

    def test_nicos_gets_every_disc_by_device_in_its_own_units(self):
        from chopcal import bifrost as plain
        from chopcal.nicos import bifrost, experiment_delay
        values = bifrost(pulse_delay=10_000, energy_min=7.0)
        settings = plain(energy_min=7.0)
        self.assertEqual(sorted(values), sorted(
            f'{disc}_chopper_{what}' for disc in ('psc1', 'psc2', 'foc1', 'foc2', 'bwc1', 'bwc2')
            for what in ('speed', 'delay')))
        # Spd_S is positive clockwise facing away from the source, chopcal anticlockwise
        for key, disc in (('ps1', 'psc1'), ('bw1', 'bwc1'), ('bw2', 'bwc2')):
            self.assertEqual(values[f'{disc}_chopper_speed'], -settings[key].speed)
        self.assertEqual(values['bwc2_chopper_speed'], 14)
        self.assertAlmostEqual(values['psc1_chopper_delay'],
                               1e9 * experiment_delay(settings['ps1'], 1e-5), places=3)

    def test_a_disc_of_several_openings_is_refused(self):
        from chopcal.lib import Chopper
        from chopcal.nicos import experiment_delay
        chopper = Chopper(14.0, 0.01, 0.0, [0.0, 10.0, 180.0, 190.0], 10.0, 0.0)
        with self.assertRaises(ValueError):
            experiment_delay(chopper)

    def test_nicos_passes_on_the_disk_sets_speed_limit(self):
        from chopcal.nicos import bifrost
        values = bifrost(pulse_delay=10_000, energy_min=7.0, pulse_shaping_speed_max=126)
        self.assertEqual(values['psc1_chopper_speed'], -126)
        self.assertEqual(values['psc2_chopper_speed'], -126)
