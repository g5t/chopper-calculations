import unittest


class BifrostTestCase(unittest.TestCase):
    def test_bifrost(self):
        from chopcal import bifrost
        d = bifrost(7.0, 0, 0.0002)
        for n, s in (('ps1', 14*14), ('ps2', 14*14), ('fo1', 14), ('fo2', 14), ('bw1', 14), ('bw2', -14)):
            self.assertTrue(n in d)
            print(f'{n}: {d[n]}')
            self.assertTrue(hasattr(d[n], 'speed'))
            self.assertTrue(hasattr(d[n], 'delay'))
            self.assertTrue(hasattr(d[n], 'angle'))
            self.assertTrue(hasattr(d[n], 'path'))
            self.assertEqual(d[n].speed, s)

    def test_delays_are_times_within_a_source_period(self):
        """A delay says when an opening is on the beam, so it is a time, in seconds.

        Every BIFROST chopper opens within one 14 Hz source period of t=0 -- the pulse
        shaping pair almost immediately, the bandwidth pair after the 78 m flight.
        """
        from chopcal import bifrost
        d = bifrost(7.0, 0, 0.0002)
        for name, chopper in d.items():
            self.assertGreater(chopper.delay, 0.0, name)
            self.assertLess(chopper.delay, 1 / 14, name)

        # further down the guide is later, whatever each disk's own speed happens to be
        self.assertLess(d['ps1'].delay, d['fo1'].delay)
        self.assertLess(d['fo1'].delay, d['fo2'].delay)
        self.assertLess(d['fo2'].delay, d['bw1'].delay)

    def test_the_pulse_shaping_pair_bursts_for_the_requested_time(self):
        """Two disks turning together, one lagging, cut a burst as short as the lag.

        They share a speed and a direction, so each is open for one slit crossing and
        both are open only for the part of it the lag leaves. Lagging by the crossing
        less the requested opening leaves exactly the requested opening -- and, because
        a delay is a time, it stays exact when the pair changes speed.
        """
        from chopcal import bifrost

        for shaping_time in (0.0002, 0.001, 0.005, 0.01):
            d = bifrost(7.0, 0, shaping_time)
            first, second = d['ps1'], d['ps2']
            self.assertEqual(first.speed, second.speed)  # co-rotating, not opposed

            # each disk is open for one slit crossing, centred on its own delay
            crossing = first.angle / 360.0 / abs(first.speed)
            opens = [(c.delay - crossing / 2, c.delay + crossing / 2)
                     for c in (first, second)]
            start = max(o[0] for o in opens)
            end = min(o[1] for o in opens)

            self.assertAlmostEqual(end - start, shaping_time, places=15,
                                   msg=f'burst at {shaping_time}')

    def test_the_burst_stays_centred_on_the_band_at_any_opening(self):
        """Widening the burst opens it symmetrically about the band's arrival.

        Whatever the requested opening, the pair is set so the part they are both open
        for straddles the same instant -- the middle of the band reaching them.
        """
        from chopcal import bifrost

        centres = []
        for shaping_time in (0.0002, 0.001, 0.002):
            d = bifrost(7.0, 0, shaping_time)
            crossing = d['ps1'].angle / 360.0 / abs(d['ps1'].speed)
            opens = [(c.delay - crossing / 2, c.delay + crossing / 2)
                     for c in (d['ps1'], d['ps2'])]
            centres.append((max(o[0] for o in opens) + min(o[1] for o in opens)) / 2)

        for centre in centres[1:]:
            self.assertAlmostEqual(centre, centres[0], places=15)

    def test_the_bandwidth_pair_shares_a_delay_while_counter_rotating(self):
        """Opposite rotation, same timing.

        Both disks are on the beam at the middle of the band; which way each turns to
        get there does not change when it is there. A delay can say that directly, where
        a phase relied on the sign being discarded further down.
        """
        from chopcal import bifrost
        d = bifrost(7.0, 0, 0.0002)

        self.assertEqual(d['bw1'].speed, -d['bw2'].speed)
        self.assertAlmostEqual(d['bw1'].delay, d['bw2'].delay, places=15)

    def test_a_delay_does_not_move_when_the_disks_speed_up(self):
        """The same opening at a different speed is the same delay.

        This is what a phase could not say: an angle has to be divided by the speed to
        become a time, so holding a phase fixed moved the opening. A shaping time too
        long for the pulse shaping pair to cut at 196 Hz drops them to a lower harmonic;
        the frame overlap and bandwidth disks are untouched and keep their timings.
        """
        from chopcal import bifrost
        slow = bifrost(7.0, 0, 0.01)    # too long for 196 Hz, so the pair runs at 42
        fast = bifrost(7.0, 0, 0.005)   # ...and at 84

        self.assertNotEqual(slow['ps1'].speed, fast['ps1'].speed)
        for name in ('fo1', 'fo2', 'bw1', 'bw2'):
            self.assertEqual(slow[name].speed, fast[name].speed, name)
            self.assertAlmostEqual(slow[name].delay, fast[name].delay, places=15, msg=name)


if __name__ == '__main__':
    unittest.main()
