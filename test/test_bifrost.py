import unittest


class BifrostTestCase(unittest.TestCase):
    def test_bifrost(self):
        from chopcal import bifrost
        d = bifrost(7.0, 0, 0.0002)
        for n, s in (('ps1', 14*14), ('ps2', 14*14), ('fo1', 14), ('fo2', 14), ('bw1', 14), ('bw2', -14)):
            self.assertTrue(n in d)
            print(f'{n}: {d[n]}')
            self.assertTrue(hasattr(d[n], 'speed'))
            self.assertTrue(hasattr(d[n], 'phase'))
            self.assertTrue(hasattr(d[n], 'angle'))
            self.assertTrue(hasattr(d[n], 'path'))
            self.assertEqual(d[n].speed, s)


if __name__ == '__main__':
    unittest.main()
