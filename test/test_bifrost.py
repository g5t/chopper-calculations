import unittest


class BifrostTestCase(unittest.TestCase):
    def test_bifrost(self):
        from chopcal import bifrost
        d = bifrost(7.0, 0, 0.0002)
        for n, s in (('ps1', 14*14), ('ps2', 14*14), ('fo1', 14), ('fo2', 14), ('bw1', 14), ('bw2', -14)):
            self.assertTrue(n+'speed' in d)
            self.assertTrue(n+'phase' in d)
            self.assertTrue(n+'delay' in d)
            self.assertEqual(d[n+'speed'], s)


if __name__ == '__main__':
    unittest.main()
