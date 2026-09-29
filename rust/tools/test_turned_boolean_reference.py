"""Independent checks of S9c.2b.1's closed forms, geometry checks and
fixture list (`generate_turned_boolean_fixtures.py`), without the slicing
reference."""
import unittest

import mpmath as mp

import generate_procedural_boolean_fixtures as procedural
import generate_turned_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, disc, group


def pair(obj, tool):
    return group('probe', obj, tool, {'common': 'solid'})[0]


class TurnedBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_crossing_common_is_the_perpendicular_one_over_sin_phi(self):
        s = mp.mpf(3)/5
        a, b = (0, 0, 1), (0, s, mp.mpf(4)/5)
        volume, moments = fixtures.crossing_common((0, 0, 0), a, 2, (0, 0, 0), b, 1.2)
        self.assertLess(abs(volume-procedural.legendre(2, 1.2)/s), mp.mpf(10)**-35*volume)
        self.assertLess(max(abs(m) for m in moments), mp.mpf(10)**-35*volume)

    def test_crossing_common_centre_lies_on_the_common_perpendicular(self):
        # Offset along e = x by 0.5, B's axis shifted along itself: the
        # centre stays on the x axis.
        b = (0, mp.mpf(3)/5, mp.mpf(4)/5)
        pb = (mp.mpf(0.5), 3*b[1], 3*b[2])
        volume, moments = fixtures.crossing_common((0, 0, 0), (0, 0, 1), 2, pb, b, 1)
        centre = [m/volume for m in moments]
        self.assertLess(abs(centre[1])+abs(centre[2]), mp.mpf(10)**-35)
        self.assertGreater(centre[0], 0)

    def test_equal_radii_node_models(self):
        # Ideal intervals equal; the stored models' ends apart within the
        # resolution for the declared nodes, equal exactly where a stored
        # axis lies along e (a double tangency: two conics).
        firsts = {c.pair_name: c for c in fixtures.cases()}
        for name in ('node_lean', 'node_r125'):
            (kind, gap), (mkind, mgap) = map(fixtures.classify, fixtures.intervals(firsts[name]))
            self.assertLess(gap, mp.mpf(10)**-30)
            self.assertTrue(0 < mgap < mp.mpf(10)**-15)
        exact = pair(([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -7.0, 7.0),
                     ([disc(0.0, 0.0, 2.0)], at('TILT', (0, 0, 0)), -9.0, 9.0))
        _, model = fixtures.intervals(exact)
        self.assertEqual(fixtures.classify(model)[1], 0)

    def test_cap_circle_on_the_other_face_is_found(self):
        # A pipe of radius 1.5 ending on a thin cylinder's axis: its cap's
        # circle crosses the thin wall at z = +-sqrt(1.25), on its face.
        crossing = pair(([disc(0.0, 0.0, 1.0)], at('XY', (0, 0, 0)), -5.0, 5.0),
                        ([disc(0.0, 0.0, 1.5)], at('SIDE', (-6, 0, 0)), 0.0, 6.0))
        out, _, count = fixtures.cap_circles(crossing)
        self.assertEqual(count, 4)
        self.assertLess(out, 0)
        clear = fixtures.cap_circles(pair(([disc(0.0, 0.0, 1.0)], at('XY', (0, 0, 0)), -5.0, 5.0),
                                          ([disc(0.0, 0.0, 0.5)], at('TILT', (0, 0, 0)), -3.0, 3.0)))
        self.assertEqual(clear[2], 0)

    def test_cases_are_turned_pairs_of_cylinders(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertTrue(12 <= len(listed) <= 16)
        self.assertEqual(sum(1 for c in listed if c.kind == 'degenerate'), 2)
        for c in listed:
            if c.kind == 'degenerate':
                self.assertEqual(c.reason, fixtures.NODE)
        # The thin pipe is the object of one pair (its cut two stubs).
        self.assertTrue(any(c.obj.boundaries[0].circle[2] < c.tool.boundaries[0].circle[2] for c in listed))


if __name__ == '__main__':
    unittest.main()
