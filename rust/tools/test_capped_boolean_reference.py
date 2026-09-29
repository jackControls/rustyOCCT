"""Independent checks of S9c.2b.2's closed forms, cap-circle test and
fixture list (`generate_capped_boolean_fixtures.py`), without the slicing
reference."""
import unittest

import mpmath as mp

import generate_capped_boolean_fixtures as fixtures
import generate_procedural_boolean_fixtures as procedural
from generate_curved_boolean_fixtures import at, disc, group


def pair(obj, tool):
    return group('probe', obj, tool, {'common': 'solid'})[0]


class CappedBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_perpendicular_pair_crossing_whole_is_legendre(self):
        P = fixtures.perpendicular(pair(([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -50.0, 50.0),
                                        ([disc(0.0, 0.0, 1.2)], at('SIDE', (-50, 0, 0)), 0.0, 100.0)))
        volume, moments, _ = P.volume()
        want = procedural.legendre(2, 1.2)
        self.assertLess(abs(volume-want), mp.mpf(10)**-35*want)
        self.assertLess(max(abs(m) for m in moments), mp.mpf(10)**-35*want)

    def test_perpendicular_pair_is_the_same_in_ideal_turned_frames(self):
        # rim_pipe turned rigidly: TILT's axis for XY's, TILTX's for SIDE's
        # (TILTX's x along x = TILT n x TILTX n, its y along TILT's axis).
        firsts = {c.pair_name: c for c in fixtures.cases()}
        exact = fixtures.perpendicular(firsts['rim_pipe'])
        turned = fixtures.perpendicular(pair(([disc(0.0, 0.0, 2.0)], at('TILT', (0, 0, 0)), -3.0, 3.0),
                                             ([disc(0.3, 2.9, 1.0)], at('TILTX', (0, 0, 0)), -5.0, 5.0)))
        a, b = exact.volume()[0], turned.volume()[0]
        self.assertLess(abs(a-b), mp.mpf(10)**-35*a)
        self.assertLess(max(abs(x-y) for x, y in zip(exact.areas(), turned.areas())), mp.mpf(10)**-30)

    def test_segment_and_band(self):
        r = mp.mpf(0.8)
        area, moment = fixtures.segment(r, 0)
        self.assertLess(abs(area-mp.pi*r*r/2), mp.mpf(10)**-35)
        self.assertLess(abs(moment+2*r**3/3), mp.mpf(10)**-35)
        self.assertLess(abs(fixtures.segment(r, r)[0]-mp.pi*r*r), mp.mpf(10)**-35)
        self.assertLess(abs(fixtures.band(2, -2, 2)-4*mp.pi), mp.mpf(10)**-35)
        self.assertLess(abs(fixtures.band(2, 0, 2)-2*mp.pi), mp.mpf(10)**-35)

    def test_hole_rim_forms_add_up(self):
        forms = fixtures.hole_rim_forms()
        box = procedural.box_with_hole(10, 4, 5, 5, 2)
        self.assertLess(abs(forms['cut'][0]+forms['common'][0]-box[0]), mp.mpf(10)**-30)
        self.assertGreater(forms['common'][0], 0)
        # The pipe's part below the top face, between x = 0 and x = 10, less
        # the hole's cylinder: more than the part outside the hole's strip.
        area, _ = fixtures.segment(0.8, mp.mpf(0.1))
        self.assertLess(forms['common'][0], 10*area)
        self.assertGreater(forms['common'][0], 6*area)

    def test_cap_circles_crossing_the_other_face_are_found(self):
        firsts = {c.pair_name: c for c in fixtures.cases()}
        crossings, _ = fixtures.cap_crossings(firsts['rim_pipe'])
        self.assertEqual(fixtures.crossed(crossings), {'obj 3': 4})
        # The turned generator's pipe ending inside: its cap's circle clear
        # of the other wall, no crossing on a face.
        blind = pair(([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -3.0, 3.0),
                     ([disc(-0.4, 0.3, 1.0)], at('SIDE', (-6, 0, 0)), 0.0, 6.5))
        crossings, clear = fixtures.cap_crossings(blind)
        self.assertEqual(fixtures.crossed(crossings), {})
        self.assertGreater(clear, 1e-3)

    def test_cases_cross_a_cap_circle_on_the_other_face(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertTrue(12 <= len(listed) <= 18)
        firsts = {}
        for c in listed:
            firsts.setdefault(c.pair_name, c)
        exact = {name for name, c in firsts.items() if fixtures.exact_pair(c)}
        self.assertTrue(exact and exact != set(firsts))
        for name, c in firsts.items():
            crossings, _ = fixtures.cap_crossings(c)
            self.assertTrue(fixtures.crossed(crossings), name)
            for _, _, _, margin, slope in crossings:
                self.assertGreater(abs(margin), fixtures.FACE_MARGIN)
                self.assertGreater(slope, fixtures.SLOPE)


if __name__ == '__main__':
    unittest.main()
