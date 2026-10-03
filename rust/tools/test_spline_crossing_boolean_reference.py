"""Checks of S9f.2b's crossing walls in `curved_boolean_reference.py` and of
its fixture list (`generate_spline_crossing_boolean_fixtures.py`): a spline
wall's meeting with a cylinder on a crossing axis at the caps' heights and
its turning points against closed forms, the slicing and the divergence
theorem against the product of chords on a perpendicular pair, and the
fixtures' classes and declared degeneracies' margins."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import curved_boolean_reference as ref
import generate_spline_crossing_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import disc


class SplineCrossingBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def dome_and_rod(self, y0, z0, r, h=3.0):
        """The dome `y = x (4 - x) / 2` (its span `(4 - 4 tau, 8 tau (1 -
        tau))`) in `XY` over `[0, h]`, and a rod along `x` (`SIDE`) about
        `(y0, z0)` of radius `r` from `x = -1` to `5`."""
        p = fixtures.profiles()
        S = fixtures.prism(p['dome'], fixtures.at('XY', (0, 0, 0)), 0.0, h, 91)
        C = fixtures.prism([disc(y0, z0, r)], fixtures.at('SIDE', (-1.0, 0.0, 0.0)), 0.0, 6.0, 92)
        A, B = ref.Prism(S), ref.Prism(C)
        span = [el for el in A.profile.elements if el.kind == 'spline'][0]
        return S, C, A, B, span, B.profile.elements[0]

    def test_the_meeting_at_a_cap_and_its_turning_points_in_closed_form(self):
        # Along the ruling at tau the rod's function is (y - y0)^2 + (w -
        # z0)^2 - r^2: A = 1, B = -z0, the discriminant r^2 - (y - y0)^2, so
        # the turning points are where 8 tau (1 - tau) = y0 +- r, and at a
        # cap's height h the meeting is where y = y0 +- sqrt(r^2 - (h -
        # z0)^2).
        y0, z0, r = 1.25, 1.0, 0.5
        _, _, A, B, span, circ = self.dome_and_rod(y0, z0, r)
        turns, D = ref.wall_cylinder_turns(A, span, B, circ)
        self.assertEqual(ref.pdeg(D), 4)
        want = sorted((1+s*mp.sqrt(1-mp.mpf(y)/2))/2 for y in (y0-r, y0+r) for s in (-1, 1))
        got = sorted(tau for tau, _ in turns)
        self.assertEqual(len(got), 4)
        for x, w in zip(got, want):
            self.assertLess(abs(x-w), mp.mpf(10)**-38)
        for tau, X in turns:
            self.assertLess(abs(X[2]-z0), mp.mpf(10)**-38)
        for h in (F(3, 4), F(5, 4)):
            q = mp.mpf(r)**2-(ref.M(h)-z0)**2
            ys = [mp.mpf(y0)-mp.sqrt(q), mp.mpf(y0)+mp.sqrt(q)]
            want = sorted((1+s*mp.sqrt(1-y/2))/2 for y in ys for s in (-1, 1))
            got = sorted(tau for tau, _ in ref.wall_cylinder_events(A, span, B, circ, A, h))
            self.assertEqual(len(got), 4)
            for x, w in zip(got, want):
                self.assertLess(abs(x-w), mp.mpf(10)**-38)
        # The rod's own cap at x = 5 misses the span (x = 4 - 4 tau); one
        # at x = 2 holds the dome's ruling at tau = 1/2 (its plane along the
        # wall's axis), met by a rod about (7/4, 1) at w = z0 +- sqrt(r^2 -
        # (2 - 7/4)^2) (the tower field's points, S9f.2b.2's).
        for x, count in ((F(5), 0), (F(2), 2)):
            C2 = fixtures.prism([disc(1.75, z0, r)], fixtures.at('SIDE', (-1.0, 0.0, 0.0)), 0.0, float(x+1), 92)
            B2 = ref.Prism(C2)
            got = ref.wall_cylinder_events(A, span, B2, B2.profile.elements[0], B2, B2.hi)
            self.assertEqual(len(got), count)
            for tau, X in got:
                self.assertLess(abs(tau-(4-ref.M(x))/4), mp.mpf(10)**-38)
                self.assertLess(abs((X[2]-z0)**2-(mp.mpf(r)**2-mp.mpf(1)/16)), mp.mpf(10)**-38)

    def test_a_perpendicular_common_is_the_product_of_chords(self):
        # The dome against a rod across its arch: the slicing, the
        # divergence theorem and the product L(y) W(y) agree.
        S, C, *_ = self.dome_and_rod(1.25, 1.0, 0.5)
        pair = ref.Pair(S, C)
        V, M = pair.volumes()['common']
        Vp, Mp = fixtures.product_common(S, C)
        self.assertLess(abs(V-Vp), mp.mpf(10)**-35)
        self.assertLess(max(abs(x-y) for x, y in zip(M, Mp)), mp.mpf(10)**-35)
        div = pair.divergence_volumes()
        for op in ref.OPS:
            self.assertLess(abs(div[op][0]-pair.volumes()[op][0]), mp.mpf(10)**-35)
        # The rod's turning points lie inside both faces: a loop on each
        # side of the arch (S9f.2b.2's).
        margin = fixtures.crossing_margins(pair)
        self.assertGreater(margin['loop'], fixtures.MARGIN)

    def test_declared_degenerate_pairs_fail_their_margin(self):
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        for name, reason in (('dome_touch', fixtures.TOUCH), ('knot_turn', fixtures.KNOT_TURN)):
            c = firsts[name]
            self.assertEqual(c.reason, reason)
            margin = fixtures.crossing_margins(ref.Pair(c.obj, c.tool))
            self.assertLess(margin[fixtures.DEGENERATE_MARGIN[reason]], mp.mpf(10)**-12)

    def test_cases_cover_the_decisions(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual({c.operation for c in listed}, {'fuse', 'cut', 'common'})
        self.assertEqual({c.reason for c in listed if c.kind == 'degenerate'},
                         {fixtures.TOUCH, fixtures.KNOT_TURN})
        self.assertEqual({c.step for c in listed}, {fixtures.STEP, fixtures.LOOPS})
        frames = {fixtures.frame_name(f) for c in listed for f in c.frames}
        self.assertTrue({'XY', 'SIDE', 'TILT', 'LEAN', 'STEEP'} <= frames)
        # The spline prism as the tool, a hole's cylinder, a stadium's arcs
        # and edges, a cylinder's cap inside the spline prism, knots crossed.
        names = {c.pair_name for c in listed}
        self.assertTrue({'steep_blob', 'ring_dome', 'bulge_stadium_tilt', 'blob_tilt_end', 'wave_lean',
                         'capsule_tilt'} <= names)
        self.assertTrue(any(fixtures.has_spline(c.tool) for c in listed))


if __name__ == '__main__':
    unittest.main()
