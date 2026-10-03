"""Checks of S9f.2a's parallel walls in `curved_boolean_reference.py` and
of its fixture list (`generate_spline_parallel_boolean_fixtures.py`): the
exact map between parallel frames, a spline span's crossings with a circle
(closed forms, a tangency's multiplicity) and with another span (crossing
and touching parabolas), containment in closed form both ways, S9a.2's
`SplinePair` on a quarter-turned pair, and the declared degeneracies'
margins."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import curved_boolean_reference as ref
import generate_spline_parallel_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import disc


def dome_ctrl():
    """The dome `y = x (4 - x) / 2` over [0, 4], run from (4, 0)."""
    return ((F(4), F(0)), (F(2), F(4)), (F(0), F(0)))


class SplineParallelBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_the_map_between_parallel_frames_is_exact(self):
        A = ref.Prism(fixtures.prism([disc(0.0, 0.0, 1.0)], fixtures.at('XY', (0, 0, 0)), 0.0, 1.0, 91))
        B = ref.Prism(fixtures.prism([disc(0.0, 0.0, 1.0)], fixtures.at('TURN', (2, 3, 0)), 0.0, 1.0, 92))
        m = ref.parallel_map(A, B)
        # A quarter turn: (u', v') to (2 - v', 3 + u').
        self.assertEqual(m, ((2, 0, -1), (3, 1, 0)))
        inv = ref.invert_map(m)
        for p in ((F(1, 3), F(-5, 7)), (F(4), F(9, 2))):
            self.assertEqual(ref.map_point(inv, ref.map_point(m, p)), p)
        # R125's rounded rotation: exact rationals, not a quarter turn.
        C = ref.Prism(fixtures.prism([disc(0.0, 0.0, 1.0)], fixtures.at('R125', (0, 0, 0)), 0.0, 1.0, 92))
        m = ref.parallel_map(A, C)
        self.assertNotIn(m[0][1], (0, 1, -1))
        self.assertEqual(ref.map_point(ref.invert_map(m), ref.map_point(m, (F(1), F(2)))), (F(1), F(2)))

    def test_a_span_meets_a_circle_at_its_exact_roots(self):
        # The dome against the circle about (2, 1) of radius 5/4: (x - 2)^2
        # = 3/2, y = 5/4; tau = (4 - x) / 4.
        rd = ref.Round((F(2), F(1)), F(5, 4), True, True, True)
        got = ref.span_circle(dome_ctrl(), rd, mp.mpf(10))
        self.assertEqual(len(got), 2)
        want = sorted([(2-mp.sqrt(mp.mpf(3)/2))/4, (2+mp.sqrt(mp.mpf(3)/2))/4])
        for (tau, phi, sine, mult), w in zip(sorted(got), want):
            self.assertEqual(mult, 1)
            self.assertLess(abs(tau-w), mp.mpf(10)**-38)
            self.assertLess(abs(mp.sin(phi)-mp.mpf(1)/5), mp.mpf(10)**-38)
            self.assertGreater(sine, mp.mpf(1)/10)
        # The circle about (2, 4) of radius 2 touches the apex: a double
        # root at tau = 1/2.
        rd = ref.Round((F(2), F(4)), F(2), True, True, True)
        got = ref.span_circle(dome_ctrl(), rd, mp.mpf(10))
        self.assertEqual(len(got), 1)
        tau, _, sine, mult = got[0]
        self.assertEqual(mult, 2)
        self.assertLess(abs(tau-mp.mpf(1)/2), mp.mpf(10)**-38)
        self.assertLess(sine, mp.mpf(10)**-30)
        # A near miss: the circle about (2, 4.5) of radius 2 passes 1/2
        # above the apex.
        rd = ref.Round((F(2), F(9, 2)), F(2), True, True, True)
        self.assertEqual(ref.span_circle(dome_ctrl(), rd, mp.mpf(10)), [])
        near = ref.span_circle_near(dome_ctrl(), rd, mp.mpf(10))
        self.assertLess(abs(min(d for _, d in near)-mp.mpf(1)/2), mp.mpf(10)**-38)

    def test_two_spans_cross_and_touch_by_subdivision_and_newton(self):
        P = ref.Curve2(dome_ctrl())
        # y = 1 + (x - 2)^2 / 2: crosses the dome at (1, 3/2) and (3, 3/2).
        Q = ref.Curve2(((F(0), F(3)), (F(2), F(-1)), (F(4), F(3))))
        got = sorted(ref.span_crossings(P, Q))
        self.assertEqual(len(got), 2)
        for (s, u, sine, dist), x in zip(got, (3, 1)):
            px, py = P.point(s)
            self.assertLess(abs(px-x)+abs(py-mp.mpf(3)/2), mp.mpf(10)**-38)
            self.assertLess(dist, mp.mpf(10)**-38)
            self.assertGreater(sine, mp.mpf(1)/2)
        # y = 2 + (x - 2)^2 / 2 touches it at the apex (Newton on the
        # distance's gradient).
        R = ref.Curve2(((F(0), F(4)), (F(2), F(0)), (F(4), F(4))))
        got = ref.span_crossings(P, R)
        self.assertEqual(len(got), 1)
        s, u, sine, dist = got[0]
        self.assertLess(abs(s-mp.mpf(1)/2)+abs(u-mp.mpf(1)/2), mp.mpf(10)**-18)
        self.assertLess(sine, mp.mpf(10)**-15)
        # One raised by 1/4 misses it by 1/4.
        S = ref.Curve2(((F(0), F(17, 4)), (F(2), F(1, 4)), (F(4), F(17, 4))))
        self.assertEqual(ref.span_crossings(P, S), [])
        near = ref.span_near_misses(P, S, F(3, 10))
        self.assertLess(abs(min(d for _, _, d in near)-mp.mpf(1)/4), mp.mpf(10)**-38)

    def test_containment_in_closed_form_both_ways(self):
        # The dome prism (area 16/3, height 5) inside a quarter-turned disc
        # prism spanning its heights: the common is the dome's prism, the
        # cut empty, the fuse the disc's prism; the divergence theorem
        # agrees.
        p = fixtures.profiles()
        A = fixtures.prism(p['dome'], fixtures.at('XY', (0, 0, 0)), 0.0, 5.0, 91)
        B = fixtures.prism([disc(0.0, 0.0, 3.0)], fixtures.at('TURN', (2, 1, -1)), 0.0, 7.0, 92)
        pair = ref.Pair(A, B)
        vols = pair.volumes()
        self.assertLess(abs(vols['common'][0]-mp.mpf(80)/3), mp.mpf(10)**-35)
        self.assertLess(abs(vols['cut'][0]), mp.mpf(10)**-35)
        self.assertLess(abs(vols['fuse'][0]-63*mp.pi), mp.mpf(10)**-35)
        div = pair.divergence_volumes()
        for op in ref.OPS:
            self.assertLess(abs(div[op][0]-vols[op][0]), mp.mpf(10)**-35)
        self.assertEqual(pair.slicing.cross2d, [])

    def test_a_quarter_turned_pair_is_s9a2s(self):
        # bulge_dome_turn's common two ways: this reference and S9a.2's
        # SplinePair with the dome turned into the bulge's frame exactly.
        c = [x for x in fixtures.cases() if x.name == 'bulge_dome_turn_common'][0]
        L, offset = fixtures.exact_map(c.obj, c.tool)
        self.assertEqual(L, ((0, -1), (1, 0)))
        n, V, S, C = fixtures.s9a2_result(c.obj, c.operation, c.tool)
        got = ref.Pair(c.obj, c.tool).result('common')
        self.assertEqual(n, got[0])
        self.assertLess(abs(got[1]-V), mp.mpf(10)**-30)
        self.assertLess(abs(got[2]-S), mp.mpf(10)**-30)
        self.assertLess(max(abs(got[3][i]-C[i]) for i in range(3)), mp.mpf(10)**-30)

    def test_declared_degenerate_pairs_fail_their_margin(self):
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        for name, reason in (('dome_disc_touch', fixtures.TOUCH_CYLINDER),
                             ('dome_flip_touch', fixtures.TOUCH_SPLINE)):
            c = firsts[name]
            self.assertEqual(c.reason, reason)
            margin = fixtures.curved_margins(ref.Pair(c.obj, c.tool))
            self.assertLess(margin[fixtures.DEGENERATE_MARGIN[reason]], mp.mpf(10)**-12)
        c = firsts['bulge_disc_r125']
        margin = fixtures.margins(ref.Pair(c.obj, c.tool))
        self.assertTrue(all(v > fixtures.MARGIN for v in margin.values()))

    def test_cases_cover_the_decisions(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual({c.operation for c in listed}, {'fuse', 'cut', 'common'})
        reasons = {c.reason for c in listed if c.kind == 'degenerate'}
        self.assertEqual(reasons, {fixtures.TOUCH_CYLINDER, fixtures.TOUCH_SPLINE})
        frames = {fixtures.frame_name(f) for c in listed for f in c.frames}
        self.assertTrue({'XY', 'TURN', 'FLIP', 'R125', 'TILT', 'TILT2'} <= frames)
        # Splines against splines, and against arcs or circles.
        both = [c for c in listed if fixtures.has_spline(c.obj) and fixtures.has_spline(c.tool)]
        self.assertTrue(len({c.pair_name for c in both}) >= 6)
        # Same-axis pairs whose offsets round.
        offsets = [c for c in listed if c.pair_name.endswith('_offset')]
        self.assertEqual(len({c.pair_name for c in offsets}), 2)
        for c in offsets:
            self.assertTrue(any(F(float(v)) != v for v in fixtures.s9f1.exact_offset(c.obj, c.tool)))


if __name__ == '__main__':
    unittest.main()
