"""Checks of S9f.1's spline walls in `curved_boolean_reference.py` and of
its fixture list (`generate_spline_any_boolean_fixtures.py`): the span root
finder, the implicit crossings' events, spline profiles' Green moments, a
spline prism between two planes parallel to its axis in closed form, both
ways, and the declared degeneracies' margins."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import curved_boolean_reference as ref
import generate_spline_any_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, square


def poly(roots):
    out = [mp.mpf(1)]
    for r in roots:
        out = [(out[k-1] if k > 0 else 0)-r*(out[k] if k < len(out) else 0) for k in range(len(out)+1)]
    return out


class SplineAnyBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_span_roots_are_the_roots_on_the_span(self):
        roots = [mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(4)/5, mp.mpf(3)]
        got = ref.span_roots(poly(roots))
        self.assertEqual(len(got), 3)
        for g, r in zip(got, roots):
            self.assertLess(abs(g-r), mp.mpf(10)**-38)

    def test_a_touching_extreme_is_a_double_root_only_with_touch(self):
        # (t - 1/2)^2 + 1e-30: no real root, a double one within the touch.
        c = poly([mp.mpf(1)/2, mp.mpf(1)/2])
        c[0] += mp.mpf(10)**-30
        self.assertEqual(ref.span_roots(c), [])
        got = ref.span_roots(c, touch=mp.mpf(10)**-25)
        self.assertEqual(len(got), 2)
        self.assertLess(abs(got[0]-mp.mpf(1)/2), mp.mpf(10)**-30)

    def test_point_events_are_exact_roots_of_degree_p(self):
        # The dome y = x (4 - x) / 2 and the point (x, 3/2): x = 1 and 3.
        dome = ref.SplineSpan([(4, 0), (2, 4), (0, 0)], True, True)
        events = dome.point_events((F(0), F(3, 2)), (F(1), F(0)))
        xs = sorted(x for _, x in events)
        self.assertEqual(len(xs), 2)
        self.assertLess(abs(xs[0]-1)+abs(xs[1]-3), mp.mpf(10)**-35)
        # Horizontal tangent at the apex only.
        taus = dome.tangent_params((F(1), F(0)))
        self.assertEqual(len(taus), 1)
        self.assertLess(abs(taus[0]-mp.mpf(1)/2), mp.mpf(10)**-38)

    def test_spline_profile_moments_are_greens(self):
        # The dome: area 16/3, centroid (2, 4/5).
        p = fixtures.profiles()
        prof = ref.Profile(p['dome'], 1e-7)
        A, Mu, Mw = prof.moments()
        self.assertLess(abs(A-mp.mpf(16)/3), mp.mpf(10)**-38)
        self.assertLess(abs(Mu/A-2)+abs(Mw/A-mp.mpf(4)/5), mp.mpf(10)**-38)
        # Every profile (the capsule's arc among its spline and lines, the
        # lens hole) against S9a.2's Green moments over its own elements.
        import boolean_reference as s9a
        import generate_boolean_fixtures as s9a_fixtures
        for name, boundaries in p.items():
            got = ref.Profile(boundaries, 1e-7).moments()
            want = s9a_fixtures.green(s9a.profile_elements(boundaries, 1e-7))
            for x, y in zip(got, want):
                self.assertLess(abs(x-y), mp.mpf(10)**-36*max(1, abs(y)), name)

    def test_dome_between_side_planes_in_closed_form_both_ways(self):
        # The dome prism (height 5) and a box x in [1, 3], z in [1, 6]: the
        # common is the dome's region over [1, 3] (area 11/3) times 4.
        p = fixtures.profiles()
        A = fixtures.prism(p['dome'], at('XY', (0, 0, 0)), 0.0, 5.0, 91)
        B = fixtures.prism([square(-1.0, 1.0, 3.0, 6.0)], at('SIDE', (1, 0, 0)), 0.0, 2.0, 92)
        pair = ref.Pair(A, B)
        vols = pair.volumes()
        self.assertLess(abs(vols['common'][0]-mp.mpf(44)/3), mp.mpf(10)**-35)
        self.assertLess(abs(vols['cut'][0]-12), mp.mpf(10)**-35)
        div = pair.divergence_volumes()
        for op in ref.OPS:
            self.assertLess(abs(div[op][0]-vols[op][0]), mp.mpf(10)**-35)
            for x, y in zip(div[op][1], vols[op][1]):
                self.assertLess(abs(x-y), mp.mpf(10)**-34)
        self.assertEqual(pair.solids(), {'fuse': 1, 'cut': 1, 'common': 1})

    def test_declared_degenerate_pairs_fail_their_margin(self):
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        for name, reason in (('dome_tangent', fixtures.TANGENT_PLANE), ('kink_knot', fixtures.KNOT_TANGENT),
                             ('blob_rounding', fixtures.AXIS_ROUNDING)):
            c = firsts[name]
            self.assertEqual(c.reason, reason)
            margin = fixtures.margins(ref.Pair(c.obj, c.tool))
            self.assertLess(margin[fixtures.DEGENERATE_MARGIN[reason]], mp.mpf(10)**-12)
        margin = fixtures.margins(ref.Pair(firsts['blob_lean'].obj, firsts['blob_lean'].tool))
        self.assertTrue(all(v > fixtures.MARGIN for v in margin.values()))

    def test_kink_knot_has_multiplicity_p_and_stays_c1_lifted(self):
        from curve_surface_reference import stored_axes
        spline = fixtures.profiles()['kink'][0].segments[1]
        self.assertEqual(spline.mults[1], spline.degree)
        left, right = spline.pieces()
        d0 = tuple(left[-1][i]-left[-2][i] for i in range(2))
        d1 = tuple(right[1][i]-right[0][i] for i in range(2))
        self.assertEqual(d0, d1)
        # R4: its poles lifted into TILT (binary64, as the kernel's
        # extrusion and the native rows lift them) stay C1 at the knot; the
        # same knot along (-2, 1), poles (7, 3), (5, 4), (3, 5), would not
        # (the kernel's extrusion refuses it: InvalidTopology("edge_not_c1")).
        o, x, y, _ = stored_axes(at('TILT', (0, 0, 0)))

        def c1(poles):
            lift = [tuple(F(o[i]+x[i]*p[0]+y[i]*p[1]) for i in range(3)) for p in poles]
            return tuple((lift[0][i]+lift[2][i])/2 for i in range(3)) == lift[1]

        self.assertTrue(c1(spline.poles[1:4]))
        self.assertFalse(c1([(7.0, 3.0), (5.0, 4.0), (3.0, 5.0)]))

    def test_r4_fixtures_take_the_rounded_knot_and_exact_parallels(self):
        from curve_surface_reference import stored_axes
        listed = fixtures.r4_cases()
        fixtures.validate(listed)
        self.assertEqual({c.pair_name for c in listed}, {'knot_tilt', 'blob_side_turned'})
        self.assertFalse({c.name for c in listed} & {c.name for c in fixtures.cases()})
        # R4's knot: multiplicity two (the degree), C1 exactly in the
        # profile, off C1 by its lifted poles in TILT (the kernel removes
        # the knot once before lifting).
        spline = fixtures.profiles()['knot'][0].segments[2]
        self.assertEqual(spline.mults[1], spline.degree)
        left, right = spline.pieces()
        d0 = tuple(left[-1][i]-left[-2][i] for i in range(2))
        d1 = tuple(right[1][i]-right[0][i] for i in range(2))
        self.assertEqual(d0, d1)
        o, x, y, _ = stored_axes(at('TILT', (0, 0, 0)))
        lift = [tuple(F(o[i]+x[i]*p[0]+y[i]*p[1]) for i in range(3)) for p in spline.poles[1:4]]
        self.assertNotEqual(tuple((lift[0][i]+lift[2][i])/2 for i in range(3)), lift[1])
        # TILT's axis against the caps of SIDE and the x planes of XY:
        # exactly parallel (generatrices); against TILTX's caps 8.9e-17 off.
        n = tuple(F(c) for c in stored_axes(at('TILT', (0, 0, 0)))[3])

        def cap(name):
            _, cx, cy, _ = (tuple(F(c) for c in v) for v in stored_axes(at(name, (0, 0, 0))))
            return ref.cross(cx, cy)

        self.assertEqual(ref.dot(n, cap('SIDE')), 0)
        self.assertEqual(ref.dot(n, (F(1), F(0), F(0))), 0)
        self.assertLess(abs(float(ref.dot(n, cap('TILTX')))), 1e-16)
        self.assertNotEqual(ref.dot(n, cap('TILTX')), 0)

    def test_cases_cover_the_decisions(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertTrue(30 <= len(listed) <= 45)
        self.assertEqual({c.operation for c in listed}, {'fuse', 'cut', 'common'})
        reasons = {c.reason for c in listed if c.kind == 'degenerate'}
        self.assertEqual(reasons, {fixtures.TANGENT_PLANE, fixtures.KNOT_TANGENT, fixtures.AXIS_ROUNDING})
        frames = {fixtures.frame_name(f) for c in listed for f in c.frames}
        self.assertTrue({'R125', 'LEAN', 'TILT', 'TILT2', 'TILTX'} <= frames)
        # A same-axis pair whose offset rounds.
        offset = [c for c in listed if c.pair_name == 'bulge_offset'][0]
        self.assertTrue(any(F(float(v)) != v for v in fixtures.exact_offset(offset.obj, offset.tool)))


if __name__ == '__main__':
    unittest.main()
