"""Checks of S9f.3b's reference (`spline_cone_boolean_reference.py`) and of
its fixture list (`generate_spline_cone_boolean_fixtures.py`): a frustum
halved by a straight spline wall through its axis (its common half the
frustum, its moments, areas and the wall's trapezoid in closed form), the
wall's area element's integral over a turn in closed form, `A`'s sign for
coaxial and crossing axes, a rim crossing the dome's generatrices in a
tower field in closed form, the declared degeneracies' margins and the
cases' classes."""
import unittest

import mpmath as mp

from identity_reference import Case
import spline_cone_boolean_reference as ref
import generate_spline_cone_boolean_fixtures as fixtures


def straight_and_cone(origin, bottom, top, height, frame='XY'):
    """S9a.2's square of straight splines `[0, 10]^2` (its right side a
    cubic with an interior knot at `y = 5`) over `[0, 5]`, and a cone."""
    p = fixtures.profiles()
    S = fixtures.make(fixtures.prism(p['straight'], fixtures.at('XY', (0, 0, 0)), 0.0, 5.0), 91)
    K = Case('', 1e-7, 92, fixtures.at(frame, origin), 0.0, 0.0, [], cone=(bottom, top, height))
    return S, K


class SplineConeBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_a_frustum_halved_by_a_straight_spline_wall_in_closed_form(self):
        # The frustum about (10, 5) from z = 1 to 4, radii 2 and 1: the
        # right side x = 10 (a cubic with collinear poles) holds its axis,
        # so the common is half of it.
        S, K = straight_and_cone((10.0, 5.0, 1.0), 2.0, 1.0, 3.0)
        pair = ref.Pair(S, K)
        self.assertEqual(pair.A, -mp.mpf(1)/9)
        pi = mp.pi
        b, t, h = mp.mpf(2), mp.mpf(1), mp.mpf(3)
        k = (t-b)/h
        rho2 = (t**3-b**3)/(3*k)
        rho3 = (t**4-b**4)/(4*k)
        V = pi*rho2/2
        # The half discs' x moments 2 rho^3 / 3 below the axis.
        mx = 10*V-2*rho3/3
        wbar = h*(b*b+2*b*t+3*t*t)/(4*(b*b+b*t+t*t))
        got = pair.volumes(pair.first)['common']
        self.assertLess(abs(got[0]-V), mp.mpf(10)**-35)
        want = (mx, 5*V, (1+wbar)*V)
        self.assertLess(max(abs(x-y) for x, y in zip(got[1], want)), mp.mpf(10)**-34)
        # Areas: half the lateral surface and the discs, and the trapezoid
        # of the axis's plane.
        lateral = pi*(b+t)*mp.sqrt(h*h+(t-b)**2)
        area = lateral/2+pi*(b*b+t*t)/2+(b+t)*h
        self.assertLess(abs(pair.area('common')-area), mp.mpf(10)**-34)
        wall_in = [c for tag, n, c, _ in pair.faces() if tag == 'K' and n == 'wall'][0]['in']
        self.assertLess(abs(wall_in-lateral/2), mp.mpf(10)**-34)
        self.assertEqual(pair.solids(), {'fuse': 1, 'cut': 1, 'common': 1})

    def test_the_wall_area_element_over_a_turn(self):
        # On a leaning slicing direction the element D(theta)^-2 integrates
        # over a turn to 2 pi |mu| / (mu^2 - R^2)^(3/2).
        S, K = straight_and_cone((5.0, 5.0, 1.0), 2.0, 1.0, 3.0)
        pair = ref.Pair(S, K)
        way = pair.second
        mu, R = way.mu, way.R
        self.assertLess(R, abs(mu))
        for theta in (mp.mpf(0), mp.mpf('0.3'), mp.mpf(-2)):
            got = way.inv_d2(theta+2*mp.pi)-way.inv_d2(theta)
            want = 2*mp.pi*abs(mu)/(mu*mu-R*R)**mp.mpf(1.5)
            self.assertLess(abs(got-want), mp.mpf(10)**-35)
        # Across a turn's end the antiderivative is continuous.
        psi = way.phi+mp.pi
        e = mp.mpf(10)**-20
        self.assertLess(abs(way.inv_d2(psi+e)-way.inv_d2(psi-e)), mp.mpf(10)**-15)

    def test_the_axis_terms_sign(self):
        # Coaxial: A = -k^2; across: A = 1 (the prism's axis z in the cone's
        # x row).
        S, K = straight_and_cone((5.0, 5.0, 1.0), 2.0, 1.0, 3.0)
        self.assertEqual(ref.Pair(S, K).A, -mp.mpf(1)/9)
        S, K = straight_and_cone((5.0, 2.0, 2.5), 1.0, 0.5, 3.0, frame='YAX')
        self.assertEqual(ref.Pair(S, K).A, 1)

    def test_a_rims_tower_points_in_closed_form(self):
        # `dome_side_cone`: the bottom rim, the circle about (17/8, 5/4, 11/4)
        # of radius 3/2 in the plane y = 5/4, meets the dome's generatrices
        # at x = 2 +- sqrt(3/2) at z = 11/4 +- sqrt(9/4 - (x - 17/8)^2).
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        c = firsts['dome_side_cone']
        pair = ref.Pair(c.obj, c.tool)
        got = sorted((X for _, tag, X in pair.first.events() if tag == 'rim_wall' and abs(X[1]-mp.mpf(5)/4) < 0.1),
                     key=lambda X: (X[0], X[2]))
        want = []
        for sx in (-1, 1):
            x = 2+sx*mp.sqrt(mp.mpf(3)/2)
            dz = mp.sqrt(mp.mpf(9)/4-(x-mp.mpf(17)/8)**2)
            want += [(x, mp.mpf(5)/4, mp.mpf(11)/4+sz*dz) for sz in (-1, 1)]
        self.assertEqual(len(got), 4)
        for X, Y in zip(got, want):
            self.assertLess(max(abs(x-y) for x, y in zip(X, Y)), mp.mpf(10)**-38)

    def test_declared_degenerate_pairs_fail_their_margin(self):
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        for name, reason in (('apex_wall', fixtures.APEX), ('ruling_tilt', fixtures.RULING),
                             ('dome_touch', fixtures.TOUCH)):
            c = firsts[name]
            self.assertEqual(c.reason, reason)
            margin = fixtures.margins(ref.Pair(c.obj, c.tool))
            self.assertLess(margin[fixtures.DEGENERATE_MARGIN[reason]], mp.mpf(10)**-12)

    def test_cases_cover_the_decisions(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual({c.operation for c in listed}, {'fuse', 'cut', 'common'})
        self.assertEqual({c.reason for c in listed if c.kind == 'degenerate'},
                         {fixtures.APEX, fixtures.RULING, fixtures.TOUCH})
        frames = {fixtures.frame_name(x.frame) for c in listed for x in (c.obj, c.tool)}
        self.assertTrue({'XY', 'DOWN', 'SIDE', 'YAX', 'TILT'} <= frames)
        # The cone as the object; cones and frustums; both signs of A.
        self.assertTrue(any(c.obj.cone is not None for c in listed))
        cones = [x.cone for c in listed for x in (c.obj, c.tool) if x.cone is not None]
        self.assertTrue(any(0.0 in k[:2] for k in cones) and any(0.0 not in k[:2] for k in cones))
        signs = set()
        for c in listed:
            if c.kind == 'solid':
                S, K = (c.obj, c.tool) if c.obj.cone is None else (c.tool, c.obj)
                signs.add(mp.sign(ref.axis_term(ref.sref.Spline3(S), ref.Cone(K))))
        self.assertEqual(signs, {-1, 1})


if __name__ == '__main__':
    unittest.main()
