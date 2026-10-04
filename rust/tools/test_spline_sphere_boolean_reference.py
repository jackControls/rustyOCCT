"""Checks of S9f.3a's reference (`spline_sphere_boolean_reference.py`) and of
its fixture list (`generate_spline_sphere_boolean_fixtures.py`): a sphere
across a straight spline wall (its common a sphere less a cap, its areas
and moments in closed form, its turning points on the wall's two spans),
a hemisphere's rim crossing the dome's generatrix in a tower field in
closed form, the declared degeneracies' margins and the cases' classes."""
import unittest

import mpmath as mp

from identity_reference import Case
import spline_sphere_boolean_reference as ref
import generate_spline_sphere_boolean_fixtures as fixtures


def straight_and_ball(centre, r, low=-ref.HALF_PI, high=ref.HALF_PI, frame='XY'):
    """S9a.2's square of straight splines `[0, 10]^2` (its right side a
    cubic with an interior knot at `y = 5`) over `[0, 5]`, and a sphere."""
    p = fixtures.profiles()
    S = fixtures.make(fixtures.prism(p['straight'], fixtures.at('XY', (0, 0, 0)), 0.0, 5.0), 91)
    B = Case('', 1e-7, 92, fixtures.at(frame, centre), 0.0, 0.0, [], sphere=(r, low, high))
    return S, B


class SplineSphereBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_a_sphere_across_a_straight_spline_wall_in_closed_form(self):
        # The sphere about (9.5, 5, 2.5) of radius 1 crosses the right side
        # x = 10 (a cubic with collinear poles): the common is the sphere
        # less the cap of height h = 1/2 beyond the wall.
        S, B = straight_and_ball((9.5, 5.0, 2.5), 1.0)
        pair = ref.Pair(S, B)
        pi, h, r = mp.pi, mp.mpf(1)/2, mp.mpf(1)
        cap = pi*h*h*(3*r-h)/3
        V = 4*pi/3-cap
        # The cap's centroid lies at r - ... from the centre: 3 (2r - h)^2 /
        # (4 (3r - h)) beyond it.
        zc = 3*(2*r-h)**2/(4*(3*r-h))
        mx = mp.mpf('9.5')*4*pi/3-(mp.mpf('9.5')+zc)*cap
        for way in (pair.first, pair.second):
            got = pair.volumes(way)['common']
            self.assertLess(abs(got[0]-V), mp.mpf(10)**-35)
            want = (mx, 5*V, mp.mpf('2.5')*V)
            self.assertLess(max(abs(x-y) for x, y in zip(got[1], want)), mp.mpf(10)**-34)
        # Areas: the sphere less its cap and the wall's disc.
        area = 4*pi-2*pi*r*h+pi*(r*r-(r-h)**2)
        self.assertLess(abs(pair.area('common')-area), mp.mpf(10)**-34)
        sphere_in = [c for _, n, c, _ in pair.faces() if n == 'sphere'][0]['in']
        self.assertLess(abs(sphere_in-(4*pi-2*pi*r*h)), mp.mpf(10)**-34)
        self.assertEqual(pair.solids(), {'fuse': 1, 'cut': 1, 'common': 1})
        # The turning points: where the ruling at y is tangent to the
        # sphere, (y - 5)^2 = 3/4, one on each of the cubic's spans.
        turns = []
        for el in pair.S.elements:
            A, Bq, C = ref.wall_quadratic(pair.S, el, pair.ball)
            D = ref.psub(ref.pmul(Bq, Bq), ref.pc(C, A))
            if ref.pdeg(D) > 0:
                turns += [el.point(t)[1] for t in ref.exact_roots(D, slack=mp.mpf(0)) if el.point(t)[0] == 10]
        want = sorted(5+s*mp.sqrt(3)/2 for s in (-1, 1))
        self.assertEqual(len(turns), 2)
        for x, y in zip(sorted(turns), want):
            self.assertLess(abs(x-y), mp.mpf(10)**-38)

    def test_a_hemispheres_rim_on_a_generatrix_in_closed_form(self):
        # A hemisphere on its side (x >= 5/2) about (5/2, 1, 3/2) of radius
        # 5/4: its rim's plane x = 5/2 holds the dome's ruling at t = 3/8 (y
        # = 15/8), met at z = 3/2 +- sqrt(25/16 - (7/8)^2).
        p = fixtures.profiles()
        S = fixtures.make(fixtures.prism(p['dome'], fixtures.at('XY', (0, 0, 0)), 0.0, 3.0), 91)
        B = Case('', 1e-7, 92, fixtures.at('SIDE', (2.5, 1.0, 1.5)), 0.0, 0.0, [], sphere=(1.25, 0.0, ref.HALF_PI))
        pair = ref.Pair(S, B)
        # (The plane meets the base's wall at y = 0 too.)
        got = [X for _, tag, X in pair.first.events() if tag == 'rim_wall' and X[1] > 1]
        self.assertEqual(len(got), 2)
        q = mp.mpf(25)/16-(mp.mpf(7)/8)**2
        for X, s in zip(sorted(got, key=lambda X: X[2]), (-1, 1)):
            self.assertLess(abs(X[0]-mp.mpf(5)/2)+abs(X[1]-mp.mpf(15)/8), mp.mpf(10)**-38)
            self.assertLess(abs(X[2]-(mp.mpf(3)/2+s*mp.sqrt(q))), mp.mpf(10)**-38)
        # The disc's classes sum to its area, the two slicings agree.
        one, two = pair.sliced(pair.first), pair.sliced(pair.second)
        for k in 'PBC':
            self.assertLess(abs(one[k][0]-two[k][0]), mp.mpf(10)**-34)
        disc = [c for _, n, c, _ in pair.faces() if n == 'disc'][0]
        self.assertLess(abs(disc['in']+disc['out']-pi_r2(1.25)), mp.mpf(10)**-34)

    def test_declared_degenerate_pairs_fail_their_margin(self):
        firsts = {}
        for c in fixtures.cases():
            firsts.setdefault(c.pair_name, c)
        for name, reason in (('dome_touch', fixtures.TOUCH), ('knot_turn', fixtures.KNOT_TURN),
                             ('cap_turn', fixtures.TURN_EDGE)):
            c = firsts[name]
            self.assertEqual(c.reason, reason)
            margin = fixtures.margins(ref.Pair(c.obj, c.tool))
            self.assertLess(margin[fixtures.DEGENERATE_MARGIN[reason]], mp.mpf(10)**-12)

    def test_cases_cover_the_decisions(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual({c.operation for c in listed}, {'fuse', 'cut', 'common'})
        self.assertEqual({c.reason for c in listed if c.kind == 'degenerate'},
                         {fixtures.TOUCH, fixtures.KNOT_TURN, fixtures.TURN_EDGE})
        frames = {fixtures.frame_name(x.frame) for c in listed for x in (c.obj, c.tool)}
        self.assertTrue({'XY', 'SIDE', 'TILT'} <= frames)
        # The sphere as the object; hemispheres on the axis and on the side.
        self.assertTrue(any(c.obj.sphere is not None for c in listed))
        halves = {c.pair_name for c in listed for x in (c.obj, c.tool)
                  if x.sphere is not None and x.sphere[1] == 0.0}
        self.assertEqual(halves, {'hemi_bulge', 'side_hemi_bulge'})


def pi_r2(r):
    return mp.pi*mp.mpf(r)**2


if __name__ == '__main__':
    unittest.main()
