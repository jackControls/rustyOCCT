"""Independent checks of S9d.3c's closed forms and two-way checks
(`generate_cones_loops_boolean_fixtures.py`): the hemisphere's sections
against a coaxial cone by polar quadrature, a cone wider than the ball
giving the hemisphere, an input's split parts on its own exact model, a
split against the whole on a simple pair, and the tangent cone's
construction."""
import math
import unittest
from fractions import Fraction as F

import mpmath as mp

import cones_boolean_reference as ref
import generate_cones_loops_boolean_fixtures as fixtures
from generate_cones_boolean_fixtures import Boolean, at, circle, cone, prism, sphere

HP = math.pi/2
O = (0, 0, 0)


class ConesLoopsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_segment_by_polar_quadrature(self):
        # A disc of radius m above the chord y = y0: the area swept by its
        # radii, r from the chord (or the centre) to the circle.
        for m, y0 in ((mp.mpf('1.5'), mp.mpf('0.5')), (mp.mpf(2), mp.mpf('-0.75'))):
            def ring(th):
                s = mp.sin(th)
                if y0 >= 0:
                    lo = y0/s if s > 0 else None
                    return (m*m-lo*lo)/2 if lo is not None and lo < m else mp.mpf(0)
                if s >= 0:
                    return m*m/2
                # Below the centre: r sin(th) >= y0 up to r = y0 / sin(th).
                return min(y0/s, m)**2/2
            ends = [mp.mpf(0), mp.pi, 2*mp.pi]
            if abs(y0) < m:
                a = mp.asin(abs(y0)/m)
                ends += [a, mp.pi-a, mp.pi+a, 2*mp.pi-a]
            want = mp.quad(ring, sorted(ends))
            self.close(fixtures.g3.segment(m, y0)[0], want, mp.mpf(10)**-25)
            self.close(fixtures.arc(m, y0), 2*mp.acos(y0/m))

    def test_a_wide_cone_holds_the_hemisphere(self):
        # A cone wider than the ball over its heights: the common is the
        # hemisphere, its measures in closed form.
        first = Boolean('w_common', 'common', sphere(2, at('TILT', O), 0.0, HP),
                        cone(10, 0, 20, at('XY', (0, 0, -5))), ('dome',))
        V, A, C = fixtures.dome_cone_forms(first)['common']
        R = mp.mpf(2)
        self.close(V, 2*mp.pi*R**3/3)
        self.close(A, 3*mp.pi*R*R)
        n = (0, mp.mpf(3)/5, mp.mpf(4)/5)
        for i in range(3):
            self.close(C[i], 3*R/8*n[i])

    def test_split_parts_are_the_inputs_model(self):
        # F cut at height 1: both parts on F's quadric exactly; a prism's
        # parts its model's heights.
        whole = ref.make_input(fixtures.g3.make(cone(2, 1, 2, at('XY', O)), 91), 'A')
        for part in fixtures.split_specs(cone(2, 1, 2, at('XY', O)), 1):
            x = ref.make_input(fixtures.g3.make(part, 91), 'A')
            self.assertEqual((x.surfaces[0].M, x.surfaces[0].m, x.surfaces[0].c),
                             (whole.surfaces[0].M, whole.surfaces[0].m, whole.surfaces[0].c))
        lo, hi = fixtures.split_specs(prism([circle(1, 0, 0.5)], at('LEAN', O), -3, 3), 0.5)
        self.assertEqual((lo[3], lo[4], hi[3], hi[4]), (-3, 0.5, 0.5, 3))

    def test_a_split_sums_to_the_whole(self):
        # F against a coaxial pipe, the pipe cut in two.
        first = Boolean('s_common', 'common', cone(2, 1, 2, at('XY', O)),
                        prism([circle(0, 0, 1.25)], at('XY', O), -1, 3), None)
        fixtures.SPLITS['s'] = ('tool', 1)
        try:
            self.assertLess(fixtures.splits(first), self.eps)
        finally:
            del fixtures.SPLITS['s']

    def test_the_tangent_cone_touches_the_rim(self):
        # Before rounding the cone's normal at the rim's point is normal to
        # the rim's tangent, and the point lies on the cone.
        (ax, ay, z0), r0, r1, h = fixtures.tangent_cone(2, -0.25, 0.75, -2.5, 4)
        p = (mp.mpf('1.2'), mp.mpf('1.28'), mp.mpf('-0.96'))
        T = (mp.mpf('-1.6'), mp.mpf('0.96'), mp.mpf('-0.72'))
        k = (mp.mpf(r1)-mp.mpf(r0))/mp.mpf(h)
        hx, hy = p[0]-mp.mpf(ax), p[1]-mp.mpf(ay)
        rho = mp.sqrt(hx*hx+hy*hy)
        self.assertLess(abs(rho-(mp.mpf(r0)+k*(p[2]-mp.mpf(z0)))), 1e-15)
        self.assertLess(abs(hx/rho*T[0]+hy/rho*T[1]-k*T[2]), 1e-15)


if __name__ == '__main__':
    unittest.main()
