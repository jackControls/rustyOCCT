"""Independent checks of S9d.1's closed forms (`generate_sphere_boolean_
fixtures.py`) against textbook formulas, and of the reference
(`sphere_boolean_reference.py`) on the simplest pairs."""
import math
import unittest

import mpmath as mp

import generate_sphere_boolean_fixtures as fixtures
import sphere_boolean_reference as ref
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
HP = math.pi/2


def box(x0, y0, x1, y1, z0, z1):
    return Case('t', 1e-7, 92, XY, float(z0), float(z1),
                [Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])])


def ball(r, low=-HP, high=HP):
    return Case('t', 1e-7, 91, XY, 0.0, 0.0, [], sphere=(float(r), low, high))


class SphereBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_cap_of_a_half_space(self):
        # Sphere r = 2 above z = 0.5: a cap of height h = 1.5.
        r, h = mp.mpf(2), mp.mpf(1.5)
        f = fixtures.box_sphere([(-5, 5), (-5, 5), (0.5, 5)], (0, 0, 0), 2)
        self.close(f['volume'], mp.pi*h*h*(3*r-h)/3)
        self.close(f['sphere_in'], 2*mp.pi*r*h)
        self.close(f['box_in'], mp.pi*(r*r-(r-h)**2))
        # The cap's centroid: 3 (2r - h)^2 / (4 (3r - h)) above the centre.
        self.close(f['moments'][2]/f['volume'], 3*(2*r-h)**2/(4*(3*r-h)))

    def test_octant_and_zone(self):
        r = mp.mpf(2)
        f = fixtures.box_sphere([(0, 5), (0, 5), (0, 5)], (0, 0, 0), 2)
        self.close(f['volume'], r**3*mp.pi/6)
        for m in f['moments']:
            self.close(m/f['volume'], 3*r/8)
        self.close(f['sphere_in'], mp.pi*r*r/2)
        self.close(f['box_in'], 3*mp.pi*r*r/4)
        # A zone between z = -1 and z = 1.5 inside a large box.
        z = fixtures.box_sphere([(-9, 9), (-9, 9), (-9, 9)], (0, 0, 0), 2, (-1, 1.5))
        V, Mz, A = fixtures.sphere_measures(2, (-1, 1.5))
        self.close(z['volume'], V)
        self.close(z['moments'][2], Mz)
        self.close(z['sphere_in']+z['ends_in'], A)
        self.close(V, mp.pi*(4*mp.mpf(2.5)-(mp.mpf(1.5)**3+1)/3))
        self.close(A, 2*mp.pi*2*mp.mpf(2.5)+mp.pi*(4-1)+mp.pi*(4-mp.mpf(2.25)))

    def test_inside(self):
        f = fixtures.box_sphere([(-3, 3), (-3, 3), (-3, 3)], (0.5, 0.25, 0), 1.5)
        self.close(f['volume'], 4*mp.pi*mp.mpf(1.5)**3/3)
        self.close(f['sphere_in'], 4*mp.pi*mp.mpf(1.5)**2)
        self.assertEqual(f['box_in'], 0)
        g = fixtures.box_sphere([(-1, 1), (-1, 1), (-0.5, 1)], (0, 0, 0), 2)
        self.close(g['volume'], 6)
        self.close(g['box_in'], 20)
        self.assertEqual(g['sphere_in'], 0)

    def test_rectangle_in_disc(self):
        rho, (x0, x1, y0, y1) = mp.mpf(2.5), (mp.mpf(1), mp.mpf(6), mp.mpf(-0.5), mp.mpf(6))
        q = lambda x: mp.sqrt(max(rho*rho-x*x, 0))
        top = lambda x: min(y1, q(x))
        bot = lambda x: max(y0, -q(x))
        ends = [x0, mp.sqrt(rho*rho-y0*y0), rho]
        A = mp.quad(lambda x: max(top(x)-bot(x), 0), ends)
        X = mp.quad(lambda x: x*max(top(x)-bot(x), 0), ends)
        Y = mp.quad(lambda x: (top(x)**2-bot(x)**2)/2, ends)
        a, mx, my = fixtures.rect_disc(x0, x1, y0, y1, rho)
        for p, w in ((a, A), (mx, X), (my, Y)):
            self.close(p, w, mp.mpf(10)**-25)
        self.close(fixtures.arc_in_rect(-3, 3, 0, 3, 2), mp.pi)

    def test_reference_cap_and_zone(self):
        # The sphere through a box's face: the common a cap of height 1.
        pair = ref.Pair(ball(2), box(-3, -3, 3, 3, 1, 4))
        n, V, A, C = pair.result('common')
        self.assertEqual(n, 1)
        self.close(V, 5*mp.pi/3)
        self.close(A, 7*mp.pi)
        self.close(C[2], mp.mpf('1.35'))
        n, V, A, C = pair.result('cut')
        self.close(V, 9*mp.pi)
        self.close(A, 15*mp.pi)
        # A zone's closed form from its stored heights.
        S = ref.Sphere(ball(2, -0.5, 1.0))
        lo, hi = ref.M(ref.F(2*math.sin(-0.5))), ref.M(ref.F(2*math.sin(1.0)))
        V, _, A = S.closed()
        self.close(V, mp.pi*((4*hi-hi**3/3)-(4*lo-lo**3/3)))
        self.close(A, 4*mp.pi*(hi-lo)+mp.pi*(8-lo*lo-hi*hi))

    def test_solids(self):
        # A thin slab cuts the sphere in two; a bar through it is cut in two.
        self.assertEqual(ref.Pair(ball(2), box(-3, -3, 3, 3, -0.25, 0.25)).solids('cut'), 2)
        bar = box(-0.5, -0.5, 0.5, 0.5, -4, 4)
        bar.operation = 91
        self.assertEqual(ref.Pair(bar, ball(2)).solids('cut'), 2)
        # A face tangent: no common, two solids fused.
        tangent = ref.Pair(ball(2), box(-3, -3, 3, 3, 2, 4))
        self.assertEqual(tangent.solids('common'), 0)
        self.assertEqual(tangent.solids('fuse'), 2)

    def test_distance_to_a_polyhedron(self):
        F = ref.F
        # The box [1, 2]^3 from the origin: its corner (1, 1, 1).
        ineqs = [((F(1), F(0), F(0)), F(1)), ((F(-1), F(0), F(0)), F(-2)),
                 ((F(0), F(1), F(0)), F(1)), ((F(0), F(-1), F(0)), F(-2)),
                 ((F(0), F(0), F(1)), F(1)), ((F(0), F(0), F(-1)), F(-2))]
        self.assertEqual(ref.dist2((F(0), F(0), F(0)), ineqs), 3)
        self.assertEqual(ref.dist2((F(3, 2), F(0), F(3, 2)), ineqs), 1)
        self.assertIsNone(ref.dist2((F(0),)*3, ineqs+[((F(1), F(0), F(0)), F(3))]))


if __name__ == '__main__':
    unittest.main()
