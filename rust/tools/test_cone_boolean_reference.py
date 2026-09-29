"""Independent checks of S9d.3a's closed forms (`generate_cone_boolean_
fixtures.py`) against textbook formulas, and of the reference
(`cone_boolean_reference.py`) on the simplest pairs."""
import unittest

import mpmath as mp

import cone_boolean_reference as ref
import generate_cone_boolean_fixtures as fixtures
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)


def box(x0, y0, x1, y1, z0, z1, op=92):
    return Case('t', 1e-7, op, XY, float(z0), float(z1),
                [Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])])


def cone(r0, r1, h, op=91):
    return Case('t', 1e-7, op, XY, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))


class ConeBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_cone_and_frustum(self):
        V, Mw, wall, total = fixtures.cone_measures(1.5, 0, 3)
        self.close(V, mp.pi*mp.mpf(2.25))
        self.close(Mw/V, mp.mpf(3)/4)
        self.close(wall, mp.pi*mp.mpf(1.5)*mp.sqrt(9+mp.mpf(2.25)))
        V, Mw, wall, total = fixtures.cone_measures(2, 1, 2)
        self.close(V, 14*mp.pi/3)
        # A frustum's centroid: h (R^2 + 2 R r + 3 r^2) / (4 (R^2 + R r + r^2)).
        self.close(Mw/V, mp.mpf(2)*(4+4+3)/(4*7))
        self.close(total, 3*mp.pi*mp.sqrt(5)+5*mp.pi)

    def test_normal_half_space(self):
        # The frustum above z = 1: radii 1.5 and 1 over a height of 1.
        f = fixtures.halfspace_form(2, 1, 2, (mp.mpf(0), mp.mpf(0), mp.mpf(1)), 1)
        self.close(f['volume'], mp.pi*(mp.mpf(2.25)+mp.mpf(1.5)+1)/3)
        self.close(f['face_in'], mp.pi*mp.mpf(2.25))
        self.close(f['cone_in'], mp.pi*mp.mpf(2.5)*mp.sqrt(1+mp.mpf(0.25))+mp.pi)

    def test_oblique_cone(self):
        # The cone of radius 1.5 and apex at height 3 above the plane w = 1 +
        # u/2: an oblique cone on the section's ellipse, a third of its area
        # times the apex's distance, its centroid three quarters of the way
        # to the ellipse's centre, the wall inside the projected ellipse
        # over sin(half angle).
        m, c = mp.mpf(1)/2, mp.mpf(1)
        u1, u2 = (3-c)/(2+m), -(3-c)/(2-m)
        a = abs(u1-u2)*mp.sqrt(1+m*m)/2
        uc = (u1+u2)/2
        wc = c+m*uc
        b = mp.sqrt(((3-wc)/2)**2-uc*uc)
        A = mp.pi*a*b
        f = fixtures.halfspace_form(1.5, 0, 3, (-m, mp.mpf(0), mp.mpf(1)), c)
        V = A*(3-c)/mp.sqrt(1+m*m)/3
        self.close(f['volume'], V)
        self.close(f['face_in'], A)
        self.close(f['cone_in'], A/mp.sqrt(1+m*m)*mp.sqrt(5))
        centre = (3*uc/4, mp.mpf(0), 3+3*(wc-3)/4)
        for i in range(3):
            self.close(f['moments'][i]/V, centre[i])

    def test_hyperbolic_segment(self):
        # The frustum (radii 2 and 1, height 2) beyond x = d = 5/4: the
        # integral over r of the circular segment, by its antiderivative
        # G(r) = r^3 acos(d/r)/3 - 2 d r q / 3 + d^3 ln(r + q) / 3, q^2 = r^2 - d^2.
        d = mp.mpf(1.25)

        def G(r):
            q = mp.sqrt(r*r-d*d)
            return r**3*mp.acos(d/r)/3-2*d*r*q/3+d**3*mp.log(r+q)/3
        want = 2*(G(mp.mpf(2))-G(d))
        f = fixtures.halfspace_form(2, 1, 2, (mp.mpf(1), mp.mpf(0), mp.mpf(0)), d)
        self.close(f['volume'], want)
        g = fixtures.box_form(2, 1, 2, [(1.25, 9), (-9, 9), (-9, 9)])
        self.close(g['volume'], want)
        # The reference's row for the same common (the `hyperbola` fixture).
        pair = ref.Pair(box(1.25, -4, 4, 4, -1, 3), cone(2, 1, 2, 92))
        n, V, _, _ = pair.result('common')
        self.assertEqual(n, 1)
        self.close(V, want)

    def test_segment(self):
        r, d = mp.mpf(2), mp.mpf('0.5')
        area, mom, angle, chord = fixtures.segment(r, d)
        self.close(area, mp.quad(lambda x: 2*mp.sqrt(r*r-x*x), [d, r]))
        self.close(mom, mp.quad(lambda x: 2*x*mp.sqrt(r*r-x*x), [d, r]))
        self.close(angle, 2*mp.acos(d/r))
        self.close(chord, 2*mp.sqrt(r*r-d*d))

    def test_half_cone(self):
        # A wall through the axis halves the cone: its centroid R / pi
        # (R = 1.5) off the axis, a quarter of the height up.
        f = fixtures.box_form(1.5, 0, 3, [(0, 5), (-5, 5), (-1, 4)])
        self.close(f['volume'], mp.pi*mp.mpf(2.25)/2)
        self.close(f['moments'][0]/f['volume'], mp.mpf(1.5)/mp.pi)
        self.close(f['moments'][2]/f['volume'], mp.mpf(3)/4)
        self.close(f['box_in'], 3*3/2)

    def test_reference(self):
        pair = ref.Pair(cone(2, 1, 2), box(-3, -3, 3, 3, 1, 5))
        n, V, A, C = pair.result('common')
        self.assertEqual(n, 1)
        self.close(V, mp.pi*(mp.mpf(2.25)+mp.mpf(1.5)+1)/3)
        self.close(A, mp.pi*(mp.mpf(2.25)+1+mp.mpf(2.5)*mp.sqrt(mp.mpf(1.25))))
        self.close(C[2], 1+(mp.mpf(2.25)+2*mp.mpf(1.5)+3)/(4*(mp.mpf(2.25)+mp.mpf(1.5)+1)))
        # A slab through the axis cuts the frustum in two.
        self.assertEqual(ref.Pair(cone(2, 1, 2), box(-0.25, -3, 0.25, 3, -1, 3)).solids('cut'), 2)
        # A box on the top disc: fused one solid, no common.
        top = ref.Pair(cone(2, 1, 2), box(-3, -3, 3, 3, 2, 3))
        self.assertEqual(top.result('common')[0], 0)
        self.assertEqual(top.solids('fuse'), 1)


if __name__ == '__main__':
    unittest.main()
