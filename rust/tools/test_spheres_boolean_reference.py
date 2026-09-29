"""Independent checks of S9d.2's closed forms (`generate_spheres_boolean_
fixtures.py`) against textbook formulas, and of the reference
(`spheres_boolean_reference.py`) on the simplest pairs."""
import math
import unittest

import mpmath as mp

import generate_spheres_boolean_fixtures as fixtures
import spheres_boolean_reference as ref
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
HP = math.pi/2


def ball(r, o=(0, 0, 0), low=-HP, high=HP):
    return Case('t', 1e-7, 91, tuple(float(c) for c in o)+XY[3:], 0.0, 0.0, [], sphere=(float(r), low, high))


def rod(cx, cy, r, z0, z1):
    return Case('t', 1e-7, 92, XY, float(z0), float(z1), [Boundary(circle=(cx, cy, r))])


class SpheresBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_spherical_ring(self):
        # A sphere less a coaxial cylinder through it: pi h^3 / 6 whatever
        # the radius, h the ring's height.
        for R, a in ((2.5, 1.5), (5, 4), (2, 1)):
            h = 2*mp.sqrt(mp.mpf(R)**2-mp.mpf(a)**2)
            f = fixtures.coax(R, a, -9, 9)
            self.close(4*mp.pi*mp.mpf(R)**3/3-f['volume'], mp.pi*h**3/6)
            # Its surface: the zone outside the cylinder and the hole's wall.
            self.close(4*mp.pi*mp.mpf(R)**2-f['sphere_in'], 2*mp.pi*R*h)
            self.close(f['wall_in'], 2*mp.pi*a*h)

    def test_cap_inside_a_cylinder(self):
        # The common of a sphere and a coaxial cylinder: the cylinder
        # between the parallels and two caps of height R - h/2.
        R, a = mp.mpf(2.5), mp.mpf(1.5)
        h = 2*mp.sqrt(R*R-a*a)
        f = fixtures.coax(R, a, -9, 9)
        cap = R-h/2
        self.close(f['volume'], mp.pi*a*a*h+2*mp.pi*cap*cap*(3*R-cap)/3)
        self.close(f['sphere_in'], 2*2*mp.pi*R*cap)
        # A cylinder ending at the equator: one cap and half the band.
        g = fixtures.coax(R, a, 0, 9)
        self.close(g['volume'], mp.pi*a*a*h/2+mp.pi*cap*cap*(3*R-cap)/3)
        self.close(g['caps_in'], mp.pi*a*a)
        # A hemisphere (a cap to the equator) against the same cylinder:
        # its disc's part inside the cylinder.
        k = fixtures.coax(R, a, -9, 9, (0, R))
        self.close(k['volume'], g['volume'])
        self.close(k['ends_in'], mp.pi*a*a)

    def test_lens_of_two_spheres(self):
        R1, R2, d = mp.mpf(2), mp.mpf(1.5), mp.mpf(2.25)
        forms = fixtures.spheres_forms(('sphere', 2, (0, 0, 0)+XY[3:], -HP, HP),
                                       ('sphere', 1.5, (2.25, 0, 0)+XY[3:], -HP, HP))
        V, A, C = forms['common']
        # The lens: pi (R1 + R2 - d)^2 (d^2 + 2 d (R1 + R2) - 3 (R1 - R2)^2) / (12 d).
        self.close(V, mp.pi*(R1+R2-d)**2*(d*d+2*d*(R1+R2)-3*(R1-R2)**2)/(12*d))
        x = (d*d+R1*R1-R2*R2)/(2*d)
        self.close(A, 2*mp.pi*R1*(R1-x)+2*mp.pi*R2*(R2-(d-x)))
        # The lens's centroid lies on the line of centres.
        self.close(C[1], 0)
        self.close(C[2], 0)
        # Apart and nested.
        apart = fixtures.spheres_forms(('sphere', 1, (0, 0, 0)+XY[3:], -HP, HP),
                                       ('sphere', 1, (3, 0, 0)+XY[3:], -HP, HP))
        self.assertIsNone(apart['common'])
        nested = fixtures.spheres_forms(('sphere', 2, (0, 0, 0)+XY[3:], -HP, HP),
                                        ('sphere', 0.5, (1, 0, 0)+XY[3:], -HP, HP))
        self.close(nested['cut'][0], 4*mp.pi*(8-mp.mpf(0.125))/3)
        self.close(nested['cut'][1], 4*mp.pi*(4+mp.mpf(0.25)))

    def test_lens_of_two_discs(self):
        # Equal discs a radius apart: 2 (pi/3 - sqrt(3)/4) r^2, centroid
        # halfway, each circle's angle inside the other 2 pi/3.
        r = mp.mpf(1)
        A, m, t1, t2 = fixtures.lens(r, r, r)
        self.close(A, 2*(mp.pi/3-mp.sqrt(3)/4))
        self.close(m/A, mp.mpf(1)/2)
        self.close(t1, 2*mp.pi/3)
        self.close(t2, 2*mp.pi/3)
        # Nested: the small disc whole.
        A, m, t1, t2 = fixtures.lens(mp.mpf(3), r, mp.mpf(1))
        self.close(A, mp.pi)
        self.close(m, mp.pi)
        self.close(t2, 2*mp.pi)

    def test_reference_ring(self):
        # The reference's slicing on the ring (its sections concentric, the
        # breakpoints where they coincide).
        pair = ref.Pair(ball(2.5), rod(0, 0, 1.5, -3, 3))
        res = pair.sliced()
        self.close(res['D-P'][0], mp.pi*64/6)
        self.close(res[('sphere', 'A')]['in'], 2*2*mp.pi*mp.mpf(2.5)*mp.mpf(0.5))
        self.assertEqual(pair.components(), {'common': 1, 'fuse': 1, 'D-P': 1, 'P-D': 2})

    def test_green_of_an_ellipse(self):
        # An ellipse's arc: area and moments of its whole region.
        C, a, b = (mp.mpf(0.5), mp.mpf(-1)), mp.mpf(2), mp.mpf(0.75)
        th = mp.mpf(0.3)
        A = (a*mp.cos(th), -b*mp.sin(th), a*mp.sin(th), b*mp.cos(th))
        g = ref.green_arc(C, A, mp.mpf(0), 2*mp.pi)
        self.close(g[0], mp.pi*a*b)
        self.close(g[1], mp.pi*a*b*C[0])
        self.close(g[2], mp.pi*a*b*C[1])

    def test_tangency_polynomial(self):
        # Two circles in a plane, the second's centre moving: tangent where
        # the distance is the radii's sum or difference.
        from fractions import Fraction as F
        c1 = ref.Conic([[F(1), F(0)], [F(0), F(1)]], ([F(0)], [F(0)]), [F(-4)])
        # |q - (s, 0)|^2 - 1: g = (-s, 0), k = s^2 - 1.
        c2 = ref.Conic([[F(1), F(0)], [F(0), F(1)]], ([F(0), F(-1)], [F(0)]), [F(-1), F(0), F(1)])
        roots = sorted(ref.real_roots(ref.tangency_poly(c1, c2)))
        for want in (-3, -1, 1, 3):
            self.assertTrue(any(abs(r-want) < self.eps for r in roots), (want, roots))


if __name__ == '__main__':
    unittest.main()
