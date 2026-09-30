"""Independent checks of S9d.2c's closed forms
(`generate_spheres_turned_boolean_fixtures.py`) against direct quadrature,
and of the reference's turned caps (`spheres_boolean_reference.py`:
`AxisSphere`, the oblique prism's faces) on the simplest pairs."""
import math
import unittest

import mpmath as mp

import generate_spheres_turned_boolean_fixtures as fixtures
import sphere_boolean_reference as s1
import spheres_boolean_reference as ref
from identity_reference import Boundary, Case
from fractions import Fraction as F

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (0.0, 0.0, 0.0, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)
HP = math.pi/2


def ball(r, frame, low=-HP, high=HP):
    return Case('t', 1e-7, 91, frame, 0.0, 0.0, [], sphere=(float(r), low, high))


def rod(cx, cy, r, z0, z1, frame=XY):
    return Case('t', 1e-7, 92, frame, float(z0), float(z1), [Boundary(circle=(cx, cy, r))])


class SpheresTurnedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_axis_sphere_planes(self):
        # In an exact frame the kernel's reading is S9d.1's affine plane,
        # tuple for tuple; in a turned frame its normal is the stored axis
        # and it passes through o + h n.
        a, b = s1.Sphere(ball(2, XY, 0.0, HP)), ref.AxisSphere(ball(2, XY, 0.0, HP))
        self.assertEqual(a.planes, b.planes)
        t = ref.AxisSphere(ball(2, TILT, -0.5, 0.25))
        _, _, _, n = (tuple(F(v) for v in w) for w in ref.stored_axes(TILT))
        for (m, rhs, _), h in zip(t.planes, t.heights):
            self.assertTrue(ref.is_zero(ref.cross(m, n)))
            self.assertEqual(ref.dot(m, ref.scale(n, h)), rhs)
        self.assertNotEqual(s1.Sphere(ball(2, TILT, -0.5, 0.25)).planes, t.planes)

    def test_segment(self):
        # A disc's part above a chord against the integral over its width.
        for rho, y0 in ((mp.mpf(2), mp.mpf('0.5')), (mp.mpf('1.5'), mp.mpf('-1.25'))):
            A, My = fixtures.segment(rho, y0)
            w = lambda y: 2*mp.sqrt(rho*rho-y*y)
            self.close(A, mp.quad(w, [y0, rho]))
            self.close(My, mp.quad(lambda y: y*w(y), [y0, rho]))
        self.close(fixtures.segment(mp.mpf(1), mp.mpf(-2))[0], mp.pi)
        self.assertEqual(fixtures.segment(mp.mpf(1), mp.mpf(2))[0], 0)

    def test_disc_and_ellipse(self):
        # A circle and a concentric ellipse: half the square of the smaller
        # polar radius, integrated round.
        R, a, b = mp.mpf(2), mp.mpf('1.75'), mp.mpf('2.1875')
        re2 = lambda t: 1/(mp.cos(t)**2/a**2+mp.sin(t)**2/b**2)
        tan2 = (1/a**2-1/R**2)/(1/R**2-1/b**2)
        phi = mp.atan(mp.sqrt(tan2))
        want = 4*mp.quad(lambda t: min(re2(t), R*R)/2, [0, phi, mp.pi/2])
        self.close(fixtures.disc_ellipse(R, a, b), want)

    def test_oblique_wall_matches_the_axial_one(self):
        # With the zone's plane normal to the cylinder the oblique wall's
        # varying bound is constant: both walls agree.
        pair = ref.Pair(ball(2, XY, 0.0, HP), rod(0.25, 0, 1, -3, 3))
        P = pair.P
        e = P.elements[0]
        axial = pair.round_wall(P, pair.D, e)
        oblique = pair.oblique_round_wall(P, pair.D, e)
        for k in ('in', 'out'):
            self.close(axial[2][k], oblique[2][k])
        self.close(axial[3], oblique[3])

    def test_halves_of_a_turned_hemisphere(self):
        # A hemisphere in TILT and its complement against a rod sum to the
        # whole sphere sliced along the rod's axis.
        first = fixtures.g2.Boolean('h_common', 'common', ('sphere', 2, TILT, 0.0, HP),
                                    ('prism', [Boundary(circle=(0.5, 0.25, 0.75))], XY, -3.0, 3.0))
        self.assertLess(fixtures.halves(first), self.eps)


if __name__ == '__main__':
    unittest.main()
