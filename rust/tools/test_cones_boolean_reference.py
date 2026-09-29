"""Independent checks of S9d.3b's closed forms (`generate_cones_boolean_
fixtures.py`) against textbook formulas, of the reference's exact algebra
(`cones_boolean_reference.py`: resultants, roots, the rulings' polynomials)
and of the reference on the simplest pairs."""
import math
import unittest
from fractions import Fraction as F

import mpmath as mp

import cones_boolean_reference as ref
import generate_cones_boolean_fixtures as fixtures
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)


def cone(r0, r1, h, frame=XY, op=91):
    return Case('t', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))


def pipe(cx, cy, r, z0, z1, frame=XY, op=92):
    return Case('t', 1e-7, op, frame, float(z0), float(z1), [Boundary(circle=(float(cx), float(cy), float(r)))])


def ball(r, frame, op=92):
    return Case('t', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), -math.pi/2, math.pi/2))


class ConesBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_coaxial_cylinder(self):
        # The frustum (radii 2 and 1, height 2) and the coaxial cylinder of
        # radius 3/2: below w = 1 the cylinder's disc, above the frustum's.
        B = fixtures.axial_of(('prism', [Boundary(circle=(0.0, 0.0, 1.5))], XY, -1.0, 3.0),
                              (mp.mpf(0),)*3, fixtures.ideal_axes('XY'))
        f = fixtures.axial_form(2, 1, 2, B)
        upper = mp.pi*(mp.mpf(2.25)+mp.mpf(1.5)+1)/3
        self.close(f['volume'], mp.pi*mp.mpf(2.25)+upper)
        # The frustum's wall above w = 1, its top disc and its base's
        # middle; the cylinder's wall below w = 1 and none of its caps.
        self.close(f['K_in'], mp.pi*mp.mpf(2.5)*mp.sqrt(mp.mpf(1.25))+mp.pi+mp.pi*mp.mpf(2.25))
        self.close(f['B_in'], 2*mp.pi*mp.mpf(1.5))

    def test_coaxial_sphere(self):
        # The cone of radius 3 and height 4 and the sphere of radius 3/4 at
        # height 11/4 inscribed in it: the common is the ball, its face
        # inside by Archimedes, the cone's wall inside nothing.
        B = fixtures.axial_of(('sphere', 0.75, (0.0, 0.0, 2.75, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0), -math.pi/2,
                               math.pi/2), (mp.mpf(0),)*3, fixtures.ideal_axes('XY'))
        f = fixtures.axial_form(3, 0, 4, B)
        self.close(f['volume'], 4*mp.pi*mp.mpf(0.75)**3/3, mp.mpf(10)**-25)
        self.close(f['B_in'], 4*mp.pi*mp.mpf(0.75)**2, mp.mpf(10)**-25)
        self.assertLess(f['K_in'], mp.mpf(10)**-20)

    def test_lens_of_cones(self):
        # Two equal cones of parallel axes 1 apart: at each height the lens
        # of two equal discs, symmetric about the midpoint.
        B = fixtures.axial_of(('cone', 1.5, 0.0, 3.0, (1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)),
                              (mp.mpf(0),)*3, fixtures.ideal_axes('XY'))
        f = fixtures.axial_form(1.5, 0, 3, B)

        def lens_area(r):
            if r <= mp.mpf(0.5):
                return mp.mpf(0)
            return 2*(r*r*mp.acos(mp.mpf(0.5)/r)-mp.mpf(0.5)*mp.sqrt(r*r-mp.mpf(0.25)))
        want = mp.quad(lambda z: lens_area(mp.mpf(1.5)-z/2), [0, 2, 3])
        self.close(f['volume'], want)
        self.close(f['moments'][0]/f['volume'], mp.mpf(0.5))
        self.close(f['K_in'], f['B_in'])

    def test_strip(self):
        # A rod of radius 1/2 across the frustum through its axis at height
        # 1: the volume by a direct double quadrature of the disc's chords.
        f = fixtures.strip_form(2, 1, 2, 0, 1, 0.5)

        def area(z):
            q = mp.sqrt(mp.mpf(0.25)-(z-1)**2)
            r = 2-z/2
            return mp.quad(lambda y: 2*mp.sqrt(r*r-y*y), [-q, q])
        with mp.workdps(30):
            want = mp.quad(area, [mp.mpf(0.5), 1, mp.mpf(1.5)])
        self.close(f['volume'], want, mp.mpf(10)**-25)
        # The rod's wall inside: its circumference times the chord at
        # each height, integrated over the rod's angle.
        wall = mp.quad(lambda t: mp.mpf(0.5)*2*mp.sqrt((2-(1+mp.sin(t)/2)/2)**2-(mp.cos(t)/2)**2), [0, 2*mp.pi])
        self.close(f['B_in'], wall, mp.mpf(10)**-25)

    def test_resultant_and_roots(self):
        # (w - t)(w - 2) and (w - 1)(w + t): a common root exactly where t
        # = 1 (w = 1), t = 0 (w = 0) or t = -2 (w = 2).
        P = ref.Pol
        a = (P([0, 2]), P([-2, -1]), P([1]))
        b = (P([0, -1]), P([-1, 1]), P([1]))
        r = ref.resultant(a, b)
        roots = ref.numeric_roots(r.poly())
        for want in (-2, 0, 1):
            self.assertTrue(any(abs(x-want) < self.eps for x in roots))
        self.assertTrue(all(min(abs(x-w) for w in (-2, 0, 1)) < mp.mpf(10)**-20 for x in roots))
        # A trigonometric polynomial's roots through tan(theta / 2), the
        # half turn included.
        t = ref.Trig({(1, 0): F(1), (0, 0): F(1)})
        self.assertEqual([mp.nstr(x, 10) for x in ref.ring_roots(t, mp.mpf(0), 2*mp.pi, 'trig')],
                         [mp.nstr(mp.pi, 10)])

    def test_ruling_family(self):
        # The frustum's rulings against the coaxial cylinder of radius 3/2:
        # every ruling meets it at w = 1.
        pair = ref.Pair(cone(2, 1, 2), pipe(0, 0, 1.5, -1, 3))
        wall = [f for f in pair.A.faces() if f.tag == ('wall',)][0]
        c0, c1, c2 = pair.B.surfaces[0].family(wall.A, wall.B, 'trig')
        for th in (mp.mpf(0), mp.mpf(1), mp.mpf(4)):
            roots = ref.solve_quadratic(c0.value(th), c1.value(th), c2.value(th))
            self.assertTrue(any(abs(w-1) < mp.mpf(10)**-35 for w in roots))

    def test_reference(self):
        pair = ref.Pair(cone(2, 1, 2), pipe(0, 0, 1.5, -1, 3))
        n, V, A, C = pair.result('common')
        self.assertEqual(n, 1)
        upper = mp.pi*(mp.mpf(2.25)+mp.mpf(1.5)+1)/3
        self.close(V, mp.pi*mp.mpf(2.25)+upper)
        self.close(A, mp.pi*mp.mpf(2.5)*mp.sqrt(mp.mpf(1.25))+mp.pi+2*mp.pi*mp.mpf(1.5)+mp.pi*mp.mpf(2.25))
        # The pipe less the frustum: its two ends.
        self.assertEqual(ref.Pair(pipe(0, 0, 1.5, -1, 3, op=91), cone(2, 1, 2, op=92)).solids('cut'), 2)
        # The frustum and a sphere on its axis: the wall inside two ways.
        pair = ref.Pair(cone(2, 1, 2), ball(1.5, (0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)))
        r = pair.sliced()
        wall = [c for w, tag, c, _ in pair.faces() if w == 'A' and tag == ('wall',)][0]
        self.close(r[('curved', 'A')]['in'], wall['in'])


if __name__ == '__main__':
    unittest.main()
