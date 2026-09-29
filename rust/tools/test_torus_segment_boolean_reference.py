"""Independent checks of S9d.4b.1's closed forms (`generate_torus_segment_
boolean_fixtures.py`) against textbook formulas (Pappus's theorems on the
meridian section), and of the reference (`torus_segment_boolean_
reference.py`) on the simplest pairs."""
import unittest

import mpmath as mp

import torus_segment_boolean_reference as ref
import generate_torus_segment_boolean_fixtures as fixtures
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
HALF_PI = fixtures.HALF_PI


def box(x0, y0, x1, y1, z0, z1, op=92):
    return Case('t', 1e-7, op, XY, float(z0), float(z1),
                [Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])])


def part(low, high, angle, R=2.5, r=1.5, op=91):
    return Case('t', 1e-7, op, XY, 0.0, 0.0, [], torus=(float(R), float(r), low, high, angle))


class TorusSegmentBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30
        self.R, self.r = mp.mpf(2.5), mp.mpf(1.5)

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_halves(self):
        # The outer half: a barrel, slices pi (R + q)^2; the inner half a
        # spool, pi (R - q)^2 (q^2 = r^2 - w^2); their walls by Pappus (the
        # half circle's length pi r, its centroid 2 r / pi off the tube's
        # centre), their end discs of radius R.
        R, r = self.R, self.r
        for shape, sign in (('outer_half', 1), ('inner_half', -1)):
            V, m, wall, ends = fixtures.Shape(shape, (2.5, 1.5)).measures()
            self.close(V, mp.pi*(2*r*R*R+sign*mp.pi*R*r*r+4*r**3/3))
            self.close(wall, 2*mp.pi*(R+sign*2*r/mp.pi)*mp.pi*r)
            self.close(ends, 2*mp.pi*R*R)
            self.close(m[2], 0)
            K = ref.Part(part(*{1: (-HALF_PI, HALF_PI), -1: (HALF_PI, fixtures.THREE_HALVES_PI)}[sign],
                              ref.TWO_PI))
            self.close(K.closed()[0], V)
            self.close(K.closed()[2], wall+ends)

    def test_wedge(self):
        # Pappus: the tube's disc swept through the turn about the axis; the
        # centroid of the disc's sweep through U at (R + r^2 / 4R) sin(U/2) /
        # (U/2) from the axis on the bisector.
        R, r = self.R, self.r
        K = ref.Part(part(0.0, ref.TWO_PI, HALF_PI))
        U = K.Um
        V, mom, A = K.closed()
        self.close(V, U*mp.pi*r*r*R)
        self.close(A, U*2*mp.pi*R*r+2*mp.pi*r*r)
        dist = (R+r*r/(4*R))*mp.sin(U/2)/(U/2)
        self.close(mom[0]/V, dist*mp.cos(U/2))
        self.close(mom[1]/V, dist*mp.sin(U/2))

    def test_band(self):
        # Between its end heights the band is the disc to R + q, above the
        # upper one the tube's annulus 4 pi R q.
        R, r = self.R, self.r
        sh = fixtures.Shape(fixtures.UB, (2.5, 1.5))
        zl, zh = sh.structure[0][:2]
        Q1 = lambda s: (s*mp.sqrt(r*r-s*s)+r*r*mp.asin(s/r))/2
        want = mp.pi*((R*R+r*r)*(zh-zl)-(zh**3-zl**3)/3)+2*mp.pi*R*(Q1(zh)-Q1(zl))+4*mp.pi*R*(Q1(r)-Q1(zh))
        self.close(sh.measures()[0], want)
        K = ref.Part(part(0.5, 2.25, ref.TWO_PI))
        self.close(K.closed()[0], want)

    def test_reference(self):
        R, r = self.R, self.r
        # Beyond a plane through the axis: half the outer half, its centroid
        # the meridian section's second moment over its first.
        pair = ref.Pair(part(-HALF_PI, HALF_PI, ref.TWO_PI), box(0, -6, 6, 6, -2, 2))
        n, V, A, C = pair.result('common')
        self.assertEqual(n, 1)
        Vh = mp.pi*(2*r*R*R+mp.pi*R*r*r+4*r**3/3)/2
        self.close(V, Vh)
        # The cut face: the meridian section on both half-planes, 2 (2 r R +
        # pi r^2 / 2).
        cut_face = 2*(2*r*R+mp.pi*r*r/2)
        wall = 2*mp.pi*(R+2*r/mp.pi)*mp.pi*r
        self.close(A, wall/2+mp.pi*R*R+cut_face)
        # A slab through the spool's waist: its cut two solids.
        slab = ref.Pair(part(HALF_PI, fixtures.THREE_HALVES_PI, ref.TWO_PI), box(-5, -5, 5, 5, -0.5, 0.5))
        self.assertEqual(slab.solids('cut'), 2)
        self.assertEqual(slab.solids('common'), 1)
        # A quarter turn and a box beyond it: no common, fused two solids.
        far = ref.Pair(part(0.0, ref.TWO_PI, HALF_PI), box(-5, -5, -1.5, 5, -2, 2))
        self.assertEqual(far.result('common')[0], 0)
        self.assertEqual(far.solids('fuse'), 2)


if __name__ == '__main__':
    unittest.main()
