"""Independent checks of S9d.4a's closed forms (`generate_torus_boolean_
fixtures.py`) against textbook formulas (Pappus's theorems), and of the
reference (`torus_boolean_reference.py`) on the simplest pairs."""
import unittest

import mpmath as mp

import torus_boolean_reference as ref
import generate_torus_boolean_fixtures as fixtures
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)


def box(x0, y0, x1, y1, z0, z1, op=92):
    return Case('t', 1e-7, op, XY, float(z0), float(z1),
                [Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])])


def torus(R, r, op=91):
    return Case('t', 1e-7, op, XY, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, ref.TWO_PI, ref.TWO_PI))


class TorusBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.eps = mp.mpf(10)**-30
        self.R, self.r = mp.mpf(2.5), mp.mpf(1.5)

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=None):
        self.assertLess(abs(a-b), (eps or self.eps)*max(1, abs(b)))

    def test_torus(self):
        # Pappus: the disc of radius r about a circle of radius R.
        V, A = fixtures.torus_measures(2.5, 1.5)
        self.close(V, 2*mp.pi*self.R*mp.pi*self.r**2)
        self.close(A, 2*mp.pi*self.R*2*mp.pi*self.r)

    def test_upper_half(self):
        # Above the equatorial plane: half the torus, its centroid the half
        # disc's (4 r / 3 pi) above it, half the wall, the annulus's face.
        R, r = self.R, self.r
        f = fixtures.halfspace_form(R, r, (0, 0, 1), 0)
        self.close(f['volume'], mp.pi**2*R*r*r)
        self.close(f['moments'][2]/f['volume'], 4*r/(3*mp.pi))
        self.close(f['torus_in'], 2*mp.pi**2*R*r)
        self.close(f['face_in'], mp.pi*((R+r)**2-(R-r)**2))

    def test_half_through_the_axis(self):
        # Beyond a plane through the axis: half the torus, its centroid by
        # Pappus's second moment (4 R^2 + r^2) / (2 pi R) off the axis, the
        # two discs of the tube on the plane.
        R, r = self.R, self.r
        f = fixtures.halfspace_form(R, r, (1, 0, 0), 0)
        self.close(f['volume'], mp.pi**2*R*r*r)
        self.close(f['moments'][0]/f['volume'], (4*R*R+r*r)/(2*mp.pi*R))
        self.close(f['torus_in'], 2*mp.pi**2*R*r)
        self.close(f['face_in'], 2*mp.pi*r*r)
        # The same by an aligned box.
        g = fixtures.box_form(R, r, [(0, 9), (-9, 9), (-9, 9)])
        self.close(g['volume'], f['volume'])
        self.close(g['moments'][0], f['moments'][0])
        self.close(g['torus_in'], f['torus_in'])
        self.close(g['box_in'], f['face_in'])

    def test_band(self):
        # |w| <= h: slices of area 4 pi R q, q^2 = r^2 - w^2; the wall's
        # latitudes |phi| <= asin(h / r) on both sides of the tube.
        R, r, h = self.R, self.r, mp.mpf(0.5)
        q = mp.sqrt(r*r-h*h)
        g = fixtures.box_form(R, r, [(-9, 9), (-9, 9), (-h, h)])
        self.close(g['volume'], 4*mp.pi*R*(h*q+r*r*mp.asin(h/r)))
        self.close(g['torus_in'], 8*mp.pi*R*r*mp.asin(h/r))
        self.close(g['box_in'], 2*4*mp.pi*R*q)

    def test_reference(self):
        R, r = self.R, self.r
        pair = ref.Pair(torus(2.5, 1.5), box(0, -5, 5, 5, -2, 2))
        n, V, A, C = pair.result('common')
        self.assertEqual(n, 1)
        self.close(V, mp.pi**2*R*r*r)
        self.close(A, 2*mp.pi**2*R*r+2*mp.pi*r*r)
        self.close(C[0], (4*R*R+r*r)/(2*mp.pi*R))
        # A slab across the tube: the torus less it two solids, the slab
        # less the torus two (the hole's disc and the outside).
        band = ref.Pair(torus(2.5, 1.5), box(-5, -5, 5, 5, -0.5, 0.5))
        self.assertEqual(band.solids('cut'), 2)
        self.assertEqual(band.solids('common'), 1)
        slab = ref.Pair(box(-5, -5, 5, 5, -0.5, 0.5), torus(2.5, 1.5, 92))
        self.assertEqual(slab.solids('cut'), 2)
        # A box in the hole: no common, fused two solids.
        hole = ref.Pair(torus(2.5, 1.5), box(-0.625, -0.625, 0.625, 0.625, -2, 2))
        self.assertEqual(hole.result('common')[0], 0)
        self.assertEqual(hole.solids('fuse'), 2)


if __name__ == '__main__':
    unittest.main()
