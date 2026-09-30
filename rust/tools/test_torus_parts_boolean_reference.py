"""Independent checks of S9d.4c's closed forms (`generate_torus_parts_boolean_
fixtures.py`) against textbook formulas (Pappus's theorems and the volumes
of revolution of the halves' meridian sections), and of the reference
(`torus_parts_boolean_reference.py`) on the simplest pairs."""
import math
import unittest

import mpmath as mp

import torus_parts_boolean_reference as ref
import generate_torus_parts_boolean_fixtures as fixtures
from brep_reference import cos_rn, sin_rn
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
T = ref.TWO_PI
HP = math.pi/2


def torus(low=0.0, high=T, angle=T, R=2.5, r=1.0, op=91):
    return Case('t', 1e-7, op, XY, 0.0, 0.0, [], torus=(float(R), float(r), low, high, angle))


def pipe(a, lo=-2.0, hi=2.0, op=92):
    return Case('t', 1e-7, op, XY, lo, hi, [Boundary(circle=(0.0, 0.0, float(a)))])


def ball(r, centre=(0.0, 0.0, 0.0), low=-HP, high=HP, op=92):
    return Case('t', 1e-7, op, tuple(float(c) for c in centre)+XY[3:], 0.0, 0.0, [], sphere=(float(r), low, high))


class TorusPartsBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.R, self.r = mp.mpf(2.5), mp.mpf(1)

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=mp.mpf(10)**-28):
        self.assertLess(abs(a-b), eps*max(1, abs(b)))

    def halves(self):
        """The outer and inner halves' volumes: `pi int (R +- sqrt(r^2 -
        w^2))^2 dw` over `|w| < r`, `pi (2 R^2 r +- pi R r^2 + 4 r^3 / 3)`."""
        R, r = self.R, self.r
        return (mp.pi*(2*R*R*r+mp.pi*R*r*r+4*r**3/3), mp.pi*(2*R*R*r-mp.pi*R*r*r+4*r**3/3))

    def test_closed_forms_of_the_parts(self):
        outer, inner = self.halves()
        spec = lambda low, high, angle=T: ('torus', 2.5, 1, XY, low, high, angle)
        self.close(fixtures.RadialTorus(spec(-HP, HP)).measures()[0], outer)
        self.close(fixtures.RadialTorus(spec(HP, 3*HP)).measures()[0], inner)
        # Their difference is the whole torus's 2 pi^2 R r^2.
        self.close(outer-inner, 2*mp.pi**2*self.R*self.r**2)
        # A quarter wedge: its turn times pi r^2 R, the turn the rounded
        # direction's angle.
        turn = mp.atan2(sin_rn(HP), cos_rn(HP))
        V, _, A = fixtures.RadialTorus(spec(0.0, T, HP)).measures()
        self.close(V, turn*mp.pi*self.r**2*self.R)
        self.close(A, turn*2*mp.pi*self.R*self.r+2*mp.pi*self.r**2)
        # The reference's own closed forms (S9d.4b.1's) agree.
        for low, high, want in ((-HP, HP, outer), (HP, 3*HP, inner)):
            self.close(ref.TorusPart(torus(low, high)).closed()[0], want)
        part = ref.TorusPart(torus(-HP, HP)).part
        rest = ref.complement_part(part)
        self.assertEqual(part.sigma, 1)
        self.assertEqual(rest.sigma, -1)
        self.close(ref.TorusPart(torus(-HP, HP), rest).closed()[0], inner)

    def test_outer_half_about_a_pipe(self):
        # A pipe of radius 2 through the outer half (its section there
        # `[0, R + q]`, beyond 2): the common is the cylinder between the
        # half's end planes, `8 pi`, its area `16 pi`.
        pair = ref.Pair(torus(-HP, HP), pipe(2))
        n, vol, area, centre = pair.result('common')
        self.assertEqual(n, 1)
        self.close(vol, 8*mp.pi)
        self.close(area, 16*mp.pi)
        outer, _ = self.halves()
        for k in (0, 1):
            self.close(pair.input_measures('A', k)[0], outer)
        self.assertEqual(pair.solids('cut'), 1)
        self.assertEqual(pair.solids('fuse'), 1)
        # The pipe less the half: its two ends.
        self.assertEqual(pair.swept_count('B-A'), 2)

    def test_hemisphere_in_the_tube(self):
        # A hemisphere of radius 1/2 about the tube's core point: inside the
        # torus, the common the hemisphere (2 pi r^3 / 3, 3 pi r^2), the cut
        # (the hemisphere less the torus) empty.
        pair = ref.Pair(ball(0.5, (2.5, 0, 0), 0.0, HP, op=91), torus(op=92))
        n, vol, area, _ = pair.result('common')
        self.assertEqual(n, 1)
        self.close(vol, 2*mp.pi*mp.mpf(0.5)**3/3)
        self.close(area, 3*mp.pi*mp.mpf(0.5)**2)
        self.assertEqual(pair.result('cut')[0], 0)

    def test_quarter_in_a_ball(self):
        # A ball holding the quarter wedge: the common is the wedge.
        pair = ref.Pair(torus(0.0, T, HP), ball(4))
        turn = mp.atan2(sin_rn(HP), cos_rn(HP))
        n, vol, area, _ = pair.result('common')
        self.assertEqual(n, 1)
        self.close(vol, turn*mp.pi*self.r**2*self.R)
        self.close(area, turn*2*mp.pi*self.R*self.r+2*mp.pi*self.r**2)


if __name__ == '__main__':
    unittest.main()
