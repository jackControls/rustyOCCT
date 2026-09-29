"""Independent checks of S9d.4b.2's closed forms (`generate_torus_curved_
boolean_fixtures.py`) against textbook formulas (Pappus's theorems on the
torus's meridian section), and of the reference (`torus_curved_boolean_
reference.py`) on the simplest pairs."""
import unittest

import mpmath as mp

import torus_curved_boolean_reference as ref
import generate_torus_curved_boolean_fixtures as fixtures
from identity_reference import Boundary, Case

XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
T = ref.TWO_PI


def torus(R=2.5, r=1.0, frame=XY, op=91):
    return Case('t', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, T, T))


def pipe(a, lo=-2.0, hi=2.0, op=92):
    return Case('t', 1e-7, op, XY, lo, hi, [Boundary(circle=(0.0, 0.0, float(a)))])


def ball(r, centre=(0.0, 0.0, 0.0), op=92):
    return Case('t', 1e-7, op, tuple(float(c) for c in centre)+XY[3:], 0.0, 0.0, [],
                sphere=(float(r), -ref.HALF_PI, ref.HALF_PI))


class TorusCurvedBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40
        self.R, self.r = mp.mpf(2.5), mp.mpf(1)

    def tearDown(self):
        mp.mp.dps = self.dps

    def close(self, a, b, eps=mp.mpf(10)**-28):
        self.assertLess(abs(a-b), eps*max(1, abs(b)))

    def pipe_pappus(self, a):
        """The torus's part inside a coaxial cylinder of radius `a` through
        the hole (`R - r < a < R`): the tube's disc less the segment beyond
        the chord `t = a`, revolved. The segment toward the axis has area
        `r^2 acos(d / r) - d sqrt(r^2 - d^2)` (`d = R - a`) and its centroid
        `2 (r^2 - d^2)^(3/2) / 3 A` from the disc's centre; the wall inside is
        the arc's `2 r acos(d / r)` at its centroid `R - r sin(al) / al`
        (`al = acos(d / r)`), the cylinder's wall inside the chord's `2 sqrt(r^2
        - d^2)` at `a`."""
        R, r = self.R, self.r
        d = R-a
        A = r*r*mp.acos(d/r)-d*mp.sqrt(r*r-d*d)
        cent = R-2*(r*r-d*d)**mp.mpf(1.5)/(3*A)
        al = mp.acos(d/r)
        wall = 2*mp.pi*(2*r*al)*(R-r*mp.sin(al)/al)
        return 2*mp.pi*A*cent, wall+2*mp.pi*a*2*mp.sqrt(r*r-d*d)

    def test_closed_forms_pappus(self):
        forms = {c.pair_name: fixtures.pair_forms(c) for c in fixtures.cases() if c.pair_name in ('pipe_hole',)}
        V, A = self.pipe_pappus(mp.mpf(1.75))
        common = forms['pipe_hole']['common']
        self.close(common[0], V)
        self.close(common[1], A)
        # The torus: Pappus's 2 pi^2 R r^2 and 4 pi^2 R r.
        VK, _, AK = fixtures.measures(('torus', 2.5, 1, XY))
        self.close(VK, 2*mp.pi**2*self.R*self.r**2)
        self.close(AK, 4*mp.pi**2*self.R*self.r)

    def test_reference_coaxial_pipe(self):
        pair = ref.Pair(torus(), pipe(1.75))
        V, A = self.pipe_pappus(mp.mpf(1.75))
        n, vol, area, centre = pair.result('common')
        self.assertEqual(n, 1)
        self.close(vol, V)
        self.close(area, A)
        for k in (0, 1):
            Vk, _, Ak = pair.input_measures('A', k)
            self.close(Vk, 2*mp.pi**2*self.R*self.r**2)
            self.close(Ak, 4*mp.pi**2*self.R*self.r)

    def test_ball_around_the_torus(self):
        # A sphere holding the whole torus: the common is the torus, the cut
        # empty, the fuse the ball (its wall alone).
        pair = ref.Pair(torus(), ball(4))
        n, vol, area, _ = pair.result('common')
        self.assertEqual(n, 1)
        self.close(vol, 2*mp.pi**2*self.R*self.r**2)
        self.close(area, 4*mp.pi**2*self.R*self.r)
        self.assertEqual(pair.result('cut')[0], 0)
        n, vol, area, _ = pair.result('fuse')
        self.assertEqual(n, 1)
        self.close(vol, 4*mp.pi*64/3)
        self.close(area, 4*mp.pi*16)

    def test_solids(self):
        # A ball in the hole: apart, the fuse two solids.
        self.assertEqual(ref.Pair(torus(), ball(1.25)).solids('fuse'), 2)
        # A rod across the torus less the torus: its two ends and its middle
        # in the hole.
        rod = Case('t', 1e-7, 91, (-4.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0), 0.0, 9.0,
                   [Boundary(circle=(0.0, 0.0, 0.5))])
        pair = ref.Pair(rod, torus(op=92))
        self.assertEqual(pair.solids('cut'), 3)
        self.assertEqual(pair.solids('common'), 2)


if __name__ == '__main__':
    unittest.main()
