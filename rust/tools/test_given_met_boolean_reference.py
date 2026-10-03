"""Independent checks of S9e.3b's triple points (`given_met_reference.py`) on
chains whose meetings are known in closed form: a sphere fused with a peg,
its `Rise` loop met by a wall through the peg's axis at `(+-1, 2, 2 sqrt
5)`; two crossing rods, their rings met by a wall through the thin rod's axis
at `(0, 3/2, +-sqrt 7 / 2)` and `(0, -1/2, +-sqrt 15 / 2)`; a frustum's
elliptic section by a tilted plane met by a coaxial-free rod (each point on
all three surfaces and the normals' determinant against a direct
computation); a declared tangency (a plane touching the loop at its lowest
point, determinant zero); a meeting of the two surfaces off the given
result's edges (where the peg is absent: no edge);
and the fixture list's classes (`generate_given_met_boolean_fixtures.py`)."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
from sphere_boolean_reference import cross
import generate_given_curved_boolean_fixtures as gc
import generate_given_met_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, disc, square
import given_met_reference as gm


def chain(specs):
    return ref.Chain([gc.solid(s, gc.OPERATIONS[k]) for k, s in enumerate(specs)])


def box(x0, y0, z0, x1, y1, z1):
    return fixtures.big_box(x0, y0, z0, x1, y1, z1)


class GivenMetReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def assert_points(self, got, want):
        self.assertEqual(len(got), len(want), [g[0] for g in got])
        for w in want:
            d = min(max(abs(g[0][i]-w[i]) for i in range(3)) for g in got)
            self.assertLess(d, mp.mpf(10)**-30, w)

    def test_rise_met_by_a_wall_through_the_pegs_axis(self):
        c = chain([fixtures.BALL5, fixtures.PEG, box(-8.0, 2.0, -8.0, 8.0, 8.0, 10.0)])
        pts = gm.Tracer(c, ('fuse', 0, 1), 2).points()
        s = 2*mp.sqrt(5)
        self.assert_points(pts, [(1, 2, s), (-1, 2, s)])
        for X, sine, kinds, found in pts:
            # Normals (x, y, z), (x, y - 2, 0), (0, 1, 0) at (+-1, 2, 2 sqrt 5).
            n1, n2, n3 = X, (X[0], X[1]-2, 0), (0, 1, 0)
            d = abs(gm.det3(gm.unit(n1), gm.unit(n2), gm.unit(n3)))
            self.assertLess(abs(sine-d), mp.mpf(10)**-30)
            self.assertEqual(sorted(kinds[:2]), ['cylinder', 'sphere'])
            self.assertEqual(kinds[2], 'plane')
            self.assertGreaterEqual(found, 2)

    def test_rings_met_by_a_wall_through_the_thin_rods_axis(self):
        c = chain([fixtures.ROD_X, fixtures.ROD_Z, box(0.0, -8.0, -8.0, 8.0, 8.0, 8.0)])
        pts = gm.Tracer(c, ('fuse', 0, 1), 2).points()
        a, b = mp.sqrt(7)/2, mp.sqrt(15)/2
        self.assert_points(pts, [(0, 1.5, a), (0, 1.5, -a), (0, -0.5, b), (0, -0.5, -b)])

    def test_points_lie_on_all_three_surfaces(self):
        c = chain([fixtures.FRUSTUM, fixtures.FRUSTUM_CUT,
                   gc.prism([disc(0.1, 2.5, 0.4)], at('XY', (0, 0, 0)), -1.0, 6.0)])
        pts = gm.Tracer(c, ('cut', 0, 1), 2).points()
        self.assertEqual(len(pts), 2)
        # The slab's plane from its stored axes (its chart's third row: the
        # stored frame's `x` and `y` crossed over its determinant) and the
        # rod's binary64 centre and radius.
        o, ax, ay, an = ([F(v) for v in w] for w in stored_axes(at('TILT', (0, 0, 0))))
        det = sum(ax[i]*cross(ay, an)[i] for i in range(3))
        n = [mp.mpf(v.numerator)/v.denominator for v in (c/det for c in cross(ax, ay))]
        xc, yc, rr = mp.mpf(0.1), mp.mpf(2.5), mp.mpf(0.4)
        for X, sine, kinds, _ in pts:
            # The cone x^2 + y^2 = (3 - z / 2)^2, the plane n . X = 2.2, the
            # rod (x - xc)^2 + (y - yc)^2 = rr^2.
            x, y, z = X
            self.assertLess(abs(x*x+y*y-(3-z/2)**2), mp.mpf(10)**-30)
            self.assertLess(abs(n[0]*x+n[1]*y+n[2]*z-mp.mpf(2.2)), mp.mpf(10)**-30)
            self.assertLess(abs((x-xc)**2+(y-yc)**2-rr*rr), mp.mpf(10)**-30)
            normals = [(x, y, (3-z/2)/2), tuple(n), (x-xc, y-yc, 0)]
            d = abs(gm.det3(*(gm.unit(v) for v in normals)))
            self.assertLess(abs(sine-d), mp.mpf(10)**-30)
            self.assertTrue(gm.in_scope(*kinds))

    def test_a_declared_tangency(self):
        c = chain([fixtures.BALL5, fixtures.PEG, box(-8.0, -8.0, 4.0, 8.0, 8.0, 12.0)])
        t = gm.Tracer(c, ('fuse', 0, 1), 2)
        values, d, edge, face = t.declared((0, 3, 4))
        self.assertLess(values, mp.mpf(10)**-30)
        self.assertLess(d, mp.mpf(10)**-30)
        self.assertTrue(edge and face)
        # The plane touches the loop without crossing it: no traced point.
        self.assertEqual(t.points(), [])

    def test_meetings_off_the_given_edges(self):
        # The sphere's and the peg's surfaces also meet below, where the peg
        # (over z in [0, 8]) is absent: the wall y = 2 crosses that meeting
        # at (+-1, 2, -2 sqrt 5), on no edge of the given result.
        c = chain([fixtures.BALL5, fixtures.PEG, box(-8.0, 2.0, -8.0, 8.0, 8.0, 10.0)])
        t = gm.Tracer(c, ('fuse', 0, 1), 2)
        s = 2*mp.sqrt(5)
        for x in (1, -1):
            X = (mp.mpf(x), mp.mpf(2), -s)
            self.assertFalse(t.on_edge(X, X, (X[0], X[1]-2, 0)))
            X = (mp.mpf(x), mp.mpf(2), s)
            self.assertTrue(t.on_edge(X, X, (X[0], X[1]-2, 0)))
            self.assertTrue(t.on_partner(X, (0, 1, 0)))

    def test_scope(self):
        self.assertTrue(gm.in_scope('sphere', 'cylinder', 'plane'))
        self.assertTrue(gm.in_scope('cylinder', 'cylinder', 'sphere'))
        self.assertFalse(gm.in_scope('sphere', 'sphere', 'plane'))
        self.assertTrue(gm.in_scope('plane', 'cone', 'cylinder'))
        self.assertFalse(gm.in_scope('plane', 'cone', 'plane'))
        self.assertFalse(gm.in_scope('plane', 'cylinder', 'sphere'))

    def test_fixture_classes(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual({c.klass for c in listed}, set(fixtures.CLASSES))
        self.assertEqual({c.partner for c in listed}, set(fixtures.PARTNERS))
        self.assertEqual(sum(c.kind == 'degenerate' for c in listed), 5)
        for c in listed:
            self.assertEqual(c.expr(), ref.chain_expr(c.op1, c.op2, c.swapped))


if __name__ == '__main__':
    unittest.main()
