"""Independent checks of S9e.4b.4b.1's fixtures
(`generate_primitive_chains_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
cup and the dome against their closed forms, the bodies' coaxial
primitives (every meeting a circle: the dome's and the pin's rims on their
ball, the bead's at `+-4`), the declared bodies' tangency and fillets, and
the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_primitive_chains_boolean_fixtures as fixtures
from curve_surface_reference import stored_axes


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def disc_quad(radius, lo, hi, breaks=()):
    """The volume of a solid of revolution of section radius `radius(z)`
    between `lo` and `hi` by quadrature of its discs."""
    return mp.quad(lambda z: mp.pi*radius(z)**2, [lo, *breaks, hi])


class PrimitiveChainsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-30
        # The dome: the frustum of radius `7/2 - z/4` inside the ball of
        # radius 5 about `(0, 0, -2)`, up to its top.
        ball = lambda z: mp.sqrt(max(mp.mpf(0), 25-(z+2)**2))
        cone = lambda z: mp.mpf(7)/2-z/4
        V, _ = fixtures.BODIES['dome'].closed()
        self.assertLess(abs(V-disc_quad(lambda z: min(ball(z), cone(z)), 0, 3, (2,))), tight)
        # The pin, the binary64 radius 1.4 inside that ball.
        r = mp.mpf(1.4)
        top = -2+mp.sqrt(25-r**2)
        V, _ = fixtures.BODIES['pin'].closed()
        self.assertLess(abs(V-disc_quad(lambda z: min(ball(z), r), 0, 3, (top,))), tight)
        # The knob: the ball of radius 2 about its centre fused with the rod
        # of the binary64 radius 1.2 from the centre to 5 above it.
        r = mp.mpf(1.2)
        knob = lambda w: max(mp.sqrt(max(mp.mpf(0), 4-w**2)), r if 0 <= w <= 5 else mp.mpf(0))
        rim = mp.sqrt(4-r**2)
        V, _ = fixtures.BODIES['knob'].closed()
        self.assertLess(abs(V-disc_quad(knob, -2, 5, (0, rim, 2))), tight)
        # The shaft's and the cup's sections in their chart, times the stored
        # axes' determinant (a unit within rounding).
        V, _ = fixtures.BODIES['shaft'].closed()
        chart = disc_quad(lambda w: 3 if w <= 4 else mp.mpf(1.5), 0, 9, (4,))
        self.assertLess(abs(V-chart*fixtures.det(fixtures.SHAFT_FRAME)), tight)
        self.assertLess(abs(fixtures.det(fixtures.SHAFT_FRAME)-1), mp.mpf(10)**-15)
        V, _ = fixtures.BODIES['cup'].closed()
        chart = mp.quad(lambda w: mp.pi*(9-(4 if w >= 1 else 0)), [0, 1, 5])
        self.assertLess(abs(V-chart*fixtures.det(fixtures.CUP_FRAME)), tight)

    def test_bodies_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 48
        try:
            for name in ('cup', 'dome'):
                b = fixtures.BODIES[name]
                chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
                V, _, _ = chain.measures(b.ref_expr)
                V0, _ = b.closed()
                self.assertLess(abs(V-V0), mp.mpf(10)**-28, name)
        finally:
            ref.SCAN = scan

    def test_coaxial_primitives(self):
        # Every body's two primitives share an axis (their frames' normals
        # and origins on one line), so each meeting is a circle.
        for b in fixtures.BODIES.values():
            if b.name == 'rounded':
                continue
            frames = [fixtures.pieces.spec_frame(s) for s in b.specs()]
            o0, _, _, n0 = (tuple(F(c) for c in a) for a in stored_axes(frames[0]))
            for fr in frames[1:]:
                o, _, _, n = (tuple(F(c) for c in a) for a in stored_axes(fr))
                cross = (n0[1]*n[2]-n0[2]*n[1], n0[2]*n[0]-n0[0]*n[2], n0[0]*n[1]-n0[1]*n[0])
                self.assertEqual(cross, (0, 0, 0), b.name)
                d = tuple(o[i]-o0[i] for i in range(3))
                off = (n0[1]*d[2]-n0[2]*d[1], n0[2]*d[0]-n0[0]*d[2], n0[0]*d[1]-n0[1]*d[0])
                self.assertEqual(off, (0, 0, 0), b.name)
        # The dome's cone meets its ball at radius 3, 2 above the base,
        # exactly; the bead's rims at `+-4`.
        self.assertEqual(F(7, 2)-F(2, 4), 3)
        self.assertEqual(3**2+(2+2)**2, 25)
        self.assertEqual(3**2+4**2, 25)
        # The pin's disc on the dome's base plane, inside the dome's disc.
        pin, dome = fixtures.BODIES['pin'], fixtures.BODIES['dome']
        self.assertEqual(pin.first[3], 0.0)
        self.assertEqual(dome.first[4][:3], (0.0, 0.0, 0.0))
        self.assertLess(pin.first[1][0].circle[2], dome.first[1])
        self.assertEqual(pin.second, dome.second)

    def test_declared_bodies(self):
        # The capsule's ball of its rod's radius about its top's centre:
        # tangent along the rim.
        rod, ball = fixtures.BODIES['capsule'].specs()
        self.assertEqual(rod[1][0].circle[2], ball[1])
        self.assertEqual(ball[2][:3], (0.0, 0.0, rod[4]))
        # The rounded square's arcs tangent to its lines (each arc's centre
        # a radius inside both lines it joins).
        square = fixtures.BODIES['rounded'].first[1][0]
        for k, seg in enumerate(square.segments):
            if seg is None:
                continue
            cx, cy, r, ccw = seg
            self.assertTrue(ccw)
            for p in (square.points[k], square.points[(k+1) % len(square.points)]):
                self.assertEqual((F(p[0])-F(cx))**2+(F(p[1])-F(cy))**2, F(r)**2)
        self.assertEqual({b.name for b in fixtures.BODIES.values() if b.kind}, {'capsule', 'rounded'})
        self.assertEqual(fixtures.BODIES['capsule'].kind, 'degenerate')
        self.assertEqual(fixtures.BODIES['rounded'].kind, 'unsupported')

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 24)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 21, 'empty': 0, 'degenerate': 0, 'unsupported': 3})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('shaft_box_fuse', 'knob_box_cut', 'dome_pin_common', 'chain_dome_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/chain_'), want, name)
            self.assertEqual(nat.count('\nbrep imported/chain_'), want, name)
        chain = next(c for c in listed if c.name == 'chain_dome_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('common', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'dome_pin_cut')
        self.assertEqual(both.expr(), ('cut', ('common', 0, 1), ('common', 2, 3)))
        swapped = next(c for c in listed if c.name == 'knob_box_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('fuse', 1, 2)))
        # Every body written once: two primitives about a Boolean.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        # The capsule given to no case (its tangency's cost), refused on
        # import by the kernel's tests.
        self.assertEqual(used, set(fixtures.BODIES)-{'capsule'})
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), [0, 1])
            self.assertGreaterEqual(len(b.curved_kinds()), 2, b.name)
        # None of the earlier steps' bodies' files.
        import generate_imported_boolean_fixtures as a
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_polyhedra_boolean_fixtures as c
        import generate_imported_pieces_boolean_fixtures as d
        import generate_one_sphere_boolean_fixtures as e
        import generate_one_sphere_incidence_boolean_fixtures as f
        import generate_piece_forms_boolean_fixtures as g
        import generate_piece_trees_boolean_fixtures as h
        import generate_imported_joints_boolean_fixtures as i
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
