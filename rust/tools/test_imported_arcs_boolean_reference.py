"""Independent checks of S9e.4b.1's fixtures
(`generate_imported_arcs_boolean_fixtures.py`): the profiles' arcs end on
their circles in the constructions (exactly, in rationals), the chained
reference on a turned slot against results known without it (the slot's own
closed form in the `TILT` frame; a box cutting it through its straight part
by a direct quadrature of its sections), the off-circle test on points known
to be on and off, and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_imported_arcs_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import square


def mpq(x):
    return mp.mpf(x.numerator)/x.denominator


class ImportedArcsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_constructions_end_their_arcs_on_their_circles(self):
        # Each body OCCT was given ends every arc on its circle exactly (its
        # profile's binary64 points): only the files' roundings leave them.
        for b in fixtures.BODIES.values():
            for p, circles in fixtures.arc_ends(b.spec):
                for cx, cy, r in circles:
                    self.assertEqual((F(p[0])-F(cx))**2+(F(p[1])-F(cy))**2, F(r)**2, (b.name, p))

    def test_slot_measures(self):
        # The slot (a stadium of straight length 6 and radius 2, height 5)
        # alone: the chained reference's volume and area its closed form,
        # `(24 + 4 pi) 5` (times the stored axes' determinant) and `2 (24 +
        # 4 pi) + (12 + 4 pi) 5`.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            slot = fixtures.construction(fixtures.BODIES['slot'].spec, 91)
            box = fixtures.construction(fixtures.prism([square(30.0, 30.0, 31.0, 31.0)], fixtures.at('XY', (0, 0, 0)),
                                                       0.0, 1.0), 92)
            chain = ref.Chain([slot, box])
            v, _, a = chain.measures(0)
            # The stored axes are binary64 roundings of an orthonormal frame:
            # volumes scale by their determinant exactly, areas within
            # rounding.
            _, x, y, n = (tuple(F(c) for c in w) for w in fixtures.stored_axes(fixtures.BODIES['slot'].spec[2]))
            det = (x[0]*(y[1]*n[2]-y[2]*n[1])-x[1]*(y[0]*n[2]-y[2]*n[0])+x[2]*(y[0]*n[1]-y[1]*n[0]))
            self.assertLess(abs(v-(24+4*mp.pi)*5*mpq(det)), mp.mpf(10)**-30)
            self.assertLess(abs(a-(2*(24+4*mp.pi)+(12+4*mp.pi)*5)), mp.mpf(10)**-14)
        finally:
            ref.SCAN = scan

    def test_slot_cut_by_its_sections(self):
        # The slot less the box `x in [7, 12]` (the arc end beyond its
        # straight part's x = 8 cut too): its volume by the slot's sections
        # normal to its `x` axis (the world's `x`), each the stadium's chord
        # height times 5, outside `[7, 12] x [0, 10] x [1, 2.5]` the full
        # section, inside it the part of the section's rectangle of `(v, w)`
        # outside the slab `1 <= -0.6 v + 0.8 w <= 2.5` within `0 <= y <= 10`
        # (every point of the slot has `y` in `[1.4, 7.6]`).
        scan, ref.SCAN = ref.SCAN, 48
        try:
            chain = ref.Chain([fixtures.construction(fixtures.BODIES['slot'].spec, 91),
                               fixtures.construction(fixtures.BOX, 92)])
            V, _, _ = chain.measures(('cut', 0, 1))
        finally:
            ref.SCAN = scan

        def chord(u):
            # The stadium's half height at u (|u| <= 5).
            return 2 if abs(u) <= 3 else mp.sqrt(4-(abs(u)-3)**2)

        def removed(h):
            # The area of {|v| <= h, 0 <= w <= 5, 1 <= -0.6 v + 0.8 w <= 2.5}:
            # w between (1 + 0.6 v) / 0.8 (or 0, below v = -5/3) and (2.5 +
            # 0.6 v) / 0.8, both within [0, 5] for |v| <= 2.
            f = lambda v: ((mp.mpf(5)/2+mp.mpf(3)/5*v)-max(0, 1+mp.mpf(3)/5*v))/(mp.mpf(4)/5)
            return mp.quad(f, [-h, -mp.mpf(5)/3, h] if h > mp.mpf(5)/3 else [-h, h])
        full = mp.quad(lambda u: 2*chord(u)*5, [-5, -3, 3, 5])
        # Breakpoints where the chord's height passes v = 5/3.
        cut = mp.quad(lambda u: removed(chord(u)), [2, 3, 3+mp.sqrt(11)/3, 5])
        self.assertLess(abs(V-(full-cut)), mp.mpf(10)**-12)

    def test_off_circle_test(self):
        # The off-circle test on the slot's construction: its own profile
        # points (on), and those points moved by an ulp (off).
        b = fixtures.BODIES['slot']
        local = fixtures.local_map(b.spec[2])
        o, x, y, n = (tuple(F(c) for c in v) for v in fixtures.stored_axes(b.spec[2]))
        for p, circles in fixtures.arc_ends(b.spec):
            world = tuple(o[i]+F(p[0])*x[i]+F(p[1])*y[i] for i in range(3))
            u = local(world)
            self.assertEqual((u[0], u[1]), (F(p[0]), F(p[1])))
        cx, cy, r = 3.0, 0.0, 2.0
        on = (F(5.0)-F(cx))**2+(F(0.0)-F(cy))**2 == F(r)**2
        off = (F(5.000000000000001)-F(cx))**2+(F(0.0)-F(cy))**2 == F(r)**2
        self.assertTrue(on)
        self.assertFalse(off)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 36)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        # The lens's 3 solid since S9e.4b.4a.
        self.assertEqual(kinds, {'solid': 30, 'empty': 0, 'degenerate': 6, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('slot_box_fuse', 'box_slot_cut', 'both_common', 'chain_slot_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.specs)
            self.assertEqual(enc.count('\nbrep imported/'), want, name)
            self.assertEqual(nat.count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_slot_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        # Every body written once, by its rows: a plane, a path and a prism.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nprism '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.specs if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of S9e.4a's bodies' names.
        import generate_imported_boolean_fixtures as earlier
        self.assertFalse(set(fixtures.BODIES) & set(earlier.BODIES))


if __name__ == '__main__':
    unittest.main()
