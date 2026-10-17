"""Independent checks of S9e.4b.4c.2b.2's fixtures
(`generate_inside_first_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
dimple against its closed form, every part's place (each ball wider than its
rod, about its axis, the dimple's through the rod's bottom rim exactly; the
flange's spigot standing in its recess, its lug on its top and inside its
disc, no pair of its plane faces facing apart with walls along them), the
partners' places, the declared case (a ball resting on the lug's top), the
bodies' native rows, and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_inside_first_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def lug_width(x):
    """The lug's section across `x` (between the rays `y = +-3x/4`, inside
    the circle of radius 5/2, right of the chord `x = 1`)."""
    if x < 1 or x > mp.mpf(5)/2:
        return mp.mpf(0)
    return 2*min(3*x/4, mp.sqrt(mp.mpf(25)/4-x**2))


class InsideFirstReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        # A dimpled rod by heights: the rod's disc `9 pi` less the ball's
        # section `pi (25 - (w - c)^2)` where the ball reaches.
        for name, c in (('dimple', -4), ('cup', mp.mpf(-9)/2)):
            V, _ = fixtures.BODIES[name].closed()
            top = c+5

            def section(w):
                ball = 25-(w-c)**2
                return mp.pi*(9-min(9, ball)) if ball > 0 else 9*mp.pi
            got = mp.quad(section, [0, top, 6])
            self.assertLess(abs(V-got), tight, name)
        # The flange by heights: the spigot's disc below its bottom, the
        # annulus outside the recess up to `1/2`, the disc to `3/2`, the lug
        # above.
        V, _ = fixtures.BODIES['flange'].closed()
        lug = mp.quad(lug_width, [1, 2, mp.mpf(5)/2])
        got = mp.pi*4*mp.mpf(1)/2+(mp.pi*4+mp.pi*7)*mp.mpf(1)/2+mp.pi*16+lug
        self.assertLess(abs(V-got), tight)

    def test_dimple_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['dimple']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            self.assertLess(abs(V-b.closed()[0]), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        # Each ball wider than its rod, about its axis, its top above the
        # rod's bottom and below its top; the dimple's through the rod's
        # bottom rim exactly.
        for name, origin, c in (('dimple', fixtures.DIMPLE_ORIGIN, -4), ('cup', fixtures.CUP_ORIGIN, F(-9, 2))):
            b = fixtures.BODIES[name]
            rod, ball = b.ref_specs
            self.assertEqual(rod[2], fixtures.at('XY', origin))
            self.assertEqual((rod[3], rod[4]), (-1.0 if name == 'dimple' else 0.0, 6.0))
            self.assertEqual(rod[1][0].circle, (0.0, 0.0, 3.0))
            self.assertEqual(ball[1], 5.0)
            self.assertEqual(tuple(F(x) for x in ball[2][:3]), (F(origin[0]), F(origin[1]), F(origin[2])+c))
            self.assertTrue(0 < c+5 < 6)
            self.assertEqual(c**2+9 == 25, name == 'dimple')
            self.assertEqual(b.ref_expr, ('cut', 0, 1))
        # The flange: its disc, recess and spigot about one axis, the recess
        # from below the disc to inside it, the spigot from below the recess
        # into the disc above it; the lug on the disc's top inside its rim.
        f = fixtures.BODIES['flange']
        disc, recess, spigot, lug = f.ref_specs
        for s in (disc, recess, spigot):
            self.assertEqual(s[2], f.frame)
            self.assertEqual(s[1][0].circle[:2], (0.0, 0.0))
        radii = [s[1][0].circle[2] for s in (disc, recess, spigot)]
        self.assertEqual(radii, [4.0, 3.0, 2.0])
        self.assertTrue(recess[3] < disc[3] < recess[4] < disc[4])
        self.assertTrue(spigot[3] < recess[3] or spigot[3] < disc[3])
        self.assertTrue(recess[4] < spigot[4] < disc[4])
        self.assertEqual(lug[3], disc[4])
        for x, y in lug[1][0].points:
            self.assertLessEqual(F(x)**2+F(y)**2, F(25, 4))
        self.assertEqual(f.ref_expr, ('fuse', ('fuse', ('cut', 0, 1), 2), 3))
        # The profile's points and the native rows: the revolve about the
        # world's `z` through the origin, then the lug fused.
        rows = f.native_rows()
        self.assertTrue(rows[0].startswith('plane 3 3 0.5 0 -1 0 1 0 0'))
        self.assertIn('revolve 3 3 0.5 0 0 1 6.283185307179586', rows)
        self.assertIn('boolean fuse', rows)

    def test_partners_placed(self):
        # The dimple's rod across its wall (and the dimple below).
        _, (c,), f, h0, h1 = fixtures.DIMPLE_ROD
        cx, cy, r = c.circle
        d = mp.sqrt(mp.mpf(cx)**2+mp.mpf(cy)**2)
        self.assertTrue(d-r < 3 < d+r and h0 < 0 and 1 < h1 < 6)
        # The flange's rod across the spigot's wall, the recess's and the
        # lug's outer arc.
        _, (c,), f, h0, h1 = fixtures.FLANGE_ROD
        cx, cy, r = c.circle
        d = mp.sqrt(mp.mpf(cx)**2+mp.mpf(cy)**2)
        self.assertTrue(d-r < 2 < d+r and d-r < 3 < d+r and d-r < 2.5 < d+r)
        self.assertTrue(h0 < -0.5 and h1 > 2.5)
        # The tool's rod across the dimple's wall between its ends.
        _, (c,), f, h0, h1 = fixtures.ROD_DIMPLE
        cx, cy, r = c.circle
        d = mp.sqrt(mp.mpf(cx)**2+mp.mpf(cy)**2)
        self.assertTrue(d-r < 3 < d+r and 1 < h0 < h1 < 6)
        # The chain's rod across the cup's wall, its box across the hole and
        # the dimple.
        _, (c,), f, h0, h1 = fixtures.CHAIN_ROD
        cx, cy, r = c.circle
        d = mp.sqrt(mp.mpf(cx)**2+mp.mpf(cy)**2)
        self.assertTrue(d-r < 3 < d+r and h0 < 0 and h1 > 6)
        _, (b,), f, h0, h1 = fixtures.CHAIN_BOX
        xs = [p[0] for p in b.points]
        self.assertTrue(min(xs) < cx-r and cx+r < max(xs) and h0 < 0 < h1)

    def test_declared(self):
        # The touching ball: radius 1 about the lug's point `(7/4, 0)` 1 above
        # its top.
        touch = fixtures.FLANGE_TOUCH
        self.assertEqual(touch[1], 1.0)
        o = fixtures.FLANGE_ORIGIN
        self.assertEqual(tuple(touch[2][:3]), (o[0]+1.75, float(o[1]), o[2]+3.5))
        self.assertEqual(fixtures.BODIES['flange'].lug[4]+touch[1], 3.5)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'flange_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 24)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 21, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('dimple_rod_fuse', 1), ('rod_dimple_cut', 1), ('flange_cyl_common', 2),
                           ('chain_cup_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_cup_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'flange_cyl_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', ('fuse', ('cut', 0, 1), 2), 3), 4))
        swapped = next(c for c in listed if c.name == 'rod_dimple_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, 2)))
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nrevolve '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean fuse'), 1)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), list(range(len(b.ref_specs))))
        # None of the earlier steps' bodies' files.
        import generate_imported_boolean_fixtures as a
        import generate_crossed_leaves_boolean_fixtures as b
        import generate_nested_leaves_boolean_fixtures as c
        import generate_turned_profiles_boolean_fixtures as d
        import generate_primitive_chains_boolean_fixtures as e
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
