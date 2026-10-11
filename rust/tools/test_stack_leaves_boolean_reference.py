"""Independent checks of S9e.4b.4b.2b.2b.3a's fixtures
(`generate_stack_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
bar against its closed form, every stack's place (the upper prism on the
lower one's top on one frame, the two caps exactly on one plane, the
outlines' crossings rational and no wall on another's circle's tangent; the
step's caps both partly exposed, the bar across the plate's arc end), the
partners' places, the declared case (a ball resting on the bar's top) and
the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_stack_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def half_width(r, u):
    """A disc's half chord at `u` from its centre (zero outside it)."""
    return mp.sqrt(max(mp.mpf(0), r**2-u**2))


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def stadium_half(x, x0, length, r):
    """A stadium's half width at `x` (half discs of radius `r` about `x0` and
    `x0 + length` on its axis)."""
    if x < x0:
        return half_width(r, x-x0)
    if x <= x0+length:
        return r
    return half_width(r, x-x0-length)


def stadium_area(x0, length, r):
    return mp.quad(lambda x: 2*stadium_half(x, mpq(x0), mpq(length), mpq(r)),
                   [mpq(x0)-mpq(r), mpq(x0), mpq(x0)+mpq(length), mpq(x0)+mpq(length)+mpq(r)])


class StackLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        # The step: the box `[-3, 3] x [-2, 2]` over `3/2` and the stadium of
        # radius 5/4 about `(1, 1)` and `(5, 1)` over `3/2`.
        V, _ = fixtures.BODIES['step'].closed()
        want = (24*mpq(F(3, 2))+mpq(F(3, 2))*stadium_area(1, 4, F(5, 4)))*det(fixtures.STEP_FRAME)
        self.assertLess(abs(V-want), tight)
        # The bar: the plate of radius 25/8 about `(0, 0)` and `(6, 0)` over 2
        # and the box `11/2 x 15/4` over 1.
        V, _ = fixtures.BODIES['bar'].closed()
        want = (2*stadium_area(0, 6, F(25, 8))+mpq(F(11, 2))*mpq(F(15, 4)))*det(fixtures.BAR_FRAME)
        self.assertLess(abs(V-want), tight)
        for frame in (fixtures.STEP_FRAME, fixtures.BAR_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)
            self.assertEqual(stored_axes(frame)[3], (0.0, 0.0, 1.0))

    def test_bar_on_the_chained_reference(self):
        # The bar alone on the chained reference (the bodies on whole cases
        # are the generator's check): the plate fused with the box standing on
        # it, its closed form `75 + 165/8 + 625 pi / 32`.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['bar']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            want = (75+mpq(F(165, 8))+625*mp.pi/32)*fixtures.det(fixtures.BAR_FRAME)
            self.assertLess(abs(V-want), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def stack(self, name):
        """A body's two prisms' native rows: frames, profiles and heights; the
        upper one's bottom the lower one's top, bit for bit, on one frame."""
        b = fixtures.BODIES[name]
        lower, upper = (s.native if hasattr(s, 'native') else s for s in (b.first, b.second))
        self.assertEqual(b.op, 'fuse')
        self.assertEqual(lower[2], upper[2])
        self.assertEqual(lower[4], upper[3])
        self.assertTrue(lower[3] < lower[4] < upper[4])
        return lower, upper

    def test_stacks_placed(self):
        # The step: the box `[-3, 3] x [-2, 2]` and the stadium of radius 5/4
        # about `(1, 1)` and `(5, 1)`: the outlines cross at `(3, -1/4)` (the
        # box's wall and the stadium's line) and `(1/4, 2)` (the box's wall `y
        # = 2` and the left arc, `(x - 1)^2 + 1 = 25/16`), its other crossing
        # with that circle `x = 7/4` right of the arc's end (inside the
        # stadium's straight part); each circle's distance from every wall of
        # the box's carrier at least 1/4 from its radius; the box's corners
        # off the stadium's outline; both caps partly exposed.
        lower, upper = self.stack('step')
        r = F(5, 4)
        for c in ((F(1), F(1)), (F(5), F(1))):
            for d in (abs(c[0]-3), abs(c[0]+3), abs(c[1]-2), abs(c[1]+2)):
                self.assertGreaterEqual(abs(d-r), F(1, 4))
        self.assertEqual((F(1, 4)-1)**2+(2-1)**2, r**2)
        self.assertEqual((F(7, 4)-1)**2+(2-1)**2, r**2)
        self.assertTrue(1 < F(7, 4) < 5 and F(-1, 4) < 2 < F(9, 4))
        for corner in ((3, 2), (3, -2), (-3, 2), (-3, -2)):
            x, y = (F(v) for v in corner)
            line = 1 <= x <= 5 and y in (F(-1, 4), F(9, 4))
            arc = min((x-1)**2+(y-1)**2, (x-5)**2+(y-1)**2) == r**2
            self.assertFalse(line or arc, corner)
        # A point of the box's top outside the stadium, of the stadium's bottom
        # outside the box.
        self.assertGreater((F(-2)-1)**2+(F(-1)-1)**2, r**2)
        self.assertTrue(F(11, 2) > 3 and (F(11, 2)-5)**2 < r**2)
        # The bar: the plate of radius 25/8 about `(0, 0)` and `(6, 0)` and the
        # box `[5, 21/2] x [-15/8, 15/8]`: the arc crossing `y = +-15/8` at `x =
        # 17/2`; the box's wall `x = 5` inside the plate, `x = 21/2` beyond its
        # end `x = 73/8`; the plate's lines `y = +-25/8` clear of the box.
        lower, upper = self.stack('bar')
        R = F(25, 8)
        self.assertEqual((F(17, 2)-6)**2+F(15, 8)**2, R**2)
        self.assertTrue(F(5) < 6 and F(21, 2) > 6+R and F(15, 8) < R)

    def test_partners_placed(self):
        def rod(spec):
            _, (c,), frame, h0, h1 = spec
            return tuple(mpq(v) for v in c.circle), frame, h0, h1
        # The step's rod inside both outlines across the shared plane, its
        # seam's middle on it: inside the box and the stadium's straight part.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.STEP_ROD)
        self.assertEqual(frame, fixtures.STEP_FRAME)
        self.assertTrue(abs(cx)+r < 3 and abs(cy)+r < 2 and 1 < cx-r and cx+r < 5)
        self.assertTrue(-mpq(F(1, 4)) < cy-r and cy+r < mpq(F(9, 4)))
        self.assertEqual((h0+h1)/2, 1.5)
        self.assertTrue(h0 < 0 and h1 > 3)
        # The wall's rod across the box's wall `x = 3` under the stadium's
        # straight part, clear of its arcs and lines.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.STEP_WALL)
        self.assertTrue(cx-r < 3 < cx+r and 1 < cx-r and cx+r < 5)
        self.assertTrue(-mpq(F(1, 4)) < cy-r and cy+r < 2)
        for c in ((1, 1), (5, 1)):
            self.assertGreater(mp.sqrt((cx-c[0])**2+(cy-c[1])**2), mpq(F(5, 4))+r)
        # The bar's rod inside both outlines (the plate's right half disc and
        # the box) across the shared plane, its seam's middle on it, through
        # the plate's bottom and the box's top.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.BAR_ROD)
        self.assertEqual(frame, fixtures.BAR_FRAME)
        self.assertLess(mp.sqrt((cx-6)**2+cy**2)+r, mpq(F(25, 8)))
        self.assertTrue(5 < cx-r and cx+r < mpq(F(21, 2)) and abs(cy)+r < mpq(F(15, 8)))
        self.assertTrue((h0+h1)/2 == 2.0 and h0 < 0 and h1 > 3)
        # The tool's rod across the box's walls `x = 21/2` and `y = 15/8`, its
        # corner inside the rod by 1/8, clear of the plate.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.ROD_BAR)
        self.assertTrue(cx-r < mpq(F(21, 2)) < cx+r and cy-r < mpq(F(15, 8)) < cy+r)
        self.assertLessEqual(mp.sqrt((cx-mpq(F(21, 2)))**2+(cy-mpq(F(15, 8)))**2)+mp.mpf(1)/8, r)
        self.assertGreater(mp.sqrt((cx-6)**2+cy**2), mpq(F(25, 8))+r)
        # The chain's rod through the step's box beside the stadium (left of its
        # left arc `x >= -1/4`); its box across the box's wall `x = 3`, the
        # shared plane, the stadium's line `y = -1/4` and its right arc (about
        # `(5, 1)`, radius 5/4: crossing `x = 6` and `y = 1/2`, clear of `y =
        # -1` and `x = 5/2`).
        (cx, cy, r), frame, h0, h1 = rod(fixtures.CHAIN_ROD)
        self.assertEqual(frame, fixtures.STEP_FRAME)
        self.assertTrue(-3 < cx-r and cx+r < -mpq(F(1, 4)) and abs(cy)+r < 2 and h0 < 0 and h1 > 3)
        _, (b,), frame, h0, h1 = fixtures.CHAIN_BOX
        x0, x1 = min(F(p[0]) for p in b.points), max(F(p[0]) for p in b.points)
        y0, y1 = min(F(p[1]) for p in b.points), max(F(p[1]) for p in b.points)
        self.assertTrue(x0 < 3 < x1 and y0 < F(-1, 4) < y1 and h0 < 1.5 < h1 < 3)
        r = F(5, 4)
        self.assertTrue(abs(x1-5) < r and abs(y1-1) < r)
        self.assertTrue(abs(y0-1) > r+F(1, 4) and abs(x0-5) > r+F(1, 4))
        self.assertLess(cx, mpq(x0))

    def test_declared(self):
        # The touching ball: radius 1 about the bar's chart point `(7, 0)` 1
        # above the box's top cap (its frame's normal the world's `z` exactly:
        # the centre's height exact), over the box's material, the plate's arcs'
        # cylinders and the box's walls clear of the ball by 1/4 or more (no
        # tangency but the cap's).
        for cx, cy, r in ((0, 0, F(25, 8)), (6, 0, F(25, 8))):
            d = mp.sqrt(mpq((7-cx)**2+cy**2))
            self.assertTrue(d+1+mp.mpf(1)/4 <= mpq(r) or d-1-mp.mpf(1)/4 >= mpq(r), (cx, cy, r))
        for wall in (7-5, F(21, 2)-7, F(15, 8)):
            self.assertGreaterEqual(wall, 1+F(1, 4))
        touch = fixtures.BAR_TOUCH
        self.assertEqual(touch[1], 1.0)
        o = stored_axes(fixtures.BAR_FRAME)[0]
        self.assertEqual(touch[2][2], o[2]+4.0)
        self.assertEqual(fixtures.BODIES['bar'].second[4]+touch[1], 4.0)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'bar_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('step_rod_fuse', 1), ('rod_bar_cut', 1), ('step_bar_common', 2),
                           ('chain_step_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_step_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'step_bar_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), ('fuse', 2, 3)))
        swapped = next(c for c in listed if c.name == 'rod_bar_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('fuse', 1, 2)))
        # Each body written once, two prisms each.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nprism '), 2*len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), list(range(len(b.ref_specs))))
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
        import generate_primitive_chains_boolean_fixtures as j
        import generate_prism_leaves_boolean_fixtures as k
        import generate_plane_parts_boolean_fixtures as m
        import generate_polyhedra_curved_boolean_fixtures as p
        import generate_deep_trees_boolean_fixtures as t
        import generate_two_leaves_boolean_fixtures as w
        import generate_rim_leaves_boolean_fixtures as r
        import generate_flush_leaves_boolean_fixtures as s
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES, t.BODIES, w.BODIES, r.BODIES, s.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
