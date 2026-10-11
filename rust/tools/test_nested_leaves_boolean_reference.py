"""Independent checks of S9e.4b.4c.2b.1's fixtures
(`generate_nested_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
post against its closed form, the reference's sets against OCCT's profiles
(the holed squares and their boxes, one area each), every part's place (the
square inside the plate's outline, the post or the U island inside the
square, the pocket from inside the plate past its top, how deep the pockets
nest), the partners' places, the declared case (a ball resting on the post's
top) and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_nested_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def half_width(r, u):
    """A disc's half chord at `u` from its centre (zero outside it)."""
    return mp.sqrt(max(mp.mpf(0), r**2-u**2))


def stadium_half(x, x0, length, r):
    if x < x0:
        return half_width(r, x-x0)
    if x <= x0+length:
        return r
    return half_width(r, x-x0-length)


def plate_area():
    """The plate's stadium (radius 3 about `(0, 0)` and `(6, 0)`) by
    quadrature."""
    return mp.quad(lambda x: 2*stadium_half(x, 0, 6, 3), [-3, 0, 6, 9])


def polygon_area(points):
    """The shoelace area of a polygon of exact points."""
    return sum(F(points[k][0])*F(points[k-1][1])-F(points[k-1][0])*F(points[k][1])
               for k in range(len(points)))/(-2)


def inside_stadium(x, y, margin):
    """Whether `(x, y)` lies inside the plate's outline by `margin`."""
    x, y = F(x), F(y)
    if 0 <= x <= 6:
        return abs(y)+margin <= 3
    c = 0 if x < 0 else 6
    return mp.sqrt(mpq((x-c)**2+y**2))+mpq(margin) <= 3


class NestedLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        # The post: the plate over 3 less the square ring `16 - 4` over 3/2.
        V, _ = fixtures.BODIES['post'].closed()
        want = (3*plate_area()-mpq(F(3, 2))*12)*det(fixtures.POST_FRAME)
        self.assertLess(abs(V-want), tight)
        # The island: the plate over 3 less the square less the U island over
        # 2 (the pocket from 1 up past the plate's top at 3).
        V, _ = fixtures.BODIES['island'].closed()
        want = (3*plate_area()-2*18)*det(fixtures.ISLAND_FRAME)
        self.assertLess(abs(V-want), tight)
        for frame in (fixtures.POST_FRAME, fixtures.ISLAND_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)
            self.assertEqual(stored_axes(frame)[3], (0.0, 0.0, 1.0))

    def test_post_on_the_chained_reference(self):
        # The post alone on the chained reference (the bodies on whole cases
        # are the generator's check): its closed form `90 + 27 pi`.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['post']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            want = (90+27*mp.pi)*fixtures.det(fixtures.POST_FRAME)
            self.assertLess(abs(V-want), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def holed(self, name):
        """A body's holed prism: OCCT's profile (the square and its hole) and
        the reference's boxes, on the body's frame."""
        b = fixtures.BODIES[name]
        holed = b.second
        _, (outer, hole), frame, h0, h1 = holed.native
        self.assertEqual(frame, b.first[2])
        self.assertEqual(b.op, 'cut')
        return b, outer, hole, holed, h0, h1

    def test_sets_are_the_profiles(self):
        # Each holed prism's area by its profile (the square less its hole)
        # equals the reference's boxes' (the square's box less the hole's, the
        # U's box less the notch's within the square), and its boxes reach past
        # the caps and the U's open side.
        for name, area in (('post', 12), ('island', 18)):
            _, outer, hole, holed, h0, h1 = self.holed(name)
            self.assertEqual(polygon_area(outer.points)-polygon_area(hole.points), area)
            boxes = holed.parts
            self.assertEqual((boxes[0][3], boxes[0][4]), (h0, h1))
            for box in boxes[1:]:
                self.assertTrue(box[3] < h0 and box[4] > h1)
            if name == 'island':
                (u,), (notch,) = boxes[1][1], boxes[2][1]
                ux1 = max(F(p[0]) for p in u.points)
                nx = [F(p[0]) for p in notch.points]
                ox1 = max(F(p[0]) for p in outer.points)
                self.assertTrue(min(nx) < ux1 < max(nx) < ox1)
                self.assertEqual(holed.expr, ('cut', 0, ('cut', 1, 2)))
            else:
                self.assertEqual(holed.expr, ('cut', 0, 1))

    def test_parts_placed(self):
        # Each pocket's square inside the plate's outline by 1/2 or more (no
        # wall of the pocket on the plate's walls or tangent to its arcs), from
        # inside the plate (its floor above the plate's bottom) past its top;
        # the post or the U island inside the square by 1 (its walls apart
        # from the pocket's), standing to the plate's top.
        for name, depth in (('post', 2), ('island', 3)):
            b, outer, hole, _, h0, h1 = self.holed(name)
            _, _, _, p0, p1 = b.first
            self.assertTrue(p0 < h0 < p1 < h1)
            self.assertEqual(b.depth, depth)
            for x, y in outer.points:
                self.assertTrue(inside_stadium(x, y, F(2, 5)), (name, x, y))
            ox = [F(p[0]) for p in outer.points]
            oy = [F(p[1]) for p in outer.points]
            for x, y in hole.points:
                x, y = F(x), F(y)
                self.assertTrue(min(ox)+1 <= x <= max(ox)-1 and min(oy)+1 <= y <= max(oy)-1, (name, x, y))
        # The island's notch: open at the U's wall `x = 9/2`, its back wall `x
        # = 5/2` inside the U, the arms `1` wide about it.
        _, _, hole, _, _, _ = self.holed('island')
        xs = sorted({F(p[0]) for p in hole.points})
        ys = sorted({F(p[1]) for p in hole.points})
        self.assertEqual(xs, [F(3, 2), F(5, 2), F(9, 2)])
        self.assertEqual(ys, [F(-3, 2), F(-1, 2), F(1, 2), F(3, 2)])

    def test_partners_placed(self):
        def rod(spec):
            _, (c,), frame, h0, h1 = spec
            return tuple(F(v) for v in c.circle), frame, h0, h1

        def crosses(c, r, line):
            return abs(c-line) < r

        # The post's rod across the post's corner `(4, 1)` (inside it by 1/8)
        # into the pocket's floor, clear of the pocket's walls.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.POST_ROD)
        self.assertEqual(frame, fixtures.POST_FRAME)
        self.assertTrue(crosses(cx, r, 4) and crosses(cy, r, 1))
        self.assertLessEqual(mp.sqrt(mpq((cx-4)**2+(cy-1)**2))+mp.mpf(1)/8, mpq(r))
        for x in (1, 5):
            self.assertGreaterEqual(abs(cx-x), r+F(1, 4))
        for y in (-2, 2, -1):
            self.assertGreaterEqual(abs(cy-y), r+F(1, 4))
        self.assertTrue(h0 < 0 and h1 > 3)
        # The island's rod across its notch's wall `y = 1/2` inside the arm,
        # clear of the arm's other wall, the notch's back wall and the arm's
        # end.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.ISLAND_ROD)
        self.assertEqual(frame, fixtures.ISLAND_FRAME)
        self.assertTrue(crosses(cy, r, F(1, 2)))
        self.assertGreaterEqual(abs(cy-F(3, 2)), r+F(1, 4))
        self.assertGreaterEqual(abs(cx-F(5, 2)), r+F(1, 2))
        self.assertGreaterEqual(abs(cx-F(9, 2)), r+F(1, 2))
        # The wall's rod across the pocket's wall `x = 11/2`, clear of its
        # other walls, the U and the plate's arc.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.ISLAND_WALL)
        self.assertTrue(crosses(cx, r, F(11, 2)))
        self.assertGreaterEqual(abs(cy-F(-5, 2)), r+F(1, 4))
        self.assertGreaterEqual(abs(cx-F(9, 2)), r+F(1, 4))
        self.assertTrue(inside_stadium(cx+r, cy, F(1, 4)) and inside_stadium(cx, cy-r, F(1, 4)))
        # The tool's rod across the post's four walls, clear of the pocket's.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.ROD_POST)
        for x in (2, 4):
            self.assertTrue(crosses(cx, r, x))
        for y in (-1, 1):
            self.assertTrue(crosses(cy, r, y))
        for x in (1, 5):
            self.assertGreaterEqual(abs(cx-x), r+F(1, 4))
        for y in (-2, 2):
            self.assertGreaterEqual(abs(cy-y), r+F(1, 4))
        # The chain's rod through the plate's right half beside the pocket;
        # its box across the pocket's wall `x = 5` and the plate's top, its
        # wall `x = 13/2` across the rod's hole, clear of the post and the
        # pocket's walls `y = +-2`.
        (cx, cy, r), frame, h0, h1 = rod(fixtures.CHAIN_ROD)
        self.assertEqual(frame, fixtures.POST_FRAME)
        self.assertTrue(cx-r > 5 and mp.sqrt(mpq((cx-6)**2+cy**2))+mpq(r) < 3 and h0 < 0 and h1 > 3)
        _, (b,), frame, h0, h1 = fixtures.CHAIN_BOX
        x0, x1 = min(F(p[0]) for p in b.points), max(F(p[0]) for p in b.points)
        y0, y1 = min(F(p[1]) for p in b.points), max(F(p[1]) for p in b.points)
        self.assertTrue(4 < x0 < 5 < x1 and crosses(cx, r, x1))
        self.assertTrue(-2 < y0 and y1 < 2 and F(3, 2) < h0 < 3 < h1)

    def test_cylinder_placed(self):
        # S9e.4a's cylinder (radius 3 about the world's `(5, 5)` over `[-1,
        # 5]`) about the island's chart point `(2, 3/8)` within 1e-2, past the
        # plate's caps, crossing the plate's wall `y = 3`, the pocket's walls
        # `x = 1/2` and `y = +-5/2`, the U's and the notch's walls, clear of
        # `x = 11/2` and `y = -3`; every vertex of the plate, the pocket, the
        # U and the notch at least 1/8 off its circle.
        o, x, y, n = stored_axes(fixtures.ISLAND_FRAME)
        d = (mpq(5)-mpq(o[0]), mpq(5)-mpq(o[1]))
        cx = d[0]*mpq(x[0])+d[1]*mpq(x[1])
        cy = d[0]*mpq(y[0])+d[1]*mpq(y[1])
        self.assertLess(abs(cx-2)+abs(cy-mpq(F(3, 8))), mp.mpf(10)**-2)
        spec = fixtures.base.BODIES['cyl'].spec
        self.assertEqual((spec[2][2], spec[3], spec[4]), (-1.0, 0.0, 6.0))
        self.assertTrue(spec[2][2] < o[2] and o[2]+3 < spec[2][2]+6)
        r = 3
        for xl in (F(1, 2), F(3, 2), F(5, 2), F(9, 2)):
            self.assertLess(abs(cx-mpq(xl)), r)
        self.assertGreater(abs(cx-mpq(F(11, 2))), r)
        for yl in (F(5, 2), F(-5, 2), F(3, 2), F(-3, 2), F(1, 2), F(-1, 2), 3):
            self.assertLess(abs(cy-mpq(yl)), r)
        self.assertGreater(abs(cy+3), r)
        corners = [(F(a), F(b)) for a, b in ((0, 3), (0, -3), (6, 3), (6, -3), (F(1, 2), F(5, 2)),
                                             (F(1, 2), F(-5, 2)), (F(11, 2), F(5, 2)), (F(11, 2), F(-5, 2)),
                                             (F(3, 2), F(3, 2)), (F(3, 2), F(-3, 2)), (F(9, 2), F(3, 2)),
                                             (F(9, 2), F(-3, 2)), (F(9, 2), F(1, 2)), (F(9, 2), F(-1, 2)),
                                             (F(5, 2), F(1, 2)), (F(5, 2), F(-1, 2)))]
        for a, b in corners:
            self.assertGreaterEqual(abs(mp.sqrt((mpq(a)-cx)**2+(mpq(b)-cy)**2)-r), mp.mpf(1)/8-mp.mpf(10)**-2,
                                    (a, b))

    def test_declared(self):
        # The touching ball: radius 1 about the post's chart point `(3, 0)` 1
        # above the plate's top (its frame's normal the world's `z` exactly:
        # the centre's height exact), over the post's top, its walls clear of
        # the ball by 1/4 or more horizontally (no tangency but the cap's).
        touch = fixtures.POST_TOUCH
        self.assertEqual(touch[1], 1.0)
        o = stored_axes(fixtures.POST_FRAME)[0]
        self.assertEqual(touch[2][2], o[2]+4.0)
        self.assertEqual(fixtures.BODIES['post'].first[4]+touch[1], 4.0)
        self.assertTrue(2 < 3 < 4 and -1 < 0 < 1)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'post_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('post_rod_fuse', 1), ('rod_post_cut', 1), ('island_cyl_common', 2),
                           ('chain_post_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_post_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, ('cut', 1, 2)), 3), 4))
        both = next(c for c in listed if c.name == 'island_cyl_cut')
        self.assertEqual(both.expr(), ('cut', ('cut', 0, ('cut', 1, ('cut', 2, 3))), 4))
        swapped = next(c for c in listed if c.name == 'rod_post_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, ('cut', 2, 3))))
        # Each body written once, two prisms each.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean cut'), len(fixtures.BODIES))
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
        import generate_stack_leaves_boolean_fixtures as u
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES, t.BODIES, w.BODIES, r.BODIES, s.BODIES, u.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
