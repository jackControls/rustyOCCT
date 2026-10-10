"""Independent checks of S9e.4b.4b.2b.2b.1's fixtures
(`generate_rim_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
well against its closed form, every part's place (the pocket's square
across the plate's rim from inside it past its top with its post inside the
plate, the ring's near wall and hole over the plate's top and its long
walls across the plate's arc end), the partners' places, the declared case
(a ball resting on the well's ring) and the case list and its protocol
rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_rim_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def half_width(r, u):
    """A disc's half chord at `u` from its centre (zero outside it)."""
    return mp.sqrt(max(mp.mpf(0), r**2-u**2))


def stadium_section(length, r):
    """A stadium's area (half discs of radius `r` about `(0, 0)` and
    `(length, 0)`) by quadrature of its widths along `x`."""
    return mp.quad(lambda x: 2*half_width(r, x) if x < 0 else (2*r if x <= length else 2*half_width(r, x-length)),
                   [-r, 0, length, length+r])


def disc_section(r):
    return mp.quad(lambda y: 2*half_width(r, y), [-r, 0, r])


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


class RimLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        r = mp.mpf(5)/2
        plate = stadium_section(6, r)
        # The post: the plate's section over `[0, 1]`, then over `[1, 5/2]`
        # the plate less the pocket, the plate's part left of `x = 2` less the
        # post's disc (inside it).
        V, _ = fixtures.BODIES['post'].closed()
        left = mp.quad(lambda x: 2*half_width(r, x) if x < 0 else 2*r, [-r, 0, 2])
        post = plate*1+(plate-(left-disc_section(1)))*mp.mpf(3)/2
        self.assertLess(abs(V-post*det(fixtures.POST_FRAME)), tight)
        # The well: the plate over `[0, 1]`, the plate and the ring over
        # `[1, 2]` less their common part (the plate's widths inside `|y| <=
        # 2` for `x` in `[4, 17/2]`, less the hole's `4`), the ring `20` over
        # `[2, 3]`.
        V, _ = fixtures.BODIES['well'].closed()

        def clipped(x):
            w = 2*r if x <= 6 else 2*half_width(r, x-6)
            return min(w, mp.mpf(4))
        common = mp.quad(clipped, [4, 6, mp.mpf(15)/2, mp.mpf(17)/2])-4
        well = plate+(plate+20-common)+20
        self.assertLess(abs(V-well*det(fixtures.WELL_FRAME)), tight)
        for frame in (fixtures.POST_FRAME, fixtures.WELL_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)
            self.assertEqual(stored_axes(frame)[3], (0.0, 0.0, 1.0))

    def test_well_on_the_chained_reference(self):
        # The well alone on the chained reference (the bodies on whole cases
        # are the generator's check): the plate fused with the ring (its box
        # less its hole's), its closed form `93 + 75 pi / 8 + 25 / 4 asin
        # 3/5`.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['well']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            want = (93+75*mp.pi/8+mp.mpf(25)/4*mp.asin(mp.mpf(3)/5))*fixtures.det(fixtures.WELL_FRAME)
            self.assertLess(abs(V-want), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        r = F(5, 2)
        # The post: both prisms on the plate's frame; the square across the
        # plate's rim (its wall `x = 2` inside the plate's straight part, its
        # other walls outside the plate), from inside the plate past its top;
        # the hole's disc inside the plate and the square, its circle apart
        # from the plate's arcs and the square's walls.
        b = fixtures.BODIES['post']
        plate, tool = b.first, b.second
        self.assertEqual(plate[2], fixtures.POST_FRAME)
        self.assertEqual(tool.native[2][3:], fixtures.POST_FRAME[3:])
        self.assertEqual(tool.native[2][:2], fixtures.POST_FRAME[:2])
        x0, y0, x1, y1 = (F(p) for p in (tool.parts[0][1][0].points[0] + tool.parts[0][1][0].points[2]))
        self.assertTrue(x0 < -r and 0 < x1 < 6 and y0 < -r and y1 > r)
        cx, cy, rad = (F(c) for c in tool.parts[1][1][0].circle)
        self.assertTrue(mp.sqrt(mpq(cx**2+cy**2))+mpq(rad) < mpq(r)-mp.mpf(1)/2)
        self.assertTrue(x0 < cx-rad and cx+rad < x1 and y0 < cy-rad and cy+rad < y1)
        self.assertTrue(x1-(cx+rad) > F(1, 2))
        self.assertTrue(plate[3] < tool.native[3] < plate[4] < tool.native[4])
        # The well: the plate and the ring on one frame; the ring's near wall
        # `x = 4` and its hole inside the plate's straight part or its arc's
        # disc by at least 1/4, its long walls `|y| = 2` across the arc (whose
        # half chord there is 3/2 at `x = 15/2`, inside the ring's `[4, 10]`),
        # from inside the plate past its top.
        b = fixtures.BODIES['well']
        plate, ring = b.first, b.second
        self.assertEqual(plate[2], fixtures.WELL_FRAME)
        self.assertEqual(ring.native[2][3:], fixtures.WELL_FRAME[3:])
        sq, hole = ring.native[1]
        xs, ys = [F(p[0]) for p in sq.points], [F(p[1]) for p in sq.points]
        self.assertEqual((min(xs), max(xs), min(ys), max(ys)), (4, 10, -2, 2))
        self.assertTrue(max(ys) < r and 6+r < max(xs))
        for x, y in hole.points:
            x, y = F(x), F(y)
            gap = mpq(r)-mp.sqrt(mpq((x-6)**2+y**2)) if x > 6 else mpq(r-abs(y))
            self.assertGreater(gap, mp.mpf(1)/4)
        self.assertTrue(plate[3] < ring.native[3] < plate[4] < ring.native[4])

    def test_partners_placed(self):
        r = F(5, 2)
        # The post's rod across the pocket's wall `x = 2` inside the plate,
        # apart from the post, through every cap.
        _, (c,), frame, h0, h1 = fixtures.POST_ROD
        cx, cy, rad = (F(v) for v in c.circle)
        self.assertEqual(frame, fixtures.POST_FRAME)
        self.assertTrue(cx-rad < 2 < cx+rad and abs(cy)+rad < r and h0 < 0 and h1 > F(5, 2))
        self.assertGreater(mp.sqrt(mpq((cx+F(1, 2))**2+(cy-F(1, 4))**2)), mpq(1+rad))
        # The tool's rod across the post's circle, inside the plate's left
        # half disc.
        _, (c,), frame, h0, h1 = fixtures.ROD_POST
        cx, cy, rad = (F(v) for v in c.circle)
        d = mp.sqrt(mpq((cx+F(1, 2))**2+(cy-F(1, 4))**2))
        self.assertTrue(abs(1-mpq(rad)) < d < 1+mpq(rad))
        self.assertLess(mp.sqrt(mpq(cx**2+cy**2))+mpq(rad), mpq(r))
        # The well's rod across the hole's wall `x = 7` inside the hole's
        # `|y| < 1` band, inside the plate.
        _, (c,), frame, h0, h1 = fixtures.WELL_ROD
        cx, cy, rad = (F(v) for v in c.circle)
        self.assertEqual(frame, fixtures.WELL_FRAME)
        self.assertTrue(cx-rad < 7 < cx+rad and abs(cy)+rad < 1 and h0 < 0 and h1 > 3)
        self.assertLess(mp.sqrt(mpq((cx-6)**2+cy**2))+mpq(rad), mpq(r))
        # The chain's rod inside the plate's right part, its box across the
        # pocket's wall `x = 2` and the plate's top.
        _, (c,), frame, _, _ = fixtures.CHAIN_ROD
        cx, cy, rad = (F(v) for v in c.circle)
        self.assertTrue(2 < cx-rad and abs(cy)+rad < r)
        _, (b,), frame, h0, h1 = fixtures.CHAIN_BOX
        xs = [F(p[0]) for p in b.points]
        self.assertTrue(min(xs) < 2 < max(xs) and h0 < F(5, 2) < h1)

    def test_declared(self):
        # The touching ball: radius 1 about the ring's chart point `(17/2,
        # 0)` 1 above its top cap (its frame's normal the world's `z`
        # exactly: the centre's height exact), over the ring's material.
        touch = fixtures.WELL_TOUCH
        self.assertEqual(touch[1], 1.0)
        o = stored_axes(fixtures.WELL_FRAME)[0]
        self.assertEqual(touch[2][2], o[2]+4.0)
        self.assertEqual(fixtures.BODIES['well'].second.native[4]+touch[1], 4.0)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'well_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 18)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 15, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('post_rod_fuse', 1), ('rod_post_cut', 1), ('post_well_common', 2),
                           ('chain_post_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_post_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, ('cut', 1, 2)), 3), 4))
        both = next(c for c in listed if c.name == 'post_well_cut')
        self.assertEqual(both.expr(), ('cut', ('cut', 0, ('cut', 1, 2)), ('fuse', 3, ('cut', 4, 5))))
        swapped = next(c for c in listed if c.name == 'rod_post_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, ('cut', 2, 3))))
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
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES, t.BODIES, w.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
