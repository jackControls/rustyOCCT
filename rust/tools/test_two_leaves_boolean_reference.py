"""Independent checks of S9e.4b.4b.2b.2a's fixtures
(`generate_two_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
cross against its closed form, every part's place (the stadium through the
box from below its bottom past its top inside its footprint, the arm through
the stadium's arc walls between its caps, the pocket inside the plate from inside it
past its top), the partners' places, the declared case (a ball resting on
the stadium's top cap), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_two_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def stadium_width(length, r, y):
    """A stadium's width at `y` about its centre line (half discs of radius
    `r` about `(0, 0)` and `(length, 0)`)."""
    return length+2*mp.sqrt(max(mp.mpf(0), r**2-y**2))


def chart(frame, point):
    """A world point's coordinates in a frame's stored axes, exactly (Cramer's
    rule)."""
    o, x, y, n = (tuple(F(c) for c in a) for a in stored_axes(frame))
    d = tuple(F(point[i])-o[i] for i in range(3))

    def det3(a, b, c):
        return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    D = det3(x, y, n)
    return det3(d, y, n)/D, det3(x, d, n)/D, det3(x, y, d)/D


def stadium_gap(length, r, u, v):
    """A point's distance inside a stadium's outline (negative outside)."""
    s = min(max(u, 0), length)
    return r-mp.sqrt((u-s)**2+v**2)


class TwoLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        half = mp.mpf(3)/2
        # The through: the box's section `48` over `[0, 2]`, the stadium's
        # alone below and above it, over its heights' outside parts.
        V, _ = fixtures.BODIES['through'].closed()
        stadium = mp.quad(lambda y: stadium_width(4, half, y), [-half, 0, half])
        through = mp.quad(lambda w: 48 if 0 <= w <= 2 else stadium, [-1, 0, 2, 3])
        self.assertLess(abs(V-through*det(fixtures.THROUGH_FRAME)), tight)
        # The cross: by heights, the stadium's section to 1, its union with
        # the arm's `24` from 1 to 2 (the arm's band `|y| <= 1` across it),
        # the stadium alone from 2 to 3.
        V, _ = fixtures.BODIES['cross'].closed()
        s = mp.quad(lambda y: stadium_width(4, 2, y), [-2, 0, 2])
        band = mp.quad(lambda y: stadium_width(4, 2, y), [-1, 0, 1])
        cross = s+(s+24-band)+s
        self.assertLess(abs(V-cross), tight)
        # The pocketed plate: the stadium's section over `[0, 5/2]` less the
        # pocket's square over `[1, 5/2]`, each in its chart.
        V, _ = fixtures.BODIES['pocketed'].closed()
        plate = mp.quad(lambda y: stadium_width(6, mp.mpf(5)/2, y), [-mp.mpf(5)/2, 0, mp.mpf(5)/2])
        want = mp.mpf(5)/2*plate*det(fixtures.POCKETED_FRAME)-9*mp.mpf(3)/2*det(fixtures.POCKETED_HOLE)
        self.assertLess(abs(V-want), tight)
        for frame in (fixtures.THROUGH_FRAME, fixtures.POCKETED_FRAME, fixtures.POCKETED_HOLE):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)

    def test_cross_on_the_chained_reference(self):
        # The cross alone on the chained reference (the bodies on whole cases
        # are the generator's check): the stadium fused with the arm, its
        # closed form `64 + 32 pi / 3 - 2 sqrt 3`.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['cross']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            self.assertLess(abs(V-(64+32*mp.pi/3-2*mp.sqrt(3))), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        # The through: the stadium on the box's axes inside its footprint
        # (its outline at least 1/2 from the box's walls), from below the box
        # past its top.
        box, stadium = fixtures.BODIES['through'].first, fixtures.BODIES['through'].second
        self.assertEqual(box[2], fixtures.THROUGH_FRAME)
        self.assertEqual(stadium[2][3:], fixtures.THROUGH_FRAME[3:])
        u, v, w = chart(fixtures.THROUGH_FRAME, stadium[2][:3])
        self.assertLess(abs(w), F(1, 10**15))
        eps = F(1, 10**12)
        self.assertTrue(u-F(3, 2) > F(1, 2)-eps and u+4+F(3, 2) < 8-F(1, 2)+eps)
        self.assertTrue(v-F(3, 2) > F(3, 2)-eps and v+F(3, 2) < 6-F(3, 2)+eps)
        self.assertTrue(stadium[3] < box[3] and box[4] < stadium[4])
        # The cross: the arm through the stadium's arc walls, past its
        # outline along `x`, between its caps, inside its line walls.
        st, arm = fixtures.BODIES['cross'].first, fixtures.BODIES['cross'].second
        self.assertEqual((st[3], st[4]), (0.0, 3.0))
        _, (b,), frame, h0, h1 = arm
        self.assertEqual(frame[3:], (1.0, 0.0, 0.0, 0.0, 1.0, 0.0))
        ys = [p[0] for p in b.points]
        zs = [p[1] for p in b.points]
        self.assertTrue(-2 < min(ys) and max(ys) < 2 and 0 < min(zs) and max(zs) < 3)
        self.assertTrue(frame[0]+h0 < -2 and frame[0]+h1 > 6)
        # The pocketed plate: the pocket's corners inside the plate's outline
        # (in its chart) by at least 1/4, from inside the plate past its top.
        plate, hole = fixtures.BODIES['pocketed'].first, fixtures.BODIES['pocketed'].second
        _, (sq,), hframe, h0, h1 = hole
        o, x, y, n = stored_axes(hframe)
        for cx, cy in sq.points:
            corner = tuple(F(o[i])+F(cx)*F(x[i])+F(cy)*F(y[i]) for i in range(3))
            u, v, w = chart(fixtures.POCKETED_FRAME, corner)
            self.assertGreater(stadium_gap(6, mp.mpf(5)/2, mp.mpf(u.numerator)/u.denominator, mp.mpf(v.numerator)/v.denominator), mp.mpf(1)/4)
            self.assertEqual(w, 0)
        self.assertTrue(plate[3] < h0 < plate[4] < h1)
        self.assertEqual(stored_axes(hframe)[3], (0.0, 0.0, 1.0))
        self.assertEqual(stored_axes(fixtures.POCKETED_FRAME)[3], (0.0, 0.0, 1.0))

    def test_partners_placed(self):
        # The through's rod across the stadium's arc end and the box's wall
        # `x = 8`, through every cap.
        _, (c,), frame, h0, h1 = fixtures.THROUGH_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame, fixtures.THROUGH_FRAME)
        d = abs(F(cx)-6)
        self.assertTrue(abs(F(3, 2)-F(r)) < d < F(3, 2)+F(r) and cx-r < 8 < cx+r and h0 < -1 and h1 > 3)
        self.assertEqual(cy, 3.0)
        # The cross's rod inside the stadium's arc end, across the arm's wall
        # `y = 1`.
        _, (c,), frame, h0, h1 = fixtures.CROSS_ROD
        cx, cy, r = c.circle
        self.assertLess(mp.sqrt((cx-4)**2+cy**2)+r, 2)
        self.assertTrue(cy-r < 1 < cy+r and h0 < 0 and h1 > 3)
        # The tool's rod across the stadium's line wall `y = 2` beside the arm.
        _, (c,), frame, h0, h1 = fixtures.ROD_CROSS
        cx, cy, r = c.circle
        self.assertTrue(cy-r < 2 < cy+r and cy-r > 1 and 0 < cx-r and cx+r < 4 and h0 < 0 and h1 > 3)
        # The pocketed plate's rod inside the plate, crossing the pocket's
        # walls.
        _, (c,), frame, h0, h1 = fixtures.POCKETED_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame, fixtures.POCKETED_FRAME)
        self.assertGreater(stadium_gap(6, mp.mpf(5)/2, mp.mpf(cx), mp.mpf(cy)), r)
        self.assertTrue(h0 < 0 and h1 > 2.5)
        # The chain's rod inside the stadium, its box across the box's wall
        # `x = 8`, its top `z = 2` and the stadium.
        _, (c,), frame, _, _ = fixtures.CHAIN_ROD
        cx, cy, r = c.circle
        self.assertTrue(2 <= cx <= 6 and abs(cy-3)+r < F(3, 2))
        _, (b,), frame, h0, h1 = fixtures.CHAIN_BOX
        xs = [p[0] for p in b.points]
        self.assertTrue(min(xs) < 7.5 < 8 < max(xs) and h0 < 2 < h1)

    def test_declared(self):
        # The touching ball: radius 1 about the stadium's centre line 1 above
        # its top cap (rounded once into the world).
        touch = fixtures.THROUGH_TOUCH
        self.assertEqual(touch[1], 1.0)
        u, v, w = chart(fixtures.THROUGH_STADIUM, touch[2][:3])
        self.assertLess(abs(w-4), F(1, 10**14))
        self.assertTrue(0 < u < 4 and abs(v) < F(1, 10**14))
        self.assertEqual(fixtures.BODIES['through'].second[4]+touch[1], 4.0)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'through_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('through_rod_fuse', 1), ('rod_cross_cut', 1), ('cross_pocketed_common', 2),
                           ('chain_through_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_through_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'cross_pocketed_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), ('cut', 2, 3)))
        swapped = next(c for c in listed if c.name == 'rod_cross_cut')
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
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES, t.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
