"""Independent checks of S9e.4b.4b.2b.2b.2's fixtures
(`generate_flush_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
open pocket against its closed form, every part's place (the diamond within
the plate's heights standing out of its line walls and arc end in two parts,
the disc with a round hole within the plate's heights across its arc end with
its post inside the plate, the square with a square hole from inside the
plate past its top across its arc end, its post flush with the plate's arc),
the partners' places, the declared case (a ball resting on the slot's top)
and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_flush_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def half_width(r, u):
    """A disc's half chord at `u` from its centre (zero outside it)."""
    return mp.sqrt(max(mp.mpf(0), r**2-u**2))


def plate_half(x):
    """The plate's half width at `x` (half discs of radius 2 about `(0, 0)`
    and `(6, 0)`)."""
    if x < 0:
        return half_width(2, x)
    if x <= 6:
        return mp.mpf(2)
    return half_width(2, x-6)


def diamond_half(x):
    """The diamond's half width at `x` (its corners `(9/2, 0)`, `(13/2, +-7/2)`
    and `(15/2, 0)`)."""
    if x <= mp.mpf(13)/2:
        return max(mp.mpf(0), mp.mpf(7)/4*(x-mp.mpf(9)/2))
    return max(mp.mpf(0), mp.mpf(7)/2*(mp.mpf(15)/2-x))


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


class FlushLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        plate = mp.quad(lambda x: 2*plate_half(x), [-2, 0, 6, 8])
        # The arms: the plate over `[0, 3]` and the diamond's part outside it
        # over `[1, 2]`, its widths beyond the plate's, the crossings found
        # by root finding.
        a = mp.findroot(lambda x: diamond_half(x)-2, mp.mpf('5.6'))
        b = mp.findroot(lambda x: diamond_half(x)-plate_half(x), mp.mpf('7.0'))
        out = mp.quad(lambda x: 2*max(mp.mpf(0), diamond_half(x)-plate_half(x)),
                      [mp.mpf(9)/2, a, 6, mp.mpf(13)/2, b, mp.mpf(15)/2])
        V, _ = fixtures.BODIES['arms'].closed()
        self.assertLess(abs(V-(3*plate+out)*det(fixtures.ARMS_FRAME)), tight)
        # The slot: the plate over `[0, 3]` less, over `[1, 2]`, its part in
        # the disc of radius 3 about `(-3/2, 0)` less the post's disc.
        c = mp.findroot(lambda x: half_width(3, x+mp.mpf(3)/2)-2, mp.mpf('0.7'))
        inside = mp.quad(lambda x: 2*min(plate_half(x), half_width(3, x+mp.mpf(3)/2)),
                         [-2, 0, c, mp.mpf(3)/2])
        post = mp.quad(lambda y: 2*half_width(1, y), [-1, 0, 1])
        V, _ = fixtures.BODIES['slot'].closed()
        self.assertLess(abs(V-(3*plate-(inside-post))*det(fixtures.SLOT_FRAME)), tight)
        # The open pocket: the plate over `[0, 5/2]` less, over `[1, 5/2]`, its
        # part left of `x = 2` less the post (the hole's square `[-5/2, 1] x
        # [-1, 1]` within the plate).
        d = mp.findroot(lambda x: plate_half(x)-1, mp.mpf('-1.7'))
        left = mp.quad(lambda x: 2*plate_half(x), [-2, 0, 2])
        post = mp.quad(lambda x: 2*min(plate_half(x), mp.mpf(1)), [-2, d, 0, 1])
        V, _ = fixtures.BODIES['open'].closed()
        want = (mp.mpf(5)/2*plate-mp.mpf(3)/2*(left-post))*det(fixtures.OPEN_FRAME)
        self.assertLess(abs(V-want), tight)
        for frame in (fixtures.ARMS_FRAME, fixtures.SLOT_FRAME, fixtures.OPEN_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)
            self.assertEqual(stored_axes(frame)[3], (0.0, 0.0, 1.0))

    def test_open_on_the_chained_reference(self):
        # The open pocket alone on the chained reference (the bodies on whole
        # cases are the generator's check): the plate less the square ring
        # (its box less its hole's), its closed form `51 + 8 pi + 3 sqrt 3 /
        # 2`.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['open']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            want = (51+8*mp.pi+3*mp.sqrt(3)/2)*fixtures.det(fixtures.OPEN_FRAME)
            self.assertLess(abs(V-want), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        # The arms: the diamond on the plate's frame within its heights, its
        # corners `(9/2, 0)` and `(15/2, 0)` inside the plate by 1/2 or more,
        # its corners `(13/2, +-7/2)` outside it: two parts, each crossing the
        # plate's line wall left of its joint `(6, +-2)` and its arc right of
        # it, every crossing at least 1/4 from the joint.
        b = fixtures.BODIES['arms']
        plate, part = b.first, b.second
        self.assertEqual(plate[2], fixtures.ARMS_FRAME)
        self.assertEqual(part[2], fixtures.ARMS_FRAME)
        self.assertTrue(plate[3] < part[3] < part[4] < plate[4])
        pts = [(F(x), F(y)) for x, y in part[1][0].points]
        inside = [p for p in pts if p[1] == 0]
        outside = [p for p in pts if p[1] != 0]
        self.assertEqual(sorted(inside), [(F(9, 2), 0), (F(15, 2), 0)])
        for x, _ in inside:
            self.assertGreaterEqual(mpq(2)-mpq(max(x-6, 0)), mp.mpf(1)/2)
        for x, y in outside:
            self.assertTrue(abs(y) > 2 and 6 < x < 8)
        line = mpq(F(9, 2))+mp.mpf(8)/7
        arc = 6+(mp.mpf(147)-2*mp.sqrt(407))/106
        self.assertGreater(6-line, mp.mpf(1)/4)
        self.assertGreater(arc-6, mp.mpf(1)/4)
        # The slot: the disc of radius 3 about `(-3/2, 0)` with its hole on the
        # plate's frame within its heights; the plate's left half disc inside
        # the disc by 1/2, the disc's circle crossing the plate's line walls at
        # `x = sqrt 5 - 3/2` (in their straight part) and its right end `x =
        # 3/2` inside the plate; the hole inside the plate and the disc.
        b = fixtures.BODIES['slot']
        plate, tool = b.first, b.second
        self.assertEqual(plate[2], fixtures.SLOT_FRAME)
        self.assertEqual(tool.native[2], fixtures.SLOT_FRAME)
        outer, hole = (tuple(F(c) for c in bd.circle) for bd in tool.native[1])
        self.assertEqual(outer, (F(-3, 2), 0, 3))
        # The half disc's points `(2 cos t, 2 sin t)`, `cos t <= 0`, lie
        # `sqrt(25/4 + 6 cos t) <= 5/2` from the disc's centre.
        far = max(mp.sqrt(mpq(F(25, 4))+6*mp.cos(t)) for t in mp.linspace(mp.pi/2, 3*mp.pi/2, 33))
        self.assertLessEqual(far, mpq(F(5, 2))+mp.mpf(10)**-30)
        self.assertTrue(0 < mp.sqrt(5)-mp.mpf(3)/2 < 6 and outer[0]+outer[2] < 6)
        hx, hy, hr = hole
        self.assertLess(mp.sqrt(mpq(hx**2+hy**2))+mpq(hr), 2-mp.mpf(1)/4)
        self.assertLess(mp.sqrt(mpq((hx-outer[0])**2+hy**2))+mpq(hr), mpq(outer[2])-mp.mpf(1)/2)
        self.assertTrue(plate[3] < tool.native[3] < tool.native[4] < plate[4])
        # The open pocket: the square ring from inside the plate past its top,
        # its square across the plate's arc end (its wall `x = 2` in the
        # plate's straight part), its hole's walls `y = +-1` crossing the arc
        # at `x = -sqrt 3`, its wall `x = 1` inside the plate.
        b = fixtures.BODIES['open']
        plate, ring = b.first, b.second
        self.assertEqual(plate[2], fixtures.OPEN_FRAME)
        self.assertEqual(ring.native[2][3:], fixtures.OPEN_FRAME[3:])
        sq, hole = ring.native[1]
        xs, ys = [F(p[0]) for p in sq.points], [F(p[1]) for p in sq.points]
        self.assertEqual((min(xs), max(xs), min(ys), max(ys)), (-4, 2, -3, 3))
        xs, ys = [F(p[0]) for p in hole.points], [F(p[1]) for p in hole.points]
        self.assertEqual((min(xs), max(xs), min(ys), max(ys)), (F(-5, 2), 1, -1, 1))
        self.assertTrue(-mpq(F(5, 2)) < -mp.sqrt(3) < 0)
        self.assertTrue(plate[3] < ring.native[3] < plate[4] < ring.native[4])

    def test_partners_placed(self):
        # The arms' rod across the upper part's wall (the diamond's edge from
        # `(9/2, 0)` to `(13/2, 7/2)`), the plate's line wall `y = 2` and its
        # arc, clear of the diamond's other edge and its top corner.
        _, (c,), frame, h0, h1 = fixtures.ARMS_ROD
        cx, cy, rad = (mpq(v) for v in c.circle)
        self.assertEqual(frame, fixtures.ARMS_FRAME)
        edge = (mpq(F(7, 2))*(cx-mpq(F(9, 2)))-2*cy)/mp.sqrt(mpq(F(65, 4)))
        self.assertLess(abs(edge), rad)
        self.assertLess(abs(cy-2), rad)
        self.assertTrue(2-rad < mp.sqrt((cx-6)**2+cy**2) < 2+rad)
        other = abs(mpq(F(7, 2))*(cx-mpq(F(15, 2)))+cy)/mp.sqrt(mpq(F(53, 4)))
        self.assertGreater(other, rad)
        self.assertTrue(h0 < 0 and h1 > 3)
        # The slot's rod across the pocket's wall (the circle of radius 3 about
        # `(-3/2, 0)`), clear of the post and the plate's walls.
        _, (c,), frame, h0, h1 = fixtures.SLOT_ROD
        cx, cy, rad = (mpq(v) for v in c.circle)
        self.assertEqual(frame, fixtures.SLOT_FRAME)
        self.assertTrue(3-rad < mp.sqrt((cx+mpq(F(3, 2)))**2+cy**2) < 3+rad)
        self.assertGreater(mp.sqrt((cx+mpq(F(1, 2)))**2+(cy-mpq(F(1, 4)))**2), 1+rad)
        self.assertLess(abs(cy)+rad, 2)
        # The tool's rod across the post's circle and the plate's arc, inside
        # the disc, clear of the plate's line walls' carriers.
        _, (c,), frame, h0, h1 = fixtures.ROD_SLOT
        cx, cy, rad = (mpq(v) for v in c.circle)
        d = mp.sqrt((cx+mpq(F(1, 2)))**2+(cy-mpq(F(1, 4)))**2)
        self.assertTrue(abs(1-rad) < d < 1+rad)
        self.assertTrue(2-rad < mp.sqrt(cx**2+cy**2) < 2+rad)
        self.assertLess(mp.sqrt((cx+mpq(F(3, 2)))**2+cy**2)+rad, 3)
        self.assertLess(abs(cy)+rad, 2)
        # The open pocket's rod across its wall `x = 2`, clear of the post and
        # the plate's line walls.
        _, (c,), frame, h0, h1 = fixtures.OPEN_ROD
        cx, cy, rad = (mpq(v) for v in c.circle)
        self.assertEqual(frame, fixtures.OPEN_FRAME)
        self.assertTrue(cx-rad < 2 < cx+rad and abs(cy)+rad < 2 and cx-rad > 1)
        # The chain's rod inside the plate's straight part, its box across the
        # upper part and the plate's top.
        _, (c,), frame, _, _ = fixtures.CHAIN_ROD
        cx, cy, rad = (mpq(v) for v in c.circle)
        self.assertTrue(cx-rad > 0 and cx+rad < 6 and abs(cy)+rad < 2)
        _, (b,), frame, h0, h1 = fixtures.CHAIN_BOX
        xs = [F(p[0]) for p in b.points]
        ys = [F(p[1]) for p in b.points]
        self.assertTrue(min(xs) < 6 < max(xs) < 7 and min(ys) < 2 < max(ys) and h0 < 2 and 3 < h1)

    def test_declared(self):
        # The touching ball: radius 1 about the slot's chart point `(4, 0)` 1
        # above its top cap (its frame's normal the world's `z` exactly: the
        # centre's height exact), over the plate's material.
        touch = fixtures.SLOT_TOUCH
        self.assertEqual(touch[1], 1.0)
        o = stored_axes(fixtures.SLOT_FRAME)[0]
        self.assertEqual(touch[2][2], o[2]+4.0)
        self.assertEqual(fixtures.BODIES['slot'].first[4]+touch[1], 4.0)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'slot_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('arms_rod_fuse', 1), ('rod_slot_cut', 1), ('arms_slot_common', 2),
                           ('chain_arms_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_arms_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'arms_slot_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), ('cut', 2, ('cut', 3, 4))))
        swapped = next(c for c in listed if c.name == 'rod_slot_cut')
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
        import generate_rim_leaves_boolean_fixtures as r
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES, t.BODIES, w.BODIES, r.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
