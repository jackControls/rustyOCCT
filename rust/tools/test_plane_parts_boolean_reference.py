"""Independent checks of S9e.4b.4b.2b.1's fixtures
(`generate_plane_parts_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
stack and the cake against their closed forms, every part's place (the
frustum's base inside the ball's half and its wall through the disc, the
hexagon within the stadium's footprint, the boss within the box's from
inside it, the slot's tooth crossing the plate's hole and each hole's seam
inside the other's disc), the declared case (a ball resting on the frustum's
top), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_plane_parts_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def stadium_width(length, r):
    """The stadium's width at height `y` (about its centre line)."""
    return lambda y: length+2*mp.sqrt(max(mp.mpf(0), r**2-y**2))


def polygon_width(points):
    """A convex polygon's width at height `y`: its edges' crossings."""
    def width(y):
        xs = []
        n = len(points)
        for i in range(n):
            (x0, y0), (x1, y1) = points[i], points[(i+1) % n]
            if y0 != y1 and min(y0, y1) <= y <= max(y0, y1):
                xs.append(x0+(x1-x0)*(y-y0)/(y1-y0))
        return max(xs)-min(xs) if xs else mp.mpf(0)
    return width


def chord(cx, cy, r, y):
    """A disc's chord at height `y` (an interval, or none)."""
    h = r**2-(y-cy)**2
    if h <= 0:
        return None
    s = mp.sqrt(h)
    return (cx-s, cx+s)


class PlanePartsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        half = mp.mpf(1)/2
        # The flat: the ball's half and the frustum's discs above the disc
        # (its height `w` in the frustum's chart, 1 within rounding).
        V, _ = fixtures.BODIES['flat'].closed()
        w = -fixtures.q(fixtures.height(fixtures.FLAT_FRAME, fixtures.FLAT_CONE[4][:3]))
        self.assertLess(abs(w-1), mp.mpf(10)**-15)
        ball = mp.quad(lambda z: mp.pi*(9-z**2), [-3, 0])
        frustum = mp.quad(lambda h: mp.pi*(mp.mpf(5)/4-h/7)**2, [w, mp.mpf(7)/2])
        self.assertLess(abs(V-(ball+frustum*det(fixtures.FLAT_FRAME))), tight)
        # The stack: the hexagon over 1, the stadium over 2.
        V, _ = fixtures.BODIES['stack'].closed()
        hexagon = mp.quad(polygon_width([(mp.mpf(x), mp.mpf(y)) for x, y in fixtures.HEXAGON]), [-1.5, 0, 1.5])
        stadium = mp.quad(stadium_width(6, mp.mpf(5)/2), [-mp.mpf(5)/2, 0, mp.mpf(5)/2])
        self.assertLess(abs(V-(hexagon+2*stadium)*det(fixtures.STACK_FRAME)), tight)
        # The cake: the box `9 x 7 x 2` and the boss above it, over `2 + h`
        # (`h` its frame's origin's height in the box's chart, 0 within
        # rounding).
        V, _ = fixtures.BODIES['cake'].closed()
        boss = mp.quad(stadium_width(4, mp.mpf(3)/2), [-mp.mpf(3)/2, 0, mp.mpf(3)/2])
        h = fixtures.boss_height()
        self.assertLess(abs(h), mp.mpf(10)**-15)
        self.assertLess(abs(V-(126+(2+h)*boss)*det(fixtures.CAKE_FRAME)), tight)
        # The slot: the plate less, over its height 1, the rectangle `[3/4,
        # 8] x [5/4, 8]` less both holes' discs, by rows.
        V, _ = fixtures.BODIES['slot'].closed()

        def free(y):
            spans = [s for s in (chord(4, 4, 2, y), chord(mp.mpf(11)/2, 4, mp.mpf(3)/2, y)) if s]
            covered = sum(s[1]-s[0] for s in spans)
            if len(spans) == 2:
                covered -= max(mp.mpf(0), min(s[1] for s in spans)-max(s[0] for s in spans))
            return mp.mpf(29)/4-covered
        # The circles cross at `x = 16/3`, `y = 4 +- sqrt(20) / 3`: the chords
        # part there.
        t = mp.sqrt(20)/3
        slot = mp.quad(free, [mp.mpf(5)/4, 2, mp.mpf(5)/2, 4-t, 4, 4+t, mp.mpf(11)/2, 6, 8])
        plate = mp.quad(lambda y: 8-2*mp.sqrt(max(mp.mpf(0), 4-(y-4)**2)), [0, 2, 4, 6, 8])
        self.assertLess(abs(V-(3*plate-slot)), mp.mpf(10)**-20)
        for frame in (fixtures.FLAT_FRAME, fixtures.STACK_FRAME, fixtures.CAKE_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)
        self.assertEqual(half, mp.mpf(1)/2)

    def test_bodies_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 48
        try:
            for name in ('stack', 'cake'):
                b = fixtures.BODIES[name]
                chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
                V, _, _ = chain.measures(b.ref_expr)
                V0, _ = b.closed()
                self.assertLess(abs(V-V0), mp.mpf(10)**-28, name)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        # The flat: a ball's half (its south pole to its equator) and a
        # frustum on its frame's axes about `(3/4, -1/2)`, from radius 5/4 at
        # -1 (its base inside the ball: its reach from the ball's axis, its
        # offset and radius, within the ball's radius there) to 3/4 at 5/2.
        flat = fixtures.BODIES['flat']
        self.assertEqual(flat.first.native[3:], (-fixtures.HP, 0.0))
        cone = fixtures.FLAT_CONE
        self.assertEqual(cone[4][3:], fixtures.FLAT_FRAME[3:])
        self.assertEqual(cone[4][:3], fixtures.point(fixtures.FLAT_FRAME, 0.75, -0.5, -1.0))
        offset = mp.sqrt(mp.mpf(13))/4
        self.assertLess((offset+mp.mpf(5)/4)**2, 9-1)
        self.assertEqual(cone[1:4], (1.25, 0.75, 3.5))
        # The stack: the hexagon inside the stadium's footprint (each corner
        # within a half disc or between the lines), the stadium on it.
        stack = fixtures.BODIES['stack']
        self.assertEqual(stack.first[2], stack.second[2])
        self.assertEqual((stack.first[3], stack.first[4], stack.second[3], stack.second[4]), (0.0, 1.0, 1.0, 3.0))
        for x, y in fixtures.HEXAGON:
            inside = (0 <= x <= 6 and abs(y) < 2.5) or min(x**2, (x-6)**2)+y**2 < 2.5**2
            self.assertTrue(inside, (x, y))
        # The cake: the boss within the box's footprint, from inside it.
        cake = fixtures.BODIES['cake']
        boss = cake.second
        self.assertEqual(boss[2][3:], cake.first[2][3:])
        self.assertEqual(boss[2][:3], fixtures.point(fixtures.CAKE_FRAME, 2.5, 3.5, 0.0))
        self.assertTrue(0 < 2.5-1.5 and 2.5+4+1.5 < 9 and 0 < 3.5-1.5 and 3.5+1.5 < 7)
        self.assertTrue(cake.first[3] < boss[3] < cake.first[4] < boss[4])
        # The slot: the copy half a turn about `z` at `(10, 10)`, its hole
        # about the world's `(11/2, 4)`, crossing the plate's hole; each
        # hole's seam (its frame's `+x`) inside the other's disc.
        slot = fixtures.BODIES['slot']
        o, x, y, n = stored_axes(fixtures.SLOT_TOOL)
        self.assertEqual((o, x, n), ((10.0, 10.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)))
        tool = slot.second.native
        u, v, r = tool[1][1].circle
        centre = (F(o[0])+F(u)*F(x[0])+F(v)*F(y[0]), F(o[1])+F(u)*F(x[1])+F(v)*F(y[1]))
        self.assertEqual(centre, (F(11, 2), F(4)))
        d2 = (centre[0]-4)**2+(centre[1]-4)**2
        self.assertTrue((2-F(r))**2 < d2 < (2+F(r))**2)
        seam_tool = (centre[0]-F(r), centre[1])
        self.assertLess((seam_tool[0]-4)**2+(seam_tool[1]-4)**2, 4)
        self.assertLess((6-centre[0])**2+(4-centre[1])**2, F(r)**2)
        # The copy's rectangle past the plate's walls `x = 8` and `y = 8`,
        # both discs inside its part in the plate.
        x0, x1 = 10-tool[1][0].points[1][0], 10-tool[1][0].points[0][0]
        self.assertEqual((x0, x1), (0.75, 8.75))
        self.assertTrue(x0 < 2 and 8 < x1)

    def test_declared(self):
        # The touching ball: radius 1 about the frustum's axis 1 above its top
        # disc (radius 3/4): resting on it at its centre.
        touch = fixtures.FLAT_TOUCH
        cone = fixtures.FLAT_CONE
        self.assertEqual(touch[1], 1.0)
        self.assertEqual(touch[2][:3], fixtures.point(cone[4], 0.0, 0.0, 4.5))
        self.assertEqual(cone[3]+touch[1], 4.5)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'flat_touch'},
                         {'a plane crossing a sphere within the resolution of tangency (S9d.1)'})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 24)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 21, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('flat_box_fuse', 'rod_flat_cut', 'cake_stack_common', 'chain_cake_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/part_'), want, name)
            self.assertEqual(nat.count('\nbrep imported/part_'), want, name)
        chain = next(c for c in listed if c.name == 'chain_cake_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'cake_stack_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), ('fuse', 2, 3)))
        swapped = next(c for c in listed if c.name == 'rod_flat_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('fuse', ('common', 1, 2), 3)))
        slot = next(c for c in listed if c.name == 'slot_rod_common')
        self.assertEqual(slot.expr(), ('common', ('cut', ('cut', 0, 1), ('cut', 2, 3)), 4))
        # Every body written once: two solids about a Boolean.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), list(range(len(b.ref_specs))))
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
        import generate_primitive_chains_boolean_fixtures as j
        import generate_prism_leaves_boolean_fixtures as k
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
