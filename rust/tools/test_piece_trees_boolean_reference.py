"""Independent checks of S9e.4b.3c.3b's fixtures
(`generate_piece_trees_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures, the chained reference on the flattened boss
against its closed form, the reference's U prisms (a box less a slot's box)
against OCCT's profiles, the frames that keep every section a circle or a
line and every seam off its body, and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_piece_trees_boolean_fixtures as fixtures
from curve_surface_reference import stored_axes


def shoelace(points):
    n = len(points)
    return sum(F(points[i][0])*F(points[(i+1) % n][1])-F(points[(i+1) % n][0])*F(points[i][1])
               for i in range(n))/2


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


class PieceTreesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # A ball beyond a plane at distance `d` from its centre by quadrature
        # of its discs, and the hemisphere's part beyond such a plane.
        for r, d in ((3, 1), (3, 2), (2.5, 0)):
            q = mp.quad(lambda s: mp.pi*(r**2-s**2), [d, r])
            self.assertLess(abs(fixtures.cap_volume(r, mp.mpf(d))-q), mp.mpf(10)**-30, (r, d))
        # The flattened boss: the half disc above the box's top face less
        # half the segment beyond the flat, by quadrature of its chords.
        r = mp.mpf(2.5)
        above = mp.quad(lambda z: mp.sqrt(r**2-z**2)+min(mp.sqrt(r**2-z**2), mp.mpf(1.5)), [0, 2, r])
        V, _ = fixtures.BODIES['dee_boss'].closed()
        self.assertLess(abs(V-(300+6*above)), mp.mpf(10)**-30)
        # The U's boss: the U's area (a box less its slot) times its height
        # and the boss's cylinder above the face.
        V, _ = fixtures.BODIES['u_boss'].closed()
        self.assertLess(abs(V-((11*8-6*2)*4+mp.pi*mp.mpf(1.5)**2*3)), mp.mpf(10)**-35)
        # The hemisphere's boss and dish on exact axes: the box less or with
        # the hemisphere's part outside it.
        frame = fixtures.at('XY', (0, 0, 0))
        for w, d in ((-1, 1), (2, -2)):
            centre = fixtures.point(frame, 4, 4, w)
            self.assertEqual(fixtures.plane_distance(frame, 0 if w < 0 else 4, centre), mp.mpf(-abs(d)))

    def test_dee_boss_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 48
        try:
            b = fixtures.BODIES['dee_boss']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
        finally:
            ref.SCAN = scan
        V0, _ = b.closed()
        self.assertLess(abs(V-V0), mp.mpf(10)**-28)

    def test_slotted_prisms(self):
        # OCCT's U prism and the reference's box less a slot's box: the same
        # profile area, the slot reaching past the open end and both caps, and
        # the U not convex.
        for name in ('u_scoop', 'u_boss', 'bites', 'u_notch', 'tooth'):
            b = fixtures.BODIES[name]
            u = next(s for s in (b.first, b.second) if isinstance(s, fixtures.Slotted))
            x0, y0, x1, y1 = u.rect
            s0, s1, t0, t1 = u.slot_rect
            profile = u.prism[1][0]
            self.assertEqual(shoelace(profile.points), F(x1-x0)*F(y1-y0)-F(s1-s0)*F(t1-t0), name)
            self.assertTrue(fixtures.reflex(profile), name)
            box, slot = u.box, u.slot
            self.assertEqual(slot[3], box[3] - 1.0)
            self.assertEqual(slot[4], box[4] + 1.0)
            self.assertGreater(slot[1][0].points[1][0], x1)
            self.assertEqual(b.ref_specs.count(u.box), 1)
            self.assertEqual(len(b.ref_specs), 3)
            self.assertEqual(sorted(leaves(b.ref_expr)), [0, 1, 2])

    def test_frames_keep_sections_circles_and_lines(self):
        # The U's groove: its ball's frame normal along the U's `y` and its
        # `x` the U's normal, exactly in the rows' integers (as the scoop's),
        # its centre on the top face within rounding.
        slab = fixtures.FRAMES['SKEW4']
        n, hint = slab[0:3], slab[3:6]
        y = (n[1]*hint[2]-n[2]*hint[1], n[2]*hint[0]-n[0]*hint[2], n[0]*hint[1]-n[1]*hint[0])
        ball = fixtures.U_SCOOP_BALL[3:]
        scale = F(y[0])/F(ball[0])
        self.assertEqual([F(c) for c in y], [scale*F(c) for c in ball[0:3]])
        self.assertEqual(ball[3:6], n)
        d = fixtures.plane_distance(fixtures.U_SCOOP_FRAME, 6, fixtures.U_SCOOP_CENTRE)
        self.assertLess(abs(d), mp.mpf(10)**-14)
        # The hemispheres on their boxes' frames: their discs and the boxes'
        # faces they meet parallels' planes; the reference's half box from
        # the disc's plane.
        for body in ('cap_boss', 'cap_pocket'):
            b = fixtures.BODIES[body]
            self.assertEqual(stored_axes(b.primitive()[2])[1:], stored_axes(b.first[2])[1:])
            self.assertEqual(b.primitive()[3:], (0.0, fixtures.HP))
            self.assertEqual(b.second.half[3], 0.0)
            self.assertEqual(b.second.half[2], b.second.ball[2])
            self.assertEqual(b.second.half[1][0].circle, (0.0, 0.0, 6.0))
        # The flat's ends on its circle exactly, the arc's frame's `x` (the
        # seam) beyond the flat; the box's top face through the axis.
        dee = fixtures.BODIES['dee_boss'].second[1][0]
        (cx, cy, r, ccw) = dee.segments[1]
        for p in dee.points:
            self.assertEqual((F(p[0])-F(cx))**2+(F(p[1])-F(cy))**2, F(r)**2)
        self.assertTrue(ccw)
        self.assertGreater(F(cx)+F(r), F(dee.points[0][0]))
        self.assertEqual(fixtures.BODIES['dee_boss'].first[4], cy)
        # The bitten cylinder's seam (its frame's `x`, at `(7, 4)`) off the
        # U; the notch's `x` reversed (its seam off the body).
        bites = fixtures.BODIES['bites']
        self.assertGreater(7.0, bites.second.rect[2])
        self.assertEqual(stored_axes(fixtures.U_NOTCH_AXIS)[1], (-1.0, 0.0, 0.0))

    def test_bodies(self):
        for b in fixtures.BODIES.values():
            self.assertTrue(b.why, b.name)
            self.assertIn(b.op, ('fuse', 'cut'))
        self.assertEqual({b.name for b in fixtures.BODIES.values() if b.kind}, {'u_notch', 'tooth'})
        self.assertEqual(fixtures.BODIES['u_notch'].kind, 'degenerate')
        self.assertEqual(fixtures.BODIES['tooth'].kind, 'unsupported')

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 30)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 24, 'empty': 0, 'degenerate': 3, 'unsupported': 3})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('u_scoop_rod_fuse', 'rod_cap_pocket_cut', 'pieces_common', 'chain_dee_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/form_'), want, name)
            self.assertEqual(nat.count('\nbrep imported/form_'), want, name)
        chain = next(c for c in listed if c.name == 'chain_dee_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'pieces_fuse')
        self.assertEqual(both.expr(), ('fuse', ('fuse', 0, ('common', 1, 2)), ('fuse', 3, 4)))
        swapped = next(c for c in listed if c.name == 'rod_cap_pocket_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, ('common', 2, 3))))
        # Every body written once: two solids about a Boolean.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of the earlier steps' bodies' files.
        import generate_imported_boolean_fixtures as a
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_polyhedra_boolean_fixtures as c
        import generate_imported_pieces_boolean_fixtures as d
        import generate_one_sphere_boolean_fixtures as e
        import generate_one_sphere_incidence_boolean_fixtures as f
        import generate_piece_forms_boolean_fixtures as g
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
