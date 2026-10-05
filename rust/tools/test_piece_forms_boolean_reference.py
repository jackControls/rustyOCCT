"""Independent checks of S9e.4b.3c.3a's fixtures
(`generate_piece_forms_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures, the chained reference on the slot against its
closed form, the frames that keep every section a circle or a line, the test
that a body is no S9e.4b.3a piece, and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_piece_forms_boolean_fixtures as fixtures
from curve_surface_reference import stored_axes


class PieceFormsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # A disc's segment beyond a chord by quadrature of its chords' lengths
        # (`2 sqrt(r^2 - s^2)` for `s` from `d` to `r`).
        for r, d in ((1.25, 0.5), (3, 1), (2, 0)):
            q = mp.quad(lambda s: 2*mp.sqrt(r**2-s**2), [d, r])
            self.assertLess(abs(fixtures.segment_area(r, d)-q), mp.mpf(10)**-30, (r, d))
        # A ball's cap beyond a plane by quadrature of its discs.
        for r, d in ((2.5, 1.5), (2.5, -1), (3.75, 0)):
            q = mp.quad(lambda s: mp.pi*(r**2-s**2), [d, r])
            self.assertLess(abs(fixtures.cap_volume(r, mp.mpf(d))-q), mp.mpf(10)**-30, (r, d))
        # G4's part: the box and the boss's cylinder, the boss's disc moved
        # from the box's face to the cap.
        V, A = fixtures.BODIES['boss'].closed()
        self.assertLess(abs(V-(10*10*8+mp.pi*2.5**2*4)), mp.mpf(10)**-35)
        self.assertLess(abs(A-(2*(100+80+80)+2*mp.pi*2.5*4)), mp.mpf(10)**-35)
        # A plane's distance on exact axes: the dimple's top face 3/2 under
        # its ball's centre.
        frame = fixtures.at('XY', (1, 2, 0.5))
        centre = fixtures.point(frame, 4, 4, 5.5)
        self.assertEqual(fixtures.plane_distance(frame, 4, centre), mp.mpf(1.5))

    def test_slot_on_the_chained_reference(self):
        # The slot (a box less a cylinder ending inside it): the chained
        # reference's volume of `B cut P` its closed form.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            b = fixtures.BODIES['slot']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.specs()])
            V, _, _ = chain.measures(('cut', 0, 1))
        finally:
            ref.SCAN = scan
        V0, _ = b.closed()
        self.assertLess(abs(V-V0), mp.mpf(10)**-28)

    def test_frames_keep_sections_circles_and_lines(self):
        # The scoop's ball: its frame's normal along the slab's `y` and its
        # `x` the slab's normal, exactly in the rows' integers, so the top
        # face is a meridian plane of the ball's frame and the side faces
        # its parallels' planes.
        slab = fixtures.FRAMES['SKEW4']
        n, hint = slab[0:3], slab[3:6]
        y = (n[1]*hint[2]-n[2]*hint[1], n[2]*hint[0]-n[0]*hint[2], n[0]*hint[1]-n[1]*hint[0])
        ball = fixtures.SCOOP_BALL[3:]
        scale = F(y[0])/F(ball[0])
        self.assertEqual([F(c) for c in y], [scale*F(c) for c in ball[0:3]])
        self.assertEqual(ball[3:6], n)
        # The ball's centre on the top face's plane within rounding.
        d = fixtures.plane_distance(fixtures.SCOOP_FRAME, 6, fixtures.SCOOP_CENTRE)
        self.assertLess(abs(d), mp.mpf(10)**-14)
        # The three-quarter frustum's seam between the removed quarter's
        # planes; the notch's seam off the body (its `x` against the world's).
        x = stored_axes(fixtures.QUARTER_CONE)[1]
        self.assertGreater(x[0], 0)
        self.assertGreater(x[1], 0)
        self.assertEqual(stored_axes(fixtures.NOTCH_AXIS)[1], (-1.0, 0.0, 0.0))
        # G4's boss and its tool along the world's `y`, their caps on planes
        # of exact coordinates.
        self.assertEqual(stored_axes(fixtures.BOSS_FRAME)[3], (0.0, 1.0, 0.0))

    def test_bodies_are_this_steps(self):
        for b in fixtures.BODIES.values():
            self.assertTrue(fixtures.this_steps(b), b.name)
        self.assertEqual({b.form for b in fixtures.BODIES.values()}, {'groove', 'bite', 'boss', None})
        # A primitive common a box is S9e.4b.3a's.
        common = fixtures.Body('c', fixtures.sphere(1.0, fixtures.at('XY', (0, 0, 0))), 'common',
                               fixtures.block(fixtures.at('XY', (0, 0, 0)), 0.0, 0.0, 2.0, 2.0, 0.0, 2.0), None)
        self.assertFalse(fixtures.this_steps(common))
        for name in ('notch', 'quarter'):
            self.assertIsNone(fixtures.BODIES[name].form)
            self.assertTrue(fixtures.BODIES[name].reason)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 39)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 33, 'empty': 0, 'degenerate': 6})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('boss_rod_fuse', 'slab_scoop_cut', 'pieces_common', 'chain_slot_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/form_'), want, name)
            self.assertEqual(nat.count('\nbrep imported/form_'), want, name)
        chain = next(c for c in listed if c.name == 'chain_slot_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'pieces_fuse')
        self.assertEqual(both.expr(), ('fuse', ('fuse', 0, 1), ('cut', 2, 3)))
        swapped = next(c for c in listed if c.name == 'slab_scoop_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, 2)))
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
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
