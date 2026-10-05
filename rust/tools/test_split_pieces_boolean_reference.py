"""Independent checks of S9e.4b.3b's fixtures
(`generate_split_pieces_boolean_fixtures.py`): the cylinder's piece's closed
form against a direct integration of its section, the chained reference on
the piece against it, the splits' classes, the half-space boxes, and the
case list and its protocol rows."""
import unittest
from fractions import Fraction as F

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_split_pieces_boolean_fixtures as fixtures


class SplitPiecesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_cylinder_cut_closed_form(self):
        # A plane normal to the axis cuts a cylinder of height h: pi r^2 h
        # and its two discs and wall. A leaning plane's piece: the volume a
        # double integral of the plane's height over the disc, its section
        # the ellipse of semi-axes r and r / cos.
        R = mp.mpf(3)
        flat = fixtures.at('XY', (5, 5, 3))
        V, A = fixtures.cylinder_cut_closed(3, (5, 5, -1), flat)
        self.assertLess(abs(V-mp.pi*R**2*4), mp.mpf(10)**-35)
        self.assertLess(abs(A-(2*mp.pi*R**2+2*mp.pi*R*4)), mp.mpf(10)**-35)
        plane = fixtures.CYL_PLANE
        V, A = fixtures.cylinder_cut_closed(3, (5, 5, -1), plane)
        o, x, y, _ = fixtures.stored_axes(plane)
        x, y = fixtures.exact(x), fixtures.exact(y)
        m = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
        p = fixtures.exact(o)

        def height(x, y):
            return p[2]-(m[0]*(x-p[0])+m[1]*(y-p[1]))/m[2]-(-1)

        def strip(x):
            w = mp.sqrt(R**2-(x-5)**2)
            # The height is linear in y: its integral the midpoint's times
            # the chord.
            return 2*w*height(x, mp.mpf(5))
        Vq = mp.quad(strip, [5-R, 5, 5+R])
        self.assertLess(abs(V-Vq), mp.mpf(10)**-30)
        # The section's area: the disc's over the normal's z-cosine.
        cos = abs(m[2])/mp.sqrt(sum(c*c for c in m))
        wall = mp.quad(lambda t: R*height(5+R*mp.cos(t), 5+R*mp.sin(t)), [0, mp.pi/2, mp.pi, 3*mp.pi/2, 2*mp.pi])
        self.assertLess(abs(A-(mp.pi*R**2+wall+mp.pi*R**2/cos)), mp.mpf(10)**-30)

    def test_cylinder_piece_on_the_chained_reference(self):
        # The cylinder below its leaning plane: the chained reference's
        # measures of `P common H` its closed form.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            piece = fixtures.PIECES['cyl_low']
            chain = ref.Chain([fixtures.construction(s, 0) for s in piece.specs()])
            V, _, A = chain.measures(('common', 0, 1))
        finally:
            ref.SCAN = scan
        V0, A0 = piece.closed()
        self.assertLess(abs(V-V0), mp.mpf(10)**-28)
        self.assertLess(abs(A-A0), mp.mpf(10)**-28)

    def test_split_classes(self):
        # Each piece is the class it claims, and the classification tells
        # them apart: the cylinder's plane normal to its axis is no oblique
        # piece, the zone's plane off its centre no half, a torus's oblique
        # plane no band.
        for p in fixtures.PIECES.values():
            self.assertEqual(fixtures.split_class(p), p.kind, p.name)
        cyl = fixtures.Piece('c', 'clipped', fixtures.CYL, fixtures.at('XY', (5, 5, 3)), 'below', 7.5)
        self.assertEqual(fixtures.split_class(cyl), 'other')
        off = fixtures.Piece('z', 'half', fixtures.ZONE, (5.0, 5.0, 5.5)+fixtures.ZONE_PLANE[3:], 'above', 5.5)
        self.assertEqual(fixtures.split_class(off), 'conic')
        tilted = fixtures.Piece('t', 'band', fixtures.TORUS, fixtures.at('TILT', (5, 5, 4.5)), 'above', 6.0)
        self.assertEqual(fixtures.split_class(tilted), 'other')

    def test_halfspace_boxes(self):
        # Each box stands on its plane (its face at offset zero through the
        # plane's origin) on the kept side, and each solid's bounding points
        # lie within four fifths of its reach from the plane's origin.
        for p in fixtures.PIECES.values():
            box = p.halfspace()
            self.assertEqual(box[2], p.plane)
            self.assertEqual(0.0 in (box[3], box[4]), True, p.name)
            self.assertEqual(box[4] > 0 if p.side == 'above' else box[3] < 0, True, p.name)
            far = fixtures.farthest(p)
            self.assertLess(far, mp.mpf(4)/5*p.reach, p.name)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 48)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 41, 'empty': 1, 'degenerate': 6, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('cyl_box_fuse', 'box_cyl_cut', 'pair_common', 'band_ball_fuse', 'chain_cyl_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'split' for s in c.items)
            self.assertEqual(enc.count('\nsplit '), want, name)
            self.assertEqual(nat.count('\nsplit '), want, name)
        # The native row: the plane's stored origin, normal and x axis.
        c = next(c for c in listed if c.name == 'cyl_box_fuse')
        row = next(r for r in c.native().split('\n') if r.startswith('split '))
        self.assertEqual(row, 'split 5.0 5.0 3.0 0.6 0.0 0.8 0.0 1.0 -0.0 below')
        chain = next(c for c in listed if c.name == 'chain_cyl_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('common', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'pair_fuse')
        self.assertEqual(both.expr(), ('fuse', ('common', 0, 1), ('common', 2, 3)))
        # A zone's piece its whole sphere common the half-space and the
        # slab between its ends' parallels.
        zone = next(c for c in listed if c.name == 'zone_box_fuse')
        self.assertEqual(zone.expr(), ('fuse', ('common', ('common', 0, 1), 2), 3))
        used = {s[1] for c in listed for s in c.items if s[0] == 'split'}
        self.assertEqual(used, set(fixtures.PIECES))


if __name__ == '__main__':
    unittest.main()
