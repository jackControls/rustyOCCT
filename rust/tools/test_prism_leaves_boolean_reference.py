"""Independent checks of S9e.4b.4b.2a's fixtures
(`generate_prism_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
dimple and the pocket against their closed forms, the prisms' profiles (the
rounded square's and the stadium's arcs tangent to their lines), every
primitive's place (the bore across the flats of two walls, the boss and the
pocket along their plates' axes, the dimple's and the dome's balls meeting
the top cap alone), the declared bodies and case (the post on a fillet's
circle, the notch through both caps at the rim, a rod tangent to the
bore), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_prism_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def rounded_width(side, r):
    """The rounded square's width at height `y` in its chart."""
    def width(y):
        d = min(y, side-y)
        if d >= r:
            return mp.mpf(side)
        return side-2*(r-mp.sqrt(r**2-(r-d)**2))
    return width


def stadium_width(length, r):
    """The stadium's width at height `y` (about its centre line)."""
    return lambda y: length+2*mp.sqrt(max(mp.mpf(0), r**2-y**2))


def disc_quad(radius, lo, hi, breaks=()):
    return mp.quad(lambda z: mp.pi*radius(z)**2, [lo, *breaks, hi])


class PrismLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-30
        square = mp.quad(rounded_width(9, 2), [0, 2, 7, 9])
        stadium = mp.quad(stadium_width(6, mp.mpf(5)/2), [-mp.mpf(5)/2, 0, mp.mpf(5)/2])
        det = lambda frame: fixtures.det(frame)
        # The bore: the square over 6 less the disc of radius 3/2 over 9 (its
        # part between the walls on the stored axes).
        V, _ = fixtures.BODIES['bore'].closed()
        bore = disc_quad(lambda w: mp.mpf(3)/2, 0, 9)
        self.assertLess(abs(fixtures.bore_volume()-bore), mp.mpf(10)**-14)
        self.assertLess(abs(V-(6*square*det(fixtures.BORE_FRAME)-fixtures.bore_volume())), tight)
        # The boss: the stadium over 2 and the boss's disc over `[2, 5]`.
        V, _ = fixtures.BODIES['boss'].closed()
        self.assertLess(abs(V-(2*stadium+disc_quad(lambda w: mp.mpf(5)/4, 2, 5))*det(fixtures.BOSS_FRAME)), tight)
        # The dimple: the ball of radius 5/2 about height 9/2 below the cap 3.
        V, _ = fixtures.BODIES['dimple'].closed()
        ball = lambda c, r: (lambda z: mp.sqrt(max(mp.mpf(0), r**2-(z-c)**2)))
        dimple = disc_quad(ball(mp.mpf(9)/2, mp.mpf(5)/2), 2, 3)
        self.assertLess(abs(V-(3*square-dimple)), tight)
        # The pocket: the frustum's radius `1/2 + (w - b) / 2` over `[b, 3]`,
        # `b` its base's height in the plate's chart.
        V, _ = fixtures.BODIES['pocket'].closed()
        b = fixtures.q(fixtures.pocket_depth())
        self.assertLess(abs(b-1), mp.mpf(10)**-15)
        pocket = disc_quad(lambda w: mp.mpf(1)/2+(w-b)/2, b, 3)
        self.assertLess(abs(V-(3*stadium-pocket)*det(fixtures.POCKET_FRAME)), tight)
        # The dome: the ball of radius 2 beyond the cap at its world distance
        # `g` from the centre (its chart height 5/2 within rounding).
        V, _ = fixtures.BODIES['dome'].closed()
        c = fixtures.q(fixtures.dome_centre())
        self.assertLess(abs(c-mp.mpf(5)/2), mp.mpf(10)**-15)
        g = fixtures.dome_gap()
        self.assertLess(abs(g-mp.mpf(1)/2), mp.mpf(10)**-15)
        dome = disc_quad(ball(0, 2), g, 2)
        self.assertLess(abs(V-(3*square*det(fixtures.DOME_FRAME)+dome)), tight)
        for frame in (fixtures.BORE_FRAME, fixtures.BOSS_FRAME, fixtures.POCKET_FRAME, fixtures.DOME_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)

    def test_bodies_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 48
        try:
            for name in ('dimple', 'pocket'):
                b = fixtures.BODIES[name]
                chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
                V, _, _ = chain.measures(b.ref_expr)
                V0, _ = b.closed()
                self.assertLess(abs(V-V0), mp.mpf(10)**-28, name)
        finally:
            ref.SCAN = scan

    def test_profiles(self):
        # Every arc tangent to the lines it joins: its ends a radius from
        # its centre, and each joining line along the radius's normal there.
        for b in fixtures.BODIES.values():
            boundary = b.first[1][0]
            pts, segs = boundary.points, boundary.segments
            for k, seg in enumerate(segs):
                if seg is None:
                    continue
                cx, cy, r, ccw = seg
                self.assertTrue(ccw, b.name)
                for j in (k, (k+1) % len(pts)):
                    p = pts[j]
                    radius = (F(p[0])-F(cx), F(p[1])-F(cy))
                    self.assertEqual(radius[0]**2+radius[1]**2, F(r)**2, b.name)
                    # The line on the other side of this end.
                    i = j-1 if j == k else j
                    a, c = pts[i % len(pts)], pts[(i+1) % len(pts)]
                    self.assertIsNone(segs[i % len(pts)], b.name)
                    line = (F(c[0])-F(a[0]), F(c[1])-F(a[1]))
                    self.assertEqual(line[0]*radius[0]+line[1]*radius[1], 0, b.name)

    def test_primitives_placed(self):
        # The bore along `y` (`ALONGY`: its profile's `(3, 9/2)` is `z = 3`,
        # `x = 9/2`) at the plate's origin: within the flats `[2, 7]` of the
        # walls it crosses, clear of both caps, through both walls.
        self.assertEqual(fixtures.FRAMES['ALONGY'], (0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
        self.assertEqual(fixtures.BORE_HOLE[:3], fixtures.BORE_FRAME[:3])
        self.assertEqual(fixtures.BORE_FRAME[3:], fixtures.FRAMES['XY'])
        hole = fixtures.BODIES['bore'].second
        cx, cy, r = hole[1][0].circle
        self.assertEqual((cx, cy, r), (3.0, 4.5, 1.5))
        self.assertTrue(2 <= cy-r and cy+r <= 7 and 0 < cx-r and cx+r < 6)
        self.assertTrue(hole[3] < 0 and hole[4] > 9)
        # The boss and the pocket on their plates' frames' axes.
        boss = fixtures.BODIES['boss']
        self.assertEqual(boss.second[2][3:], boss.first[2][3:])
        self.assertEqual(boss.second[2][:3], boss.first[2][:3])
        self.assertTrue(boss.second[3] < 2 < boss.second[4])
        pocket = fixtures.BODIES['pocket']
        self.assertEqual(pocket.second[4][3:], pocket.first[2][3:])
        # The dimple's and the dome's balls: their sections by the top cap
        # parallels, clear of the bottom cap and of the walls and fillets.
        half = F(9, 2)
        for name, centre, r, top in (('dimple', (half, half, F(9, 2)), F(5, 2), 3), ('dome', (half, half, F(5, 2)), 2, 3)):
            b = fixtures.BODIES[name]
            self.assertEqual(b.second[2][3:], b.first[2][3:], name)
            self.assertLess(abs(centre[2]-top), r, name)
            self.assertGreater(centre[2]-r, 0, name)
            # Its reach from the axis within the plate (below the top cap)
            # inside the rounded square's flats `[2, 7]` both ways.
            reach2 = r**2-max(0, centre[2]-top)**2
            for c in centre[:2]:
                self.assertTrue((c-2)**2 >= reach2 and (7-c)**2 >= reach2, name)
        self.assertEqual(fixtures.BODIES['dome'].second[2][:3],
                         fixtures.point(fixtures.DOME_FRAME, 4.5, 4.5, 2.5))

    def test_declared(self):
        # The post on the corner's fillet's circle: about its centre `(2, 2)`
        # with its radius, from inside the plate to above its top.
        post = fixtures.BODIES['post']
        cx, cy, r = post.second[1][0].circle
        self.assertEqual((cx, cy, r), (2.0, 2.0, 2.0))
        arcs = [seg for seg in post.first[1][0].segments if seg is not None]
        self.assertIn((2.0, 2.0, 2.0, True), arcs)
        self.assertEqual(post.second[2], post.first[2])
        self.assertTrue(post.first[3] < post.second[3] < post.first[4] < post.second[4])
        # The notch's rod through both caps (`[0, 2]`), ending inside.
        notch = fixtures.BODIES['notch']
        cx, cy, r = notch.second[1][0].circle
        self.assertTrue(cx-r < notch.first[3] and notch.first[4] < cx+r)
        self.assertTrue(notch.second[3] < 0 < notch.second[4] < 9)
        self.assertEqual({b.name: b.kind for b in fixtures.BODIES.values() if b.kind},
                         {'post': 'degenerate', 'notch': 'unsupported'})
        # The touching rod: radius 1 along `x` (`SIDE`: its profile's `(9/2,
        # 11/2)` is `y = 9/2`, `z = 11/2`) above the bore's axis `(9/2, *, 3)`
        # by the radii's sum, the axes crossing at right angles: tangent at
        # one point.
        self.assertEqual(fixtures.FRAMES['SIDE'], (1.0, 0.0, 0.0, 0.0, 1.0, 0.0))
        touch = fixtures.BORE_TOUCH
        self.assertEqual(touch[2][:3], fixtures.BORE_FRAME[:3])
        u, v, rr = touch[1][0].circle
        bore = fixtures.BODIES['bore'].second[1][0].circle
        self.assertEqual((u, v-rr-bore[2]), (bore[1], bore[0]))
        self.assertTrue(touch[3] < 0 and touch[4] > 9)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 24)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 21, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('bore_rod_fuse', 'dome_box_cut', 'boss_pocket_common', 'chain_dimple_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/leaf_'), want, name)
            self.assertEqual(nat.count('\nbrep imported/leaf_'), want, name)
        chain = next(c for c in listed if c.name == 'chain_dimple_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'boss_pocket_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), ('cut', 2, 3)))
        swapped = next(c for c in listed if c.name == 'dome_box_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('fuse', 1, 2)))
        touch = [c for c in listed if c.group == 'bore_touch']
        self.assertEqual({c.reason for c in touch}, {'a tangency between the inputs (S9c)'})
        # Every body written once: a prism and a primitive about a Boolean.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        # The post and the notch given to no case (the reference's cost),
        # refused on import by the kernel's tests.
        self.assertEqual(used, set(fixtures.BODIES)-{'post', 'notch'})
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), [0, 1])
            self.assertEqual(b.first[0], 'prism', b.name)
            self.assertGreaterEqual(len(b.curved_kinds()), 3, b.name)
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
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
