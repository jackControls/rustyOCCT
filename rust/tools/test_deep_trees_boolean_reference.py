"""Independent checks of S9e.4b.4c.2a's fixtures
(`generate_deep_trees_boolean_fixtures.py`): the bodies' closed forms against
direct quadratures of their sections, the chained reference on the holed
prisms against their areas, the reference's sets against OCCT's
profiles (a square with a square or a U-shaped hole), every part's place
(each pocket inside its primitive from below its top cap past it, the post
and the U island inside their pockets, the U's notch open, the pockets'
depths), the partners' places, the declared case (a ball resting on the
post's top), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_deep_trees_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def area(points):
    """A polygon's signed area by the shoelace formula, exactly."""
    p = [(F(x), F(y)) for x, y in points]
    n = len(p)
    return sum(p[i][0]*p[(i+1) % n][1]-p[(i+1) % n][0]*p[i][1] for i in range(n))/2


def reflex(points):
    """Whether a polygon has a reflex corner (counter-clockwise)."""
    p = points
    n = len(p)
    turns = [(p[(i+1) % n][0]-p[i][0])*(p[(i+2) % n][1]-p[(i+1) % n][1])
             - (p[(i+1) % n][1]-p[i][1])*(p[(i+2) % n][0]-p[(i+1) % n][0]) for i in range(n)]
    return any(t < 0 for t in turns)


def chart(frame, point):
    """A world point's coordinates in a frame's stored axes (exactly
    orthonormal to rounding: within 1e-15)."""
    o, x, y, n = stored_axes(frame)
    d = [F(point[i])-F(o[i]) for i in range(3)]
    return tuple(float(sum(d[i]*F(a[i]) for i in range(3))) for a in (x, y, n))


class DeepTreesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        # The post: the cylinder's discs less the ring `25 - 4` from 3 to 6.
        V, _ = fixtures.BODIES['post'].closed()
        post = mp.quad(lambda z: 16*mp.pi-(21 if z > 3 else 0), [0, 3, 6])
        self.assertLess(abs(V-post*det(fixtures.POST_FRAME)), tight)
        # The well: the cylinder's discs less the square `25` less the U
        # island `9 - 2` from 5/2 to 13/2.
        V, _ = fixtures.BODIES['well'].closed()
        well = mp.quad(lambda z: 16*mp.pi-(18 if z > mp.mpf(5)/2 else 0), [0, mp.mpf(5)/2, mp.mpf(13)/2])
        self.assertLess(abs(V-well*det(fixtures.WELL_FRAME)), tight)
        for frame in (fixtures.POST_FRAME, fixtures.WELL_FRAME):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)

    def test_holed_prisms_on_the_chained_reference(self):
        # The reference's holed prisms alone (the bodies' closed forms on the
        # whole chain are the generator's check): the post's ring `21` over
        # 4, the well's square less its U island `25 - 7` over 5, each in
        # its chart times the stored axes' determinant.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            for name, want in (('post', 84), ('well', 90)):
                h = fixtures.BODIES[name].holed
                chain = ref.Chain([fixtures.construction(s, 0) for s in h.parts])
                V, _, _ = chain.measures(h.expr)
                V0 = want*fixtures.det(h.native[2])
                self.assertLess(abs(V-V0), mp.mpf(10)**-28, name)
        finally:
            ref.SCAN = scan

    def test_reference_sets_are_the_profiles(self):
        # Each holed prism: OCCT's square and hole, the reference's boxes of
        # the same planes (the hole's past both caps, the notch's past the U's
        # open side), and the same area in the chart.
        for name, want in (('post', 21), ('well', 18)):
            h = fixtures.BODIES[name].holed
            outer, inner = h.native[1]
            boxes = [p[1][0].points for p in h.parts]
            self.assertEqual(outer.points, boxes[0])
            self.assertEqual(area(outer.points)-area(inner.points), want, name)
            for p in h.parts[1:]:
                self.assertLess(p[3], h.native[3])
                self.assertGreater(p[4], h.native[4])
            self.assertEqual({p[2] for p in h.parts}, {h.native[2]})
        well = fixtures.BODIES['well'].holed
        u, notch = well.native[1][1].points, well.parts[2][1][0].points
        self.assertTrue(reflex(u))
        self.assertEqual(area(u), area(well.parts[1][1][0].points)-area([(-0.5, -0.5), (1.5, -0.5), (1.5, 0.5),
                                                                         (-0.5, 0.5)]))
        # The notch's box: the notch's planes, past the U's open side.
        self.assertEqual((notch[0], notch[2][1]), ((-0.5, -0.5), 0.5))
        self.assertGreater(notch[1][0], max(x for x, _ in u))
        self.assertEqual(well.expr, ('cut', 0, ('cut', 1, 2)))
        post = fixtures.BODIES['post'].holed
        self.assertFalse(reflex(post.native[1][1].points))
        self.assertEqual(post.expr, ('cut', 0, 1))

    def test_parts_placed(self):
        # Each pocket from inside its primitive past its top cap, its square's
        # corners inside the primitive's wall at every height of the pocket.
        post, well = fixtures.BODIES['post'], fixtures.BODIES['well']
        cyl = post.first
        self.assertEqual(cyl[2], fixtures.POST_FRAME)
        self.assertEqual((cyl[3], cyl[4], cyl[1][0].circle), (0.0, 6.0, (0.0, 0.0, 4.0)))
        self.assertTrue(cyl[3] < post.holed.native[3] < cyl[4] < post.holed.native[4])
        self.assertLess(2*F(5, 2)**2, F(4)**2)
        cyl = well.first
        self.assertEqual(cyl[2], fixtures.WELL_FRAME)
        self.assertEqual((cyl[3], cyl[4], cyl[1][0].circle), (0.0, 6.5, (0.0, 0.0, 4.0)))
        self.assertTrue(cyl[3] < well.holed.native[3] < cyl[4] < well.holed.native[4])
        # The post and the U island inside their squares, apart from them.
        for b in (post, well):
            outer, inner = b.holed.native[1]
            ox = [x for x, _ in outer.points]
            ix = [x for x, _ in inner.points]
            iy = [y for _, y in inner.points]
            self.assertTrue(min(ox) < min(ix) and max(ix) < max(ox) and min(ox) < min(iy) and max(iy) < max(ox))
        # Pockets nest: the tooth's and the post's two deep, the well's three.
        self.assertEqual({b.name: getattr(b, 'depth', 2) for b in fixtures.BODIES.values()},
                         {'post': 2, 'well': 3, 'tooth': 2})
        tooth = fixtures.BODIES['tooth']
        self.assertEqual(tooth.path, 'imported/form_tooth.brep')
        self.assertEqual(tooth.ref_expr, ('cut', 0, ('cut', 1, 2)))

    def test_partners_placed(self):
        # The pin crosses the tooth's wall `x = 2` within its face.
        _, (c,), frame, h0, h1 = fixtures.TOOTH_PIN
        cx, cy, r = c.circle
        self.assertTrue(cx-r < 2 < cx+r < 3 and 3 < cy-r and cy+r < 4.5)
        self.assertTrue(h0 < 0 and h1 > 6)
        # The post's rod crosses the post's wall `x = 1` within its face, the
        # pocket's floor and both caps.
        _, (c,), frame, h0, h1 = fixtures.POST_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame, fixtures.POST_FRAME)
        self.assertTrue(cx-r < 1 < cx+r < 2.5 and -1 < cy-r and cy+r < 1 and h0 < 0 and h1 > 6)
        # The well's rod crosses the notch's side `y = 1/2` inside the notch.
        _, (c,), frame, h0, h1 = fixtures.WELL_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame, fixtures.WELL_FRAME)
        self.assertTrue(cy-r < 0.5 < cy+r < 1.5 and -0.5 < cx-r and cx+r < 1.5 and h0 < 0 and h1 > 6.5)
        # The tool's rod crosses the post's wall `x = -1` within its face.
        _, (c,), frame, h0, h1 = fixtures.ROD_POST
        cx, cy, r = c.circle
        self.assertEqual(frame, fixtures.POST_FRAME)
        self.assertTrue(-2.5 < cx-r < -1 < cx+r and -1 < cy-r and cy+r < 1 and h0 < 0 and h1 > 6)
        # The chain's first rod inside the post, its second across the
        # pocket's wall `x = 5/2` from below the floor.
        _, (c,), frame, _, _ = fixtures.CHAIN_ROD
        cx, cy, r = c.circle
        self.assertTrue(-1 < cx-r and cx+r < 1 and -1 < cy-r and cy+r < 1)
        _, (c,), frame, h0, h1 = fixtures.CHAIN_PIN
        cx, cy, r = c.circle
        self.assertTrue(1 < cx-r < 2.5 < cx+r and -2.5 < cy-r and h0 < 3 and h1 > 6)
        # S9e.4a's cylinder across the post: its axis's chart and its distance
        # from each plane of the post along it (0.2 at least from tangency,
        # crossing the planes of the post's walls and of the pocket's walls
        # `x = 5/2` and `y = 5/2`, missing the others), the post's cylinder
        # crossing it, its top above the cap's plane.
        spec = fixtures.base.BODIES['cyl'].spec
        _, (c,), frame, h0, h1 = spec
        self.assertEqual(c.circle, (0.0, 0.0, 3.0))
        u, v, w = chart(fixtures.POST_FRAME, frame[:3])
        crossed = 0
        for d in (u-2.5, u+2.5, u-1, u+1, v-2.5, v+2.5, v-1, v+1):
            self.assertGreater(abs(abs(d)-3), 0.2)
            crossed += abs(d) < 3
        self.assertEqual(crossed, 6)
        rho = mp.sqrt(u*u+v*v)
        self.assertTrue(4-3 < rho < 4+3)
        self.assertEqual((w+h0, w+h1), (0.25, 6.25))

    def test_declared(self):
        # The touching ball: radius 1 about the post's axis 1 above its top,
        # exactly (the post's frame's normal the world's `z`, its origin and
        # the height dyadic).
        touch = fixtures.POST_TOUCH
        self.assertEqual(touch[1], 1.0)
        o, _, _, n = stored_axes(fixtures.POST_FRAME)
        self.assertEqual(n, (0.0, 0.0, 1.0))
        self.assertEqual(touch[2][:3], (o[0], o[1], o[2]+7.0))
        self.assertEqual(fixtures.BODIES['post'].first[4]+touch[1], 7.0)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'post_touch'},
                         {'a tangency between the inputs (S9c)'})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('post_rod_fuse', 1), ('rod_post_cut', 1), ('post_cyl_common', 2), ('chain_post_cut', 1),
                           ('tooth_pin_cut', 1)):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            self.assertEqual(enc.count('\nbrep imported/'), want, name)
            self.assertEqual(nat.count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_post_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('cut', 0, ('cut', 1, 2)), 3), 4))
        both = next(c for c in listed if c.name == 'post_cyl_cut')
        self.assertEqual(both.expr(), ('cut', ('cut', 0, ('cut', 1, 2)), 3))
        swapped = next(c for c in listed if c.name == 'rod_post_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('cut', 1, ('cut', 2, 3))))
        well = next(c for c in listed if c.name == 'well_rod_cut')
        self.assertEqual(well.expr(), ('cut', ('cut', 0, ('cut', 1, ('cut', 2, 3))), 4))
        # This step's bodies written once each, a primitive less a prism.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.OWN))
        self.assertEqual(text.count('\nboolean cut'), len(fixtures.OWN))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), list(range(len(b.ref_specs))))
        # None of the earlier steps' bodies' files but the tooth's.
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
        paths = {b_.path for b_ in fixtures.OWN.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES, i.BODIES,
                        j.BODIES, k.BODIES, m.BODIES, p.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
