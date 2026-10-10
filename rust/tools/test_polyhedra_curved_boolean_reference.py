"""Independent checks of S9e.4b.4c.1's reference and fixtures
(`polyhedra_curved_boolean_reference.py`,
`generate_polyhedra_curved_boolean_fixtures.py`): the hull input (its faces,
closed forms against exact determinants and S9e.4b.2's constructions, its
chords' sweep against its closed form, its binary64 rays), each body the
construction S9e.4b.2 or S9e.4a gave OCCT, the cavity's inner box inside its
outer one and the slab between its faces, the partners' places (the hollow
box's ball through its wall into its cavity and tangent to no plane of it,
the inner ball inside the cavity, the touching ball on the tetrahedron's
base), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as earlier_curved
import generate_imported_polyhedra_boolean_fixtures as earlier
import generate_polyhedra_curved_boolean_fixtures as fixtures
import imported_polyhedra_boolean_reference as ipr
import polyhedra_curved_boolean_reference as ref
import polyhedral_reference as pr


def local(frame, p):
    """A point's coordinates on a frame's stored axes, exactly."""
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    d = tuple(F(p[i])-o[i] for i in range(3))
    det = pr.det3(x, y, n)
    return (pr.det3(d, y, n)/det, pr.det3(x, d, n)/det, pr.det3(x, y, d)/det)


def corners(frame, size):
    """A box's corners (`MakeBox` on the frame's stored axes), exactly."""
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    return [tuple(o[i]+F(u)*x[i]+F(v)*y[i]+F(w)*n[i] for i in range(3))
            for u in (0, size[0]) for v in (0, size[1]) for w in (0, size[2])]


class PolyhedraCurvedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_hull(self):
        # A unit cube's hull: six faces of four points, volume one, its
        # centre, area six.
        cube = ref.Hull(ref.HullCase([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]))
        self.assertEqual(sorted(len(c) for c in cube.cycles), [4]*6)
        V, mom, A = cube.closed()
        self.assertEqual((V, mom), (1, (mp.mpf(1)/2,)*3))
        self.assertLess(abs(A-6), mp.mpf(10)**-35)
        self.assertTrue(cube.contains((mp.mpf(1)/2,)*3))
        self.assertFalse(cube.contains((mp.mpf(3)/2, mp.mpf(1)/2, mp.mpf(1)/2)))
        # A tetrahedron of skew points: its volume the determinant's sixth,
        # every plane through three of them with the fourth inside.
        pts = [(0, 0, 0), (5, 1, F(1, 3)), (1, 4, 2), (F(3, 2), 1, 5)]
        tet = ref.Hull(ref.HullCase(pts))
        self.assertEqual(len(tet.cycles), 4)
        V, _, _ = tet.closed()
        a, b, c, d = (tuple(F(x) for x in p) for p in pts)
        self.assertLess(abs(V-pr.M(abs(pr.det3(pr.sub(b, a), pr.sub(c, a), pr.sub(d, a)))/6)), mp.mpf(10)**-35)
        for n, dd in tet.planes:
            inner = [p for p in tet.points if pr.dot(n, p) < dd]
            self.assertEqual(len(inner), 1)
        # Its sweep (both families of every face) its closed form.
        scan = ref.ch.SCAN
        ref.ch.SCAN = 24
        try:
            chain = ref.Chain([ref.HullCase(pts)])
            for k in (0, 1):
                v, m, area = chain.measures(0, k)
                V, M_, A = chain.inputs[0].closed()
                self.assertLess(abs(v-V), mp.mpf(10)**-30)
                self.assertLess(abs(area-A), mp.mpf(10)**-30)
                self.assertLess(max(abs(m[i]-M_[i]) for i in range(3)), mp.mpf(10)**-30)
        finally:
            ref.ch.SCAN = scan
        # Its binary64 rays: the cube's interval along `x` at its middle.
        fm = ref.Float(cube)
        self.assertEqual(fm.intervals([-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 3.0), [(1.0, 2.0)])
        self.assertTrue(fm.contains([0.5, 0.25, 0.75]))
        self.assertFalse(fm.contains([0.5, 1.25, 0.75]))

    def test_bodies(self):
        # S9e.4b.2's convex bodies: their hulls the constructions OCCT was
        # given (a wedge's corners on its frame's stored axes, a polyhedron's
        # binary64 points), their volumes S9e.4b.2's closed forms.
        for name in ('pyramid', 'truncated', 'wedge', 'tetra', 'octa'):
            b = fixtures.BODIES[name]
            self.assertEqual(b.path, earlier.BODIES[name].path)
            (kind, points), = b.specs
            self.assertEqual(kind, 'hull')
            cell, = earlier.BODIES[name].body().cells
            want = {p for c in cell for p in c}
            self.assertEqual(set(points), want, name)
            V, _, _ = ref.Hull(ref.HullCase(points)).closed()
            self.assertLess(abs(V-pr.M(F(earlier.BODIES[name].closed()))), mp.mpf(10)**-35, name)
        # The tetrahedron's and the octahedron's points as S9e.4b.2 lists
        # them.
        self.assertEqual(set(fixtures.BODIES['tetra'].specs[0][1]),
                         {tuple(F(c) for c in p) for p in earlier.TETRA})
        self.assertEqual(set(fixtures.BODIES['octa'].specs[0][1]),
                         {tuple(F(c) for c in p) for p in earlier.OCTA})
        # S9e.4b.2's Booleans of boxes: its rows' boxes.
        for name in ('notched', 'hollow'):
            rows = earlier.BODIES[name].rows
            a, op, c = fixtures.BOXES[name]
            self.assertEqual(rows, [earlier.box_row(*a), f'boolean {op}', earlier.box_row(*c)])
        self.assertEqual(fixtures.BODIES['hollow'].closed(), 1000-64)
        # S9e.4a's ball.
        self.assertEqual(fixtures.BODIES['ball'].specs, [earlier_curved.BODIES['ball'].spec])
        self.assertEqual(fixtures.BODIES['ball'].path, 'imported/ball.brep')

    def test_cavity(self):
        (fo, so), (fi, si) = fixtures.CAVITY_OUTER, fixtures.CAVITY_INNER
        # The inner box's corners inside the outer one by at least 1.35 of
        # its local coordinates.
        least = min(min(c[i], F(so[i])-c[i]) for p in corners(fi, si) for c in [local(fo, p)] for i in range(3))
        self.assertGreater(least, F(135, 100))
        # Its volume: the outer box's less the inner one's, each times its
        # stored axes' determinant (each a unit within rounding).
        V = fixtures.cavity_volume()
        for frame in (fo, fi):
            self.assertLess(abs(fixtures.det(frame)-1), mp.mpf(10)**-15)
        self.assertLess(abs(V-(480-15)), mp.mpf(10)**-12)
        hull_volume = sum((pr.volume([c for c in ipr.hull(corners(f, s))])[0]*(1 if k == 0 else -1)
                           for k, (f, s) in enumerate((fixtures.CAVITY_OUTER, fixtures.CAVITY_INNER))), F(0))
        self.assertLess(abs(V-pr.M(hull_volume)), mp.mpf(10)**-30)
        # The slab between the inner box's faces normal to `TILT`'s `n`: its
        # heights in the inner box's chart strictly inside `(0, 2)`, its
        # sides past the outer box's corners.
        slab = fixtures.CAVITY_SLAB
        frame, h0, h1 = slab[2], slab[3], slab[4]
        self.assertEqual(frame[3:], fi[3:])
        lo, hi = (local(fi, local_point(frame, 0, 0, h))[2] for h in (h0, h1))
        self.assertTrue(0 < lo < hi < 2, (lo, hi))
        x0, y0 = slab[1][0].points[0]
        x1, y1 = slab[1][0].points[2]
        for p in corners(fo, so):
            u, v, _ = local(frame, p)
            self.assertTrue(x0 < u < x1 and y0 < v < y1, (u, v))

    def test_partners(self):
        # The hollow box's ball: through its wall `x = 10` and its cavity's
        # wall `x = 7` within those faces, tangent to no plane of either box.
        r, frame = fixtures.HOLLOW_BALL[1], fixtures.HOLLOW_BALL[2]
        c = tuple(F(v) for v in frame[:3])
        r = F(r)
        for axis, values in ((0, (0, 3, 7, 10)), (1, (0, 3, 7, 10)), (2, (0, 3, 7, 10))):
            for v in values:
                self.assertNotEqual(abs(c[axis]-v), r, (axis, v))
        self.assertTrue(c[0]-r < 7 < c[0] < 10 < c[0]+r)
        for x in (7, 10):
            rho2 = r*r-(c[0]-x)**2
            # Its circle there inside the face (`[3, 7]` or `[0, 10]` in
            # `y` and `z`) by squares.
            lo, hi = (3, 7) if x == 7 else (0, 10)
            for axis in (1, 2):
                self.assertTrue((c[axis]-lo)**2 > rho2 and (hi-c[axis])**2 > rho2, (x, axis))
        # The inner ball inside the cavity `[3, 7]^3`.
        r, frame = F(fixtures.HOLLOW_INNER[1]), fixtures.HOLLOW_INNER[2]
        for v in frame[:3]:
            self.assertTrue(3 < F(v)-r and F(v)+r < 7)
        # The touching ball on the tetrahedron's base: its corners' heights
        # -1/2 (`TURN30`'s normal the world's `z`), the ball's lowest point
        # there inside the base triangle.
        base = [p for p in fixtures.BODIES['tetra'].specs[0][1] if p[2] == F(-1, 2)]
        self.assertEqual(len(base), 3)
        r, frame = F(fixtures.TETRA_TOUCH[1]), fixtures.TETRA_TOUCH[2]
        foot = (F(frame[0]), F(frame[1]), F(frame[2])-r)
        self.assertEqual(foot[2], F(-1, 2))
        n = pr.normal_of(base)
        sides = [pr.dot(pr.cross(pr.sub(base[(i+1) % 3], base[i]), pr.sub(foot, base[i])), n) for i in range(3)]
        self.assertTrue(all(s > 0 for s in sides) or all(s < 0 for s in sides))

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 34)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 30, 'empty': 0, 'degenerate': 3, 'unsupported': 1})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        chain = next(c for c in listed if c.name == 'chain_pyramid_fuse')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('fuse', ('cut', 0, 1), 2))
        swapped = next(c for c in listed if c.name == 'rod_truncated_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, 1))
        hollow = next(c for c in listed if c.name == 'hollow_ball_cut')
        self.assertEqual(hollow.expr(), ('cut', ('cut', 0, 1), 2))
        both = next(c for c in listed if c.name == 'tetra_ball_common')
        self.assertEqual(both.encode().count('\nbrep imported/'), 2)
        self.assertEqual(both.native().count('\nbrep imported/'), 2)
        # This step's body written once: two boxes about a Boolean; the
        # others the earlier steps' files.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.OWN))
        self.assertEqual(text.count('\nboolean '), len(fixtures.OWN))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        mine = {b.path for b in fixtures.OWN}
        others = {b.path for b in fixtures.BODIES.values()}-mine
        self.assertEqual(others, {earlier.BODIES[n].path for n in fixtures.BODIES if n in earlier.BODIES}
                         | {'imported/ball.brep'})
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_pieces_boolean_fixtures as d
        import generate_one_sphere_boolean_fixtures as e
        import generate_one_sphere_incidence_boolean_fixtures as f
        import generate_piece_forms_boolean_fixtures as g
        import generate_piece_trees_boolean_fixtures as h
        import generate_imported_joints_boolean_fixtures as i
        import generate_primitive_chains_boolean_fixtures as j
        import generate_prism_leaves_boolean_fixtures as k
        import generate_plane_parts_boolean_fixtures as m
        for older in (earlier_curved.BODIES, b.BODIES, earlier.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES,
                      h.BODIES, i.BODIES, j.BODIES, k.BODIES, m.BODIES):
            self.assertFalse(mine & {x.path for x in older.values()})


def local_point(frame, u, v, w):
    """The frame's point `o + u x + v y + w n` on its stored axes, exactly."""
    o, x, y, n = (tuple(F(c) for c in a) for a in stored_axes(frame))
    return tuple(o[i]+F(u)*x[i]+F(v)*y[i]+F(w)*n[i] for i in range(3))


if __name__ == '__main__':
    unittest.main()
