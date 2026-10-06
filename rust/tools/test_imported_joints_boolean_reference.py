"""Independent checks of S9e.4b.4a's fixtures
(`generate_imported_joints_boolean_fixtures.py`): every profile's arcs end on
their circles and every joint of two circles lies on both in the
constructions (exactly, in rationals), the chained reference on each body
alone against its profile's closed form (a polygon and its arcs' circular
segments), the quad less nothing but a box (their common) against a direct
quadrature of its sections, the joints' off-circle test on points known to
be on and off, and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_imported_joints_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import square


def mpq(x):
    return mp.mpf(x.numerator)/x.denominator


def profile_area(boundary):
    """A path's area in closed form: its points' polygon (the shoelace sum)
    and each arc's circular segment beyond its chord, `r^2 (t - sin t) / 2`
    for the arc's angle `t`, added where the arc turns the profile's way
    (counter-clockwise) and taken off where it turns against it."""
    pts = [tuple(mp.mpf(c) for c in p) for p in boundary.points]
    n = len(pts)
    area = sum(pts[j][0]*pts[(j+1) % n][1]-pts[(j+1) % n][0]*pts[j][1] for j in range(n))/2
    for j, s in enumerate(boundary.segments):
        if s is None:
            continue
        cx, cy, r, ccw = (mp.mpf(s[0]), mp.mpf(s[1]), mp.mpf(s[2]), s[3])
        a0 = mp.atan2(pts[j][1]-cy, pts[j][0]-cx)
        a1 = mp.atan2(pts[(j+1) % n][1]-cy, pts[(j+1) % n][0]-cx)
        t = (a1-a0) % (2*mp.pi) if ccw else (a0-a1) % (2*mp.pi)
        segment = r**2*(t-mp.sin(t))/2
        area += segment if ccw else -segment
    return area


class ImportedJointsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_joints_lie_on_both_circles(self):
        # Each body OCCT was given ends every arc on its circle and meets two
        # circles on both exactly (its profile's binary64 points): only the
        # files' roundings leave them.
        seen = 0
        for b in fixtures.BODIES.values():
            for bound in b.spec[1]:
                n = len(bound.points)
                for j, s in enumerate(bound.segments):
                    if s is None:
                        continue
                    for p in (bound.points[j], bound.points[(j+1) % n]):
                        self.assertTrue(fixtures.on_circle(p, s[:3]), (b.name, p))
            for p, circles in fixtures.joints(b.spec):
                seen += 1
                for c in circles:
                    self.assertTrue(fixtures.on_circle(p, c), (b.name, p))
        # The quad's 4, the arch's 1, the cam's 1, the blade's 1, the split
        # lens's 2.
        self.assertEqual(seen, 9)

    def test_tangent_and_crossing_joints(self):
        # The cam's joint is one circle touching the other inside (centres
        # apart by the radii's difference); the arch's and the blade's two
        # circles crossing, at asin(24 / 25) (radii along (5, 4) and (-5,
        # 4) less (+-3, 0)) and asin(36 / 325), 0.11 rad (radii along (24, 7)
        # and (12, 5)).
        def joint(name):
            (p, (c0, c1)), = fixtures.joints(fixtures.BODIES[name].spec)
            return p, c0, c1

        def angle(p, c0, c1):
            u = [(mp.mpf(p[i])-c[i])/c[2] for c in (c0, c1) for i in range(2)]
            return mp.asin(abs(u[0]*u[3]-u[1]*u[2]))
        _, c0, c1 = joint('cam')
        self.assertEqual((F(c0[0])-F(c1[0]))**2+(F(c0[1])-F(c1[1]))**2, (F(c0[2])-F(c1[2]))**2)
        self.assertLess(abs(angle(*joint('arch'))-mp.asin(mp.mpf(24)/25)), mp.mpf(10)**-35)
        self.assertLess(abs(angle(*joint('blade'))-mp.asin(mp.mpf(36)/325)), mp.mpf(10)**-35)
        # The quad's corners, radii along (4, 3) and (3, 4): asin(7 / 25).
        for p, (c0, c1) in fixtures.joints(fixtures.BODIES['quad'].spec):
            self.assertLess(abs(angle(p, c0, c1)-mp.asin(mp.mpf(7)/25)), mp.mpf(10)**-35)

    def test_bodies_alone_are_their_closed_forms(self):
        # Each body alone (beside a box far away): the chained reference's
        # volume its profile's area times its height times the stored axes'
        # determinant (the stored axes are binary64 roundings of an
        # orthonormal frame: volumes scale by it exactly).
        scan, ref.SCAN = ref.SCAN, 48
        try:
            far = fixtures.construction(fixtures.prism([square(60.0, 60.0, 61.0, 61.0)], fixtures.at('XY', (0, 0, 0)),
                                                       0.0, 1.0), 92)
            for name in ('quad', 'arch', 'cam', 'blade'):
                spec = fixtures.BODIES[name].spec
                chain = ref.Chain([fixtures.construction(spec, 91), far])
                v, _, _ = chain.measures(0)
                _, x, y, n = (tuple(F(c) for c in w) for w in fixtures.stored_axes(spec[2]))
                det = (x[0]*(y[1]*n[2]-y[2]*n[1])-x[1]*(y[0]*n[2]-y[2]*n[0])+x[2]*(y[0]*n[1]-y[1]*n[0]))
                want = profile_area(spec[1][0])*(spec[4]-spec[3])*mpq(det)
                self.assertLess(abs(v-want)/want, mp.mpf(10)**-30, name)
        finally:
            ref.SCAN = scan

    def test_quad_box_common_by_its_sections(self):
        # The quad in common with the box `[6, 12]^2 x [1, 2.5]`: 1.5 times
        # the area of the quad's section inside `[6, 12]^2`, each world line
        # `x = const` meeting every disc in a chord (the quad's frame taken as
        # the exact rotation by (12, 5) / 13: the stored axes are its
        # roundings, a relative 1e-16).
        scan, ref.SCAN = ref.SCAN, 48
        try:
            spec = fixtures.BODIES['quad'].spec
            chain = ref.Chain([fixtures.construction(spec, 91), fixtures.construction(fixtures.BOX, 92)])
            V, _, _ = chain.measures(('common', 0, 1))
        finally:
            ref.SCAN = scan
        c, s = mp.mpf(12)/13, mp.mpf(5)/13
        world = lambda u, v: (5+c*u-s*v, 5+s*u+c*v)
        discs = [world(*p) for p in ((1, 0), (-1, 0), (0, 1), (0, -1))]

        def length(x):
            lo, hi = mp.mpf(6), mp.mpf(12)
            for cx, cy in discs:
                d = 25-(x-cx)**2
                if d <= 0:
                    return mp.mpf(0)
                lo, hi = max(lo, cy-mp.sqrt(d)), min(hi, cy+mp.sqrt(d))
            return max(hi-lo, 0)
        # Breakpoints: the corners', the discs' extremes' and where a disc's
        # chord crosses y = 6 or y = 12, inside (6, 12).
        xs = [world(*p)[0] for p in ((3, 3), (3, -3), (-3, 3), (-3, -3))]
        xs += [cx+s*5 for cx, _ in discs for s in (-1, 1)]
        for cx, cy in discs:
            for y in (6, 12):
                d = 25-(y-cy)**2
                if d > 0:
                    xs += [cx-mp.sqrt(d), cx+mp.sqrt(d)]
        breaks = [mp.mpf(6)]+sorted(x for x in xs if 6 < x < 12)+[mp.mpf(12)]
        area = mp.quad(length, breaks)
        self.assertLess(abs(V-mp.mpf(3)/2*area)/V, mp.mpf(10)**-12)

    def test_off_joint_test(self):
        # The joints' test: the quad's own corners (on both circles), a corner
        # moved by an ulp (off one).
        circles = [(-1.0, 0.0, 5.0), (0.0, -1.0, 5.0)]
        self.assertTrue(all(fixtures.on_circle((3.0, 3.0), c) for c in circles))
        self.assertFalse(all(fixtures.on_circle((3.0000000000000004, 3.0), c) for c in circles))
        b = fixtures.BODIES['quad']
        local = fixtures.local_map(b.spec[2])
        o, x, y, n = (tuple(F(c) for c in v) for v in fixtures.stored_axes(b.spec[2]))
        for p, _ in fixtures.joints(b.spec):
            world = tuple(o[i]+F(p[0])*x[i]+F(p[1])*y[i] for i in range(3))
            u = local(world)
            self.assertEqual((u[0], u[1]), (F(p[0]), F(p[1])))

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 39)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 33, 'empty': 0, 'degenerate': 3, 'unsupported': 3})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('quad_box_fuse', 'box_quad_cut', 'both_common', 'chain_quad_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.specs)
            self.assertEqual(enc.count('\nbrep imported/'), want, name)
            self.assertEqual(nat.count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_quad_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        # Every body written once, by its rows: a plane, a path and a prism.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nprism '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.specs if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of an earlier step's bodies' names.
        import generate_imported_boolean_fixtures as s9e4a
        import generate_imported_arcs_boolean_fixtures as s9e4b1
        self.assertFalse(set(fixtures.BODIES) & (set(s9e4a.BODIES) | set(s9e4b1.BODIES)))


if __name__ == '__main__':
    unittest.main()
