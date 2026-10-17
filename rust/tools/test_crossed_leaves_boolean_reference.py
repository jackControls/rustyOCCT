"""Independent checks of S9e.4b.4b.2b.2b.3b's fixtures
(`generate_crossed_leaves_boolean_fixtures.py`): the bodies' closed forms
against direct quadratures of their sections, the chained reference on the
waist against its closed form, every part's place (each slab across its
plate's width and height, the joint's across the plate's line-arc joints,
the waist's across its straight middle; the plates' walls along the slabs'
axes exactly), the partners' places, the declared case (a ball resting on
the joint's top), and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
from curve_surface_reference import stored_axes
import generate_crossed_leaves_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def stadium_width(length, r, y):
    """A stadium's width at `y` about its centre line (half discs of radius
    `r` about `(0, 0)` and `(length, 0)`)."""
    return length+2*mp.sqrt(max(mp.mpf(0), r**2-y**2))


def stadium_span(length, r, x):
    """A stadium's height at `x` along its centre line (`2 sqrt(r^2 - x^2)`
    over its half discs, `2 r` between)."""
    if 0 <= x <= length:
        return 2*mp.mpf(r)
    d = x if x < 0 else x-length
    return 2*mp.sqrt(max(mp.mpf(0), r**2-d**2))


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


class CrossedLeavesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        det = fixtures.det
        half = mp.mpf(3)/2
        slab = mp.quad(lambda y: stadium_width(4, 2, y), [-2, 0, 2])
        # The joint: the plate's section by its `x` (its span `stadium_span`)
        # over `3/2`, the slab's section over its unit length, less the
        # plate's section over the slab's `[-1/4, 3/4]`.
        V, _ = fixtures.BODIES['joint'].closed()
        plate = mp.quad(lambda x: stadium_span(4, 2, x), [-2, 0, 4, 6])
        common = mp.quad(lambda x: stadium_span(4, 2, x), [-mp.mpf(1)/4, 0, mp.mpf(3)/4])
        want = half*(plate-common)*det(fixtures.JOINT_PLATE)+slab*det(fixtures.JOINT_SLAB)
        self.assertLess(abs(V-want), tight)
        # The waist: the plate of length 8 less its section over `[3, 4]`.
        V, _ = fixtures.BODIES['waist'].closed()
        plate = mp.quad(lambda x: stadium_span(8, 2, x), [-2, 0, 8, 10])
        common = mp.quad(lambda x: stadium_span(8, 2, x), [3, 4])
        want = half*(plate-common)*det(fixtures.WAIST_PLATE)+slab*det(fixtures.WAIST_SLAB)
        self.assertLess(abs(V-want), tight)
        for frame in (fixtures.JOINT_PLATE, fixtures.JOINT_SLAB, fixtures.WAIST_PLATE, fixtures.WAIST_SLAB):
            self.assertLess(abs(det(frame)-1), mp.mpf(10)**-15)

    def test_waist_on_the_chained_reference(self):
        # The waist alone on the chained reference (the bodies on whole cases
        # are the generator's check) against its closed form.
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['waist']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, _ = chain.measures(b.ref_expr)
            self.assertLess(abs(V-b.closed()[0]), mp.mpf(10)**-28)
        finally:
            ref.SCAN = scan

    def test_parts_placed(self):
        for name, plate_frame, slab_frame, length, x0, x1 in (
                ('joint', fixtures.JOINT_PLATE, fixtures.JOINT_SLAB, 4, F(-1, 4), F(3, 4)),
                ('waist', fixtures.WAIST_PLATE, fixtures.WAIST_SLAB, 8, F(3), F(4))):
            b = fixtures.BODIES[name]
            plate, slab = b.first, b.second
            self.assertEqual(plate[2], plate_frame)
            self.assertEqual(slab[2], slab_frame)
            self.assertEqual((plate[3], plate[4]), (0.0, 1.5))
            # The plate on the world's axes, the slab's normal the world's
            # `x` and its chart the plate's `(y, z)`: the plate's walls along
            # the slab's axis exactly.
            _, px, py, pn = stored_axes(plate_frame)
            _, sx, sy, sn = stored_axes(slab_frame)
            self.assertEqual((px, py, pn), ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
            self.assertEqual((sn, sx, sy), (px, py, pn))
            # The slab's range along the plate's `x` and its origin, exactly.
            u, v, w = chart(plate_frame, slab_frame[:3])
            self.assertEqual((u+F(slab[3]), u+F(slab[4])), (x0, x1))
            self.assertEqual((v, w), (-1, F(3, 4)))
            # The slab's section (half discs of radius 2 about the plate's `(y,
            # z) = (-1, 3/4)` and `(3, 3/4)`) covers the plate's `[-2, 2] x
            # [0, 3/2]` by at least 1/4.
            for y, z in ((-2, 0), (2, 0), (-2, F(3, 2)), (2, F(3, 2))):
                gap = stadium_gap(4, 2, fixtures.q(y+1), fixtures.q(z-F(3, 4)))
                self.assertGreater(gap, mp.mpf(1)/4)
            # Across the joint's line-arc joints at `x = 0` (inside by 1/4
            # either way), or within the waist's straight part.
            if name == 'joint':
                self.assertTrue(x0 < 0 < x1 and min(-x0, x1) >= F(1, 4))
            else:
                self.assertTrue(0 < x0 < x1 < length)

    def test_partners_placed(self):
        # Each rod on its plate's frame, across a slab's cap inside the
        # plate's outline, through both plate caps.
        for rod, frame, length, cap in ((fixtures.JOINT_ROD, fixtures.JOINT_PLATE, 4, F(3, 4)),
                                        (fixtures.JOINT_ARC, fixtures.JOINT_PLATE, 4, F(-1, 4)),
                                        (fixtures.WAIST_ROD, fixtures.WAIST_PLATE, 8, F(4))):
            _, (c,), f, h0, h1 = rod
            cx, cy, r = c.circle
            self.assertEqual(f, frame)
            self.assertTrue(F(cx)-F(r) < cap < F(cx)+F(r))
            self.assertGreater(stadium_gap(length, 2, mp.mpf(cx), mp.mpf(cy)), r+mp.mpf(1)/8)
            self.assertTrue(h0 < 0 and h1 > 1.5)
        # The tool's rod across both slab caps, the plate's wall `y = 2` and
        # its left arc, between the plate's caps.
        _, (c,), f, h0, h1 = fixtures.ROD_JOINT
        cx, cy, r = c.circle
        self.assertTrue(cx-r < -0.25 and cx+r > 0.75 and cy-r < 2 < cy+r and 0 < h0 < h1 < 1.5)
        self.assertLess(abs(mp.sqrt(mp.mpf(cx)**2+mp.mpf(cy)**2)-2), r)
        # The chain's rod across the waist's slab cap `x = 3`, its box across
        # the plate's right piece and the rod's hole.
        _, (c,), f, _, _ = fixtures.CHAIN_ROD
        cx, cy, r = c.circle
        self.assertTrue(cx-r < 3 < cx+r and abs(cy)+r < 2)
        _, (b,), f, h0, h1 = fixtures.CHAIN_BOX
        xs = [p[0] for p in b.points]
        self.assertTrue(min(xs) < cx < 4 < max(xs) and h0 < 1.5 < h1)

    def test_declared(self):
        # The touching ball: radius 1 about the joint plate's point `(5/2,
        # 0)` 1 above its top (rounded once into the world), off the slab.
        touch = fixtures.JOINT_TOUCH
        self.assertEqual(touch[1], 1.0)
        u, v, w = chart(fixtures.JOINT_PLATE, touch[2][:3])
        self.assertLess(abs(w-F(5, 2)), F(1, 10**14))
        self.assertTrue(abs(u-F(5, 2)) < F(1, 10**14) and abs(v) < F(1, 10**14))
        self.assertEqual(fixtures.BODIES['joint'].first[4]+touch[1], 2.5)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'joint_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('joint_rod_fuse', 1), ('rod_joint_cut', 1), ('joint_cyl_common', 2),
                           ('chain_waist_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_waist_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('fuse', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'joint_cyl_cut')
        self.assertEqual(both.expr(), ('cut', ('fuse', 0, 1), 2))
        swapped = next(c for c in listed if c.name == 'rod_joint_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, ('fuse', 1, 2)))
        # Each body written once, two prisms each, fused.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean fuse'), len(fixtures.BODIES))
        self.assertEqual(text.count('\nprism '), 2*len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(leaves(b.ref_expr)), list(range(len(b.ref_specs))))
        # None of the earlier steps' bodies' files.
        import generate_imported_boolean_fixtures as a
        import generate_two_leaves_boolean_fixtures as b
        import generate_rim_leaves_boolean_fixtures as c
        import generate_flush_leaves_boolean_fixtures as d
        import generate_stack_leaves_boolean_fixtures as e
        import generate_nested_leaves_boolean_fixtures as f
        import generate_prism_leaves_boolean_fixtures as g
        import generate_plane_parts_boolean_fixtures as h
        paths = {b_.path for b_ in fixtures.BODIES.values()}
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES, f.BODIES, g.BODIES, h.BODIES):
            self.assertFalse(paths & {x.path for x in earlier.values()})


if __name__ == '__main__':
    unittest.main()
