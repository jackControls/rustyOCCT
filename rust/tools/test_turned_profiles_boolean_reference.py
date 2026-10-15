"""Independent checks of S9e.4b.4d.2a's fixtures
(`generate_turned_profiles_boolean_fixtures.py`): the bodies' closed forms
(Pappus) against direct quadratures of their sections' areas and their
meridians' lengths, the lattice's max-min form against the profile at
heights of its own (and a wrong form caught), the waist alone on the chained
reference against its closed form, every profile's shape (dyadic, rising
from the axis to the axis, corners only, no two cones on one apex or of
parallel rulings), the partners' places, the declared case (a ball resting
on the waist's top disc) and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_profiles_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def radius(points, w):
    """The profile's radius at height `w` (linear between its points; the
    larger at a step)."""
    best = mp.mpf(0)
    for (ra, wa), (rb, wb) in zip(points, points[1:]):
        if mpq(wa) <= w <= mpq(wb) and wb > wa:
            best = max(best, mpq(ra)+(mpq(rb)-mpq(ra))*(w-mpq(wa))/(mpq(wb)-mpq(wa)))
    return best


class TurnedProfilesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # The volume `pi int r(w)^2 dw` and the area: the meridian's segments'
        # `2 pi int r ds` (a level one's an annulus's) by quadrature.
        tight = mp.mpf(10)**-25
        for b in fixtures.BODIES.values():
            pts = b.lattice.points
            V, A = b.closed()
            heights = sorted({mpq(w) for _, w in pts})
            vol = mp.pi*mp.quad(lambda w: radius(pts, w)**2, heights)
            area = mp.mpf(0)
            for (ra, wa), (rb, wb) in zip(pts, pts[1:]):
                ra, wa, rb, wb = mpq(ra), mpq(wa), mpq(rb), mpq(wb)
                length = mp.sqrt((rb-ra)**2+(wb-wa)**2)
                area += 2*mp.pi*mp.quad(lambda t: (ra+(rb-ra)*t)*length, [0, 1])
            self.assertLess(abs(V-vol), tight*V, b.name)
            self.assertLess(abs(A-area), tight*A, b.name)
        # The point's volume `pi (9 + 27 + 185/6 + 120)`, the waist's `pi (640
        # + 172 + 150)`.
        V, _ = fixtures.BODIES['point'].closed()
        self.assertLess(abs(V-mp.pi*(9+27+mpq(F(185, 6))+120)), tight)
        V, _ = fixtures.BODIES['waist'].closed()
        self.assertLess(abs(V-mp.pi*(640+172+150)), tight)

    def test_lattice_forms(self):
        # The form against the profile at heights of this test's own (a
        # hundred and one across each body's span and beyond), and a form
        # missing a term caught.
        for b in fixtures.BODIES.values():
            L = b.lattice
            lo, hi = L.w_min-2, L.w_max+2
            for k in range(101):
                w = lo+(hi-lo)*F(k, 100)+F(1, 997)
                self.assertEqual(L.form(w), L.profile(w), (b.name, w))
            self.assertEqual(sorted(leaves(L.expr())), sorted(j for t in L.terms for j in t))
            if len(L.terms) > 1:
                terms = L.terms
                L.terms = terms[1:]
                with self.assertRaises(AssertionError):
                    L.check()
                L.terms = terms
        # The point's terms, made as small as the form allows: the apex cone
        # common the first cylinder, fused with the widening cone common the
        # last cylinder (`bug28626`'s solids' form).
        self.assertEqual(fixtures.BODIES['point'].lattice.terms, [[0, 1], [2, 3]])
        # Each term holds its own face's part.
        for b in fixtures.BODIES.values():
            for i, t in enumerate(b.lattice.every):
                self.assertIn(i, t)

    def test_waist_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['waist']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, A = chain.measures(b.ref_expr)
            want, area = b.closed()
            self.assertLess(abs(V-want), mp.mpf(10)**-28*want)
            self.assertLess(abs(A-area), mp.mpf(10)**-28*area)
        finally:
            ref.SCAN = scan

    def test_profiles(self):
        for b in fixtures.BODIES.values():
            pts = b.lattice.points
            self.assertEqual((pts[0][0], pts[-1][0]), (0, 0), b.name)
            for r, w in pts:
                for x in (r, w):
                    self.assertEqual(F(float(x)), x)
            # Rising or level, never falling; no two level segments in a row
            # (a step is one annulus); corners only (consecutive rising
            # segments of different slopes).
            segs = list(zip(pts, pts[1:]))
            for (a, b_), (c, d) in zip(segs, segs[1:]):
                self.assertGreaterEqual(b_[1], a[1])
                level = (a[1] == b_[1], c[1] == d[1])
                self.assertNotEqual(level, (True, True), b.name)
                if level == (False, False):
                    self.assertNotEqual((b_[0]-a[0])/(b_[1]-a[1]), (d[0]-c[0])/(d[1]-c[1]), b.name)
            # Its cones: no two on one apex, no two of slopes `k` and `-k`
            # (rulings of one direction), none within the resolution of flat.
            slopes = [p['k'] for p in b.lattice.parts if p['k'] != 0]
            apexes = [-p['r0']/p['k'] for p in b.lattice.parts if p['k'] != 0]
            self.assertEqual(len(set(apexes)), len(apexes), b.name)
            for k in slopes:
                self.assertNotIn(-k, slopes, b.name)
                self.assertGreater(abs(k), F(1, 4))
            # Every reference solid on the axis: a cone's frame's origin on
            # it, its radii its line's.
            for spec, p in zip(b.ref_specs, b.lattice.parts):
                lo, hi = p['range']
                if spec[0] == 'cone':
                    _, r0, r1, h, frame = spec
                    self.assertEqual((F(r0), F(r1), F(h)), (b.lattice.at(p, lo), b.lattice.at(p, hi), hi-lo))
                    self.assertEqual(frame, b.frame_at(lo))
                else:
                    self.assertEqual((F(spec[3]), F(spec[4])), (lo, hi))

    def test_partners_placed(self):
        # The point's box across its first corner (radius 3 at height 3): the
        # corner's point `(3, 0, 3)` at the box's chart point `(1, 1/5,
        # 53/20)` in the `TILT` frame's exact axes `x`, `(0, 4, -3) / 5`, `(0,
        # 3, 4) / 5`, inside it by 1/5 at least.
        _, (bx,), frame, h0, h1 = fixtures.POINT_BOX
        self.assertEqual(frame[3:], fixtures.FRAMES['TILT'])
        o = [F(c) for c in frame[:3]]
        u, v, h = F(1), F(1, 5), F(53, 20)
        corner = (o[0]+u, o[1]+v*F(4, 5)+h*F(3, 5), o[2]-v*F(3, 5)+h*F(4, 5))
        self.assertEqual(corner, (3, 0, 3))
        xs = sorted({F(p[0]) for p in bx.points})
        ys = sorted({F(p[1]) for p in bx.points})
        for x, (lo, hi) in ((u, (xs[0], xs[-1])), (v, (ys[0], ys[-1])), (h, (F(h0), F(h1)))):
            self.assertTrue(lo+F(1, 5) <= x <= hi-F(1, 5))
        # The waist's rod along `x` about `(y, z) = (3/4, 12)` through its cone
        # (heights 10 to 14), off the waist's axis, through the whole body.
        _, (c,), frame, h0, h1 = fixtures.WAIST_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame[3:], fixtures.FRAMES['SIDE'])
        self.assertTrue(10 < cy-r and cy+r < 14 and cx == 0 and frame[1] == 0.75)
        self.assertTrue(frame[0]+h0 < -8 and frame[0]+h1 > 8)
        # The steps' ball across the annulus at height 5 (radii 3 to 5).
        _, r, frame, _, _ = fixtures.STEPS_BALL
        d = mp.sqrt(mpq(frame[0])**2+mpq(frame[1])**2)
        self.assertTrue(abs(mpq(frame[2])-5) < r and 3 < d < 5)

    def test_declared(self):
        # The touching ball: radius 2 about `(3/2, -1, 22)`, resting on the
        # waist's top disc (radius 5 at 20) at `(3/2, -1, 20)`, inside the
        # disc by more than its radius's quarter, clear of the cylinder below.
        _, r, frame, _, _ = fixtures.WAIST_TOUCH
        self.assertEqual(r, 2.0)
        self.assertEqual(frame[2]-r, 20.0)
        d = mp.sqrt(mpq(frame[0])**2+mpq(frame[1])**2)
        self.assertLess(d+mp.mpf(1)/4, 5)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'waist_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 21)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 18, 'empty': 0, 'degenerate': 3, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name, want in (('point_box_fuse', 1), ('box_crown_cut', 1), ('point_tilt_common', 2),
                           ('chain_waist_cut', 1)):
            c = next(c for c in listed if c.name == name)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_waist_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        n = len(fixtures.BODIES['waist'].ref_specs)
        self.assertEqual(chain.expr(), ('cut', ('cut', fixtures.BODIES['waist'].ref_expr, n), n+1))
        swapped = next(c for c in listed if c.name == 'box_crown_cut')
        self.assertEqual(swapped.expr(), ('cut', 0, fixtures.BODIES['crown'].lattice.expr(1)))
        # Each body written once, its profile turned a whole turn.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nrevolve '), len(fixtures.BODIES))
        self.assertEqual(text.count(' 6.283185307179586\n'), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual({s[1] for c in listed for s in c.items if s[0] == 'base'}, {'box_tilt'})
        self.assertEqual(used, set(fixtures.BODIES))
        for b in fixtures.BODIES.values():
            self.assertEqual(sorted(set(leaves(b.ref_expr))), list(range(len(b.ref_specs))))
        # None of the earlier turned bodies' files.
        import write_turned_bodies as w
        paths = {b.path for b in fixtures.BODIES.values()}
        self.assertFalse(paths & {f'imported/{n}.brep' for n in w.BODIES})


if __name__ == '__main__':
    unittest.main()
