"""Independent checks of the moved turned bodies' fixtures
(`generate_turned_tilted_boolean_fixtures.py`): the turned bodies' closed
forms (Pappus over the model's profile, the height's unit 3) against direct
quadratures, the pin's model profile against OCCT's, the tilted frame's
axes, the lattices' forms against the model profiles, the waist alone on the
chained reference against its closed form, the box's closed form and
corners, the partners' places (none parallel or normal to the axis) and the
case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_tilted_boolean_fixtures as fixtures


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


TURNED = [b for b in fixtures.BODIES.values() if isinstance(b, fixtures.Body)]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


class TurnedTiltedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # The volume `pi int r(h)^2 dh` over the world's heights `h = 3 w`
        # and the lateral areas `2 pi int r ds` by quadrature.
        tight = mp.mpf(10)**-25
        for b in TURNED:
            pts = b.model
            V, A = b.closed()
            heights = sorted({mpq(w) for _, w in pts})
            vol = 3*mp.pi*mp.quad(lambda w: radius(pts, w)**2, heights)
            area = mp.mpf(0)
            for (ra, wa), (rb, wb) in zip(pts, pts[1:]):
                ra, wa, rb, wb = mpq(ra), mpq(wa), mpq(rb), mpq(wb)
                length = mp.sqrt((rb-ra)**2+9*(wb-wa)**2)
                area += 2*mp.pi*mp.quad(lambda t: (ra+(rb-ra)*t)*length, [0, 1])
            self.assertLess(abs(V-vol), tight*V, b.name)
            self.assertLess(abs(A-area), tight*A, b.name)
        # The waist's volume `pi (36 * 3 + 3 (36 + 24 + 16) / 3 + 16 * 3)`.
        V, _ = fixtures.BODIES['waist'].closed()
        self.assertLess(abs(V-mp.pi*(108+76+48)), tight)
        # The box's: its edges 3, 11/4 and 13/4.
        V, A = fixtures.BODIES['skewbox'].closed()
        self.assertEqual((V, A), (mpq(F(3*11*13, 16)), mpq(2*(F(33, 4)+F(143, 16)+F(39, 4)))))

    def test_frames(self):
        # The axis frame: `x` and `(1, 2, 2) / 3` orthonormal; the box's axes
        # orthonormal, `x * y` its normal.
        X, A = fixtures.X, [c/3 for c in fixtures.AXIS]
        self.assertEqual((dot(X, X), dot(A, A), dot(X, A)), (1, 1, 0))
        x, y, n = fixtures.BoxBody.AXES
        for u, v in ((x, x), (y, y), (n, n)):
            self.assertEqual(dot(u, v), 1)
        for u, v in ((x, y), (y, n), (n, x)):
            self.assertEqual(dot(u, v), 0)
        self.assertEqual([x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]], list(n))
        # No face of the box parallel or normal to the axis.
        for u in (x, y, n):
            self.assertNotEqual(dot(u, A), 0)
            cross = [u[1]*A[2]-u[2]*A[1], u[2]*A[0]-u[0]*A[2], u[0]*A[1]-u[1]*A[0]]
            self.assertNotEqual(cross, [0, 0, 0])

    def test_model_against_occt(self):
        pin = fixtures.BODIES['pin']
        flats = [(a, b) for a, b in zip(pin.occt, pin.occt[1:]) if fixtures.near(a, b)]
        self.assertEqual(len(flats), 2)
        self.assertEqual(fixtures.BODIES['waist'].occt, fixtures.BODIES['waist'].model)
        wrong = list(pin.model)
        wrong[3] = (F(4)+F(1, 2**20), F(1))
        with self.assertRaises(AssertionError):
            fixtures.Body('wrong', pin.occt, wrong)

    def test_lattice_forms(self):
        for b in TURNED:
            L = b.lattice
            lo, hi = L.w_min-2, L.w_max+2
            for k in range(101):
                w = lo+(hi-lo)*F(k, 100)+F(1, 997)
                self.assertEqual(L.form(w), L.profile(w), (b.name, w))
            self.assertEqual(sorted(set(leaves(L.expr()))), list(range(len(L.parts))))
            # Every reference solid on the tilted frame at a dyadic point on
            # the axis, its heights three times the profile's.
            for spec, p in zip(b.ref_specs, L.parts):
                lo_, hi_ = p['range']
                frame = spec[-1] if spec[0] == 'cone' else spec[2]
                self.assertEqual(frame[3:], fixtures.GEN)
                if spec[0] == 'cone':
                    self.assertEqual(F(spec[3]), 3*(hi_-lo_))
                    self.assertEqual(tuple(F(c) for c in frame[:3]), b.place(0, lo_))
                else:
                    self.assertEqual((F(spec[3]), F(spec[4])), (3*lo_, 3*hi_))

    def test_waist_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['waist']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, A = chain.measures(b.ref_expr)
            want, area = b.closed()
            # The reference's primitives on the kernel's stored axes of the
            # tilted frame: within a rounding of the heights of one line.
            self.assertLess(abs(V-want), mp.mpf(10)**-15*want)
            self.assertLess(abs(A-area), mp.mpf(10)**-15*area)
        finally:
            ref.SCAN = scan

    def test_partners_placed(self):
        # The box about the bodies' point at radius 6 and height 1 (the world's
        # `(11/2, 15/4, -2)`), which lies inside it by a quarter at least.
        box = fixtures.BODIES['skewbox']
        corner = fixtures.BODIES['waist'].place(6, 1)
        self.assertEqual(corner, box.origin)
        x0, y0, x1, y1, h0, h1 = box.extent
        for lo, hi in ((x0, x1), (y0, y1), (h0, h1)):
            self.assertTrue(lo+F(1, 4) <= 0 <= hi-F(1, 4))
        # The ball across the waist's upper cylinder (radius 4, heights 2 to
        # 3): its centre's radius and height off the axis.
        _, r, frame, _, _ = fixtures.WAIST_BALL
        d = [mpq(frame[i])-mpq(fixtures.ORIGIN[i]) for i in range(3)]
        w = dot(d, [mpq(c) for c in fixtures.AXIS])/9
        rho = mp.sqrt(dot(d, d)-9*w*w)
        self.assertTrue(rho < 4 < rho+r)
        self.assertTrue(2 < w-mpq(F(r))/3 and w+mpq(F(r))/3 < 3)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 9)
        self.assertEqual({c.kind for c in listed}, {'solid'})
        for c in listed:
            want = 2 if c.group.endswith('_box') else 1
            self.assertEqual(c.encode().count('\nbrep imported/tilted_'), want, c.name)
            self.assertEqual(c.native().count('\nbrep imported/tilted_'), want, c.name)
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nrevolve '), len(TURNED))
        self.assertEqual(text.count('\nprism '), 1)
        self.assertEqual(text.count('\ntolerance 0.001\n'), len(fixtures.BODIES))


if __name__ == '__main__':
    unittest.main()
