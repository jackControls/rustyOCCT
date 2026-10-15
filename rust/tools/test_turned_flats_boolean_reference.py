"""Independent checks of S9e.4b.4d.2b's fixtures
(`generate_turned_flats_boolean_fixtures.py`): the bodies' closed forms
(Pappus over the model's profile) against direct quadratures, the model's
profile against OCCT's (each near-flat segment level at its wider end's
height, a wrong model caught), the near-flat segments within the resolution
of flat and every other cone far from it, the lattice's form against the
model's profile at heights of its own, the pin alone on the chained
reference against its closed form, the partners' places, the declared case
(a box resting within the resolution of the pin's step) and the case list
and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_flats_boolean_fixtures as fixtures


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


LATTICED = [b for b in fixtures.BODIES.values() if b.model]


class TurnedFlatsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        for b in LATTICED:
            pts = b.model
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
        # The pin's volume `pi (36 * 4 + 8 (16 + 12 + 9) / 3)`, the collar's
        # `pi (9 (6 - 2 D) + 36 (4 + D) + 2 (36 + 24 + 16) / 3)`.
        D = fixtures.D
        V, _ = fixtures.BODIES['pin'].closed()
        self.assertLess(abs(V-mp.pi*(144+mpq(F(8*37, 3)))), tight)
        V, _ = fixtures.BODIES['collar'].closed()
        self.assertLess(abs(V-mp.pi*mpq(9*(6-2*D)+36*(4+D)+F(2*76, 3))), tight)

    def test_model_against_occt(self):
        # Each near-flat segment rises or falls by `D` over 2 to 6: within
        # 1e-3 rad of flat and within a hundredth of the resolution, its
        # model level at its wider end; every other cone at least 0.1 rad
        # from flat and from a cylinder.
        res = mpq(F(fixtures.RESOLUTION))
        flats = 0
        for b in fixtures.BODIES.values():
            for (ra, wa), (rb, wb) in zip(b.occt, b.occt[1:]):
                if wa == wb or ra == rb:
                    continue
                angle = mp.atan(abs(mpq(wb-wa))/abs(mpq(rb-ra)))
                if fixtures.near((ra, wa), (rb, wb)):
                    flats += 1
                    self.assertLess(angle, mp.mpf('1e-5'))
                    self.assertLessEqual(abs(mpq(wb-wa)), res/100)
                else:
                    self.assertGreater(min(angle, mp.pi/2-angle), mp.mpf('0.1'), b.name)
        self.assertEqual(flats, 5)
        # A model whose step is level at its narrower end is caught.
        pin = fixtures.BODIES['pin']
        wrong = list(pin.model)
        wrong[2:4] = [(F(6), 4+fixtures.D), (4-fixtures.D/8, 4+fixtures.D)]
        with self.assertRaises(AssertionError):
            fixtures.Body('wrong', pin.occt, wrong)
        # A model whose neighbour leaves its OCCT line is caught.
        wrong = list(pin.model)
        wrong[3] = (F(4)+F(1, 2**20), F(4))
        with self.assertRaises(AssertionError):
            fixtures.Body('wrong', pin.occt, wrong)

    def test_lattice_forms(self):
        for b in LATTICED:
            L = b.lattice
            lo, hi = L.w_min-2, L.w_max+2
            for k in range(101):
                w = lo+(hi-lo)*F(k, 100)+F(1, 997)
                self.assertEqual(L.form(w), L.profile(w), (b.name, w))
            self.assertEqual(sorted(set(leaves(L.expr()))), list(range(len(L.parts))))
        # The pin: its cylinder fused with its draft cone, reaching from its
        # step's height 4 down into the cylinder by a quarter of it.
        L = fixtures.BODIES['pin'].lattice
        self.assertEqual(L.terms, [[0], [1]])
        self.assertEqual(L.parts[1]['range'], (F(3), F(12)))

    def test_pin_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            b = fixtures.BODIES['pin']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.ref_specs])
            V, _, A = chain.measures(b.ref_expr)
            want, area = b.closed()
            self.assertLess(abs(V-want), mp.mpf(10)**-28*want)
            self.assertLess(abs(A-area), mp.mpf(10)**-28*area)
        finally:
            ref.SCAN = scan

    def test_partners_placed(self):
        # The pin's ball about `(1/2, 1/4, -1)` of radius 2: across the bottom
        # plane in a circle of radius `sqrt(3)` about the dimple's apex, the
        # origin inside it, clear of the cylinder of radius 6 and the step.
        _, r, frame, _, _ = fixtures.PIN_BALL
        d = mp.sqrt(mpq(frame[0])**2+mpq(frame[1])**2)
        self.assertLess(d, mp.sqrt(r*r-mpq(frame[2])**2))
        self.assertLess(d+r, 6)
        self.assertLess(mpq(frame[2])+r, 4)
        # The collar's rod about `(y, z) = (3/4, 21/4)` of radius 5/4: across
        # the step out's plane `6 - D`, clear of the bottom and the chamfer.
        _, (c,), frame, h0, h1 = fixtures.COLLAR_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['SIDE'])
        self.assertTrue(cy-r < 6-fixtures.D < cy+r and cy-r > 1 and cy+r < 10)
        self.assertTrue(frame[0]+h0 < -6 and frame[0]+h1 > 6)

    def test_declared(self):
        # The resting box: its floor 2^-12 above the pin's step (within the
        # resolution), over the annulus between radii 4 and 6, clear of the
        # draft cone (radius at most 4 above the step).
        _, (bx,), frame, h0, h1 = fixtures.PIN_REST
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['XY'])
        gap = mpq(frame[2])+h0-4
        self.assertTrue(0 < gap <= mpq(F(fixtures.RESOLUTION))/4)
        xs = sorted({F(p[0]) for p in bx.points})
        self.assertTrue(4 < xs[0] < 6 < xs[-1])
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'pin_rest'}, {fixtures.REST})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 12)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 9, 'empty': 0, 'degenerate': 3})
        for c in listed:
            self.assertEqual(c.encode().count('\nbrep imported/flat_'), 1, c.name)
            self.assertEqual(c.native().count('\nbrep imported/flat_'), 1, c.name)
            self.assertEqual(c.expr(), (c.last, fixtures.BODIES[c.items[0][1]].lattice.expr(), len(
                fixtures.BODIES[c.items[0][1]].ref_specs)))
        # Each body written once, its profile turned a whole turn at the
        # resolution 1e-3.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nrevolve '), len(fixtures.BODIES))
        self.assertEqual(text.count('\ntolerance 0.001\n'), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, {'pin', 'collar'})
        # None of the earlier turned bodies' files.
        import write_turned_bodies as w
        import generate_turned_profiles_boolean_fixtures as profiles
        paths = {b.path for b in fixtures.BODIES.values()}
        self.assertFalse(paths & {f'imported/{n}.brep' for n in w.BODIES})
        self.assertFalse(paths & {b.path for b in profiles.BODIES.values()})


if __name__ == '__main__':
    unittest.main()
