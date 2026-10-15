"""Independent checks of S9e.4b.4d.2c.2's fixtures
(`generate_turned_crossings_boolean_fixtures.py`) and of the consistency
generator's crossings (`generate_turned_consistency_fixtures.py`): the
ribbed body's closed form (Pappus over its model profile) against a direct
quadrature, its arcs tangent to their planes and cylinders, its sharp
lattice's form against its sharp profile, its construction's membership
(the lattice with its pieces) against the profile with arcs at points of the
meridian, the reference's turned input of it (2c.1b's, with its chamfer's
cone) exactly and in binary64, and alone on the chained reference against
its closed form; every partner crossing the tangent circles it is placed
across transversally (each circle's points inside the partner's face and
the face's normal off the circle's tangent there), clear of every other
circle's touching; the declared case (a box side through a circle's
tangent); the case list and its protocol rows; and the consistency
generator's lens of two discs and its crossings' rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_crossings_boolean_fixtures as fixtures
import generate_turned_consistency_fixtures as consistency
import test_turned_fillets_boolean_reference as fillets_test

mpq = fillets_test.mpq
radius = fillets_test.radius
inside = fillets_test.inside

# Each body's tangent circles `(r, w)`: a fillet's or a round's two ends.
D = consistency.D


def circles(body):
    out = []
    pts, segs = body.model.pts, body.model.segs
    for i, seg in enumerate(segs):
        if seg[0] == 'A':
            out += [pts[i], pts[i+1]]
    return out


def triangle_inside(tri, p):
    """Whether `p` lies strictly inside the triangle `tri` (either turn)."""
    s = []
    for k in range(3):
        (x0, y0), (x1, y1) = tri[k], tri[(k+1) % 3]
        s.append((x1-x0)*(p[1]-y0)-(y1-y0)*(p[0]-x0))
    return all(v > 0 for v in s) or all(v < 0 for v in s)


class TurnedCrossingsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_ribbed_closed_form(self):
        b = fixtures.OWN['ribbed']
        m = b.model
        V, A = b.closed()
        heights = sorted({mpq(w) for _, w in m.pts})
        vol = mp.pi*mp.quad(lambda w: radius(m, w)**2, heights)
        self.assertLess(abs(V-vol), mp.mpf(10)**-25*V)
        # The area: each plane's annulus, the cylinders', the chamfer's
        # frustum, the fillet's and the round's quarter bands (`2 pi t (pi R
        # / 2 -+ t)`), exactly but for pi and the chamfer's slant.
        pi = mp.pi
        planes = pi*24**2+pi*(21**2-10**2)+pi*mpq(F(13, 2))**2
        cylinders = 2*pi*24*32+2*pi*mpq(F(19, 2))*mpq(41-F(71, 2))
        chamfer = pi*(24+21)*mp.sqrt(18)
        fillet = 2*pi*mpq(F(1, 2))*(pi*10/2-mpq(F(1, 2)))
        round_ = 2*pi*3*(pi*mpq(F(13, 2))/2+3)
        self.assertLess(abs(A-(planes+cylinders+chamfer+fillet+round_)), mp.mpf(10)**-25*A)

    def test_ribbed_arcs_and_lattice(self):
        b = fixtures.OWN['ribbed']
        m = b.model
        self.assertEqual(len(m.corners), 2)
        fillet, round_ = m.corners
        self.assertEqual((fillet['centre'], fillet['t'], fillet['round']), ((10, F(71, 2)), F(1, 2), False))
        self.assertEqual((round_['centre'], round_['t'], round_['round']), ((F(13, 2), 41), 3, True))
        self.assertEqual(fillet['corner'], (F(19, 2), 35))
        self.assertEqual(round_['corner'], (F(19, 2), 44))
        for c in m.corners:
            # The kite a rectangle of the centre, the tangent points and the
            # corner.
            k = m.kite(c)
            for p in k[1:]:
                self.assertTrue(p[0] in (c['centre'][0], c['corner'][0]) and p[1] in (c['centre'][1],
                                                                                      c['corner'][1]))
        L = b.lattice
        lo, hi = L.w_min-2, L.w_max+2
        for k in range(101):
            w = lo+(hi-lo)*F(k, 100)+F(1, 997)
            self.assertEqual(L.form(w), L.profile(w), w)
        # An arc not tangent to its plane is caught.
        bad = list(fixtures.make_ribbed())
        r, w, seg = bad[4]
        bad[4] = (r, w, fixtures.arc(seg[1], seg[2], seg[3]+F(1, 64), False))
        with self.assertRaises(AssertionError):
            fixtures.Model(bad)

    def test_ribbed_membership(self):
        # The construction (the lattice with its pieces), the profile with
        # arcs and the reference's turned input (2c.1b's, exactly and in
        # binary64) at a grid of meridian points off every surface.
        b = fixtures.OWN['ribbed']
        m = b.model
        inputs = [ref.make_input(fixtures.construction(s, 0)) for s in m.construction_specs(b.frame_at)]
        turned = ref.make_input(fixtures.TurnedCase(b))
        self.assertIsInstance(turned, fixtures.cones.TurnedInput)
        count = 0
        for i in range(-2, 61):
            for j in range(-2, 61):
                r = 25*(mp.mpf(i)+mp.mpf(1)/13)/58
                w = 46*(mp.mpf(j)+mp.mpf(5)/13)/58-1
                if r <= 0:
                    continue
                X = (r, mp.mpf(0), w)
                mem = [bool(I.contains(X)) for I in inputs]
                got = ref.evaluate(m.construction_expr(), mem)
                self.assertEqual(got, inside(m, r, w), (r, w))
                self.assertEqual(turned.contains(X), got, (r, w))
                self.assertEqual(ref.Float(turned).contains([float(c) for c in X]), got, (r, w))
                count += got
        self.assertGreater(count, 500)
        # Near each tube, along its corner's diagonal, just inside and
        # outside its arc.
        for c in m.corners:
            (cr, cw), (pr, pw), t = c['centre'], c['corner'], c['t']
            n = mp.sqrt(mpq((pr-cr)**2+(pw-cw)**2))
            for off in (F(-1, 1024), F(1, 1024)):
                r = mpq(cr)+(mpq(t)+mpq(off))*mpq(pr-cr)/n
                w = mpq(cw)+(mpq(t)+mpq(off))*mpq(pw-cw)/n
                X = (r, mp.mpf(0), w)
                self.assertEqual(turned.contains(X), inside(m, r, w), (c['centre'], off))
                self.assertEqual(ref.Float(turned).contains([float(x) for x in X]), inside(m, r, w))

    def test_ribbed_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            chain = ref.Chain([fixtures.TurnedCase(fixtures.OWN['ribbed'])])
            want, mom, area = chain.inputs[0].closed()
            for k in (0, 1):
                V, m, A = chain.measures(0, k)
                self.assertLess(abs(V-want), mp.mpf(10)**-26*want)
                self.assertLess(abs(A-area), mp.mpf(10)**-26*area)
                for i in range(3):
                    self.assertLess(abs(m[i]-mom[i]), mp.mpf(10)**-26*want*50)
        finally:
            ref.SCAN = scan

    def check_box(self, spec, crossed, others, a):
        """A box on the frame `XY`: its sides `y = +-a` across each circle of
        `crossed` (the circle's points there strictly inside the side, its
        normal `y` off the circle's tangent), and no circle of `others`
        touching a side or a cap."""
        _, (sq,), frame, h0, h1 = spec
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['XY'])
        xs = sorted({F(p[0]) for p in sq.points})
        ys = sorted({F(p[1]) for p in sq.points})
        z0, z1 = F(frame[2])+F(h0), F(frame[2])+F(h1)
        self.assertEqual(ys, [-a, a])
        for r, w in crossed:
            x = mp.sqrt(mpq(r*r-a*a))
            self.assertTrue(mpq(xs[0]) < x < mpq(xs[1]) and z0 < w < z1, (spec, r, w))
            # The tangent `(-a, x) / r` against the normal `(0, 1, 0)`.
            self.assertGreater(x, 0)
        for r, w in others:
            # No cap at a circle's height, no side x = c tangent to it.
            self.assertNotIn(w, (z0, z1))
            self.assertNotIn(r, [abs(x) for x in xs])

    def test_partners_cross(self):
        bodies = fixtures.BODIES
        ribbed = circles(bodies['ribbed'])
        self.assertEqual(ribbed, [(10, 35), (F(19, 2), F(71, 2)), (F(19, 2), 41), (F(13, 2), 44)])
        # The box across the fillet's two circles and the round's lower one;
        # the collar's across its fillet's and its upper round's lower one.
        self.check_box(fixtures.RIBBED_BOX, ribbed[:3], ribbed[3:], 1)
        collar = circles(bodies['collar'])
        self.check_box(fixtures.COLLAR_BOX, collar[2:5], collar[:2]+collar[5:], F(1, 2))
        # The rib: each crossed circle's points on y = +-1 inside its
        # triangle, its hypotenuse past the round's lower circle.
        for spec, crossed, a in ((fixtures.RIBBED_RIB, ribbed[:3], 1),
                                 (fixtures.ROUND_RIB, circles(bodies['round']), F(1, 2))):
            _, (tri,), frame, y0, y1 = spec
            self.assertEqual((F(y0), F(y1)), (-a, a))
            pts = [(mpq(F(u)), mpq(F(v))) for u, v in tri.points]
            for r, w in crossed:
                x = mp.sqrt(mpq(r)**2-mpq(a)**2)
                self.assertTrue(triangle_inside(pts, (mpq(w), x)), (r, w))
        self.assertEqual(sorted(F(v) for _, v in fixtures.RIBBED_RIB[1][0].points), [8, 8, 20])
        self.assertEqual(max(F(u) for u, _ in fixtures.RIBBED_RIB[1][0].points), 42)
        # The post: radii 9 to 41/4 from the axis, its heights 33 to 47,
        # across the fillet's two circles and the round's lower one, clear
        # of the round's upper one and of every tube's top and bottom.
        _, (disc,), frame, h0, h1 = fixtures.RIBBED_POST
        cx, cy, rho = (F(c) for c in disc.circle)
        z0, z1 = F(frame[2])+F(h0), F(frame[2])+F(h1)
        for r, w in ribbed[:3]:
            self.assertTrue(cx-rho < r < cx+rho and z0 < w < z1)
        self.assertFalse(cx-rho < ribbed[3][0] < cx+rho)
        for c in bodies['ribbed'].model.corners:
            for w in (c['centre'][1]-c['t'], c['centre'][1]+c['t']):
                self.assertNotIn(w, (z0, z1))
        # The rods along x: about (0, zc) of radius 1/4, across a circle at
        # height w (|w - zc| < 1/4).
        for spec, (r, w) in ((fixtures.STEP_RING, circles(bodies['step'])[2]),
                             (fixtures.FILLET_ROD, circles(bodies['fillet'])[1])):
            _, (disc,), frame, h0, h1 = spec
            self.assertEqual(frame[3:], fixtures.forms.FRAMES['SIDE'])
            yc, zc, rho = (F(c) for c in disc.circle)
            self.assertEqual(yc, 0)
            self.assertLess(abs(w-zc), rho)
            self.assertTrue(F(frame[0])+F(h0) < -r and F(frame[0])+F(h1) > r)
        # The skew cube: about the round's lower circle's point, no axis of
        # its frame parallel or normal to the body's, its half-diagonal
        # clear of every other circle.
        _, (sq,), frame, h0, h1 = fixtures.RIBBED_SKEW
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['SKEW'])
        self.assertEqual(tuple(F(c) for c in frame[:3]), (F(19, 2), 0, 41))
        o, x, y, n = fixtures.stored_axes(frame)
        for axis in (x, y, n):
            self.assertTrue(all(abs(c) > 0.1 for c in axis), axis)
        for r, w in ribbed[:2]+ribbed[3:]:
            self.assertGreater(abs(w-41), 3**0.5)
        self.assertEqual(circles(bodies['step'])[2], (4, 7))
        self.assertEqual(circles(bodies['fillet'])[1], (F(70, 17), F(43, 17)))

    def test_declared(self):
        # The box's side x = 9/2 through the step body's fillet's circle on
        # its step (of radius 9/2 at 3) at (9/2, 0, 3), where the circle's
        # tangent is (0, 1, 0), in the side; no cap at a tube's top or
        # bottom.
        _, (sq,), frame, h0, h1 = fixtures.STEP_TOUCH
        xs = sorted({F(p[0]) for p in sq.points})
        ys = sorted({F(p[1]) for p in sq.points})
        self.assertEqual(xs[1], F(9, 2))
        self.assertEqual(circles(fixtures.BODIES['step'])[0], (F(9, 2), 3))
        self.assertTrue(ys[0] < 0 < ys[1] and F(h0) < 3 < F(h1))
        for c in fixtures.BODIES['step'].model.corners:
            for w in (c['centre'][1]-c['t'], c['centre'][1]+c['t']):
                self.assertNotIn(w, (F(h0), F(h1)))
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.kind == 'degenerate'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 27)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 24, 'empty': 0, 'degenerate': 3})
        for c in listed:
            b = fixtures.BODIES[c.items[0][1]]
            self.assertEqual(c.encode().count(f'\nbrep {b.path}\n'), 1, c.name)
            self.assertEqual(c.native().count(f'\nbrep {b.path}\n'), 1, c.name)
            self.assertEqual(c.expr(), (c.last, b.ref_expr, len(b.ref_specs)))
        # This step's own body written once, turned a whole turn, at the
        # resolution 3e-3; the others the earlier steps' files.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), 1)
        self.assertEqual(text.count('\nrevolve '), 1)
        self.assertEqual(text.count('\ntolerance 0.003\n'), 1)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        import generate_turned_fillets_boolean_fixtures as fillets
        import generate_turned_cone_fillets_boolean_fixtures as cones
        earlier = {b.path for b in fillets.BODIES.values()} | {b.path for b in cones.BODIES.values()}
        self.assertFalse({b.path for b in fixtures.OWN.values()} & earlier)
        self.assertTrue({fixtures.BODIES[k].path for k in ('step', 'collar', 'round', 'fillet')} <= earlier)

    def test_consistency_crossings(self):
        # The lens of two discs: a disc inside the other, apart, and a
        # crossing against the chords' widths.
        lens = consistency.lens
        m = consistency.m
        self.assertEqual(lens(m(1), m(1), m(10)), mp.pi)
        self.assertEqual(lens(m(20), m(1), m(10)), 0)
        d, rho, r = m(F(77, 8)), m(F(5, 8)), m(F(19, 2))
        widths = mp.quad(lambda x: 2*min(mp.sqrt(max(m(0), rho*rho-(x-d)**2)),
                                         mp.sqrt(max(m(0), r*r-x*x))),
                         [d-rho, (r*r-rho*rho+d*d)/(2*d), r])
        self.assertLess(abs(widths-lens(d, rho, r)), mp.mpf(10)**-25)
        # The crossings' rows: the known failure's rib first, each partner
        # across the fillet and the round of the consistency body.
        rows = consistency.crossings_text().splitlines()
        self.assertTrue(rows[1].startswith('body\tribbed imported/consistency_ribbed.brep volume '))
        self.assertTrue(rows[2].startswith('case\trib_high rib 8 34.5 20 42 -1 1 volume 90.0 common '))
        self.assertTrue(rows[3].startswith('case\tbox_both box 8 -1 34 3 2 9 volume 54.0 common '))
        self.assertTrue(rows[4].startswith('case\tpost post 9.625 0 0.625 33 47 volume '))
        self.assertEqual(len(rows), 5)


if __name__ == '__main__':
    unittest.main()
