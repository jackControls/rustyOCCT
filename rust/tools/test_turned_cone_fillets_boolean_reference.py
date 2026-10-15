"""Independent checks of S9e.4b.4d.2c.1b's fixtures
(`generate_turned_cone_fillets_boolean_fixtures.py`): the bodies' closed
forms (Pappus over the model's profile with arcs) against direct quadratures
and, for two bodies, against their sharp profiles less their kites' corner
regions, each arc tangent to its cone and to its plane or cylinder along a
direction of dyadic half-angle tangent (a wrong one caught), the refitted
model against OCCT's profile, the sharp lattice's form against the sharp
profile, the construction's membership (the lattice with its pieces' slant
boxes less their tubes) against the profile with arcs at points of the
meridian, two bodies alone on the chained reference against their closed
forms, the partners' places (clear of every fillet's and round's circle),
the declared case (a ball touching the round) and the case list and its
protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_cone_fillets_boolean_fixtures as fixtures


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def radius(model, w):
    """The model profile's radius at height `w`: its rising lines' and its
    arcs' (the outer of the two points at a step)."""
    best = mp.mpf(0)
    pts, segs = model.pts, model.segs
    for i in range(len(pts)-1):
        (ra, wa), (rb, wb) = pts[i], pts[i+1]
        lo, hi = sorted((mpq(wa), mpq(wb)))
        if not lo <= w <= hi or wa == wb:
            continue
        if segs[i][0] == 'A':
            _, cr, cw, t, ccw = segs[i]
            d = mp.sqrt(max(mp.mpf(0), mpq(t)**2-(w-mpq(cw))**2))
            # The arc's side of its centre: its ends' radii.
            side = 1 if max(ra, rb) > cr else -1
            best = max(best, mpq(cr)+side*d)
        else:
            best = max(best, mpq(ra)+(mpq(rb)-mpq(ra))*(w-mpq(wa))/(mpq(wb)-mpq(wa)))
    return best


def polygon_moment(pts):
    """A polygon's `int r dA` in the meridian (Green's theorem, exact):
    `sum (r_i^2 + r_i r_j + r_j^2) (w_j - w_i) / 6` counter-clockwise."""
    out = F(0)
    for (ra, wa), (rb, wb) in zip(pts, pts[1:]+pts[:1]):
        out += (ra*ra+ra*rb+rb*rb)*(wb-wa)/6
    return out


class TurnedConeFilletsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        tight = mp.mpf(10)**-25
        for b in fixtures.BODIES.values():
            m = b.model
            V, A = b.closed()
            heights = sorted({mpq(w) for _, w in m.pts})
            vol = mp.pi*mp.quad(lambda w: radius(m, w)**2, heights)
            area = mp.mpf(0)
            for i in range(len(m.pts)-1):
                (ra, wa), (rb, wb) = [(mpq(a), mpq(c)) for a, c in m.pts[i:i+2]]
                if m.segs[i][0] == 'A':
                    _, cr, cw, t, _ = m.segs[i]
                    th0 = mp.atan2(wa-mpq(cw), ra-mpq(cr))
                    th1 = mp.atan2(wb-mpq(cw), rb-mpq(cr))
                    if th1-th0 > mp.pi:
                        th1 -= 2*mp.pi
                    if th0-th1 > mp.pi:
                        th1 += 2*mp.pi
                    area += 2*mp.pi*mpq(t)*abs(mp.quad(lambda th: mpq(cr)+mpq(t)*mp.cos(th), [th0, th1]))
                else:
                    length = mp.sqrt((rb-ra)**2+(wb-wa)**2)
                    area += 2*mp.pi*mp.quad(lambda s: (ra+(rb-ra)*s)*length, [0, 1])
            self.assertLess(abs(V-vol), tight*V, b.name)
            self.assertLess(abs(A-area), tight*A, b.name)
        # The round's and the top's volumes: `2 pi` times the sharp
        # profile's `int r dA` less the round's corner region's (its kite's,
        # exactly, less its sector's: `R t^2 a / 2 + t^3 (sin p1 - sin p0) /
        # 3` over the directions `p0` to `p1` of its ends from its centre).
        for name in ('round', 'top'):
            b = fixtures.BODIES[name]
            m = b.model
            sharp = polygon_moment(m.sharp)
            (c,) = m.corners
            (R, W), t = c['centre'], c['t']
            kite = abs(polygon_moment(m.kite(c)))
            p = sorted(mp.atan2(mpq(e[1]), mpq(e[0])) for e in c['dirs'])
            sector = mpq(R*t*t)*(p[1]-p[0])/2+mpq(t**3)*(mp.sin(p[1])-mp.sin(p[0]))/3
            V, _ = b.closed()
            want = 2*mp.pi*(mpq(sharp)-(mpq(kite)-sector))
            self.assertLess(abs(V-want), tight*V, name)

    def test_arcs_tangent(self):
        for b in fixtures.BODIES.values():
            (c,) = b.model.corners
            (cr, cw), t = c['centre'], c['t']
            for (er, ew), end in zip(c['dirs'], c['ends']):
                # The tube's points its centre's feet, a unit apart.
                self.assertEqual(er*er+ew*ew, 1)
                self.assertEqual((cr+t*er, cw+t*ew), end)
            # The direction to the cone of a half-angle tangent on the
            # multiples of 2^-40 (dyadic but for the refit's).
            s = fixtures.half_tangent(c['dirs'][c['cone']])
            self.assertEqual((s*2**40).denominator, 1, b.name)
            self.assertEqual(fixtures.direction(s), c['dirs'][c['cone']])
        # An arc not tangent to its cone is caught.
        good = fixtures.BODIES['round'].model.profile
        bad = list(good)
        r, w, seg = bad[2]
        bad[2] = (r, w, fixtures.arc(seg[1], seg[2], seg[3]+F(1, 64), True))
        with self.assertRaises(AssertionError):
            fixtures.Model(bad)
        # A round between a plane and a cylinder is the earlier part's.
        plain = [(0, 0, fixtures.L), (6, 0, fixtures.L), (6, 4, fixtures.arc(4, 4, 2, True)), (4, 6, fixtures.L),
                 (0, 6, fixtures.L)]
        with self.assertRaises(AssertionError):
            fixtures.Model(plain)

    def test_refit_against_occt(self):
        D = fixtures.D
        b = fixtures.BODIES['refit']
        # OCCT's round: about `(4, 4 + D)` through the cylinder's top `(6,
        # 4)`, its radius 2 within 3e-13; the model's tangent at `4 + D`.
        (_, _, seg) = b.occt[2]
        self.assertEqual(seg[2], 4+D)
        rho = mp.mpf(seg[3])
        self.assertLess(abs(rho-mp.sqrt(mpq(4+D*D))), mp.mpf(2)**-52*2)
        self.assertLess(abs(rho-2), mp.mpf('3e-13'))
        self.assertEqual(b.model.profile[2][1], 4+D)
        # The model's cone within 1e-11 of OCCT's at its ends.
        for k in (3, 4):
            self.assertLess(abs(mpq(b.model.profile[k][0])-mpq(b.occt[k][0])), mp.mpf('1e-11'))
        # OCCT's cone no Pythagorean direction: its direction's half-angle
        # tangent no dyadic number.
        (t1, q_) = (b.occt[3][:2], b.occt[4][:2])
        self.assertNotEqual(fixtures.kind_of(t1, q_), 'cylinder')
        # A model off OCCT's by more than 4 D is caught.
        occt, model = fixtures.refit()
        model[2] = (6, 4+5*D, fixtures.arc(4, 4+5*D, 2, True))
        with self.assertRaises(AssertionError):
            fixtures.Body('wrong', occt, model, tolerance='0.001')

    def test_lattice_forms(self):
        for b in fixtures.BODIES.values():
            L = b.lattice
            lo, hi = L.w_min-2, L.w_max+2
            for k in range(101):
                w = lo+(hi-lo)*F(k, 100)+F(1, 997)
                self.assertEqual(L.form(w), L.profile(w), (b.name, w))

    def test_construction_membership(self):
        # The lattice with its pieces against the profile with arcs at a
        # grid of meridian points (none on a surface: thirteenths off the
        # dyadic grid), exactly; the reference's turned input alike, exactly
        # and in binary64.
        for b in fixtures.BODIES.values():
            m = b.model
            turned = ref.make_input(fixtures.TurnedCase(b))
            r_max = max(r for r, _ in m.sharp)
            w_max = max(w for _, w in m.sharp)
            count = 0
            for i in range(-2, 61):
                for j in range(-2, 61):
                    r = (r_max+1)*(F(i)+F(1, 13))/58
                    w = (w_max+2)*(F(j)+F(5, 13))/58-1
                    if r <= 0:
                        continue
                    got = m.construction_inside(r, w)
                    self.assertEqual(got, m.inside(mpq(r), mpq(w)), (b.name, r, w))
                    X = (mpq(r), mp.mpf(0), mpq(w))
                    self.assertEqual(turned.contains(X), got, (b.name, r, w))
                    self.assertEqual(ref.Float(turned).contains([float(c) for c in X]), got, (b.name, r, w))
                    count += got
            self.assertGreater(count, 500, b.name)
            # Close to each tube: points along its directions to its ends
            # and between them, just inside and outside the arc.
            (c,) = m.corners
            (cr, cw), t = c['centre'], c['t']
            e0, e1 = c['dirs']
            for k in range(1, 8):
                d = (e0[0]*(8-k)+e1[0]*k, e0[1]*(8-k)+e1[1]*k)
                n = mp.sqrt(mpq(d[0]*d[0]+d[1]*d[1]))
                for off in (F(-1, 1024), F(1, 1024)):
                    r = mpq(cr)+(mpq(t)+mpq(off))*mpq(d[0])/n
                    w = mpq(cw)+(mpq(t)+mpq(off))*mpq(d[1])/n
                    rr, ww = F(str(r)), F(str(w))
                    self.assertEqual(m.construction_inside(rr, ww), m.inside(mpq(rr), mpq(ww)), (b.name, k, off))

    def test_turned_inputs_on_the_chained_reference(self):
        # Two bodies as the reference's turned inputs: both families'
        # measures their closed forms' (volume, moments, area).
        scan, ref.SCAN = ref.SCAN, 24
        try:
            for name in ('round', 'top'):
                chain = ref.Chain([fixtures.TurnedCase(fixtures.BODIES[name])])
                want, mom, area = chain.inputs[0].closed()
                for k in (0, 1):
                    V, m, A = chain.measures(0, k)
                    self.assertLess(abs(V-want), mp.mpf(10)**-28*want, name)
                    self.assertLess(abs(A-area), mp.mpf(10)**-28*area, name)
                    for i in range(3):
                        self.assertLess(abs(m[i]-mom[i]), mp.mpf(10)**-28*want*10, name)
        finally:
            ref.SCAN = scan

    def test_partners_placed(self):
        # The rod about `(y, z) = (1/2, 9/2)` of radius 1/4: between the
        # round's tangent circles' heights 4 and 84/17, across the body.
        _, (c,), frame, h0, h1 = fixtures.ROUND_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['SIDE'])
        self.assertTrue(4 < cy-r and cy+r < F(84, 17))
        self.assertTrue(frame[0]+h0 < -6 and frame[0]+h1 > 6)
        # The ball in the fillet: across its tube, clear of the step's plane
        # (at 2) and of the boss's cone (its line through `(70/17, 43/17)`
        # normal to `(-15/17, -8/17)`).
        _, r, frame, _, _ = fixtures.FILLET_BALL
        x, z = F(frame[0]), F(frame[2])
        r = F(r)
        self.assertGreater(z-r, 2)
        e = (F(-15, 17), F(-8, 17))
        dist = abs(e[0]*(x-F(70, 17))+e[1]*(z-F(43, 17)))
        self.assertGreater(dist, r)
        d = mp.sqrt(mpq((x-5)**2+(z-3)**2))
        self.assertTrue(d-mpq(r) < 1 < d+mpq(r))
        # The boxes, at their eight corners: the top's across its bottom
        # corner clear of its round's circles (heights 21/5 and 6), the
        # refit's across its top corner clear of its round's (heights 4 + D
        # and its point on the cone, about 4.89).
        refit_cone = fixtures.BODIES['refit'].model.profile[3][1]
        for spec, ws in ((fixtures.TOP_BOX, (F(21, 5), F(6))), (fixtures.REFIT_BOX, (4+fixtures.D, refit_cone))):
            _, (sq,), frame, h0, h1 = spec
            o, x, y, n = [tuple(mpq(v) for v in a) for a in fixtures.stored_axes(frame)]
            corners = []
            for u in sorted({F(p[0]) for p in sq.points}):
                for v in sorted({F(p[1]) for p in sq.points}):
                    for h in (h0, h1):
                        corners.append(tuple(o[k]+mpq(u)*x[k]+mpq(v)*y[k]+mpq(h)*n[k] for k in range(3)))
            heights = [p[2] for p in corners]
            for ww in ws:
                self.assertFalse(min(heights) < mpq(ww) < max(heights), spec)

    def test_declared(self):
        # The touching ball: its centre 25/8 from the round's centre `(4,
        # 4)` along `(24/25, 7/25)` (within the round's turn from `(1, 0)`
        # to `(15/17, 8/17)`), its radius 9/8.
        _, r, frame, _, _ = fixtures.ROUND_TOUCH
        cx, cz = F(frame[0]), F(frame[2])
        self.assertEqual((cx-4, cz-4), (F(25, 8)*F(24, 25), F(25, 8)*F(7, 25)))
        self.assertEqual(F(r), F(9, 8))
        self.assertLess(F(7, 25)/F(24, 25), F(8, 15))
        # Clear of the cylinder (radius 6 to height 4) and of the cone past
        # its tangent point `(98/17, 84/17)`.
        self.assertGreater((cx-6)**2+(cz-4)**2, F(r)**2)
        self.assertGreater((cx-F(98, 17))**2+(cz-F(84, 17))**2, F(r)**2)
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'round_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 15)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 12, 'empty': 0, 'degenerate': 3})
        for c in listed:
            self.assertEqual(c.encode().count('\nbrep imported/cone_fillet_'), 1, c.name)
            self.assertEqual(c.native().count('\nbrep imported/cone_fillet_'), 1, c.name)
            b = fixtures.BODIES[c.items[0][1]]
            self.assertEqual(c.expr(), (c.last, b.ref_expr, len(b.ref_specs)))
        # Each body written once, its profile turned a whole turn, the
        # refitted one at the resolution 1e-3.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nrevolve '), len(fixtures.BODIES))
        self.assertEqual(text.count('\ntolerance 0.001\n'), 1)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of the earlier turned bodies' files.
        import write_turned_bodies as w
        import generate_turned_fillets_boolean_fixtures as fillets
        paths = {b.path for b in fixtures.BODIES.values()}
        self.assertFalse(paths & {f'imported/{n}.brep' for n in w.BODIES})
        self.assertFalse(paths & {b.path for b in fillets.BODIES.values()})


if __name__ == '__main__':
    unittest.main()
