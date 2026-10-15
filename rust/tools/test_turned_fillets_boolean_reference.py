"""Independent checks of S9e.4b.4d.2c.1a's fixtures
(`generate_turned_fillets_boolean_fixtures.py`): the bodies' closed forms
(Pappus over the model's profile with arcs) against direct quadratures, each
arc tangent to its plane and its cylinder (a wrong one caught), the refitted
model against OCCT's profile, the sharp lattice's form against the sharp
profile, the construction's membership (the lattice with its pieces) against
the profile with arcs at points of the meridian, the step alone on the
chained reference against its closed form, the partners' places (clear of
every fillet's and round's circle), the declared case (a ball touching the
round) and the case list and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_fillets_boolean_fixtures as fixtures


def leaves(expr):
    if isinstance(expr, int):
        return [expr]
    return leaves(expr[1])+leaves(expr[2])


def mpq(x):
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def radius(model, w):
    """The model profile's radius at height `w`: its upright lines' and its
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


def inside(model, r, w):
    """Whether a meridian point lies inside the profile with arcs: inside
    the sharp polygon, less each round's corner beyond its arc, with each
    fillet's corner beyond its arc."""
    pts = [(mpq(a), mpq(b)) for a, b in model.sharp]
    wind = False
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]+pts[:1]):
        if (y0 > w) != (y1 > w) and r < x0+(w-y0)*(x1-x0)/(y1-y0):
            wind = not wind
    for c in model.corners:
        (cr, cw), (pr, pw), t = c['centre'], c['corner'], c['t']
        in_kite = mpq(min(cr, pr)) < r < mpq(max(cr, pr)) and mpq(min(cw, pw)) < w < mpq(max(cw, pw))
        beyond = (r-mpq(cr))**2+(w-mpq(cw))**2 > mpq(t)**2
        if in_kite and beyond:
            wind = not c['round']
    return wind


class TurnedFilletsReferenceTests(unittest.TestCase):
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
                    # A quarter turn either way: the short way.
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
        # The step's volume: its sharp profile's `pi (64 * 3 + 16 * 21 / 4)`
        # with its fillet's corner (the square of side 1/2 less a quarter
        # disc, its `int r dA` exactly) and less its round's.
        V, _ = fixtures.BODIES['step'].closed()

        def corner(cr, cw, t, sr):
            # int r dA over the square from the centre toward the corner
            # (side `t`, toward `sr` in r) less the quarter disc.
            square = t*t*(cr+sr*t/2)
            quarter = cr*mp.pi*t*t/4+sr*t**3/3
            return square-quarter
        sharp = mp.pi*(192+84)
        fillet = 2*mp.pi*corner(mp.mpf(9)/2, None, mp.mpf(1)/2, -1)
        rounded = 2*mp.pi*corner(mp.mpf(11)/4, None, mp.mpf(5)/4, 1)
        self.assertLess(abs(V-(sharp+fillet-rounded)), tight*V)

    def test_arcs_tangent(self):
        for b in fixtures.BODIES.values():
            for c in b.model.corners:
                (cr, cw), (pr, pw), t = c['centre'], c['corner'], c['t']
                # Its centre the tube's radius from both lines.
                self.assertEqual(abs(cr-pr), t)
                self.assertEqual(abs(cw-pw), t)
            self.assertEqual(len(b.model.corners), 3 if b.name == 'collar' else 2)
        # An arc not tangent to its cylinder is caught.
        step = list(fixtures.BODIES['step'].occt)
        step[3] = (F(9, 2), 3, fixtures.arc(F(9, 2), F(7, 2), F(1, 2)+F(1, 64), False))
        with self.assertRaises(AssertionError):
            fixtures.Model(step)
        # A round beside a cone is the next part's.
        cone = list(step)
        cone[2] = (8, F(5, 2), fixtures.L)
        with self.assertRaises(AssertionError):
            fixtures.Model(cone)

    def test_refit_against_occt(self):
        D = fixtures.D
        b = fixtures.BODIES['refit']
        # OCCT's round: about `(11/4, 7 + D)` through the shaft's top `(4,
        # 7)`, its radius 5/4 within 4e-13, the model's tangent at `7 + D`.
        (_, _, seg) = b.occt[5]
        self.assertEqual(seg[2], 7+D)
        rho = mp.mpf(seg[3])
        self.assertLess(abs(rho-mp.sqrt(mpq(F(5, 4)**2+D*D))), mp.mpf(2)**-52*2)
        self.assertLess(abs(rho-mp.mpf(5)/4), mp.mpf('4e-13'))
        self.assertEqual(b.model.profile[5][1], 7+D)
        self.assertLessEqual(D, F(1, 1000)/1000)
        # A model off OCCT's by more than 4 D is caught.
        occt, model = fixtures.refit()
        model[5] = (4, 7+5*D, fixtures.arc(F(11, 4), 7+5*D, F(5, 4), True))
        model[6] = (F(11, 4), F(33, 4)+5*D, fixtures.L)
        model[7] = (0, F(33, 4)+5*D, fixtures.L)
        with self.assertRaises(AssertionError):
            fixtures.Body('wrong', occt, model, tolerance='0.001')

    def test_lattice_forms(self):
        for b in fixtures.BODIES.values():
            L = b.lattice
            lo, hi = L.w_min-2, L.w_max+2
            for k in range(101):
                w = lo+(hi-lo)*F(k, 100)+F(1, 997)
                self.assertEqual(L.form(w), L.profile(w), (b.name, w))
            self.assertEqual(sorted(set(leaves(L.expr()))), list(range(len(L.parts))))
            # The construction: every piece's three solids after the
            # lattice's, once each.
            specs = b.model.construction_specs(b.frame_at)
            self.assertEqual(sorted(leaves(b.model.construction_expr())), list(range(len(specs))))

    def test_construction_membership(self):
        # The lattice with its pieces against the profile with arcs at a
        # grid of meridian points (none on a surface: thirteenths off the
        # dyadic grid).
        for b in fixtures.BODIES.values():
            inputs = [ref.make_input(fixtures.construction(s, 0)) for s in b.model.construction_specs(b.frame_at)]
            turned = ref.make_input(fixtures.TurnedCase(b))
            m = b.model
            r_max = max(mpq(r) for r, _ in m.sharp)
            w_max = max(mpq(w) for _, w in m.sharp)
            count = 0
            for i in range(-2, 61):
                for j in range(-2, 61):
                    r = (r_max+1)*(mp.mpf(i)+mp.mpf(1)/13)/58
                    w = (w_max+2)*(mp.mpf(j)+mp.mpf(5)/13)/58-1
                    if r <= 0:
                        continue
                    X = (r, mp.mpf(0), w)
                    mem = [bool(I.contains(X)) for I in inputs]
                    got = ref.evaluate(b.model.construction_expr(), mem)
                    self.assertEqual(got, inside(m, r, w), (b.name, r, w))
                    # The reference's turned input alike, exactly and in
                    # binary64.
                    self.assertEqual(turned.contains(X), got, (b.name, r, w))
                    self.assertEqual(ref.Float(turned).contains([float(c) for c in X]), got, (b.name, r, w))
                    count += got
            self.assertGreater(count, 500, b.name)

    def test_turned_inputs_on_the_chained_reference(self):
        # Each body as the reference's turned input: both families' measures
        # its closed form's (volume, moments, area).
        scan, ref.SCAN = ref.SCAN, 24
        try:
            for name in ('step', 'collar'):
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
        # The rod about `(y, z) = (1/2, 61/8)` of radius 1/4: between the
        # round's tangent circles' heights 7 and 33/4, across the body.
        _, (c,), frame, h0, h1 = fixtures.STEP_ROD
        cx, cy, r = c.circle
        self.assertEqual(frame[3:], fixtures.forms.FRAMES['SIDE'])
        self.assertTrue(7 < cy-r and cy+r < F(33, 4))
        self.assertTrue(frame[0]+h0 < -8 and frame[0]+h1 > 8)
        # The ball in the fillet: a quarter off the step's plane and the
        # shaft, across the tube.
        _, r, frame, _, _ = fixtures.STEP_BALL
        x, z = mpq(frame[0]), mpq(frame[2])
        self.assertTrue(z-3 > r and x-4 > r)
        d = mp.sqrt((x-mpq(F(9, 2)))**2+(z-mpq(F(7, 2)))**2)
        self.assertTrue(d-r < mpq(F(1, 2)) < d+r)
        # The boxes, at their eight corners: the collar's across its sharp
        # corner at radius 6 and height 5, clear of its fillet's circle
        # (radius 7/2 at 5) and its round's (radius 6 at 7).
        for spec, (rs, ws) in ((fixtures.COLLAR_BOX, ((F(7, 2), F(6)), (F(5), F(7)))),
                               (fixtures.REFIT_BOX, ((F(15, 4),), (F(3),)))):
            _, (sq,), frame, h0, h1 = spec
            o, x, y, n = [tuple(mpq(v) for v in a) for a in fixtures.stored_axes(frame)]
            corners = []
            for u in sorted({F(p[0]) for p in sq.points}):
                for v in sorted({F(p[1]) for p in sq.points}):
                    for h in (h0, h1):
                        corners.append(tuple(o[k]+mpq(u)*x[k]+mpq(v)*y[k]+mpq(h)*n[k] for k in range(3)))
            radii = [mp.sqrt(p[0]**2+p[1]**2) for p in corners]
            heights = [p[2] for p in corners]
            # Every smooth circle's radius below the box's least, or its
            # height outside the box's range.
            for rr, ww in zip(rs, ws):
                self.assertTrue(mpq(rr) < min(radii) or not min(heights) < mpq(ww) < max(heights), spec)

    def test_declared(self):
        # The touching ball: its centre 5/4 + 5/8 from the round's centre
        # `(11/4, 7)` along `(3/5, 4/5)`, inside the round's quarter turn.
        _, r, frame, _, _ = fixtures.STEP_TOUCH
        cx, cz = F(frame[0]), F(frame[2])
        self.assertEqual((cx-F(11, 4), cz-7), (F(15, 8)*F(3, 5), F(15, 8)*F(4, 5)))
        self.assertEqual(F(r), F(5, 8))
        # Clear of the top's disc (radius 11/4 at 33/4) and the shaft (to 7).
        self.assertTrue(cz-F(r) > 7)
        self.assertGreater(cx-F(r)*F(3, 5), F(11, 4))
        self.assertEqual({c.reason for c in fixtures.all_cases() if c.group == 'step_touch'}, {fixtures.TOUCH})

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 15)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 12, 'empty': 0, 'degenerate': 3})
        for c in listed:
            self.assertEqual(c.encode().count('\nbrep imported/fillet_'), 1, c.name)
            self.assertEqual(c.native().count('\nbrep imported/fillet_'), 1, c.name)
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
        import generate_turned_profiles_boolean_fixtures as profiles
        import generate_turned_flats_boolean_fixtures as flats
        paths = {b.path for b in fixtures.BODIES.values()}
        self.assertFalse(paths & {f'imported/{n}.brep' for n in w.BODIES})
        self.assertFalse(paths & {b.path for b in profiles.BODIES.values()})
        self.assertFalse(paths & {b.path for b in flats.BODIES.values()})


if __name__ == '__main__':
    unittest.main()
