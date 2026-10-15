"""Independent checks of S9e.4b.4d.2c.1c's fixtures
(`generate_turned_corners_boolean_fixtures.py`) and of the consistency
generator's corner body (`generate_turned_consistency_fixtures.py`): each
body's closed form (Pappus over its model profile) against a direct
quadrature of its radius, its rounds at corners (tangent to the wider
cylinder at the centre's height, the tube crossing the narrower at the
corner's exact height, a sixth of a turn round for the shaft's first two)
and its tangent round (2c.1b's), a wrong corner caught; the shaft against
`OCC496b`'s numbers; each corner's step past the corner and short of the
tube's top and of the narrower cylinder's far end; the construction (the
sharp profile's lattice, each corner's piece taken from its wider
cylinder's part) against the profile at heights across the meridian; the
reference's turned input exactly and in binary64 against the profile; the
step alone on the chained reference against its closed form; every partner
across the corner ring it is placed across, the declared ball touching a
round between its tangent circle and its corner; the case list and its
protocol rows; and the consistency generator's corner rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_turned_corners_boolean_fixtures as fixtures
import generate_turned_consistency_fixtures as consistency

val = fixtures.val


def heights(model):
    return sorted({val(w) for _, w in model.pts})


class TurnedCornersReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        for b in fixtures.BODIES.values():
            m = b.model
            V, A = b.closed()
            hs = heights(m)
            vol = mp.pi*mp.quad(lambda w: m.radius(w)**2, hs)
            self.assertLess(abs(V-vol), mp.mpf(10)**-25*V, b.name)
            # The area: each disc and annulus, each line's frustum, each
            # arc's band `2 pi int (R + T cos) T` over its angles by
            # quadrature.
            area = mp.mpf(0)
            arcs = m.arcs()
            pts = m.pts
            for i, ((ra, wa), (rb, wb)) in enumerate(zip(pts, pts[1:])):
                Ra, Rb, Wa, Wb = val(ra), val(rb), val(wa), val(wb)
                if i in arcs:
                    R, W, T, th0, th1 = arcs[i]
                    area += 2*mp.pi*abs(mp.quad(lambda th: (R+T*mp.cos(th))*T, [th0, th1]))
                elif Wa == Wb:
                    area += mp.pi*abs(Rb**2-Ra**2)
                else:
                    area += mp.pi*(Ra+Rb)*mp.sqrt((Rb-Ra)**2+(Wb-Wa)**2)
            self.assertLess(abs(A-area), mp.mpf(10)**-25*A, b.name)
        self.assertLess(abs(fixtures.BODIES['step'].closed()[0]-consistency.pappus_of(consistency.CORNER_MODEL)),
                        mp.mpf(10)**-30)

    def test_corners(self):
        shaft = fixtures.BODIES['shaft'].model
        kinds = [c['kind'] for c in shaft.corners]
        self.assertEqual(kinds, ['corner', 'corner', 'tangent'])
        for c, ry in zip(shaft.corners[:2], (fixtures.R14, fixtures.R13)):
            (cr, cw), t = c['centre'], c['t']
            self.assertEqual(t, 2)
            self.assertEqual(c['ry'], ry)
            self.assertTrue(c['up'])
            # A sixth of a turn round within 6e-7 rad.
            cos = (ry-cr)/t
            self.assertLess(abs(mp.acos(val(cos))-mp.pi/3), 6e-7)
            corner = shaft.pts[c['at']+1]
            self.assertIsInstance(corner[1], fixtures.Surd)
            self.assertEqual(corner, (ry, fixtures.Surd(cw, 1, t*t-(ry-cr)**2)))
        step = fixtures.BODIES['step'].model
        (c,) = step.corners
        self.assertEqual((c['centre'], c['t'], c['rx'], c['ry'], c['up']), ((3, 3), 1, 4, F(7, 2), True))
        self.assertEqual(step.pts[3], (F(7, 2), fixtures.Surd(3, 1, F(3, 4))))
        self.assertLess(abs(val(step.pts[3][1])-(3+mp.sqrt(3)/2)), mp.mpf(10)**-35)
        collar = fixtures.BODIES['collar'].model
        (c,) = collar.corners
        self.assertEqual((c['centre'], c['t'], c['rx'], c['ry'], c['up']), ((F(11, 4), 5), F(5, 4), 4, F(7, 2), False))
        # The collar's corner a Pythagorean direction of dyadic half-angle
        # tangent: (3/5, -4/5) from its centre, `s = -1/2`.
        self.assertEqual(collar.pts[2], (F(7, 2), 4))
        self.assertEqual(fixtures.direction(F(-1, 2)), (F(3, 5), F(-4, 5)))
        # A round whose tube misses the narrower cylinder, or a fillet at a
        # corner, is caught.
        bad = [(0, 0, fixtures.L), (4, 0, fixtures.L), (4, 3, fixtures.arc(3, 3, 1, True)),
               (F(2), 3, fixtures.L), (2, 6, fixtures.L), (0, 6, fixtures.L)]
        with self.assertRaises(AssertionError):
            fixtures.Model(bad)
        corner = fixtures.corner_height(F(3), F(5), F(1), F(9, 2), True)
        fillet = [(0, 0, fixtures.L), (4, 0, fixtures.L), (4, 3, fixtures.arc(5, 3, 1, False)),
                  (F(9, 2), corner, fixtures.L), (F(9, 2), 6, fixtures.L), (0, 6, fixtures.L)]
        with self.assertRaises(AssertionError):
            fixtures.Model(fillet)

    def test_shaft_is_occ496b(self):
        # `OCC496b`'s stored numbers, each within 5e-10 of the shaft's.
        stored = {
            'OY': (fixtures.OY, '-4.7683715820310001e-07'), 'OZ': (fixtures.OZ, '38.500011444089999'),
            'W0': (fixtures.W0, '-38.250450134280001'), 'R15': (fixtures.R15, '15.000000476837158'),
            'C1': (fixtures.C1, '-33.750469207759998'), 'R14': (fixtures.R14, '13.999999523167158'),
            'C2': (fixtures.C2, '-31.35044670105'), 'R13': (fixtures.R13, '13.000000476837158'),
            'C3': (fixtures.C3, '-20.23610687256'), 'WP': (fixtures.WP, '-3.0000028610229998'),
            'RT': (fixtures.RT, '5.2986021041871583'),
        }
        for name, (x, text) in stored.items():
            self.assertLess(abs(val(x)-mp.mpf(text)), 5e-10, name)
            self.assertEqual((F(x)*2**30).denominator, 1, name)
        # Its second round's stored edge 9 2^-19 below its centre
        # (`OCC496b`'s 1.7166e-5), OCCT's tube 7.4e-11 off that point.
        b = fixtures.BODIES['shaft']
        self.assertEqual(b.occt[4][:2], (fixtures.R14, fixtures.C2-F(9, 2**19)))
        self.assertLess(abs(F(9, 2**19)-F('1.71661376953e-05')), F(1, 10**15))
        off = mp.sqrt(4+val(F(9, 2**19))**2)-2
        self.assertTrue(7e-11 < off < 8e-11)
        # The third round tangent to the cylinder of radius 13 + 2^-21 (its
        # centre 9.5e-7 off `OCC496b`'s 8.499999523163158).
        c3 = b.model.corners[2]
        self.assertEqual(c3['centre'], (fixtures.R13-F(9, 2), fixtures.C3))
        self.assertLess(abs(val(c3['centre'][0])-mp.mpf('8.499999523163158')-mp.mpf('9.5367e-7')), 1e-11)
        # The resolution and the axis: the world's x.
        rows = b.native_rows()
        self.assertEqual(rows[-1], 'tolerance 0.003')
        self.assertTrue(rows[2].startswith('revolve 0 ') and rows[2].endswith(' 1 0 0 6.283185307179586'))

    def test_steps(self):
        for b in fixtures.BODIES.values():
            m = b.model
            for c in m.corners:
                if c['kind'] != 'corner':
                    continue
                (cr, cw), t, h = c['centre'], c['t'], c['step']
                i = c['at']
                corner = m.pts[i+1][1] if c['up'] else m.pts[i][1]
                far = m.pts[i+2][1] if c['up'] else m.pts[i-1][1]
                top = cw+t if c['up'] else cw-t
                lo, hi = sorted((val(corner), val(min(top, far) if c['up'] else max(top, far))))
                self.assertTrue(lo < val(h) < hi, (b.name, h))
                g = fixtures.power_above(t)/2**24
                self.assertEqual((h/g).denominator, 1)
        self.assertEqual(fixtures.BODIES['collar'].model.corners[0]['step'], F(31, 8))

    def test_construction(self):
        # The kernel's construction (the sharp profile's lattice, each
        # corner's piece taken from its wider cylinder's part) is the
        # profile at heights across each body, off its nodes.
        for b in fixtures.BODIES.values():
            m = b.model
            hs = heights(m)
            lo, hi = hs[0], hs[-1]
            n = 0
            # Off a tangent round's piece (2c.1a's and 2c.1b's, taken from
            # the lattice: its heights from its tangent circle on the
            # cylinder to its lines' corner).
            skip = []
            for c in m.corners:
                if c['kind'] == 'tangent':
                    i = c['at']
                    ws = [val(m.pts[i][1]), val(m.pts[i+1][1]), val(c['corner'][1])]
                    skip.append((min(ws), max(ws)))
            for k in range(400):
                w = lo+(hi-lo)*(mp.mpf(k)+mp.mpf(3)/7)/400
                if any(a <= w <= b for a, b in skip):
                    continue
                want = m.radius(w)
                got = m.construction_radius(w)
                if want is None:
                    continue
                self.assertIsNotNone(got, (b.name, w))
                self.assertLess(abs(got-want), mp.mpf(10)**-30, (b.name, w))
                n += 1
            self.assertGreater(n, 300)
            # Without the pieces the sharp profile is off the round.
            c = next(c for c in m.corners if c['kind'] == 'corner')
            (cr, cw), t = c['centre'], c['t']
            w = val(cw)+(mp.mpf(1)/4 if c['up'] else -mp.mpf(1)/4)*t
            self.assertGreater(val(m.lattice.form(F(float(w)))), m.radius(w))

    def test_membership(self):
        # The reference's turned input (exactly and in binary64) against the
        # profile's radius at a grid of meridian points.
        for b in fixtures.BODIES.values():
            m = b.model
            turned = ref.make_input(fixtures.TurnedCase(b))
            self.assertIsInstance(turned, fixtures.TurnedInput)
            fl = ref.Float(turned)
            hs = heights(m)
            r_max = max(val(r) for r, _ in m.pts)
            count = 0
            for i in range(1, 41):
                for j in range(-1, 42):
                    r = r_max*(mp.mpf(i)+mp.mpf(1)/13)/38
                    w = hs[0]+(hs[-1]-hs[0])*(mp.mpf(j)+mp.mpf(5)/13)/40
                    want = m.inside(r, w)
                    X = b.frame_point(r, w)
                    self.assertEqual(turned.contains(X), want, (b.name, r, w))
                    self.assertEqual(fl.contains([float(x) for x in X]), want, (b.name, r, w))
                    count += want
            self.assertGreater(count, 300, b.name)

    def test_step_on_the_chained_reference(self):
        scan, ref.SCAN = ref.SCAN, 24
        try:
            chain = ref.Chain([fixtures.TurnedCase(fixtures.BODIES['step'])])
            want, mom, area = chain.inputs[0].closed()
            for k in (0, 1):
                V, mm, A = chain.measures(0, k)
                self.assertLess(abs(V-want), mp.mpf(10)**-26*want)
                self.assertLess(abs(A-area), mp.mpf(10)**-26*area)
                for i in range(3):
                    self.assertLess(abs(mm[i]-mom[i]), mp.mpf(10)**-26*want*10)
        finally:
            ref.SCAN = scan

    def test_partners_cross(self):
        bodies = fixtures.BODIES

        def ring(name, k=0):
            m = bodies[name].model
            c = [c for c in m.corners if c['kind'] == 'corner'][k]
            i = c['at']
            return m.pts[i+1] if c['up'] else m.pts[i]

        # Boxes on `XY` about the world's z: sides y = +-a across the ring
        # (its points there strictly inside the side), caps off it.
        for spec, (r, w), a in ((fixtures.STEP_BOX, ring('step'), F(1, 2)),
                                (fixtures.COLLAR_BOX, ring('collar'), F(1, 2))):
            _, (sq,), frame, h0, h1 = spec
            xs = sorted({F(p[0]) for p in sq.points})
            ys = sorted({F(p[1]) for p in sq.points})
            self.assertEqual(ys, [-a, a])
            x = mp.sqrt(val(r)**2-val(a)**2)
            self.assertTrue(val(xs[0]) < x < val(xs[1]) and val(h0) < val(w) < val(h1), spec)
        # The shaft's box: its sides y = +-1 across the first ring about the
        # world's x through (0, OY, OZ); its bottom below the ring's top.
        r, w = ring('shaft')
        _, (sq,), frame, h0, h1 = fixtures.SHAFT_BOX
        xs = sorted({F(p[0]) for p in sq.points})
        ys = sorted({F(p[1]) for p in sq.points})
        self.assertTrue(val(xs[0]) < val(w) < val(xs[1]))
        for y in ys:
            z = val(fixtures.OZ)+mp.sqrt(val(r)**2-(val(y)-val(fixtures.OY))**2)
            self.assertTrue(val(h0) < z < val(h1))
        # The rods: the ring's nearest point inside, its far side outside.
        _, (disc,), frame, h0, h1 = fixtures.STEP_ROD
        yc, zc, rho = (val(F(c)) for c in disc.circle)
        r, w = ring('step')
        self.assertLess(abs(val(w)-zc), rho)
        self.assertTrue(val(F(frame[0])+F(h0)) < -val(r) and val(F(frame[0])+F(h1)) > val(r))
        _, (disc,), frame, h0, h1 = fixtures.SHAFT_ROD
        zc, xc, rho = (val(F(c)) for c in disc.circle)
        r, w = ring('shaft')
        top = val(fixtures.OZ)+val(r)
        self.assertLess(mp.sqrt((val(w)-xc)**2+(top-zc)**2), rho)
        self.assertGreater(abs(val(w)-xc), 0)
        # The skew cubes and the ball about a ring's point.
        for spec, (r, w), centre in ((fixtures.STEP_SKEW, ring('step'), (F(7, 2), 0, 4)),
                                     (fixtures.SHAFT_SKEW, ring('shaft', 1), (F(-237, 8), 0, F(103, 2)))):
            _, (sq,), frame, h0, h1 = spec
            self.assertEqual(frame[3:], fixtures.forms.FRAMES['SKEW'])
            self.assertEqual(tuple(F(c) for c in frame[:3]), centre)
        r, w = ring('shaft', 1)
        p = (val(w), val(fixtures.OY), val(fixtures.OZ)+val(r))
        self.assertLess(max(abs(p[0]+mp.mpf(237)/8), abs(p[2]-mp.mpf(103)/2)), 0.2)
        _, rho, frame = fixtures.COLLAR_BALL[:3]
        r, w = ring('collar')
        centre = tuple(F(c) for c in frame[:3])
        self.assertEqual(centre, (F(29, 8), 0, 4))
        self.assertLess(abs(centre[0]-r), F(rho))
        self.assertLess(centre[0]+F(rho), 4)
        self.assertLess(abs(val(w)-4), mp.mpf(10)**-30)

    def test_declared(self):
        # The ball touching `step`'s round at (19/5, 0, 18/5): its centre
        # 1 + 1/4 from the round's along (4/5, 3/5), between the tangent
        # circle (angle 0) and the corner (a sixth of a turn).
        _, rho, frame = fixtures.STEP_TOUCH[:3]
        centre = tuple(F(c) for c in frame[:3])
        self.assertEqual((F(rho), centre), (F(1, 4), (4, 0, F(15, 4))))
        c = fixtures.BODIES['step'].model.corners[0]
        (cr, cw), t = c['centre'], c['t']
        d = (centre[0]-cr, centre[2]-cw)
        self.assertEqual(d[0]**2+d[1]**2, (t+F(rho))**2)
        self.assertTrue(0 < mp.atan2(val(d[1]), val(d[0])) < mp.pi/3)
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
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), 3)
        self.assertEqual(text.count('\nrevolve '), 3)
        self.assertEqual(text.count('\ntolerance 0.003\n'), 1)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))

    def test_consistency_corners(self):
        rows = consistency.corners_text().splitlines()
        self.assertTrue(rows[1].startswith('body\tstep imported/corner_step.brep volume 273.39888274647'))
        self.assertTrue(rows[2].startswith('case\tstep_box box 3 -0.5 3.25 1.5 1 1.25 volume 1.875 common '))
        self.assertTrue(rows[3].startswith('case\tbox_off box 2.25 2.25 3.25 1 1 1.25 volume 1.25 common '))
        self.assertTrue(rows[4].startswith('case\tpost post 3.5 0 0.25 2.25 5 volume '))
        self.assertEqual(len(rows), 5)
        _, (sq,), frame, h0, h1 = fixtures.STEP_BOX
        xs = sorted({F(p[0]) for p in sq.points})
        ys = sorted({F(p[1]) for p in sq.points})
        self.assertEqual((xs[0], ys[0], F(h0), xs[1]-xs[0], ys[1]-ys[0], F(h1)-F(h0)),
                         consistency.CORNER_BOXES[0][1])


if __name__ == '__main__':
    unittest.main()
