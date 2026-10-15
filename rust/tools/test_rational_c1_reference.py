"""U10 (b): the rational C1 condition of spline_cell_reference.py checked
independently with mpmath.

At a knot of multiplicity equal to the degree a rational B-spline
`P = A / w` is C1 when `(A' w - A w') / w^2` agrees from both sides. The
reference decides it exactly in `Fraction`s by de Boor's recurrence (and on
a surface by the jump's polynomial on each span). Here every curve and
surface of the U10 evidence (generate_brep_fixtures.py's rational cells and
models) is evaluated again by the Cox-de Boor basis recursion in mpmath at
60 digits, from the stored binary64 data converted exactly, and differenced
one-sidedly at the knot with a step of 1e-25: the quotients' difference
must be the reference's exact jump within 1e-15 of the curve's size, and
the reference's verdict must be that jump's being zero. No math.cos, sin or
hypot is used anywhere: the quarter arcs' weight is sqrt(2)/2 correctly
rounded (mpmath), and every point is from the stored data.
"""
from fractions import Fraction as F
import unittest

import mpmath as mp

import generate_brep_fixtures as g
import spline_cell_reference as spline

STEP = mp.mpf(10)**-25
# The one-sided quotients' error is about STEP times the second derivative;
# the comparison's bound is far above that and far below any jump here.
BOUND = mp.mpf(10)**-15


def flat_knots(basis):
    """The flat knot sequence (a periodic basis's extended by the
    neighbouring period's knots, OCCT's convention) and the pole count."""
    flat, n, _ = basis.flat()
    return [mp.mpf(k.numerator)/k.denominator for k in flat], n


def basis_value(t, i, p, u, left):
    """Cox-de Boor: N_{i,p}(u), spans half-open to the right, or to the left
    when `left` (a limit from below)."""
    if p == 0:
        return 1 if ((t[i] < u <= t[i+1]) if left else (t[i] <= u < t[i+1])) else 0
    out = mp.mpf(0)
    if t[i+p] != t[i]:
        out += (u-t[i])/(t[i+p]-t[i])*basis_value(t, i, p-1, u, left)
    if t[i+p+1] != t[i+1]:
        out += (t[i+p+1]-u)/(t[i+p+1]-t[i+1])*basis_value(t, i+1, p-1, u, left)
    return out


def homogeneous_poles(poles, weights):
    return [[mp.mpf(w)*mp.mpf(x) for x in p]+[mp.mpf(w)] for p, w in zip(poles, weights)]


def curve_point(basis, h, u, left=False):
    """The rational point at u (a mpf) by the basis functions."""
    t, n = flat_knots(basis)
    count = len(t)-basis.degree-1
    acc = [mp.mpf(0)]*len(h[0])
    for i in range(count):
        b = basis_value(t, i, basis.degree, u, left)
        if b:
            pole = h[i % n if basis.periodic else i]
            acc = [a+b*x for a, x in zip(acc, pole)]
    return [a/acc[-1] for a in acc[:-1]]


def surface_point(s, u, v, left=False):
    """The rational point at (u, v), one-sided in u when `left`."""
    tu, nu = flat_knots(s.u)
    tv, nv = flat_knots(s.v)
    h = homogeneous_poles(s.poles, s.weights)
    acc = [mp.mpf(0)]*4
    for i in range(nu):
        bu = basis_value(tu, i, s.u.degree, u, left)
        if not bu:
            continue
        for j in range(nv):
            bv = basis_value(tv, j, s.v.degree, v, False)
            if bv:
                acc = [a+bu*bv*x for a, x in zip(acc, h[i*nv+j])]
    return [a/acc[-1] for a in acc[:-1]]


def mpf(x):
    return mp.mpf(x.numerator)/x.denominator


def exact_jump(basis, hpoles, ul, ur):
    """The reference's exact rational derivatives' jump at a tested knot."""
    ctrl = spline._controls(basis, hpoles)
    _, left = spline.rational_jet(*spline._one_sided(basis, ctrl, ul, True))
    _, right = spline.rational_jet(*spline._one_sided(basis, ctrl, ur, False))
    return [r-l for l, r in zip(left, right)]


class RationalC1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cells = {c.name: c for c in g.rational_c1_cells(g.base_cases())}

    def check_curve(self, curve):
        """Every tested knot's quotients against the exact jump; returns the
        reference's verdict, which must be every jump's being zero."""
        hpoles = spline.homogeneous(curve.poles, curve.weights)
        lo, hi = spline.span_of(curve)
        size = max(abs(mp.mpf(x)) for p in curve.poles for x in p)
        smooth = True
        with mp.workdps(60):
            h = homogeneous_poles(curve.poles, curve.weights)
            for ul, ur in spline._tested(curve.basis, lo, hi):
                jump = exact_jump(curve.basis, hpoles, ul, ur)
                a, b = mpf(ul), mpf(ur)
                at_l, at_r = curve_point(curve.basis, h, a, left=True), curve_point(curve.basis, h, b)
                for x, y in zip(at_l, at_r):  # C0 (a seam's two ends)
                    self.assertLess(abs(x-y), BOUND*size)
                below = curve_point(curve.basis, h, a-STEP, left=True)
                above = curve_point(curve.basis, h, b+STEP)
                for k in range(len(jump)):
                    d_left = (at_l[k]-below[k])/STEP
                    d_right = (above[k]-at_r[k])/STEP
                    self.assertLess(abs((d_right-d_left)-mpf(jump[k])), BOUND*max(size, 1))
                smooth = smooth and not any(jump)
        self.assertEqual(spline.curve_c1(curve), smooth)
        return smooth

    def test_nurbs_circle(self):
        """Quarter arcs at double knots: C1 rationally, not homogeneously,
        and on the circle within the weight's rounding; with one span twice
        as long, G1 only (a jump along the tangent)."""
        circle = g.nurbs_circle(1.5, 0.0)
        self.assertTrue(self.check_curve(circle))
        self.assertFalse(spline.homogeneous_curve_c1(circle))
        with mp.workdps(60):
            h = homogeneous_poles(circle.poles, circle.weights)
            for k in range(41):
                x, y, _ = curve_point(circle.basis, h, mp.mpf(k)/10)
                self.assertLess(abs(mp.sqrt(x*x+y*y)-mp.mpf(1.5)), mp.mpf(10)**-15)
        g1 = g.nurbs_circle(1.5, 0.0, (0.0, 1.0, 2.0, 3.0, 5.0))
        self.assertFalse(self.check_curve(g1))
        for name in ('rational_ring_nurbs_circle', 'rational_ring_nurbs_circle_g1'):
            ring = next(e.curve for e in self.cells[name].edges if e.start is None)
            self.assertEqual(self.check_curve(ring), not name.endswith('_g1'))

    def test_half_circle(self):
        """The rational cylinder's circle (STEP-b's rational_cylinder, the
        half disc's wall) and the same with spans 1 and 2."""
        for knots, smooth in (((0.0, 1.0, 2.0), True), ((0.0, 1.0, 3.0), False)):
            poles, weights, basis = g.nurbs_half_circle(5.0, knots)
            curve = spline.BSpline3(basis, [(x, y, 0.0) for x, y in poles], weights)
            self.assertEqual(self.check_curve(curve), smooth)
            self.assertFalse(spline.homogeneous_curve_c1(curve))
        for name, smooth in (('rational_half_disc', True), ('rational_half_disc_g1', False)):
            c = self.cells[name]
            curves = [e.curve for e in c.edges if isinstance(e.curve, spline.BSpline3)]
            self.assertEqual(len(curves), 2)
            for curve in curves:
                self.assertEqual(self.check_curve(curve), smooth)

    def test_line_weights(self):
        """A line spline over a double knot: the knot's weight doubled alone
        keeps the rational speed (C1, homogeneously not); with the next
        weight tripled the speed jumps along the line (G1, refused); the
        homogeneous-C1 rational edge of R4's fixtures still passes."""
        a, b = (0.0, 0.0, 0.0), (3.0, 0.0, 0.0)
        poles = g.along(a, b, g.QUADRATIC_LINE)
        c1 = spline.BSpline3(g.QUADRATIC, poles, g.WEIGHTS_C1)
        self.assertTrue(self.check_curve(c1))
        self.assertFalse(spline.homogeneous_curve_c1(c1))
        self.assertFalse(self.check_curve(spline.BSpline3(g.QUADRATIC, poles, g.WEIGHTS_G1)))
        # spline_edge_rational_c1's construction: w2 P2 the mean of w1 P1 and
        # w3 P3 with weights 1, 2, 3 about the knot.
        p1, p3 = [F(3, 4), F(0), F(0)], [F(9, 4), F(0), F(0)]
        p2 = [(x+3*y)/4 for x, y in zip(p1, p3)]
        homogeneous = spline.BSpline3(g.QUADRATIC, [a, g.exact(p1), g.exact(p2), g.exact(p3), b],
                                      [1.0, 1.0, 2.0, 3.0, 1.0])
        self.assertTrue(spline.homogeneous_curve_c1(homogeneous))
        self.assertTrue(self.check_curve(homogeneous))
        for name, smooth in (('rational_edge_weights_c1', True), ('rational_edge_weights_g1', False)):
            self.assertEqual(self.check_curve(self.cells[name].edges[0].curve), smooth)
        for name, smooth in (('rational_pcurve_weights_c1', True), ('rational_pcurve_weights_g1', False)):
            c = self.cells[name]
            pcurve = c.fins[c.loops[c.faces[1].loops[0]].fins[0]].pcurve
            self.assertEqual(self.check_curve(pcurve), smooth)

    def check_surface(self, s, smooth):
        """The one-sided quotients of S_u across the knot line u = 1.5 at
        points along it: all equal when the reference says C1, and some
        apart by more than 1e-6 otherwise."""
        self.assertEqual(spline.surface_c1(s), smooth)
        k = mp.mpf(1.5)
        worst = mp.mpf(0)
        with mp.workdps(60):
            for j in range(9):
                v = mp.mpf(2)*(2*j+1)/18  # inside the domain
                at_l, at_r = surface_point(s, k, v, left=True), surface_point(s, k, v)
                below, above = surface_point(s, k-STEP, v, left=True), surface_point(s, k+STEP, v)
                for c in range(3):
                    self.assertLess(abs(at_l[c]-at_r[c]), BOUND)
                    jump = abs((above[c]-at_r[c])/STEP-(at_l[c]-below[c])/STEP)
                    worst = max(worst, jump)
        if smooth:
            self.assertLess(worst, BOUND)
        else:
            self.assertGreater(worst, mp.mpf(10)**-6)

    def test_rational_caps(self):
        """The cap's surface: C1 rationally only; G1 only; and rows each C1
        rationally (the middle one homogeneously) whose surface is not, as
        their jumps are different multiples of their points."""
        c1 = g.rational_cap([g.WEIGHTS_C1]*3)
        self.check_surface(c1, True)
        self.assertFalse(spline.homogeneous_surface_c1(c1))
        self.check_surface(g.rational_cap([g.WEIGHTS_G1]*3), False)
        rows = g.rational_cap([g.WEIGHTS_C1, [1.0]*5, g.WEIGHTS_C1])
        self.check_surface(rows, False)
        rows_u, _ = spline._rows(rows)
        self.assertTrue(all(spline._row_rational_c1(rows.u, r) for r in rows_u))
        for name, smooth in (('rational_face_c1', True), ('rational_face_g1', False),
                             ('rational_face_rows_not_c1', False)):
            self.check_surface(self.cells[name].faces[1].surface, smooth)

    def test_half_disc_wall(self):
        """The half disc's rational cylinder wall across its knot line, by
        the same quotients along the ruling."""
        for name, smooth in (('rational_half_disc', True), ('rational_half_disc_g1', False)):
            wall = next(f.surface for f in self.cells[name].faces if isinstance(f.surface, spline.BSplineSurface))
            self.assertEqual(spline.surface_c1(wall), smooth)
            self.assertFalse(spline.homogeneous_surface_c1(wall))
            k = mp.mpf(1)
            worst = mp.mpf(0)
            with mp.workdps(60):
                for j in range(5):
                    v = mp.mpf(2*j+1)/10
                    at_l, at_r = surface_point(wall, k, v, left=True), surface_point(wall, k, v)
                    below = surface_point(wall, k-STEP, v, left=True)
                    above = surface_point(wall, k+STEP, v)
                    for c in range(3):
                        worst = max(worst, abs((above[c]-at_r[c])/STEP-(at_l[c]-below[c])/STEP))
            self.assertEqual(worst < BOUND, smooth)
            if not smooth:
                self.assertGreater(worst, mp.mpf(10)**-6)


if __name__ == '__main__':
    unittest.main()
