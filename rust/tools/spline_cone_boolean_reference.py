#!/usr/bin/env python3
"""Independent reference for S9f.3b of REVIEW_NOTES.md: Booleans of a spline
prism (a profile of lines and S8b's nonrational spline segments) against a
cone or frustum (`Solid::cone_with`: radius `bottom` at the frame's origin
and `top` at `height` along its normal, a zero radius an apex) in any
relative position, either the object: every spline wall meets the cone in
a curve that is, along the wall's ruling at the spline's parameter `t`, a
root of `A w^2 + 2 B(t) w + C(t)`, `A = q_u^2 + q_v^2 - k^2 q_w^2` (`q` the
prism's axis in the cone's frame, `k` its slope) of either sign.

An extension of `spline_sphere_boolean_reference.py`: its prism
(`Spline3`, the profile's elements as exact power polynomials, its parsing
and ray test shared with `curved_boolean_reference.py`), Green's integrals
of segments and ellipses, its loops and its solids' combinatorics are
reused; the cone, its sections, its faces' classes and the breakpoints are
this module's. Each input is its construction's exact model in rationals
from the stored binary64 data (`stored_axes`): the prism as S9f.3a's; the
cone the points `o + u x + v y + w n` of its stored frame with `0 <= w <=
h` and `u^2 + v^2 <= (b + k w)^2` (`k = (top - bottom) / h`), its frame
exact here (axes along the world's, either sense: `x * y = n`, unit), as
S9f.3a's hemispheres; any relative position comes from the prism's frame.
Nothing here uses the kernel, a surface/surface intersection or an
arrangement.

* **Slicing (volume, first moments, solids, the cone's wall).** Both
  solids are sliced by the planes `d . X = s` with `d . n_p > 0` and `d`
  within the cone's elliptic directions (`(d . n)^2 > k^2 |d x n|^2`, so
  every slice cuts the cone's quadric in an ellipse): the first `n` leaned
  toward the prism's axis `n_p` as little as that needs (`n` itself where
  `n . n_p > 0`), the second a further lean along `(2, -3, 5)`. Every point
  of a slice is projected along the prism's axis onto its profile plane, so
  the prism's section is the profile cut by the strip of its heights
  (S9f.3a's), and the cone's is the quadric restricted to the slice (in an
  orthonormal basis of the slice its quadratic form, centre and principal
  axes: an ellipse `E0 + cos(t) A1 + sin(t) A2` after the projection),
  cut by the lines of its end planes `w = 0` and `w = h`. The regions are
  taken apart by their boundaries as S9f.3a's (every curve cut at its
  crossings with the other region's curves and classified at its pieces'
  midpoints), their areas and first moments in `(u, v)` by Green's theorem
  in closed form. The cone's wall: along the slice, the ruling at the
  cone's angle `theta` has `r = G(s) / D(theta)`, `D = d . n + k d . (cos
  x + sin y)` (never zero: `d` is elliptic), `G = b d . n + k (s - d .
  o)`, and `dw = ds / D`, so its area element `r sqrt(1 + k^2) dtheta dw`
  is `sqrt(1 + k^2) |G(s)| dtheta ds / D^2` and its part inside the prism
  `sqrt(1 + k^2) |G(s)|` times the integral of `D^-2` over the section's
  arcs inside the prism's section, in closed form (`D = mu + R cos(theta -
  phi)`, `|mu| > R`: `atan` and a rational term, unwrapped across turns).
  Breakpoints: every vertex's slice and every edge's extremes along `d`
  (the prism's vertices and cap edges' extremes; the apex, real or virtual;
  the rims' extremes; the cap planes' conics' extremes, on the line where
  the quadric's gradient lies in the span of the cap's normal and `d`;
  vertical edges, cap edges and rims against the quadric, the end planes
  and the caps; the meeting's extremes, `F = 0` and `d . (F_w S' - F_t n_p)
  = 0` linear in `w`, eliminated: `A N^2 - 2 B N D + C D^2`; the end
  planes' creases' extremes). Each interval by
  `curved_boolean_reference.integrate`, as S9f.3a's.
* **Surface area.** The prism's walls swept along their generatrices:
  along each, the part inside the cone (`F <= 0`: between the roots for `A
  > 0`, outside them for `A < 0`, a half line for `A = 0`) within the end
  planes' slab and the prism's heights, times `|S'(t) x n_p|`, integrated
  over `t` between the roots of `B^2 - A C`, of `F` at the caps' heights
  and along the end planes' creases, and of the creases at the caps. The
  caps: by chords in their plane along an irrational direction, each
  chord's part inside the profile (its crossings) and inside the cone's
  conic there (a quadratic of any type) within the slab, integrated across
  between the chords through the profile's vertices, its spans' tangents
  along the chords, the conic's tangents along them, the conic's and the
  slab's lines' crossings of the profile and each other. The cone's wall
  by the slicing above; its end discs as S9f.3a's hemisphere's disc (the
  projection along `n_p` where the plane is transverse to it, else the
  profile's chords along the plane's trace times the heights).
* **Solids** (combinatorics only, at 20 digits): the pair's box cut into
  slices across the prism's axis and each slice into chords, on each chord
  both inputs' intervals exactly (the profile's chords within the prism's
  heights, the cone's quadric within its slab), the operation's from them,
  joined where they overlap on neighbouring chords and slices (`GRID` by
  `GRID`; S9f.3a's sections' overlaps in the profile's plane miss a thin
  layer under a cap that a slice leaning from the prism's axis crosses as a
  moving sliver).

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9f
references.
"""
from fractions import Fraction as F

import mpmath as mp

from curve_surface_reference import stored_axes
import curved_boolean_reference as cref
import spline_sphere_boolean_reference as sref
from boolean_reference import pder, pmul, psub, padd, peval, pdeg

mp.mp.dps = 40

OPS = sref.OPS
SECOND = sref.SECOND
M, Mv = cref.M, cref.Mv
sub, add, scale, dot, cross, apply = cref.sub, cref.add, cref.scale, cref.dot, cref.cross, cref.apply
pc, pm, mnorm, exact_roots = sref.pc, sref.pm, sref.mnorm, sref.exact_roots
vpoly_dot, wall_length = sref.vpoly_dot, sref.wall_length
Piece, seg_green, ellipse_green = sref.Piece, sref.seg_green, sref.ellipse_green
# The chords' direction in a cap's plane: an irrational slope, along no
# rational line and no conic's asymptote of these fixtures.
CHORD_ANGLE = mp.mpf('0.6180339887498948482045868343656381177203')


# ------------------------------------------------------------------ the cone

class Cone:
    """A cone or frustum on its exact model, its frame exact."""

    def __init__(self, case):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        assert cross(x, y) == n and all(dot(v, v) == 1 for v in (x, y, n)) and dot(x, y) == 0, \
            'a cone in an exact frame'
        bottom, top, height = case.cone
        self.case = case
        self.o, self.x, self.y, self.n = o, x, y, n
        self.b, self.t, self.h = F(bottom), F(top), F(height)
        assert self.b >= 0 and self.t >= 0 and self.b+self.t > 0 and self.h > 0 and self.b != self.t
        self.k = (self.t-self.b)/self.h
        self.om, self.xm, self.ym, self.nm = Mv(o), Mv(x), Mv(y), Mv(n)
        self.bm, self.tm, self.hm, self.km = M(self.b), M(self.t), M(self.h), M(self.k)
        self.N = mp.sqrt(1+self.km**2)
        # The apex: a vertex where a radius is zero, else virtual.
        self.apex_w = -self.b/self.k
        self.apex = add(o, scale(n, self.apex_w))
        self.real_apex = self.b == 0 or self.t == 0
        self.ends = [(F(0), self.b), (self.h, self.t)]

    def local(self, X):
        Y = sub(X, self.om)
        return dot(self.xm, Y), dot(self.ym, Y), dot(self.nm, Y)

    def local_exact(self, X):
        Y = sub(X, self.o)
        return dot(self.x, Y), dot(self.y, Y), dot(self.n, Y)

    def q(self, X):
        u, v, w = self.local(X)
        return u*u+v*v-(self.bm+self.km*w)**2

    def gradient(self, X):
        """Half the quadric's gradient (world)."""
        u, v, w = self.local(X)
        rho = self.bm+self.km*w
        return tuple(u*self.xm[i]+v*self.ym[i]-self.km*rho*self.nm[i] for i in range(3))

    def contains(self, X):
        u, v, w = self.local(X)
        return 0 < w < self.hm and u*u+v*v < (self.bm+self.km*w)**2

    def closed(self):
        """Volume, first moments and area in closed form."""
        b, t, h = self.bm, self.tm, self.hm
        s2 = b*b+b*t+t*t
        V = mp.pi*h*s2/3
        wbar = h*(b*b+2*b*t+3*t*t)/(4*s2)
        centre = add(self.om, scale(self.nm, wbar))
        area = mp.pi*(b+t)*mp.sqrt(h*h+(t-b)**2)+mp.pi*(b*b+t*t)
        return V, scale(centre, V), area

    def lateral(self):
        b, t, h = self.bm, self.tm, self.hm
        return mp.pi*(b+t)*mp.sqrt(h*h+(t-b)**2)


def cone_quadratic(S, el, cone):
    """The cone's function along the element's wall's rulings `P(t) + w n_p`:
    `A w^2 + 2 B(t) w + C(t)`, (A, B, C) exact, and the rows `P_u`, `P_v`,
    `P_w` (polynomials) and `q` (the axis's)."""
    P, _ = S.wall(el)
    rows = []
    for e in (cone.x, cone.y, cone.n):
        rows.append(psub(vpoly_dot(e, P), [dot(e, cone.o)]))
    q = tuple(dot(e, S.n) for e in (cone.x, cone.y, cone.n))
    k, b = cone.k, cone.b
    rho = padd([b], pc(rows[2], k))
    A = q[0]*q[0]+q[1]*q[1]-k*k*q[2]*q[2]
    B = padd(padd(pc(rows[0], q[0]), pc(rows[1], q[1])), pc(rho, -k*q[2]))
    C = psub(padd(pmul(rows[0], rows[0]), pmul(rows[1], rows[1])), pmul(rho, rho))
    return A, B, C, rows, q


def axis_term(S, cone):
    """`A` of the pair: the prism's axis in the cone's rows."""
    q = tuple(dot(e, S.n) for e in (cone.x, cone.y, cone.n))
    return q[0]*q[0]+q[1]*q[1]-cone.k*cone.k*q[2]*q[2]


def elliptic(d, cone):
    """`(d . n)^2 - k^2 |d x n|^2`, positive for a direction whose planes cut
    the cone's quadric in ellipses."""
    dn = dot(d, cone.n)
    perp = dot(d, d)-dn*dn
    return dn*dn-cone.k*cone.k*perp


def first_direction(S, cone):
    """`sigma n + lambda n_p`, the first with `d . n_p > 0` and an elliptic
    margin (its lean at most half the cone's complement's tangent)."""
    nn = dot(cone.n, S.n)
    sigma = 1 if nn >= 0 else -1
    for lam in (F(0), F(1, 8), F(1, 4), F(1, 2), F(1), F(2)):
        d = add(scale(cone.n, sigma), scale(S.n, lam))
        if dot(d, S.n) <= 0:
            continue
        dn = dot(d, cone.n)
        if dn*dn > 4*cone.k*cone.k*(dot(d, d)-dn*dn):
            return d
    raise AssertionError('no elliptic slicing direction')


def second_direction(S, cone, d1):
    for eps in (F(1, 8), F(1, 16), F(1, 32), F(1, 64)):
        for sg in (1, -1):
            d = add(d1, scale(SECOND, sg*eps))
            if dot(d, S.n) <= 0:
                continue
            dn = dot(d, cone.n)
            if dn*dn > 2*cone.k*cone.k*(dot(d, d)-dn*dn):
                return d
    raise AssertionError('no second slicing direction')


# ------------------------------------------------------------------ a slicing direction

class Way(sref.Way):
    """Slices `d . X = s` of the pair (`d . n_p > 0`, elliptic for the cone)."""

    def __init__(self, pair, d):
        super().__init__(pair, d)
        cone = pair.cone
        assert elliptic(self.d, cone) > 0, 'a slicing direction cutting the cone in ellipses'
        # D(theta) = mu + alpha cos + beta sin; G(s) = gb + gk s.
        self.mu = M(dot(self.d, cone.n))
        self.alpha = M(cone.k*dot(self.d, cone.x))
        self.beta = M(cone.k*dot(self.d, cone.y))
        self.gb = M(cone.b*dot(self.d, cone.n)-cone.k*dot(self.d, cone.o))
        self.gk = cone.km
        R = mp.sqrt(self.alpha**2+self.beta**2)
        self.R, self.phi = R, mp.atan2(self.beta, self.alpha)
        self.ia, self.ib = abs(self.mu), (R if self.mu > 0 else -R)

    def events(self):
        return events(self.pair, self)

    def inv_d2(self, theta):
        """An antiderivative of `D(theta)^-2`, continuous in `theta`: with
        `D^2 = (a + b cos psi)^2`, `psi = theta - phi`, `|b| < a`,
        `a / (a^2 - b^2) I1 - b sin(psi) / ((a^2 - b^2) (a + b cos psi))`,
        `I1 = 2 atan(sqrt((a - b) / (a + b)) tan(psi / 2)) / sqrt(a^2 - b^2)`
        unwrapped across turns."""
        a, b = self.ia, self.ib
        psi = theta-self.phi
        root = mp.sqrt(a*a-b*b)
        kappa = mp.sqrt((a-b)/(a+b))
        turns = mp.floor((psi+mp.pi)/(2*mp.pi))
        rest = psi-2*mp.pi*turns
        at = -mp.pi/2 if rest <= -mp.pi else mp.atan(kappa*mp.tan(rest/2))
        i1 = 2*(at+mp.pi*turns)/root
        return a/(a*a-b*b)*i1-b*mp.sin(psi)/((a*a-b*b)*(a+b*mp.cos(psi)))


# ------------------------------------------------------------------ one slice's regions

class Slice:
    """The pair's sections in the slice `d . X = s` as regions of `(u, v)`:
    the prism's `P` (the profile cut by the strip, unless `strip` is off)
    and the cone's `B` (its quadric's ellipse in the slice, cut by its end
    planes' lines unless `cuts` is off), their boundaries' pieces
    classified."""

    def __init__(self, way, s, strip=True, cuts=True, samples=0):
        self.way, self.s = way, s
        pair = way.pair
        S, cone = pair.S, pair.cone
        self.S, self.cone = S, cone
        w0 = way.w0(s)
        lo, hi = M(S.lo), M(S.hi)
        self.wcoef = (w0, M(way.wu), M(way.wv))
        self.strip_lines = []
        if not strip:
            self.strip_all = True
        elif way.constant:
            self.strip_all = lo < w0 < hi
        else:
            self.strip_all = None
            wu, wv = M(way.wu), M(way.wv)
            g2 = wu*wu+wv*wv
            for k, sgn in ((lo, 1), (hi, -1)):
                base = ((k-w0)*wu/g2, (k-w0)*wv/g2)
                d = (sgn*wv, -sgn*wu)
                self.strip_lines.append((base, d))
        # The ellipse: the quadric on the slice's plane `X0 + alpha e1 + beta
        # e2`.
        self.ellipse = None
        X0 = way.X0(s)
        e1, e2 = way.e
        Y0 = sub(X0, cone.om)
        kk = 1+cone.km**2
        nu = (dot(cone.nm, e1), dot(cone.nm, e2))
        nu0 = dot(cone.nm, Y0)
        bk = cone.bm*cone.km
        Mq = [[1-kk*nu[0]*nu[0], -kk*nu[0]*nu[1]], [-kk*nu[0]*nu[1], 1-kk*nu[1]*nu[1]]]
        g = [dot(Y0, e)-kk*nu0*nu[i]-bk*nu[i] for i, e in enumerate((e1, e2))]
        q0 = dot(Y0, Y0)-kk*nu0*nu0-2*bk*nu0-cone.bm**2
        det = Mq[0][0]*Mq[1][1]-Mq[0][1]*Mq[1][0]
        z = (-(Mq[1][1]*g[0]-Mq[0][1]*g[1])/det, -(-Mq[1][0]*g[0]+Mq[0][0]*g[1])/det)
        qc = q0+g[0]*z[0]+g[1]*z[1]
        centre3 = add(X0, add(scale(e1, z[0]), scale(e2, z[1])))
        rho_c = cone.bm+cone.km*cone.local(centre3)[2]
        if qc < 0 and rho_c > 0:
            # Principal axes of the 2x2 form.
            tr = Mq[0][0]+Mq[1][1]
            disc = mp.sqrt(max((Mq[0][0]-Mq[1][1])**2/4+Mq[0][1]**2, mp.mpf(0)))
            lams = (tr/2+disc, tr/2-disc)
            if Mq[0][1] == 0:
                vecs = ((mp.mpf(1), mp.mpf(0)), (mp.mpf(0), mp.mpf(1))) if Mq[0][0] >= Mq[1][1] else \
                    ((mp.mpf(0), mp.mpf(1)), (mp.mpf(1), mp.mpf(0)))
            else:
                v = (Mq[0][1], lams[0]-Mq[0][0])
                ln = mp.sqrt(v[0]**2+v[1]**2)
                v = (v[0]/ln, v[1]/ln)
                vecs = (v, (-v[1], v[0]))
            a1 = tuple(mp.sqrt(-qc/lams[0])*c for c in vecs[0])
            a2 = tuple(mp.sqrt(-qc/lams[1])*c for c in vecs[1])
            L1, L2 = way.L
            A1 = (a1[0]*L1[0]+a1[1]*L2[0], a1[0]*L1[1]+a1[1]*L2[1])
            A2 = (a2[0]*L1[0]+a2[1]*L2[0], a2[0]*L1[1]+a2[1]*L2[1])
            if A1[0]*A2[1]-A1[1]*A2[0] < 0:
                a2, A2 = (-a2[0], -a2[1]), (-A2[0], -A2[1])
            E0 = S.local(centre3)[:2]
            dt = A1[0]*A2[1]-A1[1]*A2[0]
            K = ((A2[1]/dt, -A2[0]/dt), (-A1[1]/dt, A1[0]/dt))
            ru = mp.sqrt(A1[0]**2+A2[0]**2)
            rv = mp.sqrt(A1[1]**2+A2[1]**2)
            self.ellipse = (E0, A1, A2, K, (E0[0]-ru, E0[0]+ru, E0[1]-rv, E0[1]+rv))
            # Its 3D points at `t`, for the cone's angle.
            self.space = (centre3, tuple(a1[0]*e1[i]+a1[1]*e2[i] for i in range(3)),
                          tuple(a2[0]*e1[i]+a2[1]*e2[i] for i in range(3)))
            step = (self.theta(mp.mpf(10)**-6)-self.theta(mp.mpf(0))) % (2*mp.pi)
            self.sense = 1 if step < mp.pi else -1
        # The end planes' halfplanes `h0 + hu u + hv v > 0` (`w_c` affine in
        # the slice's `(u, v)`).
        self.cuts = []
        self.cuts_all = True
        if cuts:
            nX0 = dot(cone.nm, sub(X0, cone.om))
            nu_, nv_ = dot(cone.nm, way.Xu), dot(cone.nm, way.Xv)
            parallel = dot(cone.n, add(S.x, scale(S.n, way.wu))) == 0 and \
                dot(cone.n, add(S.y, scale(S.n, way.wv))) == 0
            for h0, hu, hv in ((nX0, nu_, nv_), (cone.hm-nX0, -nu_, -nv_)):
                if parallel:
                    if not h0 > 0:
                        self.cuts_all = False
                    continue
                g2 = hu*hu+hv*hv
                self.cuts.append(((h0, hu, hv), ((-h0*hu/g2, -h0*hv/g2), (hv, -hu))))
        self.samples = samples
        self.pieces = []
        self.build()

    # ---- membership

    def in_strip(self, p):
        if self.strip_all is not None:
            return self.strip_all
        w0, wu, wv = self.wcoef
        w = w0+wu*p[0]+wv*p[1]
        return M(self.S.lo) < w < M(self.S.hi)

    def in_ellipse(self, p):
        if self.ellipse is None:
            return False
        E0, _, _, K, _ = self.ellipse
        a = (p[0]-E0[0], p[1]-E0[1])
        x = K[0][0]*a[0]+K[0][1]*a[1]
        y = K[1][0]*a[0]+K[1][1]*a[1]
        return x*x+y*y < 1

    def in_cuts(self, p, skip=None):
        if not self.cuts_all:
            return False
        for i, ((h0, hu, hv), _) in enumerate(self.cuts):
            if i != skip and not h0+hu*p[0]+hv*p[1] > 0:
                return False
        return True

    def in_B(self, p):
        return self.in_ellipse(p) and self.in_cuts(p)

    # ---- crossings (S9f.3a's)

    ellipse_angle = sref.Slice.ellipse_angle
    element_ellipse = sref.Slice.element_ellipse
    element_line = sref.Slice.element_line
    line_ellipse = sref.Slice.line_ellipse
    line_line = staticmethod(sref.Slice.line_line)

    def point_at(self, t):
        E0, A1, A2, _, _ = self.ellipse
        return (E0[0]+mp.cos(t)*A1[0]+mp.sin(t)*A2[0], E0[1]+mp.cos(t)*A1[1]+mp.sin(t)*A2[1])

    def theta(self, t):
        """The cone's angle of the ellipse's point at `t`."""
        c, a1, a2 = self.space
        X = tuple(c[i]+mp.cos(t)*a1[i]+mp.sin(t)*a2[i] for i in range(3))
        u, v, _ = self.cone.local(X)
        return mp.atan2(v, u)

    def build(self):
        els = self.S.elements
        P_present = self.strip_all is not False
        ell = self.ellipse
        lines = list(self.strip_lines) if P_present else []
        cutl = [c[1] for c in self.cuts] if ell is not None and self.cuts_all else []
        el_params = [[] for _ in els]
        el_toggles = []
        ell_params = []
        line_params = [[] for _ in lines]
        line_toggles = [[] for _ in lines]
        cut_params = [[] for _ in cutl]
        cut_toggles = [[] for _ in cutl]
        if P_present:
            for i, el in enumerate(els):
                for k, line in enumerate(lines):
                    for t, lam in self.element_line(el, line):
                        el_params[i].append(t)
                        line_params[k].append(lam)
                        line_toggles[k].append(lam)
                if ell is not None:
                    for t in self.element_ellipse(el):
                        el_params[i].append(t)
                        th = self.ellipse_angle(el.point(t))
                        ell_params.append(th)
                        el_toggles.append(th)
                for j, cl in enumerate(cutl):
                    for t, lam in self.element_line(el, cl):
                        el_params[i].append(t)
                        cut_params[j].append(lam)
                        cut_toggles[j].append(lam)
        if ell is not None:
            for k, line in enumerate(lines):
                for lam in self.line_ellipse(line):
                    line_params[k].append(lam)
                    base, d = line
                    ell_params.append(self.ellipse_angle((base[0]+lam*d[0], base[1]+lam*d[1])))
            for j, cl in enumerate(cutl):
                for lam in self.line_ellipse(cl):
                    cut_params[j].append(lam)
                    base, d = cl
                    ell_params.append(self.ellipse_angle((base[0]+lam*d[0], base[1]+lam*d[1])))
                for k, line in enumerate(lines):
                    x = self.line_line(line, cl)
                    if x is not None:
                        line_params[k].append(x[0])
                        cut_params[j].append(x[1])
        n = self.samples
        if P_present:
            for i, el in enumerate(els):
                ts = sorted(el_params[i])
                bounds = [mp.mpf(0)]+ts+[mp.mpf(1)]
                for a, b in zip(bounds, bounds[1:]):
                    if b <= a:
                        continue
                    mid = el.point((a+b)/2)
                    if not self.in_strip(mid):
                        continue
                    pts = [el.point(a+(b-a)*k/n) for k in range(n+1)] if n else None
                    if el.sign < 0 and pts:
                        pts.reverse()
                    self.pieces.append(Piece('P', el.green(a, b), self.in_B(mid), pts, ('el', i)))
            for k, line in enumerate(lines):
                self.line_pieces(line, line_params[k], line_toggles[k], 'P', ('strip', k))
        if ell is not None:
            E0, A1, A2, _, _ = ell
            angles = sorted(set(th % (2*mp.pi) for th in ell_params))
            toggles = sorted(th % (2*mp.pi) for th in el_toggles)
            if not angles:
                spans = [(mp.mpf(0), 2*mp.pi)]
            else:
                spans = [(a, b) for a, b in zip(angles, angles[1:])]+[(angles[-1], angles[0]+2*mp.pi)]
            first = None
            for a, b in spans:
                m = (a+b)/2
                p = self.point_at(m)
                if first is None:
                    first = (a, self.S.profile.inside(p[0], p[1]) if P_present else False)
                    inprof = first[1]
                else:
                    passed = sum(1 for t in toggles if first[0] < t <= a)
                    inprof = first[1] != (passed % 2 == 1)
                inP = inprof and self.in_strip(p)
                inH = self.in_cuts(p)
                pts = [self.point_at(a+(b-a)*k/n) for k in range(n+1)] if n else None
                angle = self.arc_measure(a, b) if inH else None
                self.pieces.append(Piece('B' if inH else 'E', ellipse_green(E0, A1, A2, a, b), inP, pts,
                                         ('ellipse',), angle=angle, flags=(inH, inP)))
            for j, cl in enumerate(cutl):
                self.line_pieces(cl, cut_params[j], cut_toggles[j], 'B', ('cut', j), skip=j)

    def arc_measure(self, a, b):
        """`|integral of D(theta)^-2 dtheta|` over the arc `[a, b]` of `t`
        (the cone's angle monotone along the section, which winds once about
        the axis)."""
        ta, tb = self.theta(a), self.theta(b)
        turn = 2*mp.pi
        if b-a >= turn-mp.mpf(10)**-30:
            delta = self.sense*turn
        else:
            delta = self.sense*((self.sense*(tb-ta)) % turn)
        return abs(self.way.inv_d2(ta+delta)-self.way.inv_d2(ta))

    def line_pieces(self, line, params, toggles, family, curve, skip=None):
        base, d = line
        params = sorted(params)
        toggles = sorted(toggles)
        at = lambda lam: (base[0]+lam*d[0], base[1]+lam*d[1])
        n = self.samples
        for a, b in zip(params, params[1:]):
            if b <= a:
                continue
            mid = at((a+b)/2)
            passed = sum(1 for t in toggles if t <= a)
            inprof = passed % 2 == 1
            pa, pb = at(a), at(b)
            pts = [at(a+(b-a)*k/n) for k in range(n+1)] if n else None
            if family == 'P':
                if not inprof:
                    continue
                self.pieces.append(Piece('P', seg_green(pa, pb), self.in_B(mid), pts, curve))
            else:
                if not (self.in_ellipse(mid) and self.in_cuts(mid, skip)):
                    continue
                inP = inprof and self.in_strip(mid) and self.strip_all is not False
                self.pieces.append(Piece('B', seg_green(pa, pb), inP, pts, curve))

    regions = sref.Slice.regions
    loops = sref.Slice.loops

    def angles(self):
        """`integral of D^-2` over the section's arcs on the cone's face,
        inside and outside the prism's section."""
        a_in = a_out = mp.mpf(0)
        for pc_ in self.pieces:
            if pc_.angle is None:
                continue
            if pc_.flags[1]:
                a_in += pc_.angle
            else:
                a_out += pc_.angle
        return a_in, a_out


# ------------------------------------------------------------------ events

def line_points(P0, d, cone):
    """The line `P0 + w d` (mpf) against the cone's quadric: its points."""
    Y = sub(P0, cone.om)
    k, b = cone.km, cone.bm
    du, dv, dw = dot(d, cone.xm), dot(d, cone.ym), dot(d, cone.nm)
    yu, yv, yw = dot(Y, cone.xm), dot(Y, cone.ym), dot(Y, cone.nm)
    rho = b+k*yw
    a = du*du+dv*dv-k*k*dw*dw
    bb = yu*du+yv*dv-rho*k*dw
    c = yu*yu+yv*yv-rho*rho
    out = []
    if a == 0:
        if bb != 0:
            out.append(-c/(2*bb))
    else:
        disc = bb*bb-a*c
        if disc >= 0:
            sq = mp.sqrt(disc)
            out += [(-bb-sq)/a, (-bb+sq)/a]
    return [add(P0, scale(d, w)) for w in out]


def circle_plane(c, nrm, r, m, p0):
    """The circle about `c` of radius `r` in the plane normal to `nrm`
    (unit) against the plane `m . (X - p0) = 0`: its points (mpf)."""
    L = cross(nrm, m)
    LL = dot(L, L)
    if LL == 0:
        return []
    a = dot(m, sub(p0, c))
    nn, nm_, mm = dot(nrm, nrm), dot(nrm, m), dot(m, m)
    det = nn*mm-nm_*nm_
    alpha, beta = (-nm_*a)/det, (nn*a)/det
    X0 = add(c, add(scale(nrm, alpha), scale(m, beta)))
    Q = sub(X0, c)
    QL, QQ = dot(Q, L), dot(Q, Q)-r*r
    disc = QL*QL-LL*QQ
    if disc <= 0:
        return []
    return [add(X0, scale(L, (-QL+sg*mp.sqrt(disc))/LL)) for sg in (-1, 1)]


def rim_wall_points(S, el, cone, e, rho_e, rows, q):
    """The rim at the cone's height `e` (radius `rho_e`) against the
    element's wall: [(t, w)]."""
    Pu, Pv, Pw = rows
    out = []
    if q[2] != 0:
        # The crease w = (e - P_w) / q_w in the rim's cylinder.
        wh = pc(psub([e], Pw), 1/q[2])
        fu = padd(Pu, pc(wh, q[0]))
        fv = padd(Pv, pc(wh, q[1]))
        poly = psub(padd(pmul(fu, fu), pmul(fv, fv)), [rho_e*rho_e])
        for t in exact_roots(poly):
            out.append((t, peval(pm(wh), t)))
    else:
        for t in exact_roots(psub(Pw, [e])):
            pu, pv = peval(pm(Pu), t), peval(pm(Pv), t)
            a = M(q[0]*q[0]+q[1]*q[1])
            bb = pu*M(q[0])+pv*M(q[1])
            c = pu*pu+pv*pv-M(rho_e)**2
            disc = bb*bb-a*c
            if a != 0 and disc > 0:
                for sg in (-1, 1):
                    out.append((t, (-bb+sg*mp.sqrt(disc))/a))
    return out


def events(pair, way):
    """The breakpoints of `way`'s slicing: [(s, tag, point)]."""
    S, cone = pair.S, pair.cone
    d, dm = way.d, way.dm
    out = []

    def add_point(tag, X):
        out.append((dot(dm, X), tag, X))
    los = (S.lo, S.hi)
    # The prism's vertices and cap edges' extremes.
    for v in S.profile.vertices:
        for w in los:
            add_point('vertex', S.world(M(v[0]), M(v[1]), M(w)))
    for el in S.elements:
        dp = padd(pc(el.dX, dot(d, S.x)), pc(el.dY, dot(d, S.y)))
        for t in exact_roots(dp):
            u, v = el.point(t)
            for w in los:
                add_point('cap_extreme', S.world(u, v, M(w)))
    # The apex, real or virtual.
    add_point('apex', Mv(cone.apex))
    # The rims' extremes.
    perp = sub(dm, scale(cone.nm, dot(dm, cone.nm)))
    pl = mnorm(perp)
    for e, r in cone.ends:
        c = add(cone.om, scale(cone.nm, M(e)))
        if r > 0 and pl > 0:
            for sg in (-1, 1):
                add_point('rim', add(c, scale(perp, sg*M(r)/pl)))
        elif r > 0:
            add_point('rim', c)
    # The cap planes' conics' extremes: on the cap plane, the gradient in
    # span(m, d).
    m = S.m
    ell = cross(m, d)
    if dot(ell, ell) != 0:
        kk = 1+cone.k*cone.k
        g = sub(ell, scale(cone.n, kk*dot(ell, cone.n)))
        rhs = cone.b*cone.k*dot(ell, cone.n)+dot(g, cone.o)
        L = cross(m, g)
        if dot(L, L) != 0:
            for w in los:
                p0 = add(S.o, scale(S.n, w))
                # A point on both planes m . X = m . p0, g . X = rhs.
                a1, a2 = dot(m, p0), rhs
                mm, mg, gg = dot(m, m), dot(m, g), dot(g, g)
                det = mm*gg-mg*mg
                X1 = add(scale(m, (a1*gg-a2*mg)/det), scale(g, (a2*mm-a1*mg)/det))
                for X in line_points(Mv(X1), Mv(L), cone):
                    add_point('cap_conic', X)
    # Vertical edges against the quadric and the end planes.
    qn = dot(cone.n, S.n)
    for v in S.profile.vertices:
        P0 = add(S.o, add(scale(S.x, v[0]), scale(S.y, v[1])))
        for X in line_points(Mv(P0), S.nm, cone):
            add_point('edge_cone', X)
        if qn != 0:
            for e, _ in cone.ends:
                w = (e-dot(cone.n, sub(P0, cone.o)))/qn
                add_point('edge_end', add(Mv(P0), scale(S.nm, M(w))))
    # The rims against the cap planes.
    for e, r in cone.ends:
        if r == 0:
            continue
        c = add(cone.om, scale(cone.nm, M(e)))
        for w in los:
            for X in circle_plane(c, cone.nm, M(r), Mv(S.m), Mv(add(S.o, scale(S.n, w)))):
                add_point('rim_cap', X)
    for el in S.elements:
        P, dP = S.wall(el)
        Aq, B, C, rows, q = cone_quadratic(S, el, cone)
        # Cap edges against the quadric and the end planes.
        for w in los:
            for t in exact_roots(padd(padd([Aq*w*w], pc(B, 2*w)), C)):
                u, v = el.point(t)
                add_point('cap_cone', S.world(u, v, M(w)))
            for e, _ in cone.ends:
                for t in exact_roots(padd(rows[2], [q[2]*w-e])):
                    u, v = el.point(t)
                    add_point('cap_end', S.world(u, v, M(w)))
        # The rims against the wall.
        for e, r in cone.ends:
            if r == 0:
                continue
            for t, w in rim_wall_points(S, el, cone, e, r, rows, q):
                u, v = el.point(t)
                add_point('rim_wall', S.world(u, v, w))
        # The meeting's extremes.
        dP_ = vpoly_dot(d, dP)
        dn = dot(d, S.n)
        dB, dC = pder(B), pder(C)
        N = psub(pmul(B, dP_), pc(dC, dn/2))
        Dd = psub(pc(dP_, Aq), pc(dB, dn))
        poly = padd(psub(pc(pmul(N, N), Aq), pc(pmul(pmul(B, N), Dd), 2)), pmul(C, pmul(Dd, Dd)))
        if pdeg(poly) > 0:
            for t in exact_roots(poly):
                den = peval(pm(Dd), t)
                u, v = el.point(t)
                if abs(den) > mp.mpf(10)**-20:
                    add_point('meeting_extreme', S.world(u, v, -peval(pm(N), t)/den))
                    continue
                b, cc = peval(pm(B), t), peval(pm(C), t)
                if Aq == 0:
                    if b != 0:
                        add_point('meeting_extreme', S.world(u, v, -cc/(2*b)))
                    continue
                disc = b*b-M(Aq)*cc
                if disc >= 0:
                    for sg in (-1, 1):
                        add_point('meeting_extreme', S.world(u, v, (-b+sg*mp.sqrt(disc))/M(Aq)))
        # The end planes' creases' extremes.
        if q[2] != 0:
            for e, _ in cone.ends:
                wh = pc(psub([e], rows[2]), 1/q[2])
                for t in exact_roots(padd(vpoly_dot(d, dP), pc(pder(wh), dot(d, S.n)))):
                    u, v = el.point(t)
                    add_point('crease_extreme', S.world(u, v, peval(pm(wh), t)))
    return out


# ------------------------------------------------------------------ the pair

def case_size(S, cone):
    us = [v[0] for v in S.profile.vertices]
    vs = [v[1] for v in S.profile.vertices]
    for el in S.elements:
        us += [p[0] for p in el.ctrl]
        vs += [p[1] for p in el.ctrl]
    pts = [S.world(M(u), M(v), M(w)) for u in (min(us), max(us)) for v in (min(vs), max(vs))
           for w in (S.lo, S.hi)]
    for e, r in cone.ends:
        c = add(cone.om, scale(cone.nm, M(e)))
        for sx in (-1, 1):
            for sy in (-1, 1):
                pts.append(add(c, add(scale(cone.xm, sx*M(r)), scale(cone.ym, sy*M(r)))))
    lo = [min(p[i] for p in pts) for i in range(3)]
    hi = [max(p[i] for p in pts) for i in range(3)]
    return max(h-l for l, h in zip(lo, hi)), lo, hi


class Pair:
    def __init__(self, obj, tool):
        self.obj_prism = obj.cone is None
        prism, cone = (obj, tool) if self.obj_prism else (tool, obj)
        assert prism.cone is None and cone.cone is not None
        self.S, self.cone = sref.Spline3(prism), Cone(cone)
        self.size, self.lo, self.hi = case_size(self.S, self.cone)
        d1 = first_direction(self.S, self.cone)
        self.first = Way(self, d1)
        self.second = Way(self, second_direction(self.S, self.cone, d1))
        self.A = axis_term(self.S, self.cone)
        self._sliced = {}
        self._faces = None
        self._solids = None
        self.stats = {'breaks': 0}

    breakpoints = sref.Pair.breakpoints

    def sliced(self, way):
        """The integrals over `way`'s slices: per region (P, B, C) its volume
        and world moments, and (the first way) the cone's wall's area inside
        and outside the prism."""
        key = id(way)
        if key in self._sliced:
            return self._sliced[key]
        breaks, _ = self.breakpoints(way)
        self.stats['breaks'] = max(self.stats['breaks'], len(breaks))
        wall = way is self.first
        cone = self.cone

        def f(s):
            sl = Slice(way, s)
            reg = sl.regions()
            row = []
            for k in 'PBC':
                row += way.world_moments(s, reg[k])
            if wall:
                g = abs(way.gb+way.gk*s)*cone.N
                a_in, a_out = sl.angles()
                row += [g*a_in, g*a_out]
            return row
        n = 14 if wall else 12
        total = [mp.mpf(0)]*n
        err = mp.mpf(0)
        tol = mp.mpf(10)**-33*self.size**4
        for a, b in zip(breaks, breaks[1:]):
            est, e = cref.integrate(f, a, b, tol)
            if est is not None:
                total = [x+y for x, y in zip(total, est)]
                err += e
        J = way.J
        out = {}
        for i, k in enumerate('PBC'):
            out[k] = (J*total[4*i], tuple(J*total[4*i+1+j] for j in range(3)))
        if wall:
            # `ds` along `d` per unit of the slices' parameter: the area
            # element is per `ds`, `s = d . X` (no Jacobian).
            out['wall'] = (total[12], total[13])
        out['error'] = err
        out['breaks'] = breaks
        self._sliced[key] = out
        return out

    volumes = sref.Pair.volumes

    def faces(self):
        """[(input tag 'S' or 'K', face name, {'in': area, 'out': area},
        closed area)]."""
        if self._faces is not None:
            return self._faces
        S, cone = self.S, self.cone
        out = []
        mm = mnorm(S.m)
        prof_area = S.profile.moments()[0]
        for name, w in (('cap_lo', S.lo), ('cap_hi', S.hi)):
            inside, whole = cap_classes(self, w)
            out.append(('S', name, {'in': inside*mm, 'out': (whole-inside)*mm}, prof_area*mm))
        for i, el in enumerate(S.elements):
            cls = wall_classes(self, el)
            out.append(('S', f'wall{i}', cls, M(S.hi-S.lo)*wall_length(S, el)))
        a_in, a_out = self.sliced(self.first)['wall']
        out.append(('K', 'wall', {'in': a_in, 'out': a_out}, cone.lateral()))
        for idx, (e, r) in enumerate(cone.ends):
            if r > 0:
                out.append(('K', f'disc{idx}', disc_classes(self, e, r), mp.pi*M(r)**2))
        self._faces = out
        return out

    def area(self, op):
        keep = {'fuse': ('out', 'out'), 'common': ('in', 'in')}
        if op == 'cut':
            keep['cut'] = ('out', 'in') if self.obj_prism else ('in', 'out')
        total = mp.mpf(0)
        for tag, _, cls, _ in self.faces():
            total += cls[keep[op][0 if tag == 'S' else 1]]
        return total

    def input_area(self, tag):
        return sum((c for t, _, _, c in self.faces() if t == tag), mp.mpf(0))

    def solids(self):
        if self._solids is None:
            self._solids = count_solids(self)
        return self._solids

    result = sref.Pair.result


def wall_classes(pair, el):
    """A wall's area inside and outside the cone."""
    S, cone = pair.S, pair.cone
    A, B, C, rows, q = cone_quadratic(S, el, cone)
    lo, hi = S.lo, S.hi
    polys = [psub(pmul(B, B), pc(C, A))]
    for w in (lo, hi):
        polys.append(padd(padd([A*w*w], pc(B, 2*w)), C))
    for e, r in cone.ends:
        if q[2] != 0:
            wh = pc(psub([e], rows[2]), 1/q[2])
            polys.append(padd(padd(pc(pmul(wh, wh), A), pc(pmul(B, wh), 2)), C))
            for w in (lo, hi):
                polys.append(psub(wh, [w]))
        else:
            polys.append(psub(rows[2], [e]))
    if A == 0:
        polys.append(B)
    ts = [mp.mpf(0), mp.mpf(1)]
    for p in polys:
        if pdeg(p) > 0:
            ts += exact_roots(p, slack=mp.mpf(0))
    breaks = cref.merge_breaks(ts, mp.mpf(0), mp.mpf(1), mp.mpf(10)**-30)
    Am, Bm, Cm, Pw = M(A), pm(B), pm(C), pm(rows[2])
    lom, him, h = M(lo), M(hi), M(hi-lo)
    qn, hc = M(q[2]), cone.hm

    def f(t):
        du, dv = el.deriv(t)
        T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
        g = mnorm(cross(T, S.nm))
        b, cc = peval(Bm, t), peval(Cm, t)
        # The ruling's part inside the double cone, as intervals of w.
        if Am == 0:
            if b == 0:
                ivals = [] if cc > 0 else [(lom, him)]
            else:
                w0 = -cc/(2*b)
                ivals = [(lom, w0)] if b > 0 else [(w0, him)]
        else:
            disc = b*b-Am*cc
            if disc <= 0:
                ivals = [] if Am > 0 else [(lom, him)]
            else:
                sq = mp.sqrt(disc)
                w1, w2 = sorted(((-b-sq)/Am, (-b+sq)/Am))
                ivals = [(w1, w2)] if Am > 0 else [(lom, w1), (w2, him)]
        # Within the slab and the heights.
        pw = peval(Pw, t)
        if qn == 0:
            slab = (lom, him) if 0 < pw < hc else (lom, lom)
        else:
            a1, a2 = sorted(((0-pw)/qn, (hc-pw)/qn))
            slab = (a1, a2)
        inside = mp.mpf(0)
        for a, z in ivals:
            a, z = max(a, lom, slab[0]), min(z, him, slab[1])
            inside += max(z-a, mp.mpf(0))
        return [inside*g, (h-inside)*g]
    total = [mp.mpf(0), mp.mpf(0)]
    for a, b in zip(breaks, breaks[1:]):
        est, _ = cref.integrate(f, a, b, mp.mpf(10)**-34*pair.size**2)
        if est is not None:
            total = [x+y for x, y in zip(total, est)]
    return {'in': total[0], 'out': total[1]}


def cap_classes(pair, w):
    """A cap's area in `(u, v)` inside the cone, and its whole area, both by
    chords along `CHORD_ANGLE`."""
    S, cone = pair.S, pair.cone
    dlt = (mp.cos(CHORD_ANGLE), mp.sin(CHORD_ANGLE))
    nrm = (-dlt[1], dlt[0])
    wm = M(w)
    # The cone's quadric on the cap: Q(u, v) with X = o + u x + v y + w n.
    base = add(S.om, scale(S.nm, wm))
    Yb = sub(base, cone.om)
    rowsx = (dot(S.xm, cone.xm), dot(S.xm, cone.ym), dot(S.xm, cone.nm))
    rowsy = (dot(S.ym, cone.xm), dot(S.ym, cone.ym), dot(S.ym, cone.nm))
    y0 = (dot(Yb, cone.xm), dot(Yb, cone.ym), dot(Yb, cone.nm))
    k, b, hc = cone.km, cone.bm, cone.hm

    def local(p):
        return tuple(y0[i]+p[0]*rowsx[i]+p[1]*rowsy[i] for i in range(3))

    def along(p0, dd):
        """The quadric and the cone's height along `p0 + lam dd`."""
        c0, c1 = local(p0), tuple(dd[0]*rowsx[i]+dd[1]*rowsy[i] for i in range(3))
        r0, r1 = b+k*c0[2], k*c1[2]
        a2 = c1[0]**2+c1[1]**2-r1*r1
        a1 = 2*(c0[0]*c1[0]+c0[1]*c1[1]-r0*r1)
        a0 = c0[0]**2+c0[1]**2-r0*r0
        return (a2, a1, a0), (c0[2], c1[2])

    def chord(tau):
        p0 = (tau*nrm[0], tau*nrm[1])
        prof = [(a[0], z[0]) for a, z in S.profile.chords(p0, dlt)]
        (a2, a1, a0), (h0, h1) = along(p0, dlt)
        # Inside the quadric.
        if a2 == 0:
            ivals = [] if a1 == 0 else ([(-mp.inf, -a0/a1)] if a1 > 0 else [(-a0/a1, mp.inf)])
        else:
            disc = a1*a1-4*a2*a0
            if disc <= 0:
                ivals = [] if a2 > 0 else [(-mp.inf, mp.inf)]
            else:
                sq = mp.sqrt(disc)
                l1, l2 = sorted(((-a1-sq)/(2*a2), (-a1+sq)/(2*a2)))
                ivals = [(l1, l2)] if a2 > 0 else [(-mp.inf, l1), (l2, mp.inf)]
        if h1 == 0:
            slab = (-mp.inf, mp.inf) if 0 < h0 < hc else (mp.mpf(0), mp.mpf(0))
        else:
            slab = tuple(sorted(((0-h0)/h1, (hc-h0)/h1)))
        inside = mp.mpf(0)
        whole = mp.mpf(0)
        for p, z in prof:
            whole += z-p
            for a, e in ivals:
                lo_, hi_ = max(p, a, slab[0]), min(z, e, slab[1])
                inside += max(hi_-lo_, mp.mpf(0))
        return [inside, whole]
    # Breakpoints: points of the cap's 2D arrangement and the curves'
    # tangents along the chords, by their `tau = p . nrm`.
    taus = []
    for v in S.profile.vm:
        taus.append(v[0]*nrm[0]+v[1]*nrm[1])
    for el in S.elements:
        if el.spline:
            # Where the span runs along the chords: S'(t) x delta.
            n = max(len(el.dXm), len(el.dYm))
            dX = el.dXm+[0]*(n-len(el.dXm))
            dY = el.dYm+[0]*(n-len(el.dYm))
            for t in sref.roots01([x*dlt[1]-y*dlt[0] for x, y in zip(dX, dY)]):
                u, v = el.point(t)
                taus.append(u*nrm[0]+v*nrm[1])
    # The conic: Q(p) = p^T Mq p + 2 g . p + q0 in (u, v).
    c0, cu, cv = local((0, 0)), rowsx, rowsy
    r0 = b+k*c0[2]
    Mq = [[cu[0]**2+cu[1]**2-k*k*cu[2]**2, cu[0]*cv[0]+cu[1]*cv[1]-k*k*cu[2]*cv[2]], [0, 0]]
    Mq[1] = [Mq[0][1], cv[0]**2+cv[1]**2-k*k*cv[2]**2]
    gq = (c0[0]*cu[0]+c0[1]*cu[1]-r0*k*cu[2], c0[0]*cv[0]+c0[1]*cv[1]-r0*k*cv[2])
    q0 = c0[0]**2+c0[1]**2-r0*r0
    # Tangent along the chords: the discriminant of the quadratic in lam as a
    # function of tau (quadratic).
    dd = dlt
    a2 = dd[0]*(Mq[0][0]*dd[0]+Mq[0][1]*dd[1])+dd[1]*(Mq[1][0]*dd[0]+Mq[1][1]*dd[1])
    mn = (Mq[0][0]*nrm[0]+Mq[0][1]*nrm[1], Mq[1][0]*nrm[0]+Mq[1][1]*nrm[1])
    # a1(tau)/2 = tau (d . M nrm) + g . d; a0 = tau^2 nrm M nrm + 2 tau g . nrm + q0.
    b1, b0 = dd[0]*mn[0]+dd[1]*mn[1], gq[0]*dd[0]+gq[1]*dd[1]
    c2, c1_, c0_ = nrm[0]*mn[0]+nrm[1]*mn[1], 2*(gq[0]*nrm[0]+gq[1]*nrm[1]), q0
    # (a1/2)^2 - a2 a0 = 0.
    poly = [b0*b0-a2*c0_, 2*b1*b0-a2*c1_, b1*b1-a2*c2]
    taus += [x for x in quadratic_roots(poly)]
    # The slab's lines on the cap: w_c = 0, h along the chord direction.
    slab_lines = []
    if cu[2] != 0 or cv[2] != 0:
        for e in (mp.mpf(0), hc):
            g2 = cu[2]**2+cv[2]**2
            bse = ((e-c0[2])*cu[2]/g2, (e-c0[2])*cv[2]/g2)
            slab_lines.append((bse, (cv[2], -cu[2])))
    for el in S.elements:
        # The conic's crossings of the element: Q(S(t)), degree 2 p.
        n = el.deg+1
        X = el.Xm+[0]*(n-len(el.Xm))
        Y = el.Ym+[0]*(n-len(el.Ym))
        poly = [mp.mpf(0)]*(2*n-1)
        for i in range(n):
            for j in range(n):
                poly[i+j] += X[i]*X[j]*Mq[0][0]+2*X[i]*Y[j]*Mq[0][1]+Y[i]*Y[j]*Mq[1][1]
            poly[i] += 2*(gq[0]*X[i]+gq[1]*Y[i])
        poly[0] += q0
        for t in sref.roots01(poly):
            u, v = el.point(t)
            taus.append(u*nrm[0]+v*nrm[1])
        for line in slab_lines:
            bse, dl = line
            if el.misses_line(bse, dl):
                continue
            for t in sref.roots01(el.line_poly(bse, dl)):
                u, v = el.point(t)
                taus.append(u*nrm[0]+v*nrm[1])
    for bse, dl in slab_lines:
        # The slab's line against the conic.
        aa = dl[0]*(Mq[0][0]*dl[0]+Mq[0][1]*dl[1])+dl[1]*(Mq[1][0]*dl[0]+Mq[1][1]*dl[1])
        bb = 2*(dl[0]*(Mq[0][0]*bse[0]+Mq[0][1]*bse[1])+dl[1]*(Mq[1][0]*bse[0]+Mq[1][1]*bse[1])
                + gq[0]*dl[0]+gq[1]*dl[1])
        cc = bse[0]*(Mq[0][0]*bse[0]+Mq[0][1]*bse[1])+bse[1]*(Mq[1][0]*bse[0]+Mq[1][1]*bse[1]) \
            + 2*(gq[0]*bse[0]+gq[1]*bse[1])+q0
        for lam in quadratic_roots([cc, bb, aa]):
            p = (bse[0]+lam*dl[0], bse[1]+lam*dl[1])
            taus.append(p[0]*nrm[0]+p[1]*nrm[1])
    prof_taus = [v[0]*nrm[0]+v[1]*nrm[1] for v in S.profile.vm]
    for el in S.elements:
        if el.spline:
            prof_taus += [p[0]*nrm[0]+p[1]*nrm[1] for p in el.ctrlm]
    lo_t, hi_t = min(prof_taus)-1, max(prof_taus)+1
    breaks = cref.merge_breaks(taus, lo_t, hi_t, mp.mpf(10)**-30*pair.size)
    total = [mp.mpf(0), mp.mpf(0)]
    for a, z in zip(breaks, breaks[1:]):
        est, _ = cref.integrate(chord, a, z, mp.mpf(10)**-34*pair.size**2)
        if est is not None:
            total = [x+y for x, y in zip(total, est)]
    return total[0], total[1]


def quadratic_roots(c):
    """The real roots of `c0 + c1 x + c2 x^2` (mpf), a double one once."""
    c0, c1, c2 = c
    if c2 == 0:
        return [] if c1 == 0 else [-c0/c1]
    disc = c1*c1-4*c2*c0
    if disc < 0:
        return []
    sq = mp.sqrt(disc)
    return [(-c1-sq)/(2*c2), (-c1+sq)/(2*c2)]


def disc_classes(pair, e, r):
    """An end disc (at the cone's height `e`, radius `r`) inside and outside
    the prism."""
    S, cone = pair.S, pair.cone
    n = cone.n
    c = add(cone.o, scale(n, e))
    if dot(n, S.n) != 0:
        d = n if dot(n, S.n) > 0 else scale(n, -1)
        way = Way(pair, d)
        sl = Slice(way, way.s_of(Mv(c)), cuts=False)
        reg = sl.regions()
        k = way.area_factor
        return {'in': reg['C'][0]*k, 'out': (reg['B'][0]-reg['C'][0])*k}
    # The plane holds the prism's axis: the profile's chords along its trace
    # times the heights (S9f.3a's hemisphere's disc).
    ev = cross(n, S.n)
    lc = apply(S.inv, sub(c, S.o))
    le = apply(S.inv, ev)
    chords = S.profile.chords(Mv(lc[:2]), Mv(le[:2]))
    lo, hi = M(S.lo), M(S.hi)
    wc, ez = M(lc[2]), M(le[2])
    ee, en, nn = M(dot(ev, ev)), M(dot(ev, S.n)), M(dot(S.n, S.n))
    r2 = M(r)**2
    inside = mp.mpf(0)

    def width(lam):
        disc = (en*lam)**2-nn*(ee*lam*lam-r2)
        if disc <= 0:
            return [mp.mpf(0)]
        sq = mp.sqrt(disc)
        m0, m1 = (-en*lam-sq)/nn, (-en*lam+sq)/nn
        a, b = max(m0, lo-wc-lam*ez), min(m1, hi-wc-lam*ez)
        return [max(b-a, mp.mpf(0))]
    lim = M(r)*mp.sqrt(nn/(ee*nn-en*en))
    for (t0, _), (t1, _) in chords:
        pts = [t0, t1, -lim, lim]
        for kk in (lo, hi):
            m0 = kk-wc
            qa = ee-2*en*ez+nn*ez*ez
            qb = 2*en*m0-2*nn*m0*ez
            qc = nn*m0*m0-r2
            disc = qb*qb-4*qa*qc
            if disc > 0:
                pts += [(-qb+sg*mp.sqrt(disc))/(2*qa) for sg in (-1, 1)]
        brk = cref.merge_breaks(pts, t0, t1, mp.mpf(10)**-30)
        for a, b in zip(brk, brk[1:]):
            est, _ = cref.integrate(width, a, b, mp.mpf(10)**-34*pair.size**2)
            if est is not None:
                inside += est[0]
    k = mnorm(cross(ev, S.n))
    inside *= k
    return {'in': inside, 'out': mp.pi*M(r)**2-inside}


# ------------------------------------------------------------------ solids

GRID = 320


def count_solids(pair):
    """The solids of each operation (combinatorics only, at 20 digits): the
    pair's box cut into `GRID` slices across the prism's axis (`w` the
    prism's height) and each slice into `GRID` chords along `CHORD_ANGLE` in
    the profile's plane; on each chord the prism's intervals (the profile's
    chords, within its heights) and the cone's (its quadric's chord within
    its slab) exactly, the operation's intervals from them; intervals joined
    where they overlap on neighbouring chords of a slice and on one chord of
    neighbouring slices, the solids the classes. (Unlike S9f.3a's sections,
    whose overlaps in the profile's plane miss a thin layer under a cap that
    a slice leaning from the prism's axis crosses as a moving sliver, this
    needs no section of the cone by a plane of any kind.)"""
    S, cone = pair.S, pair.cone
    counts = {}
    with mp.workdps(20):
        dlt = (mp.cos(CHORD_ANGLE), mp.sin(CHORD_ANGLE))
        nrm = (-dlt[1], dlt[0])
        # The extent of both in the profile's plane and along the axis.
        corners = [S.local(Mv(add(cone.o, add(scale(cone.n, e), add(scale(cone.x, sx*r), scale(cone.y, sy*r))))))
                   for e, r in cone.ends for sx in (-1, 1) for sy in (-1, 1)]
        taus = [p[0]*nrm[0]+p[1]*nrm[1] for p in corners]
        ws = [p[2] for p in corners]+[M(S.lo), M(S.hi)]
        for v in S.profile.vm:
            taus.append(v[0]*nrm[0]+v[1]*nrm[1])
        for el in S.elements:
            taus += [q[0]*nrm[0]+q[1]*nrm[1] for q in el.ctrlm]
        t0, t1 = min(taus), max(taus)
        w0, w1 = min(ws), max(ws)
        n = GRID
        tau = [t0+(t1-t0)*(i+mp.mpf(1)/2)/n for i in range(n)]
        hs = [w0+(w1-w0)*(j+mp.mpf(1)/2)/n for j in range(n)]
        prof = [[(a[0], z[0]) for a, z in S.profile.chords((t*nrm[0], t*nrm[1]), dlt)] for t in tau]
        lo, hi = M(S.lo), M(S.hi)
        k, b, hc = cone.km, cone.bm, cone.hm

        def cone_ivals(t, w):
            base = S.world(t*nrm[0], t*nrm[1], w)
            Y = sub(base, cone.om)
            d3 = add(scale(S.xm, dlt[0]), scale(S.ym, dlt[1]))
            c0 = (dot(Y, cone.xm), dot(Y, cone.ym), dot(Y, cone.nm))
            c1 = (dot(d3, cone.xm), dot(d3, cone.ym), dot(d3, cone.nm))
            r0, r1 = b+k*c0[2], k*c1[2]
            a2 = c1[0]**2+c1[1]**2-r1*r1
            a1 = 2*(c0[0]*c1[0]+c0[1]*c1[1]-r0*r1)
            a0 = c0[0]**2+c0[1]**2-r0*r0
            disc = a1*a1-4*a2*a0
            if a2 == 0 or disc <= 0:
                return []
            sq = mp.sqrt(disc)
            l1, l2 = sorted(((-a1-sq)/(2*a2), (-a1+sq)/(2*a2)))
            ivals = [(l1, l2)] if a2 > 0 else [(-mp.inf, l1), (l2, mp.inf)]
            if c1[2] == 0:
                slab = (-mp.inf, mp.inf) if 0 < c0[2] < hc else None
            else:
                slab = tuple(sorted(((0-c0[2])/c1[2], (hc-c0[2])/c1[2])))
            if slab is None:
                return []
            out = []
            for x, y in ivals:
                x, y = max(x, slab[0]), min(y, slab[1])
                if y > x:
                    out.append((x, y))
            return out

        def combine(op, P, K):
            """The operation's intervals from the prism's `P` and the cone's
            `K` (each sorted and disjoint)."""
            pts = sorted(set([x for i in P+K for x in i]))
            out = []
            for x, y in zip(pts, pts[1:]):
                m = (x+y)/2
                inP = any(a <= m <= z for a, z in P)
                inK = any(a <= m <= z for a, z in K)
                obj, tool = (inP, inK) if pair.obj_prism else (inK, inP)
                keep = {'fuse': obj or tool, 'cut': obj and not tool, 'common': obj and tool}[op]
                if keep:
                    if out and out[-1][1] == x:
                        out[-1] = (out[-1][0], y)
                    else:
                        out.append((x, y))
            return out
        grid = []
        for j, w in enumerate(hs):
            inside = lo < w < hi
            row = []
            for i, t in enumerate(tau):
                row.append((prof[i] if inside else [], cone_ivals(t, w)))
            grid.append(row)
        for op in OPS:
            parent = []

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            cells = []
            for j in range(n):
                cells.append([])
                for i in range(n):
                    ivs = combine(op, *grid[j][i])
                    ids = []
                    for iv in ivs:
                        parent.append(len(parent))
                        ids.append((len(parent)-1, iv))
                    cells[j].append(ids)

            def link(A, B):
                for x, (a0, a1) in A:
                    for y, (b0, b1) in B:
                        if min(a1, b1) > max(a0, b0):
                            parent[find(x)] = find(y)
            for j in range(n):
                for i in range(n):
                    if i+1 < n:
                        link(cells[j][i], cells[j][i+1])
                    if j+1 < n:
                        link(cells[j][i], cells[j+1][i])
            counts[op] = len({find(x) for x in range(len(parent))})
    return counts


# ------------------------------------------------------------------ rows

def number(x):
    return cref.number(x)


def rows(obj, operation, tool, pair=None):
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair

