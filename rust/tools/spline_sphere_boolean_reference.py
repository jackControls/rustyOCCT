#!/usr/bin/env python3
"""Independent reference for S9f.3a of REVIEW_NOTES.md: Booleans of a spline
prism (a profile of lines and S8b's nonrational spline segments) against a
sphere or a hemisphere (`Solid::sphere_with`, latitudes -pi/2, 0 and pi/2:
heights `R sin(latitude)` exact) in any relative position, either the
object: every spline wall meets the sphere in a curve that is, along the
wall's ruling at the spline's parameter `t`, a root of `A w^2 + 2 B(t) w +
C(t)` (`A = n . n`).

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal). The prism: the points `o + u x + v y + w n` with `(u, v)` in
the profile (`curved_boolean_reference.Profile`: the stored points, each
spline segment cut into its exact Bezier spans; its parsing, ray test and
chords, the quadrature and the root finders are what is shared) and `w`
between the offsets. The sphere:
`|X - c|^2 <= r^2` about its frame's stored origin, a hemisphere also
`s n . (X - c) >= 0` (`s` the side, `n` the stored normal; hemispheres only
in exact frames, where `n = x * y`). Nothing here uses the kernel, a
surface/surface intersection or an arrangement.

* **Slicing (volume, first moments, solids, the sphere's face).** Both
  solids are sliced by the planes `d . X = s` (`d . n > 0`). Every point of
  a slice is projected along the prism's axis onto its profile plane: `X =
  o + u x + v y + w(u, v) n` with `w` affine in `(u, v)`, so in `(u, v)` the
  prism's section is the profile cut by the two lines `w = lo`, `w = hi`
  (a strip; when `d` is normal to the caps, all of the profile or none),
  and the sphere's the projection of the section's circle, an ellipse
  `E0 + cos(theta) e1 + sin(theta) e2` (`theta` the circle's own angle in
  an orthonormal basis of the plane), cut by a hemisphere's line. The
  four regions are taken apart by their boundaries: each curve (profile
  elements, strip lines, the ellipse, the hemisphere's line) cut at its
  crossings with the other region's curves (the real roots of `|K (S(t) -
  E0)|^2 - 1`, degree `2 p` on a span; a line's `(S(t) - b) x d`, degree
  `p`; quadratics between lines and the ellipse) and classified at its
  pieces' midpoints (inside the ellipse and the line directly, inside the
  profile by parity: along a line from its far end, along the ellipse from
  one ray test); common = `dP n B + dB n P`, `P - B = P - common`; areas
  and first moments in `(u, v)` by Green's theorem on every piece in
  closed form (segments; spans exactly by their antiderivatives; the
  ellipse's `u dv - v du` a constant plus a first harmonic). Volumes are
  `det(x, y, n) / (d . n)` times the integral over `s`. The integrands are
  analytic between breakpoints, every one the slice of an arrangement's
  vertex or an extreme of `d . X` along one of its edges: the prism's
  vertices; its cap edges' extremes (`d . S'(t)`, degree `p - 1`); the
  sphere's poles along `d`; its sections by the cap planes and the
  hemisphere's rim (their extremes); vertical edges, cap edges and the
  rim against the sphere and the hemisphere's plane (quadratics, `F(t,
  w)` at a cap's height, degree `2 p`, the crease `w = h(t)` in `F`,
  degree `2 p`, or on a plane holding the axis its generatrices, degree
  `p`); the meeting's extremes, where `F = 0` and `d . (F_w S' - F_t n) =
  0`, linear in `w` (eliminated: `A N^2 - 2 B N D + C D^2`, degree `4 p -
  2`); the hemisphere's crease's extremes. Each interval is integrated by
  `curved_boolean_reference.integrate` (Gauss-Legendre after `s = a + (b
  - a)(1 - cos t)/2`, doubling until two estimates agree within 1e-33 of
  the case's size, else halving). Two directions: the caps' normal (the
  rows' values) and an oblique one, independently.
* **Surface area.** The prism's walls swept along their generatrices:
  along each, the part inside the sphere (the roots `w-`, `w+` of `F`,
  cut by the hemisphere's plane and the caps) in closed form, times `|S'(t)
  x n|`, integrated over `t` between the roots of `B^2 - A C` (turning
  points), of `F` at the caps' and the crease's heights and of the
  hemisphere's line. The caps: their profile against the sphere's
  section at their height (the slicing's regions). The sphere's face by
  Archimedes: its area between the slices `z` and `z + dz` along the
  caps' normal is `r dz dtheta`, so its part inside the prism is `r / |d|`
  times the integral of the angle of the section's circle inside the
  prism's section (the ellipse's classified pieces). A hemisphere's disc
  against the prism's section in its plane: the projection along `n`
  where the plane is transverse to it, else the profile's chords along the
  plane's trace times the heights, as rectangles.
* **Solids** (combinatorics only, binary64): the slices along the caps'
  normal at five places in every interval between breakpoints (two of
  them within 1e-9 of its ends), each result's region as loops of its
  kept pieces, components outer loops with their holes; components of
  consecutive slices joined when an interior point of one lies in the
  other.

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9f
references.
"""
from fractions import Fraction as F
import math

import mpmath as mp

from curve_surface_reference import stored_axes
import curved_boolean_reference as cref
from boolean_reference import pder, pint, pmul, psub, padd, peval, ptrim, pdeg

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
# Latitudes are compared with the stored binary64 values only.
HALF_PI = math.pi/2
# An oblique second slicing direction (its `n` component positive).
SECOND = (F(2), F(-3), F(5))

M, Mv = cref.M, cref.Mv
sub, add, scale, dot, cross, det3, apply = cref.sub, cref.add, cref.scale, cref.dot, cref.cross, \
    cref.det3, cref.apply


def pc(p, k):
    return ptrim([x*k for x in p]) if k else [F(0)]


def pm(p):
    return [M(c) for c in p]


def mnorm(v):
    return mp.sqrt(sum(M(c)**2 for c in v))


def roots01(c):
    """The real roots of an mpf polynomial in [0, 1] (monotone runs)."""
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    if len(c) <= 1:
        return []
    return cref._roots_in(c, mp.mpf(0), mp.mpf(1), None)


def exact_roots(p, lo=F(0), hi=F(1), slack=mp.mpf(10)**-25):
    """The real roots of an exact polynomial in [lo, hi] (mpf)."""
    if pdeg(p) <= 0:
        return []
    return [t for t in cref.real_roots(p) if M(lo)-slack <= t <= M(hi)+slack]


# ------------------------------------------------------------------ inputs

class Ball:
    """A sphere or a hemisphere on its exact model."""

    def __init__(self, case):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        radius, low, high = case.sphere
        self.case = case
        self.c, self.r = o, F(radius)
        self.r2 = self.r*self.r
        assert low in (-HALF_PI, 0.0) and high in (0.0, HALF_PI) and low < high, (low, high)
        self.half = 0 if (low, high) == (-HALF_PI, HALF_PI) else (1 if low == 0.0 else -1)
        self.ns = None
        if self.half:
            assert cross(x, y) == n and dot(n, n) == 1, 'a hemisphere in an exact frame'
            self.ns = scale(n, self.half)
        self.cm, self.rm = Mv(o), M(radius)

    def closed(self):
        """Volume, first moments and area in closed form."""
        r = self.rm
        if not self.half:
            V = 4*mp.pi*r**3/3
            return V, scale(self.cm, V), 4*mp.pi*r*r
        V = 2*mp.pi*r**3/3
        centre = add(self.cm, scale(Mv(self.ns), 3*r/8))
        return V, scale(centre, V), 3*mp.pi*r*r

    def contains(self, X):
        d = sub(X, self.cm)
        if dot(d, d) >= self.rm**2:
            return False
        return not self.half or dot(d, Mv(self.ns)) > 0


class Element:
    """A profile element of lines and spline spans, exact: `S(t) = (X(t),
    Y(t))` over [0, 1], `sign` +1 on the outer boundary (traversed with the
    material on its left), -1 on a hole."""

    def __init__(self, el):
        if el.kind == 'seg':
            X, Y = [el.p[0], el.e[0]], [el.p[1], el.e[1]]
        elif el.kind == 'spline':
            X, Y = list(el.X), list(el.Y)
        else:
            raise AssertionError("an arc: S9f.3a's profiles hold lines and splines")
        self.spline = el.kind == 'spline'
        self.ctrl = el.ctrl if self.spline else (el.p, el.q)
        self.ctrlm = [Mv(p) for p in self.ctrl]
        self.sign = 1 if el.outer else -1
        self.X, self.Y = ptrim([F(c) for c in X]), ptrim([F(c) for c in Y])
        self.dX, self.dY = pder(self.X), pder(self.Y)
        self.deg = max(pdeg(self.X), pdeg(self.Y), 1)
        self.Xm, self.Ym, self.dXm, self.dYm = pm(self.X), pm(self.Y), pm(self.dX), pm(self.dY)
        cr = psub(pmul(self.X, self.dY), pmul(self.Y, self.dX))
        self.G = [pm(pc(pint(cr), F(1, 2))), pm(pc(pint(pmul(self.X, cr)), F(1, 3))),
                  pm(pc(pint(pmul(self.Y, cr)), F(1, 3)))]

    def point(self, t):
        return peval(self.Xm, t), peval(self.Ym, t)

    def deriv(self, t):
        return peval(self.dXm, t), peval(self.dYm, t)

    def green(self, a, b):
        return [self.sign*(peval(g, b)-peval(g, a)) for g in self.G]

    def line_poly(self, base, d):
        """`(S(t) - base) x d` (mpf)."""
        n = self.deg+1
        X = self.Xm+[0]*(n-len(self.Xm))
        Y = self.Ym+[0]*(n-len(self.Ym))
        out = [X[k]*d[1]-Y[k]*d[0] for k in range(n)]
        out[0] -= base[0]*d[1]-base[1]*d[0]
        return out

    def misses_line(self, base, d):
        side = [(p[0]-base[0])*d[1]-(p[1]-base[1])*d[0] for p in self.ctrlm]
        return all(x > 0 for x in side) or all(x < 0 for x in side)

    def box(self):
        us = [p[0] for p in self.ctrlm]
        vs = [p[1] for p in self.ctrlm]
        return min(us), max(us), min(vs), max(vs)


class Spline3:
    """The prism: its exact frame, heights and profile elements."""

    def __init__(self, case):
        self.case = case
        o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(case.frame))
        self.o, self.x, self.y, self.n = o, x, y, n
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        self.profile = cref.Profile(case.boundaries, case.tolerance)
        self.elements = [Element(el) for el in self.profile.elements]
        self.det = det3(x, y, n)
        assert self.det > 0
        self.inv = tuple(scale(r, 1/self.det) for r in (cross(y, n), cross(n, x), cross(x, y)))
        self.invm = tuple(Mv(r) for r in self.inv)
        self.om, self.xm, self.ym, self.nm = Mv(o), Mv(x), Mv(y), Mv(n)
        self.m = cross(x, y)

    def world(self, u, v, w):
        return tuple(self.om[i]+u*self.xm[i]+v*self.ym[i]+w*self.nm[i] for i in range(3))

    def local(self, X):
        return apply(self.invm, sub(X, self.om))

    def contains(self, X):
        u, v, w = self.local(X)
        return M(self.lo) < w < M(self.hi) and self.profile.inside(u, v)

    def wall(self, el):
        """The wall's base `P(t) = o + x X(t) + y Y(t)` as exact polynomial
        vectors, and `P'(t)`."""
        P = tuple(padd(padd([self.o[i]], pc(el.X, self.x[i])), pc(el.Y, self.y[i])) for i in range(3))
        dP = tuple(pder(p) for p in P)
        return P, dP

    def closed(self):
        A, Mu, Mw = self.profile.moments()
        h = M(self.hi-self.lo)
        det = M(self.det)
        vol = A*h*det
        c = self.world(Mu/A, Mw/A, M(self.lo+self.hi)/2)
        mm = mnorm(self.m)
        area = 2*A*mm
        for el in self.elements:
            area += h*wall_length(self, el)
        return vol, scale(c, vol), area


def vpoly_dot(v, P):
    """`v . P(t)` for a constant vector and a polynomial vector."""
    out = [F(0)]
    for i in range(3):
        out = padd(out, pc(P[i], v[i]))
    return out


def vpoly_sq(P):
    out = [F(0)]
    for i in range(3):
        out = padd(out, pmul(P[i], P[i]))
    return out


def wall_quadratic(S, el, ball):
    """The sphere's function along the element's wall's rulings: `A w^2 + 2
    B(t) w + C(t)`: (A, B, C) exact."""
    P, _ = S.wall(el)
    Q = tuple(psub(P[i], [ball.c[i]]) for i in range(3))
    A = dot(S.n, S.n)
    B = vpoly_dot(S.n, Q)
    C = psub(vpoly_sq(Q), [ball.r2])
    return A, B, C


def wall_length(S, el):
    """`int_0^1 |S'(t) x n| dt` (world)."""
    def f(t):
        du, dv = el.deriv(t)
        T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
        return [mnorm(cross(T, S.nm))]
    if not el.spline:
        return f(mp.mpf(0))[0]
    est, _ = cref.integrate(f, mp.mpf(0), mp.mpf(1), mp.mpf(10)**-36)
    return est[0]


# ------------------------------------------------------------------ a slicing direction

class Way:
    """Slices `d . X = s` of the pair (`d . n > 0`)."""

    def __init__(self, pair, d):
        S = pair.S
        self.pair, self.d = pair, tuple(F(c) for c in d)
        dn = dot(self.d, S.n)
        assert dn > 0, 'a slicing direction along the caps'
        self.dn = dn
        self.wu, self.wv = -dot(self.d, S.x)/dn, -dot(self.d, S.y)/dn
        self.do = dot(self.d, S.o)
        self.constant = self.wu == 0 and self.wv == 0
        self.J = M(S.det/dn)
        self.dm = Mv(self.d)
        self.dd = M(dot(self.d, self.d))
        self.dlen = mp.sqrt(self.dd)
        # An orthonormal basis of the slices' plane, oriented so its image in
        # (u, v) runs counter-clockwise.
        dh = scale(self.dm, 1/self.dlen)
        seed = min(((1, 0, 0), (0, 1, 0), (0, 0, 1)), key=lambda e: abs(dot(e, dh)))
        e1 = sub(Mv(seed), scale(dh, dot(Mv(seed), dh)))
        e1 = scale(e1, 1/mp.sqrt(dot(e1, e1)))
        e2 = cross(dh, e1)
        L1, L2 = apply(S.invm, e1)[:2], apply(S.invm, e2)[:2]
        if L1[0]*L2[1]-L1[1]*L2[0] < 0:
            e2, L2 = scale(e2, -1), (-L2[0], -L2[1])
        self.e, self.L = (e1, e2), (L1, L2)
        # World points of (u, v) in a slice: X0(s) + u Xu + v Xv.
        self.Xu = tuple(S.xm[i]+M(self.wu)*S.nm[i] for i in range(3))
        self.Xv = tuple(S.ym[i]+M(self.wv)*S.nm[i] for i in range(3))
        self.area_factor = mnorm(cross(self.Xu, self.Xv))

    def w0(self, s):
        return (s-M(self.do))/M(self.dn)

    def s_of(self, X):
        return dot(self.dm, X)

    def X0(self, s):
        S = self.pair.S
        w0 = self.w0(s)
        return tuple(S.om[i]+w0*S.nm[i] for i in range(3))

    def world_moments(self, s, g):
        """A region's (A, Mu, Mv) in (u, v) as its world volume element's
        integrand: (A, the moments' three components)."""
        A, Mu, Mw = g
        X0 = self.X0(s)
        return [A]+[X0[i]*A+self.Xu[i]*Mu+self.Xv[i]*Mw for i in range(3)]

    def events(self):
        """[(s, tag, point)]: the slices of the arrangement's vertices and of
        its edges' extremes (see the module's docstring)."""
        return events(self.pair, self)


# ------------------------------------------------------------------ one slice's regions

class Piece:
    __slots__ = ('family', 'green', 'inside', 'points', 'curve', 'angle', 'flags')

    def __init__(self, family, green, inside, points, curve, angle=None, flags=None):
        self.family, self.green, self.inside, self.points = family, green, inside, points
        self.curve, self.angle, self.flags = curve, angle, flags


def seg_green(p, q):
    (x0, y0), (x1, y1) = p, q
    cr = x0*y1-x1*y0
    return [cr/2, (x0+x1)*cr/6, (y0+y1)*cr/6]


def ellipse_green(E0, e1, e2, a, b):
    """`(u dv - v du) / 2`, `u (u dv - v du) / 3`, `v (u dv - v du) / 3` over
    `E0 + cos(t) e1 + sin(t) e2` from `a` to `b`."""
    a0, a1, a2 = E0[0], e1[0], e2[0]
    b0, b1, b2 = E0[1], e1[1], e2[1]
    De = a1*b2-a2*b1
    al = a0*b2-b0*a2
    be = b0*a1-a0*b1
    Ic = mp.sin(b)-mp.sin(a)
    Is = -(mp.cos(b)-mp.cos(a))
    Icc = (b-a)/2+(mp.sin(2*b)-mp.sin(2*a))/4
    Iss = (b-a)/2-(mp.sin(2*b)-mp.sin(2*a))/4
    Ics = (mp.sin(b)**2-mp.sin(a)**2)/2
    L = b-a
    area = (De*L+al*Ic+be*Is)/2

    def mom(k0, k1, k2):
        return (k0*De*L+(k0*al+k1*De)*Ic+(k0*be+k2*De)*Is+k1*al*Icc+(k1*be+k2*al)*Ics+k2*be*Iss)/3
    return [area, mom(a0, a1, a2), mom(b0, b1, b2)]


def solve2(a, b, c, d, e, f):
    """`a x + b y = e`, `c x + d y = f`."""
    det = a*d-b*c
    return (e*d-b*f)/det, (a*f-e*c)/det


class Slice:
    """The pair's sections in the slice `d . X = s` as regions of `(u, v)`:
    the prism's `P` (the profile cut by the strip, unless `strip` is off)
    and the sphere's `B` (the projected circle, cut by the hemisphere's line
    unless `hemi` is off), their boundaries' pieces classified."""

    def __init__(self, way, s, strip=True, hemi=True, samples=0):
        self.way, self.s = way, s
        pair = way.pair
        S, ball = pair.S, pair.ball
        self.S, self.ball = S, ball
        w0 = way.w0(s)
        lo, hi = M(S.lo), M(S.hi)
        self.wcoef = (w0, M(way.wu), M(way.wv))
        # The strip.
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
        # The ellipse.
        self.ellipse = None
        dc = way.s_of(ball.cm)
        rho2 = ball.rm**2-(s-dc)**2/way.dd
        if rho2 > 0:
            rho = mp.sqrt(rho2)
            cp = add(ball.cm, scale(way.dm, (s-dc)/way.dd))
            E0 = S.local(cp)[:2]
            e1 = (rho*way.L[0][0], rho*way.L[0][1])
            e2 = (rho*way.L[1][0], rho*way.L[1][1])
            det = e1[0]*e2[1]-e1[1]*e2[0]
            K = ((e2[1]/det, -e2[0]/det), (-e1[1]/det, e1[0]/det))
            ru = mp.sqrt(e1[0]**2+e2[0]**2)
            rv = mp.sqrt(e1[1]**2+e2[1]**2)
            self.ellipse = (E0, e1, e2, K, (E0[0]-ru, E0[0]+ru, E0[1]-rv, E0[1]+rv))
        # The hemisphere's halfplane h0 + hu u + hv v > 0.
        self.hemi_all = True
        self.hemi_line = None
        if ball.half and hemi:
            ns = Mv(ball.ns)
            X0 = way.X0(s)
            h0 = dot(ns, sub(X0, ball.cm))
            hu, hv = dot(ns, way.Xu), dot(ns, way.Xv)
            self.hcoef = (h0, hu, hv)
            if dot(ball.ns, S.x)+way.wu*dot(ball.ns, S.n) == 0 and dot(ball.ns, S.y)+way.wv*dot(ball.ns, S.n) == 0:
                self.hemi_all = h0 > 0
            else:
                self.hemi_all = None
                g2 = hu*hu+hv*hv
                self.hemi_line = ((-h0*hu/g2, -h0*hv/g2), (hv, -hu))
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

    def in_hemi(self, p):
        if self.hemi_all is not None:
            return self.hemi_all
        h0, hu, hv = self.hcoef
        return h0+hu*p[0]+hv*p[1] > 0

    def in_B(self, p):
        return self.in_ellipse(p) and self.in_hemi(p)

    # ---- crossings

    def ellipse_angle(self, p):
        E0, _, _, K, _ = self.ellipse
        a = (p[0]-E0[0], p[1]-E0[1])
        return mp.atan2(K[1][0]*a[0]+K[1][1]*a[1], K[0][0]*a[0]+K[0][1]*a[1])

    def element_ellipse(self, el):
        """[t]: the element's crossings of the ellipse."""
        E0, _, _, K, (u0, u1, v0, v1) = self.ellipse
        bu0, bu1, bv0, bv1 = el.box()
        if bu1 < u0 or bu0 > u1 or bv1 < v0 or bv0 > v1:
            return []
        n = el.deg+1
        X = el.Xm+[0]*(n-len(el.Xm))
        Y = el.Ym+[0]*(n-len(el.Ym))
        X = [X[0]-E0[0]]+X[1:]
        Y = [Y[0]-E0[1]]+Y[1:]
        a = [K[0][0]*x+K[0][1]*y for x, y in zip(X, Y)]
        b = [K[1][0]*x+K[1][1]*y for x, y in zip(X, Y)]
        q = [mp.mpf(0)]*(2*n-1)
        for i in range(n):
            for j in range(n):
                q[i+j] += a[i]*a[j]+b[i]*b[j]
        q[0] -= 1
        return roots01(q)

    def element_line(self, el, line):
        """[(t, lambda)]."""
        base, d = line
        if el.misses_line(base, d):
            return []
        out = []
        dd = d[0]*d[0]+d[1]*d[1]
        for t in roots01(el.line_poly(base, d)):
            p = el.point(t)
            out.append((t, ((p[0]-base[0])*d[0]+(p[1]-base[1])*d[1])/dd))
        return out

    def line_ellipse(self, line):
        """[lambda]."""
        base, d = line
        E0, _, _, K, _ = self.ellipse
        a0 = (base[0]-E0[0], base[1]-E0[1])
        p0 = (K[0][0]*a0[0]+K[0][1]*a0[1], K[1][0]*a0[0]+K[1][1]*a0[1])
        p1 = (K[0][0]*d[0]+K[0][1]*d[1], K[1][0]*d[0]+K[1][1]*d[1])
        A = p1[0]**2+p1[1]**2
        B = p0[0]*p1[0]+p0[1]*p1[1]
        C = p0[0]**2+p0[1]**2-1
        disc = B*B-A*C
        if disc <= 0:
            return []
        sq = mp.sqrt(disc)
        return [(-B-sq)/A, (-B+sq)/A]

    @staticmethod
    def line_line(l1, l2):
        (b1, d1), (b2, d2) = l1, l2
        den = d1[0]*d2[1]-d1[1]*d2[0]
        if den == 0:
            return None
        w = (b2[0]-b1[0], b2[1]-b1[1])
        return (w[0]*d2[1]-w[1]*d2[0])/den, (w[0]*d1[1]-w[1]*d1[0])/den

    # ---- the pieces

    def build(self):
        els = self.S.elements
        P_present = self.strip_all is not False
        ell = self.ellipse
        lines = list(self.strip_lines) if P_present else []
        hemi = self.hemi_line if ell is not None else None
        # Crossing parameters per curve.
        el_params = [[] for _ in els]
        el_toggles = []  # along the ellipse: angles of element crossings
        ell_params = []
        line_params = [[] for _ in lines]
        line_toggles = [[] for _ in lines]
        hemi_params, hemi_toggles = [], []
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
                if hemi is not None:
                    for t, lam in self.element_line(el, hemi):
                        el_params[i].append(t)
                        hemi_params.append(lam)
                        hemi_toggles.append(lam)
        if ell is not None:
            for k, line in enumerate(lines):
                for lam in self.line_ellipse(line):
                    line_params[k].append(lam)
                    base, d = line
                    ell_params.append(self.ellipse_angle((base[0]+lam*d[0], base[1]+lam*d[1])))
            if hemi is not None:
                for lam in self.line_ellipse(hemi):
                    hemi_params.append(lam)
                    base, d = hemi
                    ell_params.append(self.ellipse_angle((base[0]+lam*d[0], base[1]+lam*d[1])))
                for k, line in enumerate(lines):
                    x = self.line_line(line, hemi)
                    if x is not None:
                        line_params[k].append(x[0])
                        hemi_params.append(x[1])
        n = self.samples
        # The profile's elements (dP where inside the strip).
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
            # Strip lines (dP where inside the profile).
            for k, line in enumerate(lines):
                self.line_pieces(line, line_params[k], line_toggles[k], 'P', ('strip', k))
        # The ellipse (dB where inside the hemisphere's halfplane).
        if ell is not None:
            E0, e1, e2, _, _ = ell
            angles = sorted(set(th % (2*mp.pi) for th in ell_params))
            toggles = sorted(th % (2*mp.pi) for th in el_toggles)
            if not angles:
                spans = [(mp.mpf(0), 2*mp.pi)]
            else:
                spans = [(a, b) for a, b in zip(angles, angles[1:])]+[(angles[-1], angles[0]+2*mp.pi)]
            first = None
            for a, b in spans:
                m = (a+b)/2
                p = (E0[0]+mp.cos(m)*e1[0]+mp.sin(m)*e2[0], E0[1]+mp.cos(m)*e1[1]+mp.sin(m)*e2[1])
                if first is None:
                    first = (a, self.S.profile.inside(p[0], p[1]) if P_present else False)
                    inprof = first[1]
                else:
                    # Element crossings passed since the first piece's start.
                    passed = sum(1 for t in toggles if first[0] < t <= a)
                    inprof = first[1] != (passed % 2 == 1)
                inP = inprof and self.in_strip(p)
                inH = self.in_hemi(p)
                pts = None
                if n:
                    pts = [(E0[0]+mp.cos(a+(b-a)*k/n)*e1[0]+mp.sin(a+(b-a)*k/n)*e2[0],
                            E0[1]+mp.cos(a+(b-a)*k/n)*e1[1]+mp.sin(a+(b-a)*k/n)*e2[1]) for k in range(n+1)]
                self.pieces.append(Piece('B' if inH else 'E', ellipse_green(E0, e1, e2, a, b), inP, pts,
                                         ('ellipse',), angle=b-a, flags=(inH, inP)))
            if hemi is not None:
                self.line_pieces(hemi, hemi_params, hemi_toggles, 'B', ('hemi',))

    def line_pieces(self, line, params, toggles, family, curve):
        """A line's bounded pieces: on dP where inside the profile (a strip
        line), on dB where inside the ellipse (the hemisphere's)."""
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
                if not self.in_ellipse(mid):
                    continue
                inP = inprof and self.in_strip(mid) and self.strip_all is not False
                self.pieces.append(Piece('B', seg_green(pa, pb), inP, pts, curve))

    # ---- totals

    def regions(self):
        """{'P', 'B', 'C'}: each region's (A, Mu, Mv) in (u, v)."""
        out = {k: [mp.mpf(0)]*3 for k in 'PBC'}
        for pc_ in self.pieces:
            if pc_.family == 'E':
                continue
            out[pc_.family] = [x+y for x, y in zip(out[pc_.family], pc_.green)]
            if pc_.inside:
                out['C'] = [x+y for x, y in zip(out['C'], pc_.green)]
        return out

    def angles(self):
        """The section's circle's angles inside the hemisphere's halfplane
        and inside, or outside, the prism's section."""
        a_in = a_out = mp.mpf(0)
        for pc_ in self.pieces:
            if pc_.angle is None:
                continue
            inH, inP = pc_.flags
            if inH:
                if inP:
                    a_in += pc_.angle
                else:
                    a_out += pc_.angle
        return a_in, a_out

    def loops(self, op, obj_prism):
        """The result's boundary as loops of binary64 points."""
        kept = []
        for pc_ in self.pieces:
            if pc_.family == 'E':
                continue
            pts = [(float(p[0]), float(p[1])) for p in pc_.points]
            fam, ins = pc_.family, pc_.inside
            if op == 'common':
                keep, rev = ins, False
            elif op == 'fuse':
                keep, rev = not ins, False
            else:
                first = 'P' if obj_prism else 'B'
                keep, rev = (not ins, False) if fam == first else (ins, True)
            if keep:
                kept.append(pts[::-1] if rev else pts)
        return chain(kept)


def chain(pieces):
    """Closed loops from oriented polylines joined end to start."""
    loops = []
    pool = [p for p in pieces if len(p) >= 2]
    while pool:
        cur = pool.pop()
        loop = list(cur)
        for _ in range(len(pieces)+2):
            end = loop[-1]
            if math.dist(end, loop[0]) < 1e-9 and len(loop) > 2:
                break
            best, bi = None, None
            for i, q in enumerate(pool):
                dist = math.dist(end, q[0])
                if best is None or dist < best:
                    best, bi = dist, i
            if bi is None or best > 1e-6:
                break
            loop += pool.pop(bi)[1:]
        loops.append(loop)
    return loops


# ------------------------------------------------------------------ events

def events(pair, way):
    """The breakpoints of `way`'s slicing (see the module's docstring):
    [(s, tag, point)]."""
    S, ball = pair.S, pair.ball
    d, dm = way.d, way.dm
    out = []

    def add_point(tag, X):
        out.append((dot(dm, X), tag, X))
    los = (S.lo, S.hi)
    A, c, r2 = dot(S.n, S.n), ball.c, ball.r2
    # E1: the prism's vertices.
    for v in S.profile.vertices:
        for w in los:
            add_point('vertex', S.world(M(v[0]), M(v[1]), M(w)))
    # E2: cap edges' extremes.
    for el in S.elements:
        dp = padd(pc(el.dX, dot(d, S.x)), pc(el.dY, dot(d, S.y)))
        for t in exact_roots(dp):
            u, v = el.point(t)
            for w in los:
                add_point('cap_extreme', S.world(u, v, M(w)))
    # E3: the sphere's poles along d.
    for sg in (-1, 1):
        add_point('pole', add(ball.cm, scale(dm, sg*ball.rm/way.dlen)))
    # E4: the caps' circles' extremes.
    mm = dot(S.m, S.m)
    for w in los:
        off = dot(S.m, sub(add(S.o, scale(S.n, w)), c))
        rho2 = r2-off*off/mm
        if rho2 > 0:
            cw = add(c, scale(S.m, off/mm))
            mh = Mv(S.m)
            perp = sub(dm, scale(mh, dot(dm, mh)/M(mm)))
            pl = mnorm(perp)
            if pl > 0:
                for sg in (-1, 1):
                    add_point('cap_circle', add(Mv(cw), scale(perp, sg*mp.sqrt(M(rho2))/pl)))
    # E5: vertical edges against the sphere.
    for v in S.profile.vertices:
        P0 = add(S.o, add(scale(S.x, v[0]), scale(S.y, v[1])))
        Q = sub(P0, c)
        B, C = dot(S.n, Q), dot(Q, Q)-r2
        disc = B*B-A*C
        if disc > 0:
            for sg in (-1, 1):
                w = (-M(B)+sg*mp.sqrt(M(disc)))/M(A)
                add_point('edge_sphere', add(Mv(P0), scale(S.nm, w)))
    for el in S.elements:
        P, dP = S.wall(el)
        Aq, B, C = wall_quadratic(S, el, ball)
        # E6: cap edges against the sphere.
        for w in los:
            for t in exact_roots(padd(padd([Aq*w*w], pc(B, 2*w)), C)):
                u, v = el.point(t)
                add_point('cap_sphere', S.world(u, v, M(w)))
        # E7: the meeting's extremes.
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
                # `D = 0` there (the slices' trace of the ruling's plane): the
                # meeting's points on the ruling.
                b, cc = peval(pm(B), t), peval(pm(C), t)
                disc = b*b-M(Aq)*cc
                if disc >= 0:
                    for sg in (-1, 1):
                        add_point('meeting_extreme', S.world(u, v, (-b+sg*mp.sqrt(disc))/M(Aq)))
    if ball.half:
        ns = ball.ns
        nsm = Mv(ns)
        # E8: the rim's extremes.
        perp = sub(dm, scale(nsm, dot(dm, nsm)))
        pl = mnorm(perp)
        if pl > 0:
            for sg in (-1, 1):
                add_point('rim', add(ball.cm, scale(perp, sg*ball.rm/pl)))
        else:
            add_point('disc', ball.cm)
        nsn = dot(ns, S.n)
        # E9: the hemisphere's plane against vertical edges.
        if nsn != 0:
            for v in S.profile.vertices:
                P0 = add(S.o, add(scale(S.x, v[0]), scale(S.y, v[1])))
                w = dot(ns, sub(c, P0))/nsn
                add_point('edge_disc', add(Mv(P0), scale(S.nm, M(w))))
        # E11: the rim against the caps' planes.
        for w in los:
            for X in plane_circle(S.m, add(S.o, scale(S.n, w)), ball):
                add_point('rim_cap', X)
        for el in S.elements:
            P, dP = S.wall(el)
            Aq, B, C = wall_quadratic(S, el, ball)
            hp = vpoly_dot(ns, tuple(psub(P[i], [c[i]]) for i in range(3)))
            # E10: the plane against cap edges.
            for w in los:
                for t in exact_roots(padd(hp, [nsn*w])):
                    u, v = el.point(t)
                    add_point('cap_disc', S.world(u, v, M(w)))
            if nsn != 0:
                # E12: the crease w = -hp / nsn in F; E13: its extremes.
                wh = pc(hp, -1/nsn)
                for t in exact_roots(padd(padd(pc(pmul(wh, wh), Aq), pc(pmul(B, wh), 2)), C)):
                    u, v = el.point(t)
                    add_point('rim_wall', S.world(u, v, peval(pm(wh), t)))
                for t in exact_roots(padd(vpoly_dot(d, dP), pc(pder(wh), dot(d, S.n)))):
                    u, v = el.point(t)
                    add_point('crease_extreme', S.world(u, v, peval(pm(wh), t)))
            else:
                # E12 on a plane holding the axis: its generatrices.
                for t in exact_roots(hp):
                    b, cc = peval(pm(B), t), peval(pm(C), t)
                    disc = b*b-M(Aq)*cc
                    if disc > 0:
                        for sg in (-1, 1):
                            u, v = el.point(t)
                            add_point('rim_wall', S.world(u, v, (-b+sg*mp.sqrt(disc))/M(Aq)))
    return out


def plane_circle(m, p0, ball):
    """The hemisphere's rim against the plane through `p0` normal to `m`:
    its points (mpf)."""
    ns, c = ball.ns, ball.c
    L = cross(ns, m)
    if dot(L, L) == 0:
        return []
    # A point on both planes: ns . (X - c) = 0, m . (X - p0) = 0.
    a = dot(m, sub(p0, c))
    # X = c + alpha ns + beta m with ns . (alpha ns + beta m) = 0, m . (...) = a.
    nn, nm_, mm = dot(ns, ns), dot(ns, m), dot(m, m)
    det = nn*mm-nm_*nm_
    alpha, beta = (-nm_*a)/det, (nn*a)/det
    X0 = add(c, add(scale(ns, alpha), scale(m, beta)))
    Q = sub(X0, c)
    LL, QL, QQ = dot(L, L), dot(Q, L), dot(Q, Q)-ball.r2
    disc = QL*QL-LL*QQ
    if disc <= 0:
        return []
    return [add(Mv(X0), scale(Mv(L), (-M(QL)+sg*mp.sqrt(M(disc)))/M(LL))) for sg in (-1, 1)]


# ------------------------------------------------------------------ the pair

def case_size(S, ball):
    us = [v[0] for v in S.profile.vertices]
    vs = [v[1] for v in S.profile.vertices]
    for el in S.elements:
        us += [p[0] for p in el.ctrl]
        vs += [p[1] for p in el.ctrl]
    pts = [S.world(M(u), M(v), M(w)) for u in (min(us), max(us)) for v in (min(vs), max(vs))
           for w in (S.lo, S.hi)]
    lo = [min(min(p[i] for p in pts), ball.cm[i]-ball.rm) for i in range(3)]
    hi = [max(max(p[i] for p in pts), ball.cm[i]+ball.rm) for i in range(3)]
    return max(h-l for l, h in zip(lo, hi)), lo, hi


class Pair:
    def __init__(self, obj, tool):
        self.obj_prism = obj.sphere is None
        prism, sphere = (obj, tool) if self.obj_prism else (tool, obj)
        assert prism.sphere is None and sphere.sphere is not None
        self.S, self.ball = Spline3(prism), Ball(sphere)
        self.size, self.lo, self.hi = case_size(self.S, self.ball)
        self.first = Way(self, scale(self.S.m, 1/dot(self.S.m, self.S.n)))
        self.second = Way(self, SECOND if dot(SECOND, self.S.n) > 0 else scale(SECOND, -1))
        self._sliced = {}
        self._faces = None
        self._solids = None
        self.stats = {'breaks': 0}

    def breakpoints(self, way):
        ev = way.events()
        ss = [s for s, _, _ in ev]
        lo, hi = min(ss), max(ss)
        return cref.merge_breaks(ss, lo, hi, mp.mpf(10)**-30*self.size), ev

    def sliced(self, way):
        """The integrals over `way`'s slices: per region (P, B, C) its volume
        and world moments, and the sphere's face's angles."""
        key = id(way)
        if key in self._sliced:
            return self._sliced[key]
        breaks, _ = self.breakpoints(way)
        self.stats['breaks'] = max(self.stats['breaks'], len(breaks))
        archimedes = way is self.first

        def f(s):
            sl = Slice(way, s)
            reg = sl.regions()
            row = []
            for k in 'PBC':
                row += way.world_moments(s, reg[k])
            if archimedes:
                row += list(sl.angles())
            return row
        n = 14 if archimedes else 12
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
        if archimedes:
            out['angles'] = (total[12], total[13])
        out['error'] = err
        out['breaks'] = breaks
        self._sliced[key] = out
        return out

    def volumes(self, way=None):
        """{op: (volume, moments)} by one slicing."""
        sl = self.sliced(way or self.first)
        P, B, C = sl['P'], sl['B'], sl['C']
        obj, tool = (P, B) if self.obj_prism else (B, P)
        out = {'common': C,
               'cut': (obj[0]-C[0], tuple(x-y for x, y in zip(obj[1], C[1]))),
               'fuse': (P[0]+B[0]-C[0], tuple(x+y-z for x, y, z in zip(P[1], B[1], C[1])))}
        return out

    # ---- faces

    def faces(self):
        """[(input tag 'S' or 'B', face name, {'in': area, 'out': area},
        closed area)]."""
        if self._faces is not None:
            return self._faces
        S, ball = self.S, self.ball
        out = []
        mm = mnorm(S.m)
        prof_area = S.profile.moments()[0]
        # The caps.
        for name, w in (('cap_lo', S.lo), ('cap_hi', S.hi)):
            sl = Slice(self.first, M(dot(self.first.d, S.o)+w), strip=False)
            reg = sl.regions()
            inside = reg['C'][0]*mm
            out.append(('S', name, {'in': inside, 'out': reg['P'][0]*mm-inside}, prof_area*mm))
        # The walls.
        for i, el in enumerate(S.elements):
            cls = wall_classes(self, el)
            out.append(('S', f'wall{i}', cls, M(S.hi-S.lo)*wall_length(S, el)))
        # The sphere's face by Archimedes.
        a_in, a_out = self.sliced(self.first)['angles']
        k = ball.rm/self.first.dlen
        closed = 4*mp.pi*ball.rm**2 if not ball.half else 2*mp.pi*ball.rm**2
        out.append(('B', 'sphere', {'in': k*a_in, 'out': k*a_out}, closed))
        if ball.half:
            out.append(('B', 'disc', disc_classes(self), mp.pi*ball.rm**2))
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

    # ---- solids

    def solids(self):
        if self._solids is None:
            self._solids = count_solids(self)
        return self._solids

    def result(self, op):
        vol, mom = self.volumes()[op]
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        return self.solids()[op], vol, self.area(op), tuple(m/vol for m in mom)


def wall_classes(pair, el):
    """A wall's area inside and outside the sphere."""
    S, ball = pair.S, pair.ball
    A, B, C = wall_quadratic(S, el, ball)
    lo, hi = S.lo, S.hi
    P, _ = S.wall(el)
    polys = [psub(pmul(B, B), pc(C, A))]
    for w in (lo, hi):
        polys.append(padd(padd([A*w*w], pc(B, 2*w)), C))
    hp = None
    nsn = F(0)
    if ball.half:
        ns = ball.ns
        nsn = dot(ns, S.n)
        hp = vpoly_dot(ns, tuple(psub(P[i], [ball.c[i]]) for i in range(3)))
        if nsn != 0:
            wh = pc(hp, -1/nsn)
            polys.append(padd(padd(pc(pmul(wh, wh), A), pc(pmul(B, wh), 2)), C))
            for w in (lo, hi):
                polys.append(psub(wh, [w]))
        else:
            polys.append(hp)
    ts = [mp.mpf(0), mp.mpf(1)]
    for p in polys:
        ts += exact_roots(p, slack=mp.mpf(0))
    breaks = cref.merge_breaks(ts, mp.mpf(0), mp.mpf(1), mp.mpf(10)**-30)
    Am, Bm, Cm = M(A), pm(B), pm(C)
    hpm = pm(hp) if hp is not None else None
    lom, him, h = M(lo), M(hi), M(hi-lo)

    def f(t):
        du, dv = el.deriv(t)
        T = tuple(du*S.xm[i]+dv*S.ym[i] for i in range(3))
        g = mnorm(cross(T, S.nm))
        b, cc = peval(Bm, t), peval(Cm, t)
        disc = b*b-Am*cc
        inside = mp.mpf(0)
        if disc > 0:
            sq = mp.sqrt(disc)
            a, z = (-b-sq)/Am, (-b+sq)/Am
            a, z = max(a, lom), min(z, him)
            if hpm is not None:
                hv = peval(hpm, t)
                if nsn > 0:
                    a = max(a, -hv/M(nsn))
                elif nsn < 0:
                    z = min(z, -hv/M(nsn))
                elif hv <= 0:
                    z = a
            inside = max(z-a, mp.mpf(0))
        return [inside*g, (h-inside)*g]
    total = [mp.mpf(0), mp.mpf(0)]
    for a, b in zip(breaks, breaks[1:]):
        est, _ = cref.integrate(f, a, b, mp.mpf(10)**-34*pair.size**2)
        if est is not None:
            total = [x+y for x, y in zip(total, est)]
    return {'in': total[0], 'out': total[1]}


def disc_classes(pair):
    """A hemisphere's disc inside and outside the prism."""
    S, ball = pair.S, pair.ball
    ns = ball.ns
    if dot(ns, S.n) != 0:
        way = Way(pair, ns if dot(ns, S.n) > 0 else scale(ns, -1))
        sl = Slice(way, way.s_of(ball.cm), hemi=False)
        reg = sl.regions()
        k = way.area_factor
        return {'in': reg['C'][0]*k, 'out': (reg['B'][0]-reg['C'][0])*k}
    # The plane holds the prism's axis: the profile's chords along its trace
    # times the heights, against the disc in the plane's coordinates `X = c +
    # lam e + mu n` (`e = ns x n` along the trace; the prism's height there
    # `wc + lam ez + mu`).
    e = cross(ns, S.n)
    lc = apply(S.inv, sub(ball.c, S.o))
    le = apply(S.inv, e)
    chords = S.profile.chords(Mv(lc[:2]), Mv(le[:2]))
    lo, hi = M(S.lo), M(S.hi)
    wc, ez = M(lc[2]), M(le[2])
    ee, en, nn = M(dot(e, e)), M(dot(e, S.n)), M(dot(S.n, S.n))
    r2 = ball.rm**2
    inside = mp.mpf(0)

    def width(lam):
        # The disc's mu at lam: ee lam^2 + 2 en lam mu + nn mu^2 <= r^2.
        disc = (en*lam)**2-nn*(ee*lam*lam-r2)
        if disc <= 0:
            return [mp.mpf(0)]
        sq = mp.sqrt(disc)
        m0, m1 = (-en*lam-sq)/nn, (-en*lam+sq)/nn
        a, b = max(m0, lo-wc-lam*ez), min(m1, hi-wc-lam*ez)
        return [max(b-a, mp.mpf(0))]
    lim = ball.rm*mp.sqrt(nn/(ee*nn-en*en))
    for (t0, _), (t1, _) in chords:
        pts = [t0, t1, -lim, lim]
        # Where the disc's edge crosses a cap: mu = k - wc - lam ez.
        for k in (lo, hi):
            m0 = k-wc
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
    k = mnorm(cross(e, S.n))
    inside *= k
    return {'in': inside, 'out': mp.pi*ball.rm**2-inside}


# ------------------------------------------------------------------ solids

def interior_points(loops):
    """Points strictly inside a component (its outer loop first, holes
    after), by horizontal scanlines."""
    out = []
    vs = [p[1] for p in loops[0]]
    lo, hi = min(vs), max(vs)
    for k in range(1, 8):
        y = lo+(hi-lo)*k/8
        xs = []
        for loop in loops:
            for (x0, y0), (x1, y1) in zip(loop, loop[1:]+loop[:1]):
                if (y0 <= y < y1) or (y1 <= y < y0):
                    xs.append(x0+(y-y0)*(x1-x0)/(y1-y0))
        xs.sort()
        for a, b in zip(xs[0::2], xs[1::2]):
            if b > a:
                out.append(((a+b)/2, y))
    return out


def in_polygon(p, loop):
    x, y = p
    inside = False
    for (x0, y0), (x1, y1) in zip(loop, loop[1:]+loop[:1]):
        if (y0 <= y < y1) or (y1 <= y < y0):
            if x < x0+(y-y0)*(x1-x0)/(y1-y0):
                inside = not inside
    return inside


def signed_area(loop):
    return sum(x0*y1-x1*y0 for (x0, y0), (x1, y1) in zip(loop, loop[1:]+loop[:1]))/2


def components(loops):
    """[[outer, holes...]]: outer loops (counter-clockwise) with the holes
    each lies directly in."""
    outer = [lp for lp in loops if signed_area(lp) > 0]
    holes = [lp for lp in loops if signed_area(lp) < 0]
    comps = [[lp] for lp in outer]
    for h in holes:
        best = None
        for c in comps:
            if in_polygon(h[0], c[0]) and (best is None or signed_area(c[0]) < signed_area(best[0])):
                best = c
        if best is not None:
            best.append(h)
    return comps


def in_component(p, comp):
    return in_polygon(p, comp[0]) and not any(in_polygon(p, h) for h in comp[1:])


def count_solids(pair):
    way = pair.first
    breaks = pair.sliced(way)['breaks']
    counts = {}
    with mp.workdps(20):
        slices = []
        for a, b in zip(breaks, breaks[1:]):
            for f in (1e-9, 0.25, 0.5, 0.75, 1-1e-9):
                s = mp.mpf(a)+(mp.mpf(b)-mp.mpf(a))*f
                slices.append(Slice(way, s, samples=24))
        for op in OPS:
            nodes = []
            per = []
            for sl in slices:
                comps = components(sl.loops(op, pair.obj_prism))
                comps = [c for c in comps if abs(signed_area(c[0])) > 1e-14*float(pair.size)**2]
                ids = []
                for c in comps:
                    nodes.append(len(nodes))
                    ids.append((len(nodes)-1, c, interior_points(c)))
                per.append(ids)
            parent = list(range(len(nodes)))

            def find(i):
                while parent[i] != i:
                    parent[i] = parent[parent[i]]
                    i = parent[i]
                return i
            for prev, cur in zip(per, per[1:]):
                for i, ci, pi in prev:
                    for j, cj, pj in cur:
                        if any(in_component(p, cj) for p in pi) or any(in_component(p, ci) for p in pj):
                            parent[find(i)] = find(j)
            counts[op] = len({find(i) for i in range(len(nodes))})
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
