#!/usr/bin/env python3
"""Independent reference for S9d.2 of REVIEW_NOTES.md: Booleans of a sphere
(whole, a cap or a zone, `Solid::sphere_with`) against a prism whose profile
holds lines, arcs and circles (planar and cylindrical walls) in any relative
position, and of two spheres. A cylinder and a sphere meet in parallels (the
centre on the axis) or a quartic; two spheres in a circle.

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal): the sphere `|X - o|^2 <= R^2` about its stored origin with a
cap's or zone's end planes as the kernel's `Ball` reads them, through `o +
h n` normal to the stored axis `n` (`AxisSphere`: in an exact frame
`sphere_boolean_reference.Sphere`'s affine plane, `x * y = n` exactly; in a
turned frame, S9d.2c, a plane turned from it by the axes' rounding); the
prism the points `o + u x + v y + w n` with `(u, v)` in the profile (its
stored points, arcs' centres and radii, `identity_reference.stored`; an arc
the exact circle between its end points, holes reversed so the material is
on the left of every element) and `w` between the offsets. Nothing here uses
the kernel, a surface/surface intersection or an arrangement, and nothing is
shared with another reference but `stored`, `stored_axes`, S9d.1's number
and vector helpers, its sphere model and quadrature, and the exact
polynomial helpers of `boolean_reference.py`.

* **Slices.** Both solids are sliced by the planes `d . X = s` (`d` the
  prism's cap normal `x * y`, or the first sphere's axis for two spheres,
  or a cap's or zone's stored axis where its end planes are not normal to
  the prism's axis, S9d.2c's turned caps; a cap's or zone's end planes must
  be slices; any rational `d` is accepted and a second direction checks
  the first). A slice is charted affinely by
  rationals, `X = O0 + s tau + a e1 + b e2` (for the prism's own direction
  `e1 = x`, `e2 = y`, `tau = n / (d . n)`: the chart is the profile's own
  `(u, v)`), so in `q = (a, b)` every section is exact: a sphere's an
  ellipse `q^T G q + 2 q^T g(s) + k(s) <= 0` (`G` the chart's Gram matrix, a
  circle in an exact frame), the prism's the preimage of its profile under
  the affine map `q -> (u, v)` (segments to segments, a circle `C + r (cos
  t, sin t)` to `C' + A' (cos t, sin t)`, the same parameter) cut by the
  strip between the caps (lines, or all or nothing when the caps are
  slices). Every boundary is a segment or an arc `C + A (cos t, sin t)`;
  a sphere's arc is taken with `A = sqrt(rho^2) L^-1` (`G = L^T L`), so `t`
  is the section circle's true angle. Two regions are taken apart as in
  S9d.1: each boundary cut where it crosses the other's (a segment against
  an arc a quadratic, two arcs the quartic `z^2 f(z)` in `z = e^(it)`, its
  roots on the unit circle refined by Newton), each piece classified at its
  midpoint (an ellipse's inequality, half-planes, a profile's parity along a
  ray of irrational slope): common `dD n P + dP n D`, fuse `dD \\ P + dP \\
  D`, `D - P` `dD \\ P - dP n D`, `P - D` `dP \\ D - dD n P`, their areas
  and first moments by Green's theorem over segments and arcs in closed form
  (an arc's integrands trigonometric polynomials of degree three,
  integrated exactly). The sphere's face by Archimedes: `R / |d|` times the
  integral of its section's angle inside the other. Breakpoints, all roots
  of exact polynomials in `s`: the planes tangent to a sphere; through a
  prism vertex or a zone plane; an edge's meetings with a sphere
  (quadratics); a face's plane tangent to a sphere's section (S9d.1's
  quadratic); a cap circle's meetings with a sphere (the exact quartic in
  `tan(theta / 2)`); a cap's trace tangent to its circle; and two sections
  tangent (a profile circle's preimage and a sphere's, or two spheres'):
  the discriminant in `lambda` of the cubic `det(lambda Q1(s) + Q2(s))`
  (conics tangent where the pencil has a double root), or, where it
  vanishes identically (concentric sections of proportional Gram matrices),
  where the two conics coincide. Between breakpoints Gauss-Legendre after
  `s = a + (b - a)(1 - cos t)/2`, doubling the rule (12 to 96 nodes) until
  successive estimates agree within 1e-33 of the case's size to the fourth,
  else halving (S9d.1's `integrate`).
* **Other faces.** A prism cap (a slice) and a zone's end disc (a slice)
  against the other's section there, a flat wall (its parallelogram in its
  own `(lambda, w)`) against the sphere's ellipse there, in closed form (in
  by the classified boundaries, out the rest; a coplanar face of the same or
  the opposite orientation apart). A cylindrical wall by its angle: each
  generatrix `P(theta) + w n` inside the sphere between the roots of a
  quadratic in `w`, clipped to the caps and the zone, integrated with the
  wall's element `r |(-sin x + cos y) * n|` between the breakpoints where
  the generatrix touches the sphere or a cap circle or a zone plane's
  circle meets it (quartics in `tan(theta / 2)`). A prism cut obliquely by
  the slices (S9d.2c: a turned cap's end planes the slices) has its caps
  classified in their own `(u, v)` against the sphere's ellipse and the
  zone's half-planes there, as a flat wall is, and a cylindrical wall's
  zone bound along each generatrix varies with its angle: it adds the
  angles where the bound meets the sphere's roots (the zone's rim crossing
  the wall, a quartic) or the prism's ends (a quadratic).
* **Solids** (the regularized Boolean's maximal connected regions). In each
  interval between breakpoints the result's section has fixed topology: its
  boundary pieces chained into loops, loops into components (an outer loop
  and the holes inside it), each keyed by its loops' cyclic sequences of
  source faces. Components are followed through an interval by their keys
  (at both of its ends, `delta` inside) and joined across a breakpoint when
  most of 64 points spread over the smaller of their sections on both sides
  lie inside the other: the limits there of one solid's sections nest (the
  sections change continuously, or by the other input's section added or
  taken away at a cap), those of two solids meet in no area (regions
  touching along a curve or at a point stay apart). (Sections sharing a
  wall's line make boundary classification ambiguous there; points are
  not.)

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9d.1
references. A result keeps (object `A`, tool `B`): fuse `A`'s outside and
shared same-orientation pieces and `B`'s outside, cut `A`'s outside and
shared opposite pieces and `B`'s inside, common `A`'s inside and shared
same pieces and `B`'s inside.
"""
from fractions import Fraction as F
import math

import mpmath as mp

from identity_reference import Spline, stored
from curve_surface_reference import stored_axes
from boolean_reference import padd, pdeg, pmul, psub, ptrim, squarefree
import sphere_boolean_reference as s1
from sphere_boolean_reference import M, Mv, add, compose, cross, dot, integrate, is_zero, roots2, scale, sub

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
HALF_PI = math.pi/2

# ------------------------------------------------------------------ exact polynomials in s


def pk(p, k):
    return ptrim([c*k for c in p]) if k else [F(0)]


def pzero(p):
    return pdeg(p) < 0


def plin(v0, v1):
    """The vector polynomial `v0 + s v1` (one polynomial per coordinate)."""
    return tuple(ptrim([F(a), F(b)]) for a, b in zip(v0, v1))


def pdotv(u, v):
    out = [F(0)]
    for a, b in zip(u, v):
        out = padd(out, pmul(a, b))
    return out


def pdet3(m):
    """The determinant of a 3x3 matrix of polynomials."""
    out = [F(0)]
    for (i, j, k), sign in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                            ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
        out = padd(out, pk(pmul(pmul(m[0][i], m[1][j]), m[2][k]), sign))
    return out


def real_roots(p):
    """The real roots (mpf) of an exact polynomial: each square-free factor
    by `mp.polyroots` at 60 digits, an imaginary part below 1e-12 counting
    as real (a spurious breakpoint costs nothing)."""
    p = ptrim([F(c) for c in p])
    if pdeg(p) <= 0:
        return []
    out = []
    for f, _ in squarefree(p):
        if pdeg(f) == 1:
            out.append(-M(f[0])/M(f[1]))
            continue
        if pdeg(f) == 2:
            out += roots2(f)
            continue
        with mp.workdps(60):
            desc = [M(c) for c in reversed(f)]
            try:
                found = mp.polyroots(desc, maxsteps=400, extraprec=400)
            except mp.NoConvergence:
                found = mp.polyroots(desc, maxsteps=4000, extraprec=2000)
            for r in found:
                re, im = mp.re(r), mp.im(r)
                if abs(im) <= mp.mpf(10)**-12*(1+abs(re)):
                    out.append(+re)
    return out


class Conic:
    """`q^T G q + 2 q^T g(s) + k(s)`: `G` a constant symmetric 2x2 of
    Fractions, `g` two polynomials, `k` one."""

    def __init__(self, G, g, k):
        self.G, self.g, self.k = G, g, k

    def matrix(self):
        G, g, k = self.G, self.g, self.k
        c = lambda x: [F(x)]
        return [[c(G[0][0]), c(G[0][1]), g[0]], [c(G[1][0]), c(G[1][1]), g[1]], [g[0], g[1], k]]

    def centre_poly(self):
        """`-G^-1 g(s)` and `rho(s) = g^T G^-1 g - k` (the conic
        `(q - q0)^T G (q - q0) = rho`)."""
        G, g = self.G, self.g
        det = G[0][0]*G[1][1]-G[0][1]*G[1][0]
        inv = [[G[1][1]/det, -G[0][1]/det], [-G[1][0]/det, G[0][0]/det]]
        q0 = tuple(pk(padd(pk(g[0], inv[i][0]), pk(g[1], inv[i][1])), -1) for i in range(2))
        Ginv_g = tuple(padd(pk(g[0], inv[i][0]), pk(g[1], inv[i][1])) for i in range(2))
        rho = psub(pdotv(g, Ginv_g), self.k)
        return q0, rho


def tangency_poly(c1, c2):
    """A polynomial in `s` vanishing where the two conics are tangent: the
    discriminant in `lambda` of `det(lambda Q1 + Q2)`, or, when that
    vanishes identically, where concentric conics of proportional Gram
    matrices coincide."""
    Q1, Q2 = c1.matrix(), c2.matrix()
    vals = {}
    for lam in (0, 1, -1, 2):
        m = [[padd(pk(Q1[i][j], F(lam)), Q2[i][j]) for j in range(3)] for i in range(3)]
        vals[lam] = pdet3(m)
    # Cubic in lambda through the four values (Lagrange, exact).
    lams = (0, 1, -1, 2)
    coef = [[F(0)] for _ in range(4)]
    for i, li in enumerate(lams):
        # basis polynomial prod_{j != i} (lam - lj)/(li - lj), ascending in lam
        basis = [F(1)]
        den = F(1)
        for j, lj in enumerate(lams):
            if j != i:
                basis = pmul(basis, [F(-lj), F(1)])
                den *= li-lj
        for k in range(len(basis)):
            coef[k] = padd(coef[k], pk(vals[li], basis[k]/den))
    d0, c1_, b, a = coef  # a lam^3 + b lam^2 + c lam + d
    c_, d = c1_, d0
    disc = [F(0)]
    terms = [(pmul(pmul(b, b), pmul(c_, c_)), 1), (pmul(a, pmul(c_, pmul(c_, c_))), -4),
             (pmul(pmul(b, pmul(b, b)), d), -4), (pmul(pmul(a, a), pmul(d, d)), -27),
             (pmul(pmul(a, b), pmul(c_, d)), 18)]
    for t, k in terms:
        disc = padd(disc, pk(t, F(k)))
    if not pzero(disc):
        return disc
    # Concentric sections of proportional Gram matrices.
    G1, G2 = c1.G, c2.G
    ratio = G1[0][0]/G2[0][0]
    assert all(G1[i][j] == ratio*G2[i][j] for i in range(2) for j in range(2)), 'a degenerate pencil'
    q1, r1 = c1.centre_poly()
    q2, r2 = c2.centre_poly()
    assert all(pzero(psub(q1[i], q2[i])) for i in range(2)), 'a degenerate pencil off the centres'
    return psub(r1, pk(r2, ratio))


# ------------------------------------------------------------------ trigonometric polynomials

class Trig:
    """A trigonometric polynomial in `(cos t, sin t)` with Fraction
    coefficients: {(i, j): c} for `c cos^i sin^j`."""

    def __init__(self, terms=None):
        self.t = {k: v for k, v in (terms or {}).items() if v != 0}

    @staticmethod
    def linear(c0, cc, cs):
        return Trig({(0, 0): F(c0), (1, 0): F(cc), (0, 1): F(cs)})

    def __add__(self, o):
        out = dict(self.t)
        for k, v in o.t.items():
            out[k] = out.get(k, F(0))+v
        return Trig(out)

    def scale(self, k):
        return Trig({a: v*k for a, v in self.t.items()})

    def __mul__(self, o):
        out = {}
        for (i, j), v in self.t.items():
            for (k, l), w in o.t.items():
                key = (i+k, j+l)
                out[key] = out.get(key, F(0))+v*w
        return Trig(out)

    def deriv(self):
        out = {}
        for (i, j), v in self.t.items():
            if i:
                out[(i-1, j+1)] = out.get((i-1, j+1), F(0))-i*v
            if j:
                out[(i+1, j-1)] = out.get((i+1, j-1), F(0))+j*v
        return Trig(out)

    def degree(self):
        return max((i+j for i, j in self.t), default=0)

    def tan_half(self):
        """The polynomial in `t = tan(theta / 2)` times `(1 + t^2)^deg`."""
        n = self.degree()
        out = [F(0)]
        for (i, j), v in self.t.items():
            term = [v]
            for _ in range(i):
                term = pmul(term, [F(1), F(0), F(-1)])
            for _ in range(j):
                term = pmul(term, [F(0), F(2)])
            for _ in range(n-i-j):
                term = pmul(term, [F(1), F(0), F(1)])
            out = padd(out, term)
        return out

    def value(self, th):
        c, s = mp.cos(th), mp.sin(th)
        return sum((M(v)*c**i*s**j for (i, j), v in self.t.items()), mp.mpf(0))

    def is_zero(self):
        return not self.t

    def angles(self):
        """Its real roots in `(-pi, pi]`."""
        if self.is_zero():
            return []
        out = [2*mp.atan(t) for t in real_roots(self.tan_half())]
        # theta = pi (t infinite): the value at cos = -1, sin = 0.
        if sum((v*(-1)**i for (i, j), v in self.t.items() if j == 0), F(0)) == 0:
            out.append(+mp.pi)
        return out


def trig_vec(v0, v1, v2):
    """`v0 + cos v1 + sin v2` as three Trig coordinates."""
    return tuple(Trig.linear(a, b, c) for a, b, c in zip(v0, v1, v2))


def trig_dot(u, v):
    out = Trig()
    for a, b in zip(u, v):
        out = out+a*b
    return out


# ------------------------------------------------------------------ 2D pieces



def span_param(t0, t1, phi):
    """`phi` (an angle) mapped into the open traversal range from `t0` to
    `t1` (either way), or None."""
    T = 2*mp.pi
    if t1 > t0:
        x = t0+(phi-t0) % T
        return x if t0 < x < t1 else None
    x = t0-(t0-phi) % T
    return x if t1 < x < t0 else None


class Seg:
    kind = 'seg'

    def __init__(self, p, q, tag):
        self.p, self.q, self.tag = p, q, tag

    def start(self):
        return self.p

    def end(self):
        return self.q

    def point(self, t):
        return (self.p[0]+t*(self.q[0]-self.p[0]), self.p[1]+t*(self.q[1]-self.p[1]))

    def lo_hi(self):
        return mp.mpf(0), mp.mpf(1)

    def sub(self, a, b):
        return Seg(self.point(a), self.point(b), self.tag)

    def reversed(self):
        return Seg(self.q, self.p, self.tag)

    def green(self):
        return s1.green_seg(self.p, self.q)

    def angle(self):
        return mp.mpf(0)

    def line(self):
        """Its line as a point and a direction."""
        return self.p, (self.q[0]-self.p[0], self.q[1]-self.p[1])

    def line_params(self, base, e):
        """Parameters `lam` where `base + lam e` crosses this segment (open)."""
        d = (self.q[0]-self.p[0], self.q[1]-self.p[1])
        den = e[0]*d[1]-e[1]*d[0]
        if den == 0:
            return []
        w = (self.p[0]-base[0], self.p[1]-base[1])
        lam = (w[0]*d[1]-w[1]*d[0])/den
        t = (w[0]*e[1]-w[1]*e[0])/den
        return [(lam, t)] if 0 < t < 1 else []


def inv2(A):
    a, b, c, d = A
    det = a*d-b*c
    return (d/det, -b/det, -c/det, a/det)


def mat2(A, v):
    return (A[0]*v[0]+A[1]*v[1], A[2]*v[0]+A[3]*v[1])


class Arc:
    """`C + A (cos t, sin t)` (`A = (a, b, c, d)`: `x = C0 + a cos + b sin`,
    `y = C1 + c cos + d sin`) for `t` from `t0` to `t1` (either way)."""
    kind = 'arc'

    def __init__(self, C, A, t0, t1, tag, sphere=None):
        self.C, self.A, self.t0, self.t1, self.tag, self.sphere = C, A, t0, t1, tag, sphere
        self._inv = None

    def inv(self):
        if self._inv is None:
            self._inv = inv2(self.A)
        return self._inv

    def point(self, t):
        c, s = mp.cos(t), mp.sin(t)
        a, b, cc, d = self.A
        return (self.C[0]+a*c+b*s, self.C[1]+cc*c+d*s)

    def start(self):
        return self.point(self.t0)

    def end(self):
        return self.point(self.t1)

    def lo_hi(self):
        return self.t0, self.t1

    def sub(self, a, b):
        return Arc(self.C, self.A, a, b, self.tag, self.sphere)

    def reversed(self):
        return Arc(self.C, self.A, self.t1, self.t0, self.tag, self.sphere)

    def angle(self):
        return self.t1-self.t0

    def green(self):
        return green_arc(self.C, self.A, self.t0, self.t1)

    def unit(self, X):
        """`A^-1 (X - C)`."""
        return mat2(self.inv(), (X[0]-self.C[0], X[1]-self.C[1]))

    def inside(self, X):
        u = self.unit(X)
        return u[0]*u[0]+u[1]*u[1] < 1

    def line_params(self, base, e):
        """`[(lam, t)]` where `base + lam e` crosses this arc (inside its
        span)."""
        w0 = self.unit(base)
        w1 = mat2(self.inv(), e)
        a = w1[0]**2+w1[1]**2
        b = 2*(w0[0]*w1[0]+w0[1]*w1[1])
        c = w0[0]**2+w0[1]**2-1
        disc = b*b-4*a*c
        if disc <= 0:
            return []
        sq = mp.sqrt(disc)
        out = []
        for lam in ((-b-sq)/(2*a), (-b+sq)/(2*a)):
            u = (w0[0]+lam*w1[0], w0[1]+lam*w1[1])
            t = span_param(self.t0, self.t1, mp.atan2(u[1], u[0]))
            if t is not None:
                out.append((lam, t))
        return out


def trig_ints(t0, t1):
    """Integrals over [t0, t1] of `cos^i sin^j` (i + j <= 3)."""
    c0, s0, c1, s1_ = mp.cos(t0), mp.sin(t0), mp.cos(t1), mp.sin(t1)
    D = lambda f: f(c1, s1_, t1)-f(c0, s0, t0)
    return {(0, 0): t1-t0, (1, 0): s1_-s0, (0, 1): c0-c1,
            (2, 0): D(lambda c, s, t: t/2+s*c/2), (0, 2): D(lambda c, s, t: t/2-s*c/2),
            (1, 1): D(lambda c, s, t: s*s/2), (3, 0): D(lambda c, s, t: s-s**3/3),
            (0, 3): D(lambda c, s, t: -c+c**3/3), (2, 1): D(lambda c, s, t: -c**3/3),
            (1, 2): D(lambda c, s, t: s**3/3)}


def tmul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = (i+k, j+l)
            out[key] = out.get(key, 0)+x*y
    return out


def green_arc(C, A, t0, t1):
    """Green's integrals `(x dy - y dx)/2`, `x^2 dy / 2`, `-y^2 dx / 2`
    over `C + A (cos t, sin t)`, `t` from `t0` to `t1`, exactly (the
    integrands polynomials of degree three in cos and sin)."""
    a1, b1, a2, b2 = A
    I = trig_ints(t0, t1)
    area = (C[0]*(b2*I[(1, 0)]-a2*I[(0, 1)])-C[1]*(b1*I[(1, 0)]-a1*I[(0, 1)])+(a1*b2-a2*b1)*I[(0, 0)])/2
    X = {(0, 0): C[0], (1, 0): a1, (0, 1): b1}
    Y = {(0, 0): C[1], (1, 0): a2, (0, 1): b2}
    dX = {(0, 1): -a1, (1, 0): b1}
    dY = {(0, 1): -a2, (1, 0): b2}
    mx = sum(v*I[k] for k, v in tmul(tmul(X, X), dY).items())/2
    my = -sum(v*I[k] for k, v in tmul(tmul(Y, Y), dX).items())/2
    return (area, mx, my)


def float_roots(coeffs):
    """Approximate roots of a complex polynomial (descending coefficients)
    by Aberth's iteration in binary64."""
    c = [complex(x) for x in coeffs]
    n = len(c)-1
    c = [x/c[0] for x in c]
    rad = 1+max(abs(x) for x in c[1:])
    z = [rad*0.5*complex(math.cos(2*math.pi*k/n+0.4), math.sin(2*math.pi*k/n+0.4)) for k in range(n)]
    for _ in range(500):
        moved = 0.0
        for k in range(n):
            p = dp = 0j
            for x in c:
                dp = dp*z[k]+p
                p = p*z[k]+x
            if p == 0:
                continue
            ratio = p/dp if dp != 0 else 1e-3
            corr = sum(1/(z[k]-z[j]) for j in range(n) if j != k and z[k] != z[j])
            w = ratio/(1-ratio*corr)
            z[k] -= w
            moved = max(moved, abs(w))
        if moved < 1e-15*rad:
            break
    return z


def arc_arc_params(p, o):
    """`[(t_p, t_o)]` where arc `p` crosses arc `o`: `f(t) = |w + B u(t)|^2
    - 1` on `p` in `o`'s unit coordinates, a trigonometric polynomial of
    degree two, whose roots are those of the quartic `z^2 f` on the unit
    circle (found in binary64, refined by Newton on `f`, or by
    `mp.polyroots` at 70 digits where two refine to one root)."""
    Oi = o.inv()
    w = mat2(Oi, (p.C[0]-o.C[0], p.C[1]-o.C[1]))
    a, b, c, d = p.A
    B = (Oi[0]*a+Oi[1]*c, Oi[0]*b+Oi[1]*d, Oi[2]*a+Oi[3]*c, Oi[2]*b+Oi[3]*d)
    P = B[0]**2+B[2]**2
    R = B[1]**2+B[3]**2
    Q = B[0]*B[1]+B[2]*B[3]
    al = 2*(w[0]*B[0]+w[1]*B[2])
    be = 2*(w[0]*B[1]+w[1]*B[3])
    c0 = w[0]**2+w[1]**2-1+(P+R)/2
    ga, de = (P-R)/2, Q
    f = lambda t: c0+al*mp.cos(t)+be*mp.sin(t)+ga*mp.cos(2*t)+de*mp.sin(2*t)
    df = lambda t: -al*mp.sin(t)+be*mp.cos(t)-2*ga*mp.sin(2*t)+2*de*mp.cos(2*t)
    coeffs = [mp.mpc(ga, -de)/2, mp.mpc(al, -be)/2, mp.mpc(c0), mp.mpc(al, be)/2, mp.mpc(ga, de)/2]
    scale_ = max(abs(x) for x in coeffs)
    eps = mp.mpf(10)**-32*scale_
    while coeffs and abs(coeffs[0]) <= eps:
        coeffs = coeffs[1:-1]
    if len(coeffs) <= 1:
        return []
    fs = abs(c0)+abs(al)+abs(be)+abs(ga)+abs(de)
    roots = None
    # The binary64 stage without the coefficients below 1e-10 of the
    # largest (a nearly circular pair: its far roots are not on the circle).
    fc = list(coeffs)
    while len(fc) > 2 and abs(fc[0]) <= mp.mpf(10)**-10*scale_:
        fc = fc[1:-1]
    if len(fc) >= 2:
        found = []
        ok = True
        for z in float_roots(fc):
            if abs(abs(z)-1) > 1e-4:
                continue
            t = mp.mpf(math.atan2(z.imag, z.real))
            for _ in range(40):
                dv = df(t)
                if dv == 0:
                    break
                step = f(t)/dv
                t -= step
                if abs(step) < mp.mpf(10)**-38:
                    break
            if abs(f(t)) > mp.mpf(10)**-34*fs:
                continue
            t = mp.atan2(mp.sin(t), mp.cos(t))
            if any(abs(t-u) < mp.mpf(10)**-20 for u in found):
                ok = False
                break
            found.append(t)
        if ok and len(found) % 2 == 0:
            roots = found
    if roots is None:
        with mp.workdps(70):
            try:
                zs = mp.polyroots(coeffs, maxsteps=400, extraprec=300)
            except mp.NoConvergence:
                zs = mp.polyroots(coeffs, maxsteps=4000, extraprec=1000)
        roots = []
        for z in zs:
            if abs(abs(z)-1) > mp.mpf(10)**-20:
                continue
            t = mp.arg(z)
            for _ in range(4):
                dv = df(t)
                if dv == 0:
                    break
                t -= f(t)/dv
            roots.append(t)
    out = []
    for t in roots:
        tp = span_param(p.t0, p.t1, t)
        if tp is None:
            continue
        X = p.point(t)
        u = o.unit(X)
        to = span_param(o.t0, o.t1, mp.atan2(u[1], u[0]))
        if to is not None:
            out.append((tp, to))
    return out


def crossings(x, y):
    """Parameters on piece `x` where it crosses piece `y` (open ranges)."""
    if x.kind == 'seg':
        base, e = x.line()
        return [lam for lam, _ in y.line_params(base, e) if 0 < lam < 1]
    if y.kind == 'seg':
        base, e = y.line()
        return [t for lam, t in x.line_params(base, e) if 0 < lam < 1]
    return [tp for tp, _ in arc_arc_params(x, y)]


def split(piece, params):
    a, b = piece.lo_hi()
    ts = sorted(set(params), reverse=b < a)
    marks = [a]+[t for t in ts if min(a, b) < t < max(a, b)]+[b]
    return [piece.sub(t0, t1) for t0, t1 in zip(marks, marks[1:]) if t1 != t0]


def midpoint(piece):
    a, b = piece.lo_hi()
    return piece.point((a+b)/2)


RAY = None


def ray():
    global RAY
    if RAY is None or RAY[0] != mp.mp.prec:
        a = mp.mpf('0.6180339887498948482045868343656381177203')
        RAY = (mp.mp.prec, (mp.cos(a), mp.sin(a)))
    return RAY[1]


def parity(pieces, X):
    """Whether `X` is inside the closed loops `pieces` by the parity of a
    ray's crossings (the ray of irrational slope)."""
    r = ray()
    count = 0
    for p in pieces:
        count += sum(1 for lam, _ in p.line_params(X, r) if lam > 0)
    return count % 2 == 1


# ------------------------------------------------------------------ regions

class Region:
    """A base region (closed loops of pieces, the region on their left, and
    its membership test) cut by half-planes `al x + be y > ga` (with tags);
    `empty` when it holds nothing."""

    def __init__(self, base, contains_base, halfplanes=(), empty=False):
        self.base, self.contains_base = base, contains_base
        self.halfplanes = list(halfplanes)
        self.is_empty = empty
        self._boundary = None

    def in_half(self, X, skip=None):
        return all(al*X[0]+be*X[1] > ga for k, (al, be, ga, _) in enumerate(self.halfplanes) if k != skip)

    def contains(self, X):
        if self.is_empty:
            return False
        return self.in_half(X) and self.contains_base(X)

    def boundary(self):
        if self._boundary is not None:
            return self._boundary
        out = []
        if self.is_empty:
            self._boundary = out
            return out
        lines = []
        for al, be, ga, tag in self.halfplanes:
            nn = al*al+be*be
            foot = (al*ga/nn, be*ga/nn)
            lines.append((foot, (be, -al), tag))
        for piece in self.base:
            params = []
            for foot, t, _ in lines:
                if piece.kind == 'seg':
                    params += seg_line_params(piece, foot, t)
                else:
                    params += [tt for _, tt in piece.line_params(foot, t)]
            for sub_ in split(piece, params):
                if self.in_half(midpoint(sub_)):
                    out.append(sub_)
        for k, (foot, t, tag) in enumerate(lines):
            lams = []
            for piece in self.base:
                lams += [lam for lam, _ in piece.line_params(foot, t)]
            for j, (f2, t2, _) in enumerate(lines):
                if j != k:
                    den = t[0]*t2[1]-t[1]*t2[0]
                    if den != 0:
                        w = (f2[0]-foot[0], f2[1]-foot[1])
                        lams.append((w[0]*t2[1]-w[1]*t2[0])/den)
            lams.sort()
            for a, b in zip(lams, lams[1:]):
                if b <= a:
                    continue
                m = (a+b)/2
                X = (foot[0]+m*t[0], foot[1]+m*t[1])
                if self.in_half(X, skip=k) and self.contains_base(X):
                    out.append(Seg((foot[0]+a*t[0], foot[1]+a*t[1]), (foot[0]+b*t[0], foot[1]+b*t[1]), tag))
        self._boundary = out
        return out


def seg_line_params(seg, foot, t):
    """Parameters on the segment where it crosses the line `foot + lam t`."""
    d = (seg.q[0]-seg.p[0], seg.q[1]-seg.p[1])
    den = d[0]*t[1]-d[1]*t[0]
    if den == 0:
        return []
    w = (foot[0]-seg.p[0], foot[1]-seg.p[1])
    u = (w[0]*t[1]-w[1]*t[0])/den
    return [u] if 0 < u < 1 else []


class Pieces:
    """A region given by its boundary pieces and a membership test."""

    def __init__(self, pieces, contains):
        self.pieces, self._contains = pieces, contains
        self.is_empty = not pieces

    def boundary(self):
        return self.pieces

    def contains(self, X):
        return self._contains(X)


def classify(D, P):
    """The boundary pieces of `D` with whether each is inside `P`, and of
    `P` with whether each is inside `D`."""
    bd, bp = D.boundary(), P.boundary()
    dpieces, ppieces = [], []
    for piece in bd:
        params = []
        for o in bp:
            params += crossings(piece, o)
        for sp in split(piece, params):
            dpieces.append((sp, P.contains(midpoint(sp))))
    for piece in bp:
        params = []
        for o in bd:
            params += crossings(piece, o)
        for sp in split(piece, params):
            ppieces.append((sp, D.contains(midpoint(sp))))
    return dpieces, ppieces


def zero3():
    return (mp.mpf(0), mp.mpf(0), mp.mpf(0))


def gsum(pieces):
    out = zero3()
    for p in pieces:
        out = s1.vadd(out, p.green())
    return out


def op_pieces(dpieces, ppieces, op):
    """The result's boundary pieces for `op` (common, fuse, D-P, P-D)."""
    if op == 'common':
        return [p for p, i in dpieces if i]+[p for p, i in ppieces if i]
    if op == 'fuse':
        return [p for p, i in dpieces if not i]+[p for p, i in ppieces if not i]
    if op == 'D-P':
        return [p for p, i in dpieces if not i]+[p.reversed() for p, i in ppieces if i]
    return [p for p, i in ppieces if not i]+[p.reversed() for p, i in dpieces if i]


def region_classes(X, Y):
    """Areas of `X` inside and outside `Y`."""
    xp, yp = classify(X, Y)
    xin = gsum([p for p, i in xp if i])[0]
    xout = gsum([p for p, i in xp if not i])[0]
    yin = gsum([p for p, i in yp if i])[0]
    return xin+yin, xout-yin


# ------------------------------------------------------------------ loops and components

def chain(pieces, tol):
    """The pieces joined end to start into closed loops."""
    left = list(pieces)
    loops = []
    while left:
        loop = [left.pop()]
        start = loop[0].start()
        while True:
            e = loop[-1].end()
            if abs(e[0]-start[0])+abs(e[1]-start[1]) <= tol and (len(loop) > 1 or loop[0].kind == 'arc'):
                break
            best, bd = None, None
            for k, p in enumerate(left):
                s = p.start()
                dd = abs(s[0]-e[0])+abs(s[1]-e[1])
                if bd is None or dd < bd:
                    best, bd = k, dd
            assert best is not None and bd <= tol, 'an open loop'
            loop.append(left.pop(best))
        loops.append(loop)
    return loops


def loop_key(loop):
    tags = [p.tag for p in loop]
    merged = []
    for t in tags:
        if not merged or merged[-1] != t:
            merged.append(t)
    if len(merged) > 1 and merged[0] == merged[-1]:
        merged.pop()
    n = len(merged)
    rots = [tuple(merged[i:]+merged[:i]) for i in range(n)]
    return min(rots, key=repr)


def components(pieces, tol):
    """[(key, Pieces region)] of the result's section."""
    loops = chain(pieces, tol)
    outer, holes = [], []
    for loop in loops:
        a = gsum(loop)[0]
        (outer if a > 0 else holes).append((loop, a))
    comps = [[loop, a, []] for loop, a in outer]
    for loop, _ in holes:
        X = midpoint(loop[0])
        owners = [c for c in comps if parity(c[0], X)]
        assert owners, 'a hole outside every outer loop'
        min(owners, key=lambda c: c[1])[2].append(loop)
    out = []
    for loop, _, hs in comps:
        key = tuple(sorted([repr(loop_key(loop))]+[repr(loop_key(h)) for h in hs]))
        allp = list(loop)+[p for h in hs for p in h]
        region = Pieces(allp, (lambda L, H: lambda X: parity(L, X) and not any(parity(h, X) for h in H))(loop, hs))
        out.append((key, region))
    # Components of one key (S9d.2c: a sphere's circle less an oblique
    # cylinder's ellipse crossing it, two crescents of the same faces):
    # numbered by their centres along the ray's direction, which keeps
    # their order while they deform without meeting.
    r = ray()
    by_key = {}
    for i, (key, region) in enumerate(out):
        by_key.setdefault(key, []).append(i)
    for key, idx in by_key.items():
        if len(idx) < 2:
            continue
        where = []
        for i in idx:
            a, mx, my = gsum(out[i][1].pieces)
            where.append(((mx*r[0]+my*r[1])/a, i))
        for k, (_, i) in enumerate(sorted(where)):
            out[i] = (key+(f'#{k}',), out[i][1])
    keys = [k for k, _ in out]
    assert len(keys) == len(set(keys)), 'two components with one key'
    return out


def inner_points(X, n=64):
    """`n` points of a region, uniform in its bounding box by rejection
    (a fixed linear congruential sequence)."""
    xs, ys = [], []
    for p in X.boundary():
        a, b = p.lo_hi()
        for j in range(17):
            q = p.point(a+(b-a)*j/16)
            xs.append(q[0])
            ys.append(q[1])
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    state = 12345
    out = []
    tries = 0
    while len(out) < n and tries < 200*n:
        tries += 1
        state = (state*6364136223846793005+1442695040888963407) % 2**64
        u = mp.mpf(state >> 11)/2**53
        state = (state*6364136223846793005+1442695040888963407) % 2**64
        v = mp.mpf(state >> 11)/2**53
        q = (x0+(x1-x0)*u, y0+(y1-y0)*v)
        if X.contains(q):
            out.append(q)
    return out


def joined(X, Y):
    """Whether most points of the smaller of two sections lie in the other
    (sections on both sides of a breakpoint: their limits there nest when
    they belong to one solid and meet in no area otherwise)."""
    if area_of(X) > area_of(Y):
        X, Y = Y, X
    pts = inner_points(X)
    return sum(1 for q in pts if Y.contains(q)) > len(pts)//2


def area_of(X):
    return gsum(X.boundary())[0]


# ------------------------------------------------------------------ the slicing chart

def kinv(x, y, n):
    """The inverse of the matrix of columns `x, y, n` (Fractions)."""
    det = s1.det3(x, y, n)
    rows = [cross(y, n), cross(n, x), cross(x, y)]
    return [tuple(c/det for c in r) for r in rows]


def apply3(rows, v):
    return tuple(dot(r, v) for r in rows)


class Chart:
    """The slices `d . X = s`, charted `X = O0 + s tau + a e1 + b e2`."""

    def __init__(self, d, e1=None, e2=None, tau=None, O0=None):
        d = tuple(F(v) for v in d)
        self.d = d
        self.dd = dot(d, d)
        self.dn = mp.sqrt(M(self.dd))
        if e1 is None:
            axes = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))]
            a = min(axes, key=lambda v: abs(dot(v, d)))
            e1 = cross(d, a)
            e2 = cross(d, e1)
            tau = scale(d, 1/self.dd)
            O0 = (F(0), F(0), F(0))
        self.e1, self.e2, self.tau, self.O0 = e1, e2, tau, O0
        assert dot(d, e1) == 0 and dot(d, e2) == 0 and dot(d, tau) == 1 and dot(d, O0) == 0
        self.J = abs(M(s1.det3(e1, e2, tau)))
        self.G = [[dot(e1, e1), dot(e1, e2)], [dot(e2, e1), dot(e2, e2)]]
        g = [[M(v) for v in r] for r in self.G]
        l11 = mp.sqrt(g[0][0])
        l12 = g[0][1]/l11
        l22 = mp.sqrt(g[1][1]-l12*l12)
        # G = L^T L, L upper triangular; A = sqrt(rho) L^-1.
        self.Linv = (1/l11, -l12/(l11*l22), mp.mpf(0), 1/l22)
        self.m = (Mv(e1), Mv(e2), Mv(tau), Mv(O0))

    def world(self, s, q):
        e1, e2, tau, O0 = self.m
        return tuple(O0[i]+s*tau[i]+q[0]*e1[i]+q[1]*e2[i] for i in range(3))

    def origin_poly(self):
        """`O0 + s tau` as polynomials."""
        return plin(self.O0, self.tau)


def prism_chart(P):
    d = cross(P.x, P.y)
    dn = dot(d, P.n)
    tau = scale(P.n, 1/dn)
    O0 = sub(P.o, scale(tau, dot(d, P.o)))
    return Chart(d, P.x, P.y, tau, O0)


# ------------------------------------------------------------------ inputs

class AxisSphere(s1.Sphere):
    """S9d.1's sphere, cap or zone with its end planes as the kernel's `Ball`
    reads them: `(X - o) . n = h |n|^2`, through `o + h n` normal to the
    stored axis `n` (`m = n`). In an exact frame this is S9d.1's affine plane
    (`x * y = n` exactly, the same tuples); in a turned frame the two
    differ by the stored axes' rounding (S9d.2c), and the kernel's
    enclosures hold its own."""

    def __init__(self, case):
        super().__init__(case)
        o, _, _, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        self.m = n
        self.mn = dot(n, n)
        mo = dot(n, o)
        self.planes = []
        if self.heights[0] is not None:
            self.planes.append((n, mo+self.heights[0]*self.mn, 'low'))
        if self.heights[1] is not None:
            self.planes.append((scale(n, -1), -(mo+self.heights[1]*self.mn), 'high'))


class Ball:
    """A sphere, cap or zone (`AxisSphere`'s model)."""

    def __init__(self, case, role):
        self.S = AxisSphere(case)
        self.role = role
        self.c, self.r, self.r2 = self.S.c, self.S.r, self.S.r2
        self.planes = self.S.planes

    def setup(self, chart):
        """Its zone planes as bounds on `s` (they must be slices)."""
        self.chart = chart
        self.sbounds = []
        for a, b, _ in self.planes:
            assert is_zero(cross(a, chart.d)), 'S9d.2 reference: a zone plane must be a slice'
            k = dot(a, chart.d)/chart.dd
            self.sbounds.append((k, b))
        Op = chart.origin_poly()
        w = tuple(psub(Op[i], [self.c[i]]) for i in range(3))
        e1, e2 = chart.e1, chart.e2
        g = (padd([F(0)], pdotv(w, tuple([v] for v in e1))), padd([F(0)], pdotv(w, tuple([v] for v in e2))))
        k = psub(pdotv(w, w), [self.r2])
        self.conic = Conic(chart.G, g, k)
        self.cm, self.rm = Mv(self.c), M(self.r)

    def levels(self):
        out = [b/k for k, b in self.sbounds]
        return out

    def extent_poly(self):
        dc = dot(self.chart.d, self.c)
        return [dc*dc-self.r2*self.chart.dd, -2*dc, F(1)]

    def in_zone(self, s):
        return all(M(k)*s > M(b) for k, b in self.sbounds)

    def section(self, s, tag=None, ignore_zone=False):
        """Its section at `s` as a Region (an ellipse)."""
        if not ignore_zone and not self.in_zone(s):
            return Region([], lambda X: False, empty=True)
        G = self.chart.G
        g = [s1_peval(p, s) for p in self.conic.g]
        k = s1_peval(self.conic.k, s)
        det = M(G[0][0])*M(G[1][1])-M(G[0][1])**2
        gi = ((M(G[1][1])*g[0]-M(G[0][1])*g[1])/det, (-M(G[1][0])*g[0]+M(G[0][0])*g[1])/det)
        q0 = (-gi[0], -gi[1])
        rho = g[0]*gi[0]+g[1]*gi[1]-k
        if rho <= 0:
            return Region([], lambda X: False, empty=True)
        sq = mp.sqrt(rho)
        L = self.chart.Linv
        A = tuple(sq*v for v in L)
        arc = Arc(q0, A, mp.mpf(0), 2*mp.pi, tag or (self.role, 'sphere'), sphere=self.role)
        return Region([arc], arc.inside)

    def closed(self):
        return self.S.closed()

    def z_range(self):
        return self.S.z_range()

    def face_area(self):
        lo, hi = self.z_range()
        return 2*mp.pi*self.rm*(hi-lo)


def s1_peval(p, s):
    out = mp.mpf(0)
    for c in reversed(p):
        out = out*s+M(c)
    return out


class ArcPrism:
    """A prism of lines, arcs and circles on its exact model."""

    def __init__(self, case, role):
        self.case, self.role = case, role
        o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(case.frame))
        self.o, self.x, self.y, self.n = o, x, y, n
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        self.det = s1.det3(x, y, n)
        assert self.det > 0
        self.Kinv = kinv(x, y, n)
        self.elements = []   # ('seg', p, q, k) | ('arc', C, r, t0, t1, k, full), region on the left
        self.vertices = []
        k = 0
        for bi, b in enumerate(case.boundaries):
            outer = bi == 0
            pts, _ = stored(b, case.tolerance)
            els = []
            if pts is None:
                cx, cy, r = b.circle
                els.append(['arc', (F(cx), F(cy)), F(r), mp.mpf(0), 2*mp.pi, True])
            else:
                if b.segments is None:
                    points, segs = pts, [None]*len(pts)
                else:
                    points, segs = pts
                ring = [(F(u), F(v)) for u, v in points]
                self.vertices += ring
                N = len(ring)
                for i in range(N):
                    p, q, sg = ring[i], ring[(i+1) % N], segs[i]
                    if sg is None:
                        els.append(['seg', p, q])
                    else:
                        assert not isinstance(sg, Spline), 'S9d.2: lines, arcs and circles'
                        cx, cy, r, ccw = sg
                        C, r = (F(cx), F(cy)), F(r)
                        for pt in (p, q):
                            assert (pt[0]-C[0])**2+(pt[1]-C[1])**2 == r*r, 'an arc end off its circle'
                        a0 = mp.atan2(M(p[1]-C[1]), M(p[0]-C[0]))
                        a1 = mp.atan2(M(q[1]-C[1]), M(q[0]-C[0]))
                        T = 2*mp.pi
                        t1 = a0+(a1-a0) % T if ccw else a0-(a0-a1) % T
                        els.append(['arc', C, r, a0, t1, False])
            if not outer:
                rev = []
                for e in reversed(els):
                    if e[0] == 'seg':
                        rev.append(['seg', e[2], e[1]])
                    else:
                        rev.append(['arc', e[1], e[2], e[4], e[3], e[5]])
                els = rev
            for e in els:
                self.elements.append(tuple(e)+(k, bi))
                k += 1
        # The profile's own pieces in (u, v).
        self.uv_pieces = []
        for e in self.elements:
            if e[0] == 'seg':
                self.uv_pieces.append(Seg(Mv(e[1]), Mv(e[2]), (role, 'wall', e[3])))
            else:
                _, C, r, t0, t1, full, k, _ = e
                self.uv_pieces.append(Arc(Mv(C), (M(r), mp.mpf(0), mp.mpf(0), M(r)), t0, t1, (role, 'wall', k)))

    def world(self, uv, w):
        return tuple(self.o[i]+uv[0]*self.x[i]+uv[1]*self.y[i]+w*self.n[i] for i in range(3))

    def in_profile(self, uv):
        return parity(self.uv_pieces, uv)

    def profile_moments(self):
        return gsum(self.uv_pieces)

    def points3(self):
        return [self.world(v, w) for v in self.vertices for w in (self.lo, self.hi)]

    def setup(self, chart):
        """The map `q -> (u, v, w) = a(s) + B q` of the chart's slices."""
        self.chart = chart
        K = self.Kinv
        self.a0 = apply3(K, sub(chart.O0, self.o))
        self.a1 = apply3(K, chart.tau)
        B1, B2 = apply3(K, chart.e1), apply3(K, chart.e2)
        self.B = B1, B2   # columns
        self.Buv = (B1[0], B2[0], B1[1], B2[1])
        det = self.Buv[0]*self.Buv[3]-self.Buv[1]*self.Buv[2]
        assert det != 0, 'a slice parallel to the prism axis'
        self.flip = det < 0
        self.Buv_inv = tuple(M(v) for v in inv2(self.Buv))
        self.Bw = (B1[2], B2[2])
        self.strip = self.Bw != (0, 0)
        self.slo = None if self.strip else (self.lo-self.a0[2])/self.a1[2]
        self.shi = None if self.strip else (self.hi-self.a0[2])/self.a1[2]
        if not self.strip and self.slo > self.shi:
            self.slo, self.shi = self.shi, self.slo
        # Each circle's conic in q.
        self.conics = []
        G = [[B1[0]*B1[0]+B1[1]*B1[1], B1[0]*B2[0]+B1[1]*B2[1]],
             [B2[0]*B1[0]+B2[1]*B1[1], B2[0]*B2[0]+B2[1]*B2[1]]]
        for e in self.elements:
            if e[0] == 'arc':
                C, r = e[1], e[2]
                au = (plin([self.a0[0]], [self.a1[0]])[0], plin([self.a0[1]], [self.a1[1]])[0])
                dvec = (psub(au[0], [C[0]]), psub(au[1], [C[1]]))
                g = (padd(pk(dvec[0], B1[0]), pk(dvec[1], B1[1])), padd(pk(dvec[0], B2[0]), pk(dvec[1], B2[1])))
                k = psub(padd(pmul(dvec[0], dvec[0]), pmul(dvec[1], dvec[1])), [r*r])
                self.conics.append(Conic(G, g, k))

    def uv_of(self, s, q):
        a = (M(self.a0[0])+s*M(self.a1[0]), M(self.a0[1])+s*M(self.a1[1]))
        B = self.Buv
        return (a[0]+M(B[0])*q[0]+M(B[1])*q[1], a[1]+M(B[2])*q[0]+M(B[3])*q[1])

    def section(self, s):
        if not self.strip and not (M(self.slo) < s < M(self.shi)):
            return Region([], lambda X: False, empty=True)
        a = (M(self.a0[0])+s*M(self.a1[0]), M(self.a0[1])+s*M(self.a1[1]))
        Bi = self.Buv_inv
        pre = lambda uv: mat2(Bi, (uv[0]-a[0], uv[1]-a[1]))
        pieces = []
        for piece in self.uv_pieces:
            if piece.kind == 'seg':
                p = Seg(pre(piece.p), pre(piece.q), piece.tag)
            else:
                C = pre(piece.C)
                A = (Bi[0]*piece.A[0], Bi[1]*piece.A[0], Bi[2]*piece.A[0], Bi[3]*piece.A[0])
                p = Arc(C, A, piece.t0, piece.t1, piece.tag)
            pieces.append(p.reversed() if self.flip else p)
        halfplanes = []
        if self.strip:
            aw = M(self.a0[2])+s*M(self.a1[2])
            bw = (M(self.Bw[0]), M(self.Bw[1]))
            halfplanes = [(bw[0], bw[1], M(self.lo)-aw, (self.role, 'cap', 0)),
                          (-bw[0], -bw[1], aw-M(self.hi), (self.role, 'cap', 1))]
        return Region(pieces, lambda X: self.in_profile(self.uv_of(s, X)), halfplanes)

    def cap_levels(self):
        return [] if self.strip else [self.slo, self.shi]

    def closed(self):
        A, Mu, Mv_ = self.profile_moments()
        h = M(self.hi-self.lo)
        det = M(self.det)
        V = A*h*det
        cen = tuple(M(self.o[i])+(Mu/A)*M(self.x[i])+(Mv_/A)*M(self.y[i])+M((self.lo+self.hi)/2)*M(self.n[i])
                    for i in range(3))
        return V, tuple(V*c for c in cen)


# ------------------------------------------------------------------ breakpoints

def plane_tangent_poly(ball, chart, a, b):
    """S9d.1's quadratic: the plane `a . X = b` tangent to the ball's section."""
    d = chart.d
    if is_zero(cross(a, d)):
        return None
    aa, ad = dot(a, a), dot(a, d)
    det = aa*chart.dd-ad*ad
    dc = dot(d, ball.c)
    r1 = b-dot(a, ball.c)
    return compose(aa, -2*ad*r1, chart.dd*r1*r1-ball.r2*det, -dc, F(1))


def line_sphere_poly(ball, chart, v, e):
    """The line `v + t e` meeting the ball's sphere, in `s` (None when the
    line lies in a slice)."""
    de = dot(chart.d, e)
    if de == 0:
        return None
    w = sub(v, ball.c)
    return compose(dot(e, e), 2*dot(w, e), dot(w, w)-ball.r2, -dot(chart.d, v)/de, 1/de)


def cap_circle_trig(prism, C, r, w, ball):
    """`|P(theta) - c|^2 - R^2` on the prism's circle `C, r` at height `w`."""
    V0 = sub(prism.world(C, w), ball.c)
    P = trig_vec(V0, scale(prism.x, r), scale(prism.y, r))
    return trig_dot(P, P)+Trig({(0, 0): -ball.r2})


def pair_breakpoints(inputs, chart):
    """Levels and exact polynomials whose roots are the slicing's
    breakpoints."""
    levels, polys = [], []
    balls = [x for x in inputs if isinstance(x, Ball)]
    prisms = [x for x in inputs if isinstance(x, ArcPrism)]
    for b in balls:
        levels += b.levels()
        polys.append(b.extent_poly())
    for P in prisms:
        pts = P.points3()
        levels += [dot(chart.d, p) for p in pts]
        levels += P.cap_levels()
    for b in balls:
        for P in prisms:
            N = len(P.vertices)
            # Straight edges: vertical ones and the caps' segments.
            edges = [(P.world(v, P.lo), sub(P.world(v, P.hi), P.world(v, P.lo))) for v in P.vertices]
            planes = []
            for e in P.elements:
                if e[0] == 'seg':
                    p, q = e[1], e[2]
                    for w in (P.lo, P.hi):
                        edges.append((P.world(p, w), sub(P.world(q, w), P.world(p, w))))
                    ew = add(scale(P.x, q[0]-p[0]), scale(P.y, q[1]-p[1]))
                    a = cross(ew, P.n)
                    planes.append((a, dot(a, P.world(p, P.lo))))
                else:
                    C, r = e[1], e[2]
                    for w in (P.lo, P.hi):
                        tr = cap_circle_trig(P, C, r, w, b)
                        for th in tr.angles():
                            X = tuple(M(P.o[i])+(M(C[0])+M(r)*mp.cos(th))*M(P.x[i])
                                      + (M(C[1])+M(r)*mp.sin(th))*M(P.y[i])+M(w)*M(P.n[i]) for i in range(3))
                            polys.append(('value', dot(Mv(chart.d), X)))
                        if P.strip:
                            # The cap's trace tangent to the circle: the
                            # line (d.x) u + (d.y) v = s - d.(o + w n).
                            dx, dy = dot(chart.d, P.x), dot(chart.d, P.y)
                            base = dot(chart.d, P.world((F(0), F(0)), w))+dx*C[0]+dy*C[1]
                            polys.append([base*base-r*r*(dx*dx+dy*dy), -2*base, F(1)])
            cap_normal = cross(P.x, P.y)
            for w in (P.lo, P.hi):
                planes.append((cap_normal, dot(cap_normal, P.world((F(0), F(0)), w))))
            for v, e in edges:
                p = line_sphere_poly(b, chart, v, e)
                if p is not None:
                    polys.append(p)
            for a, bb in planes:
                p = plane_tangent_poly(b, chart, a, bb)
                if p is not None:
                    polys.append(p)
            for cn in P.conics:
                polys.append(tangency_poly(cn, b.conic))
    if len(balls) == 2:
        polys.append(tangency_poly(balls[0].conic, balls[1].conic))
    pts = [M(v) for v in levels]
    for p in polys:
        if isinstance(p, tuple):
            pts.append(p[1])
        else:
            pts += real_roots(p)
    return pts


# ------------------------------------------------------------------ the pair

class Pair:
    """Two inputs (a sphere and a prism in either order, or two spheres)."""

    def __init__(self, obj, tool, d=None):
        make = lambda c, role: Ball(c, role) if c.sphere is not None else ArcPrism(c, role)
        self.A, self.B = make(obj, 'A'), make(tool, 'B')
        balls = [x for x in (self.A, self.B) if isinstance(x, Ball)]
        assert balls, 'S9d.2: a sphere against a prism or a sphere'
        if len(balls) == 2:
            self.D, self.P = self.A, self.B
        else:
            self.D = balls[0]
            self.P = self.B if self.D is self.A else self.A
        self.two_balls = len(balls) == 2
        self.size = self.case_size()
        self.chart = self.make_chart(d)
        self._res = self._faces = None
        self.quad_error = mp.mpf(0)

    def make_chart(self, d):
        oblique = [] if self.two_balls else \
            [a for a, _, _ in self.D.planes if not is_zero(cross(a, cross(self.P.x, self.P.y)))]
        if d is None and oblique:
            # A turned cap's end planes (S9d.2c): the slices along its
            # stored axis, the prism cut obliquely.
            chart = Chart(oblique[0])
        elif d is None and not self.two_balls:
            chart = prism_chart(self.P)
        else:
            if d is None:
                zoned = [b for b in (self.A, self.B) if b.planes]
                d = (zoned[0] if zoned else self.A).S.m
            chart = Chart(d)
        for x in (self.D, self.P):
            x.setup(chart)
        return chart

    def case_size(self):
        vals = [F(1)]
        for x in (self.A, self.B):
            if isinstance(x, Ball):
                vals += [abs(c)+x.r for c in x.c]
            else:
                vals += [abs(c) for p in x.points3() for c in p]
                for e in x.elements:
                    if e[0] == 'arc':
                        C, r = e[1], e[2]
                        for w in (x.lo, x.hi):
                            vals += [abs(c)+r*2 for c in x.world(C, w)]
        return M(max(vals))

    # ---- slicing

    def sections(self, s):
        return self.D.section(s), self.P.section(s)

    def integrand(self, s):
        D, P = self.sections(s)
        dp, pp = classify(D, P)
        din = gsum([p for p, i in dp if i])
        dout = gsum([p for p, i in dp if not i])
        pin = gsum([p for p, i in pp if i])
        pout = gsum([p for p, i in pp if not i])
        out = []
        for g in (din, dout, pin, pout):
            out += [g[0], s*g[0], g[1], g[2]]
        ang = [mp.mpf(0)]*4
        for k, pieces in enumerate((dp, pp)):
            for p, i in pieces:
                if p.kind == 'arc' and p.sphere is not None:
                    a = abs(p.angle())
                    if i:
                        ang[2*k] += a
                    ang[2*k+1] += a
        return out+ang

    def breakpoints(self):
        pts = pair_breakpoints([self.D, self.P], self.chart)
        lo, hi = None, None
        for x in (self.D, self.P):
            if isinstance(x, Ball):
                dc = M(dot(self.chart.d, x.c))
                a, b = dc-x.rm*self.chart.dn, dc+x.rm*self.chart.dn
            else:
                vals = [M(dot(self.chart.d, p)) for p in x.points3()]
                # A cylinder's extreme points along d (its caps' circles).
                for e in x.elements:
                    if e[0] == 'arc':
                        C, r = e[1], e[2]
                        rad = M(r)*mp.sqrt(M(dot(self.chart.d, x.x))**2+M(dot(self.chart.d, x.y))**2)
                        for w in (x.lo, x.hi):
                            c0 = M(dot(self.chart.d, x.world(C, w)))
                            vals += [c0-rad, c0+rad]
                a, b = min(vals), max(vals)
            lo = a if lo is None else min(lo, a)
            hi = b if hi is None else max(hi, b)
        pts = sorted(p for p in pts if lo <= p <= hi)
        self.raw_breaks = pts
        out = []
        gap = mp.mpf(10)**-30*self.size
        for p in [lo]+pts+[hi]:
            if not out or p-out[-1] > gap:
                out.append(p)
        return out

    def measure(self):
        self.breaks = self.breakpoints()
        tol = mp.mpf(10)**-33*self.size**4
        total = None
        for a, b in zip(self.breaks, self.breaks[1:]):
            est, diff = integrate(self.integrand, a, b, tol)
            if est is None:
                continue
            self.quad_error = max(self.quad_error, diff)
            total = est if total is None else [x+y for x, y in zip(total, est)]
        return total

    def sliced(self):
        """{name: (volume, moments)} for `D`, `P`, `common`, `fuse`, `D-P`,
        `P-D`, and the spheres' faces' areas inside and outside the other."""
        if self._res is not None:
            return self._res
        v = self.measure()
        G = {name: v[4*k:4*k+4] for k, name in enumerate(('din', 'dout', 'pin', 'pout'))}
        combos = {'D': [('din', 1), ('dout', 1)], 'P': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'D-P': [('dout', 1), ('pin', -1)], 'P-D': [('pout', 1), ('din', -1)]}
        ch = self.chart
        e1, e2, tau, O0 = ch.m
        out = {}
        for name, parts in combos.items():
            A = sum(k*G[p][0] for p, k in parts)
            sA = sum(k*G[p][1] for p, k in parts)
            ma = sum(k*G[p][2] for p, k in parts)
            mb = sum(k*G[p][3] for p, k in parts)
            vol = ch.J*A
            mom = tuple(ch.J*(O0[i]*A+tau[i]*sA+e1[i]*ma+e2[i]*mb) for i in range(3))
            out[name] = (vol, mom)
        for k, x in enumerate((self.D, self.P)):
            if isinstance(x, Ball):
                f = x.rm/ch.dn
                th_in, th_all = v[16+2*k], v[17+2*k]
                out[('sphere', x.role)] = {'in': f*th_in, 'out': f*(th_all-th_in), 'total': f*th_all}
        self._res = out
        return out

    # ---- faces

    def faces(self):
        """[(role, tag, classes, exact area)] of both inputs' faces."""
        if self._faces is not None:
            return self._faces
        out = []
        res = self.sliced()
        for x, y in ((self.D, self.P), (self.P, self.D)):
            if isinstance(x, Ball):
                sf = res[('sphere', x.role)]
                out.append((x.role, ('sphere',), {'in': sf['in'], 'out': sf['out'], 'same': mp.mpf(0),
                                                   'opp': mp.mpf(0)}, x.face_area()))
                out += self.end_discs(x, y)
            else:
                out += self.prism_faces(x, y)
        self._faces = out
        return out

    def coplanar(self, y, s):
        """`y`'s face in the slice `s` exactly: (outward normal sign along
        d, region) or None."""
        if isinstance(y, Ball):
            for (k, b), (a, bb, _) in zip(y.sbounds, y.planes):
                if b/k == s:
                    return (-1 if dot(a, self.chart.d) > 0 else 1), y.section(M(s), ignore_zone=True)
            return None
        if not y.strip:
            for lv, sign in ((y.slo, -1), (y.shi, 1)):
                if lv == s:
                    # The cap's region: the profile's preimage at that slice.
                    region = ArcPrismCap(y, M(s))
                    # Outward normal of the low cap points to -d when the
                    # map keeps the order of s and w.
                    return (sign, region)
        return None

    def end_discs(self, x, y):
        out = []
        for (k, b), (a, bb, end) in zip(x.sbounds, x.planes):
            s = b/k
            disc = x.section(M(s), ignore_zone=True)
            area = area_of(disc)*self.face_scale()
            cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
            outward = -1 if dot(a, self.chart.d) > 0 else 1
            co = self.coplanar(y, s)
            if co is not None:
                sign, region = co
                inside, outside = region_classes(disc, region)
                cls['same' if sign == outward else 'opp'] = inside*self.face_scale()
                cls['out'] = outside*self.face_scale()
            else:
                other = y.section(M(s))
                if other.is_empty:
                    cls['out'] = area
                else:
                    inside, outside = region_classes(disc, other)
                    cls['in'], cls['out'] = inside*self.face_scale(), outside*self.face_scale()
            out.append((x.role, ('end', end), cls, area))
        return out

    def face_scale(self):
        """Area in 3D per unit area of the chart."""
        ch = self.chart
        c = cross(ch.e1, ch.e2)
        return mp.sqrt(M(dot(c, c)))

    def prism_faces(self, P, ball):
        if P.strip:
            return self.oblique_prism_faces(P, ball)
        out = []
        # Caps.
        for lv, cap in ((P.slo, 0), (P.shi, 1)):
            region = ArcPrismCap(P, M(lv))
            area = area_of(region)*self.face_scale()
            cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
            co = self.coplanar(ball, lv)
            outward = -1 if cap == 0 else 1
            if co is not None:
                sign, disc = co
                inside, outside = region_classes(region, disc)
                cls['same' if sign == outward else 'opp'] = inside*self.face_scale()
                cls['out'] = outside*self.face_scale()
            else:
                other = ball.section(M(lv))
                if other.is_empty:
                    cls['out'] = area
                else:
                    inside, outside = region_classes(region, other)
                    cls['in'], cls['out'] = inside*self.face_scale(), outside*self.face_scale()
            out.append((P.role, ('cap', cap), cls, area))
        for e in P.elements:
            if e[0] == 'seg':
                out.append(self.flat_wall(P, ball, e))
            else:
                out.append(self.round_wall(P, ball, e))
        return out

    def flat_wall(self, P, ball, e):
        p, q, k = e[1], e[2], e[3]
        ew = add(scale(P.x, q[0]-p[0]), scale(P.y, q[1]-p[1]))
        X0 = P.world(p, F(0))
        cn = mp.sqrt(M(dot(cross(ew, P.n), cross(ew, P.n))))
        lo, hi = M(P.lo), M(P.hi)
        rect = [Seg((mp.mpf(0), lo), (mp.mpf(1), lo), 'r0'), Seg((mp.mpf(1), lo), (mp.mpf(1), hi), 'r1'),
                Seg((mp.mpf(1), hi), (mp.mpf(0), hi), 'r2'), Seg((mp.mpf(0), hi), (mp.mpf(0), lo), 'r3')]
        R = Region(rect, lambda X: 0 < X[0] < 1 and lo < X[1] < hi)
        area = cn*(hi-lo)
        cls = {'in': mp.mpf(0), 'out': area, 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        disc = ball_region(ball, X0, ew, P.n)
        if disc is not None:
            inside, outside = region_classes(R, disc)
            cls['in'], cls['out'] = inside*cn, outside*cn
        return (P.role, ('wall', k), cls, area)

    def oblique_prism_faces(self, P, ball):
        """A prism's faces when the slices cut it obliquely (S9d.2c): each cap
        in its own `(u, v)` against the ball's ellipse and zone half-planes
        in its plane, flat walls as `flat_wall`, cylindrical walls by
        `oblique_round_wall`."""
        out = []
        cn = cross(P.x, P.y)
        k_uv = mp.sqrt(M(dot(cn, cn)))
        profile = Region(P.uv_pieces, P.in_profile)
        parea = area_of(profile)
        for w, cap in ((P.lo, 0), (P.hi, 1)):
            X0 = P.world((F(0), F(0)), w)
            for a, b, _ in ball.planes:
                assert not (is_zero(cross(a, cn)) and dot(a, X0) == b), 'S9d.2c reference: a prism cap on a zone plane'
            area = parea*k_uv
            cls = {'in': mp.mpf(0), 'out': area, 'same': mp.mpf(0), 'opp': mp.mpf(0)}
            region = ball_region(ball, X0, P.x, P.y)
            if region is not None:
                inside, outside = region_classes(profile, region)
                cls['in'], cls['out'] = inside*k_uv, outside*k_uv
            out.append((P.role, ('cap', cap), cls, area))
        for e in P.elements:
            if e[0] == 'seg':
                out.append(self.flat_wall(P, ball, e))
            else:
                out.append(self.oblique_round_wall(P, ball, e))
        return out

    def oblique_round_wall(self, P, ball, e):
        """`round_wall` with zone planes of any direction (S9d.2c): along the
        generatrix `P(theta) + w n` the plane `a . X > b` bounds `w` by
        `(b - a . P(theta)) / (a . n)`, from below where `a . n > 0`."""
        _, C, r, t0, t1, full, k, _ = e
        x, y, n = P.x, P.y, P.n
        base = P.world(C, F(0))
        V0 = sub(base, ball.c)
        Pt = trig_vec(V0, scale(x, r), scale(y, r))
        nn = dot(n, n)
        l = trig_dot(tuple(Trig({(0, 0): c}) for c in n), Pt)
        disc = l*l+(trig_dot(Pt, Pt)+Trig({(0, 0): -ball.r2})).scale(-nn)
        events = disc.angles()
        for w in (P.lo, P.hi):
            events += cap_circle_trig(P, C, r, w, ball).angles()
        bounds = []
        for a, b, _ in ball.planes:
            an = dot(a, n)
            assert an != 0, 'S9d.2c reference: a zone plane along a generatrix'
            # a . n times the bound: b - a . base - r (a . x cos + a . y sin).
            lin = Trig.linear(b-dot(a, base), -r*dot(a, x), -r*dot(a, y))
            bounds.append((M(an), M(b-dot(a, base)), M(-r*dot(a, x)), M(-r*dot(a, y))))
            for w in (P.lo, P.hi):
                events += (lin+Trig({(0, 0): -an*w})).angles()
            # The bound on the sphere: the zone's rim crossing the wall.
            rim = trig_dot(Pt, Pt).scale(an*an)+(l*lin).scale(2*an)+(lin*lin).scale(nn) \
                + Trig({(0, 0): -ball.r2*an*an})
            events += rim.angles()
        a0, a1 = min(t0, t1), max(t0, t1)
        pts = [a0, a1]
        T = 2*mp.pi
        for th in events+[mp.pi]:
            b0 = th+T*mp.floor((a0-th)/T)
            for j in range(0, 4):
                v = b0+j*T
                if a0 < v < a1:
                    pts.append(v)
        pts = sorted(pts)
        xm, ym, nm = Mv(x), Mv(y), Mv(n)
        V0m = Mv(V0)
        rm = M(r)
        nnm = M(nn)
        wlo0, whi0 = M(P.lo), M(P.hi)

        def f(th):
            c, s_ = mp.cos(th), mp.sin(th)
            t = tuple(-s_*xm[i]+c*ym[i] for i in range(3))
            tn = cross(t, nm)
            el = rm*mp.sqrt(dot(tn, tn))
            Pv = tuple(V0m[i]+rm*(c*xm[i]+s_*ym[i]) for i in range(3))
            b = dot(nm, Pv)
            dd = b*b-nnm*(dot(Pv, Pv)-M(ball.r2))
            wlo, whi = wlo0, whi0
            for an, k0, kc, ks in bounds:
                w = (k0+kc*c+ks*s_)/an
                if an > 0:
                    wlo = max(wlo, w)
                else:
                    whi = min(whi, w)
            inside = mp.mpf(0)
            if dd > 0:
                sq = mp.sqrt(dd)
                w1, w2 = (-b-sq)/nnm, (-b+sq)/nnm
                lo_, hi_ = max(w1, wlo), min(w2, whi)
                if hi_ > lo_:
                    inside = hi_-lo_
            return [el*inside, el]
        tol = mp.mpf(10)**-33*self.size**2
        tot = [mp.mpf(0), mp.mpf(0)]
        for a, b in zip(pts, pts[1:]):
            if b-a <= mp.mpf(10)**-35:
                continue
            est, _ = integrate(f, a, b, tol)
            tot = [u+v for u, v in zip(tot, est)]
        h = M(P.hi-P.lo)
        area = tot[1]*h
        cls = {'in': tot[0], 'out': area-tot[0], 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        return (P.role, ('wall', k), cls, area)

    def round_wall(self, P, ball, e):
        _, C, r, t0, t1, full, k, _ = e
        x, y, n = P.x, P.y, P.n
        V0 = sub(P.world(C, F(0)), ball.c)
        V1, V2 = scale(x, r), scale(y, r)
        Pt = trig_vec(V0, V1, V2)
        nn = dot(n, n)
        l = trig_dot(tuple(Trig({(0, 0): c}) for c in n), Pt)
        disc = l*l+(trig_dot(Pt, Pt)+Trig({(0, 0): -ball.r2})).scale(-nn)
        # Zone bounds along the generatrix: a . (P + w n) > b, a . P constant.
        wl, wh = [M(P.lo)], [M(P.hi)]
        for a, b, _ in ball.planes:
            an = dot(a, n)
            aP = dot(a, P.world(C, F(0)))
            assert dot(a, x) == 0 and dot(a, y) == 0
            bound = (b-aP)/an
            (wl if an > 0 else wh).append(M(bound))
        wlo, whi = max(wl), min(wh)
        events = disc.angles()
        heights = [P.lo, P.hi]
        for a, b, _ in ball.planes:
            heights.append((b-dot(a, P.world(C, F(0))))/dot(a, n))
        for w in heights:
            events += cap_circle_trig(P, C, r, w, ball).angles()
        a0, a1 = min(t0, t1), max(t0, t1)
        pts = [a0, a1]
        T = 2*mp.pi
        for th in events+[mp.pi]:
            base = th+T*mp.floor((a0-th)/T)
            for j in range(0, 4):
                v = base+j*T
                if a0 < v < a1:
                    pts.append(v)
        pts = sorted(pts)
        xm, ym, nm = Mv(x), Mv(y), Mv(n)
        V0m = Mv(V0)
        rm = M(r)
        nnm = M(nn)

        def f(th):
            c, s_ = mp.cos(th), mp.sin(th)
            t = tuple(-s_*xm[i]+c*ym[i] for i in range(3))
            tn = cross(t, nm)
            el = rm*mp.sqrt(dot(tn, tn))
            Pv = tuple(V0m[i]+rm*(c*xm[i]+s_*ym[i]) for i in range(3))
            b = dot(nm, Pv)
            dd = b*b-nnm*(dot(Pv, Pv)-M(ball.r2))
            inside = mp.mpf(0)
            if dd > 0:
                sq = mp.sqrt(dd)
                w1, w2 = (-b-sq)/nnm, (-b+sq)/nnm
                lo_, hi_ = max(w1, wlo), min(w2, whi)
                if hi_ > lo_:
                    inside = hi_-lo_
            return [el*inside, el]
        tol = mp.mpf(10)**-33*self.size**2
        tot = [mp.mpf(0), mp.mpf(0)]
        for a, b in zip(pts, pts[1:]):
            if b-a <= mp.mpf(10)**-35:
                continue
            est, _ = integrate(f, a, b, tol)
            tot = [u+v for u, v in zip(tot, est)]
        h = M(P.hi-P.lo)
        area = tot[1]*h
        cls = {'in': tot[0], 'out': area-tot[0], 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        return (P.role, ('wall', k), cls, area)

    # ---- results

    def classes(self, role):
        tot = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        for w, _, cls, _ in self.faces():
            if w == role:
                for k in tot:
                    tot[k] += cls[k]
        return tot

    def area(self, op):
        ca, cb = self.classes('A'), self.classes('B')
        keep = {'fuse': (('out', 'same'), ('out',)), 'cut': (('out', 'opp'), ('in',)),
                'common': (('in', 'same'), ('in',))}[op]
        return sum(ca[k] for k in keep[0])+sum(cb[k] for k in keep[1])

    def slice_op(self, op):
        """The sliced name of `op` (the object is `A`)."""
        if op != 'cut':
            return op
        return 'D-P' if self.D is self.A else 'P-D'

    def volume(self, op):
        return self.sliced()[self.slice_op(op)]

    def solids(self, op):
        return self.components()[self.slice_op(op)]

    def components(self):
        """Solid counts of common, fuse, D-P and P-D (see the module's note)."""
        if getattr(self, '_solids', None) is not None:
            return self._solids
        breaks = self.breaks if hasattr(self, 'breaks') else self.breakpoints()
        tol = mp.mpf(10)**-24*self.size
        names = ('common', 'fuse', 'D-P', 'P-D')
        ends = []   # per interval: {op: (comps at start, comps at end)}
        for a, b in zip(breaks, breaks[1:]):
            delta = min(mp.mpf(10)**-12*self.size, (b-a)/8)
            ends.append({})
            per = {}
            for side, s in (('lo', a+delta), ('hi', b-delta)):
                D, P = self.sections(s)
                dp, pp = classify(D, P)
                for op in names:
                    pieces = op_pieces(dp, pp, op)
                    comps = components(pieces, tol) if pieces else []
                    comps = [(k, r) for k, r in comps if area_of(r) > mp.mpf(10)**-30*self.size**2]
                    per.setdefault(op, {})[side] = comps
            for op in names:
                lo_keys = sorted(k for k, _ in per[op]['lo'])
                hi_keys = sorted(k for k, _ in per[op]['hi'])
                assert lo_keys == hi_keys, ('components change inside an interval', op, lo_keys, hi_keys)
                ends[-1][op] = per[op]
        counts = {}
        for op in names:
            parent = {}

            def find(x):
                while parent[x] != x:
                    x = parent[x]
                return x
            for i, e in enumerate(ends):
                for k, _ in e[op]['lo']:
                    parent[(i, k)] = (i, k)
            for i in range(len(ends)-1):
                for ka, ra in ends[i][op]['hi']:
                    for kb, rb in ends[i+1][op]['lo']:
                        if joined(ra, rb):
                            parent[find((i, ka))] = find((i+1, kb))
            counts[op] = len({find(x) for x in parent})
        self._solids = counts
        return counts

    def result(self, op):
        vol, mom = self.volume(op)
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        n = self.solids(op)
        assert n > 0, 'a result of positive volume without solids'
        return n, vol, self.area(op), tuple(m/vol for m in mom)


def ball_region(ball, X0, e1, e2):
    """The ball in the plane `X0 + a e1 + b e2` as a Region in `(a, b)`: the
    sphere's ellipse cut by the zone's half-planes, or None where the plane
    misses the sphere."""
    G = [[dot(e1, e1), dot(e1, e2)], [dot(e2, e1), dot(e2, e2)]]
    w = sub(X0, ball.c)
    g = [dot(e1, w), dot(e2, w)]
    kk = dot(w, w)-ball.r2
    Gm = [[M(v) for v in r] for r in G]
    det = Gm[0][0]*Gm[1][1]-Gm[0][1]**2
    gm = [M(v) for v in g]
    gi = ((Gm[1][1]*gm[0]-Gm[0][1]*gm[1])/det, (-Gm[1][0]*gm[0]+Gm[0][0]*gm[1])/det)
    rho = gm[0]*gi[0]+gm[1]*gi[1]-M(kk)
    if rho <= 0:
        return None
    l11 = mp.sqrt(Gm[0][0])
    l12 = Gm[0][1]/l11
    l22 = mp.sqrt(Gm[1][1]-l12*l12)
    sq = mp.sqrt(rho)
    A = (sq/l11, -sq*l12/(l11*l22), mp.mpf(0), sq/l22)
    arc = Arc((-gi[0], -gi[1]), A, mp.mpf(0), 2*mp.pi, 'sphere')
    halfplanes = []
    for a, b, end in ball.planes:
        # a . (X0 + lam e1 + mu e2) > b
        halfplanes.append((M(dot(a, e1)), M(dot(a, e2)), M(b-dot(a, X0)), end))
    return Region([arc], arc.inside, halfplanes)


def ArcPrismCap(P, s):
    """A prism's cap at its slice `s`: the profile's preimage without the
    level test."""
    a = (M(P.a0[0])+s*M(P.a1[0]), M(P.a0[1])+s*M(P.a1[1]))
    Bi = P.Buv_inv
    pre = lambda uv: mat2(Bi, (uv[0]-a[0], uv[1]-a[1]))
    pieces = []
    for piece in P.uv_pieces:
        if piece.kind == 'seg':
            p = Seg(pre(piece.p), pre(piece.q), piece.tag)
        else:
            C = pre(piece.C)
            A = (Bi[0]*piece.A[0], Bi[1]*piece.A[0], Bi[2]*piece.A[0], Bi[3]*piece.A[0])
            p = Arc(C, A, piece.t0, piece.t1, piece.tag)
        pieces.append(p.reversed() if P.flip else p)
    return Region(pieces, lambda X: P.in_profile(P.uv_of(s, X)))


def number(x):
    return s1.number(x)


def rows(obj, operation, tool, pair=None):
    """`result N volume area cx cy cz` or `empty`, and the pair (reused
    across the three operations)."""
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair
