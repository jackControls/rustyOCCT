#!/usr/bin/env python3
"""Independent reference for S9d.3b of REVIEW_NOTES.md: Booleans of a cone
or frustum (`Solid::cone_with`) against a prism whose profile holds lines,
arcs and circles (planar and cylindrical walls), against a sphere, cap or
zone (`Solid::sphere_with`), and against another cone or frustum, in any
relative position. A cone and a cylinder, a sphere or another cone meet in
circles (a shared axis) or quartics (S7b.2's procedural curves).

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal): the cone S9d.3a's (`cone_boolean_reference.Cone`: the points
`o + u x + v y + w n` of its stored frame with `0 <= w <= h` and `u^2 + v^2
<= (bottom + k w)^2`), the sphere S9d.1's (`sphere_boolean_reference.Sphere`,
a cap's or zone's end planes at its stored heights; S9d.3c: a turned cap's
end planes as the kernel reads them, S9d.2c's `AxisSphere`, its circles and
discs spanned by vectors exactly in them), the prism S9d.2's
(`spheres_boolean_reference.ArcPrism`: its profile's exact circles and
segments, holes reversed). Every input is a set of exact surfaces, each
`Q(X) = X^T M X + 2 m . X + c` in world coordinates (a plane with `M = 0`):
the cone's quadric `(a1 . d)^2 + (a2 . d)^2 - (bottom + k a3 . d)^2` (`d = X
- o`, `a1, a2, a3` the rows of the frame's inverse) and its end planes, the
sphere and its zone planes, the prism's walls (planes and the cylinders of
its arcs, `|(u, v) - C|^2 - r^2`) and caps. Nothing here uses the kernel, a
surface/surface intersection or an arrangement; nothing is shared with
another reference but `stored`, `stored_axes`, S9d.1's number, vector and
quadrature helpers and its sphere model, S9d.2's prism model, 2D pieces
(segments and arcs `C + A (cos t, sin t)`, their classification, Green's
integrals, loops and components), slicing chart, conic pencils and
trigonometric polynomials, and the exact polynomial helpers of
`boolean_reference.py`.

* **Slices (volume, first moments, solids, the sphere's face).** Both solids
  are sliced by planes `d . X = s` in an affine rational chart `X = O0 + s
  tau + a e1 + b e2` (S9d.2's `Chart`), `d` chosen so that every quadric's
  section is an ellipse: a cone's own axis where the other input allows it
  (its chart `e1 = x`, `e2 = y`: its sections circles about the axis in
  their own angle), else a prism's, else a rational combination of the
  inputs' axes (checked exactly: the cone's form positive definite on the
  slices, no slice parallel to a prism's axis, a cap's or zone's planes
  slices). A cone's section is the ellipse `q^T G q + 2 q^T g(s) + k(s) <=
  0` (`G` its form on the slice, positive definite, so one nappe) cut by
  its end planes' lines (or all or nothing where they are slices), taken as
  `C + sqrt(rho) L^-1 (cos t, sin t)` (`G = L^T L`); a sphere's and a
  prism's as S9d.2's. The operations apart as S9d.2's (boundaries cut where
  they cross, pieces classified at their midpoints; Green's theorem over
  segments and arcs exactly), the sphere's face by Archimedes, and where
  the chart is a cone's own the cone's wall inside the other a second way
  (`r(w) N(theta)` over its circle's arcs inside, S9d.3a's element).
  Breakpoints, all roots of exact polynomials in `s` or levels of points
  found from them: every vertex (a prism's, a cone's apex), every edge of
  one input (a line: a prism's vertical edges and cap segments; a conic in
  its plane: a cap's arcs, a cone's rims, a zone's circles) meeting every
  surface of the other (a line a quadratic; a conic in its plane against
  the other surface's trace there, the resultant in one coordinate, a
  quartic), every plane tangent to every quadric's section (a quadratic),
  two sections tangent (S9d.2's pencil discriminant, or concentric sections
  coinciding), a section shrinking to a point (a cone's apex, a sphere's
  poles), and the planes that are slices. Gauss-Legendre between them as
  S9d.1's `integrate` (1e-33 of the case's size to the fourth).
* **Faces (areas inside and outside the other).** Every face but the
  sphere's is swept by a family of lines `X = A(tau) + w B(tau)` on it: a
  cone's wall by its rulings (`A = o + bottom U(theta)`, `B = n + k
  U(theta)`, `U = cos x + sin y`, the element `r(w) |U' x B|`), a
  cylindrical wall by its generatrices, a flat wall by its generatrices
  over the segment's parameter, a planar face (a prism's cap, a cone's end
  disc, a zone's disc) by parallel lines of a skew rational direction in
  its plane. Along each line every surface of both inputs is a quadratic in
  `w` whose coefficients are exact polynomials in `tau` (trigonometric for
  rulings and generatrices of circles); the line is cut at their roots and
  each piece classified at its midpoint (on the face: within the ruled
  face's range or, for a planar face, inside its own solid with its plane
  left out; then inside, outside, or on a coplanar face of the other of the
  same or the opposite orientation). Breakpoints in `tau`, roots of exact
  polynomials: a surface's quadratic of vanishing leading coefficient or
  discriminant (a line tangent to it), and two surfaces' quadratics with a
  common root (their resultant: the line through a point where the traces
  of two surfaces on the face meet, among them every vertex of the
  arrangement on the face and every point where a trace meets the face's
  boundary); a trigonometric polynomial's roots those of its polynomial in
  `tan(theta / 2)`, found at 100 digits. Between them Gauss-Legendre (1e-33
  of the case's size squared).
* **Solids** as S9d.2's: in each interval between breakpoints the result's
  section chained into loops and components keyed by their source faces,
  followed through intervals (components of one key by their centroids,
  the step halved until each is nearest its own: both ends of a rod in one
  tilted slice) and joined across a breakpoint when most of 64 points of
  the smaller section lie in the other, unless both vanish there (below
  1e-20 of the case's size squared at 1e-12 from it: two apexes touching,
  a pinch along a circle), where they meet in a point or a curve.

Two inputs of the same exact model (equal coaxial cones) have the input as
fuse and common and nothing as either cut, every face shared with the same
orientation. `rows(obj, operation, tool)` gives `result N volume area cx cy
cz` (totals over the N solids, world coordinates) or `empty`, as the S9a to
S9d.3a references. A result keeps (object `A`, tool `B`): fuse `A`'s outside
and shared same-orientation pieces and `B`'s outside, cut `A`'s outside and
shared opposite pieces and `B`'s inside, common `A`'s inside and shared same
pieces and `B`'s inside.
"""
from fractions import Fraction as F
import itertools
import math

import mpmath as mp

from curve_surface_reference import stored_axes
from boolean_reference import padd, pdeg, pmul, psub, ptrim, squarefree
import sphere_boolean_reference as s1
from sphere_boolean_reference import M, Mv, add, compose, cross, dot, integrate, is_zero, roots2, scale, sub
import spheres_boolean_reference as s2
from spheres_boolean_reference import (
    Arc, Chart, Conic, Region, Trig, area_of, classify, components, gsum, joined, op_pieces, parity, pk, prism_chart,
    tangency_poly)

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
SKEW = F(2, 7)   # a planar face's sweep direction e_b = SKEW x + y
Z = mp.mpf(0)


# ------------------------------------------------------------------ exact rings in the sweep parameter

class Pol:
    """An exact polynomial in one variable (ascending Fractions) with the
    interface of `Trig`."""

    def __init__(self, c):
        self.c = ptrim([F(v) for v in c]) if c else [F(0)]

    def __add__(self, o):
        return Pol(padd(self.c, o.c))

    def __mul__(self, o):
        return Pol(pmul(self.c, o.c))

    def scale(self, k):
        return Pol([v*k for v in self.c])

    def is_zero(self):
        return pdeg(self.c) < 0

    def value(self, x):
        out = Z
        for v in reversed(self.c):
            out = out*x+M(v)
        return out

    def poly(self):
        return self.c


def const(ring, v):
    return Trig({(0, 0): F(v)}) if ring == 'trig' else Pol([v])


def neg(x):
    return x.scale(F(-1))


def numeric_roots(p, dps=100):
    """The real roots (mpf) of an exact polynomial (ascending): its content
    removed and zero roots split off exactly, the rest by `mp.polyroots` at
    `dps` digits (a root of multiplicity two still to 50), an imaginary part
    below 1e-25 counting as real (a spurious breakpoint costs nothing)."""
    p = ptrim([F(v) for v in p])
    if pdeg(p) <= 0:
        return []
    out = []
    while p[0] == 0:
        out.append(Z)
        p = p[1:]
    if len(p) <= 1:
        return out
    if len(p) == 2:
        return out+[-M(p[0])/M(p[1])]
    if len(p) == 3:
        return out+roots2(p)
    den = 1
    for v in p:
        den = den*v.denominator//math.gcd(den, v.denominator)
    ints = [int(v*den) for v in p]
    g = 0
    for v in ints:
        g = math.gcd(g, v)
    ints = [v//g for v in ints]
    with mp.workdps(dps):
        desc = [mp.mpf(v) for v in reversed(ints)]
        try:
            found = mp.polyroots(desc, maxsteps=800, extraprec=4*dps)
        except mp.NoConvergence:
            found = mp.polyroots(desc, maxsteps=8000, extraprec=16*dps)
        for r in found:
            re, im = mp.re(r), mp.im(r)
            if abs(im) <= mp.mpf(10)**-25*(1+abs(re)):
                out.append(+re)
    return [+x for x in out]


def ring_roots(x, lo, hi, ring):
    """Roots of a ring element in `[lo, hi]` (angles mapped by whole turns)."""
    if x.is_zero():
        return []
    if ring == 'trig':
        out = [2*mp.atan(t) for t in numeric_roots(x.tan_half())]
        if sum((v*(-1)**i for (i, j), v in x.t.items() if j == 0), F(0)) == 0:
            out.append(+mp.pi)
        T = 2*mp.pi
        res = []
        for th in out:
            base = th+T*mp.floor((lo-th)/T)
            for j in range(0, 4):
                v = base+j*T
                if lo <= v <= hi:
                    res.append(v)
        return res
    return [t for t in numeric_roots(x.poly()) if lo <= t <= hi]


def resultant(a, b):
    """The resultant in `w` of two quadratics `(c0, c1, c2)` over a ring, at
    their formal degrees (None when either has no `w`)."""
    da = 2 if not a[2].is_zero() else (1 if not a[1].is_zero() else 0)
    db = 2 if not b[2].is_zero() else (1 if not b[1].is_zero() else 0)
    if da == 0 or db == 0:
        return None
    a0, a1, a2 = a
    b0, b1, b2 = b
    if da == 1 and db == 1:
        return a1*b0+neg(a0*b1)
    if da == 1:
        return b2*a0*a0+neg(b1*a0*a1)+b0*a1*a1
    if db == 1:
        return a2*b0*b0+neg(a1*b0*b1)+a0*b1*b1
    x = a2*b0+neg(a0*b2)
    y = a2*b1+neg(a1*b2)
    z = a1*b0+neg(a0*b1)
    return x*x+neg(y*z)


# ------------------------------------------------------------------ surfaces

class Surf:
    """`Q(X) = X^T M X + 2 m . X + c` (exact); a plane `M = 0` with its solid
    on `Q >= 0` side flagged by `plane`."""

    def __init__(self, Mq, m, c, tag, plane=False):
        self.M = tuple(tuple(F(v) for v in r) for r in Mq)
        self.m = tuple(F(v) for v in m)
        self.c = F(c)
        self.tag, self.plane = tag, plane
        self.Mm = tuple(tuple(M(v) for v in r) for r in self.M)
        self.mm = Mv(self.m)
        self.cm = M(self.c)

    @staticmethod
    def halfspace(a, b, tag):
        """`a . X >= b` as the plane `a . X - b` (the solid where it is
        positive)."""
        z = (F(0),)*3
        return Surf((z, z, z), tuple(F(v)/2 for v in a), -F(b), tag, plane=True)

    def a(self):
        return tuple(2*v for v in self.m)

    def value(self, X):
        """Exact or mp."""
        Mq, m = (self.M, self.m) if isinstance(X[0], F) else (self.Mm, self.mm)
        c = self.c if isinstance(X[0], F) else self.cm
        return sum(Mq[i][j]*X[i]*X[j] for i in range(3) for j in range(3))+2*sum(m[i]*X[i] for i in range(3))+c

    def line(self, P0, e):
        """Exact `(c0, c1, c2)` of `Q(P0 + t e)`."""
        Me = tuple(sum(self.M[i][j]*e[j] for j in range(3)) for i in range(3))
        c2 = dot(e, Me)
        c1 = 2*(dot(P0, Me)+dot(self.m, e))
        return self.value(P0), c1, c2

    def family(self, A, B, ring):
        """`(c0, c1, c2)` of `Q(A + w B)` over the ring (A, B vectors of ring
        elements)."""
        zero = const(ring, 0)

        def form(U, V):
            out = zero
            for i in range(3):
                for j in range(3):
                    if self.M[i][j]:
                        out = out+(U[i]*V[j]).scale(self.M[i][j])
            return out

        def lin(U):
            out = zero
            for i in range(3):
                if self.m[i]:
                    out = out+U[i].scale(self.m[i])
            return out
        c2 = form(B, B)
        c1 = (form(A, B)+lin(B)).scale(F(2))
        c0 = form(A, A)+lin(A).scale(F(2))+const(ring, self.c)
        return c0, c1, c2

    def num_line(self, P, e):
        """mp `(c0, c1, c2)` of `Q(P + t e)`."""
        Mm, mm = self.Mm, self.mm
        Me = tuple(sum(Mm[i][j]*e[j] for j in range(3)) for i in range(3))
        c2 = sum(e[i]*Me[i] for i in range(3))
        c1 = 2*(sum(P[i]*Me[i] for i in range(3))+sum(mm[i]*e[i] for i in range(3)))
        return self.value(P), c1, c2

    def conic(self, chart):
        """Its trace on the chart's slices as S9d.2's `Conic` (quadrics)."""
        e1, e2, tau, O0 = chart.e1, chart.e2, chart.tau, chart.O0
        Mv_ = lambda v: tuple(sum(self.M[i][j]*v[j] for j in range(3)) for i in range(3))
        G = [[dot(e1, Mv_(e1)), dot(e1, Mv_(e2))], [dot(e2, Mv_(e1)), dot(e2, Mv_(e2))]]
        g = tuple(ptrim([dot(e, Mv_(O0))+dot(self.m, e), dot(e, Mv_(tau))]) for e in (e1, e2))
        k = ptrim([self.value(O0), 2*(dot(O0, Mv_(tau))+dot(self.m, tau)), dot(tau, Mv_(tau))])
        return Conic(G, g, k)


def same_plane(p, q):
    """Whether two plane surfaces are one plane: (True, same outward
    orientation) or (False, None)."""
    a, b = p.a(), q.a()
    if not is_zero(cross(a, b)):
        return False, None
    i = next(j for j in range(3) if a[j] != 0)
    r = b[i]/a[i]
    if q.c != r*p.c:
        return False, None
    return True, r > 0


def solve_quadratic(c0, c1, c2):
    """The real roots (mp) of `c2 t^2 + c1 t + c0`, a vanishing leading
    coefficient (relative 1e-35) taken as linear, a discriminant within
    1e-30 (relative) of zero as no crossing."""
    scale_ = abs(c0)+abs(c1)+abs(c2)
    if scale_ == 0:
        return []
    if abs(c2) <= mp.mpf(10)**-35*scale_:
        if abs(c1) <= mp.mpf(10)**-35*scale_:
            return []
        return [-c0/c1]
    disc = c1*c1-4*c2*c0
    if disc <= mp.mpf(10)**-30*(c1*c1+abs(4*c2*c0)):
        # Within 1e-30 (relative) of tangency: a touching line crosses
        # nothing (a ruling of a cone tangent to a sphere all round).
        return []
    sq = mp.sqrt(disc)
    q = -(c1+(sq if c1 >= 0 else -sq))/2
    out = [q/c2]
    if q != 0:
        out.append(c0/q)
    return out


# ------------------------------------------------------------------ inputs

def frame_rows(case):
    o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
    det = s1.det3(x, y, n)
    assert det > 0
    inv = (scale(cross(y, n), 1/det), scale(cross(n, x), 1/det), scale(cross(x, y), 1/det))
    return o, x, y, n, inv


def outer(a, b, k=F(1)):
    return tuple(tuple(k*a[i]*b[j] for j in range(3)) for i in range(3))


def madd(*ms):
    return tuple(tuple(sum(m[i][j] for m in ms) for j in range(3)) for i in range(3))


def mvec(Mq, v):
    return tuple(sum(Mq[i][j]*v[j] for j in range(3)) for i in range(3))


def own_chart(o, x, y, n):
    """The chart of the slices normal to `x * y`, `(a, b)` the frame's own
    `(u, v)` (S9d.2's `prism_chart`)."""
    d = cross(x, y)
    tau = scale(n, 1/dot(d, n))
    O0 = sub(o, scale(tau, dot(d, o)))
    return Chart(d, x, y, tau, O0)


class Input:
    """Common parts: surfaces, containment with a surface left out, edges
    (lines and conics in their planes), vertices."""
    role = None

    def plane_surfaces(self):
        return [(i, s) for i, s in enumerate(self.surfaces) if s.plane]

    def quad_surfaces(self):
        return [(i, s) for i, s in enumerate(self.surfaces) if not s.plane]


class ConeIn(Input):
    """A cone or frustum on S9d.3a's model."""
    kind = 'cone'

    def __init__(self, case, role):
        self.case, self.role = case, role
        o, x, y, n, inv = frame_rows(case)
        r0, r1, h = (F(v) for v in case.cone)
        assert h > 0 and r0 >= 0 and r1 >= 0 and r0 != r1, 'a cone or a frustum'
        self.o, self.x, self.y, self.n, self.inv = o, x, y, n, inv
        self.r0, self.r1, self.h = r0, r1, h
        self.k = (r1-r0)/h
        a1, a2, a3 = inv
        Mq = madd(outer(a1, a1), outer(a2, a2), outer(a3, a3, -self.k*self.k))
        m = sub(scale(mvec(Mq, o), -1), scale(a3, r0*self.k))
        c = dot(o, mvec(Mq, o))+2*r0*self.k*dot(a3, o)-r0*r0
        self.surfaces = [Surf(Mq, m, c, (role, 'cone')),
                         Surf.halfspace(a3, dot(a3, o), (role, 'end', 0)),
                         Surf.halfspace(scale(a3, -1), -(dot(a3, o)+h), (role, 'end', 1))]
        self.om, self.invm = Mv(o), tuple(Mv(r) for r in inv)
        self.chart0 = own_chart(o, x, y, n)
        K = s1.det3(x, y, n)
        self.det = K
        a, b, cc = cross(y, n), cross(x, n), cross(x, y)
        self.wall_vecs = (Mv(a), Mv(b), Mv(cc))
        self.orthonormal = (dot(a, b) == 0 and dot(a, cc) == 0 and dot(b, cc) == 0 and dot(a, a) == dot(b, b))
        self.N0 = mp.sqrt(M(dot(a, a)+self.k*self.k*dot(cc, cc)))
        self.xy = mp.sqrt(M(dot(cc, cc)))

    def uvw(self, X):
        d = tuple(X[i]-self.om[i] for i in range(3))
        return tuple(sum(r[i]*d[i] for i in range(3)) for r in self.invm)

    def contains(self, X, ignore=None):
        u, v, w = self.uvw(X)
        r = M(self.r0)+M(self.k)*w
        tests = (u*u+v*v < r*r, w > 0, w < M(self.h))
        return all(t for i, t in enumerate(tests) if i != ignore)

    def rim(self, w):
        return self.r0+self.k*w

    def world(self, u, v, w):
        return tuple(self.o[i]+u*self.x[i]+v*self.y[i]+w*self.n[i] for i in range(3))

    def edges(self):
        out = []
        for idx, w in ((1, F(0)), (2, self.h)):
            if self.rim(w) > 0:
                out.append(ConicEdge(self.world(F(0), F(0), w), self.x, self.y, self.surfaces[0], self))
        return out

    def vertices(self):
        out = []
        for w in (F(0), self.h):
            if self.rim(w) == 0:
                out.append(self.world(F(0), F(0), w))
        return out

    def extent(self, d):
        dm = Mv(d)
        vals = []
        for w in (F(0), self.h):
            c = dot(dm, Mv(self.world(F(0), F(0), w)))
            rad = M(self.rim(w))*mp.sqrt(M(dot(d, self.x))**2+M(dot(d, self.y))**2)
            vals += [c-rad, c+rad]
        return min(vals), max(vals)

    # ---- slicing

    def setup(self, chart):
        self.chart = chart
        cn = self.surfaces[0].conic(chart)
        G = cn.G
        assert G[0][0] > 0 and G[0][0]*G[1][1]-G[0][1]*G[1][0] > 0, 'a section of the cone not an ellipse'
        self.conic = cn
        g = [[M(v) for v in r] for r in G]
        l11 = mp.sqrt(g[0][0])
        l12 = g[0][1]/l11
        l22 = mp.sqrt(g[1][1]-l12*l12)
        self.Linv = (1/l11, -l12/(l11*l22), Z, 1/l22)
        self.sbounds, self.halfplanes = [], []
        for idx in (1, 2):
            s = self.surfaces[idx]
            a, b = s.a(), -s.c
            ae1, ae2 = dot(a, chart.e1), dot(a, chart.e2)
            if ae1 == 0 and ae2 == 0:
                self.sbounds.append((dot(a, chart.tau), b-dot(a, chart.O0)))
            else:
                self.halfplanes.append((M(ae1), M(ae2), M(b-dot(a, chart.O0)), M(dot(a, chart.tau)), s.tag))
        self.own = chart is self.chart0

    def levels(self):
        return [b/k for k, b in self.sbounds]

    def section(self, s, ignore_bounds=False):
        if not ignore_bounds and not all(M(k)*s > M(b) for k, b in self.sbounds):
            return Region([], lambda X: False, empty=True)
        G = self.conic.G
        g = [s2.s1_peval(p, s) for p in self.conic.g]
        kk = s2.s1_peval(self.conic.k, s)
        det = M(G[0][0])*M(G[1][1])-M(G[0][1])**2
        gi = ((M(G[1][1])*g[0]-M(G[0][1])*g[1])/det, (-M(G[1][0])*g[0]+M(G[0][0])*g[1])/det)
        q0 = (-gi[0], -gi[1])
        rho = g[0]*gi[0]+g[1]*gi[1]-kk
        if rho <= 0:
            return Region([], lambda X: False, empty=True)
        sq = mp.sqrt(rho)
        A = tuple(sq*v for v in self.Linv)
        arc = Arc(q0, A, Z, 2*mp.pi, (self.role, 'cone'))
        hp = [(al, be, ga-s*gt, tag) for al, be, ga, gt, tag in self.halfplanes]
        return Region([arc], arc.inside, hp)

    def section_w(self, s):
        """The level `w` of a slice of its own chart."""
        d = self.chart0.d
        return (s-M(dot(d, self.o)))/M(dot(d, self.n))

    def weight(self, t0, t1):
        """S9d.3a's integral of `N(theta)` over `[t0, t1]`."""
        if self.orthonormal:
            return self.N0*(t1-t0)
        a, b, c = self.wall_vecs
        k = M(self.k)
        total = Z
        pieces = max(1, int(mp.ceil(abs(t1-t0)/(mp.pi/4))))
        step = (t1-t0)/pieces
        for i in range(pieces):
            lo = t0+i*step
            for xi, wi in s1.nodes(4):
                t = lo+step*(xi+1)/2
                v = tuple(mp.cos(t)*a[j]-mp.sin(t)*b[j]-k*c[j] for j in range(3))
                total += wi*step/2*mp.sqrt(dot(v, v))
        return total

    # ---- measures

    def closed(self):
        """Volume, moments and area (world) in closed form (S9d.3a's)."""
        r0, r1, h = M(self.r0), M(self.r1), M(self.h)
        V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
        Mw = mp.pi*h*h*(r0*r0+2*r0*r1+3*r1*r1)/12
        det = M(self.det)
        Vw = det*V
        nm = Mv(self.n)
        mom = tuple(self.om[i]*Vw+det*Mw*nm[i] for i in range(3))
        return Vw, mom, self.wall_area()+mp.pi*(r0*r0+r1*r1)*self.xy

    def wall_area(self):
        return M(self.h)*(M(self.r0)+M(self.r1))/2*self.weight(Z, 2*mp.pi)

    # ---- faces

    def faces(self):
        out = [Ruled(self, ('wall',), 'trig', trig_vec(add(self.o, scale(self.x, 0)), scale(self.x, self.r0),
                                                      scale(self.y, self.r0)),
                     trig_vec(self.n, scale(self.x, self.k), scale(self.y, self.k)), Z, 2*mp.pi,
                     (F(0), self.h), (self.r0, self.k), self.ruling_element, self.wall_area(), own=0)]
        for idx, w in ((1, F(0)), (2, self.h)):
            r = self.rim(w)
            if r == 0:
                continue
            P0 = self.world(F(0), F(0), w)
            span = (-r*(1+SKEW)-1, r*(1+SKEW)+1)
            out.append(PlanarFace(self, ('end', idx-1), idx, P0, self.x, add(scale(self.x, SKEW), self.y), span,
                                  mp.pi*M(r*r)*self.xy))
        return out

    def ruling_element(self, th):
        c, s_ = mp.cos(th), mp.sin(th)
        x, y, n = Mv(self.x), Mv(self.y), Mv(self.n)
        k = M(self.k)
        dU = tuple(-s_*x[i]+c*y[i] for i in range(3))
        B = tuple(n[i]+k*(c*x[i]+s_*y[i]) for i in range(3))
        v = cross(dU, B)
        return mp.sqrt(dot(v, v))


class BallIn(Input):
    """A sphere, cap or zone on S9d.1's model (S9d.2's `Ball` for slices)."""
    kind = 'sphere'

    def __init__(self, case, role):
        self.case, self.role = case, role
        self.ball = s2.Ball(case, role)
        S = self.ball.S
        self.S = S
        self.c, self.r2 = S.c, S.r2
        I3 = ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1)))
        self.surfaces = [Surf(I3, scale(self.c, -1), dot(self.c, self.c)-self.r2, (role, 'sphere'))]
        for a, b, end in S.planes:
            self.surfaces.append(Surf.halfspace(a, b, (role, 'zone', end)))
        o, x, y, n, _ = frame_rows(case)
        self.o, self.x, self.y, self.n = o, x, y, n
        self.cm = Mv(self.c)

    def contains(self, X, ignore=None):
        tests = [sum((X[i]-self.cm[i])**2 for i in range(3)) < M(self.r2)]
        for s in self.surfaces[1:]:
            tests.append(s.value(X) > 0)
        return all(t for i, t in enumerate(tests) if i != ignore)

    def zone_point(self, height):
        return tuple(self.o[i]+height*self.n[i] for i in range(3))

    def zone_heights(self):
        """[(surface index, height along n)]."""
        out = []
        idx = 1
        for h, end in zip(self.S.heights, ('low', 'high')):
            if h is not None:
                out.append((idx, h))
                idx += 1
        return out

    def zone_basis(self, idx):
        """Two vectors spanning a zone's plane: the frame's `x` and `y` in
        an exact frame, else (S9d.3c, a turned cap: its stored axes not
        exactly normal to its plane, whose normal is the stored axis) `y x a`
        and `a x (y x a)` for the plane's normal `a`, exactly in it."""
        a = self.surfaces[idx].a()
        if dot(a, self.x) == 0 and dot(a, self.y) == 0:
            return self.x, self.y
        e1 = cross(self.y, a)
        return e1, cross(a, e1)

    def edges(self):
        return [ConicEdge(self.zone_point(h), *self.zone_basis(idx), self.surfaces[0], self)
                for idx, h in self.zone_heights()]

    def vertices(self):
        return []

    def extent(self, d):
        dc = M(dot(d, self.c))
        rr = mp.sqrt(M(self.r2))*mp.sqrt(M(dot(d, d)))
        return dc-rr, dc+rr

    def setup(self, chart):
        self.ball.setup(chart)
        self.chart = chart
        self.conic = self.ball.conic

    def levels(self):
        return self.ball.levels()

    def section(self, s):
        return self.ball.section(s)

    def closed(self):
        return self.S.closed()

    def face_area(self):
        return self.ball.face_area()

    def faces(self):
        out = []
        R = mp.sqrt(M(self.r2))
        for idx, h in self.zone_heights():
            P0 = self.zone_point(h)
            w = self.surfaces[idx]
            dist = w.value(Mv(self.c))/mp.sqrt(M(dot(w.a(), w.a())))
            area = mp.pi*(R*R-dist*dist)
            rr = F(math.ceil(float(R)))+1
            e1, e2 = self.zone_basis(idx)
            out.append(PlanarFace(self, ('zone', idx), idx, P0, e1, add(scale(e1, SKEW), e2),
                                  (-rr*(1+SKEW), rr*(1+SKEW)), area))
        return out


class PrismIn(Input):
    """A prism of lines, arcs and circles on S9d.2's model."""
    kind = 'prism'

    def __init__(self, case, role):
        self.case, self.role = case, role
        P = s2.ArcPrism(case, role)
        self.P = P
        o, x, y, n = P.o, P.x, P.y, P.n
        self.o, self.x, self.y, self.n = o, x, y, n
        a1, a2, a3 = P.Kinv
        self.inv = (a1, a2, a3)
        self.surfaces = []
        self.el_surf = []
        for e in P.elements:
            k = e[-2]
            if e[0] == 'seg':
                p, q = e[1], e[2]
                ew = add(scale(x, q[0]-p[0]), scale(y, q[1]-p[1]))
                a = cross(ew, n)
                # The material is on the left of the element: inward is
                # n x ew, so a . X >= b with a = n x ew.
                a = scale(a, -1)
                self.surfaces.append(Surf.halfspace(a, dot(a, P.world(p, F(0))), (role, 'wall', k)))
            else:
                C, r = e[1], e[2]
                Mq = madd(outer(a1, a1), outer(a2, a2))
                Cw = add(scale(a1, C[0]), scale(a2, C[1]))
                m = sub(scale(mvec(Mq, o), -1), Cw)
                c = dot(o, mvec(Mq, o))+2*dot(Cw, o)+C[0]*C[0]+C[1]*C[1]-r*r
                self.surfaces.append(Surf(Mq, m, c, (role, 'wall', k)))
            self.el_surf.append(len(self.surfaces)-1)
        self.cap_idx = (len(self.surfaces), len(self.surfaces)+1)
        self.surfaces.append(Surf.halfspace(a3, dot(a3, o)+P.lo, (role, 'cap', 0)))
        self.surfaces.append(Surf.halfspace(scale(a3, -1), -(dot(a3, o)+P.hi), (role, 'cap', 1)))
        self.om = Mv(o)
        self.invm = tuple(Mv(r) for r in self.inv)

    def uvw(self, X):
        d = tuple(X[i]-self.om[i] for i in range(3))
        return tuple(sum(r[i]*d[i] for i in range(3)) for r in self.invm)

    def contains(self, X, ignore=None):
        u, v, w = self.uvw(X)
        if ignore != self.cap_idx[0] and not w > M(self.P.lo):
            return False
        if ignore != self.cap_idx[1] and not w < M(self.P.hi):
            return False
        return self.P.in_profile((u, v))

    def edges(self):
        P = self.P
        out = []
        for v in P.vertices:
            out.append(LineEdge(P.world(v, P.lo), sub(P.world(v, P.hi), P.world(v, P.lo))))
        for e in P.elements:
            for w in (P.lo, P.hi):
                if e[0] == 'seg':
                    p, q = e[1], e[2]
                    out.append(LineEdge(P.world(p, w), sub(P.world(q, w), P.world(p, w))))
                else:
                    _, C, r, t0, t1, full, k, _ = e
                    surf = self.surfaces[self.el_surf[k]]
                    out.append(ConicEdge(P.world(C, w), self.x, self.y, surf, self,
                                         arc=None if full else (C, t0, t1)))
        return out

    def vertices(self):
        return self.P.points3()

    def extent(self, d):
        x = self.P
        vals = [M(dot(d, p)) for p in x.points3()]
        for e in x.elements:
            if e[0] == 'arc':
                C, r = e[1], e[2]
                rad = M(r)*mp.sqrt(M(dot(d, x.x))**2+M(dot(d, x.y))**2)
                for w in (x.lo, x.hi):
                    c0 = M(dot(d, x.world(C, w)))
                    vals += [c0-rad, c0+rad]
        return min(vals), max(vals)

    def setup(self, chart):
        self.P.setup(chart)
        self.chart = chart
        self.conics = self.P.conics

    def levels(self):
        return self.P.cap_levels()

    def section(self, s):
        return self.P.section(s)

    def closed(self):
        V, mom = self.P.closed()
        return V, mom, self.area()

    def cap_area(self):
        A = self.P.profile_moments()[0]
        xy = cross(self.x, self.y)
        return A*mp.sqrt(M(dot(xy, xy)))

    def area(self):
        return 2*self.cap_area()+sum(f.area for f in self.faces() if f.tag[0] == 'wall')

    def faces(self):
        P = self.P
        out = []
        us, vs = [], []
        for e in P.elements:
            if e[0] == 'seg':
                us += [e[1][0], e[2][0]]
                vs += [e[1][1], e[2][1]]
            else:
                C, r = e[1], e[2]
                us += [C[0]-r, C[0]+r]
                vs += [C[1]-r, C[1]+r]
        span = (min(us)-SKEW*max(vs)-1, max(us)-SKEW*min(vs)+1)
        for cap, w in ((0, P.lo), (1, P.hi)):
            out.append(PlanarFace(self, ('cap', cap), self.cap_idx[cap], P.world((F(0), F(0)), w), self.x,
                                  add(scale(self.x, SKEW), self.y), span, self.cap_area()))
        h = P.hi-P.lo
        xm, ym, nm = Mv(self.x), Mv(self.y), Mv(self.n)
        for e in P.elements:
            k = e[-2]
            if e[0] == 'seg':
                p, q = e[1], e[2]
                ew = add(scale(self.x, q[0]-p[0]), scale(self.y, q[1]-p[1]))
                c = cross(ew, self.n)
                el = mp.sqrt(M(dot(c, c)))
                A0 = P.world(p, F(0))
                out.append(Ruled(self, ('wall', k), 'pol', tuple(Pol([A0[i], ew[i]]) for i in range(3)),
                                 tuple(Pol([v]) for v in self.n), Z, mp.mpf(1), (P.lo, P.hi), (F(1), F(0)),
                                 (lambda el_: lambda t: el_)(el), el*M(h), own=self.el_surf[k]))
            else:
                _, C, r, t0, t1, full, _, _ = e
                a0, a1_ = min(t0, t1), max(t0, t1)

                def el(th, r_=r):
                    t = tuple(-mp.sin(th)*xm[i]+mp.cos(th)*ym[i] for i in range(3))
                    c = cross(t, nm)
                    return M(r_)*mp.sqrt(dot(c, c))
                area = M(h)*mp.quad(el, [a0, a1_])
                out.append(Ruled(self, ('wall', k), 'trig', trig_vec(P.world(C, F(0)), scale(self.x, r),
                                                                    scale(self.y, r)),
                                 trig_vec(self.n, (F(0),)*3, (F(0),)*3), a0, a1_, (P.lo, P.hi), (F(1), F(0)), el,
                                 area, own=self.el_surf[k]))
        return out


def trig_vec(v0, v1, v2):
    return s2.trig_vec(v0, v1, v2)


def make_input(case, role):
    if case.cone is not None:
        return ConeIn(case, role)
    if case.sphere is not None:
        return BallIn(case, role)
    return PrismIn(case, role)


# ------------------------------------------------------------------ edges

class LineEdge:
    def __init__(self, P0, e):
        self.P0, self.e = P0, e

    def points(self, surf):
        """Its points on a surface (t in [0, 1] within 1e-20), mp."""
        c0, c1, c2 = surf.line(self.P0, self.e)
        if c0 == 0 and c1 == 0 and c2 == 0:
            return []
        out = []
        for t in roots2([c0, c1, c2]):
            if -mp.mpf(10)**-20 <= t <= 1+mp.mpf(10)**-20:
                out.append(tuple(M(self.P0[i])+t*M(self.e[i]) for i in range(3)))
        return out


class ConicEdge:
    """The trace of `surf` on the plane `P + alpha E1 + beta E2` (a circle
    of its input; `arc`: a prism's arc `(C, t0, t1)` in its `(u, v)`, the
    plane's origin its centre)."""

    def __init__(self, P, E1, E2, surf, owner, arc=None):
        self.P, self.E1, self.E2, self.surf, self.owner, self.arc = P, E1, E2, surf, owner, arc

    def points(self, other):
        A = tuple(Pol([self.P[i], self.E1[i]]) for i in range(3))
        B = tuple(Pol([v]) for v in self.E2)
        f = self.surf.family(A, B, 'pol')
        g = other.family(A, B, 'pol')
        if all(x.is_zero() for x in g):
            return []
        dg = 2 if not g[2].is_zero() else (1 if not g[1].is_zero() else 0)
        alphas = []
        if dg == 0:
            alphas = numeric_roots(g[0].poly())
        else:
            res = resultant(f, g)
            if res is None or res.is_zero():
                return []
            alphas = numeric_roots(res.poly())
        out = []
        Pm, E1m, E2m = Mv(self.P), Mv(self.E1), Mv(self.E2)
        for al in alphas:
            for be in solve_quadratic(*(x.value(al) for x in f)):
                X = tuple(Pm[i]+al*E1m[i]+be*E2m[i] for i in range(3))
                gv = other.value(X)
                scale_ = 1+sum(abs(v) for v in X)**2
                if abs(gv) > mp.mpf(10)**-20*scale_:
                    continue
                if self.arc is not None:
                    C, t0, t1 = self.arc
                    ang = mp.atan2(be, al)
                    near = min(abs((ang-t+mp.pi) % (2*mp.pi)-mp.pi) for t in (t0, t1))
                    if s2.span_param(t0, t1, ang) is None and near > mp.mpf(10)**-20:
                        continue
                out.append(X)
        return out


# ------------------------------------------------------------------ faces by lines

class Ruled:
    """A ruled face: lines `A(tau) + w B(tau)`, `tau` in `[t0, t1]`, `w` in
    `wr`, the element `element(tau) * (p0 + p1 w)`."""

    def __init__(self, owner, tag, ring, A, B, t0, t1, wr, weight, element, area, own):
        self.owner, self.tag, self.ring, self.A, self.B = owner, tag, ring, A, B
        self.t0, self.t1, self.wr, self.weight, self.element, self.area = t0, t1, wr, weight, element, area
        self.own = own
        self.planar = None

    def num(self, tau):
        if self.ring == 'trig':
            A = tuple(x.value(tau) for x in self.A)
            B = tuple(x.value(tau) for x in self.B)
        else:
            A = tuple(x.value(tau) for x in self.A)
            B = tuple(x.value(tau) for x in self.B)
        return A, B


class PlanarFace:
    """A planar face: the points of the plane `P0 + alpha E1 + beta E2`
    inside its solid with that plane left out, swept by the lines of
    constant `alpha` in `span`."""

    def __init__(self, owner, tag, idx, P0, E1, E2, span, area):
        self.owner, self.tag, self.idx = owner, tag, idx
        self.ring = 'pol'
        self.A = tuple(Pol([P0[i], E1[i]]) for i in range(3))
        self.B = tuple(Pol([v]) for v in E2)
        self.t0, self.t1 = M(span[0]), M(span[1])
        c = cross(E1, E2)
        self.el = mp.sqrt(M(dot(c, c)))
        self.element = lambda t: self.el
        self.weight = (F(1), F(0))
        self.wr = None
        self.area = area
        self.own = idx
        self.planar = owner.surfaces[idx]

    def num(self, tau):
        return tuple(x.value(tau) for x in self.A), tuple(x.value(tau) for x in self.B)


def sweep(face, other, size):
    """{in, out, same, opp, total} (world areas) of a face against the other
    input (see the module's note)."""
    owner = face.owner
    ring = face.ring
    fams = []
    for j, s in enumerate(owner.surfaces):
        if j == face.own:
            continue
        fams.append(('own', j, s, s.family(face.A, face.B, ring)))
    coplanar = None
    for j, s in enumerate(other.surfaces):
        cf = s.family(face.A, face.B, ring)
        if all(x.is_zero() for x in cf):
            if face.planar is not None and s.plane:
                same, orient = same_plane(face.planar, s)
                assert same
                # Outward normals: -a of each half-space.
                coplanar = (j, orient)
            continue
        fams.append(('other', j, s, cf))
    lo, hi = face.t0, face.t1
    events = [lo, hi]
    for _, _, _, cf in fams:
        c0, c1, c2 = cf
        if not c2.is_zero():
            events += ring_roots(c2, lo, hi, ring)
            events += ring_roots(c1*c1+neg((c0*c2).scale(F(4))), lo, hi, ring)
        elif not c1.is_zero():
            events += ring_roots(c1, lo, hi, ring)
        else:
            events += ring_roots(c0, lo, hi, ring)
    for (_, _, _, a), (_, _, _, b) in itertools.combinations(fams, 2):
        r = resultant(a, b)
        if r is not None:
            events += ring_roots(r, lo, hi, ring)
    if face.wr is not None:
        # The range's ends as planes: their resultants are each surface's
        # value at the end.
        for w in face.wr:
            for _, _, _, (c0, c1, c2) in fams:
                v = c0+c1.scale(F(w))+c2.scale(F(w)*F(w))
                events += ring_roots(v, lo, hi, ring)
    events = sorted(events)
    brk = []
    for p in events:
        if not brk or p-brk[-1] > mp.mpf(10)**-30*size:
            brk.append(p)
    num = [(kind, j, s) for kind, j, s, _ in fams]
    p0, p1 = M(face.weight[0]), M(face.weight[1])

    def wint(a, b):
        return p0*(b-a)+p1*(b*b-a*a)/2

    def f(tau):
        A, B = face.num(tau)
        ts = []
        for kind, j, s in num:
            ts += solve_quadratic(*s.num_line(A, B))
        if face.wr is not None:
            wl, wh = M(face.wr[0]), M(face.wr[1])
            ts = [t for t in ts if wl < t < wh]+[wl, wh]
        ts.sort()
        acc = {'in': Z, 'out': Z, 'same': Z, 'opp': Z}
        for a, b in zip(ts, ts[1:]):
            if b-a <= 0:
                continue
            m = (a+b)/2
            X = tuple(A[i]+m*B[i] for i in range(3))
            if face.wr is None and not owner.contains(X, ignore=face.own):
                continue
            if coplanar is not None:
                j, orient = coplanar
                if other.contains(X, ignore=j):
                    acc['same' if orient else 'opp'] += wint(a, b)
                else:
                    acc['out'] += wint(a, b)
            elif other.contains(X):
                acc['in'] += wint(a, b)
            else:
                acc['out'] += wint(a, b)
        el = face.element(tau)
        return [el*acc[k] for k in ('in', 'out', 'same', 'opp')]
    tol = mp.mpf(10)**-33*size**2
    tot = [Z]*4
    for a, b in zip(brk, brk[1:]):
        est, _ = integrate(f, a, b, tol)
        if est is not None:
            tot = [u+v for u, v in zip(tot, est)]
    return dict(zip(('in', 'out', 'same', 'opp'), tot)), len(brk)


# ------------------------------------------------------------------ charts

def chart_ok(chart, inputs):
    try:
        for x in inputs:
            x.setup(chart)
    except AssertionError:
        return False
    return True


def axis_of(x):
    if x.kind == 'cone':
        return cross(x.x, x.y)
    if x.kind == 'prism':
        return cross(x.x, x.y)
    return x.ball.S.m


def candidate_charts(inputs):
    """Charts in order of preference (see the module's note)."""
    out = []
    for x in inputs:
        if x.kind == 'cone':
            out.append(x.chart0)
    for x in inputs:
        if x.kind == 'prism':
            out.append(prism_chart(x.P))
        if x.kind == 'sphere' and x.ball.S.planes:
            out.append(Chart(x.ball.S.m))
    axes = [axis_of(x) for x in inputs]
    ks = (F(1, 4), F(1, 2), F(1), F(2), F(4))
    for a, b in itertools.permutations(axes, 2):
        for k in ks:
            for sg in (1, -1):
                d = add(a, scale(b, sg*k))
                if not is_zero(d):
                    out.append(Chart(d))
    generic = [(F(1), F(-2), F(0)), (F(2), F(1), F(-1)), (F(-1), F(3), F(2))]
    for a in axes:
        for g in generic:
            for k in (F(1, 8), F(1, 4), F(1, 2)):
                d = add(a, scale(g, k*max(1, round(math.sqrt(float(dot(a, a)))))))
                out.append(Chart(d))
    return out


def pick_chart(inputs, d=None, avoid=None):
    if d is not None:
        ch = Chart(d)
        assert chart_ok(ch, inputs), 'an unusable slicing direction'
        return ch
    for ch in candidate_charts(inputs):
        if avoid is not None and is_zero(cross(ch.d, avoid)):
            continue
        if chart_ok(ch, inputs):
            return ch
    return None


# ------------------------------------------------------------------ slicing breakpoints

def plane_quad_poly(plane, quad, d):
    """The plane's trace tangent to the quadric's section by `d . X = s`
    (S9d.3a's discriminant along the plane's line in the slice)."""
    a, b = plane.a(), -plane.c
    l = cross(d, a)
    if is_zero(l):
        return None
    X0 = s1.solve([d, a, l], [F(0), b, F(0)])
    X1 = s1.solve([d, a, l], [F(1), F(0), F(0)])
    Ml = mvec(quad.M, l)
    A = dot(l, Ml)
    B0 = 2*(dot(X0, Ml)+dot(quad.m, l))
    B1 = 2*dot(X1, Ml)
    C0 = quad.value(X0)
    MX1 = mvec(quad.M, X1)
    C1 = 2*(dot(X0, MX1)+dot(quad.m, X1))
    C2 = dot(X1, MX1)
    return [B0*B0-4*A*C0, 2*B0*B1-4*A*C1, B1*B1-4*A*C2]


def conics_of(x):
    if x.kind == 'prism':
        return list(x.conics)
    return [x.conic]


def slicing_breakpoints(inputs, chart):
    d = chart.d
    dm = Mv(d)
    pts = []
    for x in inputs:
        pts += [M(v) for v in x.levels()]
        for v in x.vertices():
            pts.append(M(dot(d, v)))
        if x.kind in ('cone', 'sphere'):
            _, rho = x.conic.centre_poly()
            pts += numeric_roots(rho)
    for x, y in itertools.permutations(inputs, 2):
        for e in x.edges():
            for s in y.surfaces:
                for X in e.points(s):
                    pts.append(dot(dm, X))
    planes = [s for x in inputs for s in x.surfaces if s.plane]
    quads = [s for x in inputs for s in x.surfaces if not s.plane]
    for p in planes:
        for q in quads:
            c = plane_quad_poly(p, q, d)
            if c is not None:
                pts += roots2(c)
    ca, cb = conics_of(inputs[0]), conics_of(inputs[1])
    for c1 in ca:
        for c2 in cb:
            t = tangency_poly(c1, c2)
            if pdeg(t) > 0:
                for f, _ in squarefree(t):
                    pts += numeric_roots(f)
    return pts


COMP_OPS = ('common', 'fuse', 'D-P', 'P-D')


def section_components(pieces, tol):
    """S9d.2's `components` ([(key, region)] of a section), two components
    allowed one key (both ends of a rod in one oblique slice)."""
    loops = s2.chain(pieces, tol)
    outer, holes = [], []
    for loop in loops:
        a = gsum(loop)[0]
        (outer if a > 0 else holes).append((loop, a))
    comps = [[loop, a, []] for loop, a in outer]
    for loop, _ in holes:
        X = s2.midpoint(loop[0])
        owners = [c for c in comps if parity(c[0], X)]
        assert owners, 'a hole outside every outer loop'
        min(owners, key=lambda c: c[1])[2].append(loop)
    out = []
    for loop, _, hs in comps:
        key = tuple(sorted([repr(s2.loop_key(loop))]+[repr(s2.loop_key(h)) for h in hs]))
        allp = list(loop)+[p for h in hs for p in h]
        region = s2.Pieces(allp, (lambda L, H: lambda X: parity(L, X) and not any(parity(h, X) for h in H))(loop, hs))
        out.append((key, region))
    return out


# ------------------------------------------------------------------ the pair

def same_model(a, b):
    if a.kind != b.kind or a.kind != 'cone':
        return False
    return (a.o, a.x, a.y, a.n, a.r0, a.r1, a.h) == (b.o, b.x, b.y, b.n, b.r0, b.r1, b.h)


class Pair:
    """Two inputs, the first the object (`A`), at least one a cone."""

    def __init__(self, obj, tool, d=None, avoid=None):
        self.A, self.B = make_input(obj, 'A'), make_input(tool, 'B')
        assert 'cone' in (self.A.kind, self.B.kind), 'S9d.3b: a cone against a curved input'
        self.inputs = (self.A, self.B)
        self.identical = same_model(self.A, self.B)
        self.size = self.case_size()
        self.chart = None if self.identical else pick_chart(self.inputs, d, avoid)
        assert self.identical or self.chart is not None, 'no slicing direction gives ellipses'
        self._res = self._faces = self._solids = None
        self.quad_error = Z
        self.sweep_breaks = 0

    def case_size(self):
        vals = [F(1)]
        for x in self.inputs:
            if x.kind == 'cone':
                vals += [abs(c)+max(x.r0, x.r1)+x.h for c in x.o]
            elif x.kind == 'sphere':
                vals += [abs(c)+x.S.r for c in x.c]
            else:
                vals += [abs(c) for p in x.P.points3() for c in p]
                for e in x.P.elements:
                    if e[0] == 'arc':
                        C, r = e[1], e[2]
                        for w in (x.P.lo, x.P.hi):
                            vals += [abs(c)+r*2 for c in x.P.world(C, w)]
        return M(max(vals))

    # ---- slicing

    def sections(self, s):
        return self.A.section(s), self.B.section(s)

    def integrand(self, s):
        D, P = self.sections(s)
        dp, pp = classify(D, P)
        out = []
        for g in (gsum([p for p, i in dp if i]), gsum([p for p, i in dp if not i]),
                  gsum([p for p, i in pp if i]), gsum([p for p, i in pp if not i])):
            out += [g[0], s*g[0], g[1], g[2]]
        extra = [Z]*4
        for k, (x, pieces) in enumerate(((self.A, dp), (self.B, pp))):
            for p, inside in pieces:
                if p.kind != 'arc':
                    continue
                if p.tag == (x.role, 'sphere'):
                    a = abs(p.angle())
                elif p.tag == (x.role, 'cone') and x.own:
                    t0, t1 = sorted((p.t0, p.t1))
                    a = x.weight(t0, t1)*(M(x.r0)+M(x.k)*x.section_w(s))
                else:
                    continue
                if inside:
                    extra[2*k] += a
                extra[2*k+1] += a
        return out+extra

    def breakpoints(self):
        pts = slicing_breakpoints(self.inputs, self.chart)
        lo = min(x.extent(self.chart.d)[0] for x in self.inputs)
        hi = max(x.extent(self.chart.d)[1] for x in self.inputs)
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
        """{name: (volume, moments)} for `A`, `B`, `common`, `fuse`, `A-B`,
        `B-A`; the spheres' faces and (in a cone's own chart) the cone's
        wall inside and outside the other."""
        if self._res is not None:
            return self._res
        if self.identical:
            V, mom, _ = self.A.closed()
            zero = (Z, (Z, Z, Z))
            self._res = {'A': (V, mom), 'B': (V, mom), 'common': (V, mom), 'fuse': (V, mom), 'A-B': zero,
                         'B-A': zero}
            return self._res
        v = self.measure()
        G = {name: v[4*k:4*k+4] for k, name in enumerate(('din', 'dout', 'pin', 'pout'))}
        combos = {'A': [('din', 1), ('dout', 1)], 'B': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'A-B': [('dout', 1), ('pin', -1)], 'B-A': [('pout', 1), ('din', -1)]}
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
        for k, x in enumerate(self.inputs):
            th_in, th_all = v[16+2*k], v[17+2*k]
            if x.kind == 'sphere':
                f = x.ball.rm/ch.dn
            elif x.kind == 'cone' and x.own:
                f = 1/M(dot(ch.d, x.n))
            else:
                continue
            out[('curved', x.role)] = {'in': f*th_in, 'out': f*(th_all-th_in), 'total': f*th_all}
        self._res = out
        return out

    # ---- faces

    def faces(self):
        """[(role, tag, classes, exact area)] of both inputs' faces."""
        if self._faces is not None:
            return self._faces
        out = []
        for x, y in ((self.A, self.B), (self.B, self.A)):
            if x.kind == 'sphere':
                if self.identical:
                    raise AssertionError
                sf = self.sliced()[('curved', x.role)]
                out.append((x.role, ('sphere',), {'in': sf['in'], 'out': sf['out'], 'same': Z, 'opp': Z},
                            x.face_area()))
            for face in x.faces():
                if self.identical:
                    cls = {'in': Z, 'out': Z, 'same': face.area, 'opp': Z}
                else:
                    cls, nb = sweep(face, y, self.size)
                    self.sweep_breaks += nb
                out.append((x.role, face.tag, cls, face.area))
        self._faces = out
        return out

    def classes(self, role):
        tot = {'in': Z, 'out': Z, 'same': Z, 'opp': Z}
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
        return {'cut': 'A-B'}.get(op, op)

    def volume(self, op):
        return self.sliced()[self.slice_op(op)]

    def solids(self, op):
        if self.identical:
            return 1
        return self.components()[self.slice_op(op)]

    def section_comps(self, s):
        """{op: [(key, region)]} of the result's section at `s`."""
        tol = mp.mpf(10)**-24*self.size
        D, P = self.sections(s)
        dp, pp = classify(D, P)
        out = {}
        for op in COMP_OPS:
            pieces = op_pieces(dp, pp, op)
            comps = section_components(pieces, tol) if pieces else []
            out[op] = [(k, r) for k, r in comps if area_of(r) > mp.mpf(10)**-30*self.size**2]
        return out

    def link(self, op, s0, c0, s1, c1, depth=0):
        """[(i, j)]: component `i` at `s0` is component `j` at `s1` (no
        breakpoint between): by their keys where each key is one
        component's, else among those of one key by their centroids (each
        nearest its own by less than half the next), halving the step
        until that holds."""
        k0, k1 = [k for k, _ in c0], [k for k, _ in c1]
        assert sorted(k0) == sorted(k1), ('components change inside an interval', op, mp.nstr(s0, 12),
                                          mp.nstr(s1, 12), k0, k1)
        if len(set(k0)) == len(k0):
            return [(i, k1.index(k)) for i, k in enumerate(k0)]

        def centroid(r):
            a, mx, my = gsum(r.boundary())
            return (mx/a, my/a)
        g0, g1 = [centroid(r) for _, r in c0], [centroid(r) for _, r in c1]
        pairs, ok = [], True
        for i, ka in enumerate(k0):
            cands = sorted((abs(g0[i][0]-g1[j][0])+abs(g0[i][1]-g1[j][1]), j) for j, kb in enumerate(k1) if kb == ka)
            if len(cands) > 1 and not cands[0][0] < cands[1][0]/2:
                ok = False
                break
            pairs.append((i, cands[0][1]))
        if ok and sorted(j for _, j in pairs) == list(range(len(c1))):
            return pairs
        assert depth < 24, 'components not followed through an interval'
        sm = (s0+s1)/2
        cm = self.section_comps(sm)[op]
        first = dict(self.link(op, s0, c0, sm, cm, depth+1))
        second = dict(self.link(op, sm, cm, s1, c1, depth+1))
        return [(i, second[m]) for i, m in first.items()]

    def components(self):
        """Solid counts of common, fuse, A-B and B-A (S9d.2's, components of
        one key followed by overlap)."""
        if self._solids is not None:
            return self._solids
        breaks = self.breaks if hasattr(self, 'breaks') else self.breakpoints()
        ends = []
        for a, b in zip(breaks, breaks[1:]):
            delta = min(mp.mpf(10)**-12*self.size, (b-a)/8)
            lo, hi = self.section_comps(a+delta), self.section_comps(b-delta)
            per = {}
            for op in COMP_OPS:
                per[op] = (lo[op], hi[op], self.link(op, a+delta, lo[op], b-delta, hi[op]))
            ends.append(per)
        counts = {}
        pinch = mp.mpf(10)**-20*self.size**2
        for op in COMP_OPS:
            parent = {}

            def find(x):
                while parent[x] != x:
                    x = parent[x]
                return x

            def union(x, y):
                parent[find(x)] = find(y)
            for i, e in enumerate(ends):
                lo, hi, pairs = e[op]
                for j in range(len(lo)):
                    parent[(i, 'lo', j)] = (i, 'lo', j)
                for j in range(len(hi)):
                    parent[(i, 'hi', j)] = (i, 'hi', j)
                for j0, j1 in pairs:
                    union((i, 'lo', j0), (i, 'hi', j1))
            for i in range(len(ends)-1):
                for ja, (_, ra) in enumerate(ends[i][op][1]):
                    for jb, (_, rb) in enumerate(ends[i+1][op][0]):
                        # Sections vanishing on both sides meet in a point
                        # or a curve there (apexes touching, a pinch).
                        if max(area_of(ra), area_of(rb)) < pinch:
                            continue
                        if joined(ra, rb):
                            union((i, 'hi', ja), (i+1, 'lo', jb))
            counts[op] = len({find(x) for x in parent})
        self._solids = {'common': counts['common'], 'fuse': counts['fuse'], 'A-B': counts['D-P'],
                        'B-A': counts['P-D']}
        return self._solids

    def result(self, op):
        vol, mom = self.volume(op)
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, Z, Z, None
        n = self.solids(op)
        assert n > 0, 'a result of positive volume without solids'
        return n, vol, self.area(op), tuple(m/vol for m in mom)


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
