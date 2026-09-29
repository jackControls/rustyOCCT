#!/usr/bin/env python3
"""Independent reference for S9d.3a of REVIEW_NOTES.md: Booleans of a cone
or frustum (`Solid::cone_with`: radius `bottom` at the frame's origin and
`top` at `height` along its normal, a zero radius an apex) against a
polyhedral prism (a profile of lines, every face a plane) in any relative
position, every pair of faces meeting in a line (two planes) or a conic (a
plane and the cone's wall: an ellipse, a parabola or a hyperbola, clipped
by the end planes' traces).

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal). The prism: S9d.1's (`sphere_boolean_reference.Prism`, a
nonconvex profile ear-clipped into convex pieces). The cone: the points `o
+ u x + v y + w n` of its stored frame with `0 <= w <= h` and `u^2 + v^2 <=
r(w)^2`, `r(w) = bottom + k w`, `k = (top - bottom) / h` (the decisions'
model, affine where the stored axes are not orthonormal). Everything is
done in the cone's chart `(u, v, w)`: the prism's vertices mapped there
exactly, the cone the quadric `q(p) = u^2 + v^2 - (bottom + k w)^2 <= 0`
between the end planes (within them `r >= 0`, so one nappe). Volumes and
moments are carried back by the chart's Jacobian `det(x, y, n)`, areas of
planes by Nanson's formula (`|a_u (y x n) + a_v (n x x) + a_w (x x y)| /
|a|` for a chart normal `a`), the wall by its own area element. Nothing
here uses the kernel, a surface/surface intersection or an arrangement.

* **Slicing (volume, first moments, the wall).** Both solids are sliced by
  the chart planes `d . p = s`, `d = (0, 0, 1)` (normal to the axis: each
  slice of the cone a disc of radius `r(s)`, the end planes slices) or any
  rational `d` tilted from the axis by less than the cone's complement
  (`k |d_uv| < |d_w|`: every slice an ellipse, cut by the end planes'
  lines), in coordinates `p = s e + alpha g + beta h` (`d . e = 1`, `g`, `h`
  spanning the slice and orthonormal for the quadric's form there, a
  Cholesky factor), where the cone's section is a disc `|(alpha, beta) -
  C(s)|^2 <= rho^2(s)` cut by the end planes' lines and each prism piece's a
  convex polygon (its edges' crossings). The regions are bounded by
  segments and arcs of one circle, and the four operations are taken apart
  as S9d.1's (`classify`: the disc's boundary cut at the pieces' lines and
  the pieces' edges at the circle, every piece classified at its midpoint;
  common, fuse and both cuts from the four classes, the common a second way
  by clipping each piece, its volume per piece kept for the solids; a line
  within 1e-30 of tangency to the circle taken as tangent, since a face
  tangent along a ruling is tangent to every slice's circle). The
  wall: on `(theta, w)` its area element is `r(w) N(theta) dtheta dw`, `N =
  |cos theta (y x n) - sin theta (x x n) - k (x x y)|` (`sqrt(1 + k^2)` in
  an orthonormal frame, else integrated along each arc by Gauss-Legendre),
  so its part inside the prism is the integral over `w` of `r(w)` times the
  weighted angle of the slice's circle inside the pieces (normal slices,
  where `theta` is the circle's own angle). Breakpoints, the roots of exact
  polynomials of degree at most two in `s`: a piece's vertex; an edge's
  meetings with the cone's quadric and with its end planes; a face's line
  tangent to the section (the discriminant of the quadric along the line
  of the face and the slice); the end circles' extreme levels; the apex; a
  face's and an end plane's line meeting the quadric. Each interval by
  Gauss-Legendre after `s = a + (b - a)(1 - cos t)/2`, doubling the rule
  until successive estimates agree within 1e-33 of the case's size to the
  fourth, else halving (`sphere_boolean_reference.integrate`).
* **Planar faces.** Each face of the prism (a wall's parallelogram, a
  cap's convex pieces) in the chart, by its own parameters: a face normal
  to the axis against the disc of its level in closed form (Green's
  theorem; inside by clipping, outside by the boundaries classified), or on
  an end plane with the same or the opposite orientation; any other face
  sliced by the lines of constant `w` across it, the chord's part inside
  the level's disc a quadratic's roots, integrated over `w` (the spacing
  `|a| / |a_uv|`) between the breakpoints (the face's vertices, the end
  planes, the chord tangent to the circle, the face's edges meeting the
  quadric). Each end disc of the cone against the prism's section by its
  plane (every piece crossing it cut into a polygon, a piece's face on it
  with its orientation) in closed form.
* **Solids** (combinatorics; the cone and every piece convex). Common: the
  pieces whose part of the common has positive volume, joined across a
  shared face whose part inside the cone has positive area. Fuse: one solid
  when the inputs overlap or share a face of opposite orientation, else two.
  `K - P` (a convex prism): the union of the open convex sets `int K n H_f`
  (`H_f` the open outer half-space of `P`'s face `f`), its components those
  of the graph joining two when `int K n H_f n H_g` is not empty (a volume
  sliced as above, positive). `P - K`: as S9d.1's, the components of `dP -
  K`, each face minus the cone's convex section one component per maximal
  run of its boundary outside it (edges cut where they enter the cone),
  pieces joined through their shared edges. Regions meeting along a curve
  or at a point are separate solids (the regularized Boolean).

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9d.2
references. A result keeps (object `A`, tool `B`): fuse `A`'s outside and
shared same-orientation pieces and `B`'s outside, cut `A`'s outside and
shared opposite pieces and `B`'s inside, common `A`'s inside and shared
same pieces and `B`'s inside. Nothing is shared with another reference but
`stored`, `stored_axes`, S9d.1's prism and its 2D regions and quadrature.
"""
from fractions import Fraction as F
import itertools

import mpmath as mp

from curve_surface_reference import stored_axes
from sphere_boolean_reference import (
    M, Mv, Disc2, Polys2, Prism, add, clip_poly, compose, cross, cross2, det3, dot, green_arc, green_seg, hull2,
    integrate, is_zero, line_of, nodes, roots2, scale, solve, sub, tau, vadd, zero3)

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
AXIS = (F(0), F(0), F(1))


# ------------------------------------------------------------------ 2D regions, tangencies snapped

# S9d.1's classification of a disc (cut by lines) against convex polygons,
# with a line within 1e-30 (relative) of tangency to the circle taken as
# tangent: a face tangent to the cone along a ruling is tangent to every
# slice's circle, where rounding would split a boundary at a point.
SNAP = mp.mpf(10)**-30


def circle_line_angles(C, rho, line):
    al, be, ga = line
    nn = mp.sqrt(al*al+be*be)
    if nn == 0 or rho == 0:
        return []
    val = (ga-al*C[0]-be*C[1])/(rho*nn)
    if abs(val) >= 1-SNAP:
        return []
    phi = mp.atan2(be, al)
    w = mp.acos(val)
    return [phi-w, phi+w]


def seg_circle_params(p, q, C, rho2):
    """Parameters `t` where `p + t (q - p)` crosses the circle."""
    e = (q[0]-p[0], q[1]-p[1])
    w = (p[0]-C[0], p[1]-C[1])
    A = e[0]*e[0]+e[1]*e[1]
    B = 2*(w[0]*e[0]+w[1]*e[1])
    Cc = w[0]*w[0]+w[1]*w[1]-rho2
    disc = B*B-4*A*Cc
    if A == 0 or disc <= SNAP*(B*B+abs(4*A*Cc)):
        return []
    sq = mp.sqrt(disc)
    return [(-B-sq)/(2*A), (-B+sq)/(2*A)]


def split_classify(curve, circle, other, other_circle):
    """Pieces of `curve` cut at `other`'s lines (and circle), each with
    whether its midpoint is inside `other`: [(green triple, inside, t0, t1)]."""
    out = []
    if curve[0] == 'seg':
        _, p, q = curve
        ts = [mp.mpf(0), mp.mpf(1)]
        for al, be, ga in other.split_lines():
            fp, fq = al*p[0]+be*p[1]-ga, al*q[0]+be*q[1]-ga
            if (fp < 0 < fq) or (fq < 0 < fp):
                ts.append(fp/(fp-fq))
        if other_circle is not None and other_circle[1] > 0:
            ts += [t for t in seg_circle_params(p, q, other_circle[0], other_circle[1]) if 0 < t < 1]
        ts.sort()
        for t0, t1 in zip(ts, ts[1:]):
            if t1-t0 <= 0:
                continue
            a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
            b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
            mid = ((a[0]+b[0])/2, (a[1]+b[1])/2)
            out.append((green_seg(a, b), other.contains(mid), t0, t1))
        return out
    _, a0, a1 = curve
    C, rho = circle
    T = tau()
    angs = [a0, a1]
    for ln in other.split_lines():
        for t in circle_line_angles(C, rho, ln):
            t = t+(mp.floor((a0-t)/T)+1)*T
            while t < a1:
                if t > a0:
                    angs.append(t)
                t += T
    angs.sort()
    for t0, t1 in zip(angs, angs[1:]):
        if t1-t0 <= 0:
            continue
        m = (t0+t1)/2
        mid = (C[0]+rho*mp.cos(m), C[1]+rho*mp.sin(m))
        out.append((green_arc(C, rho, t0, t1), other.contains(mid), t0, t1))
    return out


def classify(D, P):
    """Green's triples of `dD n P`, `dD - P`, `dP n D`, `dP - D`, and the
    arcs of `D`'s circle with whether each is inside `P`."""
    din, dout, pin, pout = zero3(), zero3(), zero3(), zero3()
    arcs = []
    circ = (D.C, D.rho)
    for curve in D.boundary():
        for g, inside, t0, t1 in split_classify(curve, circ, P, None):
            if inside:
                din = vadd(din, g)
            else:
                dout = vadd(dout, g)
            if curve[0] == 'arc':
                arcs.append((t0, t1, inside))
    dcirc = (D.C, D.rho2) if not D.empty() else None
    for curve in P.boundary():
        for g, inside, _, _ in split_classify(curve, None, D, dcirc):
            if inside and not D.empty():
                pin = vadd(pin, g)
            else:
                pout = vadd(pout, g)
    return din, dout, pin, pout, arcs


def convex_disc(poly, D):
    """Green's triple of a convex polygon (counter-clockwise) inside `D`: the
    polygon clipped by `D`'s lines, then its edges' chords inside the circle
    and the circle's arcs inside it."""
    if D.empty():
        return zero3()
    for ln in D.lines:
        poly = clip_poly(poly, ln)
        if not poly:
            return zero3()
    lines = [line_of(poly[i], poly[(i+1) % len(poly)]) for i in range(len(poly))]
    inside = lambda q: all(al*q[0]+be*q[1] > ga for al, be, ga in lines)
    total = zero3()
    angs = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i+1) % n]
        ts = seg_circle_params(p, q, D.C, D.rho2)
        if len(ts) == 2:
            t0, t1 = max(ts[0], mp.mpf(0)), min(ts[1], mp.mpf(1))
        else:
            m = ((p[0]+q[0])/2, (p[1]+q[1])/2)
            pin = (m[0]-D.C[0])**2+(m[1]-D.C[1])**2 < D.rho2
            t0, t1 = (mp.mpf(0), mp.mpf(1)) if pin else (mp.mpf(1), mp.mpf(0))
        if t1 > t0:
            a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
            b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
            total = vadd(total, green_seg(a, b))
        for t in ts:
            if 0 < t < 1:
                X = (p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1]))
                angs.append(mp.atan2(X[1]-D.C[1], X[0]-D.C[0]))
    T = tau()
    if not angs:
        if inside(D.point(mp.mpf(0))):
            total = vadd(total, green_arc(D.C, D.rho, mp.mpf(0), T))
        return total
    angs = sorted(a % T for a in angs)
    for i in range(len(angs)):
        a = angs[i]
        b = angs[i+1] if i+1 < len(angs) else angs[0]+T
        if b > a and inside(D.point((a+b)/2)):
            total = vadd(total, green_arc(D.C, D.rho, a, b))
    return total


def region_classes(D, polys):
    """(inside by clipping, outside by the boundaries classified) of a union
    of convex polygons against a disc."""
    inside = zero3()
    for p in polys:
        inside = vadd(inside, convex_disc(p, D))
    din, dout, pin, pout, _ = classify(D, Polys2(polys))
    return inside[0], vadd(pout, din, -1)[0]


# ------------------------------------------------------------------ the cone

class Cone:
    """A cone or frustum on its exact model, in its chart `(u, v, w)`."""

    def __init__(self, case):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        r0, r1, h = (F(v) for v in case.cone)
        assert h > 0 and r0 >= 0 and r1 >= 0 and r0 != r1, 'S9d.3a: a cone or a frustum'
        self.case = case
        self.o, self.x, self.y, self.n = o, x, y, n
        self.r0, self.r1, self.h = r0, r1, h
        self.k = (r1-r0)/h
        self.det = det3(x, y, n)
        assert self.det > 0
        self.inv = [scale(cross(y, n), 1/self.det), scale(cross(n, x), 1/self.det),
                    scale(cross(x, y), 1/self.det)]
        # The quadric q(p) = p^T Q p + 2 c . p + kappa.
        self.Q = (F(1), F(1), -self.k*self.k)
        self.c = (F(0), F(0), -r0*self.k)
        self.kappa = -r0*r0
        self.apex = (F(0), F(0), -r0/self.k)
        a, b, c = cross(y, n), cross(x, n), cross(x, y)
        self.wall = (Mv(a), Mv(b), Mv(c))
        self.orthonormal = (dot(a, b) == 0 and dot(a, c) == 0 and dot(b, c) == 0 and dot(a, a) == dot(b, b))
        self.N0 = mp.sqrt(M(dot(a, a)+self.k*self.k*dot(c, c)))
        self.xy = mp.sqrt(M(dot(c, c)))
        self.om = Mv(o)

    def chart(self, X):
        d = sub(X, self.o)
        return tuple(dot(row, d) for row in self.inv)

    def to_world_moments(self, V, m):
        """World volume and moments from chart ones (mpf)."""
        det = M(self.det)
        Vw = det*V
        axes = (Mv(self.x), Mv(self.y), Mv(self.n))
        return Vw, tuple(self.om[i]*Vw+det*sum(m[j]*axes[j][i] for j in range(3)) for i in range(3))

    def qf(self, a, b):
        return sum(Qi*ai*bi for Qi, ai, bi in zip(self.Q, a, b))

    def q(self, p):
        return self.qf(p, p)+2*dot(self.c, p)+self.kappa

    def r(self, w):
        return M(self.r0)+M(self.k)*w

    def line_quadratic(self, X0, l):
        """`q(X0 + t l)` as ascending coefficients in `t` (exact)."""
        return [self.q(X0), 2*(self.qf(X0, l)+dot(self.c, l)), self.qf(l, l)]

    def weight(self, t0, t1):
        """The integral of `N(theta)` over `[t0, t1]`."""
        if self.orthonormal:
            return self.N0*(t1-t0)
        a, b, c = self.wall
        k = M(self.k)
        total = mp.mpf(0)
        pieces = max(1, int(mp.ceil((t1-t0)/(mp.pi/4))))
        step = (t1-t0)/pieces
        for i in range(pieces):
            lo = t0+i*step
            for xi, wi in nodes(4):
                t = lo+step*(xi+1)/2
                v = tuple(mp.cos(t)*a[j]-mp.sin(t)*b[j]-k*c[j] for j in range(3))
                total += wi*step/2*mp.sqrt(dot(v, v))
        return total

    def closed(self):
        """Volume, first moments and area (world) in closed form."""
        r0, r1, h = M(self.r0), M(self.r1), M(self.h)
        V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
        Mw = mp.pi*h*h*(r0*r0+2*r0*r1+3*r1*r1)/12
        Vw, mom = self.to_world_moments(V, (0, 0, Mw))
        area = h*(r0+r1)/2*self.weight(mp.mpf(0), 2*mp.pi)+mp.pi*(r0*r0+r1*r1)*self.xy
        return Vw, mom, area

    def face_factor(self, a):
        """World area per chart area of a plane of chart normal `a`."""
        y, n, x = self.y, self.n, self.x
        v = add(add(scale(cross(y, n), a[0]), scale(cross(n, x), a[1])), scale(cross(x, y), a[2]))
        return mp.sqrt(M(dot(v, v)))/mp.sqrt(M(dot(a, a)))


# ------------------------------------------------------------------ slicing

def basis(K, d):
    """(g, h, e, jacobian) of the slices `d . p = s`: `d . e = 1`, `g`, `h`
    in the slice, orthonormal for the quadric's form there."""
    if d[0] == 0 and d[1] == 0:
        e = (F(0), F(0), 1/d[2])
        return (mp.mpf(1), mp.mpf(0), mp.mpf(0)), (mp.mpf(0), mp.mpf(1), mp.mpf(0)), Mv(e), abs(M(e[2]))
    g0 = cross(d, (F(0), F(0), F(1)))
    h0 = cross(d, g0)
    G11, G12, G22 = K.qf(g0, g0), K.qf(g0, h0), K.qf(h0, h0)
    assert G11 > 0 and G11*G22-G12*G12 > 0, 'slices tilted beyond an ellipse'
    l11 = mp.sqrt(M(G11))
    l21 = M(G12)/l11
    l22 = mp.sqrt(M(G22)-l21*l21)
    gm, hm = Mv(g0), Mv(h0)
    g = scale(gm, 1/l11)
    h = tuple(-l21/(l11*l22)*gm[i]+hm[i]/l22 for i in range(3))
    e = Mv(scale(d, 1/dot(d, d)))
    return g, h, e, abs(det3(g, h, e))


class Slicing:
    """Both solids sliced by the chart planes `d . p = s`."""

    def __init__(self, K, P, d, size, wall=False):
        self.K, self.P, self.d = K, P, tuple(F(v) for v in d)
        self.size, self.wall = size, wall
        self.pc = {key: K.chart(X) for key, X in P.P.items()}
        self.pm = {key: Mv(p) for key, p in self.pc.items()}
        self.level = {key: dot(self.d, p) for key, p in self.pc.items()}
        self.g, self.h, self.e, self.jac = basis(K, self.d)
        assert not wall or self.d == AXIS
        g, h, e = self.g, self.h, self.e
        dt = det3(g, h, e)
        self.dual = (scale(cross(h, e), 1/dt), scale(cross(e, g), 1/dt))
        Qm = tuple(M(v) for v in K.Q)
        cm = Mv(K.c)
        Qe = tuple(Qm[i]*e[i] for i in range(3))
        self.mg = (dot(g, Qe), dot(g, cm))
        self.mh = (dot(h, Qe), dot(h, cm))
        self.n2 = (dot(e, Qe), 2*dot(cm, e), M(K.kappa))
        self.strip = not (g[2] == 0 and h[2] == 0)
        self.quad_error = mp.mpf(0)
        self.breaks = None

    def to2(self, X):
        return (dot(self.dual[0], X), dot(self.dual[1], X))

    def disc(self, s):
        K = self.K
        if not self.strip:
            w = s*self.e[2]
            if w <= 0 or w >= M(K.h):
                return Disc2((mp.mpf(0), mp.mpf(0)), mp.mpf(-1))
        m = (s*self.mg[0]+self.mg[1], s*self.mh[0]+self.mh[1])
        rho2 = m[0]*m[0]+m[1]*m[1]-(self.n2[0]*s*s+self.n2[1]*s+self.n2[2])
        lines = []
        if self.strip:
            gw, hw, ew = self.g[2], self.h[2], self.e[2]
            lines = [(gw, hw, -s*ew), (-gw, -hw, s*ew-M(K.h))]
        return Disc2((-m[0], -m[1]), rho2, lines)

    def sections(self, s):
        out = []
        for piece in self.P.pieces:
            pts = []
            for i, j in piece.edges:
                li, lj = M(self.level[i]), M(self.level[j])
                if (li < s < lj) or (lj < s < li):
                    t = (s-li)/(lj-li)
                    X = add(self.pm[i], scale(sub(self.pm[j], self.pm[i]), t))
                    pts.append(self.to2(X))
            out.append(hull2(pts))
        return out

    def integrand(self, s):
        D = self.disc(s)
        polys = self.sections(s)
        P = Polys2([p for p in polys if p])
        din, dout, pin, pout, arcs = classify(D, P)
        out = []
        for g in (din, dout, pin, pout):
            out += [g[0], s*g[0], g[1], g[2]]
        for poly in polys:
            g = convex_disc(poly, D) if poly else zero3()
            out += [g[0], s*g[0], g[1], g[2]]
        if self.wall:
            r = self.K.r(s)
            th_in = th_all = mp.mpf(0)
            if not D.empty():
                for t0, t1, inside in arcs:
                    wgt = self.K.weight(t0, t1)
                    th_all += wgt
                    if inside:
                        th_in += wgt
            out += [r*th_in, r*th_all]
        return out

    def breakpoints(self):
        K, P, d = self.K, self.P, self.d
        levels = set(self.level.values())
        polys = []
        edges = set(e for piece in P.pieces for e in piece.edges)
        planes = []
        for piece in P.pieces:
            for f in piece.faces:
                p0, p1, p2 = (self.pc[k] for k in f['loop'][:3])
                a = cross(sub(p1, p0), sub(p2, p0))
                planes.append((a, dot(a, p0)))
        ends = [((F(0), F(0), F(1)), F(0)), ((F(0), F(0), F(1)), K.h)]
        levels.add(dot(d, K.apex))
        for (a, b) in ends:
            # The end circle's extreme levels.
            w0 = b
            r = K.r0+K.k*w0
            duv = d[0]*d[0]+d[1]*d[1]
            polys.append([(d[2]*w0)**2-r*r*duv, -2*d[2]*w0, F(1)])
        for i, j in edges:
            vi, vj = self.pc[i], self.pc[j]
            e = sub(vj, vi)
            de = dot(d, e)
            for a, b in ends:
                ae = dot(a, e)
                if ae != 0:
                    t = (b-dot(a, vi))/ae
                    if 0 <= t <= 1:
                        levels.add(dot(d, vi)+t*de)
            if de != 0:
                A0, A1, A2 = K.line_quadratic(vi, e)
                polys.append(compose(A2, A1, A0, -dot(d, vi)/de, 1/de))
        for a, b in planes+ends:
            l = cross(d, a)
            if is_zero(l):
                continue
            X0 = solve([d, a, l], [F(0), b, F(0)])
            X1 = solve([d, a, l], [F(1), F(0), F(0)])
            # q(X0 + s X1 + t l): its discriminant in t, a quadratic in s.
            A = K.qf(l, l)
            B0, B1 = 2*(K.qf(X0, l)+dot(K.c, l)), 2*K.qf(X1, l)
            C0 = K.q(X0)
            C1 = 2*(K.qf(X0, X1)+dot(K.c, X1))
            C2 = K.qf(X1, X1)
            polys.append([B0*B0-4*A*C0, 2*B0*B1-4*A*C1, B1*B1-4*A*C2])
        for (az, bz), (af, bf) in itertools.product(ends, planes):
            l = cross(az, af)
            if is_zero(l):
                continue
            X0 = solve([az, af, l], [bz, bf, F(0)])
            dl = dot(d, l)
            if dl == 0:
                levels.add(dot(d, X0))
                continue
            A0, A1, A2 = K.line_quadratic(X0, l)
            polys.append(compose(A2, A1, A0, -dot(d, X0)/dl, 1/dl))
        pts = [M(v) for v in levels]
        for p in polys:
            pts += roots2(p)
        lo_s = min(M(v) for v in self.level.values())
        hi_s = max(M(v) for v in self.level.values())
        cone_levels = [M(dot(d, (F(0), F(0), w))) for w in (F(0), K.h)]
        for w in (F(0), K.h):
            r = M(K.r0+K.k*w)
            cone_levels += [M(d[2]*w)+sg*r*mp.sqrt(M(d[0]*d[0]+d[1]*d[1])) for sg in (-1, 1)]
        a, b = min([lo_s]+cone_levels), max([hi_s]+cone_levels)
        pts = sorted(p for p in pts if a <= p <= b)
        out = []
        gap = mp.mpf(10)**-30*self.size
        for p in [a]+pts+[b]:
            if not out or p-out[-1] > gap:
                out.append(p)
        self.raw_breaks = pts
        return out

    def measure(self):
        if self.breaks is None:
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

    def results(self):
        """{name: (world volume, world moments)} for `K`, `P`, `common`,
        `fuse`, `K-P`, `P-K`, `common2` (by clipping) and each piece's part
        of the common (`pieces`); the wall's area inside the prism and in
        all (`wall`, normal slices only)."""
        v = self.measure()
        G = {}
        for k, name in enumerate(('din', 'dout', 'pin', 'pout')):
            G[name] = v[4*k:4*k+4]
        npieces = len(self.P.pieces)
        per = [v[16+4*k:20+4*k] for k in range(npieces)]
        G['conv'] = [sum(p[i] for p in per) for i in range(4)]
        combos = {'K': [('din', 1), ('dout', 1)], 'P': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'K-P': [('dout', 1), ('pin', -1)], 'P-K': [('pout', 1), ('din', -1)],
                  'common2': [('conv', 1)]}

        def world(A, sA, ma, mb):
            V = self.jac*A
            mom = tuple(self.jac*(self.e[i]*sA+self.g[i]*ma+self.h[i]*mb) for i in range(3))
            return self.K.to_world_moments(V, mom)
        out = {}
        for name, parts in combos.items():
            vals = [sum(k*G[p][i] for p, k in parts) for i in range(4)]
            out[name] = world(*vals)
        out['pieces'] = [world(*p)[0] for p in per]
        if self.wall:
            th_in, th_all = v[16+4*npieces], v[17+4*npieces]
            out['wall'] = {'in': th_in, 'out': th_all-th_in, 'total': th_all}
        return out


# ------------------------------------------------------------------ planar faces

def ccw2(poly):
    A = sum(cross2(poly[i-1], poly[i]) for i in range(len(poly)))
    return poly if A > 0 else list(reversed(poly))


def slanted_face(K, polys, a, b, size):
    """(inside, all) chart areas of convex chart polygons in the plane `a . p
    = b` (not normal to the axis), sliced by the lines of constant `w`."""
    au, av, aw = a
    nuv2 = au*au+av*av
    spacing = mp.sqrt(M(dot(a, a)))/mp.sqrt(M(nuv2))
    levels = set()
    quads = [[b*b-nuv2*K.r0*K.r0, -2*b*aw-2*nuv2*K.r0*K.k, aw*aw-nuv2*K.k*K.k]]
    edges = []
    for poly in polys:
        for i in range(len(poly)):
            p, q = poly[i], poly[(i+1) % len(poly)]
            levels.add(p[2])
            edges.append((p, q))
            e = sub(q, p)
            if e[2] != 0:
                A0, A1, A2 = K.line_quadratic(p, e)
                quads.append(compose(A2, A1, A0, -p[2]/e[2], 1/e[2]))
    lo, hi = min(levels), max(levels)
    pts = [M(v) for v in levels | {F(0), K.h} if lo <= v <= hi]
    for c in quads:
        pts += [t for t in roots2(c) if M(lo) <= t <= M(hi)]
    pts.sort()
    brk = []
    for p in pts:
        if not brk or p-brk[-1] > mp.mpf(10)**-30*size:
            brk.append(p)
    pm = [[Mv(p) for p in poly] for poly in polys]
    hm = M(K.h)

    def f(s):
        Lin = Lall = mp.mpf(0)
        for poly in pm:
            ends = []
            for i in range(len(poly)):
                p, q = poly[i], poly[(i+1) % len(poly)]
                if (p[2] < s < q[2]) or (q[2] < s < p[2]):
                    t = (s-p[2])/(q[2]-p[2])
                    ends.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
            if len(ends) != 2:
                continue
            P1, P2 = ends
            L = mp.sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2)
            Lall += L
            if not 0 < s < hm:
                continue
            r2 = K.r(s)**2
            ts = seg_circle_params(P1, P2, (mp.mpf(0), mp.mpf(0)), r2)
            if len(ts) == 2:
                Lin += L*max(mp.mpf(0), min(ts[1], mp.mpf(1))-max(ts[0], mp.mpf(0)))
        return [Lin, Lall]
    tol = mp.mpf(10)**-33*size**3
    tin = tall = mp.mpf(0)
    for x0, x1 in zip(brk, brk[1:]):
        est, _ = integrate(f, x0, x1, tol)
        if est is not None:
            tin += est[0]
            tall += est[1]
    return spacing*tin, spacing*tall


def face_classes(K, polys, outward, size):
    """{in, out, same, opp} (world areas) of a planar face of the prism
    (convex chart polygons, its chart outward normal) against the cone."""
    p0, p1, p2 = polys[0][0], polys[0][1], polys[0][2]
    a = cross(sub(p1, p0), sub(p2, p0))
    if dot(a, outward) < 0:
        a = scale(a, -1)
    b = dot(a, p0)
    factor = K.face_factor(a)
    cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
    if a[0] == 0 and a[1] == 0:
        w0 = b/a[2]
        p2d = [ccw2([(M(p[0]), M(p[1])) for p in poly]) for poly in polys]
        total = sum(abs(sum(cross2(p[i-1], p[i]) for i in range(len(p)))/2) for p in p2d)
        r = K.r0+K.k*w0
        if w0 < 0 or w0 > K.h or r == 0:
            cls['out'] = total*factor
            return cls
        D = Disc2((mp.mpf(0), mp.mpf(0)), M(r*r))
        inside, outside = region_classes(D, p2d)
        if 0 < w0 < K.h:
            cls['in'] = inside*factor
        else:
            # On an end plane: the end's outward normal is -w at 0, +w at h.
            end_out = -1 if w0 == 0 else 1
            cls['same' if (a[2] > 0) == (end_out > 0) else 'opp'] = inside*factor
        cls['out'] = outside*factor
        return cls
    inside, total = slanted_face(K, polys, a, b, size)
    cls['in'] = inside*factor
    cls['out'] = (total-inside)*factor
    return cls


def prism_face_classes(K, prism, size):
    """[(tag, {in, out, same, opp}, exact area)] for the prism's faces."""
    out = []
    for tag, polys, piece, area in prism.faces():
        cp = [[K.chart(X) for X in poly] for poly in polys]
        p0, p1, p2 = cp[0][0], cp[0][1], cp[0][2]
        a = cross(sub(p1, p0), sub(p2, p0))
        if dot(a, K.chart(piece.centroid)) > dot(a, p0):
            a = scale(a, -1)
        out.append((tag, face_classes(K, cp, a, size), area))
    return out


def end_disc_classes(K, prism):
    """[(end, {in, out, same, opp}, exact area)] for the cone's end discs
    (none at an apex)."""
    out = []
    for w0, sign, end in ((F(0), 1, 'bottom'), (K.h, -1, 'top')):
        r = K.r0+K.k*w0
        if r == 0:
            continue
        D = Disc2((mp.mpf(0), mp.mpf(0)), M(r*r))
        area = mp.pi*M(r*r)*K.xy
        groups = {'in': [], 'same': [], 'opp': []}
        for piece in prism.pieces:
            ch = {k: K.chart(prism.P[k]) for k in piece.keys}
            vals = {k: sign*(p[2]-w0) for k, p in ch.items()}
            pos = any(v > 0 for v in vals.values())
            neg = any(v < 0 for v in vals.values())
            zeros = [k for k, v in vals.items() if v == 0]
            if pos and neg:
                pts = [(M(ch[k][0]), M(ch[k][1])) for k in zeros]
                for i, j in piece.edges:
                    if (vals[i] < 0 < vals[j]) or (vals[j] < 0 < vals[i]):
                        t = vals[i]/(vals[i]-vals[j])
                        X = add(ch[i], scale(sub(ch[j], ch[i]), t))
                        pts.append((M(X[0]), M(X[1])))
                poly = hull2(pts)
                if poly:
                    groups['in'].append(poly)
            elif len(zeros) >= 3:
                poly = hull2([(M(ch[k][0]), M(ch[k][1])) for k in zeros])
                # A piece on the cone's side has its face's outward normal
                # the disc's; on the other side, opposite.
                groups['same' if pos else 'opp'].append(poly)
        cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        allp = groups['in']+groups['same']+groups['opp']
        for key, polys in groups.items():
            for p in polys:
                cls[key] += convex_disc(p, D)[0]*K.xy
        if allp:
            din, dout, pin, pout, _ = classify(D, Polys2(allp))
            cls['out'] = (dout[0]-pin[0])*K.xy
        else:
            cls['out'] = area
        out.append((end, cls, area))
    return out


# ------------------------------------------------------------------ solids

def region_area(D):
    total = mp.mpf(0)
    for curve in D.boundary():
        if curve[0] == 'arc':
            total += green_arc(D.C, D.rho, curve[1], curve[2])[0]
        else:
            total += green_seg(curve[1], curve[2])[0]
    return total


def cut_volume(K, halfspaces, size):
    """The chart volume of the open cone cut by open half-spaces `a . p > b`
    (chart), by normal slices."""
    lines, bounds = [], []
    levels = {F(0), K.h}
    polys = []
    for a, b in halfspaces:
        if a[0] == 0 and a[1] == 0:
            bounds.append((a[2], b))
            levels.add(b/a[2])
            continue
        lines.append((a, b))
        nuv2 = a[0]*a[0]+a[1]*a[1]
        polys.append([b*b-nuv2*K.r0*K.r0, -2*b*a[2]-2*nuv2*K.r0*K.k, a[2]*a[2]-nuv2*K.k*K.k])
    for (a1, b1), (a2, b2) in itertools.combinations(halfspaces, 2):
        l = cross(a1, a2)
        if is_zero(l):
            continue
        X0 = solve([a1, a2, l], [b1, b2, F(0)])
        if l[2] == 0:
            levels.add(X0[2])
            continue
        A0, A1, A2 = K.line_quadratic(X0, l)
        polys.append(compose(A2, A1, A0, -X0[2]/l[2], 1/l[2]))
    pts = sorted(set([M(v) for v in levels if 0 <= v <= K.h]
                     + [t for c in polys for t in roots2(c) if 0 <= t <= M(K.h)]))
    hm = M(K.h)

    def f(s):
        if any(M(aw)*s <= M(b) for aw, b in bounds) or not 0 < s < hm:
            return [mp.mpf(0)]
        D = Disc2((mp.mpf(0), mp.mpf(0)), K.r(s)**2,
                  [(M(a[0]), M(a[1]), M(b)-M(a[2])*s) for a, b in lines])
        return [region_area(D)]
    tol = mp.mpf(10)**-30*size**3
    total = mp.mpf(0)
    for x0, x1 in zip(pts, pts[1:]):
        est, _ = integrate(f, x0, x1, tol)
        if est is not None:
            total += est[0]
    return total


class Solids:
    """Solid counts by convexity (see the module's note)."""

    def __init__(self, pair):
        self.pair, self.K, self.P = pair, pair.K, pair.P
        self.eps = mp.mpf(10)**-25*pair.size**3
        self.aeps = mp.mpf(10)**-25*pair.size**2

    def common(self):
        vols = self.pair.sliced()['pieces']
        pieces = [k for k, v in enumerate(vols) if v > self.eps]
        parent = {k: k for k in pieces}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k
        for i, j in itertools.combinations(pieces, 2):
            pi, pj = self.P.pieces[i], self.P.pieces[j]
            for f in pi.faces:
                if f['tag'][0] != 'internal' or not any(g['tag'] == f['tag'] for g in pj.faces):
                    continue
                poly = [self.K.chart(self.P.P[k]) for k in f['loop']]
                # The world covector -a (outward) in the chart.
                outward = tuple(-dot(f['a'], col) for col in (self.K.x, self.K.y, self.K.n))
                if face_classes(self.K, [poly], outward, self.pair.size)['in'] > self.aeps:
                    parent[find(i)] = find(j)
        return len({find(k) for k in pieces})

    def fuse(self, opp_area):
        if self.pair.sliced()['common'][0] > self.eps or opp_area > self.aeps:
            return 1
        return 2

    def cone_minus(self):
        """Components of `K - P`, `P` convex."""
        assert len(self.P.pieces) == 1, 'S9d.3a reference: K - P for a convex prism only'
        piece = self.P.pieces[0]
        if self.pair.sliced()['common'][0] <= self.eps:
            return 1
        outer = []
        for f in piece.faces:
            # The chart form of the world half-space a . X < b.
            a, b = f['a'], f['b']
            ac = tuple(-sum(a[i]*col[i] for i in range(3)) for col in (self.K.x, self.K.y, self.K.n))
            outer.append((ac, -(b-dot(a, self.K.o))))
        size = self.pair.size
        nodes_ = [k for k, hs in enumerate(outer) if cut_volume(self.K, [hs], size) > self.eps]
        parent = {k: k for k in nodes_}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k
        for i, j in itertools.combinations(nodes_, 2):
            if cut_volume(self.K, [outer[i], outer[j]], size) > self.eps:
                parent[find(i)] = find(j)
        return len({find(k) for k in nodes_})

    def edge_outside(self, i, j):
        """The parts of the edge from key `i` to key `j` outside the cone, as
        S9d.1's."""
        K = self.K
        vi, vj = K.chart(self.P.P[i]), K.chart(self.P.P[j])
        e = sub(vj, vi)
        ts = [mp.mpf(0), mp.mpf(1)]
        ts += [t for t in roots2(K.line_quadratic(vi, e)) if 0 < t < 1]
        if e[2] != 0:
            for w0 in (F(0), K.h):
                t = (w0-vi[2])/e[2]
                if 0 < t < 1:
                    ts.append(M(t))
        ts.sort()
        vm, em = Mv(vi), Mv(e)
        inside = []
        for t0, t1 in zip(ts, ts[1:]):
            if t1-t0 <= mp.mpf(10)**-30:
                continue
            t = (t0+t1)/2
            p = tuple(vm[k]+t*em[k] for k in range(3))
            if 0 < p[2] < M(K.h) and p[0]**2+p[1]**2 < K.r(p[2])**2:
                inside.append((t0, t1))
        if not inside:
            return ('whole', i, j), ('whole', i, j)
        t0, t1 = inside[0][0], inside[-1][1]
        return (('start', i, j) if t0 > 0 else None), (('end', i, j) if t1 < 1 else None)

    def prism_minus(self):
        """Components of `P - K`."""
        ends = {}
        for piece in self.P.pieces:
            for i, j in piece.edges:
                if (i, j) not in ends:
                    ends[(i, j)] = self.edge_outside(i, j)
        nodes_ = set(x for v in ends.values() for x in v if x is not None)
        parent = {k: k for k in nodes_}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k

        def at(edge_from, edge_to, key):
            e = tuple(sorted((edge_from, edge_to)))
            pa, pb = ends[e]
            return pa if key == e[0] else pb
        for piece in self.P.pieces:
            for f in piece.faces:
                L = f['loop']
                n = len(L)
                for t in range(n):
                    prev, cur, nxt = L[t-1], L[t], L[(t+1) % n]
                    a, b = at(prev, cur, cur), at(cur, nxt, cur)
                    if a is not None and b is not None:
                        parent[find(a)] = find(b)
        return len({find(k) for k in nodes_})


# ------------------------------------------------------------------ the pair

def case_size(K, P):
    vals = [F(1)]
    vals += [abs(x) for p in P.P.values() for x in p]
    vals += [abs(x)+max(K.r0, K.r1)+K.h for x in K.o]
    return M(max(vals))


class Pair:
    """A cone and a prism, in either order (`obj`, `tool` Boolean cases)."""

    def __init__(self, obj, tool):
        self.cone_first = obj.cone is not None
        k_case, p_case = (obj, tool) if self.cone_first else (tool, obj)
        assert k_case.cone is not None and p_case.cone is None and p_case.sphere is None, \
            'S9d.3a: a cone and a prism'
        self.K, self.P = Cone(k_case), Prism(p_case)
        self.size = case_size(self.K, self.P)
        self.slicing = Slicing(self.K, self.P, AXIS, self.size, wall=True)
        self._res = self._faces = None

    def sliced(self):
        if self._res is None:
            self._res = self.slicing.results()
        return self._res

    def faces(self):
        """[(input 'K' or 'P', tag, classes, exact area)]."""
        if self._faces is None:
            out = [('P', tag, cls, area) for tag, cls, area in prism_face_classes(self.K, self.P, self.size)]
            out += [('K', ('end', end), cls, area) for end, cls, area in end_disc_classes(self.K, self.P)]
            wall = self.sliced()['wall']
            K = self.K
            total = M(K.h)*(M(K.r0)+M(K.r1))/2*K.weight(mp.mpf(0), 2*mp.pi)
            out.append(('K', ('wall',), {'in': wall['in'], 'out': wall['out'], 'same': mp.mpf(0),
                                         'opp': mp.mpf(0)}, total))
            self._faces = out
        return self._faces

    def classes(self, which):
        tot = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        for w, _, cls, _ in self.faces():
            if w == which:
                for k in tot:
                    tot[k] += cls[k]
        return tot

    def roles(self):
        return ('K', 'P') if self.cone_first else ('P', 'K')

    def area(self, op):
        a, b = self.roles()
        ca, cb = self.classes(a), self.classes(b)
        keep = {'fuse': (('out', 'same'), ('out',)), 'cut': (('out', 'opp'), ('in',)),
                'common': (('in', 'same'), ('in',))}[op]
        return sum(ca[k] for k in keep[0])+sum(cb[k] for k in keep[1])

    def volume(self, op):
        r = self.sliced()
        key = {'fuse': 'fuse', 'common': 'common', 'cut': 'K-P' if self.cone_first else 'P-K'}[op]
        return r[key]

    def solids(self, op):
        sol = Solids(self)
        if op == 'common':
            return sol.common()
        if op == 'fuse':
            return sol.fuse(self.classes('K')['opp'])
        return sol.cone_minus() if self.cone_first else sol.prism_minus()

    def result(self, op):
        vol, mom = self.volume(op)
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        n = self.solids(op)
        assert n > 0, 'a result of positive volume without solids'
        return n, vol, self.area(op), tuple(m/vol for m in mom)


def number(x):
    x = M(x)
    return '0.0' if abs(x) < M(10)**-30 else mp.nstr(x, 25, min_fixed=-5, max_fixed=5)


def rows(obj, operation, tool, pair=None):
    """`result N volume area cx cy cz` or `empty`, and the pair (reused
    across the three operations)."""
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair
