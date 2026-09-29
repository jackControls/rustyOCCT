#!/usr/bin/env python3
"""Independent reference for S9d.4a of REVIEW_NOTES.md: Booleans of a whole
torus (`Solid::torus_with` with the full tube and turn) against a polyhedral
prism (a profile of lines, every face a plane) in any relative position,
every pair of faces meeting in a line (two planes) or a spiric section (a
plane and the torus's wall).

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal). The prism: S9d.1's (`sphere_boolean_reference.Prism`, a
nonconvex profile ear-clipped into convex pieces). The torus: the points
`o + u x + v y + w n` of its stored frame with `(|p|^2 + R^2 - r^2)^2 <= 4
R^2 (u^2 + v^2)`, `p = (u, v, w)` (the decisions' model, affine where the
stored axes are not orthonormal), `R > r > 0`. Everything is done in the
torus's chart `(u, v, w)`: the prism's vertices mapped there exactly;
volumes and moments carried back by the chart's Jacobian `det(x, y, n)`,
areas of planes by Nanson's formula (as S9d.3a's), the wall by its own area
element: at `p(theta, phi) = ((R + r cos phi) cos theta, (R + r cos phi) sin
theta, r sin phi)` it is `r (R + r cos phi) |N_u (y x n) + N_v (n x x) + N_w
(x x y)| dtheta dphi`, `N = (cos phi cos theta, cos phi sin theta, sin phi)`
(`r (R + r cos phi)` in an orthonormal frame; otherwise integrated by
Gauss-Legendre along each arc, the factor within 1e-15 of one). Nothing
here uses the kernel, a surface/surface intersection or an arrangement.

* **Normal slices (volume, first moments; the wall a second way).** Both
  solids are sliced by the chart planes `w = s`: the torus's slice is the
  annulus between the circles of radii `R - q` and `R + q` about the axis,
  `q = sqrt(r^2 - s^2)`, and each prism piece's a convex polygon. S9d.1's
  classification (`cone_boolean_reference.classify`: boundaries cut and
  classified, Green's theorem over segments and arcs) is run twice, the
  outer disc and the inner one against the pieces, and the operations are
  their combinations: common `(o_in + p_in(o)) - (i_in + p_in(i))`, fuse
  `o_out + p_out(o) - i_out + p_in(i)`, `K - P` `(o_out - p_in(o)) - (i_out
  - p_in(i))`, `P - K` `(p_out(o) - o_in) + (i_in + p_in(i))`; the common a
  second way by clipping each piece by both discs. The wall: at level `s`
  its element is `r rho / q` times the circle's angle (`rho` either
  radius), so its part inside the prism is the integral of the weighted
  angles of both circles inside the pieces. Breakpoints: the pieces'
  vertices' levels; `s = +-r`, where the circles meet; an edge's meetings
  with the torus (the quartic `F(X0 + t e)`); a face's line tangent to a
  circle (`4 R^2 N D^2 = (D^2 + N (R^2 - r^2 + s^2))^2`, `D = b - a_w s`,
  `N = a_u^2 + a_v^2`, a quartic in `s`). Real roots of the quartics by
  mpmath's `polyroots` at 40 digits, a near-real pair taken as a root (a
  breakpoint more is harmless).
* **Meridian half-planes (the wall in its own parameters; volume and
  moments a second way; solids).** Both solids are sliced by the half-planes
  `theta = const` about the axis, in coordinates `(t, w)` (`t >= 0` the
  distance from the axis in the chart): the torus's section is the disc of
  radius `r` about `(R, 0)` for every `theta`, the prism's a convex polygon
  (a piece's section by the plane, clipped at the axis). The volume element
  is `t dt dw dtheta`, so volume and moments are Green's integrals of `t`,
  `t^2` and `t w` over the classified boundaries (closed forms: arcs by
  their antiderivatives in `phi`, segments exactly), and the wall's part
  inside is the integral over `theta` of `r (R + r cos phi)` over the arcs
  of the disc's circle inside the section (`phi` its angle about `(R, 0)`,
  the torus's own). Breakpoints in `theta`: the prism's vertices' angles
  (and opposite); an edge's and any two faces' line's meetings with the
  torus; a face's line tangent to the tube's circle (`(alpha R - b)^2 = r^2
  (alpha^2 + a_w^2)`, `alpha = a_u cos theta + a_v sin theta`: a quadratic
  in `alpha`, the decisions' spiric form) and parallel to the axis (`alpha
  = 0`).
* **Planar faces.** Each face of the prism, by its own parameters as
  S9d.3a's: a face normal to the axis against the annulus of its level in
  closed form (both discs classified), any other face sliced by the lines of
  constant `w` across it, the chord's part inside the annulus (inside the
  outer circle less inside the inner) integrated over `w` between the
  breakpoints (the face's vertices, `w = +-r`, the chord tangent to a
  circle, the face's edges meeting the torus).
* **Solids** (a convex prism: one piece). In each meridian half-plane the
  sections are `D` (the disc) and `C` (a convex polygon). The common's
  section `D n C` is convex; `D - C` is the union of the convex sets `D n
  H_f` (`H_f` a face's open outer half-plane), its components those of the
  graph joining two where `D n H_f n H_g` has area; `C - D` has one
  component per maximal run of `C`'s boundary outside `D` (one if none of
  it is inside), each holding at least one of `C`'s vertices (an edge's
  crossing, or a face's meeting with the axis). Between breakpoints this
  structure is fixed; components are followed around the axis and joined
  across a breakpoint: the common's where its section has area there, `D -
  C`'s where one face's `D n H_f` does, `C - D`'s where a vertex of one
  side and one of the other meet there outside `D` (a vertex of the prism
  on the plane there: every crossing near it). Sections vanishing at a
  breakpoint are not joined (regions meeting along a curve or at a point
  are separate solids, the regularized Boolean). Fuse: one solid when the
  inputs overlap, else two (the torus has no planar face to share).

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9d.3
references. A result keeps (object `A`, tool `B`): fuse `A`'s outside and
`B`'s outside, cut `A`'s outside and `B`'s inside, common `A`'s inside and
`B`'s inside. Nothing is shared with another reference but `stored`,
`stored_axes`, S9d.1's prism, 2D regions and quadrature and S9d.3a's
classification.
"""
from fractions import Fraction as F
import itertools
import math

import mpmath as mp

from curve_surface_reference import stored_axes
from sphere_boolean_reference import (
    M, Mv, Disc2, Polys2, Prism, add, clip_poly, cross, cross2, det3, dot, hull2, integrate, is_zero, nodes,
    roots2, scale, solve, sub, tau, vadd, zero3)
from boolean_reference import squarefree
from cone_boolean_reference import ccw2, classify, convex_disc, region_area, region_classes, seg_circle_params

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
TWO_PI = 6.283185307179586


# ------------------------------------------------------------------ polynomials

def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0)+(q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def pscale(p, k):
    return [c*k for c in p]


def real_roots(coeffs):
    """The real roots (mpf) of an exact polynomial (ascending
    coefficients): each factor of its square-free factorization (Yun's, in
    Fractions: a tangency's double root is a simple root of a factor) of
    degree two or less exactly, else by mpmath's `polyroots` at 40 digits
    with a near-real pair (imaginary part within 1e-15) taken as a real
    root, each polished by Newton's method on the exact coefficients."""
    c = [F(v) for v in coeffs]
    while c and c[-1] == 0:
        c.pop()
    if len(c) <= 3:
        return roots2(c)
    out = []
    for factor, _ in squarefree(c):
        out += simple_roots(factor)
    return sorted(out)


def simple_roots(c):
    """The real roots of a square-free exact polynomial."""
    if len(c) <= 3:
        return roots2(c)
    cm = [M(v) for v in c]
    scale_ = max(abs(v) for v in cm)
    desc = [v/scale_ for v in reversed(cm)]
    try:
        rts = mp.polyroots(desc, maxsteps=200, extraprec=200)
    except mp.NoConvergence:
        rts = mp.polyroots(desc, maxsteps=800, extraprec=400, error=False)
    out = []
    for z in rts:
        z = mp.mpc(z)
        if abs(z.imag) <= mp.mpf(10)**-15*(1+abs(z)):
            x = z.real
            for _ in range(4):
                f = mp.polyval(desc, x)
                d = mp.polyval([v*(len(desc)-1-k) for k, v in enumerate(desc[:-1])], x)
                if d == 0:
                    break
                step = f/d
                if abs(step) > mp.mpf(10)**-10*(1+abs(x)):
                    break
                x -= step
            out.append(x)
    return sorted(out)


# ------------------------------------------------------------------ the torus

def whole_turn(case):
    _, _, low, high, angle = case.torus
    return high-low == TWO_PI and angle == TWO_PI


class Torus:
    """A whole torus on its exact model, in its chart `(u, v, w)`."""

    def __init__(self, case):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        R, r, _, _, _ = case.torus
        assert whole_turn(case), 'S9d.4a: a whole torus (the full tube and turn)'
        self.R, self.r = F(R), F(r)
        assert self.R > self.r > 0, 'S9d.4a: a ring torus'
        self.case = case
        self.o, self.x, self.y, self.n = o, x, y, n
        self.det = det3(x, y, n)
        assert self.det > 0
        self.inv = [scale(cross(y, n), 1/self.det), scale(cross(n, x), 1/self.det),
                    scale(cross(x, y), 1/self.det)]
        a, b, c = cross(y, n), cross(n, x), cross(x, y)
        self.cof = (Mv(a), Mv(b), Mv(c))
        self.orthonormal = (dot(a, b) == 0 and dot(a, c) == 0 and dot(b, c) == 0
                            and dot(a, a) == dot(b, b) == dot(c, c))
        self.N0 = mp.sqrt(M(dot(a, a)))
        self.Rm, self.rm = M(self.R), M(self.r)
        self.om = Mv(o)

    def chart(self, X):
        d = sub(X, self.o)
        return tuple(dot(row, d) for row in self.inv)

    def covector(self, a, b):
        """The world half-space `a . X >= b` as the chart's `ac . p >= bc`."""
        ac = tuple(dot(a, col) for col in (self.x, self.y, self.n))
        return ac, b-dot(a, self.o)

    def to_world_moments(self, V, m):
        det = M(self.det)
        Vw = det*V
        axes = (Mv(self.x), Mv(self.y), Mv(self.n))
        return Vw, tuple(self.om[i]*Vw+det*sum(m[j]*axes[j][i] for j in range(3)) for i in range(3))

    def F(self, p):
        """The torus's function (negative inside), exact."""
        R2, r2 = self.R*self.R, self.r*self.r
        s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]
        return (s+R2-r2)**2-4*R2*(p[0]*p[0]+p[1]*p[1])

    def line_quartic(self, X0, l):
        """`F(X0 + t l)` as ascending coefficients in `t` (exact)."""
        R2, r2 = self.R*self.R, self.r*self.r
        full = [dot(X0, X0)+R2-r2, 2*dot(X0, l), dot(l, l)]
        hz = [X0[0]*X0[0]+X0[1]*X0[1], 2*(X0[0]*l[0]+X0[1]*l[1]), l[0]*l[0]+l[1]*l[1]]
        return padd(pmul(full, full), pscale(hz, -4*R2))

    def factor(self, theta, phi):
        """`|cof(A) N|` at `(theta, phi)`."""
        if self.orthonormal:
            return self.N0
        c = mp.cos(phi)
        N = (c*mp.cos(theta), c*mp.sin(theta), mp.sin(phi))
        a, b, cc = self.cof
        v = tuple(N[0]*a[i]+N[1]*b[i]+N[2]*cc[i] for i in range(3))
        return mp.sqrt(dot(v, v))

    def _gl(self, f, t0, t1):
        total = mp.mpf(0)
        pieces = max(1, int(mp.ceil((t1-t0)/(mp.pi/4))))
        step = (t1-t0)/pieces
        for i in range(pieces):
            lo = t0+i*step
            for xi, wi in nodes(4):
                total += wi*step/2*f(lo+step*(xi+1)/2)
        return total

    def meridian_weight(self, theta, p0, p1):
        """The wall's area over the tube's arc `phi` in `[p0, p1]` at
        `theta`, per `dtheta`."""
        R, r = self.Rm, self.rm
        if self.orthonormal:
            return self.N0*r*(R*(p1-p0)+r*(mp.sin(p1)-mp.sin(p0)))
        return self._gl(lambda p: r*(R+r*mp.cos(p))*self.factor(theta, p), p0, p1)

    def parallel_weight(self, phi, t0, t1):
        """The integral of `|cof(A) N|` over `theta` in `[t0, t1]` at `phi`."""
        if self.orthonormal:
            return self.N0*(t1-t0)
        return self._gl(lambda t: self.factor(t, phi), t0, t1)

    def wall_area(self):
        """The whole wall's area: `4 pi^2 R r` times the factor, or the
        trapezoidal rule on the periodic integrand (exponentially
        accurate)."""
        R, r = self.Rm, self.rm
        if self.orthonormal:
            return self.N0*4*mp.pi**2*R*r
        m = 48
        h = 2*mp.pi/m
        total = mp.mpf(0)
        for i in range(m):
            th = i*h
            for j in range(m):
                ph = j*h
                total += r*(R+r*mp.cos(ph))*self.factor(th, ph)
        return total*h*h

    def closed(self):
        """Volume, first moments and area (world) in closed form."""
        V = 2*mp.pi**2*self.Rm*self.rm**2
        Vw, mom = self.to_world_moments(V, (0, 0, 0))
        return Vw, mom, self.wall_area()

    def face_factor(self, a):
        """World area per chart area of a plane of chart normal `a`."""
        v = add(add(scale(cross(self.y, self.n), a[0]), scale(cross(self.n, self.x), a[1])),
                scale(cross(self.x, self.y), a[2]))
        return mp.sqrt(M(dot(v, v)))/mp.sqrt(M(dot(a, a)))

    def contains(self, p):
        """Strictly inside (mpf chart point)."""
        R2, r2 = self.Rm**2, self.rm**2
        s = p[0]**2+p[1]**2+p[2]**2
        return (s+R2-r2)**2 < 4*R2*(p[0]**2+p[1]**2)


def circle_tangency_quartic(T, a, b):
    """The levels `s` where the line `a_u u + a_v v = b - a_w s` is tangent to
    a circle of the torus's slice (either radius): ascending coefficients."""
    au, av, aw = a
    N = au*au+av*av
    R2, r2 = T.R*T.R, T.r*T.r
    D2 = [b*b, -2*b*aw, aw*aw]
    E = padd(D2, [N*(R2-r2), F(0), N])
    return padd(pmul(E, E), pscale(D2, -4*R2*N))


# ------------------------------------------------------------------ normal slices

class Normal:
    """Both solids sliced by the chart planes `w = s`."""

    def __init__(self, T, P, size):
        self.T, self.P, self.size = T, P, size
        self.pc = {key: T.chart(X) for key, X in P.P.items()}
        self.pm = {key: Mv(p) for key, p in self.pc.items()}
        self.level = {key: p[2] for key, p in self.pc.items()}
        self.quad_error = mp.mpf(0)
        self.breaks = None
        self._res = None

    def discs(self, s):
        T = self.T
        q2 = T.rm**2-s*s
        if q2 <= 0:
            empty = Disc2((mp.mpf(0), mp.mpf(0)), mp.mpf(-1))
            return empty, empty, mp.mpf(0)
        q = mp.sqrt(q2)
        o = Disc2((mp.mpf(0), mp.mpf(0)), (T.Rm+q)**2)
        i = Disc2((mp.mpf(0), mp.mpf(0)), (T.Rm-q)**2)
        return o, i, q

    def sections(self, s):
        out = []
        for piece in self.P.pieces:
            pts = []
            for i, j in piece.edges:
                li, lj = M(self.level[i]), M(self.level[j])
                if (li < s < lj) or (lj < s < li):
                    t = (s-li)/(lj-li)
                    X = add(self.pm[i], scale(sub(self.pm[j], self.pm[i]), t))
                    pts.append((X[0], X[1]))
            out.append(hull2(pts))
        return out

    def integrand(self, s):
        Do, Di, q = self.discs(s)
        polys = self.sections(s)
        P2 = Polys2([p for p in polys if p])
        out = []
        co = classify(Do, P2)
        ci = classify(Di, P2)
        for c in (co, ci):
            for g in c[:4]:
                out += [g[0], s*g[0], g[1], g[2]]
        conv = zero3()
        for poly in polys:
            if poly:
                conv = vadd(conv, vadd(convex_disc(poly, Do), convex_disc(poly, Di), -1))
        out += [conv[0], s*conv[0], conv[1], conv[2]]
        w_in = w_all = mp.mpf(0)
        if q > 0:
            T = self.T
            for c, sign in ((co, 1), (ci, -1)):
                rho = T.Rm+sign*q
                phi = mp.atan2(s, sign*q)
                for t0, t1, inside in c[4]:
                    wgt = T.rm*rho/q*T.parallel_weight(phi, t0, t1)
                    w_all += wgt
                    if inside:
                        w_in += wgt
        out += [w_in, w_all]
        return out

    def breakpoints(self):
        T, P = self.T, self.P
        levels = set(self.level.values())
        levels |= {-T.r, T.r}
        polys = []
        edges = set(e for piece in P.pieces for e in piece.edges)
        for i, j in edges:
            vi, vj = self.pc[i], self.pc[j]
            e = sub(vj, vi)
            if e[2] == 0:
                continue
            # F(vi + t e) with t = (s - vi_w) / e_w.
            quart = T.line_quartic(vi, e)
            p0, p1 = -vi[2]/e[2], 1/e[2]
            lin = [p0, p1]
            poly = [F(0)]
            powk = [F(1)]
            for c in quart:
                poly = padd(poly, pscale(powk, c))
                powk = pmul(powk, lin)
            polys.append(poly)
        for piece in P.pieces:
            for f in piece.faces:
                p0, p1, p2 = (self.pc[k] for k in f['loop'][:3])
                a = cross(sub(p1, p0), sub(p2, p0))
                if a[0] == 0 and a[1] == 0:
                    continue
                polys.append(circle_tangency_quartic(T, a, dot(a, p0)))
        pts = [M(v) for v in levels]
        for p in polys:
            pts += real_roots(p)
        lo = min([M(v) for v in self.level.values()]+[-T.rm])
        hi = max([M(v) for v in self.level.values()]+[T.rm])
        pts = sorted(p for p in pts if lo <= p <= hi)
        out = []
        gap = mp.mpf(10)**-30*self.size
        for p in [lo]+pts+[hi]:
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
        `fuse`, `K-P`, `P-K`, `common2` (by clipping), and the wall's area
        inside the prism and outside (`wall`)."""
        if self._res is not None:
            return self._res
        v = self.measure()
        G = {}
        names = ('oin', 'oout', 'opin', 'opout', 'iin', 'iout', 'ipin', 'ipout', 'conv')
        for k, name in enumerate(names):
            G[name] = v[4*k:4*k+4]
        combos = {'K': [('oin', 1), ('oout', 1), ('iin', -1), ('iout', -1)],
                  'P': [('opin', 1), ('opout', 1)],
                  'common': [('oin', 1), ('opin', 1), ('iin', -1), ('ipin', -1)],
                  'fuse': [('oout', 1), ('opout', 1), ('iout', -1), ('ipin', 1)],
                  'K-P': [('oout', 1), ('opin', -1), ('iout', -1), ('ipin', 1)],
                  'P-K': [('opout', 1), ('oin', -1), ('iin', 1), ('ipin', 1)],
                  'common2': [('conv', 1)]}
        out = {}
        for name, parts in combos.items():
            A, sA, ma, mb = [sum(k*G[p][i] for p, k in parts) for i in range(4)]
            out[name] = self.T.to_world_moments(A, (ma, mb, sA))
        w_in, w_all = v[36], v[37]
        out['wall'] = {'in': w_in, 'out': w_all-w_in, 'total': w_all}
        self._res = out
        return out


# ------------------------------------------------------------------ meridian half-planes

def green_mer_seg(p, q):
    """`(t^2 dw / 2, t^3 dw / 3, t^2 w dw / 2)` over the segment from `p` to
    `q` in `(t, w)` (exact: three-point Gauss-Legendre on cubics)."""
    dt, dw = q[0]-p[0], q[1]-p[1]
    out = [mp.mpf(0)]*3
    k = mp.sqrt(mp.mpf(3)/5)
    for x, wgt in ((-k, mp.mpf(5)/9), (mp.mpf(0), mp.mpf(8)/9), (k, mp.mpf(5)/9)):
        lam = (x+1)/2
        t, w = p[0]+lam*dt, p[1]+lam*dw
        wgt = wgt/2*dw
        out[0] += wgt*t*t/2
        out[1] += wgt*t**3/3
        out[2] += wgt*t*t*w/2
    return tuple(out)


def green_mer_arc(R, r, a, b):
    """The same over the circle `(R + r cos phi, r sin phi)`, `phi` from `a` to
    `b` (closed-form antiderivatives)."""
    def I(p):
        s, c = mp.sin(p), mp.cos(p)
        ic = s
        ic2 = p/2+mp.sin(2*p)/4
        ic3 = s-s**3/3
        ic4 = 3*p/8+mp.sin(2*p)/4+mp.sin(4*p)/32
        isc = -c*c/2
        isc2 = -c**3/3
        isc3 = -c**4/4
        f1 = r/2*(R*R*ic+2*R*r*ic2+r*r*ic3)
        f2 = r/3*(R**3*ic+3*R*R*r*ic2+3*R*r*r*ic3+r**3*ic4)
        f3 = r*r/2*(R*R*isc+2*R*r*isc2+r*r*isc3)
        return (f1, f2, f3)
    ib, ia = I(b), I(a)
    return tuple(x-y for x, y in zip(ib, ia))


class Meridian:
    """Both solids sliced by the half-planes `theta = const` about the axis."""

    def __init__(self, T, P, size):
        self.T, self.P, self.size = T, P, size
        self.pc = {key: T.chart(X) for key, X in P.P.items()}
        self.pm = {key: Mv(p) for key, p in self.pc.items()}
        # Faces as chart half-spaces a . p >= b holding the piece.
        self.faces = []
        for pi, piece in enumerate(P.pieces):
            for f in piece.faces:
                ac, bc = T.covector(f['a'], f['b'])
                self.faces.append({'piece': pi, 'tag': f['tag'], 'a': ac, 'b': bc,
                                   'am': Mv(ac), 'bm': M(bc)})
        self.C = (T.Rm, mp.mpf(0))
        self.quad_error = mp.mpf(0)
        self.breaks = None
        self._res = None

    # -- sections

    def piece_section(self, pi, theta):
        """The piece's section by the half-plane: [((t, w), id)] counter-
        clockwise, `id` ('edge', key pair) or ('axis', face index)."""
        piece = self.P.pieces[pi]
        c, s = mp.cos(theta), mp.sin(theta)
        sig = {k: -s*self.pm[k][0]+c*self.pm[k][1] for k in piece.keys}
        # A vertex on the plane (on the axis: every plane's) is a vertex of
        # the section.
        pts = [((c*self.pm[k][0]+s*self.pm[k][1], self.pm[k][2]), ('vertex', k)) for k in piece.keys if sig[k] == 0]
        for i, j in piece.edges:
            si, sj = sig[i], sig[j]
            if (si < 0 < sj) or (sj < 0 < si):
                lam = si/(si-sj)
                X = add(self.pm[i], scale(sub(self.pm[j], self.pm[i]), lam))
                pts.append(((c*X[0]+s*X[1], X[2]), ('edge', (i, j))))
        if len(pts) < 3:
            return []
        mt = sum(p[0][0] for p in pts)/len(pts)
        mw = sum(p[0][1] for p in pts)/len(pts)
        pts.sort(key=lambda p: mp.atan2(p[0][1]-mw, p[0][0]-mt))
        # Clip at the axis (t >= 0), the new vertices named by the face
        # shared by the two edges whose crossings they join.
        out = []
        n = len(pts)
        for k in range(n):
            (p, ip), (q, iq) = pts[k], pts[(k+1) % n]
            if p[0] >= 0:
                out.append((p, ip))
            if (p[0] < 0 < q[0]) or (q[0] < 0 < p[0]):
                lam = p[0]/(p[0]-q[0])
                w = p[1]+lam*(q[1]-p[1])
                out.append(((mp.mpf(0), w), ('axis', self.shared_face(pi, ip, iq))))
        return out if len(out) >= 3 else []

    def shared_face(self, pi, a, b):
        """The index (in `faces`) of the piece's face holding both section
        vertices (ids `('edge', e)` or `('vertex', k)`)."""
        piece = self.P.pieces[pi]
        for f in piece.faces:
            L = f['loop']
            es = {tuple(sorted((L[t], L[(t+1) % len(L)]))) for t in range(len(L))}
            holds = lambda v: v[1] in es if v[0] == 'edge' else v[1] in L
            if holds(a) and holds(b):
                return next(k for k, g in enumerate(self.faces) if g['piece'] == pi and g['tag'] == f['tag'])
        raise AssertionError('section vertices without a shared face')

    def sections(self, theta):
        return [self.piece_section(pi, theta) for pi in range(len(self.P.pieces))]

    # -- classification in the half-plane

    def pieces2d(self, polys):
        """Boundary pieces of the disc and of the polygons, classified:
        (disc arcs [(p0, p1, inside)], polygon segments [(p, q, inside)])."""
        R, r = self.T.Rm, self.T.rm
        plines = []
        for poly in polys:
            ls = []
            n = len(poly)
            for k in range(n):
                p, q = poly[k][0], poly[(k+1) % n][0]
                al, be = -(q[1]-p[1]), q[0]-p[0]
                ls.append((al, be, al*p[0]+be*p[1]))
            plines.append(ls)
        lines = [ln for ls in plines for ln in ls]
        in_poly = lambda X: any(all(al*X[0]+be*X[1] > ga for al, be, ga in ls) for ls in plines)
        in_disc = lambda X: (X[0]-R)**2+X[1]**2 < r*r
        T = tau()
        angs = []
        for al, be, ga in lines:
            nn = mp.sqrt(al*al+be*be)
            if nn == 0:
                continue
            val = (ga-al*R)/(r*nn)
            if abs(val) < 1:
                phi, w = mp.atan2(be, al), mp.acos(val)
                angs += [(phi-w) % T, (phi+w) % T]
        arcs = []
        if not angs:
            arcs.append((mp.mpf(0), T, in_poly((R+r, mp.mpf(0)))))
        else:
            angs.sort()
            for k in range(len(angs)):
                a = angs[k]
                b = angs[k+1] if k+1 < len(angs) else angs[0]+T
                if b-a <= 0:
                    continue
                m = (a+b)/2
                arcs.append((a, b, in_poly((R+r*mp.cos(m), r*mp.sin(m)))))
        segs = []
        for poly in polys:
            n = len(poly)
            for k in range(n):
                p, q = poly[k][0], poly[(k+1) % n][0]
                ts = [mp.mpf(0), mp.mpf(1)]
                ts += [t for t in seg_circle_params(p, q, self.C, r*r) if 0 < t < 1]
                ts.sort()
                for t0, t1 in zip(ts, ts[1:]):
                    if t1-t0 <= 0:
                        continue
                    a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
                    b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
                    segs.append((a, b, in_disc(((a[0]+b[0])/2, (a[1]+b[1])/2))))
        return arcs, segs

    def integrand(self, theta):
        T = self.T
        R, r = T.Rm, T.rm
        polys = [p for p in self.sections(theta) if p]
        arcs, segs = self.pieces2d(polys)
        c, s = mp.cos(theta), mp.sin(theta)
        acc = {k: [mp.mpf(0)]*3 for k in ('din', 'dout', 'pin', 'pout')}
        w_in = w_all = mp.mpf(0)
        for a, b, inside in arcs:
            g = green_mer_arc(R, r, a, b)
            key = 'din' if inside else 'dout'
            acc[key] = [x+y for x, y in zip(acc[key], g)]
            wgt = T.meridian_weight(theta, a, b)
            w_all += wgt
            if inside:
                w_in += wgt
        for p, q, inside in segs:
            g = green_mer_seg(p, q)
            key = 'pin' if inside else 'pout'
            acc[key] = [x+y for x, y in zip(acc[key], g)]
        out = []
        for key in ('din', 'dout', 'pin', 'pout'):
            m1, m2, mtw = acc[key]
            out += [m1, c*m2, s*m2, mtw]
        out += [w_in, w_all]
        return out

    # -- breakpoints

    def line_meetings(self, X0, l, t_range=None):
        """Angles of the points where the line `X0 + t l` meets the torus
        (within `t_range` when given)."""
        out = []
        for t in real_roots(self.T.line_quartic(X0, l)):
            if t_range is not None and not (t_range[0]-mp.mpf(10)**-30 <= t <= t_range[1]+mp.mpf(10)**-30):
                continue
            X = tuple(M(X0[k])+t*M(l[k]) for k in range(3))
            if X[0] != 0 or X[1] != 0:
                out.append(mp.atan2(X[1], X[0]))
        return out

    def breakpoints(self):
        T, P = self.T, self.P
        pi = mp.pi
        angs = []
        for key, p in self.pc.items():
            if p[0] != 0 or p[1] != 0:
                a = mp.atan2(M(p[1]), M(p[0]))
                angs += [a, a+pi]
        edges = set(e for piece in P.pieces for e in piece.edges)
        for i, j in edges:
            vi, vj = self.pc[i], self.pc[j]
            angs += self.line_meetings(vi, sub(vj, vi), (0, 1))
        R, r = T.R, T.r
        for f in self.faces:
            au, av, aw = f['a']
            b = f['b']
            if au == 0 and av == 0:
                continue
            rho = mp.sqrt(M(au*au+av*av))
            phi = mp.atan2(M(av), M(au))
            angs += [phi+pi/2, phi-pi/2]
            for alpha in roots2([b*b-r*r*aw*aw, -2*R*b, R*R-r*r]):
                if abs(alpha) <= rho:
                    w = mp.acos(alpha/rho)
                    angs += [phi+w, phi-w]
        for f, g in itertools.combinations(self.faces, 2):
            l = cross(f['a'], g['a'])
            if is_zero(l):
                continue
            X0 = solve([f['a'], g['a'], l], [f['b'], g['b'], F(0)])
            angs += self.line_meetings(X0, l)
        T2 = 2*pi
        norm = sorted(((a+pi) % T2)-pi for a in angs)
        pts = [-pi]+norm+[pi]
        out = []
        gap = mp.mpf(10)**-30
        for p in sorted(pts):
            if not out or p-out[-1] > gap:
                out.append(p)
        if pi-out[-1] <= gap:
            out[-1] = pi
        else:
            out.append(pi)
        self.raw_breaks = norm
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
        """As `Normal.results` (no `common2`): the volumes and moments by the
        cylindrical element, the wall by its own parameters."""
        if self._res is not None:
            return self._res
        v = self.measure()
        G = {name: v[4*k:4*k+4] for k, name in enumerate(('din', 'dout', 'pin', 'pout'))}
        combos = {'K': [('din', 1), ('dout', 1)], 'P': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'K-P': [('dout', 1), ('pin', -1)], 'P-K': [('pout', 1), ('din', -1)]}
        out = {}
        for name, parts in combos.items():
            V, mu, mv, mw = [sum(k*G[p][i] for p, k in parts) for i in range(4)]
            out[name] = self.T.to_world_moments(V, (mu, mv, mw))
        out['wall'] = {'in': v[16], 'out': v[17]-v[16], 'total': v[17]}
        self._res = out
        return out

    # -- solids

    def half_lines(self, theta, which):
        """The faces' lines in the half-plane at `theta` as `(al, be, ga)`
        (`al t + be w > ga`): `inner` (the piece's side) or `outer` for
        each face; None for a face whose plane holds the half-plane's plane
        (the condition constant: True or False)."""
        c, s = mp.cos(theta), mp.sin(theta)
        out = []
        for f in self.faces:
            au, av, aw = f['am']
            alpha = au*c+av*s
            line = (alpha, aw, f['bm']) if which == 'inner' else (-alpha, -aw, -f['bm'])
            if aw == 0 and abs(alpha) <= mp.mpf(10)**-30*(abs(au)+abs(av)):
                # The face's plane holds the half-plane's: the condition
                # `0 > ga` is constant.
                out.append(0 > line[2])
            else:
                out.append(line)
        return out

    def disc_area(self, lines):
        ls = []
        for ln in lines:
            if ln is True:
                continue
            if ln is False:
                return mp.mpf(0)
            ls.append(ln)
        r = self.T.rm
        crossing = False
        for al, be, ga in ls:
            nn = mp.sqrt(al*al+be*be)
            if nn > 0 and abs((ga-al*self.C[0]-be*self.C[1])/(r*nn)) < 1:
                crossing = True
        if not crossing:
            # No line crosses the circle (a tangent one touches it at a
            # point): the whole disc or nothing, by its centre.
            inside = all(al*self.C[0]+be*self.C[1] > ga for al, be, ga in ls)
            return mp.pi*r*r if inside else mp.mpf(0)
        return region_area(Disc2(self.C, r*r, ls))

    def structure(self, theta):
        """The components of the common (whether it has area), of `D - C`
        (sets of faces) and of `C - D` (sets of vertex ids) at `theta`."""
        eps = mp.mpf(10)**-25*self.size**2
        inner = self.half_lines(theta, 'inner')
        outer = self.half_lines(theta, 'outer')
        common = self.disc_area(inner) > eps
        nf = len(self.faces)
        present = [f for f in range(nf) if self.disc_area([outer[f]]) > eps]
        parent = {f: f for f in present}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k
        for f, g in itertools.combinations(present, 2):
            if self.disc_area([outer[f], outer[g]]) > eps:
                parent[find(f)] = find(g)
        groups = {}
        for f in present:
            groups.setdefault(find(f), set()).add(f)
        k_minus = [frozenset(g) for g in groups.values()]
        return common, k_minus, self.runs(theta)

    def runs(self, theta):
        """`C - D`'s components as sets of the section's vertex ids."""
        R, r = self.T.Rm, self.T.rm
        poly = self.piece_section(0, theta)
        if not poly:
            return []
        n = len(poly)
        items = []   # (inside, id of the start vertex or None)
        for k in range(n):
            (p, ip), (q, _) = poly[k], poly[(k+1) % n]
            ts = [mp.mpf(0), mp.mpf(1)]
            ts += [t for t in seg_circle_params(p, q, self.C, r*r) if 0 < t < 1]
            ts.sort()
            first = True
            for t0, t1 in zip(ts, ts[1:]):
                if t1-t0 <= 0:
                    continue
                tm = (t0+t1)/2
                X = (p[0]+tm*(q[0]-p[0]), p[1]+tm*(q[1]-p[1]))
                inside = (X[0]-R)**2+X[1]**2 < r*r
                items.append((inside, ip if first else None))
                first = False
        if all(i for i, _ in items):
            return []
        if not any(i for i, _ in items):
            return [frozenset(ip for _, ip in poly)]
        k0 = next(k for k, (i, _) in enumerate(items) if i)
        items = items[k0+1:]+items[:k0+1]
        out, cur = [], None
        for inside, vid in items:
            if inside:
                if cur is not None:
                    out.append(frozenset(cur))
                    cur = None
                continue
            if cur is None:
                cur = set()
            if vid is not None:
                cur.add(vid)
        if cur is not None:
            out.append(frozenset(cur))
        assert all(out), 'a run of the section outside the disc without a vertex'
        return out

    def vertex_position(self, vid, theta):
        c, s = mp.cos(theta), mp.sin(theta)
        if vid[0] == 'axis':
            f = self.faces[vid[1]]
            au, av, aw = f['am']
            return (mp.mpf(0), f['bm']/aw)
        if vid[0] == 'vertex':
            X = self.pm[vid[1]]
            return (c*X[0]+s*X[1], X[2])
        i, j = vid[1]
        sig = lambda k: -s*self.pm[k][0]+c*self.pm[k][1]
        si, sj = sig(i), sig(j)
        lam = si/(si-sj)
        X = add(self.pm[i], scale(sub(self.pm[j], self.pm[i]), lam))
        return (c*X[0]+s*X[1], X[2])

    def solids(self):
        """{'common', 'K-P', 'P-K'}: component counts (see the module's
        note)."""
        assert len(self.P.pieces) == 1, 'S9d.4a reference: solids of a convex prism only'
        if self.breaks is None:
            self.breaks = self.breakpoints()
        br = self.breaks
        ivs = list(zip(br, br[1:]))
        n = len(ivs)
        R, r = self.T.Rm, self.T.rm
        eps = mp.mpf(10)**-25*self.size**2
        tol = mp.mpf(10)**-12*self.size
        per = [self.structure((a+b)/2) for a, b in ivs]
        parent = {}

        def node(x):
            parent.setdefault(x, x)
            return x

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x

        def join(x, y):
            parent[find(x)] = find(y)
        for i, (common, km, pk) in enumerate(per):
            if common:
                node(('common', i, 0))
            for k in range(len(km)):
                node(('K-P', i, k))
            for k in range(len(pk)):
                node(('P-K', i, k))
        for i in range(n):
            j = (i+1) % n
            tb = ivs[i][1]
            delta = min(mp.mpf(10)**-20, (ivs[i][1]-ivs[i][0])/4, (ivs[j][1]-ivs[j][0])/4)
            if per[i][0] and per[j][0] and self.disc_area(self.half_lines(tb, 'inner')) > eps:
                join(('common', i, 0), ('common', j, 0))
            outer = self.half_lines(tb, 'outer')
            for ka, A in enumerate(per[i][1]):
                for kb, B in enumerate(per[j][1]):
                    if any(self.disc_area([outer[f]]) > eps for f in A & B):
                        join(('K-P', i, ka), ('K-P', j, kb))
            tl, tr = tb-delta, ivs[j][0]+delta
            left = [(ka, self.vertex_position(v, tl)) for ka, A in enumerate(per[i][2]) for v in A]
            right = [(kb, self.vertex_position(v, tr)) for kb, B in enumerate(per[j][2]) for v in B]
            for ka, X in left:
                if (X[0]-R)**2+X[1]**2 <= (r+tol)**2:
                    continue
                for kb, Y in right:
                    if abs(X[0]-Y[0])+abs(X[1]-Y[1]) < tol:
                        join(('P-K', i, ka), ('P-K', j, kb))
        counts = {}
        for op in ('common', 'K-P', 'P-K'):
            counts[op] = len({find(x) for x in parent if x[0] == op})
        return counts


# ------------------------------------------------------------------ planar faces

def slanted_face(T, polys, a, b, size):
    """(inside, all) chart areas of convex chart polygons in the plane `a . p
    = b` (not normal to the axis), sliced by the lines of constant `w`."""
    au, av, aw = a
    nuv2 = au*au+av*av
    spacing = mp.sqrt(M(dot(a, a)))/mp.sqrt(M(nuv2))
    levels = {-T.r, T.r}
    polys_ = [circle_tangency_quartic(T, a, b)]
    for poly in polys:
        for i in range(len(poly)):
            p, q = poly[i], poly[(i+1) % len(poly)]
            levels.add(p[2])
            e = sub(q, p)
            if e[2] != 0:
                quart = T.line_quartic(p, e)
                lin = [-p[2]/e[2], 1/e[2]]
                acc, powk = [F(0)], [F(1)]
                for c in quart:
                    acc = padd(acc, pscale(powk, c))
                    powk = pmul(powk, lin)
                polys_.append(acc)
    lo = min(p[2] for poly in polys for p in poly)
    hi = max(p[2] for poly in polys for p in poly)
    pts = [M(v) for v in levels if lo <= v <= hi]
    for c in polys_:
        pts += [t for t in real_roots(c) if M(lo) <= t <= M(hi)]
    pts.sort()
    brk = []
    for p in pts:
        if not brk or p-brk[-1] > mp.mpf(10)**-30*size:
            brk.append(p)
    pm = [[Mv(p) for p in poly] for poly in polys]
    rm, Rm = T.rm, T.Rm
    zero = (mp.mpf(0), mp.mpf(0))

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
            q2 = rm*rm-s*s
            if q2 <= 0:
                continue
            q = mp.sqrt(q2)
            for rho, sign in ((Rm+q, 1), (Rm-q, -1)):
                ts = seg_circle_params(P1, P2, zero, rho*rho)
                if len(ts) == 2:
                    Lin += sign*L*max(mp.mpf(0), min(ts[1], mp.mpf(1))-max(ts[0], mp.mpf(0)))
        return [Lin, Lall]
    tol = mp.mpf(10)**-33*size**3
    tin = tall = mp.mpf(0)
    for x0, x1 in zip(brk, brk[1:]):
        est, _ = integrate(f, x0, x1, tol)
        if est is not None:
            tin += est[0]
            tall += est[1]
    return spacing*tin, spacing*tall


def face_classes(T, polys, outward, size):
    """{in, out, same, opp} (world areas) of a planar face of the prism
    (convex chart polygons, its chart outward normal) against the torus."""
    p0, p1, p2 = polys[0][0], polys[0][1], polys[0][2]
    a = cross(sub(p1, p0), sub(p2, p0))
    if dot(a, outward) < 0:
        a = scale(a, -1)
    b = dot(a, p0)
    factor = T.face_factor(a)
    cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
    if a[0] == 0 and a[1] == 0:
        w0 = b/a[2]
        p2d = [ccw2([(M(p[0]), M(p[1])) for p in poly]) for poly in polys]
        total = sum(abs(sum(cross2(p[i-1], p[i]) for i in range(len(p)))/2) for p in p2d)
        q2 = T.r*T.r-w0*w0
        if q2 <= 0:
            cls['out'] = total*factor
            return cls
        q = mp.sqrt(M(q2))
        o_in, o_out = region_classes(Disc2((mp.mpf(0), mp.mpf(0)), (T.Rm+q)**2), p2d)
        i_in, _ = region_classes(Disc2((mp.mpf(0), mp.mpf(0)), (T.Rm-q)**2), p2d)
        cls['in'] = (o_in-i_in)*factor
        cls['out'] = (o_out+i_in)*factor
        return cls
    inside, total = slanted_face(T, polys, a, b, size)
    cls['in'] = inside*factor
    cls['out'] = (total-inside)*factor
    return cls


def prism_face_classes(T, prism, size):
    """[(tag, {in, out, same, opp}, exact area)] for the prism's faces."""
    out = []
    for tag, polys, piece, area in prism.faces():
        cp = [[T.chart(X) for X in poly] for poly in polys]
        p0, p1, p2 = cp[0][0], cp[0][1], cp[0][2]
        a = cross(sub(p1, p0), sub(p2, p0))
        if dot(a, T.chart(piece.centroid)) > dot(a, p0):
            a = scale(a, -1)
        out.append((tag, face_classes(T, cp, a, size), area))
    return out


# ------------------------------------------------------------------ the pair

def case_size(T, P):
    vals = [F(1)]
    vals += [abs(x) for p in P.P.values() for x in p]
    vals += [abs(x)+T.R+T.r for x in T.o]
    return M(max(vals))


class Pair:
    """A torus and a prism, in either order (`obj`, `tool` Boolean cases)."""

    def __init__(self, obj, tool):
        self.torus_first = obj.torus is not None
        t_case, p_case = (obj, tool) if self.torus_first else (tool, obj)
        assert t_case.torus is not None and p_case.torus is None and p_case.sphere is None \
            and p_case.cone is None, 'S9d.4a: a torus and a prism'
        self.T, self.P = Torus(t_case), Prism(p_case)
        self.size = case_size(self.T, self.P)
        self.normal = Normal(self.T, self.P, self.size)
        self.meridian = Meridian(self.T, self.P, self.size)
        self._faces = self._solids = None

    def sliced(self):
        return self.normal.results()

    def swept(self):
        return self.meridian.results()

    def faces(self):
        """[(input 'K' or 'P', tag, classes, exact area)]."""
        if self._faces is None:
            out = [('P', tag, cls, area) for tag, cls, area in prism_face_classes(self.T, self.P, self.size)]
            wall = self.swept()['wall']
            out.append(('K', ('wall',), {'in': wall['in'], 'out': wall['out'], 'same': mp.mpf(0),
                                         'opp': mp.mpf(0)}, self.T.wall_area()))
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
        return ('K', 'P') if self.torus_first else ('P', 'K')

    def area(self, op):
        a, b = self.roles()
        ca, cb = self.classes(a), self.classes(b)
        keep = {'fuse': (('out', 'same'), ('out',)), 'cut': (('out', 'opp'), ('in',)),
                'common': (('in', 'same'), ('in',))}[op]
        return sum(ca[k] for k in keep[0])+sum(cb[k] for k in keep[1])

    def volume(self, op):
        r = self.sliced()
        key = {'fuse': 'fuse', 'common': 'common', 'cut': 'K-P' if self.torus_first else 'P-K'}[op]
        return r[key]

    def solids(self, op):
        if self._solids is None:
            self._solids = self.meridian.solids()
        if op == 'fuse':
            return 1 if self.sliced()['common'][0] > mp.mpf(10)**-25*self.size**3 else 2
        if op == 'common':
            return self._solids['common']
        return self._solids['K-P' if self.torus_first else 'P-K']

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
