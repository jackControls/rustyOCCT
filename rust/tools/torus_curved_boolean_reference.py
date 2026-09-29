#!/usr/bin/env python3
"""Independent reference for S9d.4b.2 of REVIEW_NOTES.md: Booleans of a whole
torus (`Solid::torus_with` with the full tube and turn) against a prism with
arcs and circles (planar and cylindrical walls), a whole sphere
(`Solid::sphere_with`), a cone or frustum (`Solid::cone_with`) and another
whole torus, in any relative position. A torus and a quadric of revolution
about another axis meet in quartic curves of the torus's tube angle (a circle
of the tube against the quadric: a trigonometric polynomial of degree two in
the tube's angle); coaxial pairs in circles; two tori about different axes in
curves of degree eight in the tube's half-angle tangent.

**The inputs' exact models**, in rationals from the stored binary64 data
(`stored_axes`, the kernel's `Frame3::new`, not exactly orthonormal), each
the image of an ideal solid under its stored frame `X = o + u x + v y + w n`
(affine where the axes are not orthonormal): the torus S9d.4a's (`(|p|^2 +
R^2 - r^2)^2 <= 4 R^2 (u^2 + v^2)`), the cone S9d.3a's (`u^2 + v^2 <= (bottom
+ k w)^2`, `0 <= w <= h`), the sphere S9d.1's (`|X - o|^2 <= R^2` in the
world: its frame's axes do not enter), the prism S9d.2's (a convex profile of
segments and counter-clockwise arcs, or a circle, between its stored start
and end heights). Each input is a set of exact surfaces, each a polynomial
`f(X)` of degree one (planes), two (quadrics) or four (tori) in the world
point, and a membership test: the torus, sphere and cone inside where every
`f < 0`, the prism where its height lies between its caps and its profile
point inside the profile (its winding: a profile of arcs and segments is not
the intersection of their discs and half-planes). Nothing here uses the
kernel, a surface/surface intersection or an arrangement; shared with other
references are `stored`, `stored_axes`, S9d.1's vector helpers, Gauss-Legendre
nodes and adaptive quadrature, and `number`.

**Every face swept two ways.** Each face of both inputs is covered by two
families of curves on it, each curve a circle or ellipse `X(b) = X0 + C cos
b + S sin b` or a line `X(b) = X0 + b D` over the family's outer parameter
`a`: the tori's walls by the tube's circles (`a` the turn, `b` the tube's
angle) and by the parallels (`a` the tube's angle, `b` the turn); the sphere
by its meridians (half circles, `b` the latitude) and its parallels; a
cylindrical wall (an arc of the profile swept) by its generatrices (lines)
and its sections (arcs); a flat wall by its generatrices and its lines at
constant height; a cone's wall by its rulings and its parallels; every
planar face (the prism's caps, the cone's discs) by chords of two rational
directions in its plane (their ends where the line leaves the profile).
Along one curve every surface of the other input is a polynomial in `b` (a
line) or a trigonometric polynomial of the surface's degree (a circle: its
coefficients by the discrete Fourier transform of `2 d + 1` samples, exact
for that degree); its real roots are the unit-circle roots of `z^d f`
(`mp.polyroots` at 40 digits where roots cluster, else floating-point roots
polished by Newton's method in 40 digits), a pair within 1e-15 of the circle
taken as real (a spurious boundary costs nothing), two roots of one surface
closer than 1e-15 of the curve's range dropped (a tangency's double root
split by rounding: no sign change between them). The curve is cut at every
root and each piece classified at a point a golden fraction along it (not
its midpoint, where a symmetric tangency may sit) by the other input's
membership (`in` or `out`; the surfaces of two inputs never share a face
here). Along each piece the area element `|X_a x X_b|`, the volume `(X -
c0) . N / 3` and the first moments `(X_i - c0_i)^2 N_i / 2` (`N = X_a x
X_b` oriented outward: the divergence theorem) are integrated by 24-point
Gauss-Legendre on pieces of at most a quarter turn of a circle (exact for a
line's polynomials; the area element of an affine frame's curved face by
its own analytic square root). The results are the classes' sums: common
the `in` pieces of both, fuse the `out` pieces of both, cut the object's
`out` pieces and the tool's `in` pieces reversed; the volume and moments
come out of the divergence theorem, the areas directly.

The outer parameter is integrated between its events by S9d.1's adaptive
Gauss-Legendre after `a = a0 + (a1 - a0)(1 - cos t)/2` (square-root ends
smooth), to 1e-32 of the case's size to the fourth. An event is where the
curve's structure changes: the sequence of its pieces' classes (adjacent
equal classes merged) with the surface whose root bounds each, cyclic on a
closed curve. It changes where two roots of one surface merge (the curve
tangent to the section: a turning point of the meeting curve in that
family), where roots of two surfaces cross (the curve through an edge of the
other input), and where a root reaches a curve's end (a planar face's chord,
a meridian's pole, an arc's end). Events are found on a scan of 360 curves
per family and located by bisection on the structure to 1e-35 of the range;
the family's own kinks (a profile's vertices across a chord family, an
arc's ends) are breakpoints too. Every quadrature node's structure is
checked against its interval's (`missed`): a missing event is an error, not
a slow quadrature.

**Checks the generator makes.** Every operation's volume, moments and area
two ways (the two families of every face, independent roots and events),
both inputs' measures from their faces against their closed forms,
inclusion and exclusion (`fuse = A + B - common`, the cut's), the area
identity `area(fuse) + area(common) = area(A) + area(B)`.

**Solids** by a sweep of the torus's normal slices (S9d.4a's: the chart
planes `w = s`, the torus's section the annulus between `R - q` and `R + q`,
`q = sqrt(r^2 - s^2)`), the other input's section along 720 rays from the
axis in each of 240 slices: each ray's intervals of the operation's section
(the annulus's and the other's, from the roots of its surfaces along the ray
in binary64, classified at midpoints), intervals joined where they overlap
on adjacent rays of a slice, at the axis, and on the same ray of adjacent
slices; the connected components counted. Regions meeting along a curve or
at a point stay apart (overlaps of positive length only), the regularized
Boolean. The fuse is one solid where the inputs overlap, two otherwise; the
sweep's must agree.

**Margins.** The least sine between the two inputs' surfaces where they
meet (at the scanned curves' roots on the torus's tube circles), the least
distance between surfaces that come close without meeting (the other's
function's critical values along every scanned curve, `|f| / |grad f|`,
where no root is near), the least sine at which an edge (a prism's or a
cone's rim) crosses the torus and its least distance from tangency with it,
and a cone's apex's distance from the other input.

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9d.4b.1
references.
"""
from fractions import Fraction as F
import cmath
import math

import mpmath as mp

from curve_surface_reference import stored_axes
from identity_reference import stored
from sphere_boolean_reference import M, Mv, add, cross, dot, integrate, nodes, scale, sub
from torus_boolean_reference import number, OPS, TWO_PI

mp.mp.dps = 40

Z = mp.mpf(0)
HALF_PI = math.pi/2
SCAN = 360
INNER = 4          # 24-point Gauss-Legendre along a piece
QUARTER = None     # pi/4 at the working precision


def quarter():
    global QUARTER
    if QUARTER is None or QUARTER[0] != mp.mp.prec:
        QUARTER = (mp.mp.prec, mp.pi/2)
    return QUARTER[1]


def lin(*terms):
    """`sum c_i v_i` of mpf vectors."""
    out = [Z, Z, Z]
    for c, v in terms:
        out[0] += c*v[0]
        out[1] += c*v[1]
        out[2] += c*v[2]
    return tuple(out)


def norm(v):
    return mp.sqrt(dot(v, v))


# ------------------------------------------------------------------ frames

class Frame:
    """A stored frame: exact axes and their mpf images, the chart map."""

    def __init__(self, frame):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(frame))
        self.o, self.x, self.y, self.n = o, x, y, n
        self.det = dot(x, cross(y, n))
        assert self.det > 0
        inv = [scale(cross(y, n), 1/self.det), scale(cross(n, x), 1/self.det), scale(cross(x, y), 1/self.det)]
        self.inv = inv
        self.om, self.xm, self.ym, self.nm = Mv(o), Mv(x), Mv(y), Mv(n)
        self.invm = [Mv(row) for row in inv]
        a, b, c = cross(y, n), cross(n, x), cross(x, y)
        self.orthonormal = (dot(a, b) == 0 and dot(a, c) == 0 and dot(b, c) == 0
                            and dot(a, a) == dot(b, b) == dot(c, c))

    def world(self, u, v, w):
        o, x, y, n = self.om, self.xm, self.ym, self.nm
        return (o[0]+u*x[0]+v*y[0]+w*n[0], o[1]+u*x[1]+v*y[1]+w*n[1], o[2]+u*x[2]+v*y[2]+w*n[2])

    def vec(self, u, v, w):
        x, y, n = self.xm, self.ym, self.nm
        return (u*x[0]+v*y[0]+w*n[0], u*x[1]+v*y[1]+w*n[1], u*x[2]+v*y[2]+w*n[2])

    def chart(self, X):
        d = (X[0]-self.om[0], X[1]-self.om[1], X[2]-self.om[2])
        return tuple(r[0]*d[0]+r[1]*d[1]+r[2]*d[2] for r in self.invm)

    def covector(self, g):
        """World gradient from a chart gradient."""
        i = self.invm
        return tuple(g[0]*i[0][k]+g[1]*i[1][k]+g[2]*i[2][k] for k in range(3))


# ------------------------------------------------------------------ surfaces

class Surf:
    """A surface of an input: `f(X)` (world, mpf), its degree and gradient."""

    def __init__(self, label, degree, f, grad):
        self.label, self.degree, self.f, self.grad = label, degree, f, grad


def plane_surf(label, fr, a, b):
    """The chart plane `a . p = b` as `f = a . p - b` (world)."""
    am = Mv(a)
    bm = M(b)
    g = fr.covector(am)
    return Surf(label, 1, lambda X: dot(am, fr.chart(X))-bm, lambda X: g)


# ------------------------------------------------------------------ curves and families

class Curve:
    """One curve of a family at outer parameter `a`: a circle `X0 + C cos b + S
    sin b` (`kind` 'circle', `rng` None for the whole circle or `(lo, hi)`)
    or a line `X0 + b D` over `rng`; `dX0`, `dC`, `dS` (`dD`) the derivatives
    in `a`."""

    def __init__(self, kind, X0, C, S, dX0, dC, dS, rng):
        self.kind, self.X0, self.C, self.S = kind, X0, C, S
        self.dX0, self.dC, self.dS, self.rng = dX0, dC, dS, rng

    def point(self, b):
        if self.kind == 'line':
            return lin((1, self.X0), (b, self.C))
        return lin((1, self.X0), (mp.cos(b), self.C), (mp.sin(b), self.S))

    def frame_at(self, b, c=None, s=None):
        """(X, X_a, X_b) at `b`."""
        if self.kind == 'line':
            X = lin((1, self.X0), (b, self.C))
            Xa = lin((1, self.dX0), (b, self.dC))
            return X, Xa, self.C
        c = mp.cos(b) if c is None else c
        s = mp.sin(b) if s is None else s
        X = lin((1, self.X0), (c, self.C), (s, self.S))
        Xa = lin((1, self.dX0), (c, self.dC), (s, self.dS))
        Xb = lin((-s, self.C), (c, self.S))
        return X, Xa, Xb

    def degenerate(self):
        """A circle of radius zero or a line of no length."""
        if self.kind == 'line':
            return self.rng[1] <= self.rng[0]
        return dot(self.C, self.C)+dot(self.S, self.S) == 0


class Family:
    """A family of curves covering a face: outer range `[a0, a1]`
    (`periodic`), `curve(a)`, the orientation `sign` (`sign * X_a x X_b`
    outward) and the family's own kinks."""

    def __init__(self, name, a0, a1, periodic, curve, sign, own=()):
        self.name, self.a0, self.a1, self.periodic = name, a0, a1, periodic
        self.curve, self.sign, self.own = curve, sign, tuple(own)


class Face:
    def __init__(self, name, families):
        self.name, self.families = name, families


def circle_curve(X0, C, S, dX0, dC, dS, rng=None):
    return Curve('circle', X0, C, S, dX0, dC, dS, rng)


def line_curve(X0, D, dX0, dD, b0, b1):
    return Curve('line', X0, D, None, dX0, dD, None, (b0, b1))


# ------------------------------------------------------------------ inputs

class Torus:
    """A whole torus (S9d.4a's model)."""
    kind = 'torus'

    def __init__(self, case):
        R, r, low, high, angle = case.torus
        assert high-low == TWO_PI and angle == TWO_PI, 'S9d.4b.2: a whole torus'
        self.case = case
        self.fr = fr = Frame(case.frame)
        self.R, self.r = F(R), F(r)
        assert self.R > self.r > 0
        self.Rm, self.rm = M(self.R), M(self.r)
        R2, k = self.Rm**2, self.Rm**2-self.rm**2

        def f(X):
            p = fr.chart(X)
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]
            return (s+k)**2-4*R2*(p[0]*p[0]+p[1]*p[1])

        def grad(X):
            p = fr.chart(X)
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]
            g = 4*(s+k)
            return fr.covector((g*p[0]-8*R2*p[0], g*p[1]-8*R2*p[1], g*p[2]))
        self.surfs = [Surf('torus', 4, f, grad)]
        self.faces = [self.wall()]

    def contains(self, X):
        return self.surfs[0].f(X) < 0

    def e(self, t):
        fr = self.fr
        c, s = mp.cos(t), mp.sin(t)
        return lin((c, fr.xm), (s, fr.ym)), lin((-s, fr.xm), (c, fr.ym))

    def wall(self):
        fr, R, r = self.fr, self.Rm, self.rm
        zero = (Z, Z, Z)

        def tube(th):
            e, de = self.e(th)
            return circle_curve(lin((1, fr.om), (R, e)), scale(e, r), scale(fr.nm, r),
                                scale(de, R), scale(de, r), zero)

        def parallel(ph):
            c, s = mp.cos(ph), mp.sin(ph)
            rho = R+r*c
            return circle_curve(lin((1, fr.om), (r*s, fr.nm)), scale(fr.xm, rho), scale(fr.ym, rho),
                                scale(fr.nm, r*c), scale(fr.xm, -r*s), scale(fr.ym, -r*s))
        return Face('wall', [Family('tube circles', Z, 2*mp.pi, True, tube, 1),
                             Family('parallels', Z, 2*mp.pi, True, parallel, -1)])

    def closed(self):
        """Volume, world moments and area in closed form (Pappus; the area of
        an affine frame's wall by the periodic trapezoidal rule on its own
        factor, exponentially accurate)."""
        R, r = self.Rm, self.rm
        det = M(self.fr.det)
        V = 2*mp.pi**2*R*r*r*det
        mom = tuple(c*V for c in self.fr.om)
        if self.fr.orthonormal:
            A = 4*mp.pi**2*R*r*norm(Mv(cross(self.fr.x, self.fr.y)))
        else:
            m = 64
            h = 2*mp.pi/m
            A = Z
            for i in range(m):
                for j in range(m):
                    cur = self.faces[0].families[0].curve(i*h)
                    _, Xa, Xb = cur.frame_at(j*h)
                    A += norm(cross(Xa, Xb))
            A *= h*h
        return V, mom, A

    def box(self):
        R, r = float(self.R), float(self.r)
        return frame_box(self.fr, [(u, v, w) for u in (-R-r, R+r) for v in (-R-r, R+r) for w in (-r, r)])

    def ray_intervals_chart(self, s):
        """The annulus at chart level `s` (float) as radii."""
        R, r = float(self.R), float(self.r)
        if abs(s) >= r:
            return None
        q = math.sqrt(r*r-s*s)
        return (R-q, R+q)


class Sphere:
    """A whole sphere (S9d.1's model: `|X - o|^2 <= R^2` in the world)."""
    kind = 'sphere'

    def __init__(self, case):
        radius, low, high = case.sphere
        assert low == -HALF_PI and high == HALF_PI, 'S9d.4b.2: a whole sphere'
        self.case = case
        self.fr = Frame(case.frame)
        self.c, self.R = self.fr.o, F(radius)
        self.cm, self.Rm = Mv(self.c), M(self.R)
        cm, R2 = self.cm, self.Rm**2
        self.surfs = [Surf('sphere', 2, lambda X: (X[0]-cm[0])**2+(X[1]-cm[1])**2+(X[2]-cm[2])**2-R2,
                           lambda X: (2*(X[0]-cm[0]), 2*(X[1]-cm[1]), 2*(X[2]-cm[2])))]
        self.faces = [self.wall()]

    def contains(self, X):
        return self.surfs[0].f(X) < 0

    def wall(self):
        c, R = self.cm, self.Rm
        e1, e2, e3 = (mp.mpf(1), Z, Z), (Z, mp.mpf(1), Z), (Z, Z, mp.mpf(1))
        zero = (Z, Z, Z)

        def meridian(lam):
            cl, sl = mp.cos(lam), mp.sin(lam)
            return circle_curve(c, lin((R*cl, e1), (R*sl, e2)), scale(e3, R), zero, lin((-R*sl, e1), (R*cl, e2)),
                                zero, (-mp.pi/2, mp.pi/2))

        def parallel(beta):
            cb, sb = mp.cos(beta), mp.sin(beta)
            return circle_curve(lin((1, c), (R*sb, e3)), scale(e1, R*cb), scale(e2, R*cb), scale(e3, R*cb),
                                scale(e1, -R*sb), scale(e2, -R*sb))
        return Face('wall', [Family('meridians', Z, 2*mp.pi, True, meridian, 1),
                             Family('parallels', -mp.pi/2, mp.pi/2, False, parallel, -1)])

    def closed(self):
        R = self.Rm
        V = 4*mp.pi*R**3/3
        return V, tuple(c*V for c in self.cm), 4*mp.pi*R*R

    def box(self):
        c, R = [float(v) for v in self.c], float(self.R)
        return [v-R for v in c], [v+R for v in c]


def frame_box(fr, pts):
    ws = [[float(v) for v in fr.world(M(u), M(v), M(w))] for u, v, w in pts]
    return [min(p[i] for p in ws) for i in range(3)], [max(p[i] for p in ws) for i in range(3)]


class Profile:
    """A convex profile of segments and counter-clockwise arcs, or a circle
    (S9d.2's, exact): its elements `('seg', p, q)` or `('arc', C, a, t0, t1)`
    (`t1 > t0`, a circle's `t1 = t0 + 2 pi`), counter-clockwise."""

    def __init__(self, boundary, tol):
        pts, _ = stored(boundary, tol)
        self.els = []
        self.vertices = []
        if pts is None:
            cx, cy, a = boundary.circle
            self.els.append(('arc', (F(cx), F(cy)), F(a), Z, 2*mp.pi))
            self.circle = True
        else:
            self.circle = False
            if boundary.segments is None:
                points, segs = pts, [None]*len(pts)
            else:
                points, segs = pts
            ring = [(F(u), F(v)) for u, v in points]
            self.vertices = ring
            n = len(ring)
            for i in range(n):
                p, q, sg = ring[i], ring[(i+1) % n], segs[i]
                if sg is None:
                    self.els.append(('seg', p, q))
                else:
                    cx, cy, a, ccw = sg
                    assert ccw, 'S9d.4b.2: a convex profile (counter-clockwise arcs)'
                    C, a = (F(cx), F(cy)), F(a)
                    for pt in (p, q):
                        assert (pt[0]-C[0])**2+(pt[1]-C[1])**2 == a*a, 'an arc end off its circle'
                    t0 = mp.atan2(M(p[1]-C[1]), M(p[0]-C[0]))
                    t1 = mp.atan2(M(q[1]-C[1]), M(q[0]-C[0]))
                    t1 = t0+(t1-t0) % (2*mp.pi)
                    self.els.append(('arc', C, a, t0, t1))
            # Convex: every vertex turns left (exactly for segments; arcs
            # bulge outward).
            for i in range(n):
                a, b, c = ring[i-1], ring[i], ring[(i+1) % n]
                assert (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]) >= 0, 'S9d.4b.2: a convex profile'

    def contains(self, u, v):
        """Strictly inside (mpf): the winding of the boundary about the point
        (convex: left of every segment, and inside every arc's circle on the
        arc's side of its chord)."""
        for e in self.els:
            if e[0] == 'seg':
                p, q = Mv(e[1]), Mv(e[2])
                if (q[0]-p[0])*(v-p[1])-(q[1]-p[1])*(u-p[0]) <= 0:
                    return False
            else:
                _, C, a, t0, t1 = e
                Cm, am = Mv(C), M(a)
                if t1-t0 >= 2*mp.pi-mp.mpf(10)**-30:
                    if (u-Cm[0])**2+(v-Cm[1])**2 >= am*am:
                        return False
                    continue
                p = (Cm[0]+am*mp.cos(t0), Cm[1]+am*mp.sin(t0))
                q = (Cm[0]+am*mp.cos(t1), Cm[1]+am*mp.sin(t1))
                side = (q[0]-p[0])*(v-p[1])-(q[1]-p[1])*(u-p[0])
                if side <= 0 and (u-Cm[0])**2+(v-Cm[1])**2 >= am*am:
                    return False
        return True

    def chord(self, P, d):
        """The parameters where the line `P + t d` (profile coordinates, mpf)
        is inside the convex profile: `(t0, t1)` or None."""
        ts = []
        for e in self.els:
            if e[0] == 'seg':
                p, q = Mv(e[1]), Mv(e[2])
                e_ = (q[0]-p[0], q[1]-p[1])
                den = d[0]*e_[1]-d[1]*e_[0]
                if den == 0:
                    continue
                w = (p[0]-P[0], p[1]-P[1])
                t = (w[0]*e_[1]-w[1]*e_[0])/den
                s = (w[0]*d[1]-w[1]*d[0])/den
                if 0 <= s <= 1:
                    ts.append(t)
            else:
                _, C, a, t0, t1 = e
                Cm, am = Mv(C), M(a)
                w = (P[0]-Cm[0], P[1]-Cm[1])
                A = d[0]*d[0]+d[1]*d[1]
                B = 2*(w[0]*d[0]+w[1]*d[1])
                Cc = w[0]*w[0]+w[1]*w[1]-am*am
                D = B*B-4*A*Cc
                if D <= 0:
                    continue
                sq = mp.sqrt(D)
                for t in ((-B-sq)/(2*A), (-B+sq)/(2*A)):
                    ang = mp.atan2(w[1]+t*d[1], w[0]+t*d[0])
                    if t1-t0 >= 2*mp.pi-mp.mpf(10)**-30 or ((ang-t0) % (2*mp.pi)) <= t1-t0:
                        ts.append(t)
        if len(ts) < 2:
            return None
        lo, hi = min(ts), max(ts)
        return (lo, hi) if hi > lo else None

    def support(self, d):
        """`(min, max)` of `d . p` over the profile, and the kinks of a chord
        family across `d` (vertices' offsets)."""
        vals = [M(p[0])*d[0]+M(p[1])*d[1] for p in self.vertices]
        kinks = list(vals)
        for e in self.els:
            if e[0] == 'arc':
                _, C, a, t0, t1 = e
                Cm, am = Mv(C), M(a)
                phi = mp.atan2(d[1], d[0])
                dn = norm((d[0], d[1], Z))
                for ang in (phi, phi+mp.pi):
                    if t1-t0 >= 2*mp.pi-mp.mpf(10)**-30 or ((ang-t0) % (2*mp.pi)) <= t1-t0:
                        vals.append(Cm[0]*d[0]+Cm[1]*d[1]+am*dn*mp.cos(ang-phi))
        return min(vals), max(vals), kinks


class Prism:
    """A prism of a convex profile of segments and arcs, or a circle (S9d.2's
    model)."""
    kind = 'prism'

    def __init__(self, case):
        assert len(case.boundaries) == 1, 'S9d.4b.2: a prism of one convex profile'
        self.case = case
        self.fr = fr = Frame(case.frame)
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        self.lom, self.him = M(self.lo), M(self.hi)
        self.profile = Profile(case.boundaries[0], case.tolerance)
        self.surfs = []
        for k, e in enumerate(self.profile.els):
            if e[0] == 'seg':
                p, q = e[1], e[2]
                a = (q[1]-p[1], -(q[0]-p[0]), F(0))
                self.surfs.append(plane_surf(('wall', k), fr, a, a[0]*p[0]+a[1]*p[1]))
            else:
                self.surfs.append(self.cylinder_surf(('wall', k), e[1], e[2]))
        self.surfs.append(plane_surf('start', fr, (F(0), F(0), F(-1)), -self.lo))
        self.surfs.append(plane_surf('end', fr, (F(0), F(0), F(1)), self.hi))
        self.faces = self.make_faces()

    def cylinder_surf(self, label, C, a):
        fr = self.fr
        Cm, a2 = Mv(C), M(a)**2

        def f(X):
            p = fr.chart(X)
            return (p[0]-Cm[0])**2+(p[1]-Cm[1])**2-a2

        def grad(X):
            p = fr.chart(X)
            return fr.covector((2*(p[0]-Cm[0]), 2*(p[1]-Cm[1]), Z))
        return Surf(label, 2, f, grad)

    def contains(self, X):
        p = self.fr.chart(X)
        return self.lom < p[2] < self.him and self.profile.contains(p[0], p[1])

    def make_faces(self):
        fr = self.fr
        lo, hi = self.lom, self.him
        zero = (Z, Z, Z)
        faces = []
        for k, e in enumerate(self.profile.els):
            if e[0] == 'seg':
                p, q = Mv(e[1]), Mv(e[2])
                d = fr.vec(q[0]-p[0], q[1]-p[1], Z)

                def gen(s, p=p, d=d):
                    return line_curve(lin((1, fr.world(p[0], p[1], Z)), (s, d)), fr.nm, d, zero, lo, hi)

                def row(h, p=p, d=d):
                    return line_curve(fr.world(p[0], p[1], h), d, fr.nm, zero, Z, mp.mpf(1))
                faces.append(Face(('wall', k), [Family('generatrices', Z, mp.mpf(1), False, gen, 1),
                                                Family('rows', lo, hi, False, row, -1)]))
            else:
                _, C, a, t0, t1 = e
                Cm, am = Mv(C), M(a)
                full = t1-t0 >= 2*mp.pi-mp.mpf(10)**-30

                def gen(t, Cm=Cm, am=am):
                    c, s = mp.cos(t), mp.sin(t)
                    return line_curve(fr.world(Cm[0]+am*c, Cm[1]+am*s, Z), fr.nm, fr.vec(-am*s, am*c, Z), zero,
                                      lo, hi)

                def section(h, Cm=Cm, am=am, t0=t0, t1=t1, full=full):
                    return circle_curve(fr.world(Cm[0], Cm[1], h), fr.vec(am, Z, Z), fr.vec(Z, am, Z),
                                        fr.nm, zero, zero, None if full else (t0, t1))
                faces.append(Face(('wall', k), [Family('generatrices', t0, t1, full, gen, 1),
                                                Family('sections', lo, hi, False, section, -1)]))
        for h, name, out in ((lo, 'start', -1), (hi, 'end', 1)):
            faces.append(Face(name, [self.chords(h, (F(1), F(0)), out, 'chords x'),
                                     self.chords(h, (F(2), F(7)), out, 'chords skew')]))
        return faces

    def chords(self, h, direction, out, name):
        """Chords of the cap at height `h` along the profile direction `d`,
        the outer parameter the offset along `d^perp`."""
        fr = self.fr
        d = Mv(direction)
        dn = norm((d[0], d[1], Z))
        d = (d[0]/dn, d[1]/dn)
        nperp = (-d[1], d[0])
        lo, hi, kinks = self.profile.support(nperp)
        D = fr.vec(d[0], d[1], Z)
        dP = fr.vec(nperp[0], nperp[1], Z)
        zero = (Z, Z, Z)

        def chord(t):
            P = (t*nperp[0], t*nperp[1])
            iv = self.profile.chord(P, d)
            b0, b1 = iv if iv is not None else (Z, Z)
            return line_curve(fr.world(P[0], P[1], h), D, dP, zero, b0, b1)
        # X_t x X_b = (nperp x d) n = -n: outward `out n` needs sign -out.
        own = [k for k in kinks if lo < k < hi]
        return Family(name, lo, hi, False, chord, -out, own)

    def closed(self):
        """Volume, world moments and area in closed form: the profile's area
        and moments by Green's theorem on its segments and arcs, the walls'
        and caps' areas through the frame's factors (the walls' by 24-point
        Gauss-Legendre on each arc's factor where the frame is not
        orthonormal)."""
        fr = self.fr
        A0 = mx = my = per = Z
        wall = Z
        for e in self.profile.els:
            if e[0] == 'seg':
                p, q = Mv(e[1]), Mv(e[2])
                cr = p[0]*q[1]-q[0]*p[1]
                A0 += cr/2
                mx += (p[0]+q[0])*cr/6
                my += (p[1]+q[1])*cr/6
                wall += norm(cross(fr.vec(q[0]-p[0], q[1]-p[1], Z), fr.nm))
            else:
                _, C, a, t0, t1 = e
                Cm, am = Mv(C), M(a)
                # Green: A = 1/2 int (x dy - y dx), moments int x^2 dy / 2, -int y^2 dx / 2.
                th = t1-t0
                s0, s1, c0, c1 = mp.sin(t0), mp.sin(t1), mp.cos(t0), mp.cos(t1)
                A0 += (am*am*th+am*(Cm[0]*(s1-s0)-Cm[1]*(c1-c0)))/2
                # int x^2 dy with x = Cx + a c, dy = a c dt.
                ix = Cm[0]**2*am*(s1-s0)+2*Cm[0]*am*am*(th/2+(mp.sin(2*t1)-mp.sin(2*t0))/4) \
                    + am**3*((s1-s1**3/3)-(s0-s0**3/3))
                iy = Cm[1]**2*am*(c1-c0)*(-1)+2*Cm[1]*am*am*(th/2-(mp.sin(2*t1)-mp.sin(2*t0))/4) \
                    + am**3*((-c1+c1**3/3)-(-c0+c0**3/3))
                mx += ix/2
                my += iy/2
                if fr.orthonormal:
                    wall += am*th*norm(Mv(cross(fr.x, fr.y)))
                else:
                    pieces = max(1, int(mp.ceil(th/quarter())))
                    step = th/pieces
                    for i in range(pieces):
                        for xi, wi in nodes(INNER):
                            t = t0+step*i+step*(xi+1)/2
                            wall += wi*step/2*norm(cross(fr.vec(-am*mp.sin(t), am*mp.cos(t), Z), fr.nm))
        h = self.him-self.lom
        det = M(fr.det)
        V = A0*h*det
        cu, cv, cw = mx/A0, my/A0, (self.lom+self.him)/2
        mom = tuple(c*V for c in fr.world(cu, cv, cw))
        cap = A0*norm(Mv(cross(fr.x, fr.y)))
        return V, mom, wall*h+2*cap

    def box(self):
        us, vs = [], []
        for e in self.profile.els:
            if e[0] == 'seg':
                us += [float(e[1][0])]
                vs += [float(e[1][1])]
            else:
                C, a = e[1], e[2]
                us += [float(C[0]-a), float(C[0]+a)]
                vs += [float(C[1]-a), float(C[1]+a)]
        return frame_box(self.fr, [(u, v, w) for u in (min(us), max(us)) for v in (min(vs), max(vs))
                                   for w in (self.lo, self.hi)])


class Cone:
    """A cone or frustum (S9d.3a's model)."""
    kind = 'cone'

    def __init__(self, case):
        self.case = case
        self.fr = fr = Frame(case.frame)
        r0, r1, h = (F(v) for v in case.cone)
        assert h > 0 and r0 >= 0 and r1 >= 0 and r0 != r1
        self.r0, self.r1, self.h = r0, r1, h
        self.k = (r1-r0)/h
        r0m, km = M(r0), M(self.k)
        self.r0m, self.r1m, self.hm, self.km = r0m, M(r1), M(h), km

        def f(X):
            p = fr.chart(X)
            return p[0]*p[0]+p[1]*p[1]-(r0m+km*p[2])**2

        def grad(X):
            p = fr.chart(X)
            return fr.covector((2*p[0], 2*p[1], -2*km*(r0m+km*p[2])))
        self.surfs = [Surf('cone', 2, f, grad), plane_surf('start', fr, (F(0), F(0), F(-1)), F(0)),
                      plane_surf('end', fr, (F(0), F(0), F(1)), h)]
        self.faces = self.make_faces()

    def contains(self, X):
        p = self.fr.chart(X)
        return 0 < p[2] < self.hm and p[0]*p[0]+p[1]*p[1] < (self.r0m+self.km*p[2])**2

    def apex(self):
        if self.r1 == 0:
            return self.fr.world(Z, Z, self.hm)
        if self.r0 == 0:
            return self.fr.world(Z, Z, Z)
        return None

    def make_faces(self):
        fr = self.fr
        r0, k, h = self.r0m, self.km, self.hm
        zero = (Z, Z, Z)

        def ruling(t):
            c, s = mp.cos(t), mp.sin(t)
            e, de = fr.vec(c, s, Z), fr.vec(-s, c, Z)
            return line_curve(lin((1, fr.om), (r0, e)), lin((k, e), (1, fr.nm)), scale(de, r0), scale(de, k), Z, h)

        def parallel(w):
            rho = r0+k*w
            return circle_curve(fr.world(Z, Z, w), fr.vec(rho, Z, Z), fr.vec(Z, rho, Z), fr.nm, fr.vec(k, Z, Z),
                                fr.vec(Z, k, Z))
        faces = [Face('wall', [Family('rulings', Z, 2*mp.pi, True, ruling, 1),
                               Family('parallels', Z, h, False, parallel, -1)])]
        for w, rad, name, out in ((Z, self.r0m, 'start', -1), (h, self.r1m, 'end', 1)):
            if rad > 0:
                faces.append(Face(name, [disc_chords(fr, w, rad, (F(1), F(0)), out, 'chords x'),
                                         disc_chords(fr, w, rad, (F(2), F(7)), out, 'chords skew')]))
        return faces

    def closed(self):
        """Volume, world moments and area in closed form (S9d.3a's)."""
        r0, r1, h = self.r0m, self.r1m, self.hm
        fr = self.fr
        det = M(fr.det)
        V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3*det
        zc = h*(r0*r0+2*r0*r1+3*r1*r1)/(4*(r0*r0+r0*r1+r1*r1))
        mom = tuple(c*V for c in fr.world(Z, Z, zc))
        cap = norm(Mv(cross(fr.x, fr.y)))
        area = mp.pi*(r0*r0+r1*r1)*cap
        k = self.km
        if fr.orthonormal:
            area += mp.pi*(r0+r1)*mp.sqrt(h*h+(r1-r0)**2)*cap
        else:
            # The wall: (r0 + k w) |U'(t) x (k U + n)| over t and w.
            tot = Z
            pieces = 8
            step = 2*mp.pi/pieces
            for i in range(pieces):
                for xi, wi in nodes(INNER):
                    t = step*i+step*(xi+1)/2
                    c, s = mp.cos(t), mp.sin(t)
                    tot += wi*step/2*norm(cross(fr.vec(-s, c, Z), lin((k, fr.vec(c, s, Z)), (1, fr.nm))))
            area += tot*h*(r0+r1)/2
        return V, mom, area

    def box(self):
        R = float(max(self.r0, self.r1))
        return frame_box(self.fr, [(u, v, w) for u in (-R, R) for v in (-R, R) for w in (0, float(self.h))])


def disc_chords(fr, w, rad, direction, out, name):
    """Chords of the disc of radius `rad` about the axis at height `w`."""
    d = Mv(direction)
    dn = norm((d[0], d[1], Z))
    d = (d[0]/dn, d[1]/dn)
    nperp = (-d[1], d[0])
    D = fr.vec(d[0], d[1], Z)
    dP = fr.vec(nperp[0], nperp[1], Z)
    zero = (Z, Z, Z)

    def chord(t):
        q = mp.sqrt(max(rad*rad-t*t, Z))
        return line_curve(fr.world(t*nperp[0], t*nperp[1], w), D, dP, zero, -q, q)
    return Family(name, -rad, rad, False, chord, -out)


def make_input(case):
    if case.torus is not None:
        return Torus(case)
    if case.sphere is not None:
        return Sphere(case)
    if case.cone is not None:
        return Cone(case)
    return Prism(case)


# ------------------------------------------------------------------ roots along a curve

_TABLES = {}


def dft_table(N):
    key = (N, mp.mp.prec)
    if key not in _TABLES:
        bs = [2*mp.pi*j/N for j in range(N)]
        _TABLES[key] = (bs, [mp.cos(b) for b in bs], [mp.sin(b) for b in bs],
                        [[mp.cos(k*b) for b in bs] for k in range(N)], [[mp.sin(k*b) for b in bs] for k in range(N)])
    return _TABLES[key]


def trig_coeffs(cur, surf):
    """The trigonometric polynomial `f(X(b)) = a0 + sum a_k cos kb + b_k sin
    kb` of a surface along a circle, exact for its degree `d` from `2d + 1`
    samples."""
    d = surf.degree
    N = 2*d+1
    bs, cs, ss, ck, sk = dft_table(N)
    vals = [surf.f(lin((1, cur.X0), (cs[j], cur.C), (ss[j], cur.S))) for j in range(N)]
    a = [sum(vals)/N]
    b = [Z]
    for k in range(1, d+1):
        a.append(2*sum(v*c for v, c in zip(vals, ck[k]))/N)
        b.append(2*sum(v*s for v, s in zip(vals, sk[k]))/N)
    return a, b


def trig_value(a, b, t):
    out = a[0]
    for k in range(1, len(a)):
        out += a[k]*mp.cos(k*t)+b[k]*mp.sin(k*t)
    return out


def trig_deriv(a, b, t):
    out = Z
    for k in range(1, len(a)):
        out += k*(b[k]*mp.cos(k*t)-a[k]*mp.sin(k*t))
    return out


def float_roots(desc):
    """Complex roots of a polynomial (descending complex coefficients) in
    binary64 by the Aberth-Ehrlich iteration; None when it does not settle."""
    n = len(desc)-1
    lead = desc[0]
    c = [x/lead for x in desc]
    rad = 1+max(abs(x) for x in c[1:])
    zs = [rad*0.5*cmath.exp(1j*(2*math.pi*k/n+0.4)) for k in range(n)]
    for _ in range(200):
        done = True
        for i in range(n):
            z = zs[i]
            p = dp = 0
            for x in c:
                dp = dp*z+p
                p = p*z+x
            if p == 0:
                continue
            ratio = p/dp if dp != 0 else 1e-3
            s = sum(1/(z-zs[j]) for j in range(n) if j != i and z != zs[j])
            step = ratio/(1-ratio*s)
            zs[i] = z-step
            if abs(step) > 1e-15*(1+abs(z)):
                done = False
        if done:
            return zs
    return None


def mp_roots(desc):
    """Complex roots (mpc) of a polynomial with mp coefficients (descending):
    floating-point Aberth roots polished by Newton's method in the working
    precision, or `mp.polyroots` where roots cluster or polishing fails."""
    while len(desc) > 1 and desc[0] == 0:
        desc = desc[1:]
    n = len(desc)-1
    if n <= 0:
        return []
    fz = None
    try:
        fz = float_roots([complex(x) for x in desc])
    except (OverflowError, ZeroDivisionError):
        fz = None
    if fz is not None:
        close = any(abs(fz[i]-fz[j]) < 1e-5*(1+abs(fz[i])) for i in range(n) for j in range(i))
        if not close:
            out = []
            ok = True
            for z0 in fz:
                z = mp.mpc(z0)
                for _ in range(12):
                    p = dp = mp.mpc(0)
                    for x in desc:
                        dp = dp*z+p
                        p = p*z+x
                    if dp == 0:
                        ok = False
                        break
                    step = p/dp
                    z -= step
                    if abs(step) <= mp.mpf(10)**(-mp.mp.dps+2)*(1+abs(z)):
                        break
                else:
                    ok = False
                if not ok:
                    break
                out.append(z)
            if ok and not any(abs(out[i]-out[j]) < mp.mpf(10)**-8*(1+abs(out[i])) for i in range(n) for j in range(i)):
                return out
    try:
        return list(mp.polyroots(desc, maxsteps=400, extraprec=3*mp.mp.prec))
    except mp.NoConvergence:
        return list(mp.polyroots(desc, maxsteps=4000, extraprec=8*mp.mp.prec, error=False))


REAL_TOL = mp.mpf(10)**-15
GOLDEN = (3-mp.sqrt(5))/2


def circle_roots(a, b):
    """The real roots in `[0, 2 pi)` of a trigonometric polynomial: the unit
    circle's roots of `z^d f(z)`, a root within 1e-15 of the circle taken as
    real, polished by Newton's method on `f`."""
    d = len(a)-1
    scale_ = max([abs(x) for x in a]+[abs(x) for x in b])
    if scale_ == 0:
        return []
    while d > 0 and abs(a[d])+abs(b[d]) <= mp.mpf(10)**-34*scale_:
        d -= 1
    if d == 0:
        return []
    # c_k = (a_k - i b_k)/2, c_-k = conj; z^d f = sum c_{m-d} z^m.
    coef = []
    for m in range(2*d+1):
        k = m-d
        if k == 0:
            coef.append(mp.mpc(a[0]))
        elif k > 0:
            coef.append(mp.mpc(a[k], -b[k])/2)
        else:
            coef.append(mp.mpc(a[-k], b[-k])/2)
    zs = mp_roots(list(reversed(coef)))
    out = []
    T = 2*mp.pi
    for z in zs:
        if z == 0:
            continue
        if abs(mp.log(abs(z))) > REAL_TOL:
            continue
        t = mp.arg(z) % T
        for _ in range(6):
            dv = trig_deriv(a, b, t)
            if dv == 0:
                break
            step = trig_value(a, b, t)/dv
            if abs(step) > mp.mpf(10)**-12:
                break
            t -= step
            if abs(step) < mp.mpf(10)**(-mp.mp.dps+3):
                break
        out.append(t % T)
    return sorted(out)


def line_coeffs(cur, surf):
    """`f(X0 + b D)` as ascending coefficients in `tau`, `b = m + h tau` over
    the line's range, from `d + 1` Chebyshev samples."""
    d = surf.degree
    b0, b1 = cur.rng
    m, h = (b0+b1)/2, (b1-b0)/2
    taus = [mp.cos(mp.pi*(2*j+1)/(2*(d+1))) for j in range(d+1)]
    vals = [surf.f(lin((1, cur.X0), (m+h*t, cur.C))) for t in taus]
    A = [[t**k for k in range(d+1)] for t in taus]
    sol = mp.lu_solve(mp.matrix(A), mp.matrix(vals))
    return [sol[i] for i in range(d+1)], m, h


def line_roots(cur, surf):
    """The real roots of a surface along a line strictly inside its range."""
    c, m, h = line_coeffs(cur, surf)
    scale_ = max(abs(x) for x in c)
    if scale_ == 0:
        return []
    while len(c) > 1 and abs(c[-1]) <= mp.mpf(10)**-34*scale_:
        c.pop()
    if len(c) <= 1:
        return []
    out = []
    for z in mp_roots(list(reversed(c))):
        if abs(mp.im(z)) > REAL_TOL*(1+abs(z)):
            continue
        t = mp.re(z)
        if -1 < t < 1:
            for _ in range(6):
                p = dp = Z
                for x in reversed(c):
                    dp = dp*t+p
                    p = p*t+x
                if dp == 0:
                    break
                step = p/dp
                if abs(step) > mp.mpf(10)**-12:
                    break
                t -= step
                if abs(step) < mp.mpf(10)**(-mp.mp.dps+3):
                    break
            if -1 < t < 1:
                out.append(m+h*t)
    return sorted(out)


def curve_roots(cur, surf):
    """Real roots of a surface along a curve, inside the curve's range."""
    if cur.kind == 'line':
        return line_roots(cur, surf)
    a, b = trig_coeffs(cur, surf)
    rts = circle_roots(a, b)
    if cur.rng is None:
        return rts
    lo, hi = cur.rng
    T = 2*mp.pi
    out = []
    for t in rts:
        t = lo+(t-lo) % T
        if lo < t < hi:
            out.append(t)
    return sorted(out)


# ------------------------------------------------------------------ pieces of a curve

def curve_pieces(cur, other):
    """The curve cut at every root of the other input's surfaces, each piece
    `(b0, b1, class, label0, label1)` classified inside it."""
    if cur.degenerate():
        return []
    roots = []
    closed = cur.kind == 'circle' and cur.rng is None
    span = 2*mp.pi if closed else cur.rng[1]-cur.rng[0]
    for j, s in enumerate(other.surfs):
        rts = curve_roots(cur, s)
        # A double root (a tangency) split by rounding: both dropped, no
        # sign change between them; a pair of simple roots closer than
        # 1e-15 of the range lies within 1e-30 of an event, where no
        # quadrature node falls.
        keep = [True]*len(rts)
        n = len(rts)
        for i in range(n):
            k = (i+1) % n
            if (k == 0 and not closed) or n < 2 or not keep[i] or not keep[k]:
                continue
            gap = rts[k]-rts[i]+(span if k == 0 else 0)
            if gap < mp.mpf(10)**-15*span:
                keep[i] = keep[k] = False
        roots += [(t, j) for t, ok in zip(rts, keep) if ok]
    roots.sort(key=lambda x: x[0])
    # Each piece classified at a point a golden fraction along it, not its
    # midpoint: a symmetric tangency (dropped above) may sit at the middle.
    cls = lambda t0, t1: 'in' if other.contains(cur.point(t0+(t1-t0)*GOLDEN)) else 'out'
    out = []
    if cur.kind == 'circle' and cur.rng is None:
        T = 2*mp.pi
        if not roots:
            return [(Z, T, cls(Z, T), None, None)]
        for i, (t, j) in enumerate(roots):
            t1, j1 = roots[(i+1) % len(roots)]
            if i+1 == len(roots):
                t1 += T
            if t1 > t:
                out.append((t, t1, cls(t, t1), j, j1))
        return out
    lo, hi = cur.rng
    pts = [(lo, 'end')]+[r for r in roots if lo < r[0] < hi]+[(hi, 'end')]
    for (t0, j0), (t1, j1) in zip(pts, pts[1:]):
        if t1 > t0:
            out.append((t0, t1, cls(t0, t1), j0, j1))
    return out


def signature(cur, pieces):
    """The curve's structure: its pieces' classes (adjacent equal classes
    merged) with the surface bounding each, cyclic on a closed curve."""
    if not pieces:
        return ('none',)
    merged = []
    for p in pieces:
        if merged and merged[-1][2] == p[2]:
            merged[-1] = (merged[-1][0], p[1], p[2], merged[-1][3], p[4])
        else:
            merged.append(p)
    closed = cur.kind == 'circle' and cur.rng is None
    if closed and len(merged) > 1 and merged[0][2] == merged[-1][2]:
        last = merged.pop()
        merged[0] = (last[0], merged[0][1], last[2], last[3], merged[0][4])
    seq = tuple((p[2], p[3]) for p in merged)
    if not closed:
        return seq
    if len(seq) == 1:
        return (seq[0][0],)
    rots = [seq[i:]+seq[:i] for i in range(len(seq))]
    return min(rots, key=repr)


def normal_jet(cur, sign):
    """`sign * X_a x X_b` along the curve: on a circle a trigonometric
    polynomial of degree two (`N0 + Nc cos b + Ns sin b + N2c cos 2b + N2s
    sin 2b`), on a line `N0 + b N1`."""
    k = 1 if sign > 0 else -1
    sc = lambda v, f: (v[0]*f, v[1]*f, v[2]*f)
    if cur.kind == 'line':
        return (sc(cross(cur.dX0, cur.C), k), sc(cross(cur.dC, cur.C), k))
    A0, A1, A2, C, S = cur.dX0, cur.dC, cur.dS, cur.C, cur.S
    A1S, A2C, A2S, A1C = cross(A1, S), cross(A2, C), cross(A2, S), cross(A1, C)
    h = k*mp.mpf(0.5)
    return (sc(sub(A1S, A2C), h), sc(cross(A0, S), k), sc(cross(A0, C), -k), sc(add(A1S, A2C), h),
            sc(sub(A2S, A1C), h))


def fmul(A, B):
    """The product of two trigonometric series `(a, b)` (`a0 + sum a_k cos
    kb + b_k sin kb`), exactly by the product-to-sum rules."""
    n = len(A[0])+len(B[0])-1
    a, b = [Z]*n, [Z]*n
    for m in range(len(A[0])):
        am, bm = A[0][m], A[1][m]
        for k in range(len(B[0])):
            ak, bk = B[0][k], B[1][k]
            p, q = m+k, abs(m-k)
            sg = 1 if m >= k else -1
            # cos m cos k, sin m sin k, sin m cos k, cos m sin k
            h = (am*ak)/2
            a[p] += h
            a[q] += h
            h = (bm*bk)/2
            a[q] += h
            a[p] -= h
            h = (bm*ak)/2
            b[p] += h
            b[q] += sg*h
            h = (am*bk)/2
            b[p] += h
            b[q] -= sg*h
    b[0] = Z
    return a, b


def fint(A, cs):
    """`int_{b0}^{b1}` of a series, `cs = [(cos kb, sin kb)] at b1 minus at
    b0` (k >= 1) and `b1 - b0` first."""
    a, b = A
    out = a[0]*cs[0]
    for k in range(1, len(a)):
        dc, ds = cs[k]
        out += (a[k]*ds-b[k]*dc)/k
    return out


def multiples(t, K):
    c, s = mp.cos(t), mp.sin(t)
    out = [(mp.mpf(1), Z), (c, s)]
    for _ in range(2, K+1):
        pc, ps = out[-1]
        out.append((pc*c-ps*s, ps*c+pc*s))
    return out


def integrate_pieces(cur, pieces, sign, c0):
    """Per class: area, volume and first moments (about `c0`) of the
    pieces, by the divergence theorem with the outward `sign * X_a x X_b`:
    the volume's and moments' integrands in closed form (on a circle
    trigonometric polynomials of degree at most four multiplied out exactly,
    on a line polynomials), the area element `|X_a x X_b|` by 24-point
    Gauss-Legendre on pieces of at most a quarter turn."""
    out = {'in': [Z]*5, 'out': [Z]*5}
    jet = normal_jet(cur, sign)
    X0 = sub(cur.X0, c0)
    C, S = cur.C, cur.S
    if cur.kind == 'circle':
        N0, Nc, Ns, N2c, N2s = jet
        Ds = [([X0[i], C[i]], [Z, S[i]]) for i in range(3)]
        Nsr = [([N0[i], Nc[i], N2c[i]], [Z, Ns[i], N2s[i]]) for i in range(3)]
        vol = None
        mom = []
        for i in range(3):
            p = fmul(Ds[i], Nsr[i])
            vol = p if vol is None else ([x+y for x, y in zip(vol[0], p[0])], [x+y for x, y in zip(vol[1], p[1])])
            mom.append(fmul(Ds[i], p))
    else:
        N0, N1 = jet
    for b0, b1, cls, _, _ in pieces:
        acc = out[cls]
        if cur.kind == 'circle':
            e0, e1 = multiples(b0, 4), multiples(b1, 4)
            cs = [b1-b0]+[(e1[k][0]-e0[k][0], e1[k][1]-e0[k][1]) for k in range(1, 5)]
            acc[1] += fint(vol, cs)/3
            for i in range(3):
                acc[2+i] += fint(mom[i], cs)/2
            n = max(1, int(mp.ceil((b1-b0)/quarter())))
        else:
            # d = X0 + b C, N = N0 + b N1: the integrands' polynomials in b.
            def pint(coef):
                return sum(c*(b1**(k+1)-b0**(k+1))/(k+1) for k, c in enumerate(coef))
            vol = [Z, Z, Z]
            for i in range(3):
                x0, x1, n0, n1 = X0[i], C[i], N0[i], N1[i]
                vol[0] += x0*n0
                vol[1] += x0*n1+x1*n0
                vol[2] += x1*n1
                sq = (x0*x0, 2*x0*x1, x1*x1)
                acc[2+i] += pint([sq[0]*n0, sq[0]*n1+sq[1]*n0, sq[1]*n1+sq[2]*n0, sq[2]*n1])/2
            acc[1] += pint(vol)/3
            n = 1
        step = (b1-b0)/n
        a0 = Z
        for i in range(n):
            lo = b0+i*step
            for xi, wi in nodes(INNER):
                b = lo+step*(xi+1)/2
                if cur.kind == 'line':
                    N = (N0[0]+b*N1[0], N0[1]+b*N1[1], N0[2]+b*N1[2])
                else:
                    c, s = mp.cos(b), mp.sin(b)
                    c2, s2 = c*c-s*s, 2*c*s
                    N = tuple(N0[j]+c*Nc[j]+s*Ns[j]+c2*N2c[j]+s2*N2s[j] for j in range(3))
                a0 += wi*mp.sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2])
        acc[0] += a0*step/2
    return out


# ------------------------------------------------------------------ a family swept

class Sweep:
    """A face's family against the other input: its events (the scan and
    bisection on the curves' structure), the family's own kinks, and the
    classes' measures integrated between them."""

    def __init__(self, fam, other, c0, size):
        self.fam, self.other, self.c0, self.size = fam, other, c0, size
        self.missed = 0
        self.quad_error = Z
        self.events = None
        self.samples = []   # (a, curve, pieces) of the scan, for the margins

    def structure(self, a):
        cur = self.fam.curve(a)
        pcs = curve_pieces(cur, self.other)
        return signature(cur, pcs), cur, pcs

    def find_events(self, keep_samples=False):
        fam = self.fam
        a0, a1 = fam.a0, fam.a1
        span = a1-a0
        eps = span*mp.mpf(10)**-35
        if fam.periodic:
            grid = [a0+span*(j+mp.mpf(0.3819660112501051))/SCAN for j in range(SCAN)]
        else:
            delta = span*mp.mpf(10)**-25
            grid = [a0+delta]+[a0+span*j/SCAN for j in range(1, SCAN)]+[a1-delta]
        for k in fam.own:
            grid += [k-span*mp.mpf(10)**-25, k+span*mp.mpf(10)**-25]
        grid.sort()
        sigs = []
        for a in grid:
            s, cur, pcs = self.structure(a)
            sigs.append(s)
            if keep_samples:
                self.samples.append((a, cur, pcs))
        found = []

        def bisect(lo, slo, hi, shi):
            if hi-lo <= eps:
                found.append((lo+hi)/2)
                return
            m = (lo+hi)/2
            sm = self.structure(m)[0]
            if sm != slo and sm != shi and hi-lo <= span*mp.mpf(10)**-20:
                # Two surfaces tangent along an edge (a profile's segment
                # meeting its arc tangentially): their roots coincide to
                # second order there, which bounds the piece rounding's
                # to call; the event is placed within 1e-20.
                found.append(m)
                return
            if sm != slo:
                bisect(lo, slo, m, sm)
            if sm != shi:
                bisect(m, sm, hi, shi)
        n = len(grid)
        pairs = list(zip(range(n-1), range(1, n)))
        if fam.periodic:
            pairs.append((n-1, 0))
        for i, j in pairs:
            if sigs[i] != sigs[j]:
                lo, hi = grid[i], grid[j]
                if j == 0:
                    hi += span
                bisect(lo, sigs[i], hi, sigs[j])
        out = []
        for e in found:
            if fam.periodic:
                e = a0+(e-a0) % span
            out.append(e)
        self.events = sorted(out)
        return self.events

    def breakpoints(self):
        fam = self.fam
        pts = [fam.a0, fam.a1]+list(self.events)+list(fam.own)
        pts = sorted(set(pts))
        return [p for p in pts if fam.a0 <= p <= fam.a1]

    def measure(self):
        """{class: [area, volume, moments about c0]} of the face."""
        if self.events is None:
            self.find_events()
        tol = mp.mpf(10)**-32*self.size**4
        tot = [Z]*10
        pts = self.breakpoints()
        for lo, hi in zip(pts, pts[1:]):
            if hi-lo <= (self.fam.a1-self.fam.a0)*mp.mpf(10)**-34:
                continue
            ref = self.structure((lo+hi)/2)[0]

            def f(a, ref=ref):
                s, cur, pcs = self.structure(a)
                if s != ref:
                    self.missed += 1
                m = integrate_pieces(cur, pcs, self.fam.sign, self.c0)
                return m['in']+m['out']
            est, err = integrate(f, lo, hi, tol)
            self.quad_error += err
            tot = [x+y for x, y in zip(tot, est)]
        return {'in': tot[:5], 'out': tot[5:]}


# ------------------------------------------------------------------ the pair

def box_size(boxes):
    return max([1.0]+[abs(v) for lo, hi in boxes for v in lo+hi])


class Pair:
    """A whole torus and a prism with arcs, a sphere, a cone or another whole
    torus, in either order: `A` the object, `B` the tool."""

    def __init__(self, obj, tool):
        self.A, self.B = make_input(obj), make_input(tool)
        assert 'torus' in (self.A.kind, self.B.kind), 'S9d.4b.2: a whole torus against a curved solid'
        self.K = self.A if self.A.kind == 'torus' else self.B
        self.P = self.B if self.K is self.A else self.A
        self.size = mp.mpf(box_size([self.A.box(), self.B.box()]))
        self.c0 = self.K.fr.om
        self._sweeps = None
        self._counts = {}

    def sweeps(self):
        """{(family index, role, face name): Sweep} measured."""
        if self._sweeps is None:
            out = {}
            for k in (0, 1):
                for role, S, O in (('A', self.A, self.B), ('B', self.B, self.A)):
                    for face in S.faces:
                        sw = Sweep(face.families[k], O, self.c0, self.size)
                        sw.find_events(keep_samples=(k == 0))
                        sw.result = sw.measure()
                        out[(k, role, face.name)] = sw
            self._sweeps = out
        return self._sweeps

    def classes(self, k, role):
        """{class: [area, volume, moments about c0]} of an input's faces in
        family `k`."""
        tot = {'in': [Z]*5, 'out': [Z]*5}
        for (kk, rr, _), sw in self.sweeps().items():
            if kk == k and rr == role:
                for c in tot:
                    tot[c] = [x+y for x, y in zip(tot[c], sw.result[c])]
        return tot

    def to_world(self, vec):
        """(volume, world moments, area) from [area, volume, moments about c0]."""
        V = vec[1]
        return V, tuple(self.c0[i]*V+vec[2+i] for i in range(3)), vec[0]

    def op_vec(self, k, op):
        a, b = self.classes(k, 'A'), self.classes(k, 'B')
        if op == 'common':
            return [x+y for x, y in zip(a['in'], b['in'])]
        if op == 'fuse':
            return [x+y for x, y in zip(a['out'], b['out'])]
        inv = b['in']
        return [a['out'][0]+inv[0]]+[x-y for x, y in zip(a['out'][1:], inv[1:])]

    def measures(self, op, k=0):
        """(volume, world moments, area) of an operation by family `k`."""
        return self.to_world(self.op_vec(k, op))

    def input_measures(self, role, k=0):
        c = self.classes(k, role)
        return self.to_world([x+y for x, y in zip(c['in'], c['out'])])

    def volume(self, op, k=0):
        V, m, _ = self.measures(op, k)
        return V, m

    def area(self, op, k=0):
        return self.measures(op, k)[2]

    def result(self, op):
        """(solids, volume, area, centre)."""
        V, m, A = self.measures(op)
        n = self.solids(op)
        if n == 0:
            return 0, Z, Z, (Z, Z, Z)
        return n, V, A, tuple(x/V for x in m)

    def solids(self, op):
        key = {'common': 'common', 'fuse': 'fuse', 'cut': 'A-B'}[op]
        return self.swept_count(key)

    def swept_count(self, key):
        if key not in self._counts:
            self._counts.update(sweep_counts(self))
        return self._counts[key]

    def missed(self):
        return sum(sw.missed for sw in self.sweeps().values())

    def quad_error(self):
        return max(sw.quad_error for sw in self.sweeps().values())


# ------------------------------------------------------------------ solids: a sweep of the torus's slices

def fpoly_mul(p, q):
    out = [0.0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def fpoly_add(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0.0)+(q[i] if i < len(q) else 0.0) for i in range(n)]


class FloatModel:
    """An input in binary64 for the sweep: its surfaces' polynomials along a
    ray `P + t D` (ascending in `t`) and its membership."""

    def __init__(self, S):
        self.S = S
        fr = S.fr
        self.o = [float(v) for v in fr.om]
        self.inv = [[float(v) for v in row] for row in fr.invm]
        self.kind = S.kind
        if S.kind == 'torus':
            self.R, self.r = float(S.R), float(S.r)
        elif S.kind == 'sphere':
            self.c, self.R2 = [float(v) for v in S.cm], float(S.R)**2
        elif S.kind == 'cone':
            self.r0, self.k, self.h = float(S.r0), float(S.k), float(S.h)
        else:
            self.lo, self.hi = float(S.lo), float(S.hi)
            self.els = []
            for e in S.profile.els:
                if e[0] == 'seg':
                    self.els.append(('seg', [float(v) for v in e[1]], [float(v) for v in e[2]]))
                else:
                    self.els.append(('arc', [float(v) for v in e[1]], float(e[2]), float(e[3]), float(e[4])))

    def chart(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        return [sum(r[i]*d[i] for i in range(3)) for r in self.inv]

    def chart_vec(self, D):
        return [sum(r[i]*D[i] for i in range(3)) for r in self.inv]

    def polys(self, P, D):
        if self.kind == 'sphere':
            q = [P[i]-self.c[i] for i in range(3)]
            return [[sum(x*x for x in q)-self.R2, 2*sum(q[i]*D[i] for i in range(3)), sum(x*x for x in D)]]
        p0, p1 = self.chart(P), self.chart_vec(D)
        lin_ = lambda i: [p0[i], p1[i]]
        sq = lambda i: fpoly_mul(lin_(i), lin_(i))
        if self.kind == 'torus':
            s = fpoly_add(fpoly_add(sq(0), sq(1)), sq(2))
            s[0] += self.R*self.R-self.r*self.r
            h = fpoly_add(sq(0), sq(1))
            return [fpoly_add(fpoly_mul(s, s), [-4*self.R*self.R*c for c in h])]
        if self.kind == 'cone':
            rad = [self.r0+self.k*p0[2], self.k*p1[2]]
            return [fpoly_add(fpoly_add(sq(0), sq(1)), [-c for c in fpoly_mul(rad, rad)]), lin_(2),
                    [p0[2]-self.h, p1[2]]]
        out = [lin_(2), [p0[2]-self.hi, p1[2]]]
        for e in self.els:
            if e[0] == 'seg':
                p, q = e[1], e[2]
                a = (q[1]-p[1], -(q[0]-p[0]))
                out.append([a[0]*(p0[0]-p[0])+a[1]*(p0[1]-p[1]), a[0]*p1[0]+a[1]*p1[1]])
            else:
                C, a = e[1], e[2]
                u, v = [p0[0]-C[0], p1[0]], [p0[1]-C[1], p1[1]]
                out.append(fpoly_add(fpoly_add(fpoly_mul(u, u), fpoly_mul(v, v)), [-a*a]))
        return out

    def contains(self, X):
        if self.kind == 'sphere':
            return sum((X[i]-self.c[i])**2 for i in range(3)) < self.R2
        p = self.chart(X)
        if self.kind == 'torus':
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+self.R*self.R-self.r*self.r
            return s*s < 4*self.R*self.R*(p[0]*p[0]+p[1]*p[1])
        if self.kind == 'cone':
            return 0 < p[2] < self.h and p[0]*p[0]+p[1]*p[1] < (self.r0+self.k*p[2])**2
        if not self.lo < p[2] < self.hi:
            return False
        u, v = p[0], p[1]
        for e in self.els:
            if e[0] == 'seg':
                a, b = e[1], e[2]
                if (b[0]-a[0])*(v-a[1])-(b[1]-a[1])*(u-a[0]) <= 0:
                    return False
            else:
                C, a, t0, t1 = e[1], e[2], e[3], e[4]
                if t1-t0 >= 2*math.pi-1e-12:
                    if (u-C[0])**2+(v-C[1])**2 >= a*a:
                        return False
                    continue
                pa = (C[0]+a*math.cos(t0), C[1]+a*math.sin(t0))
                pb = (C[0]+a*math.cos(t1), C[1]+a*math.sin(t1))
                side = (pb[0]-pa[0])*(v-pa[1])-(pb[1]-pa[1])*(u-pa[0])
                if side <= 0 and (u-C[0])**2+(v-C[1])**2 >= a*a:
                    return False
        return True

    def intervals(self, P, D, tmax):
        """The ray's parameters inside, `[(t0, t1)]` within `[0, tmax]`."""
        ts = []
        for p in self.polys(P, D):
            while len(p) > 1 and abs(p[-1]) <= 1e-14*max(abs(c) for c in p):
                p = p[:-1]
            if len(p) == 2:
                ts.append(-p[0]/p[1])
            elif len(p) == 3:
                a, b, c = p[2], p[1], p[0]
                disc = b*b-4*a*c
                if disc > 0:
                    sq = math.sqrt(disc)
                    q = -(b+math.copysign(sq, b))/2
                    ts += [q/a, c/q] if q != 0 else [-b/(2*a)]
            elif len(p) > 3:
                zs = float_roots(list(reversed([complex(c) for c in p])))
                if zs is None:
                    with mp.workdps(30):
                        zs = [complex(z) for z in mp.polyroots(list(reversed(p)), maxsteps=400, extraprec=60,
                                                               error=False)]
                ts += [z.real for z in zs if abs(z.imag) <= 1e-7*(1+abs(z))]
        ts = sorted(t for t in ts if 0 < t < tmax)
        pts = [0.0]+ts+[tmax]
        out = []
        for a, b in zip(pts, pts[1:]):
            if b > a and self.contains([P[i]+(a+b)/2*D[i] for i in range(3)]):
                if out and abs(out[-1][1]-a) <= 1e-12:
                    out[-1] = (out[-1][0], b)
                else:
                    out.append((a, b))
        return out


def iv_ops(A, B):
    """Intersection, union and both differences of two sorted interval lists."""
    pts = sorted(set([x for iv in A+B for x in iv]))
    inside = lambda ivs, t: any(a < t < b for a, b in ivs)
    out = {'common': [], 'fuse': [], 'A-B': [], 'B-A': []}
    for a, b in zip(pts, pts[1:]):
        m = (a+b)/2
        ia, ib = inside(A, m), inside(B, m)
        for key, keep in (('common', ia and ib), ('fuse', ia or ib), ('A-B', ia and not ib),
                          ('B-A', ib and not ia)):
            if keep:
                lst = out[key]
                if lst and lst[-1][1] == a:
                    lst[-1] = (lst[-1][0], b)
                else:
                    lst.append((a, b))
    return out


def sweep_counts(pair, levels=200, rays=512):
    """Solids of every operation by a sweep of the torus's normal slices: in
    each slice the operation's intervals along rays from the axis, joined
    where they overlap on adjacent rays, at the axis, and on the same ray of
    the adjacent slice; the connected components counted."""
    K = pair.K
    fK = FloatModel(K)
    fA, fB = FloatModel(pair.A), FloatModel(pair.B)
    o = fK.o
    x, y, n = ([float(v) for v in K.fr.xm], [float(v) for v in K.fr.ym], [float(v) for v in K.fr.nm])
    corners = []
    for lo, hi in (pair.A.box(), pair.B.box()):
        corners += [fK.chart([(lo, hi)[i >> k & 1][k] for k in range(3)]) for i in range(8)]
    wmin, wmax = min(c[2] for c in corners)-0.01, max(c[2] for c in corners)+0.01
    tmax = max(math.hypot(c[0], c[1]) for c in corners)*1.5+1
    R, r = float(K.R), float(K.r)
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    keys = ('common', 'fuse', 'A-B', 'B-A')
    prev = None
    for l in range(levels):
        s = wmin+(wmax-wmin)*(l+0.5)/levels
        row = []
        for i in range(rays):
            th = 2*math.pi*(i+0.25)/rays
            c, sn = math.cos(th), math.sin(th)
            P = [o[k]+s*n[k] for k in range(3)]
            D = [c*x[k]+sn*y[k] for k in range(3)]
            iv = {}
            for role, fm in (('A', fA), ('B', fB)):
                if fm.S is K:
                    iv[role] = [(R-math.sqrt(r*r-s*s), R+math.sqrt(r*r-s*s))] if abs(s) < r else []
                else:
                    iv[role] = fm.intervals(P, D, tmax)
            ops = iv_ops(iv['A'], iv['B'])
            for key in keys:
                for j, _ in enumerate(ops[key]):
                    parent[(key, l, i, j)] = (key, l, i, j)
            row.append(ops)
        for key in keys:
            axis = None
            for i in range(rays):
                cur, nxt = row[i][key], row[(i+1) % rays][key]
                for j, (a, b) in enumerate(cur):
                    if a <= 1e-9:
                        if axis is not None:
                            union(axis, (key, l, i, j))
                        axis = (key, l, i, j)
                    for jj, (a2, b2) in enumerate(nxt):
                        if min(b, b2)-max(a, a2) > 1e-9:
                            union((key, l, i, j), (key, l, (i+1) % rays, jj))
                    if prev is not None:
                        for jj, (a2, b2) in enumerate(prev[i][key]):
                            if min(b, b2)-max(a, a2) > 1e-9:
                                union((key, l, i, j), (key, l-1, i, jj))
        prev = row
    counts = {}
    for key in keys:
        counts[key] = len({find(a) for a in parent if a[0] == key})
    return counts


# ------------------------------------------------------------------ margins

def prism_edges(P):
    """A prism's edges as curves: every element's bottom and top edges, and
    the vertical edge at each vertex."""
    out = []
    for face in P.faces:
        if face.name in ('start', 'end'):
            continue
        gen, row = face.families
        out += [row.curve(row.a0), row.curve(row.a1)]
        if not gen.periodic:
            out.append(gen.curve(gen.a0))
    return out


def edges_of(S):
    if S.kind == 'prism':
        return prism_edges(S)
    if S.kind == 'cone':
        par = S.faces[0].families[1]
        return [par.curve(w) for w, rad in ((par.a0, S.r0m), (par.a1, S.r1m)) if rad > 0]
    return []


def vertices_of(S):
    if S.kind == 'prism':
        return [S.fr.world(M(u), M(v), h) for u, v in S.profile.vertices for h in (S.lom, S.him)]
    if S.kind == 'cone':
        a = S.apex()
        return [] if a is None else [a]
    return []


def critical_points(cur, surf):
    """Where `f` is stationary along the curve (inside its range)."""
    if cur.kind == 'line':
        c, m, h = line_coeffs(cur, surf)
        dc = [k*c[k] for k in range(1, len(c))]
        while len(dc) > 1 and dc[-1] == 0:
            dc.pop()
        if not dc or all(x == 0 for x in dc):
            return []
        out = []
        for z in mp_roots(list(reversed(dc))):
            if abs(mp.im(z)) <= REAL_TOL*(1+abs(z)) and -1 < mp.re(z) < 1:
                out.append(m+h*mp.re(z))
        return out
    a, b = trig_coeffs(cur, surf)
    da = [Z]+[k*b[k] for k in range(1, len(a))]
    db = [Z]+[-k*a[k] for k in range(1, len(a))]
    rts = circle_roots(da, db)
    if cur.rng is None:
        return rts
    lo, hi = cur.rng
    return [lo+(t-lo) % (2*mp.pi) for t in rts if lo < lo+(t-lo) % (2*mp.pi) < hi]


def distance_estimate(surf, X):
    g = surf.grad(X)
    gn = norm(g)
    return abs(surf.f(X))/gn if gn > 0 else mp.inf


def on_face(O, surf, X, size):
    """Whether the foot of `X` on the surface (Newton's steps along the
    gradient) lies on a face of `O` (inside `O` just behind it, outside just
    in front): a surface's extension beyond its face does not count."""
    foot = X
    for _ in range(30):
        g = surf.grad(foot)
        g2 = dot(g, g)
        if g2 == 0:
            return False
        k = surf.f(foot)/g2
        foot = (foot[0]-k*g[0], foot[1]-k*g[1], foot[2]-k*g[2])
        if abs(k)*mp.sqrt(g2) < size*mp.mpf(10)**-12:
            break
    g = surf.grad(foot)
    e = size*mp.mpf(10)**-6/mp.sqrt(dot(g, g))
    return O.contains(sub(foot, scale(g, e))) and not O.contains(add(foot, scale(g, e)))


def sine(u, v):
    nu, nv = norm(u), norm(v)
    if nu == 0 or nv == 0:
        return mp.mpf(1)
    return norm(cross(u, v))/(nu*nv)


def margins(pair):
    """The least sine between the inputs' surfaces where they meet (at the
    first family's scanned curves' class boundaries), the least distance
    (first order, `|f| / |grad f|`) from a face to another input's surface it
    never crosses (the surface's critical values along the scanned curves),
    the least sine at which an edge crosses the torus and its least distance
    from it where stationary, and the vertices' distances: relative to the
    case's size."""
    size = pair.size
    out = {'face_crossing': mp.inf, 'face_gap': mp.inf, 'edge_crossing': mp.inf, 'edge_gap': mp.inf,
           'vertex': mp.inf}
    for (k, role, name), sw in pair.sweeps().items():
        if k != 0:
            continue
        O = pair.B if role == 'A' else pair.A
        crossed = set()
        for a, cur, pcs in sw.samples:
            n = len(pcs)
            for i in range(n):
                p, q = pcs[i], pcs[(i+1) % n]
                closed = cur.kind == 'circle' and cur.rng is None
                if i+1 == n and not closed:
                    continue
                if p[2] == q[2] or p[4] in (None, 'end'):
                    continue
                j = p[4]
                crossed.add(j)
                X, Xa, Xb = cur.frame_at(p[1])
                out['face_crossing'] = min(out['face_crossing'], sine(cross(Xa, Xb), O.surfs[j].grad(X)))
        for j, surf in enumerate(O.surfs):
            if j in crossed:
                continue
            for a, cur, pcs in sw.samples[::4]:
                if cur.degenerate():
                    continue
                for b in critical_points(cur, surf):
                    X = cur.point(b)
                    if on_face(O, surf, X, size):
                        out['face_gap'] = min(out['face_gap'], distance_estimate(surf, X)/size)
    K = pair.K
    for S in (pair.A, pair.B):
        if S is K:
            continue
        O = pair.B if S is pair.A else pair.A
        for cur in edges_of(S):
            if cur.degenerate():
                continue
            for surf in O.surfs:
                for b in curve_roots(cur, surf):
                    X, _, Xb = cur.frame_at(b)
                    g = surf.grad(X)
                    out['edge_crossing'] = min(out['edge_crossing'], abs(dot(Xb, g))/(norm(Xb)*norm(g)))
                for b in critical_points(cur, surf):
                    out['edge_gap'] = min(out['edge_gap'], distance_estimate(surf, cur.point(b))/size)
        for X in vertices_of(S):
            for surf in O.surfs:
                out['vertex'] = min(out['vertex'], distance_estimate(surf, X)/size)
    return out


def near_coincidences(pair):
    """A face's class, or a result's volume, positive but thinner than 1e-9
    of the case's size; events of one family closer than that."""
    size = pair.size
    out = []
    for (k, role, name), sw in pair.sweeps().items():
        for c, vec in sw.result.items():
            if mp.mpf(10)**-30*size**2 < abs(vec[0]) < mp.mpf(10)**-9*size**2:
                out.append(f'{role} {name} family {k} {c} {mp.nstr(vec[0], 3)}')
        span = sw.fam.a1-sw.fam.a0
        ev = sorted(sw.events)
        for p, q in zip(ev, ev[1:]):
            if mp.mpf(10)**-30*span < q-p < mp.mpf(10)**-9*span:
                out.append(f'{role} {name} family {k} events {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    for op in OPS:
        v = pair.volume(op)[0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    return out


def rows(obj, operation, tool, pair=None):
    """`result N volume area cx cy cz` or `empty`, and the pair (reused
    across the three operations)."""
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair

