#!/usr/bin/env python3
"""Independent reference for S7b.1 of REVIEW_NOTES.md: the procedural
intersection curves of two cylinders with crossing axes and of a cylinder
and a sphere off its axis.

Surfaces are the exact point sets of their stored data, as in S7a. The
curve's class is decided by exact rational predicates: for two cylinders
the squared distance of the axes `d^2` against `(r1 + r2)^2` and
`(r1 - r2)^2`; for a cylinder and a sphere the distance `e` of the centre
from the axis against `|e - r|` and `e + r` (compared through squares). The
curve is parameterised on the ruled surface (the thinner cylinder; the
cylinder of a cylinder/sphere pair) by the angle `u` of its ruling in the
canonical frame: `x` the unit common normal of the axes pointing towards
the other axis (or the unit direction from the axis towards the sphere's
centre), `y = a x x` with `a` the unit axis; of two equal cylinders the
ruled one comes first by its stored normal, then origin. A ruling `o + r (cos u x + sin u y) + v a` meets the other
quadric where `A v^2 + 2 B(u) v + C(u) = 0`; the curve is
`v = (-B +- sqrt(D)) / A`, `D = B^2 - A C`. The reference finds the roots of
`D` by dense 80-digit sampling and `findroot`, checks their number against
the class, and prints canonical rows:

* `empty`;
* `point x y z` (a tangent point);
* `loop u0 u1 p0 p1 m+ m-`: the loop's parameter range (u1 > u0, possibly
  beyond pi), its two ends and its points at the middle parameter on the
  `+` and `-` branches;
* `rings p+ p- q+ q-`: two rings around the ruled surface, by their points
  at `u = 0` and `u = pi` on each branch;
* `figure_eight node p+ p-`: two rings touching at the node, and the points
  opposite it (`u` of the node plus pi) on each branch.

The kernel uses the same parameterisation and frame, so every row is
directly comparable; it evaluates in certified intervals and finds roots by
sign changes, not by sampling.

S7b.2 adds cones (`ConeCurve`, `cone_form`); S7b.3a a torus with a plane or
a sphere on its meridians (`TorusCurve`: `D` from the other surface
evaluated on the meridian circle at three angles, classes from its own exact
predicates on the quadratic `P(m cos phi)` found by exact interpolation) and
the special and coaxial torus pairs as `circle` rows (centre, unit normal,
radius), `same`, or `not_conic` for a sphere containing a meridian circle.
"""
from fractions import Fraction as F

import mpmath as mp

import analytic_intersection_reference as ana
from analytic_intersection_reference import cross, dot, mpf, scale, sub, zero
from identity_reference import frame_axes

mp.mp.dps = 80


def unit(v):
    m = [mpf(x) for x in v]
    n = mp.sqrt(sum(x*x for x in m))
    return [x/n for x in m]


def frame(o, a, toward):
    """The canonical frame of a ruled cylinder: unit axis, x along the
    rational direction `toward` (perpendicular to the axis), y = a x x."""
    au = unit(a)
    xu = unit(toward)
    yu = [au[1]*xu[2]-au[2]*xu[1], au[2]*xu[0]-au[0]*xu[2], au[0]*xu[1]-au[1]*xu[0]]
    return [mpf(v) for v in o], au, xu, yu


class Curve:
    def __init__(self, fr, r, other):
        self.o, self.a, self.x, self.y = fr
        self.r = mpf(r)
        self.other = other     # f(p) and its quadratic form on the axis
        self.A = other['form'](self.a)

    def ruling(self, u):
        c, s = mp.cos(u), mp.sin(u)
        return [o+self.r*(c*x+s*y) for o, x, y in zip(self.o, self.x, self.y)]

    def coefficients(self, u):
        """(B, C): f(P + v a) = A v^2 + 2 B v + C."""
        p = self.ruling(u)
        f = self.other['f']
        c0 = f(p)
        f1 = f([pi+ai for pi, ai in zip(p, self.a)])
        # f(P + a) = A + 2 B + C.
        return (f1-self.A-c0)/2, c0

    def D(self, u):
        B, C = self.coefficients(u)
        return B*B-self.A*C

    def point(self, u, sign):
        B, C = self.coefficients(u)
        disc = B*B-self.A*C
        root = mp.sqrt(disc) if disc > 0 else mp.mpf(0)
        v = (-B+sign*root)/self.A
        p = self.ruling(u)
        return [pi+v*ai for pi, ai in zip(p, self.a)]

    def roots(self, samples=2880):
        """Sign changes of D over one turn, refined; u in (-pi, pi]."""
        us = [-mp.pi+2*mp.pi*(k+0.5)/samples for k in range(samples)]
        ds = [self.D(u) for u in us]
        out = []
        for k in range(samples):
            j = (k+1) % samples
            if (ds[k] > 0) != (ds[j] > 0):
                lo, hi = us[k], us[j] if j else us[j]+2*mp.pi
                root = mp.findroot(self.D, (lo, hi), solver='anderson')
                root = mp.atan2(mp.sin(root), mp.cos(root))
                out.append(root)
        return sorted(out)


def cylinder_form(o, a, r):
    """f(p) = |p - o|^2 - ((p - o) . a)^2 / |a|^2 - r^2 and its quadratic
    form on a direction."""
    om, am, rm = [mpf(v) for v in o], [mpf(v) for v in a], mpf(r)
    aa = sum(v*v for v in am)

    def f(p):
        rel = [x-y for x, y in zip(p, om)]
        h = sum(x*y for x, y in zip(rel, am))
        return sum(x*x for x in rel)-h*h/aa-rm*rm

    def form(d):
        h = sum(x*y for x, y in zip(d, am))
        return sum(x*x for x in d)-h*h/aa
    return {'f': f, 'form': form}


def sphere_form(c, R):
    cm, Rm = [mpf(v) for v in c], mpf(R)
    return {'f': lambda p: sum((x-y)**2 for x, y in zip(p, cm))-Rm*Rm,
            'form': lambda d: sum(x*x for x in d)}


def cone_form(o, a, r, angle):
    """A double cone: f(p) = rho^2 cos^2 - (r cos + h sin)^2 with h along
    the unit axis from o, and its quadratic form on a direction."""
    om, am = [mpf(v) for v in o], unit(a)
    ca, sa, rm = mp.cos(mpf(angle)), mp.sin(mpf(angle)), mpf(r)

    def f(p):
        rel = [x-y for x, y in zip(p, om)]
        h = sum(x*y for x, y in zip(rel, am))
        rho2 = sum(x*x for x in rel)-h*h
        return rho2*ca*ca-(rm*ca+h*sa)**2

    def form(d):
        h = sum(x*y for x, y in zip(d, am))
        return sum(x*x for x in d)*ca*ca-h*h
    return {'f': f, 'form': form}


class ConeCurve:
    """A sphere seen from a cone's apex: rulings V + v d(u) with
    d(u) = cos a axis + sin a (cos u x + sin u y), x towards the centre."""
    def __init__(self, cone, sphere):
        o, a = cone.axes()
        c, _ = sphere.axes()
        self.a = unit(a)
        ca, sa = mp.cos(mpf(cone.angle)), mp.sin(mpf(cone.angle))
        self.ca, self.sa = ca, sa
        back = mpf(cone.radius)/mp.tan(mpf(cone.angle))
        self.V = [mpf(x)-back*y for x, y in zip(o, self.a)]
        cm = [mpf(x) for x in c]
        rel = [x-y for x, y in zip(cm, self.V)]
        h = sum(x*y for x, y in zip(rel, self.a))
        perp = [x-h*y for x, y in zip(rel, self.a)]
        n = mp.sqrt(sum(x*x for x in perp))
        self.x = [x/n for x in perp]
        self.y = [self.a[1]*self.x[2]-self.a[2]*self.x[1], self.a[2]*self.x[0]-self.a[0]*self.x[2],
                  self.a[0]*self.x[1]-self.a[1]*self.x[0]]
        self.c, self.R = cm, mpf(sphere.radius)
        self.C = sum(x*x for x in rel)-self.R**2

    def d(self, u):
        c, s = mp.cos(u), mp.sin(u)
        return [self.ca*a+self.sa*(c*x+s*y) for a, x, y in zip(self.a, self.x, self.y)]

    def B(self, u):
        return sum((v-c)*d for v, c, d in zip(self.V, self.c, self.d(u)))

    def D(self, u):
        return self.B(u)**2-self.C

    def point(self, u, sign):
        disc = self.D(u)
        root = mp.sqrt(disc) if disc > 0 else mp.mpf(0)
        v = -self.B(u)+sign*root
        return [p+v*q for p, q in zip(self.V, self.d(u))]

    def roots(self, samples=2880):
        return Curve.roots(self, samples)


def classify_cylinders(c1, c2):
    o1, a1 = c1.axes()
    o2, a2 = c2.axes()
    r1, r2 = F(c1.radius), F(c2.radius)
    m = cross(a1, a2)
    d2 = dot(sub(o2, o1), m)**2/dot(m, m)
    if d2 > (r1+r2)**2:
        return 'empty'
    if d2 == (r1+r2)**2:
        return 'point'
    if r1 != r2 and d2 == (r1-r2)**2:
        return 'figure_eight'
    if d2 < (r1-r2)**2:
        return 'rings'
    return 'loop'


def classify_sphere(c, s):
    o, a = c.axes()
    ctr, _ = s.axes()
    r, R = F(c.radius), F(s.radius)
    w = sub(ctr, o)
    w = sub(w, scale(a, dot(w, a)/dot(a, a)))
    e2 = dot(w, w)
    # R against |e - r| and e + r, e = sqrt(e2): compare through squares.
    def cmp_sum(k):
        # sign of R - (e + k r) for k = +1, -1 (e + k r may be negative).
        t = R-k*r            # R - k r against e
        if t < 0:
            return -1
        return (t*t > e2)-(t*t < e2)
    def cmp_abs():
        # R against |e - r|: R^2 against e2 - 2 e r + r^2, i.e. 2 e r
        # against e2 + r^2 - R^2.
        rhs = e2+r*r-R*R
        if rhs < 0:
            return 1       # 2 e r >= 0 > rhs: R > |e - r|
        lhs2 = 4*e2*r*r
        return (lhs2 > rhs*rhs)-(lhs2 < rhs*rhs)
    low = cmp_abs()        # sign of R - |e - r|
    if low < 0:
        return 'empty'
    if low == 0:
        return 'point'
    high = cmp_sum(1)      # sign of R - (e + r)
    if high < 0:
        return 'loop'
    if high == 0:
        return 'figure_eight'
    return 'rings'


def curve(s1, s2):
    """(class, Curve) of a cylinder pair with crossing axes or a cylinder
    and a sphere off its axis; None for other pairs."""
    if 'torus' in (s1.kind, s2.kind):
        torus, other = (s1, s2) if s1.kind == 'torus' else (s2, s1)
        if torus_special(torus, other) == 'general':
            return 'torus', torus_curve(torus, other)
        return None
    if ana.KINDS.index(s1.kind) > ana.KINDS.index(s2.kind):
        s1, s2 = s2, s1
    if (s1.kind, s2.kind) == ('cylinder', 'cylinder'):
        o1, a1 = s1.axes()
        o2, a2 = s2.axes()
        if zero(cross(a1, a2)):
            return None
        if F(s1.radius) == F(s2.radius) and dot(sub(o2, o1), cross(a1, a2)) == 0:
            return None        # S7a's ellipses
        cls = classify_cylinders(s1, s2)
        # The thinner cylinder; of equal ones the first by its stored normal,
        # then origin (independent of argument order).
        key = lambda s: (*frame_axes(s.frame)[3], *frame_axes(s.frame)[0])
        if F(s1.radius) != F(s2.radius):
            first = F(s1.radius) < F(s2.radius)
        else:
            first = key(s1) <= key(s2)
        ruled, other = (s1, s2) if first else (s2, s1)
        o, a = ruled.axes()
        oo, ao = other.axes()
        m = cross(a, ao)
        # x towards the other axis along the common normal.
        if dot(sub(oo, o), m) < 0:
            m = scale(m, -1)
        fr = frame(o, a, m)
        return cls, Curve(fr, ruled.radius, cylinder_form(oo, ao, other.radius))
    if (s1.kind, s2.kind) == ('cylinder', 'cone'):
        o, a = s1.axes()
        ok, ak = s2.axes()
        if zero(cross(a, ak)) and zero(cross(sub(ok, o), a)):
            return None        # coaxial: S7a
        m = cross(a, ak)
        if zero(m):
            # Parallel axes: x towards the cone's axis.
            w = sub(ok, o)
            m = sub(w, scale(a, dot(w, a)/dot(a, a)))
        elif dot(sub(ok, o), m) < 0:
            m = scale(m, -1)
        cv = Curve(frame(o, a, m), s1.radius, cone_form(ok, ak, s2.radius, s2.angle))
        return 'by_roots', cv
    if (s1.kind, s2.kind) == ('cone', 'sphere'):
        o, a = s1.axes()
        c, _ = s2.axes()
        w = sub(c, o)
        w = sub(w, scale(a, dot(w, a)/dot(a, a)))
        if zero(w):
            return None        # coaxial: S7a
        cv = ConeCurve(s1, s2)
        if F(s1.radius) == 0 and dot(sub(c, o), sub(c, o)) == F(s2.radius)**2:
            return None        # the apex on the sphere: not yet parameterised
        return 'by_roots', cv
    if (s1.kind, s2.kind) == ('cylinder', 'sphere'):
        o, a = s1.axes()
        c, _ = s2.axes()
        w = sub(c, o)
        w = sub(w, scale(a, dot(w, a)/dot(a, a)))
        if zero(w):
            return None        # coaxial: S7a
        cls = classify_sphere(s1, s2)
        return cls, Curve(frame(o, a, w), s1.radius, sphere_form(c, s2.radius))
    return None


def rows(s1, s2):
    if 'torus' in (s1.kind, s2.kind):
        if s1.kind == 'torus' and s2.kind != 'torus':
            s1, s2 = s2, s1
        return torus_rows(s2, s1)
    found = curve(s1, s2)
    if found is None:
        return None
    cls, cv = found
    if cls == 'empty':
        return ['empty']
    roots = cv.roots()
    if cls == 'by_roots':
        # Components from the roots of D (all simple): loops where D > 0
        # between consecutive roots, two rings where D > 0 throughout.
        if not roots:
            if cv.D(mp.mpf(0)) < 0:
                return ['empty']
            return [('rings', cv.point(mp.mpf(0), 1), cv.point(mp.mpf(0), -1),
                     cv.point(mp.pi, 1), cv.point(mp.pi, -1))]
        out = []
        n = len(roots)
        for k in range(n):
            u0, u1 = roots[k], roots[(k+1) % n]
            if k+1 == n:
                u1 += 2*mp.pi
            mid = (u0+u1)/2
            if cv.D(mid) > 0:
                out.append(('loop', u0, u1, cv.point(u0, 1), cv.point(u1, 1), cv.point(mid, 1),
                            cv.point(mid, -1)))
        # Canonical order: a loop across pi first by its start below -pi
        # would repeat; normalise starts to (-pi, pi].
        return sorted(out, key=lambda r: float(r[1]))
    if cls == 'point':
        # D touches zero at u = 0 (towards the other surface).
        assert not roots, roots
        return [('point', cv.point(mp.mpf(0), 1))]
    if cls in ('rings', 'figure_eight'):
        assert not roots, roots
        if cls == 'rings':
            return [('rings', cv.point(mp.mpf(0), 1), cv.point(mp.mpf(0), -1),
                     cv.point(mp.pi, 1), cv.point(mp.pi, -1))]
        # The node where D touches zero: at u = 0 or pi, whichever has the
        # smaller D (the other is positive).
        node_u = mp.mpf(0) if abs(cv.D(mp.mpf(0))) < abs(cv.D(mp.pi)) else mp.pi
        return [('figure_eight', cv.point(node_u, 1), cv.point(node_u+mp.pi, 1),
                 cv.point(node_u+mp.pi, -1))]
    assert len(roots) == 2, (cls, roots)
    # The loop's range: where D >= 0, from one root to the other.
    u0, u1 = roots
    if cv.D((u0+u1)/2) < 0:
        u0, u1 = u1, u0+2*mp.pi
    mid = (u0+u1)/2
    return [('loop', u0, u1, cv.point(u0, 1), cv.point(u1, 1), cv.point(mid, 1), cv.point(mid, -1))]


def number(x):
    # Exact zeros print as zero, not as the 80-digit rounding around them.
    x = mpf(x)
    return '0.0' if abs(x) < mp.mpf(10)**-60 else mp.nstr(x, 30, min_fixed=-5, max_fixed=5)


def text(row):
    if isinstance(row, str):
        return row
    words = [row[0]]
    for part in row[1:]:
        words += [number(v) for v in part] if isinstance(part, (list, tuple)) else [number(part)]
    return ' '.join(words)


# ------------------------------------------------------------ S7b.3a: tori

def surd_sign(u, v, q):
    """The sign of u + v sqrt(q) (q >= 0, u, v, q rational), exactly."""
    su = (u > 0)-(u < 0)
    sv = ((v > 0)-(v < 0)) if q != 0 else 0
    if sv == 0 or su == sv:
        return su or sv
    if su == 0:
        return sv
    lhs, rhs = u*u, v*v*q
    return su if lhs > rhs else sv if lhs < rhs else 0


class TorusCurve:
    """A torus's meridians: C(phi) + r (cos t e + sin t a), C = o + R e,
    e(phi) = cos phi x + sin phi y. The other surface restricted to a
    meridian circle is f0 + alpha cos t + beta sin t (a plane or a sphere),
    found here by evaluating it at t = 0, pi/2 and pi."""
    def __init__(self, torus, f, toward):
        o, a = torus.axes()
        self.o, self.a, self.x, self.y = frame(o, a, toward)
        self.R, self.r = mpf(torus.radius), mpf(torus.minor)
        self.f = f

    def e(self, phi):
        c, s = mp.cos(phi), mp.sin(phi)
        return [c*x+s*y for x, y in zip(self.x, self.y)]

    def meridian(self, phi, t):
        e = self.e(phi)
        return [o+(self.R+self.r*mp.cos(t))*ei+self.r*mp.sin(t)*ai for o, ei, ai in zip(self.o, e, self.a)]

    def coefficients(self, phi):
        g0, g1, g2 = (self.f(self.meridian(phi, t)) for t in (mp.mpf(0), mp.pi/2, mp.pi))
        f0 = (g0+g2)/2
        return f0, (g0-g2)/2, g1-f0

    def D(self, phi):
        f0, alpha, beta = self.coefficients(phi)
        return alpha*alpha+beta*beta-f0*f0

    def point(self, phi, sign):
        f0, alpha, beta = self.coefficients(phi)
        rho = mp.sqrt(alpha*alpha+beta*beta)
        k = max(-1, min(1, -f0/rho))
        return self.meridian(phi, mp.atan2(beta, alpha)+sign*mp.acos(k))

    def roots(self, samples=2880):
        return Curve.roots(self, samples)


def torus_quadratic(torus, other):
    """D as an exact quadratic (P2, P1, P0) in c = m cos phi, with m^2 = q:
    the plane's normal, or the direction to the sphere's centre, has the
    component of length m normal to the axis; also that rational
    direction. Found by evaluating the exact D(c) at c = -1, 0, 1."""
    ot, A = torus.axes()
    R, r = F(torus.radius), F(torus.minor)
    AA = dot(A, A)
    if other.kind == 'plane':
        op, N = other.axes()
        d0 = dot(N, sub(ot, op))
        beta2 = r*r*dot(N, A)**2/AA
        q = dot(N, N)-dot(N, A)**2/AA
        toward = sub(N, scale(A, dot(N, A)/AA))
        D = lambda c: r*r*c*c+beta2-(d0+R*c)**2
    else:
        s, _ = other.axes()
        w = sub(s, ot)
        wa2 = dot(w, A)**2/AA
        q = dot(w, w)-wa2
        toward = sub(w, scale(A, dot(w, A)/AA))
        rho2 = F(other.radius)**2
        D = lambda c: (2*r*(R-c))**2+4*r*r*wa2-(R*R-2*R*c+dot(w, w)+r*r-rho2)**2
    P0 = D(F(0))
    P2 = (D(F(1))+D(F(-1)))/2-P0
    P1 = (D(F(1))-D(F(-1)))/2
    return (P2, P1, P0), q, toward


def torus_class(P, q):
    """The components as parameter ranges of phi, from exact signs:
    ('empty',), ('points', [phi]), ('loops', [(u0, u1)]), ('rings',),
    ('figure_eight', node). P's leading coefficient is negative; its
    discriminant is positive except for a sphere centred in the equatorial
    plane with K = 2 R^2, where P = P2 (c - R)^2: empty when R > m, the
    excluded meridian circle otherwise."""
    P2, P1, P0 = P
    assert P2 < 0, P
    if P1*P1-4*P2*P0 == 0:
        v = -P1/(2*P2)
        assert v*v > q, P
        return ('empty',)
    assert P1*P1-4*P2*P0 > 0, P
    at_plus = surd_sign(P2*q+P0, P1, q)     # D at phi = 0
    at_minus = surd_sign(P2*q+P0, -P1, q)   # D at phi = pi
    v = -P1/(2*P2)                           # the vertex
    v_plus = surd_sign(v, -1, q)            # sign of v - m
    v_minus = surd_sign(v, 1, q)            # sign of v + m
    m = mp.sqrt(mpf(q))
    sq = mp.sqrt(mpf(P1*P1-4*P2*P0))
    lo, hi = (-mpf(P1)+sq)/(2*mpf(P2)), (-mpf(P1)-sq)/(2*mpf(P2))
    t = lambda c: mp.acos(c/m)
    pi, zero_ = mp.pi, mp.mpf(0)
    if at_plus > 0 and at_minus > 0:
        return ('rings',)
    if at_plus > 0 and at_minus < 0:
        return ('loops', [(-t(lo), t(lo))])
    if at_plus < 0 and at_minus > 0:
        return ('loops', [(t(hi), 2*pi-t(hi))])
    if at_plus > 0 and at_minus == 0:
        return ('figure_eight', pi)
    if at_plus == 0 and at_minus > 0:
        return ('figure_eight', zero_)
    if at_plus == 0 and at_minus == 0:
        return ('loops', [(zero_, pi), (pi, 2*pi)])
    if at_plus < 0 and at_minus < 0:
        if v_plus < 0 and v_minus > 0:
            return ('loops', [(-t(lo), -t(hi)), (t(hi), t(lo))])
        return ('empty',)
    if at_plus == 0:
        # A root at m: the upper one (two loops touching at 0) or the lower
        # (a tangent point at 0).
        return ('loops', [(-t(lo), zero_), (zero_, t(lo))]) if v_plus < 0 else ('points', [zero_])
    return ('loops', [(t(hi), pi), (pi, 2*pi-t(hi))]) if v_minus > 0 else ('points', [pi])


def circle_circle(c1, r1, c2, r2, d2):
    """Two circles in a plane (centres as (rho, z) in mp, radii and the
    squared distance of the centres exact): their meeting points, by exact
    classes (tangent: one point)."""
    s, f = (r1+r2)**2, (r1-r2)**2
    if d2 > s or d2 < f:
        return []
    c1, c2 = [mpf(x) for x in c1], [mpf(x) for x in c2]
    d = mp.sqrt(mpf(d2))
    along = (mpf(d2)+mpf(r1)**2-mpf(r2)**2)/(2*d)
    u = [(b-a)/d for a, b in zip(c1, c2)]
    foot = [a+along*x for a, x in zip(c1, u)]
    if d2 == s or d2 == f:
        return [foot]
    h = mp.sqrt(mpf(r1)**2-along*along)
    n = [-u[1], u[0]]
    return [[a+h*x for a, x in zip(foot, n)], [a-h*x for a, x in zip(foot, n)]]


def torus_rows(torus, other):
    """Canonical rows of a torus and a plane or sphere (any position), or a
    coaxial cylinder, cone or torus; None for S7b.3b's pairs."""
    special = torus_special(torus, other)
    if special != 'general':
        return special
    P, q, toward = torus_quadratic(torus, other)
    return torus_components(torus_curve(torus, other), torus_class(P, q))


def torus_special(torus, other):
    """The rows of the special and coaxial cases, None for S7b.3b's pairs,
    'general' for a plane or sphere in general position."""
    ot, A = torus.axes()
    R, r = F(torus.radius), F(torus.minor)
    AA = dot(A, A)
    la = mp.sqrt(mpf(AA))
    am = [mpf(v)/la for v in A]
    om = [mpf(v) for v in ot]
    at_height = lambda z: tuple(o+z*a for o, a in zip(om, am))
    # (rho, z) meridian points to circles about the axis.
    circles = lambda pts: ana.canonical([ana.circle(at_height(z), A, rho) for rho, z in pts])
    meridian = [R, 0]
    oo, ao = other.axes()
    rel = sub(oo, ot)
    on_axis = zero(cross(rel, A))
    if other.kind == 'plane':
        N = ao
        if zero(cross(N, A)):
            # Normal to the axis: sin t = z / r at the plane's height z.
            z2 = dot(N, rel)**2*AA/dot(N, A)**2
            zm = mpf(dot(N, rel)*AA/dot(N, A))/la
            if z2 > r*r:
                return ['empty']
            if z2 == r*r:
                return circles([(mpf(R), zm)])
            h = mp.sqrt(mpf(r*r-z2))
            return circles([(mpf(R)+h, zm), (mpf(R)-h, zm)])
        if dot(N, A) == 0 and dot(N, rel) == 0:
            # The plane contains the axis: two meridian circles.
            u = unit(cross(A, N))
            return ana.canonical([ana.circle(tuple(o+s*mpf(R)*x for o, x in zip(om, u)), N, r)
                                  for s in (1, -1)])
    elif other.kind == 'sphere':
        if on_axis:
            za2, za = dot(rel, A)**2/AA, mpf(dot(rel, A))/la
            pts = circle_circle(meridian, r, [0, za], F(other.radius), R*R+za2)
            return circles(pts) if pts else ['empty']
        w2 = dot(rel, rel)
        if dot(rel, A) == 0 and F(other.radius)**2 == r*r+w2-R*R and w2 >= R*R:
            return ['not_conic']       # the sphere contains a meridian circle
    else:
        if not (zero(cross(ao, A)) and on_axis):
            return None                # S7b.3b
        za2, za = dot(rel, A)**2/AA, mpf(dot(rel, A))/la
        if other.kind == 'cylinder':
            k = F(other.radius)-R
            if k*k > r*r:
                return ['empty']
            if k*k == r*r:
                return circles([(mpf(other.radius), mp.mpf(0))])
            h = mp.sqrt(mpf(r*r-k*k))
            return circles([(mpf(other.radius), h), (mpf(other.radius), -h)])
        if other.kind == 'torus':
            R2, r2 = F(other.radius), F(other.minor)
            if R2 == R and r2 == r and za2 == 0:
                return ['same']
            if R2 == R and za2 == 0:
                return ['empty']
            pts = circle_circle(meridian, r, [R2, za], r2, (R2-R)**2+za2)
            return circles(pts) if pts else ['empty']
        # A coaxial cone: rho = s (rc + g (z - za) tan h) on each nappe s,
        # g = +-1 as the cone's axis runs along the torus's or against it.
        tn = mp.tan(mpf(other.angle))
        g = 1 if dot(ao, A) > 0 else -1
        pts = []
        for s in (1, -1):
            k0, k1 = s*(mpf(other.radius)-g*za*tn), s*g*tn
            a2, a1, a0 = k1*k1+1, 2*k1*(k0-mpf(R)), (k0-mpf(R))**2-mpf(r)**2
            disc = a1*a1-4*a2*a0
            assert abs(disc) > mp.mpf(10)**-40, 'a cone is never exactly tangent'
            if disc > 0:
                for z in ((-a1+mp.sqrt(disc))/(2*a2), (-a1-mp.sqrt(disc))/(2*a2)):
                    pts.append((k0+k1*z, z))
        return circles(pts) if pts else ['empty']
    return 'general'


def plane_function(pl):
    o, n = pl.axes()
    om, nm = [mpf(v) for v in o], [mpf(v) for v in n]
    return lambda p: sum((x-y)*z for x, y, z in zip(p, om, nm))


def torus_curve(torus, other):
    """The TorusCurve of a torus and a plane or sphere in general position."""
    _, q, toward = torus_quadratic(torus, other)
    oo, _ = other.axes()
    f = plane_function(other) if other.kind == 'plane' else sphere_form(oo, other.radius)['f']
    return TorusCurve(torus, f, toward)


def torus_components(cv, cls):
    """Rows of a torus curve's components; the ends of loops that are not
    nodes are checked against D's sign changes."""
    kind = cls[0]
    if kind == 'empty':
        return ['empty']
    sampled = cv.roots()
    if kind == 'points':
        assert not sampled, sampled
        return sorted([('point', cv.point(u, 1)) for u in cls[1]], key=lambda row: [float(x) for x in row[1]])
    if kind == 'rings':
        assert not sampled, sampled
        return [('rings', cv.point(mp.mpf(0), 1), cv.point(mp.mpf(0), -1), cv.point(mp.pi, 1),
                 cv.point(mp.pi, -1))]
    if kind == 'figure_eight':
        assert not sampled, sampled
        node = cls[1]
        return [('figure_eight', cv.point(node, 1), cv.point(node+mp.pi, 1), cv.point(node+mp.pi, -1))]
    # An end shared by two loops (modulo a turn) is a node, where D touches
    # zero without a sign change; the others are simple roots.
    ends = [u for span in cls[1] for u in span]
    wrap = lambda u: mp.atan2(mp.sin(u), mp.cos(u))
    simple = [u for u in ends if sum(abs(wrap(u-w)) < mp.mpf(10)**-60 for w in ends) == 1]
    # Each simple end is a sign change of the independently evaluated D
    # (tiny loops escape the sampling), and every sampled root is an end.
    eps = mp.mpf(10)**-30
    for u in simple:
        assert cv.D(u-eps)*cv.D(u+eps) < 0, u
    for s_ in sampled:
        assert any(abs(wrap(u-s_)) < mp.mpf(10)**-40 for u in simple), (s_, simple)
    out = []
    for u0, u1 in cls[1]:
        mid = (u0+u1)/2
        out.append(('loop', u0, u1, cv.point(u0, 1), cv.point(u1, 1), cv.point(mid, 1), cv.point(mid, -1)))
    return sorted(out, key=lambda row: float(row[1]))
