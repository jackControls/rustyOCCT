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
        words += [number(v) for v in part] if isinstance(part, list) else [number(part)]
    return ' '.join(words)
