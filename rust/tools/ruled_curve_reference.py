#!/usr/bin/env python3
"""Independent reference for S7b.4 of REVIEW_NOTES.md: two cones (not coaxial,
without a common apex), and a cone whose rational apex lies on a sphere, a
cylinder or the other cone.

The curve is parameterised on the first cone's rulings through its apex `V`
(the cone first by stored data when both are cones), `V + v d(u)`,
`d(u) = cos h a + sin h (cos u x + sin u y)`, `x` the other axis's component
normal to `a` (the direction to the other's origin or centre when parallel).
The other surface restricted to a ruling is `A(u) v^2 + 2 B(u) v + C`,
found here by evaluating it at `v = -1, 0, 1`. With `v = tan(t / 2)` the curve
and its points at infinity (`t = pi`) are the zero set of
`G(u, t) = ((A + C) + (C - A) cos t) / 2 + B sin t` on the torus
`(u, t)` mod `2 pi`: the scan, critical meridians and components of
`torus_curve_reference.py` apply unchanged. A component crossing `t = pi` is
unbounded.

When the apex lies on the other surface (`C = 0` exactly, a rational apex),
`G = sin(t / 2) (A sin(t / 2) + 2 B cos(t / 2))`: the apex for every ruling,
and the curve `tan psi = -2 B / A` (`psi = t / 2`), one point per ruling, a
ring through the apex where `B` vanishes. Two cones with parallel axes and
equal half-angles have `A = 0` identically: `G = cos(t / 2) (C cos(t / 2) +
2 B sin(t / 2))`, the common circle at infinity and the conic
`tan psi = -C / (2 B)` in their radical plane.

Rows, in canonical order:

* `empty`;
* `fold u t x y z` for each fold (sorted by `u` in `(-pi, pi]`, then `t`);
* `component folds w_u w_t infinite` for each component (the torus's
  windings, and its crossings of `t = pi`, the points at infinity: an
  unbounded component crosses them and may return, winding 0);
* `ring u t x y z` for each component without folds, by its point at
  `u = 0`;
* `line p d` rows (S7a's form) for two cones with one rational apex: their
  common generatrices, or `point` the apex alone;
* `apex x y z crossing|isolated`, and the apex curve's
  `component 0 1 w infinite`
  and `ring 1 psi x y z` (its point at `u = 1`, `psi` in `(-pi / 2, pi / 2]`)
  when the apex lies on the other surface.
"""
from fractions import Fraction as F

import mpmath as mp

import analytic_intersection_reference as ana
from analytic_intersection_reference import dot, mpf, scale, sub, zero
from procedural_intersection_reference import cone_form, cylinder_form, frame, sphere_form
import torus_curve_reference as tc

mp.mp.dps = 80


class RulingChart:
    deg = 1

    def __init__(self, cone, other):
        o, a = cone.axes()
        oo, ao = other.axes()
        toward = sub(ao, scale(a, dot(ao, a)/dot(a, a))) if other.kind != 'sphere' else (0, 0, 0)
        if zero(toward):
            w = sub(oo, o)
            toward = sub(w, scale(a, dot(w, a)/dot(a, a)))
        _, self.a, self.x, self.y = frame(o, a, toward)
        h = mpf(cone.angle)
        self.ch, self.sh = mp.cos(h), mp.sin(h)
        back = mpf(cone.radius)/mp.tan(h)
        self.V = [mpf(p)-back*q for p, q in zip(o, self.a)]
        if other.kind == 'sphere':
            self.f = sphere_form(oo, other.radius)['f']
        elif other.kind == 'cylinder':
            self.f = cylinder_form(oo, ao, other.radius)['f']
        else:
            self.f = cone_form(oo, ao, other.radius, other.angle)['f']
        self.C = self.f(self.V)

    def d(self, u):
        c, s = mp.cos(u), mp.sin(u)
        return [self.ch*a+self.sh*(c*x+s*y) for a, x, y in zip(self.a, self.x, self.y)]

    def AB(self, u):
        d = self.d(u)
        fp = self.f([v+e for v, e in zip(self.V, d)])
        fm = self.f([v-e for v, e in zip(self.V, d)])
        return (fp+fm)/2-self.C, (fp-fm)/4

    def G(self, u, t):
        A, B = self.AB(u)
        return ((A+self.C)+(self.C-A)*mp.cos(t))/2+B*mp.sin(t)

    def fourier(self, u):
        A, B = self.AB(u)
        return [(A+self.C)/2, (self.C-A)/2], [mp.mpf(0), B]

    def roots(self, u, loose=False):
        coeffs = tc.zpoly(*self.fourier(u))
        while coeffs and abs(coeffs[0]) < mp.mpf(10)**-60:
            coeffs = coeffs[1:]
        out = []
        for z in mp.polyroots(coeffs, maxsteps=200, extraprec=200):
            if abs(abs(z)-1) < (mp.mpf(10)**-12 if loose else mp.mpf(10)**-25):
                out.append(mp.atan2(mp.im(z), mp.re(z)))
        return sorted(out)

    def fold(self, u, t):
        g = lambda p, s: [self.G(p, s), mp.diff(lambda q: self.G(p, q), s)]
        p, s = mp.findroot(g, (u, t))
        return tc.wrap(p), tc.wrap(s)

    def point(self, u, t):
        v = mp.tan(t/2)
        return [p+v*q for p, q in zip(self.V, self.d(u))]


def order(s1, s2):
    """The chart's cone first: the only cone, or the first by stored data."""
    if s1.kind != 'cone':
        return s2, s1
    if s2.kind != 'cone':
        return s1, s2
    from identity_reference import frame_axes
    key = lambda s: (*frame_axes(s.frame)[3], *frame_axes(s.frame)[0])
    return (s1, s2) if key(s1) <= key(s2) else (s2, s1)


def apex_on(cone, other):
    """Whether the cone's apex lies exactly on the other surface: a rational
    apex (radius zero at the origin) on a sphere or a cylinder, by exact
    rationals; never on another cone (its cos^2 is transcendental, and a
    common apex is S7a's)."""
    if F(cone.radius) != 0 or other.kind == 'cone':
        return False
    V, _ = cone.axes()
    o, a = other.axes()
    w = sub(V, o)
    if other.kind == 'sphere':
        return dot(w, w) == F(other.radius)**2
    return dot(w, w)-dot(w, a)**2/dot(a, a) == F(other.radius)**2


def twins(cone, other):
    """Parallel axes and equal half-angles: A vanishes identically (every
    ruling of the first is parallel to a generatrix of the second), the
    circle at infinity is common, and the rest is the conic in the cones'
    radical plane, one point per ruling."""
    if other.kind != 'cone' or F(cone.angle) != F(other.angle):
        return False
    _, a1 = cone.axes()
    _, a2 = other.axes()
    return zero(ana.cross(a1, a2))


def common_apex(cone, other):
    """Two cones with one rational apex: their common generatrices."""
    if other.kind != 'cone' or F(cone.radius) != 0 or F(other.radius) != 0:
        return False
    return cone.axes()[0] == other.axes()[0]


def apex_lines(cv):
    """The rulings where A vanishes (a generatrix of both cones), as lines
    through the apex, or the apex alone."""
    n = 2880
    us = [-mp.pi+2*mp.pi*(k+mp.mpf(1)/3)/n for k in range(n)]
    As = [cv.AB(u)[0] for u in us]
    out = []
    for k in range(n):
        j = (k+1) % n
        if (As[k] > 0) != (As[j] > 0):
            hi = us[j] if j else us[j]+2*mp.pi
            u = mp.findroot(lambda x: cv.AB(x)[0], (us[k], hi), solver='anderson')
            out.append(ana.line(tuple(cv.V), tuple(cv.d(u))))
    return ana.canonical(out) if out else [('point', tuple(cv.V))]


def rows(s1, s2):
    cone, other = order(s1, s2)
    cv = RulingChart(cone, other)
    if common_apex(cone, other):
        return apex_lines(cv)
    if apex_on(cone, other):
        return apex_rows(cv)
    if twins(cone, other):
        return factor_rows(cv, lambda u: (2*cv.AB(u)[1], cv.C), None)
    phis, roots, folds, edges = tc.scan(cv)
    cycles = tc.components(phis, roots, edges) if edges else []
    if not cycles:
        return ['empty']
    out = []
    for u, t in sorted(folds, key=lambda x: (float(x[0]), float(x[1]))):
        assert abs(tc.wrap(t-mp.pi)) > mp.mpf(10)**-6, 'a fold at infinity'
        out.append(('fold', u, t, cv.point(u, t)))
    comps, rings = [], []
    for w, cycle in cycles:
        verts = set(cycle)
        nf = sum(1 for x, y in edges if x in verts and x[0] == y[0])
        # Crossings of t = pi (the points at infinity) along the cycle: a
        # step whose direct difference and shortest turn differ.
        ts = [roots[k][j] for k, j in cycle]
        infinite = sum(1 for a, b in zip(ts, ts[1:]+ts[:1]) if abs(b-a) > mp.pi)
        comps.append(('component', nf, w[0], w[1], infinite))
        if nf == 0:
            k0 = min(range(len(phis)), key=lambda k: abs(phis[k]))
            (j0,) = [j for k, j in cycle if k == k0]
            t0 = min(cv.roots(mp.mpf(0)), key=lambda t: tc.dist(t, roots[k0][j0]))
            rings.append(('ring', mp.mpf(0), t0, cv.point(mp.mpf(0), t0)))
    return out+sorted(comps, key=lambda r: r[1:])+sorted(rings, key=lambda r: float(r[2]))


def apex_rows(cv):
    """The apex on the other surface: `A sin psi + 2 B cos psi = 0`, one point
    per ruling; the apex is crossed where B vanishes."""
    n = 2880
    us = [-mp.pi+2*mp.pi*(k+mp.mpf(1)/3)/n for k in range(n)]
    bs = [cv.AB(u)[1] for u in us]
    crossings = sum(1 for k in range(n) if (bs[k] > 0) != (bs[(k+1) % n] > 0))
    kind = 'crossing' if crossings == 2 else 'isolated'
    assert crossings in (0, 2), crossings
    return factor_rows(cv, lambda u: cv.AB(u)[0:1]+(2*cv.AB(u)[1],), ('apex', cv.V, kind))


def factor_rows(cv, ab, head):
    """A curve `alpha(u) sin psi + beta(u) cos psi = 0` (`v = tan psi`), one
    point per ruling: `psi = atan2(-beta, alpha)` mod pi. Its winding in psi
    over a turn, its crossings of infinity (where alpha vanishes) and its point
    at `u = 1`."""
    n = 2880
    us = [-mp.pi+2*mp.pi*(k+mp.mpf(1)/3)/n for k in range(n)]
    psi = lambda u: mp.atan2(-ab(u)[1], ab(u)[0])
    total = mp.mpf(0)
    prev = psi(us[0])
    for u in us[1:]+[us[0]+2*mp.pi]:
        cur = psi(u)
        step = cur-prev
        step -= mp.pi*mp.nint(step/mp.pi)
        total += step
        prev = cur
    w = int(mp.nint(total/mp.pi))
    # The ring's point at u = 1 (u = 0 is often where it crosses the apex).
    u1 = mp.mpf(1)
    p0 = psi(u1)
    p0 = p0-mp.pi*mp.floor((p0+mp.pi/2)/mp.pi)
    if p0 <= -mp.pi/2:
        p0 += mp.pi
    point = [v+mp.tan(p0)*e for v, e in zip(cv.V, cv.d(u1))]
    # The curve meets infinity where alpha vanishes (psi = pi / 2).
    as_ = [ab(u)[0] for u in us]
    infinite = sum(1 for k in range(n) if (as_[k] > 0) != (as_[(k+1) % n] > 0))
    rows_ = [('component', 0, 1, w, infinite), ('ring', u1, p0, point)]
    return ([head] if head else [])+rows_


def curve_samples(s1, s2, every=6, reach=20):
    """Float points along each component of the reference curve within
    `reach` of the apex (unbounded components are cut there), and the apex
    when it is an isolated point, for comparisons with native lines."""
    cone, other = order(s1, s2)
    cv = RulingChart(cone, other)
    near = lambda p: sum((x-v)**2 for x, v in zip(p, cv.V)) < reach**2
    if common_apex(cone, other):
        rows_ = apex_lines(cv)
        if rows_[0][0] == 'point':
            return [], [[float(x) for x in cv.V]]
        comps = [[[float(p+k*d) for p, d in zip(row[1], row[2])] for k in range(-10, 11) if k] for row in rows_]
        return comps, []
    if apex_on(cone, other) or twins(cone, other):
        ab = (lambda u: (cv.AB(u)[0], 2*cv.AB(u)[1])) if apex_on(cone, other) else \
            (lambda u: (2*cv.AB(u)[1], cv.C))
        pts = []
        for k in range(720):
            u = -mp.pi+2*mp.pi*(k+mp.mpf(1)/3)/720
            a, b = ab(u)
            if abs(a) < mp.mpf(10)**-30:
                continue
            p = [v-(b/a)*e for v, e in zip(cv.V, cv.d(u))]
            if near(p):
                pts.append([float(x) for x in p])
        isolated = []
        if apex_on(cone, other) and apex_rows(cv)[0][2] == 'isolated':
            isolated.append([float(x) for x in cv.V])
        return [pts[::every]], isolated
    phis, roots, folds, edges = tc.scan(cv)
    cycles = tc.components(phis, roots, edges) if edges else []
    comps = []
    for _, cycle in cycles:
        pts = [cv.point(phis[k], roots[k][j]) for k, j in cycle[::every]]
        comps.append([[float(x) for x in p] for p in pts if near(p)])
    return comps, []
