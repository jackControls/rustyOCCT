#!/usr/bin/env python3
"""Independent reference for S7a of REVIEW_NOTES.md: intersections of
analytic surfaces whose result is empty, the same surface, points, lines or
a conic.

A surface is the exact point set its stored binary64 data define (the
frame as `Frame3::new` stores it, from `identity_reference.frame_axes`):
a plane through the origin with the normal; a cylinder, cone or sphere
about the line through the origin along the normal, normalised exactly,
with the stored radius and half-angle. Rational configurations are decided
in exact rationals (`Fraction`); a cone's cosine and sine are transcendental,
so its decisions use mpmath at 80 digits, and its degenerate cases (a plane
through a rational apex, containing the axis, perpendicular to it) are
decided by exact predicates.

A plane and a quadric are intersected by substituting a rational basis of
the plane, `p = P0 + s e1 + t e2` with `e1` along the quadric axis's
projection and `e2 = n x e1`, into the quadric's implicit equation: the
quadratic in `(s, t)` has no cross term, and completing its squares
classifies the conic and gives its canonical parameters. The kernel instead
uses each pair's closed form, so the two share no formula. Other pairs use
their own constructions below.

Canonical results (`items`): `point`, `line` (the point nearest the global
origin and a unit direction), `circle` (centre, unit normal, radius),
`ellipse` (centre, unit normal, unit major direction, semi-major, semi-minor)
and `hyperbola` (centre, unit normal, unit transverse direction,
semi-transverse, semi-conjugate; both branches). A unit direction's first
nonzero coordinate is positive. `empty`, `same` and `not_conic` stand alone.
"""
from fractions import Fraction as F
from dataclasses import dataclass

import mpmath as mp

import brep_reference  # noqa: F401 (sets mp.dps to 40 on import; imported first)
from identity_reference import frame_axes

mp.mp.dps = 80

KINDS = ('plane', 'cylinder', 'cone', 'sphere', 'torus')


@dataclass
class Surface:
    kind: str
    frame: tuple        # origin, normal, x hint: nine binary64 values
    radius: float = 0.0  # a torus's major radius
    angle: float = 0.0  # a cone's half-angle
    minor: float = 0.0  # a torus's minor radius

    def axes(self):
        """(origin, normal) as exact rationals of the stored frame."""
        o, _, _, n = frame_axes(self.frame)
        return tuple(F(v) for v in o), tuple(F(v) for v in n)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, k):
    return tuple(x*k for x in a)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def zero(v):
    return all(x == 0 for x in v)


def mpf(x):
    if isinstance(x, F):
        return mp.mpf(x.numerator)/x.denominator
    return mp.mpf(x)


def unit(v):
    """A canonical unit vector: the first nonzero coordinate positive. The
    sign is read from the exact vector when it is rational."""
    lead = next(x for x in v if x != 0)
    sign = 1 if lead > 0 else -1
    m = [mpf(x) for x in v]
    norm = mp.sqrt(sum(x*x for x in m))
    return tuple(sign*x/norm for x in m)


def nearest(point, direction):
    """The point of a line nearest the global origin."""
    p = [mpf(x) for x in point]
    d = [mpf(x) for x in direction]
    k = sum(a*b for a, b in zip(p, d))/sum(x*x for x in d)
    return tuple(a-k*b for a, b in zip(p, d))


def line(point, direction):
    return ('line', nearest(point, direction), unit(direction))


def circle(centre, normal, radius):
    return ('circle', tuple(mpf(x) for x in centre), unit(normal), mpf(radius))


def basis(n, axis):
    """A rational basis of the plane with normal `n`: `e1` along the
    projection of `axis` (any rational in-plane direction when the axis is
    normal to the plane), `e2 = n x e1`."""
    e1 = sub(axis, scale(n, dot(axis, n)/dot(n, n))) if axis is not None else (0, 0, 0)
    if zero(e1):
        for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            e1 = cross(n, e)
            if not zero(e1):
                break
    return e1, cross(n, e1)


def quadric(q):
    """f(p) = 0 for a quadric: a function of a point (rationals or mpf)
    returning its value, exact for a sphere or cylinder."""
    o, a = q.axes()
    if q.kind == 'sphere':
        r2 = F(q.radius)**2
        return lambda p: dot(sub(p, o), sub(p, o))-r2
    if q.kind == 'cylinder':
        r2 = F(q.radius)**2
        aa = dot(a, a)
        return lambda p: dot(sub(p, o), sub(p, o))-dot(sub(p, o), a)**2/aa-r2
    # Both nappes: rho^2 cos^2 - (r cos + h sin)^2, h along the unit axis.
    ca, sa = mp.cos(mpf(q.angle)), mp.sin(mpf(q.angle))
    r = mpf(q.radius)
    la = mp.sqrt(mpf(dot(a, a)))
    am = [mpf(x) for x in a]
    om = [mpf(x) for x in o]

    def f(p):
        rel = [mpf(x)-y for x, y in zip(p, om)]
        h = sum(x*y for x, y in zip(rel, am))/la
        rho2 = sum(x*x for x in rel)-h*h
        return rho2*ca*ca-(r*ca+h*sa)**2
    return f


def plane_quadric(pl, q):
    """A plane and a sphere, cylinder or cone, by substitution."""
    p0o, n = pl.axes()
    qo, qa = q.axes()
    axis = None if q.kind == 'sphere' else qa
    e1, e2 = basis(n, axis)
    # P0: where the axis meets the plane, else the foot of the quadric's
    # origin on the plane (rational either way).
    if axis is not None and dot(qa, n) != 0:
        lam = dot(sub(p0o, qo), n)/dot(qa, n)
        P0 = add(qo, scale(qa, lam))
    else:
        P0 = sub(qo, scale(n, dot(sub(qo, p0o), n)/dot(n, n)))
    f = quadric(q)
    exact = q.kind != 'cone'
    # g(s, t) = A s^2 + C t^2 + 2 D s + 2 E t + G from five samples (the
    # cross term vanishes by the choice of e1 and e2).
    at = lambda s, t: f(add(P0, add(scale(e1, s), scale(e2, t)))) if exact else \
        f(tuple(mpf(x)+mpf(s)*mpf(y)+mpf(t)*mpf(z) for x, y, z in zip(P0, e1, e2)))
    one = F(1) if exact else mp.mpf(1)
    G = at(0*one, 0*one)
    A = (at(one, 0*one)+at(-one, 0*one))/2-G
    C = (at(0*one, one)+at(0*one, -one))/2-G
    D = (at(one, 0*one)-at(-one, 0*one))/4
    E = (at(0*one, one)-at(0*one, -one))/4
    B = (at(one, one)-at(one, 0*one)-at(0*one, one)+G)/2
    assert (B == 0) if exact else abs(B) < mp.mpf(10)**-60, 'a cross term survived'
    # Exact predicates: a circle when the section is normal to the axis;
    # a section through the apex only through a rational one (a cone of
    # radius zero at its origin) or when the plane contains the axis.
    circular = axis is None or zero(cross(qa, n))
    apex = q.kind == 'cone' and dot(sub(qo, p0o), n) == 0 and (F(q.radius) == 0 or dot(qa, n) == 0)
    return conic(pl, P0, e1, e2, A, C, D, E, G, exact, circular, apex)


def is0(x, exact):
    return x == 0 if exact else abs(x) < mp.mpf(10)**-60


def conic(pl, P0, e1, e2, A, C, D, E, G, exact, circular, apex):
    """Classify A s^2 + C t^2 + 2 D s + 2 E t + G = 0 in the plane's basis;
    `circular` and `apex` are the exact predicates of the pair."""
    _, n = pl.axes()
    pt = lambda s, t: tuple(mpf(x)+mpf(s)*mpf(y)+mpf(t)*mpf(z) for x, y, z in zip(P0, e1, e2))
    if is0(A, exact) and is0(C, exact):
        # Linear: a line or nothing (a plane/cylinder pair cannot be the same).
        raise AssertionError('a linear section of a quadric')
    if is0(A, exact) or is0(C, exact):
        # One direction free: parallel lines along it (a cylinder parallel
        # to the plane); the parabola cannot occur exactly (REVIEW_NOTES S7).
        if is0(A, exact):
            free, k, lin, other, e_free, e_k = 'e1', C, D, E, e1, e2
        else:
            free, k, lin, other, e_free, e_k = 'e2', A, E, D, e2, e1
        assert is0(lin, exact), 'a parabola'
        # k t^2 + 2 other t + G = 0 along e_k.
        disc = other*other-k*G
        if (disc < 0) if exact else disc < -mp.mpf(10)**-60:
            return ['empty']
        if is0(disc, exact):
            t0 = -other/k
            base = pt(0, t0) if free == 'e1' else pt(t0, 0)
            return [line(base, e_free)]
        root = mp.sqrt(mpf(disc))
        out = []
        for sign in (1, -1):
            t0 = (-mpf(other)+sign*root)/mpf(k)
            base = pt(0, t0) if free == 'e1' else pt(t0, 0)
            out.append(line(base, e_free))
        return out
    s0, t0 = -D/A, -E/C
    K = G-D*D/A-E*E/C
    # A (s - s0)^2 + C (t - t0)^2 + K = 0.
    centre = pt(s0, t0)
    # K vanishes exactly for a rational section, and for a cone exactly when
    # the plane passes through its apex.
    vanishes = K == 0 if exact else apex
    if not exact:
        assert (abs(K) < mp.mpf(10)**-60) == apex, 'the apex predicate disagrees'
    if (A > 0) == (C > 0):
        if vanishes:
            return [('point', centre)]
        if (K > 0) == (A > 0):
            return ['empty']
        a2, c2 = -K/A, -K/C
        l1, l2 = mpf(dot(e1, e1)), mpf(dot(e2, e2))
        # Semi-axes in length: sqrt(-K/A) along e1 scaled by |e1|, likewise e2.
        s1, s2 = mp.sqrt(mpf(a2)*l1), mp.sqrt(mpf(c2)*l2)
        if circular:
            return [circle(centre, n, (s1+s2)/2)]
        if s1 > s2:
            return [('ellipse', centre, unit(n), unit(e1), s1, s2)]
        return [('ellipse', centre, unit(n), unit(e2), s2, s1)]
    if vanishes:
        # Two lines through the centre: A s^2 = -C t^2.
        k = mp.sqrt(-mpf(A)/mpf(C))
        out = []
        for sign in (1, -1):
            d = tuple(mpf(x)+sign*k*mpf(y) for x, y in zip(e1, e2))
            out.append(('line', nearest(centre, d), unit_mp(d)))
        return out
    # Hyperbola: the transverse axis is the one whose coefficient has the
    # sign opposite to K.
    l1, l2 = mpf(dot(e1, e1)), mpf(dot(e2, e2))
    ax1, ax2 = mp.sqrt(abs(mpf(-K/A))*l1), mp.sqrt(abs(mpf(-K/C))*l2)
    if (mpf(-K/A) > 0):
        return [('hyperbola', centre, unit(n), unit(e1), ax1, ax2)]
    return [('hyperbola', centre, unit(n), unit(e2), ax2, ax1)]


def unit_mp(v):
    """A canonical unit vector of an irrational direction."""
    lead = next(x for x in v if abs(x) > mp.mpf(10)**-60)
    sign = 1 if lead > 0 else -1
    norm = mp.sqrt(sum(x*x for x in v))
    return tuple(sign*x/norm for x in v)


def plane_plane(p1, p2):
    o1, n1 = p1.axes()
    o2, n2 = p2.axes()
    d = cross(n1, n2)
    if zero(d):
        return ['same'] if dot(sub(o2, o1), n1) == 0 else ['empty']
    # The point alpha n1 + beta n2 on both planes.
    a11, a12, a22 = dot(n1, n1), dot(n1, n2), dot(n2, n2)
    b1, b2 = dot(o1, n1), dot(o2, n2)
    det = a11*a22-a12*a12
    alpha, beta = (b1*a22-b2*a12)/det, (a11*b2-a12*b1)/det
    return [line(add(scale(n1, alpha), scale(n2, beta)), d)]


def sphere_sphere(s1, s2):
    c1, _ = s1.axes()
    c2, _ = s2.axes()
    r1, r2 = F(s1.radius), F(s2.radius)
    d = sub(c2, c1)
    dd = dot(d, d)
    if dd == 0:
        return ['same'] if r1 == r2 else ['empty']
    if dd > (r1+r2)**2 or dd < (r1-r2)**2:
        return ['empty']
    # The radical plane: centre c1 + d k, k = (dd + r1^2 - r2^2) / (2 dd).
    k = (dd+r1*r1-r2*r2)/(2*dd)
    centre = add(c1, scale(d, k))
    rad2 = r1*r1-k*k*dd
    if rad2 == 0:
        return [('point', tuple(mpf(x) for x in centre))]
    return [circle(centre, d, mp.sqrt(mpf(rad2)))]


def parallel_cylinders(c1, c2):
    o1, a1 = c1.axes()
    o2, a2 = c2.axes()
    r1, r2 = F(c1.radius), F(c2.radius)
    # The offset between the axes, perpendicular to them.
    w = sub(o2, o1)
    w = sub(w, scale(a1, dot(w, a1)/dot(a1, a1)))
    dd = dot(w, w)
    if dd == 0:
        return ['same'] if r1 == r2 else ['empty']
    if dd > (r1+r2)**2 or dd < (r1-r2)**2:
        return ['empty']
    # In the plane across the axes: as two circles meeting.
    k = (dd+r1*r1-r2*r2)/(2*dd)
    foot = add(o1, scale(w, k))
    h2 = r1*r1-k*k*dd
    if h2 == 0:
        return [line(foot, a1)]
    v = cross(a1, w)
    h = mp.sqrt(mpf(h2)/mpf(dot(v, v)))
    return [line(tuple(mpf(x)+s*h*mpf(y) for x, y in zip(foot, v)), a1) for s in (1, -1)]


def crossing_equal_cylinders(c1, c2):
    """Equal radii, axes crossing at one point: two ellipses in the planes
    through the crossing point that bisect the axes."""
    o1, a1 = c1.axes()
    o2, a2 = c2.axes()
    # The crossing point o1 + a1 u = o2 + a2 v (exact: the axes are coplanar).
    w = sub(o2, o1)
    m = cross(a1, a2)
    u = dot(cross(w, a2), m)/dot(m, m)
    x = add(o1, scale(a1, u))
    r = mpf(c1.radius)
    l1, l2 = mp.sqrt(mpf(dot(a1, a1))), mp.sqrt(mpf(dot(a2, a2)))
    u1 = [mpf(v)/l1 for v in a1]
    u2 = [mpf(v)/l2 for v in a2]
    out = []
    for sign in (1, -1):
        # Each section lies in the plane through the crossing point spanned
        # by the common normal m and a bisector b = u1 + sign u2; its normal
        # is u1 - sign u2. Semi-minor r along m; semi-major r over the
        # cosine between an axis and that normal, along b.
        b = [p+sign*q for p, q in zip(u1, u2)]
        normal = [p-sign*q for p, q in zip(u1, u2)]
        cos_axis = abs(sum(p*q for p, q in zip(u1, normal)))/mp.sqrt(sum(v*v for v in normal))
        out.append(('ellipse', tuple(mpf(v) for v in x), unit_mp(normal), unit_mp(b), r/cos_axis, r))
    return out


def coaxial(q1, q2):
    """Coaxial cylinders, cones and spheres (a sphere centred on the
    other's axis): circles where their radius functions of the axial
    coordinate agree."""
    o, a = (q1 if q1.kind != 'sphere' else q2).axes()
    la = mp.sqrt(mpf(dot(a, a)))
    am = [mpf(v)/la for v in a]
    om = [mpf(v) for v in o]

    def radius2(q):
        """rho^2 as a polynomial in h (along the unit axis from o):
        coefficients (c0, c1, c2)."""
        qo, _ = q.axes()
        h0 = sum((mpf(x)-y)*z for x, y, z in zip(qo, om, am))
        if q.kind == 'cylinder':
            return (mpf(q.radius)**2, mp.mpf(0), mp.mpf(0))
        if q.kind == 'sphere':
            # rho^2 = R^2 - (h - h0)^2.
            return (mpf(q.radius)**2-h0*h0, 2*h0, mp.mpf(-1))
        # Cone: rho = r + (h - h0) tan a (both nappes: rho^2).
        t = mp.tan(mpf(q.angle))
        c = mpf(q.radius)-h0*t
        return (c*c, 2*c*t, t*t)
    p1, p2 = radius2(q1), radius2(q2)
    diff = [x-y for x, y in zip(p1, p2)]
    tiny = mp.mpf(10)**-60
    if all(abs(x) < tiny for x in diff):
        return ['same']
    roots = []
    if abs(diff[2]) < tiny:
        if abs(diff[1]) < tiny:
            return ['empty']
        roots = [-diff[0]/diff[1]]
    else:
        disc = diff[1]**2-4*diff[2]*diff[0]
        if disc < -tiny:
            return ['empty']
        if abs(disc) < tiny:
            roots = [-diff[1]/(2*diff[2])]
        else:
            sq = mp.sqrt(disc)
            roots = [(-diff[1]+sq)/(2*diff[2]), (-diff[1]-sq)/(2*diff[2])]
    out = []
    for h in roots:
        rho2 = p1[0]+p1[1]*h+p1[2]*h*h
        centre = tuple(x+h*y for x, y in zip(om, am))
        if rho2 < -tiny:
            continue
        if abs(rho2) < tiny:
            out.append(('point', centre))
        else:
            out.append(('circle', centre, unit(a), mp.sqrt(rho2)))
    return out or ['empty']


def on_axis(point, o, a):
    return zero(cross(sub(point, o), a))


def intersect(s1, s2):
    """The canonical items of two analytic surfaces, or ['not_conic']."""
    if KINDS.index(s1.kind) > KINDS.index(s2.kind):
        s1, s2 = s2, s1
    k = (s1.kind, s2.kind)
    if k == ('plane', 'plane'):
        return plane_plane(s1, s2)
    if s1.kind == 'plane':
        return plane_quadric(s1, s2)
    if k == ('sphere', 'sphere'):
        return sphere_sphere(s1, s2)
    o1, a1 = s1.axes()
    o2, a2 = s2.axes()
    if s2.kind == 'sphere':
        # A cylinder or cone and a sphere centred on its axis.
        return coaxial(s1, s2) if on_axis(o2, o1, a1) else ['not_conic']
    parallel = zero(cross(a1, a2))
    if parallel and on_axis(o2, o1, a1):
        return coaxial(s1, s2)
    if k == ('cylinder', 'cylinder'):
        if parallel:
            return parallel_cylinders(s1, s2)
        coplanar = dot(sub(o2, o1), cross(a1, a2)) == 0
        if coplanar and F(s1.radius) == F(s2.radius):
            return crossing_equal_cylinders(s1, s2)
    return ['not_conic']


def number(x):
    return mp.nstr(mpf(x), 30, min_fixed=-5, max_fixed=5) if not isinstance(x, str) else x


def text(item):
    if isinstance(item, str):
        return item
    words = [item[0]]
    for part in item[1:]:
        if isinstance(part, tuple):
            words += [number(v) for v in part]
        else:
            words.append(number(part))
    return ' '.join(words)


def canonical(items):
    """Items in a canonical order: by kind, then by their numbers."""
    def key(item):
        if isinstance(item, str):
            return (item, ())
        flat = []
        for part in item[1:]:
            flat += list(part) if isinstance(part, tuple) else [part]
        return (item[0], tuple(float(v) for v in flat))
    return sorted(items, key=key)


def parse(block):
    """A case block: `case NAME`, then two `surface KIND o(3) n(3) x(3)
    [radius [angle]]` rows (a torus: `major minor`)."""
    lines = [l.split() for l in block.strip().splitlines()]
    name = lines[0][1]
    surfaces = []
    for w in lines[1:]:
        assert w[0] == 'surface'
        v = [float(x) for x in w[2:]]
        s = Surface(w[1], tuple(v[:9]))
        if w[1] in ('cylinder', 'sphere', 'cone'):
            s.radius = v[9]
        if w[1] == 'cone':
            s.angle = v[10]
        if w[1] == 'torus':
            s.radius, s.minor = v[9], v[10]
        surfaces.append(s)
    return name, surfaces
