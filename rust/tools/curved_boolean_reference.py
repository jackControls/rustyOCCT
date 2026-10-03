#!/usr/bin/env python3
"""Independent reference for S9c.1 of REVIEW_NOTES.md: Booleans of two
prisms of line, arc and circle profiles (planar and cylindrical walls) in
any relative position, where every pair of faces meets in lines, circles or
ellipses: the prisms' axes parallel (so every pair of walls is parallel or
coaxial) or crossing, cylinders of crossing axes of equal radii.

Each prism is its construction's exact model: a point is `o + u x + v y + w
n` in rationals from the frame's stored binary64 origin and axes
(`stored_axes`, the kernel's `Frame3::new`, not exactly orthonormal), the
profile's binary64 points, centres and radii and the offsets' values (an
arc is the exact circle `|(u, v) - c| = r` between its end points). Nothing
here uses the kernel, a surface/surface intersection or a 3D arrangement.

* **Slicing (volume, first moments, solids).** Both solids are sliced by
  the planes `d . X = s` with `d = n_A x n_B` (crossing axes) or `d = n x e`
  for a world axis `e` (parallel axes), so every slice plane contains both
  axes' directions: in the plane's coordinates `X = o_A + s d' + lambda g +
  mu h` (`d . d' = 1`; `g = n_A` and `h = n_B`, or `g = n` and `h = n x
  d`), each prism's section is a union of parallelograms, one per chord of
  its profile along the plane's trace in its `(u, v)` (a cylinder cut along
  its generatrices), between its two cap lines: with the slicing direction
  so chosen no section holds an ellipse arc (slices in another direction
  would: a cylinder cut obliquely; here ellipses bound the planar faces'
  regions of the area part instead). A chord's ends are exact
  functions of `s`: `L(s) + k sqrt(Q(s))` (a line's crossing linear in `s`,
  a circle's a quadratic surd). Each slice's common is the union of the
  convex intersections of the two sections' parallelograms, its cut each
  object parallelogram less the tool's (successive half-plane complements)
  and its fuse the cut and the tool's; their areas and first moments by
  Green's theorem over the convex pieces' edges (segments: closed form).
  The section's structure changes only where three of the plane's boundary
  lines are concurrent or two parallel ones coincide, where the trace
  passes a profile vertex or touches a circle: all found as the real roots
  of exact polynomials (the surds squared away; squaring's spurious roots
  are filtered numerically: the meeting must lie on both sections' closed
  boundaries). Between consecutive breakpoints the integrands are analytic
  in `s` and the surds whose roots end the interval, so each interval is
  integrated by Gauss-Legendre after `s = a + (b - a)(1 - cos t)/2` (which
  makes end-point square roots analytic), doubling the rule (12 to 96
  nodes) until successive estimates agree within 1e-33 of the case's size
  to the fourth, else halving the interval. Solids (combinatorics only, in
  binary64): in each interval the plane's lines at its middle cut the plane
  into convex faces, each inside or outside each section; faces of the
  result are joined when adjacent across one line (they share a segment for
  the whole interval: a face of the solid), and across a breakpoint when
  their limit polygons there overlap in a positive area. Regions meeting
  along a curve or at a point are separate solids (the regularized
  Boolean).
* **Surface area.** Every face of each input (caps: the profile in `(u,
  v)`; flat walls: segment position and height; cylinder walls: angle and
  height) is swept by lines of its own parameterisation (`v`, or the
  height); along each line its crossings of the other solid's boundary
  (the other's caps, flat and cylindrical walls, in the other's
  coordinates) split it into pieces classified at their midpoints: inside
  or outside the other solid, or on a face of it (exactly coplanar planes
  or coincident cylinders, decided in rationals first) with the same or
  the opposite orientation. Coincident cylinders must come from frames with
  equal stored axes or exactly orthonormal ones: otherwise the two exact
  models are different elliptic cylinders within rounding of each other,
  refused. The class lengths are integrated over the
  face's parameter with the surface's own area element (`|x * y|`, `|e *
  n|`, `r |t(theta) * n|`), with breakpoints where two crossings (or a
  crossing and the face's own boundary) coincide, where a line passes the
  other profile's vertex or touches its circle: roots of exact polynomials
  in the parameter (in `tan(theta / 2)` on cylinders), filtered like the
  slices'. A result keeps: fuse the object's outside pieces, the tool's
  outside and the shared ones of the same orientation (once); cut the
  object's outside, the tool's inside and the shared of opposite
  orientation; common the object's and tool's inside and the shared of the
  same orientation.

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a and S9b
references. The exact polynomial helpers (Fractions, Yun's square-free
factorization) are `boolean_reference.py`'s; nothing else is shared with
another reference.

**S9f.1: spline walls.** A profile may also hold S8b's nonrational spline
segments (`identity_reference.Spline`), each cut into its exact Bezier
spans (blossoms, Fractions; `SplineSpan`, one element per span, its joins
the profile's vertices), when the other prism's profile has lines only (a
spline prism against a polyhedral one, S9f.1; either may be the object). A
span's wall is `o + X(tau) x + Y(tau) y + w n` over `tau` in [0, 1]
(`Face` kind `spl`, a polynomial family swept by its generatrices, whose
crossings of the other's planes are exact polynomials in `tau` of degree
`p`, as a cylinder's are in `tan(theta / 2)`). A line meets a span at the
real roots in [0, 1] of the degree-`p` polynomial `(S(tau) - b) x d`
(`span_roots`: monotone runs between the derivative's roots, then
safeguarded Newton steps at the working precision; never a resultant),
each tagged by its order along the span, which holds between breakpoints.
Where the S9c.1 events are coincidences of surds, a spline crossing is
implicit: an event between it and a line's crossing moving linearly with
the family parameter, `Q0 + x Q1 = S(tau)`, is the root of the exact
polynomial `(S(tau) - Q0) x Q1` (degree `p`), `x` then read off along
`Q1`; the family's tangencies to a span are the roots of `S'(tau) x d`
(degree `p - 1`), and its passages through a span's ends (knots) are the
profile's vertices'. A profile with splines takes its area and first
moments by the symmetric Green forms (`(x dy - y dx) / 2`, `x (x dy - y
dx) / 3`, `y (x dy - y dx) / 3`) on every element, arcs included, the
spans' exactly. For the solids' union-find, a spline chord at an
interval's end is the crossing there nearest to the chord just inside the
interval (its order along the span may change at a tangency or a knot).

With splines (`Pair(divergence=True)`, the default when either profile has
one) each face sweep also integrates, per class, `X . N / 3` and `x_i^2
N_i / 2` over its pieces (`N dA` the family's `X0'(x) x E dx dt` with the
face's outward sense), so each operation's volume and first moments are
computed a second way, by the divergence theorem over the kept pieces,
independently of the slicing (`Pair.divergence_volumes`).

**S9f.2a: parallel walls.** With exactly parallel axes the other prism may
hold arcs, circles or splines too: every wall is swept along the common
axis, so a spline's chord meets the other's arc or spline chord, and a
curved wall's generatrices change class, only where the slice's trace or
the generatrix passes a 2D crossing of the two profiles' curves projected
along the axis (`parallel_crossings`, in the exact affine map between the
frames' `(u, v)`): a span against a circle at the real roots of the exact
polynomial `|S(tau) - c|^2 - r^2` of degree `2 p` (its square-free factors:
a tangency is a root of even multiplicity), two spans by subdivision of
their Bezier forms in fractions while their control boxes meet, then
Newton's method at 40 digits (on the distance's gradient where the
crossing is tangent), as S9a.2's reference finds its meetings, never a
resultant. Each crossing's slice parameter and each face's parameter there
are breakpoints of the slicing and of the face sweeps.

**S9f.2b: crossing walls.** With crossing axes the other prism may hold
arcs and circles against a spline prism (no second spline). The slicing's
planes hold both axes, so a spline wall's section is still its chords'
parallelograms and a cylinder's its chords' (surds), and an event where a
spline chord, a cylinder's chord and a height line are concurrent is a
point of the two walls' meeting at a cap's height: along the spline wall's
ruling `X = o + S(tau) + w n` the cylinder's equation in its own frame is
`A w^2 + 2 B(tau) w + C(tau)` (degrees 0, `p`, `2 p`), so at the spline
prism's cap (`w` fixed) or the cylinder prism's (`w` linear in `S(tau)`)
the points are the real roots of an exact polynomial of degree `2 p` in
`tau` (`wall_cylinder_events`), `s` read from the point. A cylinder wall's
generatrices swept against a spline wall cross it implicitly; their events
are the same meetings at the face's own heights and the other's caps, and
the generatrices whose projection along the spline's axis touches a span
(at the span's points where it runs along the cylinder's axis's
projection, the cylinder met on that ruling at 40 digits). A cap's sweep
line meets the other's cylinder (a surd) on a spline cap edge, or a
spline wall on a cylinder prism's cap circle, at the same polynomials'
roots. The meeting's turning points (`wall_cylinder_turns`: the roots of
`B^2 - A C`, degree `2 p`) need no event of their own: the slices' combinatorics
change only at the trace's tangencies, which they are.
"""
from fractions import Fraction as F
import itertools

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

from identity_reference import Spline, stored
from curve_surface_reference import stored_axes
from boolean_reference import padd, pdeg, pder, peval, pgcd, pint, pmul, power_of, psub, ptrim, squarefree

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
SET = {
    'fuse': lambda a, b: a or b,
    'cut': lambda a, b: a and not b,
    'common': lambda a, b: a and b,
}

# ------------------------------------------------------------------ numbers and vectors


def M(v):
    if isinstance(v, F):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def Mv(v):
    return tuple(M(c) for c in v)


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, s):
    return tuple(x*s for x in a)


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), 0)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def cross2(a, b):
    return a[0]*b[1]-a[1]*b[0]


def det3(a, b, c):
    return dot(a, cross(b, c))


def apply(rows, v):
    return tuple(dot(r, v) for r in rows)


# ------------------------------------------------------------------ exact polynomials

def P(*c):
    return ptrim([F(x) for x in c]) if c else [F(0)]


def pscale(p, k):
    return ptrim([x*k for x in p]) if k else [F(0)]


def pzero(p):
    return pdeg(p) < 0


def pvec(vectors):
    """A polynomial vector from ascending vector coefficients."""
    return tuple(ptrim([F(v[i]) for v in vectors]) for i in range(len(vectors[0])))


def pvapply(rows, polys):
    """A constant matrix times a polynomial vector."""
    out = []
    for r in rows:
        acc = [F(0)]
        for k, p in zip(r, polys):
            acc = padd(acc, pscale(p, k))
        out.append(acc)
    return tuple(out)


def sq_poly(L, terms):
    """A polynomial vanishing wherever `L + sum(a sqrt(Q))` does, for every
    choice of the roots' signs (at most two square roots)."""
    terms = [(a, Q) for a, Q in terms if a != 0 and not pzero(Q)]
    if not terms:
        return ptrim(L)
    if len(terms) == 1:
        a, Q = terms[0]
        return psub(pmul(L, L), pscale(Q, a*a))
    assert len(terms) == 2, 'three square roots in one equation'
    (a, Q1), (b, Q2) = terms
    S = padd(psub(pmul(L, L), pscale(Q2, b*b)), pscale(Q1, a*a))
    return psub(pmul(S, S), pscale(pmul(pmul(L, L), Q1), 4*a*a))


ROOTS = {'polys': 0, 'roots': 0}


def real_roots(p):
    """The real roots (mpf) of an exact polynomial: each square-free factor
    by `mp.polyroots` at 60 digits; roots with an imaginary part below
    1e-12 count as real (a spurious breakpoint costs nothing)."""
    p = ptrim([F(c) for c in p])
    if pdeg(p) <= 0:
        return []
    ROOTS['polys'] += 1
    out = []
    for f, _ in squarefree(p):
        if pdeg(f) == 1:
            out.append(-M(f[0])/M(f[1]))
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
    ROOTS['roots'] += len(out)
    return out


class Desc:
    """A crossing parameter `(L(x) + s k sqrt(Q(x))) / den(x)` on a family
    of lines (`s` the branch: -1 or +1 when `k` is nonzero, else 0)."""

    def __init__(self, L, k=F(0), Q=None, key=None):
        self.L, self.k, self.Q, self.key = ptrim(L), F(k), ptrim(Q) if Q is not None else [F(0)], key
        if self.k == 0:
            self.Q = [F(0)]

    def branches(self):
        return (-1, 1) if self.k else (0,)

    def value(self, x, s, den):
        q = peval([M(c) for c in self.Q], x) if self.k else 0
        if self.k and q < 0:
            q = mp.mpf(0)
        return (peval([M(c) for c in self.L], x)+s*M(self.k)*mp.sqrt(q))/peval([M(c) for c in den], x)


def element_desc(el, Nb, den, d, key):
    """The crossing parameter of the line `Nb(x)/den(x) + t d` (uv) with a
    profile element, or None when it runs parallel to a segment (or the
    element is a spline span: its crossings are implicit, S9f.1)."""
    if el.kind == 'spline':
        return None
    if el.kind == 'seg':
        cd = cross2(d, el.e)
        if cd == 0:
            return None
        w = (psub(pscale(den, el.p[0]), Nb[0]), psub(pscale(den, el.p[1]), Nb[1]))
        L = psub(pscale(w[0], el.e[1]), pscale(w[1], el.e[0]))
        return Desc(pscale(L, 1/cd), key=key)
    a = d[0]*d[0]+d[1]*d[1]
    if a == 0:
        return None
    W = (psub(Nb[0], pscale(den, el.c[0])), psub(Nb[1], pscale(den, el.c[1])))
    bW = padd(pscale(W[0], d[0]), pscale(W[1], d[1]))
    Q = psub(pmul(bW, bW), pscale(psub(padd(pmul(W[0], W[0]), pmul(W[1], W[1])),
                                       pscale(pmul(den, den), el.r*el.r)), a))
    return Desc(pscale(bW, -1/a), 1/a, Q, key=key)


# ------------------------------------------------------------------ profiles

def angle_rel(phi, start):
    tau = 2*mp.pi
    return (phi-start) % tau


class Seg:
    kind = 'seg'

    def __init__(self, p, q, left, outer):
        self.p, self.q, self.left, self.outer = p, q, left, outer
        self.e = (q[0]-p[0], q[1]-p[1])
        self.pm, self.qm, self.em = Mv(p), Mv(q), Mv(self.e)
        self.len = mp.sqrt(self.em[0]**2+self.em[1]**2)


class Round:
    """An arc or a whole circle: its angles `start + [0, sweep]` (turning
    counter-clockwise), traversed from `p` to `q` (`ccw`)."""

    def __init__(self, c, r, full, radial_out, outer, p=None, q=None, ccw=True):
        self.kind = 'circle' if full else 'arc'
        self.c, self.r, self.radial_out, self.outer = c, r, radial_out, outer
        self.cm, self.rm = Mv(c), M(r)
        self.p, self.q, self.ccw = p, q, ccw
        if full:
            self.start, self.sweep = mp.mpf(0), 2*mp.pi
            self.a0, self.a1 = mp.mpf(0), 2*mp.pi
        else:
            for pt in (p, q):
                d2 = (pt[0]-c[0])**2+(pt[1]-c[1])**2
                assert d2 == r*r, f'arc end {pt} not on its circle exactly'
            a0 = mp.atan2(M(p[1]-c[1]), M(p[0]-c[0]))
            a1 = mp.atan2(M(q[1]-c[1]), M(q[0]-c[0]))
            tau = 2*mp.pi
            if ccw:
                self.start, self.sweep = a0, (a1-a0) % tau
                self.a0, self.a1 = a0, a0+self.sweep
            else:
                self.start, self.sweep = a1, (a0-a1) % tau
                self.a0, self.a1 = a0, a0-self.sweep
            assert self.sweep > 0

    def range_bad(self, phi):
        if self.kind == 'circle':
            return mp.mpf(0)
        rel = angle_rel(phi, self.start)
        if rel <= self.sweep:
            return mp.mpf(0)
        return min(rel-self.sweep, 2*mp.pi-rel)*self.rm

    def point(self, phi):
        return (self.cm[0]+self.rm*mp.cos(phi), self.cm[1]+self.rm*mp.sin(phi))


# ------------------------------------------------------------------ S9f.1: spline spans

SPAN_SLACK = mp.mpf(10)**-3


def _mp_trim(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def _mp_der(c):
    return [k*c[k] for k in range(1, len(c))] or [mp.mpf(0)]


def _newton(c, dc, a, b, fa):
    """The root of `c` in [a, b] where it changes sign (`fa` its sign at
    `a`), by Newton steps kept inside the shrinking bracket."""
    lo, hi = a, b
    t = (a+b)/2
    eps = mp.mpf(2)**(8-mp.mp.prec)
    for _ in range(400):
        f = peval(c, t)
        if f == 0:
            return t
        if (f > 0) == (fa > 0):
            lo = t
        else:
            hi = t
        d = peval(dc, t)
        nt = t-f/d if d != 0 else (lo+hi)/2
        if not lo < nt < hi:
            nt = (lo+hi)/2
        if abs(nt-t) <= eps*(1+abs(t)) or hi-lo <= eps*(1+abs(t)):
            return nt
        t = nt
    return t


def _roots_in(c, a, b, touch):
    """The real roots of the mpf polynomial `c` in [a, b], sorted: between
    consecutive roots of its derivative it is monotone. With `touch`, a
    local extreme within `touch` of zero (no sign change) is a double root,
    listed twice."""
    c = _mp_trim(c)
    n = len(c)-1
    if n <= 0:
        return []
    if n == 1:
        r = -c[0]/c[1]
        return [r] if a <= r <= b else []
    dc = _mp_der(c)
    crit = [t for t in _roots_in(dc, a, b, None) if a < t < b]
    pts = [a]+sorted(crit)+[b]
    vals = [peval(c, t) for t in pts]
    out = []
    for i in range(len(pts)-1):
        fa, fb = vals[i], vals[i+1]
        if fa == 0 and (not out or out[-1] != pts[i]):
            out.append(pts[i])
        if (fa < 0 < fb) or (fb < 0 < fa):
            out.append(_newton(c, dc, pts[i], pts[i+1], fa))
    if vals[-1] == 0 and (not out or out[-1] != pts[-1]):
        out.append(pts[-1])
    if touch is not None:
        for k in range(1, len(pts)-1):
            f = vals[k]
            if f != 0 and abs(f) <= touch and (vals[k-1]-f)*(vals[k+1]-f) > 0 and \
                    (vals[k-1] > 0) == (f > 0) == (vals[k+1] > 0):
                out += [pts[k], pts[k]]
    return sorted(out)


def span_roots(c, touch=None):
    """The real roots of an mpf polynomial on a span's parameter, [0, 1]
    widened by `SPAN_SLACK` (roots outside [0, 1] miss the span)."""
    return _roots_in(c, -SPAN_SLACK, 1+SPAN_SLACK, touch)


class SplineSpan:
    """S9f.1: one Bezier span of a profile's spline segment, `S(tau) = (X(tau),
    Y(tau))` on [0, 1] with exact (Fraction) power coefficients, traversed
    from `p` to `q` in the stored ring's (counter-clockwise) order; `left`
    as a segment's (the material on the left of the traversal)."""
    kind = 'spline'

    def __init__(self, ctrl, left, outer):
        self.ctrl = tuple((F(x), F(y)) for x, y in ctrl)
        self.p, self.q = self.ctrl[0], self.ctrl[-1]
        self.left, self.outer = left, outer
        self.deg = len(ctrl)-1
        self.X = ptrim(power_of([c[0] for c in self.ctrl]))
        self.Y = ptrim(power_of([c[1] for c in self.ctrl]))
        self.dX, self.dY = pder(self.X), pder(self.Y)
        self.Xm, self.Ym = [M(c) for c in self.X], [M(c) for c in self.Y]
        self.dXm, self.dYm = [M(c) for c in self.dX], [M(c) for c in self.dY]
        self.pm, self.qm = Mv(self.p), Mv(self.q)

    def point(self, tau):
        return peval(self.Xm, tau), peval(self.Ym, tau)

    def tangent(self, tau):
        return peval(self.dXm, tau), peval(self.dYm, tau)

    def outside(self, tau):
        """How far the curve's point at `tau` lies from the span (0 on it)."""
        if 0 <= tau <= 1:
            return mp.mpf(0)
        end = self.pm if tau < 0 else self.qm
        x, y = self.point(tau)
        return mp.sqrt((x-end[0])**2+(y-end[1])**2)

    def line_poly(self, base, d):
        """`(S(tau) - base) x d` (mpf coefficients)."""
        out = [self.Xm[k]*d[1] if k < len(self.Xm) else 0 for k in range(self.deg+1)]
        for k in range(len(self.Ym)):
            out[k] -= self.Ym[k]*d[0]
        out[0] -= base[0]*d[1]-base[1]*d[0]
        return out

    def line_roots(self, base, d, touch=None):
        """[(tau, t)]: the line `base + t d`'s crossings, sorted by `tau`."""
        a = d[0]*d[0]+d[1]*d[1]
        out = []
        for tau in span_roots(self.line_poly(base, d), touch):
            x, y = self.point(tau)
            out.append((tau, ((x-base[0])*d[0]+(y-base[1])*d[1])/a))
        return out

    def exact_line_poly(self, Q0, Q1):
        """`(S(tau) - Q0) x Q1` exactly (Fractions)."""
        return psub(pscale(psub(self.X, [F(Q0[0])]), F(Q1[1])), pscale(psub(self.Y, [F(Q0[1])]), F(Q1[0])))

    def point_events(self, Q0, Q1):
        """[(tau, x)] with `tau` in the span where the point `Q0 + x Q1`
        (exact, moving) lies on the curve: the real roots of the exact
        degree-`p` polynomial `(S(tau) - Q0) x Q1`, `x` read off along
        `Q1`. A fixed point (`Q1 = 0`) gives none."""
        Q0, Q1 = (F(Q0[0]), F(Q0[1])), (F(Q1[0]), F(Q1[1]))
        if Q1 == (0, 0):
            return []
        poly = self.exact_line_poly(Q0, Q1)
        assert not pzero(poly), 'a moving point along a straight span (OutOfDomain)'
        a = M(Q1[0]*Q1[0]+Q1[1]*Q1[1])
        out = []
        for tau in real_roots(poly):
            if -SPAN_SLACK <= tau <= 1+SPAN_SLACK:
                x, y = self.point(tau)
                out.append((tau, ((x-M(Q0[0]))*M(Q1[0])+(y-M(Q0[1]))*M(Q1[1]))/a))
        return out

    def tangent_params(self, d):
        """The `tau` in the span where the curve runs along `d` (exact
        direction): roots of `S'(tau) x d`, degree `p - 1`."""
        poly = psub(pscale(self.dX, F(d[1])), pscale(self.dY, F(d[0])))
        if pzero(poly):
            return []
        return [t for t in real_roots(poly) if -SPAN_SLACK <= t <= 1+SPAN_SLACK]

    def distance(self, u, v):
        """The distance from `(u, v)` to the span."""
        P = [(self.Xm[k] if k < len(self.Xm) else 0) for k in range(self.deg+1)]
        Q = [(self.Ym[k] if k < len(self.Ym) else 0) for k in range(self.deg+1)]
        P[0] -= u
        Q[0] -= v
        dP, dQ = _mp_der(P), _mp_der(Q)
        poly = [mp.mpf(0)]*(2*self.deg)
        for i, x in enumerate(P):
            for j, y in enumerate(dP):
                poly[i+j] += x*y
        for i, x in enumerate(Q):
            for j, y in enumerate(dQ):
                poly[i+j] += x*y
        cand = [mp.mpf(0), mp.mpf(1)]+[t for t in _roots_in(poly, mp.mpf(0), mp.mpf(1), None)]
        best = None
        for t in cand:
            x, y = self.point(t)
            dist = mp.sqrt((x-u)**2+(y-v)**2)
            best = dist if best is None else min(best, dist)
        return best

    def green(self):
        """Exact `(x dy - y dx) / 2`, `x (x dy - y dx) / 3`, `y (x dy - y dx)
        / 3` over the span."""
        cr = psub(pmul(self.X, self.dY), pmul(self.Y, self.dX))
        whole = lambda p: sum(pint(p), F(0))
        return (whole(cr)/2, whole(pmul(self.X, cr))/3, whole(pmul(self.Y, cr))/3)

    def speed(self, tau):
        dx, dy = self.tangent(tau)
        return mp.sqrt(dx*dx+dy*dy)


def arc_green_symmetric(el):
    """An arc's or circle's `(x dy - y dx) / 2`, `x (x dy - y dx) / 3` and `y
    (x dy - y dx) / 3` in closed form (from `a0` to `a1`)."""
    a, b = (mp.mpf(0), 2*mp.pi) if el.kind == 'circle' else (el.a0, el.a1)
    cu, cv, r = el.cm[0], el.cm[1], el.rm
    S1 = mp.sin(b)-mp.sin(a)
    S2 = -(mp.cos(b)-mp.cos(a))
    C2 = (b-a)/2+(mp.sin(2*b)-mp.sin(2*a))/4
    SS = (b-a)/2-(mp.sin(2*b)-mp.sin(2*a))/4
    CS = (mp.sin(b)**2-mp.sin(a)**2)/2
    area = r*(cu*S1+cv*S2+r*(b-a))/2
    mu = r*(cu*cu*S1+cu*cv*S2+cu*r*(b-a)+r*cu*C2+r*cv*CS+r*r*S1)/3
    mw = r*(cv*cu*S1+cv*cv*S2+cv*r*(b-a)+r*cu*CS+r*cv*SS+r*r*S2)/3
    return area, mu, mw


RAY = None


class Profile:
    """A prism's profile: its stored boundaries (every ring counter-clockwise,
    `identity_reference.stored`), the outer one's material on its left, the
    holes' on their right."""

    def __init__(self, boundaries, tolerance):
        self.elements, self.vertices = [], []
        self.splines = False
        for k, b in enumerate(boundaries):
            outer = k == 0
            pts, _ = stored(b, tolerance)
            if pts is None:
                cx, cy, r = b.circle
                self.elements.append(Round((F(cx), F(cy)), F(r), True, outer, outer))
                continue
            if b.segments is None:
                points, segs = pts, [None]*len(pts)
            else:
                points, segs = pts
            ring = [(F(x), F(y)) for x, y in points]
            n = len(ring)
            joins = []
            for i in range(n):
                p, q, s = ring[i], ring[(i+1) % n], segs[i]
                if s is None:
                    self.elements.append(Seg(p, q, outer, outer))
                elif isinstance(s, Spline):
                    # S9f.1: the spline's Bezier spans (blossoms, exact), their
                    # joins vertices of the profile.
                    s.check(points[i], points[(i+1) % n])
                    pieces = s.pieces()
                    assert pieces[0][0] == p and pieces[-1][-1] == q
                    for piece in pieces:
                        self.elements.append(SplineSpan(piece, outer, outer))
                    joins += [piece[0] for piece in pieces[1:]]
                    self.splines = True
                else:
                    cx, cy, r, ccw = s
                    self.elements.append(Round((F(cx), F(cy)), F(r), False, ccw == outer, outer, p, q, ccw))
            self.vertices += ring+joins
        self.vm = [Mv(v) for v in self.vertices]

    # ---- numeric predicates

    def inside(self, u, v):
        """Strictly inside, by the parity of a ray's crossings (a fixed
        direction of irrational slope; the points asked are never on the
        boundary)."""
        global RAY
        if RAY is None or RAY[0] != mp.mp.prec:
            a = mp.mpf('0.6180339887498948482045868343656381177203')
            RAY = (mp.mp.prec, (mp.cos(a), mp.sin(a)))
        d = RAY[1]
        count = 0
        for el in self.elements:
            if el.kind == 'seg':
                den = cross2(d, el.em)
                if den == 0:
                    continue
                w = (el.pm[0]-u, el.pm[1]-v)
                t = cross2(w, el.em)/den
                s = cross2(w, d)/den
                if t > 0 and 0 <= s < 1:
                    count += 1
            elif el.kind == 'spline':
                for tau, t in el.line_roots((u, v), d):
                    if t > 0 and 0 <= tau < 1:
                        count += 1
            else:
                W = (u-el.cm[0], v-el.cm[1])
                b = d[0]*W[0]+d[1]*W[1]
                cc = W[0]**2+W[1]**2-el.rm**2
                disc = b*b-cc
                if disc <= 0:
                    continue
                sq = mp.sqrt(disc)
                for t in (-b-sq, -b+sq):
                    if t > 0:
                        phi = mp.atan2(v+t*d[1]-el.cm[1], u+t*d[0]-el.cm[0])
                        if el.range_bad(phi) == 0:
                            count += 1
        return count % 2 == 1

    def boundary_distance(self, u, v):
        best = None
        for el in self.elements:
            if el.kind == 'seg':
                w = (u-el.pm[0], v-el.pm[1])
                t = (w[0]*el.em[0]+w[1]*el.em[1])/el.len**2
                t = min(max(t, 0), 1)
                dist = mp.sqrt((w[0]-t*el.em[0])**2+(w[1]-t*el.em[1])**2)
            elif el.kind == 'spline':
                dist = el.distance(u, v)
            else:
                W = (u-el.cm[0], v-el.cm[1])
                rho = mp.sqrt(W[0]**2+W[1]**2)
                phi = mp.atan2(W[1], W[0])
                if el.range_bad(phi) == 0:
                    dist = abs(rho-el.rm)
                else:
                    dist = min(mp.sqrt((u-a)**2+(v-b)**2) for a, b in (el.point(el.a0), el.point(el.a1)))
            best = dist if best is None else min(best, dist)
        return best

    def closure_bad(self, u, v):
        """0 inside, else the distance to the boundary."""
        if self.inside(u, v):
            return mp.mpf(0)
        return self.boundary_distance(u, v)

    def crossings(self, base, d, tol=None, skip=()):
        """The line `base + t d`'s crossings of every element: [(t, (element,
        branch), bad)], `bad` how far the crossing misses the element (0 on
        it); a circle's tangency within `tol` counts as a double crossing."""
        out = []
        for i, el in enumerate(self.elements):
            if i in skip:
                continue
            if el.kind == 'seg':
                den = cross2(d, el.em)
                if den == 0:
                    continue
                w = (el.pm[0]-base[0], el.pm[1]-base[1])
                t = cross2(w, el.em)/den
                s = cross2(w, d)/den
                bad = max(0, -s, s-1)*el.len
                out.append((t, (i, 0), bad))
            elif el.kind == 'spline':
                # Tagged by their order along the span among those on it
                # (fixed between breakpoints); a crossing off the span by
                # its distance from it (and a negative tag).
                touch = None if tol is None else tol*mp.sqrt(d[0]**2+d[1]**2)
                k = j = 0
                for tau, t in el.line_roots(base, d, touch):
                    bad = el.outside(tau)
                    if bad == 0:
                        out.append((t, (i, k), bad))
                        k += 1
                    else:
                        j += 1
                        out.append((t, (i, -j), bad))
            else:
                a = d[0]**2+d[1]**2
                W = (base[0]-el.cm[0], base[1]-el.cm[1])
                b = d[0]*W[0]+d[1]*W[1]
                cc = W[0]**2+W[1]**2-el.rm**2
                disc = b*b-a*cc
                if disc < 0:
                    if tol is None or disc < -2*a*el.rm*tol:
                        continue
                    disc = mp.mpf(0)
                sq = mp.sqrt(disc)
                for s in (-1, 1):
                    t = (-b+s*sq)/a
                    phi = mp.atan2(base[1]+t*d[1]-el.cm[1], base[0]+t*d[0]-el.cm[0])
                    out.append((t, (i, s), el.range_bad(phi)))
        out.sort(key=lambda c: c[0])
        return out

    def chords(self, base, d):
        """The line's intervals inside the profile: [((t0, tag0), (t1, tag1))]."""
        hits = [(t, tag) for t, tag, bad in self.crossings(base, d) if bad == 0]
        assert len(hits) % 2 == 0, 'an odd chord (a line through a vertex or tangent)'
        return [(hits[i], hits[i+1]) for i in range(0, len(hits), 2)]

    # ---- exact data

    def moments(self):
        """Area and first moments of the region by Green's theorem (mpf)."""
        A, Mu, Mw = mp.mpf(0), mp.mpf(0), mp.mpf(0)
        if self.splines:
            # S9f.1: the symmetric forms on every element (spans exactly).
            for el in self.elements:
                sign = 1 if el.outer else -1
                if el.kind == 'seg':
                    (x0, y0), (x1, y1) = el.pm, el.qm
                    cr = x0*y1-x1*y0
                    g = (cr/2, (x0+x1)*cr/6, (y0+y1)*cr/6)
                elif el.kind == 'spline':
                    g = tuple(M(x) for x in el.green())
                else:
                    g = arc_green_symmetric(el)
                A, Mu, Mw = A+sign*g[0], Mu+sign*g[1], Mw+sign*g[2]
            return A, Mu, Mw
        for el in self.elements:
            sign = 1 if el.outer else -1
            if el.kind == 'seg':
                (x0, y0), (x1, y1) = el.pm, el.qm
                cr = x0*y1-x1*y0
                A += sign*cr/2
                Mu += sign*(x0+x1)*cr/6
                Mw += sign*(y0+y1)*cr/6
            else:
                a, b = el.a0, el.a1
                if el.kind == 'circle':
                    a, b = mp.mpf(0), 2*mp.pi
                cu, cv, r = el.cm[0], el.cm[1], el.rm
                A += sign*(r*cu*(mp.sin(b)-mp.sin(a))-r*cv*(mp.cos(b)-mp.cos(a))+r*r*(b-a))/2
                ic = lambda t: mp.sin(t)
                ic2 = lambda t: t/2+mp.sin(2*t)/4
                ic3 = lambda t: mp.sin(t)-mp.sin(t)**3/3
                is_ = lambda t: -mp.cos(t)
                is2 = lambda t: t/2-mp.sin(2*t)/4
                is3 = lambda t: -mp.cos(t)+mp.cos(t)**3/3
                Mu += sign*(cu*cu*r*(ic(b)-ic(a))+2*cu*r*r*(ic2(b)-ic2(a))+r**3*(ic3(b)-ic3(a)))/2
                Mw += sign*(cv*cv*r*(is_(b)-is_(a))+2*cv*r*r*(is2(b)-is2(a))+r**3*(is3(b)-is3(a)))/2
        return A, Mu, Mw

    def perimeter(self):
        total = mp.mpf(0)
        for el in self.elements:
            if el.kind == 'spline':
                total += integrate(lambda t, el=el: [el.speed(t)], mp.mpf(0), mp.mpf(1), mp.mpf(10)**-36)[0][0]
                continue
            total += el.len if el.kind == 'seg' else el.rm*el.sweep
        return total

    def bounds(self):
        """A rational box holding the profile."""
        us = [v[0] for v in self.vertices]
        vs = [v[1] for v in self.vertices]
        for el in self.elements:
            if el.kind == 'spline':
                us += [c[0] for c in el.ctrl]
                vs += [c[1] for c in el.ctrl]
            elif el.kind != 'seg':
                us += [el.c[0]-el.r, el.c[0]+el.r]
                vs += [el.c[1]-el.r, el.c[1]+el.r]
        return min(us), max(us), min(vs), max(vs)


# ------------------------------------------------------------------ prisms

class Prism:
    def __init__(self, case):
        self.case = case
        o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(case.frame))
        self.o, self.x, self.y, self.n = o, x, y, n
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        assert self.lo < self.hi
        self.profile = Profile(case.boundaries, case.tolerance)
        self.det = det3(x, y, n)
        assert self.det > 0
        self.inv = tuple(scale(r, 1/self.det) for r in (cross(y, n), cross(n, x), cross(x, y)))
        self.om, self.xm, self.ym, self.nm = Mv(o), Mv(x), Mv(y), Mv(n)
        self.invm = tuple(Mv(r) for r in self.inv)
        self.orthonormal = (dot(x, x) == 1 and dot(y, y) == 1 and dot(n, n) == 1 and dot(x, y) == 0
                            and dot(x, n) == 0 and dot(y, n) == 0)
        self.faces = self.make_faces()

    def local(self, X):
        return apply(self.invm, sub(X, self.om))

    def world(self, u, v, w):
        return tuple(self.om[i]+u*self.xm[i]+v*self.ym[i]+w*self.nm[i] for i in range(3))

    def world_exact(self, u, v, w):
        return tuple(self.o[i]+u*self.x[i]+v*self.y[i]+w*self.n[i] for i in range(3))

    def contains(self, X):
        u, v, w = self.local(X)
        return M(self.lo) < w < M(self.hi) and self.profile.inside(u, v)

    def closure_bad(self, X):
        u, v, w = self.local(X)
        hb = max(0, M(self.lo)-w, w-M(self.hi))
        return hb+self.profile.closure_bad(u, v)

    def measures(self):
        """Closed-form volume, first moments (world) and boundary area."""
        A, Mu, Mw = self.profile.moments()
        h = M(self.hi-self.lo)
        det = M(self.det)
        vol = A*h*det
        c = self.world(Mu/A, Mw/A, M(self.lo+self.hi)/2)
        area = sum((f.closed_area() for f in self.faces), mp.mpf(0))
        return vol, scale(c, vol), area

    def bounds(self):
        """A rational world box holding the prism."""
        u0, u1, v0, v1 = self.profile.bounds()
        pts = [self.world_exact(u, v, w) for u in (u0, u1) for v in (v0, v1) for w in (self.lo, self.hi)]
        return [min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)]

    def make_faces(self):
        faces = [Face(self, 'cap', h=self.lo), Face(self, 'cap', h=self.hi)]
        for i, el in enumerate(self.profile.elements):
            kind = {'seg': 'wall', 'spline': 'spl'}.get(el.kind, 'cyl')
            faces.append(Face(self, kind, element=i))
        return faces


class Face:
    """A face of a prism with its parameterisation: a family of parallel
    lines `X(x, t) = X0(x) + t E` (caps: `x = u`, `t = v`; flat walls: `x`
    along the segment from 0 to 1, `t = w`; cylinder walls: `x = theta`, `t =
    w`)."""

    def __init__(self, prism, kind, h=None, element=None):
        self.prism, self.kind, self.h, self.index = prism, kind, h, element
        p = prism
        if kind == 'cap':
            self.sign = 1 if h == p.hi else -1
            self.normal = scale(cross(p.x, p.y), self.sign)
            self.point = p.world_exact(0, 0, h)
            self.E = p.y
            self.family = 'affine'
            self.X0 = (p.world_exact(0, 0, h), p.x)
            u0, u1, _, _ = p.profile.bounds()
            self.domain = (u0, u1)
            self.metric_const = mp.sqrt(M(dot(cross(p.x, p.y), cross(p.x, p.y))))
        elif kind == 'wall':
            el = p.profile.elements[element]
            self.el = el
            E = add(scale(p.x, el.e[0]), scale(p.y, el.e[1]))
            self.normal = cross(E, p.n) if el.left else cross(p.n, E)
            self.point = p.world_exact(el.p[0], el.p[1], p.lo)
            self.E = p.n
            self.family = 'affine'
            self.X0 = (p.world_exact(el.p[0], el.p[1], 0), E)
            self.domain = (F(0), F(1))
            self.metric_const = mp.sqrt(M(dot(cross(E, p.n), cross(E, p.n))))
        elif kind == 'spl':
            # S9f.1: a spline span's wall, swept by its generatrices.
            el = p.profile.elements[element]
            self.el = el
            self.E = p.n
            self.family = 'poly'
            self.X0 = tuple(add(scale(p.o, 1 if k == 0 else 0),
                                add(scale(p.x, el.X[k] if k < len(el.X) else 0),
                                    scale(p.y, el.Y[k] if k < len(el.Y) else 0)))
                            for k in range(el.deg+1))
            self.domain = (F(0), F(1))
        else:
            el = p.profile.elements[element]
            self.el = el
            self.E = p.n
            self.family = 'trig'
            self.X0 = (p.world_exact(el.c[0], el.c[1], 0), scale(p.x, el.r), scale(p.y, el.r))
            self.domain = (el.start, el.start+el.sweep)
            self.radial_out = el.radial_out
        if self.family == 'trig':
            A0, A1, A2 = self.X0
            self.den = P(1, 0, 1)
            self.num = pvec([add(A0, A1), scale(A2, 2), sub(A0, A1)])
        elif self.family == 'poly':
            self.den = P(1)
            self.num = pvec(list(self.X0))
        else:
            A0, A1 = self.X0
            self.den = P(1)
            self.num = pvec([A0, A1])
        self.X0m = tuple(Mv(v) for v in self.X0)
        self.Em = Mv(self.E)

    def base(self, x):
        if self.family == 'trig':
            A0, A1, A2 = self.X0m
            c, s = mp.cos(x), mp.sin(x)
            return tuple(A0[i]+c*A1[i]+s*A2[i] for i in range(3))
        if self.family == 'poly':
            return tuple(peval([C[i] for C in self.X0m], x) for i in range(3))
        A0, A1 = self.X0m
        return tuple(A0[i]+x*A1[i] for i in range(3))

    def derivative(self, x):
        """`X0'(x)`."""
        if self.family == 'trig':
            _, A1, A2 = self.X0m
            c, s = mp.cos(x), mp.sin(x)
            return tuple(-s*A1[i]+c*A2[i] for i in range(3))
        if self.family == 'poly':
            return tuple(peval(_mp_der([C[i] for C in self.X0m]), x) for i in range(3))
        return self.X0m[1]

    def orient(self):
        """+1 when `X0'(x) x E` points out of the prism, else -1."""
        if self.kind == 'cap':
            return self.sign
        if self.kind == 'cyl':
            return 1 if self.radial_out else -1
        return 1 if self.el.left else -1

    def metric(self, x):
        if self.family == 'trig':
            p = self.prism
            r = self.el.rm
            t = tuple(-mp.sin(x)*p.xm[i]+mp.cos(x)*p.ym[i] for i in range(3))
            c = cross(t, p.nm)
            return r*mp.sqrt(dot(c, c))
        if self.family == 'poly':
            c = cross(self.derivative(x), self.Em)
            return mp.sqrt(dot(c, c))
        return self.metric_const

    def closed_area(self):
        p = self.prism
        if self.kind == 'cap':
            return p.profile.moments()[0]*self.metric_const
        if self.kind == 'wall':
            return self.metric_const*M(p.hi-p.lo)
        if self.kind == 'spl':
            return integrate(lambda x: [self.metric(x)], mp.mpf(0), mp.mpf(1),
                             mp.mpf(10)**-36)[0][0]*M(p.hi-p.lo)
        if p.orthonormal:
            return self.el.rm*self.el.sweep*M(p.hi-p.lo)
        return integrate(lambda x: [self.metric(x)], M(self.domain[0]), M(self.domain[1]),
                         mp.mpf(10)**-36)[0][0]*M(p.hi-p.lo)

    def own_intervals(self, x, tol=None):
        """The line's own extent: [((t0, tag0), (t1, tag1))]."""
        p = self.prism
        if self.kind == 'cap':
            return [((a, ('own',)+ta), (b, ('own',)+tb)) for (a, ta), (b, tb) in
                    p.profile.chords((x, mp.mpf(0)), (mp.mpf(0), mp.mpf(1)))]
        return [((M(p.lo), ('own', 'lo')), (M(p.hi), ('own', 'hi')))]

    def own_descs(self):
        p = self.prism
        if self.kind == 'cap':
            out = []
            for i, el in enumerate(p.profile.elements):
                d = element_desc(el, (P(0, 1), P(0)), P(1), (F(0), F(1)), ('own', i))
                if d is not None:
                    out.append(d)
            return out
        return [Desc(pscale(self.den, p.lo), key=('own', 'lo')), Desc(pscale(self.den, p.hi), key=('own', 'hi'))]

    def own_structure(self):
        """Parameters where the own extent changes (caps: vertices and circle
        extremes)."""
        p = self.prism
        if self.kind != 'cap':
            return []
        out = [M(v[0]) for v in p.profile.vertices]
        for el in p.profile.elements:
            if el.kind == 'spline':
                # Where the span runs along the sweep's lines (`u` extremes).
                out += [el.point(t)[0] for t in el.tangent_params((F(0), F(1)))]
            elif el.kind != 'seg':
                out += [M(el.c[0]-el.r), M(el.c[0]+el.r)]
        return out

    def own_bad(self, x, t):
        """How far the line's point `t` lies outside the face (0 on it)."""
        p = self.prism
        if self.kind == 'cap':
            return p.profile.closure_bad(x, t)
        return max(0, M(p.lo)-t, t-M(p.hi))


# ------------------------------------------------------------------ integration

_NODES = {}


def nodes(degree):
    key = (degree, mp.mp.prec)
    if key not in _NODES:
        _NODES[key] = GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)
    return _NODES[key]


QUAD = {'intervals': 0, 'evaluations': 0, 'splits': 0}


def integrate(f, a, b, tol, depth=0):
    """The integral over [a, b] of the vector function `f` by Gauss-Legendre
    after `s = a + (b - a)(1 - cos t)/2`, doubling the rule (12, 24, 48, 96
    nodes) until two successive estimates agree within `tol`, else halving
    the interval (each half within 3/5 of `tol`: only a singularity just
    outside the interval keeps one half splitting): (estimate, largest
    difference)."""
    QUAD['intervals'] += 1
    if b <= a:
        return None, mp.mpf(0)
    prev = None
    half = (b-a)/2
    for degree in (3, 4, 5, 6):
        est = None
        for xi, wi in nodes(degree):
            th = (xi+1)*mp.pi/2
            s = a+half*(1-mp.cos(th))
            jac = wi*half*mp.sin(th)*mp.pi/2
            v = f(s)
            QUAD['evaluations'] += 1
            est = [jac*c for c in v] if est is None else [e+jac*c for e, c in zip(est, v)]
        if prev is not None:
            diff = max(abs(x-y) for x, y in zip(est, prev))
            if diff <= tol:
                return est, diff
        prev = est
    assert depth < 60, 'quadrature did not converge'
    QUAD['splits'] += 1
    m = (a+b)/2
    e1, d1 = integrate(f, a, m, tol*3/5, depth+1)
    e2, d2 = integrate(f, m, b, tol*3/5, depth+1)
    return [x+y for x, y in zip(e1, e2)], d1+d2


def merge_breaks(points, lo, hi, gap):
    pts = sorted(p for p in points if lo < p < hi)
    out = [lo]
    for p in pts:
        if p-out[-1] > gap:
            out.append(p)
    if hi-out[-1] <= gap and len(out) > 1:
        out.pop()
    out.append(hi)
    return out


# ------------------------------------------------------------------ 2D convex pieces

def clip(poly, a, b, c):
    """The part of a convex polygon where `a x + b y <= c`."""
    if not poly:
        return []
    out = []
    n = len(poly)
    vals = [a*p[0]+b*p[1]-c for p in poly]
    for i in range(n):
        p, q, fp, fq = poly[i], poly[(i+1) % n], vals[i], vals[(i+1) % n]
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp/(fp-fq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return out if len(out) >= 3 else []


def poly_moments(poly):
    """Area and first moments (Green's theorem over the edges)."""
    A = X = Y = 0
    n = len(poly)
    for i in range(n):
        (x0, y0), (x1, y1) = poly[i], poly[(i+1) % n]
        cr = x0*y1-x1*y0
        A += cr
        X += (x0+x1)*cr
        Y += (y0+y1)*cr
    return A/2, X/6, Y/6


def corner(l1, l2):
    (a1, b1, c1), (a2, b2, c2) = l1, l2
    det = a1*b2-a2*b1
    return ((c1*b2-c2*b1)/det, (a1*c2-a2*c1)/det)


def ccw(poly):
    if poly and poly_moments(poly)[0] < 0:
        return list(reversed(poly))
    return poly


# ------------------------------------------------------------------ slicing

class SplineChord:
    """S9f.1: a spline span's chord lines in a section's plane, `chord_ab .
    (lam, mu) = t(s)` with `t(s)` a crossing of the trace `U0 + s U1 + t U2`
    with the span (implicit: a root of degree `p`; branch `k` the `k`-th
    along the span)."""
    implicit = True
    k = F(0)

    def __init__(self, section, index, el, key):
        self.section, self.index, self.el, self.key = section, index, el, key

    def branches(self):
        return range(self.el.deg)

    def value_at(self, s, br):
        base, d = self.section.trace(s)
        on = [t for tau, t in self.el.line_roots(base, d) if 0 <= tau <= 1]
        assert br < len(on), 'a spline chord missing between breakpoints'
        return on[br]

    def events(self, T):
        """[(tau, s)] where the trace's point at `t = T(s)` (exact, at most
        linear in `s`) lies on the span."""
        sec = self.section
        T = ptrim(T)+[F(0)]
        assert pdeg(T) <= 1, 'a chord event of degree above one'
        Q0 = (sec.U0[0]+T[0]*sec.U2[0], sec.U0[1]+T[0]*sec.U2[1])
        Q1 = (sec.U1[0]+T[1]*sec.U2[0], sec.U1[1]+T[1]*sec.U2[1])
        return self.el.point_events(Q0, Q1)


def implicit(d):
    return getattr(d, 'implicit', False)


class Section:
    """A prism's sections by the planes `X = P0 + s d' + lam g + mu h`: its
    uv trace `U0 + s U1 + t U2` (`t` the chord variable: `mu`, or `lam`)
    and its height `w = W0 + s W1 + cl lam + cm mu`."""

    def __init__(self, prism, P0, dp, g, h, tag):
        self.prism, self.tag = prism, tag
        Lg, Lh = apply(prism.inv, g), apply(prism.inv, h)
        L0, L1 = apply(prism.inv, sub(P0, prism.o)), apply(prism.inv, dp)
        if Lg[0] == 0 and Lg[1] == 0:
            self.var, U2 = 'mu', Lh[:2]
        else:
            assert Lh[0] == 0 and Lh[1] == 0, 'a plane direction off the axis with a trace in both'
            self.var, U2 = 'lam', Lg[:2]
        self.U0, self.U1, self.U2 = L0[:2], L1[:2], U2
        self.W0, self.W1, self.cl, self.cm = L0[2], L1[2], Lg[2], Lh[2]
        self.U0m, self.U1m, self.U2m = Mv(self.U0), Mv(self.U1), Mv(self.U2)
        self.W0m, self.W1m, self.clm, self.cmm = M(self.W0), M(self.W1), M(self.cl), M(self.cm)
        self.lom, self.him = M(prism.lo), M(prism.hi)
        self.chord_ab = (F(0), F(1)) if self.var == 'mu' else (F(1), F(0))
        # Exact lines of the plane: `a lam + b mu = c(s)`.
        prof = prism.profile
        Nb = ((self.U0[0], self.U1[0]), (self.U0[1], self.U1[1]))
        Nb = (P(*Nb[0]), P(*Nb[1]))
        self.lines = []
        # S9f.1: parameters `s` where the trace touches a spline span.
        self.structure_s = []
        for i, el in enumerate(prof.elements):
            if el.kind == 'spline':
                self.lines.append((self.chord_ab, SplineChord(self, i, el, (tag, 'chord', i))))
                den = cross2(self.U1, self.U2)
                if den != 0:
                    for tau in el.tangent_params(self.U2):
                        X = el.point(tau)
                        self.structure_s.append(((X[0]-self.U0m[0])*self.U2m[1]-(X[1]-self.U0m[1])*self.U2m[0])
                                                / M(den))
                continue
            d = element_desc(el, Nb, P(1), self.U2, (tag, 'chord', i))
            if d is not None:
                self.lines.append((self.chord_ab, d))
        for name, b in (('lo', prism.lo), ('hi', prism.hi)):
            self.lines.append(((self.cl, self.cm), Desc(P(b-self.W0, -self.W1), key=(tag, 'height', name))))
        self.structure = []
        for v in prof.vertices:
            w = (P(v[0]-self.U0[0], -self.U1[0]), P(v[1]-self.U0[1], -self.U1[1]))
            self.structure.append(psub(pscale(w[0], self.U2[1]), pscale(w[1], self.U2[0])))
        for ab, d in self.lines:
            if d.k:
                self.structure.append(d.Q)

    def trace(self, s):
        return ((self.U0m[0]+s*self.U1m[0], self.U0m[1]+s*self.U1m[1]), self.U2m)

    def chords(self, s):
        base, d = self.trace(s)
        return self.prism.profile.chords(base, d)

    def heights(self, s):
        W = self.W0m+s*self.W1m
        return self.lom-W, self.him-W

    def halfplanes(self, s, t0, t1):
        lo, hi = self.heights(s)
        cl, cm = self.clm, self.cmm
        if self.var == 'mu':
            chord = [(0, -1, -t0), (0, 1, t1)]
        else:
            chord = [(-1, 0, -t0), (1, 0, t1)]
        return chord+[(-cl, -cm, -lo), (cl, cm, hi)]

    def parallelograms(self, s):
        out = []
        lo, hi = self.heights(s)
        cl, cm = self.clm, self.cmm
        for (t0, _), (t1, _) in self.chords(s):
            if self.var == 'mu':
                ls = [(0, 1, t0), (0, 1, t1)]
            else:
                ls = [(1, 0, t0), (1, 0, t1)]
            hs = [(cl, cm, lo), (cl, cm, hi)]
            poly = [corner(ls[0], hs[0]), corner(ls[0], hs[1]), corner(ls[1], hs[1]), corner(ls[1], hs[0])]
            out.append((ccw(poly), self.halfplanes(s, t0, t1)))
        return out

    def contains(self, s, lam, mu):
        t = mu if self.var == 'mu' else lam
        lo, hi = self.heights(s)
        w = self.clm*lam+self.cmm*mu
        if not (lo < w < hi):
            return False
        return any(a < t < b for (a, _), (b, _) in self.chords(s))

    def closure_bad(self, s, lam, mu):
        t = mu if self.var == 'mu' else lam
        lo, hi = self.heights(s)
        w = self.clm*lam+self.cmm*mu
        hb = max(0, lo-w, w-hi)/mp.sqrt(self.clm**2+self.cmm**2)
        base, d = self.trace(s)
        u, v = base[0]+t*d[0], base[1]+t*d[1]
        return hb+self.prism.profile.closure_bad(u, v)

    def instances(self, s, tol):
        """The plane's lines at `s` (lenient): [(key, (a, b, c), bad)]."""
        out = []
        base, d = self.trace(s)
        for t, (i, br), bad in self.prism.profile.crossings(base, d, tol):
            out.append(((self.tag, 'chord', i, br), tuple(M(x) for x in self.chord_ab)+(t,), bad))
        lo, hi = self.heights(s)
        out.append(((self.tag, 'height', 'lo', 0), (self.clm, self.cmm, lo), mp.mpf(0)))
        out.append(((self.tag, 'height', 'hi', 0), (self.clm, self.cmm, hi), mp.mpf(0)))
        return out


def difference(poly, halfplanes):
    """A convex polygon less a convex region (the intersection of half-planes)
    as convex pieces."""
    pieces, rest = [], poly
    for a, b, c in halfplanes:
        outside = clip(rest, -a, -b, -c)
        if outside:
            pieces.append(outside)
        rest = clip(rest, a, b, c)
        if not rest:
            break
    return pieces


class Slicing:
    def __init__(self, A, B, axis=None):
        self.A, self.B = A, B
        c = cross(A.n, B.n)
        if c == (0, 0, 0):
            self.parallel = True
            axes = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))]
            e = axis if axis is not None else min(axes, key=lambda e: abs(dot(A.n, e)))
            d = cross(A.n, e)
            g, h = A.n, cross(A.n, d)
        else:
            self.parallel = False
            d, g, h = c, A.n, B.n
        self.d, self.g, self.h = d, g, h
        self.dp = scale(d, 1/dot(d, d))
        self.P0 = A.o
        self.J = abs(det3(self.dp, g, h))
        self.sA = Section(A, self.P0, self.dp, g, h, 'A')
        self.sB = Section(B, self.P0, self.dp, g, h, 'B')
        self.implicit = any(implicit(d) for _, d in self.sA.lines+self.sB.lines)
        # S9f.2a: spline chords against the other's arc or spline chords
        # coincide where the trace passes a 2D crossing of the profiles.
        self.cross2d = parallel_crossings(A, B) if self.parallel else []
        lo, hi = [], []
        for prism in (A, B):
            l, u = self.s_range(prism)
            lo.append(l)
            hi.append(u)
        self.ranges = list(zip(lo, hi))
        self.size = max(max(abs(x) for x in bnd) for p in (A, B) for bnd in p.bounds())+1
        # (lam, mu) box holding both prisms' points.
        corners = []
        for p in (A, B):
            b0, b1 = p.bounds()
            for i in range(8):
                X = tuple(b1[k] if (i >> k) & 1 else b0[k] for k in range(3))
                corners.append(self.plane_coords(X))
        self.box = (min(c[1] for c in corners), max(c[1] for c in corners),
                    min(c[2] for c in corners), max(c[2] for c in corners))

    def plane_coords(self, X):
        """(s, lam, mu) of an exact point."""
        R = sub(X, self.P0)
        det = det3(self.dp, self.g, self.h)
        s = det3(R, self.g, self.h)/det
        lam = det3(self.dp, R, self.h)/det
        mu = det3(self.dp, self.g, R)/det
        return s, lam, mu

    def s_range(self, prism):
        """A rational range of `s = d . (X - P0)` over the prism."""
        p = prism.profile
        dx, dy = dot(self.d, prism.x), dot(self.d, prism.y)
        base = dot(self.d, sub(prism.o, self.P0))
        vals = [base+v[0]*dx+v[1]*dy for v in p.vertices]
        for el in p.elements:
            if el.kind == 'spline':
                vals += [base+c[0]*dx+c[1]*dy for c in el.ctrl]
            elif el.kind != 'seg':
                c = base+el.c[0]*dx+el.c[1]*dy
                vals += [c-el.r*(abs(dx)+abs(dy)), c+el.r*(abs(dx)+abs(dy))]
        return min(vals), max(vals)

    # ---- the integrand

    def pieces(self, s):
        """Each operation's convex pieces of the slice."""
        pa, pb = self.sA.parallelograms(s), self.sB.parallelograms(s)
        common, cut = [], []
        for poly, _ in pa:
            for _, hb in pb:
                q = poly
                for h in hb:
                    q = clip(q, *h)
                    if not q:
                        break
                if q:
                    common.append(q)
            rest = [poly]
            for _, hb in pb:
                rest = [piece for q in rest for piece in difference(q, hb)]
            cut += rest
        fuse = cut+[poly for poly, _ in pb]
        return {'fuse': fuse, 'cut': cut, 'common': common}

    def integrand(self, s):
        out = []
        pieces = self.pieces(s)
        for op in OPS:
            A = X = Y = 0
            for q in pieces[op]:
                a, x, y = poly_moments(q)
                A, X, Y = A+a, X+x, Y+y
            out += [A, s*A, X, Y]
        return out

    # ---- breakpoints

    def breakpoints(self):
        lo = min(r[0] for r in self.ranges)
        hi = max(r[1] for r in self.ranges)
        cand = [M(r[0]) for r in self.ranges]+[M(r[1]) for r in self.ranges]
        for sec in (self.sA, self.sB):
            for p in sec.structure:
                cand += real_roots(p)
            cand += sec.structure_s
        # S9f.2a: the traces through the profiles' 2D crossings (a spline's
        # chord meeting the other's arc or spline chord).
        sec = self.sA
        den = M(cross2(sec.U1, sec.U2))
        if den != 0:
            for c in self.cross2d:
                X = c.on(self.A)[2]
                cand.append(((X[0]-sec.U0m[0])*sec.U2m[1]-(X[1]-sec.U0m[1])*sec.U2m[0])/den)
        tol = mp.mpf(10)**-20*self.size
        lines = self.sA.lines+self.sB.lines
        found = []
        # Two parallel lines of different sections coinciding.
        for i, (ab1, d1) in enumerate(self.sA.lines):
            for ab2, d2 in self.sB.lines:
                if cross2(ab1, ab2) != 0:
                    continue
                kappa = ab1[0]/ab2[0] if ab2[0] else ab1[1]/ab2[1]
                if implicit(d1) or implicit(d2):
                    # S9f.1: a spline chord `ab . X = t(s)` meets the other's
                    # parallel line `ab' . X = L(s)` where `t(s) = kappa L(s)`
                    # (`ab = kappa ab'`).
                    if implicit(d1) and implicit(d2):
                        # S9f.2a: two spline chords, at the 2D crossings.
                        assert self.parallel, 'spline walls on both prisms with crossing axes'
                        continue
                    if implicit(d1):
                        sc, other, factor = d1, d2, kappa
                    else:
                        sc, other, factor = d2, d1, 1/kappa
                    if other.k != 0:
                        # S9f.2a: a spline chord against an arc's, at the 2D
                        # crossings.
                        assert self.parallel, 'a spline chord against a surd line'
                        continue
                    for tau, s in sc.events(pscale(other.L, factor)):
                        if M(lo)-tol <= s <= M(hi)+tol and self.relevant(s, (d1.key, d2.key), tol):
                            found.append(s)
                    continue
                L = psub(d1.L, pscale(d2.L, kappa))
                poly = sq_poly(L, [(d1.k, d1.Q), (-kappa*d2.k, d2.Q)])
                if pzero(poly):
                    continue
                for s in real_roots(poly):
                    if M(lo)-tol <= s <= M(hi)+tol and self.relevant(s, (d1.key, d2.key), tol):
                        found.append(s)
        # Three pairwise non-parallel lines, not all of one section, concurrent.
        n = len(lines)
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    trip = (lines[i], lines[j], lines[k])
                    tags = {l[1].key[0] for l in trip}
                    if len(tags) < 2:
                        continue
                    abs_ = [l[0] for l in trip]
                    if any(cross2(abs_[p], abs_[q]) == 0 for p, q in ((0, 1), (0, 2), (1, 2))):
                        continue
                    imp = [l for l in trip if implicit(l[1])]
                    if imp:
                        # S9f.1: the other two lines' meeting, linear in `s`,
                        # on the spline chord.
                        assert len(imp) == 1, 'two spline chords in a triple'
                        (ab_s, sc), = imp
                        (p1, e1), (p2, e2) = [l for l in trip if not implicit(l[1])]
                        if e1.k != 0 or e2.k != 0:
                            # S9f.2b: a spline chord, a cylinder's chord and a
                            # height line concurrent (crossing axes): the
                            # meeting of the two walls at that height.
                            assert not self.parallel, 'a spline chord against a surd line'
                            surd, height = (e1, e2) if e1.k != 0 else (e2, e1)
                            assert height.k == 0, 'two surd lines with a spline chord'
                            for s in self.crossing_events(sc, surd.key, height.key):
                                if M(lo)-tol <= s <= M(hi)+tol and \
                                        self.relevant(s, tuple(l[1].key for l in trip), tol):
                                    found.append(s)
                            continue
                        det = cross2(p1, p2)
                        X0 = psub(pscale(e1.L, p2[1]/det), pscale(e2.L, p1[1]/det))
                        X1 = psub(pscale(e2.L, p1[0]/det), pscale(e1.L, p2[0]/det))
                        T = padd(pscale(X0, ab_s[0]), pscale(X1, ab_s[1]))
                        for tau, s in sc.events(T):
                            if M(lo)-tol <= s <= M(hi)+tol and \
                                    self.relevant(s, tuple(l[1].key for l in trip), tol):
                                found.append(s)
                        continue
                    (a1, b1), (a2, b2), (a3, b3) = abs_
                    cof = (a2*b3-a3*b2, -(a1*b3-a3*b1), a1*b2-a2*b1)
                    L = [F(0)]
                    terms = []
                    for cf, (_, d) in zip(cof, trip):
                        L = padd(L, pscale(d.L, cf))
                        if d.k:
                            terms.append((cf*d.k, d.Q))
                    poly = sq_poly(L, terms)
                    if pzero(poly):
                        continue
                    for s in real_roots(poly):
                        if M(lo)-tol <= s <= M(hi)+tol and \
                                self.relevant(s, tuple(l[1].key for l in trip), tol):
                            found.append(s)
        self.events = len(found)
        return merge_breaks(cand+found, M(lo), M(hi), mp.mpf(10)**-30*self.size)

    def prism_of(self, tag):
        return self.A if tag == 'A' else self.B

    def crossing_events(self, sc, surd_key, height_key):
        """S9f.2b: the slice parameters where the spline chord `sc`, the
        cylinder chord `surd_key` and the height line `height_key` are
        concurrent: the walls' meeting at that height (`wall_cylinder_events`),
        `s = d . (X - P0)`."""
        S = sc.section.prism
        C = self.prism_of(surd_key[0])
        circ = C.profile.elements[surd_key[2]]
        H = self.prism_of(height_key[0])
        h = H.lo if height_key[2] == 'lo' else H.hi
        P0, d = Mv(self.P0), Mv(self.d)
        return [sum(((X[i]-P0[i])*d[i] for i in range(3)), mp.mpf(0))
                for _, X in wall_cylinder_events(S, sc.el, C, circ, H, h)]

    def relevant(self, s, keys, tol):
        """Whether the lines `keys` (some branch of each) meet at `s` on both
        sections' closed boundaries."""
        inst = {}
        for sec in (self.sA, self.sB):
            for key, line, bad in sec.instances(s, tol):
                if bad <= tol:
                    inst.setdefault(key[:3], []).append((line, sec))
        options = []
        for key in keys:
            if key not in inst:
                return False
            options.append(inst[key])

        for combo in itertools.product(*options):
            ls = [c[0] for c in combo]
            if len(ls) == 2:
                (a1, b1, c1), (a2, b2, c2) = ls
                kappa = a1/a2 if abs(a2) > abs(b2) else b1/b2
                if abs(c1-kappa*c2) > tol*(1+abs(kappa)):
                    continue
                # The coinciding lines' portions on the two boundaries overlap.
                if self.portions_overlap(s, combo, tol):
                    return True
                continue
            pair = None
            for p, q in ((0, 1), (0, 2), (1, 2)):
                if abs(ls[p][0]*ls[q][1]-ls[q][0]*ls[p][1]) > 0:
                    pair = (p, q)
                    break
            X = corner(ls[pair[0]], ls[pair[1]])
            r = [l for idx, l in enumerate(ls) if idx not in pair][0]
            if abs(r[0]*X[0]+r[1]*X[1]-r[2]) > tol*mp.sqrt(r[0]**2+r[1]**2):
                continue
            if self.sA.closure_bad(s, *X) <= tol and self.sB.closure_bad(s, *X) <= tol:
                return True
        return False

    def portions_overlap(self, s, combo, tol):
        """Two coinciding lines: some point of their common line (a meeting
        with another of the plane's lines, or the middle between two) lies
        in the closures of both sections."""
        (a, b, c), _ = combo[0]
        norm = a*a+b*b
        foot = (a*c/norm, b*c/norm)
        direction = (-b, a)
        ts = []
        for sec in (self.sA, self.sB):
            for key, (a2, b2, c2), bad in sec.instances(s, tol):
                det = a*b2-a2*b
                if abs(det) <= tol*mp.sqrt(norm*(a2*a2+b2*b2)) or bad > tol:
                    continue
                X = corner((a, b, c), (a2, b2, c2))
                ts.append(((X[0]-foot[0])*direction[0]+(X[1]-foot[1])*direction[1])/norm)
        ts.sort()
        ts = ts+[(p+q)/2 for p, q in zip(ts, ts[1:])]
        for t in ts:
            X = (foot[0]+t*direction[0], foot[1]+t*direction[1])
            if self.sA.closure_bad(s, *X) <= tol and self.sB.closure_bad(s, *X) <= tol:
                return True
        return False

    # ---- measures

    def measure(self):
        breaks = self.breakpoints()
        self.breaks = breaks
        tol = mp.mpf(10)**-33*self.size**4
        total = [mp.mpf(0)]*12
        err = mp.mpf(0)
        for a, b in zip(breaks, breaks[1:]):
            est, diff = integrate(self.integrand, a, b, tol)
            if est is None:
                continue
            total = [x+y for x, y in zip(total, est)]
            err = max(err, diff)
        self.quad_error = err
        out = {}
        J = M(self.J)
        P0, dp, g, h = Mv(self.P0), Mv(self.dp), Mv(self.g), Mv(self.h)
        for k, op in enumerate(OPS):
            A, SA, X, Y = total[4*k:4*k+4]
            vol = J*A
            mom = tuple(J*(P0[i]*A+dp[i]*SA+g[i]*X+h[i]*Y) for i in range(3))
            out[op] = (vol, mom)
        return out

    # ---- solids

    def solids(self):
        """Each operation's number of solids (floats: combinatorics only)."""
        breaks = [float(b) for b in self.breaks]
        size = float(self.size)
        box = [float(x) for x in self.box]
        margin = max(box[1]-box[0], box[3]-box[2])+1.0
        bx = [(box[0]-margin, box[2]-margin), (box[1]+margin, box[2]-margin),
              (box[1]+margin, box[3]+margin), (box[0]-margin, box[3]+margin)]
        thr = 1e-11*margin*margin
        canon = self.canonical_lines()
        intervals = []
        for a, b in zip(self.breaks, self.breaks[1:]):
            if b-a <= mp.mpf(10)**-9*self.size:
                intervals.append(None)
                continue
            intervals.append(self.faces_in(a, b, bx, thr, canon))
        count = {}
        for op in OPS:
            parent = {}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            def union(x, y):
                parent[find(x)] = find(y)

            live = []
            for k, iv in enumerate(intervals):
                if iv is None:
                    continue
                faces = [f for f in iv['faces'] if SET[op](f['A'], f['B'])]
                for f in faces:
                    parent[(k, f['id'])] = (k, f['id'])
                for x in range(len(faces)):
                    for y in range(x+1, len(faces)):
                        sx, sy = faces[x]['signs'], faces[y]['signs']
                        if sum(1 for key in sx if sx[key] != sy.get(key)) == 1 and set(sx) == set(sy):
                            union((k, faces[x]['id']), (k, faces[y]['id']))
                live.append((k, faces))
            for (k1, f1), (k2, f2) in zip(live, live[1:]):
                for f in f1:
                    if f['end'] is None:
                        continue
                    for g in f2:
                        if g['start'] is None:
                            continue
                        if overlap_area(f['end'], g['start']) > thr:
                            union((k1, f['id']), (k2, g['id']))
            count[op] = len({find(x) for x in parent})
        return count

    def canonical_lines(self):
        """Keys of lines equal as functions of `s` mapped to one of them
        (per branch)."""
        items = []
        for sec in (self.sA, self.sB):
            for ab, d in sec.lines:
                for br in d.branches():
                    key = d.key+(br,) if d.key[1] == 'chord' else d.key+(0,)
                    items.append((key, ab, d, br))
        canon = {}
        for i, (key, ab, d, br) in enumerate(items):
            canon.setdefault(key, key)
            for key2, ab2, d2, br2 in items[:i]:
                if canon.get(key2) != key2:
                    continue
                if same_line(ab, d, br, ab2, d2, br2):
                    canon[key] = key2
                    break
        self.line_info = {key: (ab, d, br) for key, ab, d, br in items}
        return canon

    def line_at(self, key, s):
        ab, d, br = self.line_info[key]
        if implicit(d):
            return (float(ab[0]), float(ab[1]), float(d.value_at(s, br)))
        return (float(ab[0]), float(ab[1]), float(d.value(s, br, P(1))))

    def line_end(self, key, s, a, b):
        """A line at an end `s` of the interval [a, b]. S9f.1: a spline
        chord's order may change there (a tangency, a knot), so it is the
        crossing at `s` nearest to the chord just inside the interval (a
        double root at a tangency counted)."""
        ab, d, br = self.line_info[key]
        if not implicit(d):
            return self.line_at(key, s)
        near = a+(b-a)*mp.mpf(10)**-15 if s == a else b-(b-a)*mp.mpf(10)**-15
        want = d.value_at(near, br)
        base, dd = d.section.trace(s)
        touch = mp.mpf(10)**-25*self.size*mp.sqrt(dd[0]**2+dd[1]**2)
        roots = [t for _, t in d.el.line_roots(base, dd, touch)]
        return (float(ab[0]), float(ab[1]), float(min(roots, key=lambda t: abs(t-want))))

    def faces_in(self, a, b, bx, thr, canon):
        s = (a+b)/2
        lines = {}
        for sec in (self.sA, self.sB):
            for (lo_, ta), (hi_, tb) in sec.chords(s):
                for t, tag in ((lo_, ta), (hi_, tb)):
                    key = (sec.tag, 'chord', tag[0], tag[1])
                    lines[canon[key]] = None
            lines[canon[(sec.tag, 'height', 'lo', 0)]] = None
            lines[canon[(sec.tag, 'height', 'hi', 0)]] = None
        keys = sorted(lines)
        vals = {key: self.line_at(key, s) for key in keys}
        faces = [(bx, {})]
        for key in keys:
            la, lb, lc = vals[key]
            nxt = []
            for poly, signs in faces:
                below = clip_f(poly, la, lb, lc)
                above = clip_f(poly, -la, -lb, -lc)
                for part, sg in ((below, -1), (above, 1)):
                    if part and abs(area_f(part)) > thr:
                        ns = dict(signs)
                        ns[key] = sg
                        nxt.append((part, ns))
            faces = nxt
        out = []
        for idx, (poly, signs) in enumerate(faces):
            cx = sum(p[0] for p in poly)/len(poly)
            cy = sum(p[1] for p in poly)/len(poly)
            inA = self.sA.contains(s, mp.mpf(cx), mp.mpf(cy))
            inB = self.sB.contains(s, mp.mpf(cx), mp.mpf(cy))
            ends = []
            for t in (a, b):
                q = bx
                for key in keys:
                    la, lb, lc = self.line_end(key, t, a, b)
                    sg = signs[key]
                    q = clip_f(q, la, lb, lc) if sg < 0 else clip_f(q, -la, -lb, -lc)
                    if not q:
                        break
                ends.append(q if q and abs(area_f(q)) > thr else None)
            out.append({'id': idx, 'signs': signs, 'A': inA, 'B': inB, 'start': ends[0], 'end': ends[1]})
        return {'faces': out}


def same_line(ab, d, br, ab2, d2, br2):
    if implicit(d) or implicit(d2):
        return False
    if cross2(ab, ab2) != 0:
        return False
    kappa = ab[0]/ab2[0] if ab2[0] else ab[1]/ab2[1]
    if ptrim(d.L) != ptrim(pscale(d2.L, kappa)):
        return False
    if (d.k == 0) != (d2.k == 0):
        return False
    if d.k == 0:
        return True
    if ptrim(pscale(d.Q, d.k*d.k)) != ptrim(pscale(d2.Q, kappa*kappa*d2.k*d2.k)):
        return False
    return br*(1 if d.k > 0 else -1) == br2*(1 if kappa*d2.k > 0 else -1)


def clip_f(poly, a, b, c):
    out = []
    n = len(poly)
    vals = [a*p[0]+b*p[1]-c for p in poly]
    for i in range(n):
        p, q, fp, fq = poly[i], poly[(i+1) % n], vals[i], vals[(i+1) % n]
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp/(fp-fq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return out if len(out) >= 3 else []


def area_f(poly):
    return sum(p[0]*q[1]-q[0]*p[1] for p, q in zip(poly, poly[1:]+poly[:1]))/2


def overlap_area(p, q):
    """Area of the intersection of two convex float polygons."""
    q = q if area_f(q) > 0 else list(reversed(q))
    out = p
    for (x0, y0), (x1, y1) in zip(q, q[1:]+q[:1]):
        # Inside: left of the edge.
        a, b = (y1-y0), -(x1-x0)
        c = a*x0+b*y0
        out = clip_f(out, a, b, c)
        if not out:
            return 0.0
    return abs(area_f(out))


# ------------------------------------------------------------------ S9f.2a: parallel walls
#
# Two prisms whose axes are exactly parallel: every wall is swept along the
# common axis, so a spline wall meets the other's arc or spline walls in
# generatrices over the 2D crossings of the profiles' curves projected along
# the axis. The second prism's `(u', v')` are exact affine functions of the
# first's `(u, v)` (`parallel_map`); a crossing of a spline span with a
# circle is a root of the exact polynomial `|S(tau) - c|^2 - r^2` (degree
# `2 p`) of the span mapped into the circle's frame (`real_roots`: Yun's
# factors, `mp.polyroots` at 60 digits, a tangency a root of even
# multiplicity); two spline spans (the second mapped exactly into the
# first's frame) cross where subdivision of their Bezier forms in fractions
# keeps their control boxes meeting down to 2^-26, then Newton's method on
# `C(s) = D(u)` at 40 digits, or on the gradient of `|C(s) - D(u)|^2 / 2`
# where the crossing is tangent (S9a.2's reference's method, no
# resultant). The slicing and the face sweeps take every crossing's
# parameter as a breakpoint (a spurious one costs nothing).


def parallel_map(A, B):
    """The exact affine map from B's `(u', v')` to A's `(u, v)` (axes
    exactly parallel): `((a0, a1, a2), (b0, b1, b2))`, `u = a0 + a1 u' + a2
    v'`."""
    assert cross(A.n, B.n) == (0, 0, 0), 'S9f.2a: axes exactly parallel'
    o = apply(A.inv, sub(B.o, A.o))
    x = apply(A.inv, B.x)
    y = apply(A.inv, B.y)
    return ((o[0], x[0], y[0]), (o[1], x[1], y[1]))


def invert_map(m):
    (a0, a1, a2), (b0, b1, b2) = m
    det = a1*b2-a2*b1
    assert det != 0
    i11, i12, i21, i22 = b2/det, -a2/det, -b1/det, a1/det
    return ((-(i11*a0+i12*b0), i11, i12), (-(i21*a0+i22*b0), i21, i22))


def map_point(m, p):
    return (m[0][0]+m[0][1]*p[0]+m[0][2]*p[1], m[1][0]+m[1][1]*p[0]+m[1][2]*p[1])


def _bez_halves(ctrl):
    cur, left, right = list(ctrl), [ctrl[0]], [ctrl[-1]]
    for _ in range(len(ctrl)-1):
        cur = [((x[0]+y[0])/2, (x[1]+y[1])/2) for x, y in zip(cur, cur[1:])]
        left.append(cur[0])
        right.append(cur[-1])
    return left, right[::-1]


def _box(ctrl):
    xs, ys = [c[0] for c in ctrl], [c[1] for c in ctrl]
    return min(xs), max(xs), min(ys), max(ys)


class Curve2:
    """A Bezier span's exact power coefficients and 40-digit jets."""

    def __init__(self, ctrl):
        self.ctrl = tuple((F(x), F(y)) for x, y in ctrl)
        self.X = ptrim(power_of([c[0] for c in self.ctrl]))
        self.Y = ptrim(power_of([c[1] for c in self.ctrl]))
        self.Xm, self.Ym = [M(c) for c in self.X], [M(c) for c in self.Y]
        self.dXm, self.dYm = [M(c) for c in pder(self.X)], [M(c) for c in pder(self.Y)]
        self.ddXm, self.ddYm = [M(c) for c in pder(pder(self.X))], [M(c) for c in pder(pder(self.Y))]

    def point(self, t):
        return peval(self.Xm, t), peval(self.Ym, t)

    def d1(self, t):
        return peval(self.dXm, t), peval(self.dYm, t)

    def d2(self, t):
        return peval(self.ddXm, t), peval(self.ddYm, t)


def span_span_meets(P, Q, reach=F(0), limit=F(1, 2**26)):
    """Candidate parameter boxes `(s0, s1, u0, u1)` where two spans'
    control boxes meet (within `reach`) down to `limit`."""
    stack = [(P.ctrl, F(0), F(1), Q.ctrl, F(0), F(1))]
    out, steps = [], 0
    while stack:
        cp, s0, s1, cq, u0, u1 = stack.pop()
        steps += 1
        assert steps < 2000000, 'spline/spline subdivision did not settle (an overlap?)'
        a, b = _box(cp), _box(cq)
        if a[0] > b[1]+reach or b[0] > a[1]+reach or a[2] > b[3]+reach or b[2] > a[3]+reach:
            continue
        if s1-s0 <= limit and u1-u0 <= limit:
            out.append((s0, s1, u0, u1))
            continue
        if s1-s0 >= u1-u0:
            L, R = _bez_halves(cp)
            m = (s0+s1)/2
            stack += [(L, s0, m, cq, u0, u1), (R, m, s1, cq, u0, u1)]
        else:
            L, R = _bez_halves(cq)
            m = (u0+u1)/2
            stack += [(cp, s0, s1, L, u0, m), (cp, s0, s1, R, m, u1)]
    return out


def _newton_meet(P, Q, s, u):
    """Newton's method on `P(s) = Q(u)`: (s, u), or None where singular."""
    for _ in range(100):
        (x1, y1), (x2, y2) = P.point(s), Q.point(u)
        fx, fy = x1-x2, y1-y2
        (a, b), (c, d) = P.d1(s), Q.d1(u)
        det = b*c-a*d
        if abs(det) <= mp.mpf(10)**-30*(a*a+b*b+c*c+d*d):
            return None
        ds, du = (fx*d-c*fy)/det, (b*fx-a*fy)/det
        s, u = s+ds, u+du
        if abs(ds)+abs(du) < mp.mpf(10)**-45:
            break
    return s, u


def _newton_near(P, Q, s, u):
    """Newton's method on the gradient of `|P(s) - Q(u)|^2 / 2` (a local
    minimum of the distance: a touch or a near miss)."""
    for _ in range(200):
        (x1, y1), (x2, y2) = P.point(s), Q.point(u)
        e = (x1-x2, y1-y2)
        p1, q1 = P.d1(s), Q.d1(u)
        p2, q2 = P.d2(s), Q.d2(u)
        g = (e[0]*p1[0]+e[1]*p1[1], -(e[0]*q1[0]+e[1]*q1[1]))
        h11 = p1[0]**2+p1[1]**2+e[0]*p2[0]+e[1]*p2[1]
        h22 = q1[0]**2+q1[1]**2-(e[0]*q2[0]+e[1]*q2[1])
        h12 = -(p1[0]*q1[0]+p1[1]*q1[1])
        det = h11*h22-h12*h12
        if det == 0:
            return None
        ds, du = -(g[0]*h22-g[1]*h12)/det, -(h11*g[1]-h12*g[0])/det
        s, u = s+ds, u+du
        if abs(ds)+abs(du) < mp.mpf(10)**-42:
            break
    return s, u


def span_crossings(P, Q):
    """[(s, u, sine, distance)] where spans P and Q (one frame) cross or
    touch, `s` and `u` in [0, 1]."""
    out = []
    pad = mp.mpf(2)**-20
    tiny = mp.mpf(10)**-28
    for s0, s1, u0, u1 in span_span_meets(P, Q):
        s, u = M(s0+s1)/2, M(u0+u1)/2
        if any(abs(s-a)+abs(u-b) < mp.mpf(2)**-18 for a, b, _, _ in out):
            continue
        r = _newton_meet(P, Q, s, u)
        if r is None or not (M(s0)-pad <= r[0] <= M(s1)+pad and M(u0)-pad <= r[1] <= M(u1)+pad):
            r = _newton_near(P, Q, s, u)
        if r is None:
            continue
        s, u = r
        if not (-tiny <= s <= 1+tiny and -tiny <= u <= 1+tiny):
            continue
        s, u = min(max(s, M(0)), M(1)), min(max(u, M(0)), M(1))
        (x1, y1), (x2, y2) = P.point(s), Q.point(u)
        dist = mp.sqrt((x1-x2)**2+(y1-y2)**2)
        if dist > mp.mpf(10)**-30:
            continue
        (a, b), (c, d) = P.d1(s), Q.d1(u)
        sine = abs(a*d-b*c)/mp.sqrt((a*a+b*b)*(c*c+d*d))
        if not any(abs(s-a2)+abs(u-b2) < mp.mpf(10)**-25 for a2, b2, _, _ in out):
            out.append((s, u, sine, dist))
    return out


def span_near_misses(P, Q, reach):
    """Local minima of the distance between two spans within `reach` that
    are not crossings: [(s, u, distance)]."""
    out = []
    tiny = mp.mpf(10)**-28
    for s0, s1, u0, u1 in span_span_meets(P, Q, reach, F(1, 2**10)):
        r = _newton_near(P, Q, M(s0+s1)/2, M(u0+u1)/2)
        if r is None:
            continue
        s, u = r
        if not (-tiny <= s <= 1+tiny and -tiny <= u <= 1+tiny):
            continue
        (x1, y1), (x2, y2) = P.point(s), Q.point(u)
        dist = mp.sqrt((x1-x2)**2+(y1-y2)**2)
        if dist <= mp.mpf(10)**-30:
            continue
        if not any(abs(s-a)+abs(u-b) < mp.mpf(10)**-20 for a, b, _ in out):
            out.append((s, u, dist))
    return out


def circle_poly(ctrl, c, r):
    """`|S(tau) - c|^2 - r^2` of a span `S` (exact)."""
    X = psub(ptrim(power_of([p[0] for p in ctrl])), [c[0]])
    Y = psub(ptrim(power_of([p[1] for p in ctrl])), [c[1]])
    return psub(padd(pmul(X, X), pmul(Y, Y)), [r*r])


def span_circle(ctrl, rd, size):
    """[(tau, phi, sine, multiplicity)] where a span (in the circle's frame)
    meets a round element `rd` within its range."""
    poly = ptrim(circle_poly(ctrl, rd.c, rd.r))
    assert not pzero(poly), 'a spline span on a circle'
    C = Curve2(ctrl)
    out = []
    for f, mult in squarefree(poly):
        for tau in real_roots(f):
            if not (-mp.mpf(10)**-30 <= tau <= 1+mp.mpf(10)**-30):
                continue
            tau = min(max(tau, M(0)), M(1))
            x, y = C.point(tau)
            phi = mp.atan2(y-rd.cm[1], x-rd.cm[0])
            if rd.range_bad(phi) > mp.mpf(10)**-25*size:
                continue
            dx, dy = C.d1(tau)
            sine = abs(dx*mp.cos(phi)+dy*mp.sin(phi))/mp.sqrt(dx*dx+dy*dy)
            out.append((tau, phi, sine, mult))
    return out


def span_circle_near(ctrl, rd, size):
    """Near misses of a span and a round element: [(tau, distance)] at the
    critical points of `|S - c|^2` where it does not vanish, within the
    element's range."""
    poly = ptrim(circle_poly(ctrl, rd.c, rd.r))
    C = Curve2(ctrl)
    out = []
    for tau in real_roots(pder(poly)):
        if not (0 <= tau <= 1):
            continue
        x, y = C.point(tau)
        phi = mp.atan2(y-rd.cm[1], x-rd.cm[0])
        if rd.range_bad(phi) > mp.mpf(10)**-25*size:
            continue
        dist = abs(mp.sqrt((x-rd.cm[0])**2+(y-rd.cm[1])**2)-rd.rm)
        if dist > mp.mpf(10)**-30:
            out.append((tau, dist))
    return out


class Crossing:
    """A 2D crossing (or touch) of a curved element of each prism: per
    prism its element's index, its parameter (`tau` on a span, the angle on
    an arc or circle) and the point in its `(u, v)`."""

    def __init__(self, A, B, ia, ta, XA, ib, tb, XB, sine, touch):
        self.at = {id(A): (ia, ta, XA), id(B): (ib, tb, XB)}
        self.sine, self.touch = sine, touch

    def on(self, prism):
        return self.at[id(prism)]


def parallel_crossings(A, B):
    """Every crossing of a spline span of either prism with an arc, circle
    or span of the other (S9f.2a; axes exactly parallel)."""
    if not (A.profile.splines or B.profile.splines):
        return []
    curved = lambda p: [(i, el) for i, el in enumerate(p.profile.elements) if el.kind != 'seg']
    m = parallel_map(A, B)
    inv = invert_map(m)
    size = max([M(1)]+[abs(M(c)) for p in (A, B) for v in p.profile.vertices for c in v])
    out = []
    for ia, ea in curved(A):
        for ib, eb in curved(B):
            if ea.kind != 'spline' and eb.kind != 'spline':
                continue
            if ea.kind == 'spline' and eb.kind == 'spline':
                P = Curve2(ea.ctrl)
                Q = Curve2(tuple(map_point(m, c) for c in eb.ctrl))
                Qb = Curve2(eb.ctrl)
                for s, u, sine, _ in span_crossings(P, Q):
                    out.append(Crossing(A, B, ia, s, P.point(s), ib, u, Qb.point(u), sine,
                                        sine < mp.mpf(10)**-20))
            elif ea.kind == 'spline':
                ctrl = tuple(map_point(inv, c) for c in ea.ctrl)
                P = Curve2(ea.ctrl)
                for tau, phi, sine, mult in span_circle(ctrl, eb, size):
                    XB = (eb.cm[0]+eb.rm*mp.cos(phi), eb.cm[1]+eb.rm*mp.sin(phi))
                    out.append(Crossing(A, B, ia, tau, P.point(tau), ib, phi, XB, sine, mult > 1))
            else:
                ctrl = tuple(map_point(m, c) for c in eb.ctrl)
                Q = Curve2(eb.ctrl)
                for tau, phi, sine, mult in span_circle(ctrl, ea, size):
                    XA = (ea.cm[0]+ea.rm*mp.cos(phi), ea.cm[1]+ea.rm*mp.sin(phi))
                    out.append(Crossing(A, B, ia, phi, XA, ib, tau, Q.point(tau), sine, mult > 1))
    return out


# ------------------------------------------------------------------ S9f.2b: crossing axes
#
# A spline span's wall `X = o + S(tau) + w n` of one prism against a
# cylinder wall of another whose axis crosses it: the cylinder's equation
# in its own frame along the wall's ruling is quadratic in `w` with
# coefficients of degrees `p` and `2 p` in `tau`, so the meeting's points at
# a fixed height of either prism (a cap of the spline prism: `w` constant;
# a cap of the cylinder prism: `w` linear in `S(tau)`) are the real roots of
# an exact polynomial of degree `2 p` in `tau` (`real_roots`: Yun's factors,
# `mp.polyroots` at 60 digits), and the turning points of the meeting over
# `tau` (the ruling tangent to the cylinder) the roots of its discriminant,
# of degree `2 p` too. No resultant.


def _wall_local(S, el, C):
    """The cylinder prism `C`'s local coordinates of the wall's point at
    `(tau, w)`: [(P_k(tau), q_k)] with coordinate `k = P_k + w q_k`."""
    X = [padd(padd(pscale(el.X, S.x[i]), pscale(el.Y, S.y[i])), [S.o[i]-C.o[i]]) for i in range(3)]
    out = []
    for row in C.inv:
        P = [F(0)]
        for i in range(3):
            P = padd(P, pscale(X[i], row[i]))
        out.append((ptrim(P), dot(row, S.n)))
    return out


def wall_cylinder_events(S, el, C, circ, H, h):
    """[(tau, X)]: the points of the spline span `el`'s wall (prism `S`, `tau`
    in the span widened by `SPAN_SLACK`) on the cylinder of the circle or
    arc `circ` of prism `C` (crossing axes) at the height `h` of prism `H`
    (`S` or `C`), `X` in world coordinates (mpf). On a cap of `C` holding
    the wall's axis direction (its ruling along the cap plane) the points
    lie on the generatrices at the roots of `C`'s height (degree `p`), each
    met by the cylinder at the roots of a quadratic, at 40 digits (a tower
    field: S9f.2b.2's)."""
    L = _wall_local(S, el, C)
    (P0, q0), (P1, q1), (P2, q2) = L
    cu, cv, r = circ.c[0], circ.c[1], circ.r
    if H is S:
        U = padd(P0, [q0*h])
        V = padd(P1, [q1*h])
    elif q2 != 0:
        # w = (h - P2) / q2.
        wpoly = pscale(psub([F(h)], P2), 1/q2)
        U = padd(P0, pscale(wpoly, q0))
        V = padd(P1, pscale(wpoly, q1))
    else:
        out = []
        for tau in real_roots(psub(P2, [F(h)])):
            if not -SPAN_SLACK <= tau <= 1+SPAN_SLACK:
                continue
            u0 = peval([M(c) for c in P0], tau)-M(cu)
            v0 = peval([M(c) for c in P1], tau)-M(cv)
            a = M(q0)**2+M(q1)**2
            b = u0*M(q0)+v0*M(q1)
            cc = u0*u0+v0*v0-M(r)**2
            disc = b*b-a*cc
            if disc < 0:
                continue
            for s in (-1, 1):
                w = (-b+s*mp.sqrt(disc))/a
                x, y = el.point(tau)
                out.append((tau, tuple(S.om[i]+x*S.xm[i]+y*S.ym[i]+w*S.nm[i] for i in range(3))))
        return out
    E = psub(padd(pmul(psub(U, [cu]), psub(U, [cu])), pmul(psub(V, [cv]), psub(V, [cv]))), [r*r])
    assert not pzero(E), 'a spline wall on a cylinder (S9f.2b)'
    out = []
    for tau in real_roots(E):
        if not -SPAN_SLACK <= tau <= 1+SPAN_SLACK:
            continue
        x, y = el.point(tau)
        if H is S:
            w = M(h)
        else:
            w = (M(h)-peval([M(c) for c in P2], tau))/M(q2)
        out.append((tau, tuple(S.om[i]+x*S.xm[i]+y*S.ym[i]+w*S.nm[i] for i in range(3))))
    return out


def wall_cylinder_turns(S, el, C, circ):
    """The turning points of a spline span's wall's meeting with a cylinder
    (crossing axes): [(tau, X)] at the real roots in the span of the
    discriminant `B^2 - A C` (degree `2 p`) of the cylinder's equation along
    the ruling, `A w^2 + 2 B w + C`, and the discriminant itself (exact)."""
    (P0, q0), (P1, q1), _ = _wall_local(S, el, C)
    U, V = psub(P0, [circ.c[0]]), psub(P1, [circ.c[1]])
    A = q0*q0+q1*q1
    B = padd(pscale(U, q0), pscale(V, q1))
    Cc = psub(padd(pmul(U, U), pmul(V, V)), [circ.r*circ.r])
    D = psub(pmul(B, B), pscale(Cc, A))
    out = []
    if pzero(D):
        return out, D
    for tau in real_roots(D):
        if not -SPAN_SLACK <= tau <= 1+SPAN_SLACK:
            continue
        x, y = el.point(tau)
        w = -peval([M(c) for c in B], tau)/M(A)
        out.append((tau, tuple(S.om[i]+x*S.xm[i]+y*S.ym[i]+w*S.nm[i] for i in range(3))))
    return out, D


def generatrix_tangencies(S, el, C, circ):
    """The points where a generatrix of the cylinder of `circ` (prism `C`),
    projected along the spline prism `S`'s axis, touches the span `el`: on
    the wall's rulings at the span's points whose tangent runs along the
    projection of `C`'s axis, the cylinder's points there (40 digits)."""
    d = apply(S.inv, C.n)[:2]
    if d == (0, 0):
        return []
    (P0, q0), (P1, q1), _ = _wall_local(S, el, C)
    out = []
    for tau in el.tangent_params(d):
        u0 = peval([M(c) for c in P0], tau)-M(circ.c[0])
        v0 = peval([M(c) for c in P1], tau)-M(circ.c[1])
        a = M(q0)**2+M(q1)**2
        b = u0*M(q0)+v0*M(q1)
        cc = u0*u0+v0*v0-M(circ.r)**2
        disc = b*b-a*cc
        if disc < 0:
            continue
        x, y = el.point(tau)
        for s in (-1, 1):
            w = (-b+s*mp.sqrt(disc))/a
            out.append((tau, tuple(S.om[i]+x*S.xm[i]+y*S.ym[i]+w*S.nm[i] for i in range(3))))
    return out


# ------------------------------------------------------------------ face areas

class FaceSweep:
    """One face of `prism` against the other solid `other`: the areas of its
    pieces inside, outside and on the other's faces (same or opposite
    orientation)."""

    def __init__(self, face, other, size, divergence=False, crossings=()):
        self.face, self.other, self.size = face, other, size
        self.divergence = divergence
        self.crossings = crossings
        self.parallel = cross(face.prism.n, other.n) == (0, 0, 0)
        f, o = face, other
        # The face's lines in the other's coordinates: numerators over den.
        rel = tuple(psub(f.num[i], pscale(f.den, o.o[i])) for i in range(3))
        self.loc = pvapply(o.inv, rel)
        self.dir = apply(o.inv, f.E)
        self.dirm = Mv(self.dir)
        self.g = self.dir[2]
        self.duv = self.dir[:2]
        self.coincident = coincident_faces(f, o)
        self.tol = mp.mpf(10)**-20*size
        self.eps = mp.mpf(10)**-30*size
        self.groups, self.structure = self.descs()

    def numeric_line(self, x):
        X0 = self.face.base(x)
        return self.other.local(X0), X0

    def instances(self, x, tol=None):
        """The line's crossings of the other's boundary: [(t, tag, bad)]
        (`bad` how far it misses the other's closed boundary). Strict (no
        `tol`): the other's element hit exactly, its height and cap within
        `eps` (a line in a cap's or a wall's plane meets the rim), `bad` 0
        when accepted, else positive."""
        o, f = self.other, self.face
        base, X0 = self.numeric_line(x)
        strict = tol is None
        eps = self.eps
        out = []
        g = M(self.g)
        lo, hi = M(o.lo), M(o.hi)

        def height_bad(w):
            hb = max(0, lo-w, w-hi)
            return (0 if hb <= eps else hb) if strict else hb

        if self.g != 0:
            for name, b in (('lo', o.lo), ('hi', o.hi)):
                t = (M(b)-base[2])/g
                u, v = base[0]+t*self.dirm[0], base[1]+t*self.dirm[1]
                bad = o.profile.closure_bad(u, v)
                if strict:
                    bad = 0 if bad <= eps else bad
                out.append((t, ('Oh', name), bad))
        if self.duv != (0, 0):
            skip = self.collinear_elements
            for t, (i, br), bad in o.profile.crossings(base[:2], self.dirm[:2], tol, skip):
                if self.collinear_vertices:
                    u, v = base[0]+t*self.dirm[0], base[1]+t*self.dirm[1]
                    if any(abs(u-o.profile.vm[vi][0])+abs(v-o.profile.vm[vi][1]) <= self.tol
                           for vi in self.collinear_vertices):
                        continue
                out.append((t, ('Oe', i, br), bad+height_bad(base[2]+t*g)))
            for vi in self.collinear_vertices:
                vtx = o.profile.vm[vi]
                d = self.dirm[:2]
                t = ((vtx[0]-base[0])*d[0]+(vtx[1]-base[1])*d[1])/(d[0]**2+d[1]**2)
                out.append((t, ('Ov', vi), height_bad(base[2]+t*g)))
        return out

    # ---- events

    def descs(self):
        """Exact crossing parameters by group."""
        o, f = self.other, self.face
        den = f.den
        out = {'own': f.own_descs(), 'Oh': [], 'Oe': [], 'Ov': []}
        structure = []
        Nb = self.loc[:2]
        Nw = self.loc[2]
        # S9f.1: implicit crossings with spline spans (the other's, or the
        # cap's own), and family parameters where a line touches a span.
        self.implicit_oe, self.implicit_own, self.structure_x = [], [], []
        if f.kind == 'cap':
            self.implicit_own = [(i, el) for i, el in enumerate(f.prism.profile.elements) if el.kind == 'spline']
        splines = [(i, el) for i, el in enumerate(o.profile.elements) if el.kind == 'spline']
        # S9f.2a: a curved wall's generatrices against parallel spline walls
        # change class at the 2D crossings only (`crossing_params`).
        parallel_curved = f.family != 'affine' and self.duv == (0, 0)
        # S9f.2b: a cylinder wall's generatrices against spline walls on a
        # crossing axis meet them at implicit crossings whose events are the
        # walls' meetings (`crossing_breaks`).
        self.crossing_curved = bool(splines) and f.family == 'trig' and not self.parallel
        if splines and not parallel_curved and not self.crossing_curved:
            assert f.family == 'affine' and den == P(1), 'a curved face against a spline wall'
            assert all(pdeg(c) <= 1 for c in Nb), 'a family of degree above one'
            Nb0 = (ptrim(Nb[0])[0], ptrim(Nb[1])[0])
            Nb1 = tuple((ptrim(c)+[F(0)])[1] for c in Nb)
            self.line_uv = (Nb0, Nb1)
        if self.g != 0:
            for name, b in (('lo', o.lo), ('hi', o.hi)):
                out['Oh'].append(Desc(pscale(psub(pscale(den, b), Nw), 1/self.g), key=('Oh', name)))
        else:
            for b in (o.lo, o.hi):
                structure.append(psub(Nw, pscale(den, b)))
        self.collinear_vertices, self.collinear_elements = [], set()
        if self.duv != (0, 0):
            d = self.duv
            for i, el in enumerate(o.profile.elements):
                if el.kind == 'spline':
                    self.implicit_oe.append((i, el))
                    if self.crossing_curved:
                        continue
                    c = cross2(Nb1, d)
                    if c != 0:
                        for tau in el.tangent_params(d):
                            X = el.point(tau)
                            self.structure_x.append(((X[0]-M(Nb0[0]))*M(d[1])-(X[1]-M(Nb0[1]))*M(d[0]))/M(c))
                    continue
                de = element_desc(el, Nb, den, d, ('Oe', i))
                if de is None:
                    # Parallel to a segment: collinear when its vertex polynomial vanishes.
                    continue
                out['Oe'].append(de)
                if de.k:
                    structure.append(de.Q)
            for vi, v in enumerate(o.profile.vertices):
                w = (psub(pscale(den, v[0]), Nb[0]), psub(pscale(den, v[1]), Nb[1]))
                poly = psub(pscale(w[0], d[1]), pscale(w[1], d[0]))
                if pzero(poly):
                    self.collinear_vertices.append(vi)
                    a = d[0]*d[0]+d[1]*d[1]
                    out['Ov'].append(Desc(pscale(padd(pscale(w[0], d[0]), pscale(w[1], d[1])), 1/a),
                                          key=('Ov', vi)))
                else:
                    structure.append(poly)
            for i, el in enumerate(o.profile.elements):
                if el.kind == 'seg' and cross2(d, el.e) == 0:
                    # A segment along the line: collinear when both its ends are.
                    ends = [o.profile.vertices.index(el.p), o.profile.vertices.index(el.q)]
                    if all(e in self.collinear_vertices for e in ends):
                        self.collinear_elements.add(i)
        else:
            for i, el in enumerate(o.profile.elements):
                if el.kind == 'seg':
                    w = (psub(Nb[0], pscale(den, el.p[0])), psub(Nb[1], pscale(den, el.p[1])))
                    structure.append(psub(pscale(w[0], el.e[1]), pscale(w[1], el.e[0])))
                elif el.kind == 'spline':
                    if not parallel_curved:
                        self.structure_x += [x for _, x in el.point_events(Nb0, Nb1)]
                else:
                    W = (psub(Nb[0], pscale(den, el.c[0])), psub(Nb[1], pscale(den, el.c[1])))
                    structure.append(psub(padd(pmul(W[0], W[0]), pmul(W[1], W[1])),
                                          pscale(pmul(den, den), el.r*el.r)))
            for v in o.profile.vertices:
                a = psub(Nb[0], pscale(den, v[0]))
                b = psub(Nb[1], pscale(den, v[1]))
                if pzero(a) and pzero(b):
                    continue
                structure.append(b if pzero(a) else a if pzero(b) else pgcd(a, b))
        return out, [p for p in structure if not pzero(p)]

    def params(self, poly):
        """The family parameters of a polynomial's real roots."""
        roots = real_roots(poly)
        if self.face.family != 'trig':
            return roots
        return [2*mp.atan(t) for t in roots]

    def breakpoints(self):
        f = self.face
        groups, structure = self.groups, self.structure
        lo, hi = M(f.domain[0]), M(f.domain[1])
        cand = list(f.own_structure())
        for poly in structure:
            cand += self.params(poly)
        cand += self.structure_x
        cand += self.crossing_params()
        found = []
        pairs = [('own', 'Oh'), ('own', 'Oe'), ('own', 'Ov'), ('Oh', 'Oe'), ('Oh', 'Ov')]
        tol = self.tol
        # S9f.1: an implicit spline crossing meets another crossing where
        # that crossing's point (linear in the family's parameter) lies on
        # the span.
        for i, el in self.implicit_oe:
            if self.crossing_curved:
                found += self.crossing_breaks(i, el)
                cand += [self.face_param(X) for _, X in generatrix_tangencies(self.other, el, f.prism, f.el)]
                continue
            Nb0, Nb1 = self.line_uv
            d = self.duv
            for g in ('own', 'Oh'):
                for de in groups[g]:
                    if de.k != 0 and self.parallel:
                        # S9f.2a: an arc's crossing, at the 2D crossings.
                        continue
                    if de.k != 0 and g == 'own':
                        # S9f.2b: a cylinder prism's cap circle against a
                        # spline wall on a crossing axis.
                        circ = f.prism.profile.elements[de.key[1]]
                        for _, X in wall_cylinder_events(self.other, el, f.prism, circ, f.prism, f.h):
                            for xx in self.in_domain(self.face_param(X)):
                                if self.relevant(xx, de.key, ('Oe', i), tol):
                                    found.append(xx)
                        continue
                    assert de.k == 0 and pdeg(de.L) <= 1, 'a spline crossing against a surd'
                    L = ptrim(de.L)+[F(0)]
                    Q0 = (Nb0[0]+L[0]*d[0], Nb0[1]+L[0]*d[1])
                    Q1 = (Nb1[0]+L[1]*d[0], Nb1[1]+L[1]*d[1])
                    for _, x in el.point_events(Q0, Q1):
                        for xx in self.in_domain(x):
                            if self.relevant(xx, de.key, ('Oe', i), tol):
                                found.append(xx)
        for i, el in self.implicit_own:
            for g in ('Oh', 'Oe', 'Ov'):
                for de in groups[g]:
                    if de.k != 0 and self.parallel:
                        # S9f.2a: an arc's crossing, at the 2D crossings.
                        continue
                    if de.k != 0 and g == 'Oe':
                        # S9f.2b: the spline prism's cap edge against the
                        # other's cylinder on a crossing axis.
                        circ = self.other.profile.elements[de.key[1]]
                        for _, X in wall_cylinder_events(f.prism, el, self.other, circ, f.prism, f.h):
                            for xx in self.in_domain(self.face_param(X)):
                                if self.relevant(xx, ('own', i), de.key, tol):
                                    found.append(xx)
                        continue
                    assert de.k == 0 and pdeg(de.L) <= 1, 'a spline crossing against a surd'
                    L = ptrim(de.L)+[F(0)]
                    for _, x in el.point_events((F(0), L[0]), (F(1), L[1])):
                        for xx in self.in_domain(x):
                            if self.relevant(xx, ('own', i), de.key, tol):
                                found.append(xx)
        for g1, g2 in pairs:
            for d1 in groups[g1]:
                for d2 in groups[g2]:
                    poly = sq_poly(psub(d1.L, d2.L), [(d1.k, d1.Q), (-d2.k, d2.Q)])
                    if pzero(poly):
                        continue
                    for x in self.params(poly):
                        for xx in self.in_domain(x):
                            if self.relevant(xx, d1.key, d2.key, tol):
                                found.append(xx)
        if f.family == 'trig':
            cand += [xx for xx in self.in_domain(mp.pi)]
            cand = [xx for x in cand for xx in self.in_domain(x)]
        self.events = len(found)
        return merge_breaks(cand+found, lo, hi, mp.mpf(10)**-30*max(1, hi-lo))

    def face_param(self, X):
        """The face's family parameter of a point on it: a cap's `u`, a
        cylinder wall's angle (S9f.2b's events)."""
        f = self.face
        u, v, _ = f.prism.local(X)
        if f.kind == 'cap':
            return u
        assert f.kind == 'cyl', 'a crossing event on a face other than a cap or a cylinder'
        return mp.atan2(v-f.el.cm[1], u-f.el.cm[0])

    def crossing_breaks(self, i, el):
        """S9f.2b: the cylinder wall's generatrices against the other's
        spline span `el` (element `i`) on a crossing axis: their implicit
        crossing meets the face's own heights and the other's caps at the
        walls' meetings at those heights."""
        f, o = self.face, self.other
        found = []
        for H, key in ((f.prism, 'own'), (o, 'Oh')):
            for name, h in (('lo', H.lo), ('hi', H.hi)):
                for _, X in wall_cylinder_events(o, el, f.prism, f.el, H, h):
                    for xx in self.in_domain(self.face_param(X)):
                        if self.relevant(xx, (key, name), ('Oe', i), self.tol):
                            found.append(xx)
        return found

    def crossing_params(self):
        """S9f.2a: the face's parameters at the profiles' 2D crossings (a
        cap's sweep line through one, a curved wall's generatrix on it)."""
        f = self.face
        out = []
        for c in self.crossings:
            i, t, X = c.on(f.prism)
            if f.kind == 'cap':
                out.append(X[0])
            elif f.kind in ('cyl', 'spl') and i == f.index:
                out.append(t)
        return out

    def in_domain(self, x):
        lo, hi = M(self.face.domain[0]), M(self.face.domain[1])
        eps = mp.mpf(10)**-25*max(1, hi-lo)
        if self.face.family != 'trig':
            return [x] if lo-eps <= x <= hi+eps else []
        out = []
        for k in range(-3, 4):
            y = x+2*k*mp.pi
            if lo-eps <= y <= hi+eps:
                out.append(y)
        return out

    def relevant(self, x, k1, k2, tol):
        inst = self.instances(x, tol)
        own = []
        if self.face.kind == 'cap':
            for t, tag, bad in self.face.prism.profile.crossings((x, mp.mpf(0)), (mp.mpf(0), mp.mpf(1)), tol):
                own.append((t, ('own',)+tag, bad))
        else:
            p = self.face.prism
            own = [(M(p.lo), ('own', 'lo'), 0), (M(p.hi), ('own', 'hi'), 0)]
        allinst = own+inst

        def match(key, tag):
            return tag[:len(key)] == key

        c1 = [(t, bad) for t, tag, bad in allinst if match(k1, tag) and bad <= tol]
        c2 = [(t, bad) for t, tag, bad in allinst if match(k2, tag) and bad <= tol]
        for t1, _ in c1:
            for t2, _ in c2:
                if abs(t1-t2) <= tol and self.face.own_bad(x, t1) <= tol:
                    return True
        return False

    # ---- classification

    def classify(self, X):
        for G, same in self.coincident:
            if G.region_contains(X):
                return 'same' if same else 'opp'
        return 'in' if self.other.contains(X) else 'out'

    def structure_at(self, x):
        """[(class, tag0, tag1)] of the line's pieces."""
        f = self.face
        inst = [(t, tag) for t, tag, bad in self.instances(x) if bad == 0]
        X0 = f.base(x)
        out = []
        for (a, ta), (b, tb) in f.own_intervals(x):
            inner = sorted((t, tag) for t, tag in inst if a < t < b)
            pts = [(a, ta)]+inner+[(b, tb)]
            for (t0, g0), (t1, g1) in zip(pts, pts[1:]):
                tm = (t0+t1)/2
                X = tuple(X0[i]+tm*f.Em[i] for i in range(3))
                cls = self.classify(X)
                if out and out[-1][0] == cls and out[-1][2] == g0:
                    out[-1] = (cls, out[-1][1], g1)
                else:
                    out.append((cls, g0, g1))
        return out

    def values_at(self, x):
        f = self.face
        vals = {}
        for t, tag, bad in self.instances(x):
            if bad == 0:
                vals[tag] = t
        for (a, ta), (b, tb) in f.own_intervals(x):
            vals[ta], vals[tb] = a, b
        return vals

    def areas(self):
        """{class: area} over the face (with `divergence`, also
        `self.flux`: {class: [X . N / 3, x_i^2 N_i / 2 for i]} integrated,
        `N` the prism's outward normal)."""
        f = self.face
        breaks = self.breakpoints()
        self.breaks = breaks
        classes = ('in', 'out', 'same', 'opp')
        total = {c: mp.mpf(0) for c in classes}
        width = 5 if self.divergence else 1
        flux = {c: [mp.mpf(0)]*4 for c in classes}
        self.quad_error = mp.mpf(0)
        tol = mp.mpf(10)**-33*self.size**2
        E = f.Em
        sense = f.orient() if self.divergence else 1
        # The flux terms scaled to the areas' size (the quadrature's tolerance
        # is one for every component).
        s1, s2 = 1/M(self.size), 1/M(self.size)**2
        for a, b in zip(breaks, breaks[1:]):
            mid = (a+b)/2
            struct = self.structure_at(mid)
            if not struct:
                continue

            def fn(x, struct=struct):
                vals = self.values_at(x)
                out = [mp.mpf(0)]*(4*width)
                m = f.metric(x)
                if self.divergence:
                    X0 = f.base(x)
                    Nv = scale(cross(f.derivative(x), E), sense)
                    XN = dot(X0, Nv)
                for cls, t0, t1 in struct:
                    if t0 not in vals or t1 not in vals:
                        raise AssertionError(f'a missed breakpoint on {describe(f)} at {mp.nstr(x, 15)} '
                                             f'({mp.nstr(a, 15)}, {mp.nstr(b, 15)}): {t0} {t1}')
                    k = classes.index(cls)*width
                    a_, b_ = vals[t0], vals[t1]
                    out[k] += (b_-a_)*m
                    if self.divergence:
                        # E . N = 0: X . N is constant along the line.
                        out[k+1] += XN*(b_-a_)/3*s1
                        for i in range(3):
                            out[k+2+i] += Nv[i]*(X0[i]**2*(b_-a_)+X0[i]*E[i]*(b_**2-a_**2)
                                                 + E[i]**2*(b_**3-a_**3)/3)/2*s2
                return out

            est, diff = integrate(fn, a, b, tol)
            self.quad_error = max(self.quad_error, diff)
            for k, c in enumerate(classes):
                total[c] += est[k*width]
                if self.divergence:
                    scaled = est[k*width+1:k*width+5]
                    scaled = [scaled[0]/s1]+[v/s2 for v in scaled[1:]]
                    flux[c] = [x+y for x, y in zip(flux[c], scaled)]
        self.flux = flux
        return total


def describe(face):
    return f'{face.kind} {face.h if face.kind == "cap" else face.index}'


# ------------------------------------------------------------------ coincident faces

class Region:
    """A face of the other solid lying on the swept face's surface."""

    def __init__(self, prism, face):
        self.prism, self.face = prism, face

    def region_contains(self, X):
        p, f = self.prism, self.face
        u, v, w = p.local(X)
        if f.kind == 'cap':
            return p.profile.inside(u, v)
        if not (M(p.lo) < w < M(p.hi)):
            return False
        el = f.el
        if f.kind == 'wall':
            s = ((u-el.pm[0])*el.em[0]+(v-el.pm[1])*el.em[1])/el.len**2
            return 0 < s < 1
        return el.range_bad(mp.atan2(v-el.cm[1], u-el.cm[0])) == 0


def planar(face):
    return face.kind in ('cap', 'wall')


def coincident_faces(f, other):
    """The other solid's faces on the face's surface (exactly), with
    whether their orientations agree."""
    out = []
    p = f.prism
    for G in other.faces:
        if planar(f) and planar(G):
            if cross(f.normal, G.normal) == (0, 0, 0) and dot(f.normal, sub(G.point, f.point)) == 0:
                out.append((Region(other, G), dot(f.normal, G.normal) > 0))
        elif f.kind == 'cyl' and G.kind == 'cyl':
            if cross(p.n, other.n) != (0, 0, 0) or f.el.r != G.el.r:
                continue
            a = p.world_exact(f.el.c[0], f.el.c[1], 0)
            b = other.world_exact(G.el.c[0], G.el.c[1], 0)
            if cross(sub(b, a), p.n) != (0, 0, 0):
                continue
            same_axes = (p.x, p.y, p.n) == (other.x, other.y, other.n)
            assert same_axes or (p.orthonormal and other.orthonormal), \
                'coincident cylinders in frames whose stored axes differ (not exactly circular)'
            out.append((Region(other, G), f.el.radial_out == G.el.radial_out))
    return out


# ------------------------------------------------------------------ the pair

class Pair:
    def __init__(self, obj, tool, axis=None, divergence=None):
        self.A, self.B = Prism(obj), Prism(tool)
        self.slicing = Slicing(self.A, self.B, axis)
        self.size = self.slicing.size
        if divergence is None:
            divergence = self.A.profile.splines or self.B.profile.splines
        self.divergence = divergence
        self._volumes = None
        self._areas = None
        self._solids = None

    def volumes(self):
        if self._volumes is None:
            self._volumes = self.slicing.measure()
        return self._volumes

    def face_areas(self):
        """[(prism tag, face, {class: area})]."""
        if self._areas is None:
            out = []
            for tag, p, o in (('A', self.A, self.B), ('B', self.B, self.A)):
                for f in p.faces:
                    sweep = FaceSweep(f, o, self.size, self.divergence, self.slicing.cross2d)
                    out.append((tag, f, sweep.areas(), sweep))
            self._areas = out
        return self._areas

    def area(self, op):
        keep = {'fuse': ({'out', 'same'}, {'out'}), 'cut': ({'out', 'opp'}, {'in'}),
                'common': ({'in', 'same'}, {'in'})}[op]
        total = mp.mpf(0)
        for tag, f, cls, _ in self.face_areas():
            for c in keep[0 if tag == 'A' else 1]:
                total += cls[c]
        return total

    def solids(self):
        if self._solids is None:
            self.volumes()
            self._solids = self.slicing.solids()
        return self._solids

    def divergence_volumes(self):
        """S9f.1's second way: each operation's volume and first moments by
        the divergence theorem over its kept face pieces (the tool's pieces
        inside the object bound a cut with their normals reversed)."""
        assert self.divergence
        keep = {'fuse': ({'out', 'same'}, {'out'}, 1), 'cut': ({'out', 'opp'}, {'in'}, -1),
                'common': ({'in', 'same'}, {'in'}, 1)}
        out = {}
        for op in OPS:
            total = [mp.mpf(0)]*4
            for tag, f, cls, sweep in self.face_areas():
                classes = keep[op][0 if tag == 'A' else 1]
                sign = 1 if tag == 'A' else keep[op][2]
                for c in classes:
                    total = [x+sign*y for x, y in zip(total, sweep.flux[c])]
            out[op] = (total[0], tuple(total[1:]))
        return out

    def result(self, op):
        vol, mom = self.volumes()[op]
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        return self.solids()[op], vol, self.area(op), tuple(m/vol for m in mom)


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
