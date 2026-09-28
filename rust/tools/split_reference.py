#!/usr/bin/env python3
"""Independent reference for S8a of REVIEW_NOTES.md: prisms of line and arc
profiles split by a plane.

A prism is its profile region `Omega` (the stored boundaries of
`identity_reference.stored`: an outer path and holes, lines and circular
arcs) in the stored frame axes `(o, x, y, n)` (`stored_axes`: the kernel's `Frame3::new` step by step), swept along `n`
between the offsets. The plane is `m . (p - q) = 0`; in frame coordinates
`m . (p - q) = a u + b v + c w + d` with exact rationals from the stored
binary64 data (the frame's axes taken as stored, not as exactly orthonormal:
`a = m . x`, and so on).

Each piece is the prism's points on one side: its volume and first moments
are `int_Omega` of the length (and `w` moment) of the clipped segment
`{w in [w0, w1] : side(a u + b v + c w + d)}` of the axis over the profile,
computed by slicing `Omega` at each `u` (the boundaries' crossings exact, the
integrand in `v` a polynomial of degree two between breaks, so Simpson's rule
is exact there) and integrating over `u` with `mp.quad` between every break
(vertices, arcs' extremes, the clip lines' crossings of the boundary). Its
area is its caps' parts on its side, its walls' clipped strips (a segment's
in closed form, an arc's by the antiderivative of `A + B cos t + C sin t`),
and the cut face: the region of `Omega` where the plane lies strictly
between the caps, times `|(a, b, c)| / |c|`, or the plane's vertical chord
of `Omega` times the height when `c = 0`.

S8b: a path segment may be a nonrational B-spline (`identity_reference.
Spline`). Its element is `('B', p, q, pieces)`: its Bezier pieces from exact
knot insertion (Boehm's, in Fractions, each interior knot raised to the
degree), each piece's coordinates polynomials in `t` on `[0, 1]` with exact
coefficients. The vertical line `u = x` meets a piece where `x(t) = x` on one
of its `x`-monotone runs (between the roots of `x'(t)`), found by a bracketed
Newton iteration; the breaks add each piece's ends and `x`-extremes and the
points where a clip line's `alpha x(t) + beta y(t) - gamma` vanishes (roots
of a polynomial by the same bracketing on its own monotone runs,
recursively). A spline wall's clipped strip is `mp.quad` of the clipped
height times `|C'(t)|` between those roots. `green_moments` gives a profile
of lines and splines its area and first moments by Green's theorem in exact
Fractions (Bernstein product integrals over Bezier pieces from blossoms, not
from knot insertion), a check of the slicing.

S8e: a face or wire body (a case with `make`, S6) lies on its frame's plane
`w = 0`, and the plane's trace there is `a u + b v + d = 0` (`Planar`). A
sheet's side is a prism of height one cut parallel to its axis (the same
slicing: area and first moments), its perimeter the boundary strictly on
the side plus the trace where it bounds the side; a wire's side is its
boundary's pieces on it (the elements cut exactly where the trace crosses or
touches them, a piece along the trace taking the side of the piece before it
in stored order), lengths and moments in closed form for lines and arcs and
by quadrature for splines. Rows: `side S area perimeter cx cy cz` for a
sheet, `side S length 0 cx cy cz` for a wire (a curve's boundary, its ends,
has no length).

A side holding no solid (the plane missing, touching a vertex, an edge or a
ruling, or lying in a cap: every stored point certainly on the other side or
on the plane) gives no piece: one `whole` row. A side may hold several pieces
(a U's two prongs): its row gives their totals (volume, area) and their
common centre. Rows: `side below|above|whole volume area cx cy cz` (`below`
the side against the plane's normal).
"""
from fractions import Fraction as F
import math

import mpmath as mp

from identity_reference import Spline, arc_sweep, bernstein_product, stored
from curve_surface_reference import stored_axes

mp.mp.dps = 40


def _Q(v):
    return F(v)


def M(v):
    """An mpf of a Fraction, a float or an mpf."""
    if isinstance(v, F):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def _horner(c, t):
    out = M(0)
    for k in reversed(c):
        out = out*t+k
    return out


def _derivative(c):
    return [k*c[k] for k in range(1, len(c))]


def _trim(c):
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return c


def _bracket(f, df, a, b, fa):
    """The root of `f` in `[a, b]`, where it changes sign (`fa = f(a)`):
    Newton steps kept inside the shrinking bracket, bisection otherwise."""
    t = (a+b)/2
    for _ in range(400):
        ft = f(t)
        if ft == 0:
            return t
        if (ft < 0) == (fa < 0):
            a, fa = t, ft
        else:
            b = t
        d = df(t)
        nt = t-ft/d if d != 0 else None
        if nt is None or not (a < nt < b):
            nt = (a+b)/2
        if abs(nt-t) <= M(10)**-39 or b-a <= M(10)**-39:
            return nt
        t = nt
    raise ArithmeticError('bracketed root did not converge')


def _runs(c):
    """0, the roots of `c'` in (0, 1) and 1: the ends of the polynomial
    `c`'s monotone runs on [0, 1]."""
    return [M(0)]+_roots(_derivative(_trim(c)), touching=True)+[M(1)]


def _roots(c, touching=False):
    """The roots of the polynomial `c` (ascending mpf coefficients) in the
    open (0, 1): where it changes sign on one of its monotone runs
    (recursively through its derivatives), and with `touching` also the
    runs' interior ends where it vanishes (to 1e-35 of its coefficients)."""
    c = _trim(c)
    if len(c) <= 1:
        return []
    if len(c) == 2:
        t = -c[0]/c[1]
        return [t] if 0 < t < 1 else []
    ends = _runs(c)
    dc = _derivative(c)
    f, df = (lambda t: _horner(c, t)), (lambda t: _horner(dc, t))
    out = []
    for a, b in zip(ends, ends[1:]):
        fa, fb = f(a), f(b)
        if fa != 0 and fb != 0 and (fa < 0) != (fb < 0):
            out.append(_bracket(f, df, a, b, fa))
    if touching:
        scale = max(abs(k) for k in c)
        out += [t for t in ends[1:-1] if abs(f(t)) <= scale*M(10)**-35]
    return sorted(out)


def _power(ctrl):
    """Power-basis coefficients (ascending, exact) of a Bezier coordinate
    with Fraction controls."""
    from math import comb
    n = len(ctrl)-1
    return [comb(n, k)*sum((-1)**(k-i)*comb(k, i)*ctrl[i] for i in range(k+1)) for k in range(n+1)]


class Bezier:
    """One Bezier piece of a spline element on `t` in [0, 1]: its exact
    controls, 40-digit power coefficients, its `x`-monotone runs, and its
    crossings of vertical lines cached by abscissa."""

    def __init__(self, ctrl):
        self.ctrl = tuple((F(x), F(y)) for x, y in ctrl)
        self.cx = [M(k) for k in _power([p[0] for p in self.ctrl])]
        self.cy = [M(k) for k in _power([p[1] for p in self.ctrl])]
        self.dx, self.dy = _derivative(self.cx), _derivative(self.cy)
        ts = _runs(self.cx)
        self.runs = [(a, b, self.point(a)[0], self.point(b)[0]) for a, b in zip(ts, ts[1:])]
        self.cache = {}

    def reversed(self):
        return Bezier(tuple(reversed(self.ctrl)))

    def point(self, t):
        """The point at `t`, its ends exactly the controls (so neighbouring
        pieces agree there)."""
        if t == 0:
            return M(self.ctrl[0][0]), M(self.ctrl[0][1])
        if t == 1:
            return M(self.ctrl[-1][0]), M(self.ctrl[-1][1])
        return _horner(self.cx, t), _horner(self.cy, t)

    def speed(self, t):
        return mp.sqrt(_horner(self.dx, t)**2+_horner(self.dy, t)**2)

    def ys_at(self, x):
        """The `v` where the piece crosses `u = x` (never at an end or an
        `x`-extreme: those are breaks)."""
        if x not in self.cache:
            ys = []
            for a, b, xa, xb in self.runs:
                if (xa < x < xb) or (xb < x < xa):
                    t = _bracket(lambda t: _horner(self.cx, t)-x, lambda t: _horner(self.dx, t),
                                 a, b, xa-x)
                    ys.append(_horner(self.cy, t))
            self.cache[x] = ys
        return self.cache[x]

    def line(self, alpha, beta, gamma):
        """`alpha x(t) + beta y(t) - gamma` as ascending coefficients."""
        c = [alpha*i+beta*j for i, j in zip(self.cx, self.cy)]
        c[0] -= gamma
        return c

    def line_roots(self, alpha, beta, gamma):
        """Every `t` in (0, 1) where the line `alpha x + beta y = gamma`
        meets the piece, touching ones included."""
        return _roots(self.line(M(alpha), M(beta), M(gamma)), touching=True)

    def line_sign_changes(self, alpha, beta, gamma):
        """The points where `alpha x + beta y - gamma` passes between `< 0`
        and `>= 0` along the piece (the convention of a segment's
        crossings; a touch crosses nothing), one per monotone run."""
        alpha, beta, gamma = M(alpha), M(beta), M(gamma)
        c = self.line(alpha, beta, gamma)
        dc = _derivative(_trim(c))
        g = lambda t: alpha*self.point(t)[0]+beta*self.point(t)[1]-gamma
        ends = _runs(c)
        out = []
        for a, b in zip(ends, ends[1:]):
            fa, fb = g(a), g(b)
            if (fa < 0) != (fb < 0):
                if fa == 0 or fb == 0:
                    out.append(a if fa == 0 else b)
                else:
                    out.append(_bracket(lambda t: _horner(c, t), lambda t: _horner(dc, t), a, b, fa))
        return [self.point(t) for t in out]


def bezier_pieces(spline):
    """The spline's Bezier control points per knot span, by Boehm's knot
    insertion in Fractions: each interior knot inserted until its
    multiplicity is the degree."""
    d = spline.degree
    U = spline.flat_knots()
    P = [(F(x), F(y)) for x, y in spline.poles]
    for k, m in zip(spline.knots[1:-1], spline.mults[1:-1]):
        t = F(k)
        for _ in range(d-m):
            s = max(i for i in range(len(U)-1) if U[i] <= t < U[i+1])
            Q = P[:s-d+1]
            for i in range(s-d+1, s+1):
                a = (t-U[i])/(U[i+d]-U[i])
                Q.append(tuple((1-a)*P[i-1][c]+a*P[i][c] for c in range(2)))
            Q += P[s:]
            P, U = Q, U[:s+1]+[t]+U[s+1:]
    count = len(spline.knots)-1
    assert len(P) == count*d+1
    return [tuple(P[i*d:i*d+d+1]) for i in range(count)]


def elements(boundaries, tolerance):
    """The stored boundaries' elements: ('L', p, q), ('A', p, q, (cx, cy,
    r), sweep) or ('B', p, q, Bezier pieces) (S8b), counter-clockwise outer,
    clockwise holes (as stored: holes reversed so the region is on the
    left)."""
    out = []
    for k, b in enumerate(boundaries):
        pts, _ = stored(b, tolerance)
        items = []
        if pts is None:
            cx, cy, r = b.circle
            p = (cx+r, cy)
            items.append(('A', p, p, (cx, cy, r), 2*mp.pi))
        elif b.segments is not None:
            points, segs = pts
            n = len(points)
            for i in range(n):
                p, q = points[i], points[(i+1) % n]
                if segs[i] is None:
                    items.append(('L', p, q))
                elif isinstance(segs[i], Spline):
                    segs[i].check(p, q)
                    items.append(('B', p, q, tuple(Bezier(c) for c in bezier_pieces(segs[i]))))
                else:
                    items.append(('A', p, q, segs[i][:3], arc_sweep(p, q, segs[i])))
        else:
            n = len(pts)
            for i in range(n):
                items.append(('L', pts[i], pts[(i+1) % n]))
        if k > 0:
            # A hole: the region lies outside it; reverse its direction.
            rev = []
            for it in reversed(items):
                if it[0] == 'L':
                    rev.append(('L', it[2], it[1]))
                elif it[0] == 'B':
                    rev.append(('B', it[2], it[1], tuple(b.reversed() for b in reversed(it[3]))))
                else:
                    rev.append(('A', it[2], it[1], it[3], -it[4]))
            items = rev
        out += items
    return out


def arc_angles(e):
    _, p, _, (cx, cy, r), sweep = e
    a0 = mp.atan2(M(p[1])-cy, M(p[0])-cx)
    return a0, sweep


def crossings(els, x):
    """The `v` values where the boundary crosses the vertical line `u = x`
    (`x` not at a vertex or an arc's extreme)."""
    ys = []
    for e in els:
        if e[0] == 'L':
            (x0, y0), (x1, y1) = e[1], e[2]
            if (x0 < x < x1) or (x1 < x < x0):
                t = (x-M(x0))/(M(x1)-x0)
                ys.append(M(y0)+t*(M(y1)-y0))
        elif e[0] == 'B':
            for piece in e[3]:
                ys += piece.ys_at(x)
        else:
            _, _, _, (cx, cy, r), sweep = e
            dx = x-M(cx)
            if abs(dx) >= r:
                continue
            h = mp.sqrt(M(r)**2-dx*dx)
            a0, sw = arc_angles(e)
            for y in (cy+h, cy-h):
                t = mp.atan2(y-cy, dx)
                rel = (t-a0) % (2*mp.pi) if sw > 0 else (a0-t) % (2*mp.pi)
                if 0 < rel < abs(sw) or abs(abs(sw)-2*mp.pi) < M(10)**-30:
                    ys.append(y)
    return sorted(ys)


def x_breaks(els, lines):
    """Vertices, arcs' leftmost and rightmost points, and where the clip
    lines `alpha u + beta v = gamma` cross the boundary or are vertical."""
    xs = set()
    for e in els:
        xs.add(M(e[1][0]))
        xs.add(M(e[2][0]))
        if e[0] == 'B':
            for piece in e[3]:
                for _, _, xa, xb in piece.runs:
                    xs.update((xa, xb))
        if e[0] == 'A':
            _, _, _, (cx, cy, r), sweep = e
            a0, sw = arc_angles(e)
            for t in (0, mp.pi):
                rel = (t-a0) % (2*mp.pi) if sw > 0 else (a0-t) % (2*mp.pi)
                if rel <= abs(sw):
                    xs.add(M(cx)+r*mp.cos(t))
    for alpha, beta, gamma in lines:
        if beta == 0:
            if alpha != 0:
                xs.add(M(gamma)/alpha)
            continue
        for e in els:
            if e[0] == 'L':
                (x0, y0), (x1, y1) = e[1], e[2]
                # alpha x + beta y = gamma along p + t (q - p)
                f0 = alpha*M(x0)+beta*M(y0)-gamma
                f1 = alpha*M(x1)+beta*M(y1)-gamma
                if f0 != f1:
                    t = f0/(f0-f1)
                    if 0 <= t <= 1:
                        xs.add(M(x0)+t*(M(x1)-x0))
            elif e[0] == 'B':
                for piece in e[3]:
                    for t in piece.line_roots(alpha, beta, gamma):
                        xs.add(piece.point(t)[0])
            else:
                _, _, _, (cx, cy, r), sweep = e
                # alpha (cx + r cos) + beta (cy + r sin) = gamma
                R = mp.sqrt(alpha**2+beta**2)*r
                k = gamma-alpha*M(cx)-beta*M(cy)
                if abs(k) <= R:
                    phi = mp.atan2(beta, alpha)
                    for t in (phi+mp.acos(k/R), phi-mp.acos(k/R)):
                        xs.add(M(cx)+r*mp.cos(t))
    return sorted(xs)


class Prism:
    def __init__(self, case, plane):
        self.o, self.x, self.y, self.n = (tuple(_Q(v) for v in w) for w in stored_axes(case.frame))
        self.w0, self.w1 = sorted((_Q(case.start), _Q(case.end)))
        q, m = tuple(_Q(v) for v in plane[:3]), tuple(_Q(v) for v in plane[3:6])
        dot = lambda u, v: sum(i*j for i, j in zip(u, v))
        self.a, self.b, self.c = dot(m, self.x), dot(m, self.y), dot(m, self.n)
        self.d = dot(m, tuple(i-j for i, j in zip(self.o, q)))
        self.m = m
        self.els = elements(case.boundaries, case.tolerance)
        self.case = case

    # The plane's value at a point of frame coordinates.
    def g(self, u, v, w):
        return self.a*u+self.b*v+self.c*w+self.d

    def extremes(self):
        """min and max of `a u + b v` over the profile (exact at vertices;
        arcs' support points in 60 digits)."""
        vals = []
        for e in self.els:
            vals.append(M(self.a)*e[1][0]+M(self.b)*e[1][1])
            if e[0] == 'B' and (self.a or self.b):
                # Each piece's end and its support points.
                for piece in e[3]:
                    c = piece.line(M(self.a), M(self.b), M(0))
                    for t in [M(1)]+_roots(_derivative(_trim(c)), touching=True):
                        x, y = piece.point(t)
                        vals.append(M(self.a)*x+M(self.b)*y)
            if e[0] == 'A' and (self.a or self.b):
                _, _, _, (cx, cy, r), sweep = e
                a0, sw = arc_angles(e)
                phi = mp.atan2(self.b, self.a)
                for t in (phi, phi+mp.pi):
                    rel = (t-a0) % (2*mp.pi) if sw > 0 else (a0-t) % (2*mp.pi)
                    if rel <= abs(sw):
                        vals.append(M(self.a)*(cx+r*mp.cos(t))+M(self.b)*(cy+r*mp.sin(t)))
        return min(vals), max(vals)

    def sides(self):
        """Which sides hold solid: below (g < 0) and above (g > 0)."""
        lo, hi = self.extremes()
        ws = [M(self.c)*self.w0, M(self.c)*self.w1]
        gmin = lo+min(ws)+self.d
        gmax = hi+max(ws)+self.d
        eps = M(10)**-35
        return gmin < -eps, gmax > eps

    def regime(self, u, v, below):
        """Which clipped interval holds at `(u, v)`: None (empty), or the
        pair of flags (lower end is the plane, upper end is the plane);
        `below=None` for the whole solid."""
        if below is None:
            return (False, False)
        a, b, c, d = (M(t) for t in (self.a, self.b, self.c, self.d))
        w0, w1 = M(self.w0), M(self.w1)
        base = a*u+b*v+d
        if c == 0:
            inside = base < 0 if below else base > 0
            return (False, False) if inside else None
        h = -base/c
        # g < 0: c w < -base: w < h if c > 0, w > h if c < 0.
        if (c > 0) == below:
            return (False, h < w1) if h > w0 else None
        return (h > w0, False) if h < w1 else None

    def clip(self, u, v, below, regime=None):
        """The clipped interval of `w` at `(u, v)` on one side, as (lo, hi)
        or None; with a regime fixed elsewhere (a sub-interval's midpoint),
        its formula extended to `(u, v)`."""
        r = self.regime(u, v, below) if regime is None else regime
        if r is None:
            return None
        w0, w1 = M(self.w0), M(self.w1)
        c = M(self.c)
        h = -(M(self.a)*u+M(self.b)*v+M(self.d))/c if c != 0 else None
        return (h if r[0] else w0, h if r[1] else w1)

    def volume_moments(self, below):
        a, b, c = (M(t) for t in (self.a, self.b, self.c))
        # Clip lines in (u, v): h = w0, h = w1 (c != 0), or a u + b v + d = 0.
        if c != 0:
            lines = [(a, b, -M(self.d)-c*M(w)) for w in (self.w0, self.w1)]
        else:
            lines = [(a, b, -M(self.d))]
        xs = x_breaks(self.els, lines)
        def inner(x):
            ys = crossings(self.els, x)
            out = [M(0)]*4
            for i in range(0, len(ys)-1, 2):
                y0, y1 = ys[i], ys[i+1]
                # Breaks in v where the clip changes form.
                br = [y0, y1]
                if b != 0:
                    for aa, bb, gg in lines:
                        yv = (gg-aa*x)/bb
                        if y0 < yv < y1:
                            br.append(yv)
                br.sort()
                for s0, s1 in zip(br, br[1:]):
                    r = self.regime(x, (s0+s1)/2, below)
                    def f(y):
                        if r is None:
                            return (0, 0, 0, 0)
                        iv = self.clip(x, y, below, r)
                        L = iv[1]-iv[0]
                        return (L, x*L, y*L, (iv[1]**2-iv[0]**2)/2)
                    fm, f0, f1 = f((s0+s1)/2), f(s0), f(s1)
                    for k in range(4):
                        out[k] += (s1-s0)/6*(f0[k]+4*fm[k]+f1[k])
            return out
        total = [M(0)]*4
        for x0, x1 in zip(xs, xs[1:]):
            if x1-x0 < M(10)**-35:
                continue
            for k in range(4):
                total[k] += mp.quad(lambda x: inner(x)[k], [x0, x1])
        return total

    def area(self, below):
        a, b, c, d = (M(t) for t in (self.a, self.b, self.c, self.d))
        w0, w1 = M(self.w0), M(self.w1)
        # |x cross y| for the caps, |x|, |y| scale the profile's measure; the
        # stored axes are unit up to rounding, taken as exact here.
        area = M(0)
        # Caps: the part of the profile on this side at w0 and at w1 (the
        # whole profile for the whole solid).
        for w in (w0, w1):
            if below is None:
                area += self.region_measure([], lambda u, v: True)
                continue
            lines = [(a, b, -d-c*w)]
            area += self.region_measure(lines, lambda u, v: ((a*u+b*v+c*w+d < 0) == below))
        # Walls: each element's strip, clipped.
        for e in self.els:
            area += self.wall(e, below)
        # The cut face.
        if below is None:
            pass
        elif c != 0:
            lines = [(a, b, -d-c*w) for w in (w0, w1)]
            inside = lambda u, v: w0 < -(a*u+b*v+d)/c < w1
            area += self.region_measure(lines, inside)*mp.sqrt(a*a+b*b+c*c)/abs(c)
        elif a or b:
            area += self.chord((a, b, -d))*(w1-w0)
        return area

    def region_measure(self, lines, inside):
        xs = x_breaks(self.els, lines)
        def inner(x):
            ys = crossings(self.els, x)
            out = M(0)
            for i in range(0, len(ys)-1, 2):
                y0, y1 = ys[i], ys[i+1]
                br = [y0, y1]
                for aa, bb, gg in lines:
                    if bb != 0:
                        yv = (gg-aa*x)/bb
                        if y0 < yv < y1:
                            br.append(yv)
                br.sort()
                for s0, s1 in zip(br, br[1:]):
                    if inside(x, (s0+s1)/2):
                        out += s1-s0
            return out
        total = M(0)
        for x0, x1 in zip(xs, xs[1:]):
            if x1-x0 < M(10)**-35:
                continue
            total += mp.quad(inner, [x0, x1])
        return total

    def chord(self, line):
        """The length of the line `alpha u + beta v = gamma` inside the
        profile."""
        alpha, beta, gamma = line
        norm = mp.sqrt(alpha*alpha+beta*beta)
        # Parameterise the line by arc length s from its foot point.
        foot = (alpha*gamma/norm**2, beta*gamma/norm**2)
        dirv = (-beta/norm, alpha/norm)
        ts = []
        for e in self.els:
            if e[0] == 'L':
                (x0, y0), (x1, y1) = e[1], e[2]
                f0 = alpha*x0+beta*y0-gamma
                f1 = alpha*x1+beta*y1-gamma
                if (f0 < 0) != (f1 < 0) and f0 != f1:
                    t = f0/(f0-f1)
                    px, py = x0+t*(x1-x0), y0+t*(y1-y0)
                    ts.append((px-foot[0])*dirv[0]+(py-foot[1])*dirv[1])
            elif e[0] == 'B':
                for piece in e[3]:
                    for px, py in piece.line_sign_changes(alpha, beta, gamma):
                        ts.append((px-foot[0])*dirv[0]+(py-foot[1])*dirv[1])
            else:
                _, _, _, (cx, cy, r), sweep = e
                a0, sw = arc_angles(e)
                R = norm*r
                k = gamma-alpha*cx-beta*cy
                if abs(k) < R:
                    phi = mp.atan2(beta, alpha)
                    for t in (phi+mp.acos(k/R), phi-mp.acos(k/R)):
                        rel = (t-a0) % (2*mp.pi) if sw > 0 else (a0-t) % (2*mp.pi)
                        if 0 < rel < abs(sw) or abs(abs(sw)-2*mp.pi) < M(10)**-30:
                            px, py = cx+r*mp.cos(t), cy+r*mp.sin(t)
                            ts.append((px-foot[0])*dirv[0]+(py-foot[1])*dirv[1])
        ts.sort()
        return sum(ts[i+1]-ts[i] for i in range(0, len(ts)-1, 2))

    def wall(self, e, below):
        a, b, c, d = (M(t) for t in (self.a, self.b, self.c, self.d))
        def length(u, v):
            iv = self.clip(u, v, below)
            return iv[1]-iv[0] if iv else M(0)
        if e[0] == 'L':
            (x0, y0), (x1, y1) = e[1], e[2]
            seg = mp.sqrt((M(x1)-x0)**2+(M(y1)-y0)**2)
            f = lambda t: length(x0+t*(M(x1)-x0), y0+t*(M(y1)-y0))
            br = [M(0), M(1)]
            if c != 0:
                for w in (self.w0, self.w1):
                    g0 = a*x0+b*y0+d+c*w
                    g1 = a*x1+b*y1+d+c*w
                    if g0 != g1:
                        t = g0/(g0-g1)
                        if 0 < t < 1:
                            br.append(t)
            else:
                g0, g1 = a*x0+b*y0+d, a*x1+b*y1+d
                if g0 != g1:
                    t = g0/(g0-g1)
                    if 0 < t < 1:
                        br.append(t)
            br.sort()
            total = M(0)
            for t0, t1 in zip(br, br[1:]):
                total += (t1-t0)/2*(f(t0+(t1-t0)*M(10)**-30)+f(t1-(t1-t0)*M(10)**-30)) \
                    if c != 0 else (t1-t0)*f((t0+t1)/2)
            return total*seg
        if e[0] == 'B':
            # The clipped height times the speed, between the clip levels'
            # roots on each piece.
            levels = [(a, b, -d-c*M(w)) for w in (self.w0, self.w1)] if c != 0 else [(a, b, -d)]
            total = M(0)
            for piece in e[3]:
                br = {M(0), M(1)}
                for alpha, beta, gamma in levels:
                    if alpha or beta:
                        br.update(piece.line_roots(alpha, beta, gamma))
                br = sorted(br)
                f = lambda t: length(*piece.point(t))*piece.speed(t)
                total += sum(mp.quad(f, [t0, t1]) for t0, t1 in zip(br, br[1:]))
            return total
        _, _, _, (cx, cy, r), sweep = e
        a0, sw = arc_angles(e)
        f = lambda t: length(cx+r*mp.cos(t), cy+r*mp.sin(t))
        lo, hi = (a0, a0+sw) if sw > 0 else (a0+sw, a0)
        br = [lo, hi]
        R = mp.sqrt(a*a+b*b)*r
        phi = mp.atan2(b, a)
        levels = [-d-c*w-a*cx-b*cy for w in ((self.w0, self.w1) if c != 0 else (0,))] if c != 0 \
            else [-d-a*cx-b*cy]
        # The support angles too: quadrature nodes never land on a tangency.
        candidates = [phi, phi+mp.pi] if R > 0 else []
        for k in levels:
            if R > 0 and abs(k) < R:
                candidates += [phi+mp.acos(k/R), phi-mp.acos(k/R)]
        for t in candidates:
            for j in range(-3, 4):
                tt = t+2*mp.pi*j
                if lo < tt < hi:
                    br.append(tt)
        br.sort()
        return sum(mp.quad(f, [t0, t1]) for t0, t1 in zip(br, br[1:]))*r


def rows(case, plane):
    p = Prism(case, plane)
    has_below, has_above = p.sides()
    out = []
    world = lambda u, v, w: tuple(M(p.o[i])+u*M(p.x[i])+v*M(p.y[i])+w*M(p.n[i])
                                  for i in range(3))
    sides = [('below', True), ('above', False)]
    if not (has_below and has_above):
        V, mu, mv, mw = p.volume_moments(None)
        centre = world(mu/V, mv/V, mw/V)
        return [('whole', V, p.area(None), centre)]
    for name, below in sides:
        V, mu, mv, mw = p.volume_moments(below)
        out.append((name, V, p.area(below), world(mu/V, mv/V, mw/V)))
    return out


def number(x):
    x = M(x)
    return '0.0' if abs(x) < M(10)**-30 else mp.nstr(x, 25, min_fixed=-5, max_fixed=5)


def text(row):
    side, V, A, c = row
    return ' '.join(['side', side, number(V), number(A)]+[number(v) for v in c])


# ------------------------------------------------------------------ S8b

def _triple(n, i, j, m, k):
    """The integral over [0, 1] of B(n, i) B(n, j) B(m, k)."""
    from math import comb
    return F(comb(n, i)*comb(n, j)*comb(m, k), comb(2*n+m, i+j+k)*(2*n+m+1))


def green_moments(case):
    """(area, integral of u, integral of v) of a profile of lines and
    splines in its frame's coordinates, by Green's theorem in exact
    Fractions: `A = 1/2 int x dy - y dx`, `int u dA = int x^2/2 dy`, `int v
    dA = -int y^2/2 dx` along the stored boundaries (holes reversed), each
    line a degree-1 Bezier and each spline its pieces from blossoms
    (`Spline.pieces`, checked equal to the slicing's knot insertion), with
    Bernstein product integrals."""
    area = mx = my = F(0)
    for k, b in enumerate(case.boundaries):
        pts, _ = stored(b, case.tolerance)
        assert pts is not None, 'green_moments: lines and splines only'
        points, segs = pts if b.segments is not None else (pts, [None]*len(pts))
        pieces, n = [], len(points)
        for i in range(n):
            p, q = points[i], points[(i+1) % n]
            if segs[i] is None:
                pieces.append(((F(p[0]), F(p[1])), (F(q[0]), F(q[1]))))
            else:
                assert isinstance(segs[i], Spline), 'green_moments: lines and splines only'
                blossomed = segs[i].pieces()
                assert blossomed == bezier_pieces(segs[i]), 'blossoms and knot insertion disagree'
                pieces += blossomed
        if k > 0:
            pieces = [tuple(reversed(c)) for c in reversed(pieces)]
        for c in pieces:
            d = len(c)-1
            dx = [d*(c[j+1][0]-c[j][0]) for j in range(d)]
            dy = [d*(c[j+1][1]-c[j][1]) for j in range(d)]
            for i in range(d+1):
                for j in range(d):
                    area += bernstein_product(d, i, d-1, j)*(c[i][0]*dy[j]-c[i][1]*dx[j])/2
                for j in range(d+1):
                    for l in range(d):
                        w = _triple(d, i, j, d-1, l)
                        mx += w*c[i][0]*c[j][0]*dy[l]/2
                        my -= w*c[i][1]*c[j][1]*dx[l]/2
    return area, mx, my


# ------------------------------------------------------------------ S8c

class Revolved:
    """A cone, frustum or sphere zone (S8c) about its stored frame's normal:
    the points `u x + v y + w n` (stored axes, `stored_axes`) with `u^2 +
    v^2 <= rho(w)^2`, `w` between its ends; a cone `rho = r1 + (r2 - r1) w /
    h` on `[0, h]`, a zone `rho = sqrt(R^2 - w^2)` on `[R sin a1, R sin a2]`.
    The plane `a u + b v + c w + d` has exact coefficients from the stored
    data, as for prisms.

    Each slice `w` is a disc cut by the line `a u + b v = -(c w + d)`: its
    part below is the disc less the circular segment beyond `s = -(c w +
    d) / |(a, b)|`, with area `pi rho^2 - (rho^2 acos(s/rho) - s sqrt(rho^2
    - s^2))` and moment `-(2/3)(rho^2 - s^2)^(3/2)` along `(a, b)`. Volumes
    and moments integrate these over `w` (`mp.quad` between the ends and
    every `w` where `|s| = rho`). Areas: the ends' discs' parts, the lateral
    surface's `rho sqrt(1 + rho'^2)` times the angle below (`2 pi - 2
    acos(s/rho)`), and the cut face, `|m| / |(a, b)|` times the integral of
    its chords `2 sqrt(rho^2 - s^2)` (or the disc at `w = -d/c` when the
    plane is normal to the axis)."""

    def __init__(self, kind, frame, params, plane):
        self.kind = kind
        o, x, y, n = stored_axes(frame)
        self.o, self.x, self.y, self.n = o, x, y, n
        if kind == 'cone':
            r1, r2, h = params
            self.w0, self.w1 = M(0), M(h)
            self.rho = lambda w: M(r1)+(M(r2)-M(r1))*w/M(h)
            self.slope = (M(r2)-M(r1))/M(h)
        else:
            R, a1, a2 = params
            self.R = M(R)
            lo = -self.R if a1 == -math.pi/2 else self.R*mp.sin(M(a1))
            hi = self.R if a2 == math.pi/2 else self.R*mp.sin(M(a2))
            self.w0, self.w1 = lo, hi
            self.rho = lambda w: mp.sqrt(max(self.R**2-w**2, M(0)))
        q, m = [_Q(t) for t in plane[:3]], [_Q(t) for t in plane[3:]]
        dot = lambda u, v: sum(_Q(u[i])*_Q(v[i]) for i in range(3))
        self.a, self.b, self.c = dot(m, x), dot(m, y), dot(m, n)
        self.d = sum(m[i]*(_Q(o[i])-q[i]) for i in range(3))
        self.norm = mp.sqrt(sum(M(t)**2 for t in m))
        self.ab = mp.sqrt(M(self.a)**2+M(self.b)**2)

    def s(self, w):
        return -(M(self.c)*w+M(self.d))/self.ab

    def breaks(self):
        """The ends, and every w where the slice's line touches its disc."""
        out = [self.w0, self.w1]
        if self.ab == 0:
            if self.c != 0:
                ws = -M(self.d)/M(self.c)
                if self.w0 < ws < self.w1:
                    out.append(ws)
            return sorted(out)
        a, b, c, d, ab = M(self.a), M(self.b), M(self.c), M(self.d), self.ab
        roots = []
        if self.kind == 'cone':
            # -(c w + d)/ab = +-(r1 + k w).
            r1, k = self.rho(M(0)), self.slope
            for sg in (1, -1):
                den = c/ab+sg*k
                if den != 0:
                    roots.append((-d/ab-sg*r1)/den)
        else:
            # (c w + d)^2 = ab^2 (R^2 - w^2).
            A, B, C = c*c+ab*ab, 2*c*d, d*d-ab*ab*self.R**2
            disc = B*B-4*A*C
            if disc >= 0:
                roots += [(-B+sg*mp.sqrt(disc))/(2*A) for sg in (1, -1)]
        out += [r for r in roots if self.w0 < r < self.w1]
        return sorted(out)

    def slice(self, w, below):
        """(area, moment along (a, b)/|(a, b)|) of the slice's part on a side."""
        rho = self.rho(w)
        full = mp.pi*rho**2
        if self.ab == 0:
            inside = (M(self.c)*w+M(self.d) < 0) == below
            return (full if inside else M(0)), M(0)
        s = self.s(w)
        if s >= rho:
            part, mom = full, M(0)
        elif s <= -rho:
            part, mom = M(0), M(0)
        else:
            h2 = rho**2-s**2
            part = full-(rho**2*mp.acos(s/rho)-s*mp.sqrt(h2))
            mom = -M(2)/3*h2**M(1.5)
        return (part, mom) if below else (full-part, -mom)

    def quad(self, f):
        br = self.breaks()
        return sum(mp.quad(f, [t0, t1]) for t0, t1 in zip(br, br[1:]))

    def volume_moments(self, below):
        a, b = M(self.a), M(self.b)
        e = (a/self.ab, b/self.ab) if self.ab != 0 else (M(0), M(0))
        V = self.quad(lambda w: self.slice(w, below)[0])
        mx = self.quad(lambda w: self.slice(w, below)[1])
        mw = self.quad(lambda w: w*self.slice(w, below)[0])
        return V, mx*e[0], mx*e[1], mw

    def angle(self, w, below):
        """The angle of the slice's circle on a side."""
        rho = self.rho(w)
        if self.ab == 0:
            return 2*mp.pi if (M(self.c)*w+M(self.d) < 0) == below else M(0)
        s = self.s(w)
        if rho == 0:
            return M(0)
        t = 2*mp.pi if s >= rho else (M(0) if s <= -rho else 2*mp.pi-2*mp.acos(s/rho))
        return t if below else 2*mp.pi-t

    def area(self, below):
        ends = sum(self.slice(w, below)[0] for w in (self.w0, self.w1) if self.rho(w) > 0)
        if self.kind == 'cone':
            factor = mp.sqrt(1+self.slope**2)
            lateral = self.quad(lambda w: self.rho(w)*factor*self.angle(w, below))
        else:
            lateral = self.quad(lambda w: self.R*self.angle(w, below))
        if self.ab == 0:
            ws = -M(self.d)/M(self.c)
            cut = mp.pi*self.rho(ws)**2 if self.w0 < ws < self.w1 else M(0)
        else:
            chord = lambda w: 2*mp.sqrt(max(self.rho(w)**2-self.s(w)**2, M(0)))
            cut = self.norm/self.ab*self.quad(chord)
        return ends+lateral+cut

    def world(self, u, v, w):
        return tuple(M(self.o[i])+u*M(self.x[i])+v*M(self.y[i])+w*M(self.n[i]) for i in range(3))


def revolved_rows(kind, frame, params, plane):
    p = Revolved(kind, frame, params, plane)
    rows, total = [], p.volume_moments(True)[0]+p.volume_moments(False)[0]
    sides = []
    for name, below in (('below', True), ('above', False)):
        V, mu, mv, mw = p.volume_moments(below)
        if V > total*M(10)**-25:
            sides.append((name, V, p.area(below), p.world(mu/V, mv/V, mw/V)))
    if len(sides) < 2:
        V, mu, mv, mw = [x+y for x, y in zip(p.volume_moments(True), p.volume_moments(False))]
        return [('whole', V, p.area(True)+p.area(False)-2*_cut(p), p.world(mu/V, mv/V, mw/V))]
    return sides


def _cut(p):
    """The cut face's area (counted on both sides of a whole solid)."""
    if p.ab == 0:
        ws = -M(p.d)/M(p.c)
        return mp.pi*p.rho(ws)**2 if p.w0 < ws < p.w1 else M(0)
    chord = lambda w: 2*mp.sqrt(max(p.rho(w)**2-p.s(w)**2, M(0)))
    return p.norm/p.ab*p.quad(chord)


# ------------------------------------------------------------------ S8d

class Torus(Revolved):
    """A whole ring torus (S8d) about its stored frame's normal: the points
    with `(sqrt(u^2 + v^2) - R)^2 + w^2 <= r^2`. Each slice `|w| < r` is an
    annulus between `R -+ sqrt(r^2 - w^2)`: the outer disc less the inner
    one, both cut by the plane's line. The lateral area integrates the
    angle below over the tube's angle (the area element `r (R + r cos t)`),
    and the cut face's chords are the outer disc's less the inner's."""

    def __init__(self, frame, params, plane):
        R, r = params
        Revolved.__init__(self, 'cone', frame, (1.0, 2.0, 1.0), plane)
        self.kind = 'torus'
        self.R, self.r = M(R), M(r)
        self.w0, self.w1 = -self.r, self.r
        half = lambda w: mp.sqrt(max(self.r**2-w**2, M(0)))
        self.outer = lambda w: self.R+half(w)
        self.inner = lambda w: self.R-half(w)
        self.rho = self.outer

    def breaks(self):
        out = [self.w0, self.w1]
        if self.ab == 0:
            if self.c != 0:
                ws = -M(self.d)/M(self.c)
                if self.w0 < ws < self.w1:
                    out.append(ws)
            return sorted(out)
        # Where |s| meets either radius: sign changes on a fine grid,
        # refined.
        n = 4000
        grid = [self.w0+(self.w1-self.w0)*k/n for k in range(n+1)]
        for f in (lambda w: abs(self.s(w))-self.outer(w), lambda w: abs(self.s(w))-self.inner(w)):
            for a, b in zip(grid, grid[1:]):
                fa, fb = f(a), f(b)
                if fa == 0:
                    out.append(a)
                elif fa*fb < 0:
                    out.append(mp.findroot(f, (a, b), solver='anderson'))
        return sorted(set(out))

    def disc(self, rho, w, below):
        """(area, moment) of a disc of radius rho's part on a side."""
        full = mp.pi*rho**2
        if self.ab == 0:
            inside = (M(self.c)*w+M(self.d) < 0) == below
            return (full if inside else M(0)), M(0)
        s = self.s(w)
        if s >= rho:
            part, mom = full, M(0)
        elif s <= -rho:
            part, mom = M(0), M(0)
        else:
            h2 = rho**2-s**2
            part = full-(rho**2*mp.acos(s/rho)-s*mp.sqrt(h2))
            mom = -M(2)/3*h2**M(1.5)
        return (part, mom) if below else (full-part, -mom)

    def slice(self, w, below):
        a, ma = self.disc(self.outer(w), w, below)
        b, mb = self.disc(self.inner(w), w, below)
        return a-b, ma-mb

    def area(self, below):
        # The tube's angle t: the circle of radius R + r cos t at height r
        # sin t, its angle below the plane.
        def lateral(t):
            rho, w = self.R+self.r*mp.cos(t), self.r*mp.sin(t)
            if self.ab == 0:
                ang = 2*mp.pi if (M(self.c)*w+M(self.d) < 0) == below else M(0)
            else:
                s = self.s(w)
                a = 2*mp.pi if s >= rho else (M(0) if s <= -rho else 2*mp.pi-2*mp.acos(s/rho))
                ang = a if below else 2*mp.pi-a
            return self.r*rho*ang
        ts = sorted(set([M(0), 2*mp.pi]+[t % (2*mp.pi) for w in self.breaks()
                                          for t in (mp.asin(max(min(w/self.r, M(1)), M(-1))),
                                                    mp.pi-mp.asin(max(min(w/self.r, M(1)), M(-1))))]))
        lat = sum(mp.quad(lateral, [t0, t1]) for t0, t1 in zip(ts, ts[1:]))
        return lat+self.cut_area()

    def cut_area(self):
        if self.ab == 0:
            ws = -M(self.d)/M(self.c)
            return mp.pi*(self.outer(ws)**2-self.inner(ws)**2) if self.w0 < ws < self.w1 else M(0)
        def chord(w):
            s = self.s(w)
            o = 2*mp.sqrt(max(self.outer(w)**2-s**2, M(0)))
            i = 2*mp.sqrt(max(self.inner(w)**2-s**2, M(0)))
            return o-i
        return self.norm/self.ab*self.quad(chord)


def torus_rows(frame, params, plane):
    p = Torus(frame, params, plane)
    total = p.volume_moments(True)[0]+p.volume_moments(False)[0]
    sides = []
    for name, below in (('below', True), ('above', False)):
        V, mu, mv, mw = p.volume_moments(below)
        if V > total*M(10)**-25:
            sides.append((name, V, p.area(below), p.world(mu/V, mv/V, mw/V)))
    if len(sides) < 2:
        V, mu, mv, mw = [x+y for x, y in zip(p.volume_moments(True), p.volume_moments(False))]
        return [('whole', V, p.area(True)+p.area(False)-2*p.cut_area(), p.world(mu/V, mv/V, mw/V))]
    return sides


# ------------------------------------------------------------------ S8e

def _sign(x):
    return (x > 0)-(x < 0)


class Planar:
    """A face or wire body (S6: a case with `make`) on its stored frame's
    plane `w = 0`. The plane's trace there is `g = a u + b v + d = 0`, the
    exact rationals of `Prism` (`c` plays no part); `below` is `g < 0`.

    A sheet's side is the profile's section by the trace: its area and first
    moments are a prism's of height one by a plane parallel to its axis (the
    same slicing), its perimeter the length of its boundary strictly on the
    side plus the parts of the trace bounding it (between consecutive
    contacts of the trace with the boundary, where a point just off the
    trace on that side is inside the profile). A wire's pieces are its first
    boundary's elements cut where `g` vanishes, in stored order: lines
    exactly (the crossing a rational point), arcs by exact tangency (`(a cx +
    b cy + d)^2` against `r^2 (a^2 + b^2)`) and their crossing angles,
    splines at their Bezier pieces' roots; a piece lying along the trace
    takes the side of the piece before it in stored order (S8e's
    decisions). Lengths and moments: lines and arcs in closed form, splines
    by `mp.quad` of the speed."""

    def __init__(self, case, plane):
        import dataclasses
        assert case.make in ('face', 'wire')
        boundaries = case.boundaries if case.make == 'face' else case.boundaries[:1]
        solid = dataclasses.replace(case, make=None, start=0.0, end=1.0, boundaries=boundaries)
        self.prism = Prism(solid, plane)
        self.prism.c = F(0)
        self.make = case.make
        self.a, self.b, self.d = self.prism.a, self.prism.b, self.prism.d
        self.els = self.prism.els

    def g(self, u, v):
        return M(self.a)*u+M(self.b)*v+M(self.d)

    def world(self, u, v):
        p = self.prism
        return tuple(M(p.o[i])+u*M(p.x[i])+v*M(p.y[i]) for i in range(3))

    def pieces(self, e):
        """The element's pieces in its direction, each `(sign, length, mu,
        mv, start, end)`: cut where `g` changes sign and where it touches
        zero (an arc's tangency, a spline's even root), `sign` 0 for a line
        lying along the trace."""
        a, b, d = self.a, self.b, self.d
        if e[0] == 'L':
            p, q = [tuple(F(c) for c in pt) for pt in (e[1], e[2])]
            f0, f1 = a*p[0]+b*p[1]+d, a*q[0]+b*q[1]+d
            parts = [(p, q)]
            if f0*f1 < 0:
                t = f0/(f0-f1)
                r = (p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1]))
                parts = [(p, r), (r, q)]
            out = []
            for s, t in parts:
                fm = a*(s[0]+t[0])/2+b*(s[1]+t[1])/2+d
                length = mp.sqrt(M((t[0]-s[0])**2+(t[1]-s[1])**2))
                mid = (M(s[0]+t[0])/2, M(s[1]+t[1])/2)
                out.append((_sign(fm), length, length*mid[0], length*mid[1],
                            (M(s[0]), M(s[1])), (M(t[0]), M(t[1]))))
            return out
        if e[0] == 'A':
            _, _, _, (cx, cy, r), _ = e
            a0, sw = arc_angles(e)
            k = -(a*F(cx)+b*F(cy)+d)
            R2 = F(r)**2*(a*a+b*b)
            rels = []
            if a or b:
                phi = mp.atan2(M(b), M(a))
                if k*k < R2:
                    alpha = mp.acos(M(k)/mp.sqrt(M(R2)))
                    roots = [phi+alpha, phi-alpha]
                elif k*k == R2:
                    roots = [phi if k > 0 else phi+mp.pi]
                else:
                    roots = []
                for t in roots:
                    rel = (t-a0) % (2*mp.pi) if sw > 0 else (a0-t) % (2*mp.pi)
                    if M(10)**-30 < rel < abs(sw)-M(10)**-30:
                        rels.append(rel)
            rels = [M(0)]+sorted(rels)+[abs(sw)]
            direction = 1 if sw > 0 else -1
            point = lambda t: (M(cx)+r*mp.cos(t), M(cy)+r*mp.sin(t))
            out = []
            for r0, r1 in zip(rels, rels[1:]):
                t0, t1 = a0+direction*r0, a0+direction*r1
                lo, hi = min(t0, t1), max(t0, t1)
                sign = _sign(self.g(*point((t0+t1)/2)))
                length = M(r)*(hi-lo)
                mu = M(r)*(M(cx)*(hi-lo)+M(r)*(mp.sin(hi)-mp.sin(lo)))
                mv = M(r)*(M(cy)*(hi-lo)-M(r)*(mp.cos(hi)-mp.cos(lo)))
                # The element's own ends exactly: its stored points.
                start = (M(e[1][0]), M(e[1][1])) if r0 == 0 else point(t0)
                end = (M(e[2][0]), M(e[2][1])) if r1 == abs(sw) else point(t1)
                out.append((sign, length, mu, mv, start, end))
            return out
        out = []
        for piece in e[3]:
            ts = [M(0)]+(piece.line_roots(a, b, -d) if a or b else [])+[M(1)]
            for t0, t1 in zip(ts, ts[1:]):
                sign = _sign(self.g(*piece.point((t0+t1)/2)))
                assert sign != 0, 'a spline along the trace'
                length = mp.quad(piece.speed, [t0, t1])
                mu = mp.quad(lambda t: piece.point(t)[0]*piece.speed(t), [t0, t1])
                mv = mp.quad(lambda t: piece.point(t)[1]*piece.speed(t), [t0, t1])
                out.append((sign, length, mu, mv, piece.point(t0), piece.point(t1)))
        return out

    def trace_length(self, sign):
        """The length of the trace bounding the side `sign` of the profile:
        between consecutive contacts of the trace with the boundary, where a
        point just off the trace's midpoint on that side is inside."""
        a, b, d = M(self.a), M(self.b), M(self.d)
        norm = mp.sqrt(a*a+b*b)
        foot = (-a*d/norm**2, -b*d/norm**2)
        dirv = (-b/norm, a/norm)
        ts = []
        for e in self.els:
            for piece in self.pieces(e):
                for pt in piece[4:6]:
                    if abs(self.g(*pt)) <= M(10)**-28:
                        ts.append((pt[0]-foot[0])*dirv[0]+(pt[1]-foot[1])*dirv[1])
        ts.sort()
        kept = [t for i, t in enumerate(ts) if i == 0 or t-ts[i-1] > M(10)**-30]
        delta = M(10)**-25
        total = M(0)
        for t0, t1 in zip(kept, kept[1:]):
            tm = (t0+t1)/2
            x = foot[0]+tm*dirv[0]+sign*delta*a/norm
            y = foot[1]+tm*dirv[1]+sign*delta*b/norm
            if sum(1 for v in crossings(self.els, x) if v < y) % 2 == 1:
                total += t1-t0
        return total

    def perimeter(self, sign=None):
        """The boundary's length on a side (strictly), or all of it."""
        return sum(q[1] for e in self.els for q in self.pieces(e) if sign is None or q[0] == sign)

    def sheet_rows(self):
        p = self.prism
        has_below, has_above = p.sides()
        if not (has_below and has_above):
            V, mu, mv, _ = p.volume_moments(None)
            return [('whole', V, self.perimeter(), self.world(mu/V, mv/V))]
        out = []
        for name, below, s in (('below', True, -1), ('above', False, 1)):
            V, mu, mv, _ = p.volume_moments(below)
            out.append((name, V, self.perimeter(s)+self.trace_length(s), self.world(mu/V, mv/V)))
        return out

    def wire_pieces(self):
        """The wire's pieces in stored order, each along the trace given the
        side of the piece before it (cyclically); every piece `-1` (below)
        when the whole wire lies along the trace."""
        seq = [q for e in self.els for q in self.pieces(e)]
        if all(q[0] == 0 for q in seq):
            return [(-1,)+q[1:] for q in seq]
        out = []
        for i, q in enumerate(seq):
            j = i
            while seq[j][0] == 0:
                j = (j-1) % len(seq)
            out.append((seq[j][0],)+q[1:])
        return out

    def wire_runs(self):
        """{side: its number of open wires}: the maximal cyclic runs of
        pieces on one side, or {'whole': 1}."""
        seq = [q[0] for q in self.wire_pieces()]
        if len(set(seq)) == 1:
            return {'whole': 1}
        runs = {}
        for i, s in enumerate(seq):
            if seq[i-1] != s:
                key = 'below' if s < 0 else 'above'
                runs[key] = runs.get(key, 0)+1
        return runs

    def wire_rows(self):
        pieces = self.wire_pieces()
        def row(name, chosen):
            L = sum(q[1] for q in chosen)
            return (name, L, M(0), self.world(sum(q[2] for q in chosen)/L, sum(q[3] for q in chosen)/L))
        if len(set(q[0] for q in pieces)) == 1:
            return [row('whole', pieces)]
        return [row('below', [q for q in pieces if q[0] < 0]), row('above', [q for q in pieces if q[0] > 0])]

    def rows(self):
        return self.sheet_rows() if self.make == 'face' else self.wire_rows()


def planar_rows(case, plane):
    """S8e: a face body's sides (`side S area perimeter cx cy cz`) or a wire
    body's (`side S length 0 cx cy cz`); a trace missing or touching the
    body, or a plane parallel to it, gives one `whole` row."""
    return Planar(case, plane).rows()
