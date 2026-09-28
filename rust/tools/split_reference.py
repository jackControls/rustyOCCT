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

from identity_reference import stored, arc_sweep
from curve_surface_reference import stored_axes

mp.mp.dps = 40


def _Q(v):
    return F(v)


def M(v):
    """An mpf of a Fraction, a float or an mpf."""
    if isinstance(v, F):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def elements(boundaries, tolerance):
    """The stored boundaries' elements: ('L', p, q) or ('A', p, q, (cx, cy,
    r), sweep), counter-clockwise outer, clockwise holes (as stored: holes
    reversed so the region is on the left)."""
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
