#!/usr/bin/env python3
"""Independent reference for S9a of REVIEW_NOTES.md: Booleans of two prisms
in one frame (fuse, cut, common).

Both prisms' frames have bitwise-equal stored axes (`stored_axes`, the
kernel's `Frame3::new`), and the tool's origin is the object's plus `a x +
b y + c n`, `(a, b, c)` solved exactly from the stored binary64 data (the
axes as stored, rationals, not exactly orthonormal) and required to be
binary64 (S9a's domain). In the object's frame coordinates `(u, v, w)` the
tool is its profile translated by `(a, b)` exactly (Fractions) and its
offsets moved by `c`. Nothing here uses the kernel or a 2D Boolean library.

The Boolean is a stack of height slabs between the four caps' heights; a
slab's region is the 2D Boolean of the profiles present in it (fuse: `A u
B`, cut: `A \\ B`, common: `A n B`; a profile absent from the slab is
empty). Every such region is a union of the three atoms `AB` (inside both
profiles), `A` (inside the object's only) and `B` (the tool's only), so
everything reduces to the atoms' areas and first moments and to the lengths
of the profiles' boundary pieces by their classes:

* **Atoms, by slicing in `v`.** The breakpoints are every vertex's `v`,
  every circle's `cy +- r`, and every meeting of the underlying curves of
  the two boundaries (lines exactly in Fractions; a line and a circle, and
  two circles by their radical line, from exact discriminants). Between two
  breakpoints the order of the boundaries' crossings of `v = const` is
  fixed; at the band's middle each elementary interval between consecutive
  crossings belongs to one atom (by the even-odd counts of each profile's
  crossings) and its ends to known curves: a line `u = alpha v + beta`
  (canonical Fractions, so coincident edges of the two profiles are the same
  curve) or a circle's branch `u = cx +- sqrt(r^2 - (v - cy)^2)`. The
  band's contribution to an atom's area and first moments is the closed-form
  integral over the band of the ends' `u`, `u^2 / 2` and `u v` (lines as
  polynomials, circles by `asin` and `(r^2 - t^2)^(3/2)`).
* **Boundary pieces.** Each element (line or arc) of each profile is cut at
  its meetings with the other profile's elements and at the other's
  vertices lying on it; each piece is classified at its midpoint against the
  other profile: `on` it (within 1e-25 of the case's size), with the same or
  the opposite direction (boundaries oriented with their region on the
  left: outer boundaries counter-clockwise, holes clockwise), else inside or
  outside by the parity of a ray's crossings (a direction of transcendental
  slope, retried if it grazes a vertex or a tangent). A piece of `A` outside
  `B` bounds the atom `A` against nothing, one inside `B` bounds `AB`
  against `B`, and likewise for `B`'s; a shared piece of the same direction
  bounds `AB` against nothing, one of opposite directions `A` against `B`.
  A region's boundary length is the sum of the classes whose two sides it
  holds one of.
* **The result.** Volume and first moments: each slab's region's area and
  moments times its height. Surface area: each slab's region's boundary
  length times its height (walls), and at each slab boundary the area of
  the symmetric difference of the regions below and above (caps, including
  the bottom and the top). Centre: the moments over the volume, mapped to
  world coordinates by the object's stored axes. Solids: the maximal
  connected regions, by union-find over each slab's intervals in each
  band (intervals of adjacent bands joined when their limits at the shared
  `v` overlap in a positive length; of adjacent slabs, when they overlap at
  the band's middle): regions meeting only at a point or along an edge are
  separate solids (the regularized Boolean). A result of zero volume is
  empty: no solids.

Consecutive slabs with equal regions (a zero symmetric difference) are
merged, empty ones dropped: the result is one prism per solid (S9a.1) when
one run of slabs remains, otherwise a stack (S9a.2).
"""
from fractions import Fraction as F
import math

import mpmath as mp

from identity_reference import arc_sweep, stored
from curve_surface_reference import stored_axes

mp.mp.dps = 40

ATOMS = ('AB', 'A', 'B')


def M(v):
    if isinstance(v, F):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def region(operation, a_present, b_present):
    """The atoms of a slab's region."""
    if a_present and b_present:
        return {'fuse': {'AB', 'A', 'B'}, 'cut': {'A'}, 'common': {'AB'}}[operation]
    if a_present:
        return {'AB', 'A'} if operation in ('fuse', 'cut') else set()
    if b_present:
        return {'AB', 'B'} if operation == 'fuse' else set()
    return set()


# ------------------------------------------------------------------ elements

class Line:
    kind = 'L'

    def __init__(self, p, q):
        self.p, self.q = p, q
        self.d = (q[0]-p[0], q[1]-p[1])
        assert self.d != (0, 0)

    def key(self):
        """The canonical curve `u = alpha v + beta` (None when horizontal)."""
        if self.d[1] == 0:
            return None
        alpha = self.d[0]/self.d[1]
        return ('L', alpha, self.p[0]-alpha*self.p[1])

    def length(self):
        return mp.sqrt(M(self.d[0]**2+self.d[1]**2))

    def point(self, t):
        return (M(self.p[0])+t*M(self.d[0]), M(self.p[1])+t*M(self.d[1]))

    def tangent(self, t):
        return (M(self.d[0]), M(self.d[1]))

    def param(self, pt):
        dd = M(self.d[0]**2+self.d[1]**2)
        return ((pt[0]-M(self.p[0]))*M(self.d[0])+(pt[1]-M(self.p[1]))*M(self.d[1]))/dd

    def crossings(self, v):
        lo, hi = sorted((M(self.p[1]), M(self.q[1])))
        if lo < v < hi:
            k = self.key()
            return [(x_at(k, v), k)]
        return []

    def reversed(self):
        return Line(self.q, self.p)


class Arc:
    kind = 'A'

    def __init__(self, c, r, p, q, sweep, full):
        self.c, self.r, self.p, self.q, self.sweep, self.full = c, r, p, q, sweep, full
        self.theta0 = mp.atan2(M(p[1]-c[1]), M(p[0]-c[0]))

    def length(self):
        return M(self.r)*abs(self.sweep)

    def angle(self, t):
        return self.theta0+t*self.sweep

    def point(self, t):
        th = self.angle(t)
        return (M(self.c[0])+M(self.r)*mp.cos(th), M(self.c[1])+M(self.r)*mp.sin(th))

    def tangent(self, t):
        th, s = self.angle(t), (1 if self.sweep > 0 else -1)
        return (-s*mp.sin(th), s*mp.cos(th))

    def param_of_angle(self, phi):
        """The parameter of the angle `phi` on the arc, in [0, 2 pi / |sweep|);
        values near the full turn wrap to just below 0."""
        two_pi = 2*mp.pi
        rel = (phi-self.theta0) % two_pi if self.sweep > 0 else (self.theta0-phi) % two_pi
        if rel > two_pi-M(10)**-28:
            rel -= two_pi
        return rel/abs(self.sweep)

    def param(self, pt):
        return self.param_of_angle(mp.atan2(pt[1]-M(self.c[1]), pt[0]-M(self.c[0])))

    def crossings(self, v):
        t = v-M(self.c[1])
        r = M(self.r)
        if abs(t) >= r:
            return []
        s = mp.sqrt(r*r-t*t)
        out = []
        for sigma in (1, -1):
            if self.full or 0 < self.param_of_angle(mp.atan2(t, sigma*s)) < 1:
                k = ('C', self.c[0], self.c[1], self.r, sigma)
                out.append((x_at(k, v), k))
        return out

    def reversed(self):
        return Arc(self.c, self.r, self.q, self.p, -self.sweep, self.full)


def profile_elements(boundaries, tolerance, du=F(0), dv=F(0)):
    """The stored boundaries (`identity_reference.stored`) translated by
    `(du, dv)` exactly: counter-clockwise outer boundaries, holes reversed
    (clockwise) so the region lies on the left of every element."""
    out = []
    for k, b in enumerate(boundaries):
        pts, _ = stored(b, tolerance)
        items = []
        move = lambda p: (F(p[0])+du, F(p[1])+dv)
        if pts is None:
            cx, cy, r = b.circle
            c = move((cx, cy))
            p = (c[0]+F(r), c[1])
            items.append(Arc(c, F(r), p, p, 2*mp.pi, True))
        elif b.segments is not None:
            points, segs = pts
            n = len(points)
            for i in range(n):
                p, q = points[i], points[(i+1) % n]
                if segs[i] is None:
                    items.append(Line(move(p), move(q)))
                else:
                    assert isinstance(segs[i], tuple), 'S9a: lines, arcs and circles'
                    cx, cy, r, _ = segs[i]
                    items.append(Arc(move((cx, cy)), F(r), move(p), move(q), arc_sweep(p, q, segs[i]), False))
        else:
            n = len(pts)
            for i in range(n):
                items.append(Line(move(pts[i]), move(pts[(i+1) % n])))
        if k > 0:
            items = [e.reversed() for e in reversed(items)]
        out += items
    return out


# ------------------------------------------------------------------ curves

def x_at(key, v):
    if key[0] == 'L':
        return M(key[1])*v+M(key[2])
    _, cx, cy, r, sigma = key
    t = v-M(cy)
    return M(cx)+sigma*mp.sqrt(max(M(0), M(r)**2-t*t))


def integrals(key, v0, v1):
    """The integrals over [v0, v1] of u, u^2 / 2 and u v along the curve."""
    if key[0] == 'L':
        a, b = M(key[1]), M(key[2])
        d1, d2, d3 = v1-v0, (v1**2-v0**2)/2, (v1**3-v0**3)/3
        return (a*d2+b*d1, (a*a*d3+2*a*b*d2+b*b*d1)/2, a*d3+b*d2)
    _, cx, cy, r, sigma = key
    cx, cy, r = M(cx), M(cy), M(r)

    def j0(t):
        t = min(max(t, -r), r)
        return (t*mp.sqrt(max(M(0), r*r-t*t))+r*r*mp.asin(t/r))/2

    def j1(t):
        t = min(max(t, -r), r)
        return -max(M(0), r*r-t*t)**mp.mpf(1.5)/3
    t0, t1 = v0-cy, v1-cy
    J0, J1 = j0(t1)-j0(t0), j1(t1)-j1(t0)
    dt, dt2, dt3 = t1-t0, (t1**2-t0**2)/2, (t1**3-t0**3)/3
    return (cx*dt+sigma*J0, ((cx*cx+r*r)*dt-dt3)/2+sigma*cx*J0, cx*(dt2+cy*dt)+sigma*(J1+cy*J0))


# ------------------------------------------------------------------ meetings

def _cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def meet(e, f):
    """The points (mpf pairs) where the underlying curves of `e` and `f`
    meet (none for parallel lines or concentric circles: coincident pieces
    are cut at the other's vertices instead)."""
    if e.kind == 'L' and f.kind == 'L':
        den = _cross(e.d, f.d)
        if den == 0:
            return []
        t = _cross((f.p[0]-e.p[0], f.p[1]-e.p[1]), f.d)/den
        return [(M(e.p[0]+t*e.d[0]), M(e.p[1]+t*e.d[1]))]
    if e.kind == 'A' and f.kind == 'L':
        return meet(f, e)
    if e.kind == 'L':
        c, r = f.c, f.r
        w = (e.p[0]-c[0], e.p[1]-c[1])
        a = e.d[0]**2+e.d[1]**2
        b = 2*(e.d[0]*w[0]+e.d[1]*w[1])
        cc = w[0]**2+w[1]**2-r*r
        disc = b*b-4*a*cc
        if disc < 0:
            return []
        roots = [M(-b)/(2*M(a))] if disc == 0 else \
            [(M(-b)+s*mp.sqrt(M(disc)))/(2*M(a)) for s in (1, -1)]
        return [e.point(t) for t in roots]
    (c1, r1), (c2, r2) = (e.c, e.r), (f.c, f.r)
    dx, dy = c2[0]-c1[0], c2[1]-c1[1]
    d2 = dx*dx+dy*dy
    if d2 == 0:
        return []
    k = (d2+r1*r1-r2*r2)/(2*d2)
    m = r1*r1/d2-k*k
    if m < 0:
        return []
    base = (M(c1[0]+k*dx), M(c1[1]+k*dy))
    if m == 0:
        return [base]
    s = mp.sqrt(M(m))
    return [(base[0]-s*M(dy), base[1]+s*M(dx)), (base[0]+s*M(dy), base[1]-s*M(dx))]


def distance(pt, f):
    if f.kind == 'L':
        t = min(max(f.param(pt), M(0)), M(1))
        q = f.point(t)
        return mp.sqrt((pt[0]-q[0])**2+(pt[1]-q[1])**2)
    t = f.param(pt)
    if f.full or 0 <= t <= 1:
        return abs(mp.sqrt((pt[0]-M(f.c[0]))**2+(pt[1]-M(f.c[1]))**2)-M(f.r))
    return min(mp.sqrt((pt[0]-q[0])**2+(pt[1]-q[1])**2) for q in (f.point(0), f.point(1)))


def within(e, pt, eps):
    t = e.param(pt)
    if e.kind == 'A' and e.full:
        return True
    return -eps <= t <= 1+eps


# ------------------------------------------------------------------ the pair

class Pair:
    """The object's and the tool's profiles in the object's frame, their
    atoms and boundary classes."""

    def __init__(self, obj, tool):
        self.axes = tuple(tuple(F(c) for c in v) for v in stored_axes(obj.frame))
        taxes = tuple(tuple(F(c) for c in v) for v in stored_axes(tool.frame))
        assert self.axes[1:] == taxes[1:], 'S9a: frames with bitwise-equal axes'
        o, x, y, n = self.axes
        D = tuple(taxes[0][i]-o[i] for i in range(3))
        det = lambda a, b, c: (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
                               + a[2]*(b[0]*c[1]-b[1]*c[0]))
        # D = a x + b y + c n, by Cramer's rule in Fractions.
        den = det(x, y, n)
        a, b, c = det(D, y, n)/den, det(x, D, n)/den, det(x, y, D)/den
        for value in (a, b, c):
            assert F(float(value)) == value, 'S9a: the origins\' offset must be binary64 in the frame'
        self.offset = (a, b, c)
        self.A = profile_elements(obj.boundaries, obj.tolerance)
        self.B = profile_elements(tool.boundaries, tool.tolerance, a, b)
        self.heights = {'A': tuple(sorted((F(obj.start), F(obj.end)))),
                        'B': tuple(sorted((F(tool.start)+c, F(tool.end)+c)))}
        coords = [abs(M(v)) for e in self.A+self.B for v in (*e.p, *e.q)]
        coords += [abs(M(e.c[i]))+M(e.r) for e in self.A+self.B if e.kind == 'A' for i in range(2)]
        self.scale = max([M(1)]+coords)
        self.eps = M(10)**-25*self.scale
        self._bands()
        self.classes = self._classes()

    # -------------------------------------------------------------- slicing

    def _breaks(self):
        vs = set()
        for e in self.A+self.B:
            vs.add(M(e.p[1]))
            vs.add(M(e.q[1]))
            if e.kind == 'A':
                vs.add(M(e.c[1])+M(e.r))
                vs.add(M(e.c[1])-M(e.r))
        for e in self.A:
            for f in self.B:
                for pt in meet(e, f):
                    vs.add(pt[1])
        out = []
        for v in sorted(vs):
            if not out or v-out[-1] > M(10)**-30*self.scale:
                out.append(v)
        return out

    def _bands(self):
        """Per band: (v0, v1, [(x_left, key_left, x_right, key_right, atom)])."""
        self.bands = []
        breaks = self._breaks()
        for v0, v1 in zip(breaks, breaks[1:]):
            vm = (v0+v1)/2
            xa = [c for e in self.A for c in e.crossings(vm)]
            xb = [c for e in self.B for c in e.crossings(vm)]
            assert len(xa) % 2 == 0 and len(xb) % 2 == 0, 'odd crossings'
            xs = sorted(xa+xb, key=lambda c: c[0])
            pieces = []
            for (xl, kl), (xr, kr) in zip(xs, xs[1:]):
                if xr == xl:
                    continue
                mid = (xl+xr)/2
                in_a = sum(1 for x, _ in xa if x < mid) % 2 == 1
                in_b = sum(1 for x, _ in xb if x < mid) % 2 == 1
                atom = 'AB' if in_a and in_b else 'A' if in_a else 'B' if in_b else None
                if atom:
                    pieces.append((xl, kl, xr, kr, atom))
            self.bands.append((v0, v1, pieces))
        self.atoms = {a: [M(0), M(0), M(0)] for a in ATOMS}
        for v0, v1, pieces in self.bands:
            for _, kl, _, kr, atom in pieces:
                il, ir = integrals(kl, v0, v1), integrals(kr, v0, v1)
                for i in range(3):
                    self.atoms[atom][i] += ir[i]-il[i]

    def measure(self, atoms):
        """Area and first moments (u, v) of a union of atoms."""
        return tuple(sum((self.atoms[a][i] for a in atoms), M(0)) for i in range(3))

    # -------------------------------------------------------------- boundary

    def _classify(self, pt, tangent, Y):
        """'on_same', 'on_opposite', 'inside' or 'outside' of the profile Y."""
        near = min(Y, key=lambda f: distance(pt, f))
        if distance(pt, near) <= self.eps:
            other = near.tangent(near.param(pt))
            dot = tangent[0]*other[0]+tangent[1]*other[1]
            assert abs(dot) > M(10)**-10, 'a shared piece crossing the other boundary'
            return 'on_same' if dot > 0 else 'on_opposite'
        for attempt in range(8):
            phi = M(0.7)+attempt*M(0.37)
            u = (mp.cos(phi), mp.sin(phi))
            count, ambiguous = 0, False
            for f in Y:
                if f.kind == 'L':
                    d = (M(f.d[0]), M(f.d[1]))
                    den = u[0]*d[1]-u[1]*d[0]
                    if den == 0:
                        continue
                    w = (M(f.p[0])-pt[0], M(f.p[1])-pt[1])
                    s = (w[0]*d[1]-w[1]*d[0])/den
                    t = (w[0]*u[1]-w[1]*u[0])/den
                    if s > 0 and -M(10)**-25 < t < 1+M(10)**-25:
                        if abs(t) < M(10)**-25 or abs(t-1) < M(10)**-25:
                            ambiguous = True
                        elif 0 < t < 1:
                            count += 1
                else:
                    w = (pt[0]-M(f.c[0]), pt[1]-M(f.c[1]))
                    b = u[0]*w[0]+u[1]*w[1]
                    disc = b*b-(w[0]**2+w[1]**2-M(f.r)**2)
                    if disc <= 0:
                        ambiguous = ambiguous or abs(disc) < M(10)**-25*self.scale**2
                        continue
                    for s in (-b+mp.sqrt(disc), -b-mp.sqrt(disc)):
                        if s <= 0:
                            continue
                        if f.full:
                            count += 1
                            continue
                        t = f.param((pt[0]+s*u[0], pt[1]+s*u[1]))
                        if abs(t) < M(10)**-25 or abs(t-1) < M(10)**-25:
                            ambiguous = True
                        elif 0 < t < 1:
                            count += 1
            if not ambiguous:
                return 'inside' if count % 2 else 'outside'
        raise ArithmeticError('no unambiguous ray')

    def pieces(self, X, Y):
        """Each element of X cut at its meetings with Y's elements and Y's
        vertices on it: (element, t0, t1, length, class)."""
        out = []
        for e in X:
            ts = [M(0), M(1)]
            for f in Y:
                for pt in meet(e, f):
                    if within(e, pt, M(10)**-28) and within(f, pt, M(10)**-28):
                        ts.append(e.param(pt))
                v = (M(f.p[0]), M(f.p[1]))
                if distance(v, e) <= self.eps:
                    ts.append(e.param(v))
            ts = sorted(min(max(t, M(0)), M(1)) for t in ts)
            cuts = [ts[0]]
            for t in ts[1:]:
                if t-cuts[-1] > M(10)**-28:
                    cuts.append(t)
            cuts[-1] = M(1)
            for t0, t1 in zip(cuts, cuts[1:]):
                tm = (t0+t1)/2
                cls = self._classify(e.point(tm), e.tangent(tm), Y)
                out.append((e, t0, t1, e.length()*(t1-t0), cls))
        return out

    def _classes(self):
        out = {}
        for name, X, Y in (('A', self.A, self.B), ('B', self.B, self.A)):
            for _, _, _, length, cls in self.pieces(X, Y):
                out[(name, cls)] = out.get((name, cls), M(0))+length
        return out

    def length(self, name, cls):
        return self.classes.get((name, cls), M(0))

    def perimeter(self, atoms):
        """The boundary length of a union of atoms: each class of pieces
        whose two sides the union holds exactly one of."""
        has = lambda a: a in atoms
        total = M(0)
        for name, other in (('A', 'B'), ('B', 'A')):
            if has(name):
                total += self.length(name, 'outside')
            if has('AB') != has(other):
                total += self.length(name, 'inside')
        if has('AB'):
            total += self.length('A', 'on_same')
        if has('A') != has('B'):
            total += self.length('A', 'on_opposite')
        return total

    # -------------------------------------------------------------- solids

    def intervals(self, atoms, band):
        """The region's maximal intervals at a band's middle: (key_left,
        key_right, x_left, x_right)."""
        out = []
        for xl, kl, xr, kr, atom in self.bands[band][2]:
            if atom not in atoms:
                continue
            if out and out[-1][3] == xl:
                out[-1] = (out[-1][0], kr, out[-1][2], xr)
            else:
                out.append((kl, kr, xl, xr))
        return out

    def components(self, regions):
        """The number of maximal connected solid regions of a stack of slab
        regions (a list of atom sets, bottom to top)."""
        parent = {}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            parent[find(a)] = find(b)
        nodes = {}
        for s, atoms in enumerate(regions):
            for k in range(len(self.bands)):
                nodes[(s, k)] = self.intervals(atoms, k)
                for i in range(len(nodes[(s, k)])):
                    parent[(s, k, i)] = (s, k, i)
        overlap = lambda l0, r0, l1, r1: min(r0, r1)-max(l0, l1) > self.eps
        for s in range(len(regions)):
            for k in range(len(self.bands)-1):
                v = self.bands[k][1]
                assert v == self.bands[k+1][0]
                for i, (kl, kr, _, _) in enumerate(nodes[(s, k)]):
                    for j, (ml, mr, _, _) in enumerate(nodes[(s, k+1)]):
                        if overlap(x_at(kl, v), x_at(kr, v), x_at(ml, v), x_at(mr, v)):
                            union((s, k, i), (s, k+1, j))
            if s+1 < len(regions):
                for k in range(len(self.bands)):
                    for i, (_, _, l0, r0) in enumerate(nodes[(s, k)]):
                        for j, (_, _, l1, r1) in enumerate(nodes[(s+1, k)]):
                            if overlap(l0, r0, l1, r1):
                                union((s, k, i), (s+1, k, j))
        return len({find(a) for a in parent})

    # -------------------------------------------------------------- result

    def slabs(self, operation):
        """[(w0, w1, atoms)] bottom to top over every slab between the caps."""
        ws = sorted(set(self.heights['A'] + self.heights['B']))
        out = []
        for w0, w1 in zip(ws, ws[1:]):
            present = [lo <= w0 and w1 <= hi for lo, hi in (self.heights['A'], self.heights['B'])]
            out.append((w0, w1, region(operation, *present)))
        return out

    def xor_area(self, below, above):
        return self.measure(set(below) ^ set(above))[0]

    def result(self, operation):
        """(solids, volume, area, centre in world coordinates, merged runs
        [(w0, w1, area, perimeter)])."""
        slabs = self.slabs(operation)
        V, mu, mv, mw = M(0), M(0), M(0), M(0)
        walls = M(0)
        for w0, w1, atoms in slabs:
            h = M(w1-w0)
            area, su, sv = self.measure(atoms)
            V, mu, mv, mw = V+h*area, mu+h*su, mv+h*sv, mw+h*area*M(w0+w1)/2
            walls += h*self.perimeter(atoms)
        caps = M(0)
        regions = [set()]+[atoms for _, _, atoms in slabs]+[set()]
        for below, above in zip(regions, regions[1:]):
            caps += self.xor_area(below, above)
        tiny = M(10)**-30*self.scale**3
        runs = []
        for w0, w1, atoms in slabs:
            if self.measure(atoms)[0] <= M(10)**-30*self.scale**2:
                continue
            if runs and runs[-1][1] == w0 and self.xor_area(runs[-1][2], atoms) <= M(10)**-30*self.scale**2:
                runs[-1] = (runs[-1][0], w1, runs[-1][2])
            else:
                runs.append((w0, w1, atoms))
        if V <= tiny:
            return 0, M(0), M(0), None, []
        solids = self.components([atoms for _, _, atoms in slabs])
        o, x, y, n = self.axes
        u, v, w = mu/V, mv/V, mw/V
        centre = tuple(M(o[i])+u*M(x[i])+v*M(y[i])+w*M(n[i]) for i in range(3))
        merged = [(w0, w1, self.measure(atoms)[0], self.perimeter(atoms)) for w0, w1, atoms in runs]
        return solids, V, walls+caps, centre, merged


def number(x):
    x = M(x)
    return '0.0' if abs(x) < M(10)**-30 else mp.nstr(x, 25, min_fixed=-5, max_fixed=5)


def rows(obj, operation, tool):
    """The expected rows of a Boolean: `result N volume area cx cy cz` (the
    totals over the result's N solids and their common centre, world
    coordinates) or `empty`, then per merged run of slabs `slab w0 w1 area
    perimeter` (the heights in the object's frame, the run's profile's area
    and boundary length)."""
    pair = Pair(obj, tool)
    solids, V, A, centre, runs = pair.result(operation)
    if solids == 0:
        return ['empty'], pair
    out = [' '.join(['result', str(solids), number(V), number(A)]+[number(c) for c in centre])]
    for w0, w1, area, perimeter in runs:
        out.append(' '.join(['slab', number(w0), number(w1), number(area), number(perimeter)]))
    return out, pair


def step(pair, operation):
    """S9a.1 when the result is empty or one run of slabs (each solid one
    prism), S9a.2 otherwise (a stack)."""
    runs = pair.result(operation)[4]
    return 'S9a.1' if len(runs) <= 1 else 'S9a.2'
