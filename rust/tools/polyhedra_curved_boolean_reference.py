#!/usr/bin/env python3
"""Independent reference for S9e.4b.4c.1 of REVIEW_NOTES.md: imported
polyhedra (the constructions OCCT was given for them) against curved faces,
their results given to further Booleans, and polyhedra with a cavity.

It is S9e.3a's chained reference (`chained_curved_boolean_reference.py`:
every face of every input swept by two families of curves, each curve cut at
the roots of the other inputs' surfaces, each piece classified by the
chain's set function and measured by the divergence theorem) with one more
kind of input, a convex hull of exact points (`Hull`): a wedge's corners on
its frame's stored axes, a polyhedron's binary64 points
(`imported_polyhedra_boolean_reference.hull`: each plane through three of
the points with every point on its inner side, its points ordered about its
outward normal). A body that is not convex is a Boolean of hulls and boxes
(the construction OCCT computed, or a union of its convex cells), one
expression of the chain. Nothing here uses the kernel, a first result or an
arrangement.

* **The hull's surfaces** are its faces' planes `n . X = d` (`n` twice the
  face's vector area, exact rationals, the function `n . X - d` in mpf), its
  membership every plane's strict inner side.
* **Its faces** are swept by chords of two directions in each face's plane
  (a unit direction across the face and one turned from it by `atan(7/2)`,
  as the prisms' caps: `(1, 0)` and `(2, 7)` in an orthonormal basis of the
  plane made from the plane's normal in mpf), the outer parameter the offset
  across the chord's direction, each chord's range where its line lies in
  the convex polygon (its corners' offsets the family's kinks).
* **Its closed form**: the exact volume and first moments of its cells' face
  cycles (`polyhedral_reference.volume`), the area each face's twice vector
  area's norm halved.
* **Solids by rays**: the binary64 hull, each ray's interval where every
  plane's value is negative (`Float`, the chained reference's binary64
  inputs with the hull added, which `count_solids` and `monte_carlo` use).
"""
from fractions import Fraction as F

import mpmath as mp

import chained_curved_boolean_reference as ch
import imported_polyhedra_boolean_reference as ipr
import polyhedral_reference as pr
import torus_curved_boolean_reference as tc
from torus_curved_boolean_reference import Z, Face, Family, Surf, line_curve, lin, norm
from sphere_boolean_reference import M, cross, dot, scale

mp.mp.dps = 40


class HullCase:
    """A convex hull's construction: its exact points (each a binary64 point
    or a rational one, the construction's own)."""

    def __init__(self, points):
        self.points = [tuple(F(c) for c in p) for p in points]
        self.hull = True


def unit(v):
    n = norm(v)
    return tuple(x/n for x in v)


class Hull:
    """The convex hull of exact points."""
    kind = 'hull'

    def __init__(self, case):
        self.case = case
        self.cycles = ipr.hull(case.points)
        self.points = []
        for c in self.cycles:
            for p in c:
                if p not in self.points:
                    self.points.append(p)
        # Each face's plane: its outward normal (twice its vector area) and
        # offset, exactly.
        self.planes = []
        for c in self.cycles:
            n = pr.normal_of(c)
            self.planes.append((n, dot(n, c[0])))
        self.surfs = []
        for k, (n, d) in enumerate(self.planes):
            nm, dm = tuple(M(x) for x in n), M(d)
            self.surfs.append(Surf(('face', k), 1, lambda X, nm=nm, dm=dm: dot(nm, X)-dm,
                                   lambda X, nm=nm: nm))
        self.faces = [Face(('face', k), [self.chords(k, 0), self.chords(k, 1)])
                      for k in range(len(self.cycles))]

    def contains(self, X):
        return all(s.f(X) < 0 for s in self.surfs)

    def chords(self, k, which):
        """Chords of face `k` along its first (`which` 0) or second
        direction, the outer parameter the offset across it."""
        cyc = [tuple(M(x) for x in p) for p in self.cycles[k]]
        n = unit(tuple(M(x) for x in self.planes[k][0]))
        # An axis least along the normal, its part in the plane the first
        # direction.
        a = min(range(3), key=lambda i: abs(n[i]))
        e = [Z, Z, Z]
        e[a] = mp.mpf(1)
        u = unit(cross(n, tuple(e)))
        w = cross(n, u)
        if which == 1:
            u, w = unit(lin((2, u), (7, w))), None
            w = cross(n, u)
        # In the face's plane: `X = c0 + s u + t w`; the polygon in `(s, t)`.
        c0 = cyc[0]
        st = [(dot(tuple(p[i]-c0[i] for i in range(3)), u), dot(tuple(p[i]-c0[i] for i in range(3)), w))
              for p in cyc]
        ts = [q[1] for q in st]
        lo, hi = min(ts), max(ts)
        zero = (Z, Z, Z)

        def chord(t):
            ss = []
            m = len(st)
            for i in range(m):
                (s0, t0), (s1, t1) = st[i], st[(i+1) % m]
                if t0 == t1:
                    continue
                f = (t-t0)/(t1-t0)
                if 0 <= f <= 1:
                    ss.append(s0+f*(s1-s0))
            b0, b1 = (min(ss), max(ss)) if len(ss) >= 2 else (Z, Z)
            return line_curve(lin((1, c0), (t, w)), u, w, zero, b0, b1)
        # `X_t x X_b = w x u = -n`: outward needs the sign -1.
        own = [t for t in ts if lo < t < hi]
        return Family(('chords', which), lo, hi, False, chord, -1, own)

    def closed(self):
        """Volume, world moments and area, exactly but the area's roots."""
        V, mom = pr.volume(self.cycles)
        area = sum((pr.area(c) for c in self.cycles), Z)
        return M(V), tuple(M(x) for x in mom), area

    def box(self):
        lo = [min(float(p[i]) for p in self.points) for i in range(3)]
        hi = [max(float(p[i]) for p in self.points) for i in range(3)]
        return lo, hi


def make_input(case):
    if getattr(case, 'hull', None) is True:
        return Hull(case)
    return tc.make_input(case)


class Chain(ch.Chain):
    """The chained reference's chain, its inputs hulls too."""

    def __init__(self, cases):
        self.inputs = [make_input(c) for c in cases]
        boxes = [I.box() for I in self.inputs]
        self.size = mp.mpf(tc.box_size(boxes))
        lo = [min(b[0][i] for b in boxes) for i in range(3)]
        hi = [max(b[1][i] for b in boxes) for i in range(3)]
        self.c0 = tuple(mp.mpf((lo[i]+hi[i])/2) for i in range(3))
        self._sweeps = {}


class Float(ch.Float):
    """The chained reference's binary64 input, a hull's planes too."""

    def __init__(self, S):
        if S.kind != 'hull':
            super().__init__(S)
            return
        self.S = S
        self.kind = 'hull'
        self.planes = [([float(x) for x in n], float(d)) for n, d in S.planes]

    def polys(self, P, D):
        if self.kind != 'hull':
            return super().polys(P, D)
        return [[sum(n[i]*P[i] for i in range(3))-d, sum(n[i]*D[i] for i in range(3))] for n, d in self.planes]

    def contains(self, X):
        if self.kind != 'hull':
            return super().contains(X)
        return all(sum(n[i]*X[i] for i in range(3))-d < 0 for n, d in self.planes)


# The chained reference's rays and Monte Carlo take hulls too.
ch.Float = Float

rows = ch.rows
monte_carlo = ch.monte_carlo
count_solids = ch.count_solids
leaves = ch.leaves
evaluate = ch.evaluate
chain_expr = ch.chain_expr

__all__ = ['Chain', 'Float', 'Hull', 'HullCase', 'chain_expr', 'count_solids', 'evaluate', 'leaves',
           'make_input', 'monte_carlo', 'rows', 'scale']
