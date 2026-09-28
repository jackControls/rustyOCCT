#!/usr/bin/env python3
"""Independent reference for S9b of REVIEW_NOTES.md: Booleans of two
polyhedral prisms (line profiles, polygon holes) in any relative position.

Each prism is its construction's exact model: a point is `o + u x + v y + w
n` in Fractions from the frame's stored binary64 origin and axes
(`stored_axes`, the kernel's `Frame3::new`; not exactly orthonormal) and the
profile's and offsets' binary64 values. Nothing here uses the kernel, a
plane arrangement or a 2D Boolean:

* **Convex cells.** Each profile (outer polygon and polygon holes) is cut
  into trapezoids by vertical lines through its vertices (in its `(u, v)`:
  within a strip the edges spanning it are ordered at its middle and paired
  by parity), each trapezoid swept over the offsets a convex polyhedron, so
  each prism is a union of convex cells with disjoint interiors.
* **Clipping.** A convex polyhedron (faces as exact point cycles, outward)
  is clipped by a half-space by clipping each face (Sutherland-Hodgman) and
  closing the cut with the points on the plane, ordered by exact angles
  about their centroid in the plane's projection. The common is the union of
  every pair of cells' intersections (disjoint interiors), so its volume and
  first moments are exact sums of signed tetrahedra; the fuse and the cut
  follow by inclusion and exclusion (`fuse = A + B - common`, `cut = A -
  common`).
* **Surface area.** Every boundary face of each prism (its caps cut into
  trapezoids, its walls) is split by every plane of the other prism's cells
  into cells no such plane crosses; each face cell's centroid, pushed an
  infinitesimal `e` along the face's outward normal and against it,
  classifies exactly against the other prism (inside a cell when on the
  inner side of all its planes, a zero decided by the push's sign). A cell
  of the object with the tool's membership `t+` in front and `t-` behind lies
  on the result's boundary when the result holds `(0, t+)` and not `(1, t-)`
  or the reverse; a cell of the tool counts likewise unless the object's
  membership is the same on both sides' opposite: only when the object's
  front and back memberships are equal (a cell on both boundaries is the
  object's). Areas are exact squared norms, their roots in mpmath (40
  digits).
* **Solids.** The result's convex cells (the common's pair intersections;
  a cut's each object cell less every tool cell, as successive half-space
  complements; a fuse's cut cells and the tool's cells) are joined when two
  of their faces lie on one plane, face each other and overlap in a positive
  area (exact 2D clipping): regions meeting at an edge or a point are
  separate solids (the regularized Boolean). A result of zero volume is
  empty.
"""
from fractions import Fraction as F
from functools import cmp_to_key

import mpmath as mp

from identity_reference import stored
from curve_surface_reference import stored_axes

mp.mp.dps = 40


# ------------------------------------------------------------------ vectors

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def scale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def det3(a, b, c):
    return dot(a, cross(b, c))


def M(v):
    if isinstance(v, F):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def normal_of(face):
    """Twice the vector area of a planar cycle (its outward normal)."""
    total = (F(0), F(0), F(0))
    r = face[0]
    for a, b in zip(face[1:], face[2:]):
        total = add(total, cross(sub(a, r), sub(b, r)))
    return total


def area(face):
    n = normal_of(face)
    return mp.sqrt(M(dot(n, n)))/2


# ------------------------------------------------------------------ prisms

class Prism:
    """A prism's exact model: its frame's stored axes, its profile's
    polygons (outer counter-clockwise, holes as stored) and its heights."""

    def __init__(self, case):
        self.axes = tuple(tuple(F(c) for c in v) for v in stored_axes(case.frame))
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        self.polygons = []
        for b in case.boundaries:
            pts, _ = stored(b, case.tolerance)
            assert pts is not None, 'S9b: polygon boundaries'
            if b.segments is not None:
                points, segs = pts
                assert all(s is None for s in segs), 'S9b: line segments only'
                pts = points
            self.polygons.append([(F(p[0]), F(p[1])) for p in pts])
        self.cells = [self.cell(t) for t in trapezoids(self.polygons)]
        self.faces = self.boundary_faces()

    def point(self, u, v, w):
        o, x, y, n = self.axes
        return tuple(o[i]+u*x[i]+v*y[i]+w*n[i] for i in range(3))

    def orientation(self):
        """+1 when (x, y, n) is right-handed (the stored axes are)."""
        _, x, y, n = self.axes
        return 1 if det3(x, y, n) > 0 else -1

    def cell(self, poly):
        """A convex polygon of the profile (counter-clockwise) swept: faces
        as outward cycles."""
        assert self.orientation() > 0
        bottom = [self.point(u, v, self.lo) for u, v in reversed(poly)]
        top = [self.point(u, v, self.hi) for u, v in poly]
        faces = [bottom, top]
        k = len(poly)
        for i in range(k):
            (u0, v0), (u1, v1) = poly[i], poly[(i+1) % k]
            faces.append([self.point(u0, v0, self.lo), self.point(u1, v1, self.lo),
                          self.point(u1, v1, self.hi), self.point(u0, v0, self.hi)])
        return faces

    def boundary_faces(self):
        """The prism's boundary as convex outward cycles: its caps' trapezoids
        and each wall (every boundary edge, holes' walls facing into the
        hole)."""
        faces = []
        for t in trapezoids(self.polygons):
            faces.append([self.point(u, v, self.lo) for u, v in reversed(t)])
            faces.append([self.point(u, v, self.hi) for u, v in t])
        for k, poly in enumerate(self.polygons):
            ring = poly if signed_area2(poly) > 0 else list(reversed(poly))
            if k > 0:
                ring = list(reversed(ring))
            m = len(ring)
            for i in range(m):
                (u0, v0), (u1, v1) = ring[i], ring[(i+1) % m]
                faces.append([self.point(u0, v0, self.lo), self.point(u1, v1, self.lo),
                              self.point(u1, v1, self.hi), self.point(u0, v0, self.hi)])
        return faces


def signed_area2(poly):
    return sum(p[0]*q[1]-q[0]*p[1] for p, q in zip(poly, poly[1:]+poly[:1]))


def trapezoids(polygons):
    """The region (outer polygon less holes) as convex counter-clockwise
    polygons: strips between consecutive vertex `u`s, the spanning edges
    ordered at the strip's middle and paired by parity."""
    edges = []
    for poly in polygons:
        m = len(poly)
        for i in range(m):
            p, q = poly[i], poly[(i+1) % m]
            if p[0] != q[0]:
                edges.append((p, q) if p[0] < q[0] else (q, p))
    us = sorted({p[0] for poly in polygons for p in poly})
    out = []
    for u0, u1 in zip(us, us[1:]):
        mid = (u0+u1)/2
        at = lambda e, u: e[0][1]+(e[1][1]-e[0][1])*(u-e[0][0])/(e[1][0]-e[0][0])
        span = [e for e in edges if e[0][0] <= u0 and e[1][0] >= u1]
        span.sort(key=lambda e: at(e, mid))
        assert len(span) % 2 == 0, 'an odd strip'
        for lo, hi in zip(span[0::2], span[1::2]):
            pts = [(u0, at(lo, u0)), (u1, at(lo, u1)), (u1, at(hi, u1)), (u0, at(hi, u0))]
            ring = []
            for p in pts:
                if not ring or ring[-1] != p:
                    ring.append(p)
            if len(ring) > 1 and ring[0] == ring[-1]:
                ring.pop()
            if len(ring) >= 3:
                out.append(ring)
    return out


# ------------------------------------------------------------------ clipping

def planes(cell):
    """Each face's plane `(n, a)`: inside where `n . (X - a) <= 0`."""
    return [(normal_of(f), f[0]) for f in cell if len(f) >= 3]


def clip_polygon(poly, n, a, keep):
    """The part of a planar cycle where `n . (X - a)` has `keep`'s sign or
    vanishes (`keep = -1`: `<= 0`; `+1`: `>= 0`), and its points on the
    plane."""
    out, on = [], []
    s = [-dot(n, sub(p, a))*keep for p in poly]
    m = len(poly)
    for i in range(m):
        p, q, sp, sq = poly[i], poly[(i+1) % m], s[i], s[(i+1) % m]
        if sp <= 0:
            out.append(p)
            if sp == 0:
                on.append(p)
        if (sp < 0 < sq) or (sq < 0 < sp):
            t = sp/(sp-sq)
            x = add(p, scale(sub(q, p), t))
            out.append(x)
            on.append(x)
    return out, on


def order_cycle(points, n):
    """Distinct coplanar points of a convex polygon ordered counter-clockwise
    about `n`."""
    pts = []
    for p in points:
        if p not in pts:
            pts.append(p)
    if len(pts) < 3:
        return []
    k = max(range(3), key=lambda i: abs(n[i]))
    i, j = [(1, 2), (2, 0), (0, 1)][k]
    sign = 1 if n[k] > 0 else -1
    cu = sum(p[i] for p in pts)/len(pts)
    cv = sum(p[j] for p in pts)/len(pts)

    def half(p):
        du, dv = p[i]-cu, p[j]-cv
        return 0 if (dv > 0 or (dv == 0 and du > 0)) else 1

    def cmp(p, q):
        hp, hq = half(p), half(q)
        if hp != hq:
            return hp-hq
        c = (p[i]-cu)*(q[j]-cv)-(p[j]-cv)*(q[i]-cu)
        return -1 if c > 0 else 1 if c < 0 else 0

    ordered = sorted(pts, key=cmp_to_key(cmp))
    if sign < 0:
        ordered.reverse()
    # Drop collinear points (the cycle's corners only).
    out = []
    for idx, p in enumerate(ordered):
        a, b = ordered[idx-1], ordered[(idx+1) % len(ordered)]
        if cross(sub(p, a), sub(b, p)) != (0, 0, 0):
            out.append(p)
    return out if len(out) >= 3 else []


def clip_cell(cell, n, a, keep=-1):
    """A convex polyhedron clipped by a closed half-space."""
    faces, cap, flat = [], [], False
    for f in cell:
        g, on = clip_polygon(f, n, a, keep)
        cap += on
        g = order_cycle(g, normal_of(f)) if len(g) >= 3 else []
        if g:
            faces.append(g)
            # A kept face in the plane already closes the cut.
            flat = flat or all(dot(n, sub(p, a)) == 0 for p in g)
    cyc = [] if flat else order_cycle(cap, scale(n, -keep))
    if cyc:
        faces.append(cyc)
    if not faces or volume(faces)[0] == 0:
        return []
    return faces


def intersect(p, q):
    out = p
    for n, a in planes(q):
        out = clip_cell(out, n, a)
        if not out:
            return []
    return out


def difference(p, q):
    """`p` less `q` as convex cells: the successive complements of `q`'s
    half-spaces."""
    pieces, rest = [], p
    for n, a in planes(q):
        outside = clip_cell(rest, n, a, keep=1)
        if outside:
            pieces.append(outside)
        rest = clip_cell(rest, n, a)
        if not rest:
            break
    return pieces


def volume(cell):
    """Exact volume and first moments."""
    r = cell[0][0]
    vol, mom = F(0), (F(0), F(0), F(0))
    for f in cell:
        for a, b in zip(f[1:], f[2:]):
            v = det3(sub(f[0], r), sub(a, r), sub(b, r))/6
            vol += v
            c = scale(add(add(r, f[0]), add(a, b)), F(1, 4))
            mom = add(mom, scale(c, v))
    return vol, mom


def inside(cells, x, d):
    """Whether `x + e d` lies inside a union of cells for every small `e > 0`."""
    for cell in cells:
        ok = True
        for n, a in planes(cell):
            s = dot(n, sub(x, a))
            if s == 0:
                s = dot(n, d)
                assert s != 0, 'a push along a plane'
            if s > 0:
                ok = False
                break
        if ok:
            return True
    return False


def split_face(face, cells):
    """A convex face cycle cut by every plane of the cells that crosses it."""
    parts = [face]
    for cell in cells:
        for n, a in planes(cell):
            nxt = []
            for part in parts:
                s = [dot(n, sub(p, a)) for p in part]
                if min(s) < 0 < max(s):
                    fn = normal_of(part)
                    for keep in (-1, 1):
                        g, _ = clip_polygon(part, n, a, keep)
                        g = order_cycle(g, fn)
                        if g:
                            nxt.append(g)
                else:
                    nxt.append(part)
            parts = nxt
    return parts


def centroid(poly):
    return scale(tuple(sum(p[i] for p in poly) for i in range(3)), F(1, len(poly)))


# ------------------------------------------------------------------ the Boolean

OPS = {
    'fuse': lambda a, b: a or b,
    'cut': lambda a, b: a and not b,
    'common': lambda a, b: a and b,
}


class Pair:
    def __init__(self, obj, tool):
        self.A, self.B = Prism(obj), Prism(tool)
        self.common_cells = [c for p in self.A.cells for q in self.B.cells for c in [intersect(p, q)] if c]

    def measure(self, cells):
        vol, mom = F(0), (F(0), F(0), F(0))
        for c in cells:
            v, m = volume(c)
            vol += v
            mom = add(mom, m)
        return vol, mom

    def boundary_area(self, operation):
        f = OPS[operation]
        total = mp.mpf(0)
        for x, y, first in ((self.A, self.B, True), (self.B, self.A, False)):
            for face in x.faces:
                n = normal_of(face)
                for part in split_face(face, y.cells):
                    c = centroid(part)
                    front, back = inside(y.cells, c, n), inside(y.cells, c, scale(n, -1))
                    if first:
                        counts = f(False, front) != f(True, back)
                    else:
                        counts = front == back and f(front, False) != f(back, True)
                    if counts:
                        total += area(part)
        return total

    def cells(self, operation):
        if operation == 'common':
            return list(self.common_cells)
        cut = []
        for p in self.A.cells:
            pieces = [p]
            for q in self.B.cells:
                pieces = [d for piece in pieces for d in difference(piece, q)]
            cut += pieces
        if operation == 'cut':
            return cut
        return cut+list(self.B.cells)

    def result(self, operation):
        """(solids, volume, area, centre) in world coordinates, exact volume
        and centre as Fractions; area in mpmath."""
        va, ma = self.measure(self.A.cells)
        vb, mb = self.measure(self.B.cells)
        vc, mc = self.measure(self.common_cells)
        if operation == 'common':
            vol, mom = vc, mc
        elif operation == 'cut':
            vol, mom = va-vc, sub(ma, mc)
        else:
            vol, mom = va+vb-vc, sub(add(ma, mb), mc)
        if vol == 0:
            return 0, F(0), mp.mpf(0), None
        cells = self.cells(operation)
        # The cells' measures must add up to the inclusion-exclusion's.
        cv, cm = self.measure(cells)
        assert (cv, cm) == (vol, mom), 'the cells miss the inclusion-exclusion measures'
        return solids(cells), vol, self.boundary_area(operation), scale(mom, 1/vol)


def solids(cells):
    parent = list(range(len(cells)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(cells)):
        for j in range(i+1, len(cells)):
            if find(i) != find(j) and touching(cells[i], cells[j]):
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(cells))})


def touching(c, d):
    """Two faces on one plane, facing each other, overlapping in area."""
    for f in c:
        nf = normal_of(f)
        for g in d:
            ng = normal_of(g)
            if cross(nf, ng) != (0, 0, 0) or dot(nf, ng) >= 0:
                continue
            if dot(nf, sub(g[0], f[0])) != 0:
                continue
            overlap = f
            for a, b in zip(g, g[1:]+g[:1]):
                # g's edges, g counter-clockwise about ng: inside on the left
                # seen along ng, the right seen along nf.
                edge_normal = cross(sub(b, a), ng)
                overlap, _ = clip_polygon(overlap, edge_normal, a, -1)
                overlap = order_cycle(overlap, nf) if len(overlap) >= 3 else []
                if not overlap:
                    break
            if overlap and dot(normal_of(overlap), normal_of(overlap)) > 0:
                return True
    return False


def number(x):
    x = M(x)
    return '0.0' if abs(x) < M(10)**-30 else mp.nstr(x, 25, min_fixed=-5, max_fixed=5)


def rows(obj, operation, tool):
    """`result N volume area cx cy cz` (totals over the N solids, world
    coordinates) or `empty`."""
    pair = Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair
