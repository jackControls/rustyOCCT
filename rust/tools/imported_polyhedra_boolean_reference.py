#!/usr/bin/env python3
"""Independent reference for S9e.4b.2 of REVIEW_NOTES.md: Booleans of
polyhedra in any position, each the construction OCCT was given for an
imported body (a wedge, a polyhedron of given points, a Boolean of boxes)
or a prism of lines the kernel builds from rows.

Every solid is a union of convex cells with disjoint interiors in exact
Fractions: a prism's are S9b's (`polyhedral_reference.Prism`, its frame's
stored axes and its profile's trapezoids); a convex body's one cell is the
convex hull of its exact points (`hull`: each plane through three of them
with every point on one side, its points on it ordered about its normal);
a Boolean's result's cells are S9b's construction (the common every pair
of cells' intersection, the cut each object cell less every tool cell as
successive half-space complements, the fuse the cut's cells and the tool's
cells). Nothing here uses the kernel. Per Boolean:

* **Volumes three ways.** By inclusion and exclusion of the inputs' and
  the common's cells (`fuse = A + B - common`, `cut = A - common`), by the
  result's own cells, and by the divergence theorem over the result's
  boundary faces (signed tetrahedra from a point): all exact, all equal.
  First moments likewise (the first two ways).
* **Areas two ways.** The inputs' boundary faces split by every plane of
  the other's cells and classified by an infinitesimal push either way
  (`polyhedral_reference.Pair.boundary_area`'s rule: a part of the object
  counts where the operation's set function differs across it, a part of
  the tool only where the object's membership is the same on both sides),
  and the result's own boundary (`Body.faces`: every cell's face split by
  the other cells' planes, a part kept where the union holds the side
  behind it and not the side in front). Areas are exact squared norms,
  their roots in mpmath (40 digits).
* **Solids.** The result's cells joined where two of their faces lie on
  one plane, face each other and overlap in a positive area
  (`polyhedral_reference.solids`): bodies meeting at an edge or a point
  are separate solids. A result of zero volume is empty.
* **Margins.** `clearance`: each input's boundary vertices' least
  distance from the other's boundary faces where it is not zero (a vertex
  exactly on a face is an exact contact, counted), relative to the size;
  `sine`: the least sine between two faces of different inputs whose
  planes each cross the other face.
"""
from fractions import Fraction as F
from itertools import combinations

import mpmath as mp

import polyhedral_reference as pr
from polyhedral_reference import add, cross, dot, normal_of, scale, sub, volume
from curve_surface_reference import stored_axes

mp.mp.dps = 40
DEBUG = False

OPS = pr.OPS


# ------------------------------------------------------------------ cells

def hull(points):
    """The convex hull of exact points as outward face cycles (each plane
    through three points with every point on its inner side or on it)."""
    pts = []
    for p in points:
        p = tuple(F(c) for c in p)
        if p not in pts:
            pts.append(p)
    faces, seen = [], set()
    for a, b, c in combinations(pts, 3):
        n = cross(sub(b, a), sub(c, a))
        if n == (0, 0, 0):
            continue
        s = [dot(n, sub(p, a)) for p in pts]
        if all(x <= 0 for x in s):
            pass
        elif all(x >= 0 for x in s):
            n = scale(n, -1)
        else:
            continue
        # The plane's canonical key: its normal scaled to a first
        # coordinate of +-1 (oriented), and its offset.
        k = next(x for x in n if x != 0)
        key = (tuple(x/abs(k) for x in n), dot(n, a)/abs(k))
        if key in seen:
            continue
        seen.add(key)
        on = [p for p in pts if dot(n, sub(p, a)) == 0]
        cyc = pr.order_cycle(on, n)
        assert cyc, 'a hull face'
        faces.append(cyc)
    assert faces and volume(faces)[0] > 0, 'a hull of positive volume'
    return faces


def measure(cells):
    vol, mom = F(0), (F(0), F(0), F(0))
    for c in cells:
        v, m = volume(c)
        vol += v
        mom = add(mom, m)
    return vol, mom


def boundary(cells):
    """The union's boundary as convex outward cycles: each cell's face cut
    by every other cell's planes, a part kept where the union holds the side
    behind it and not the side in front."""
    out = []
    for i, cell in enumerate(cells):
        others = cells[:i]+cells[i+1:]
        for face in cell:
            n = normal_of(face)
            for part in pr.split_face(face, others):
                c = pr.centroid(part)
                if pr.inside(cells, c, scale(n, -1)) and not pr.inside(cells, c, n):
                    out.append(part)
    return out


class Body:
    """A solid as convex cells with disjoint interiors, and its boundary."""

    def __init__(self, cells):
        self.cells = [c for c in cells if c]
        self.faces = boundary(self.cells)

    def measure(self):
        return measure(self.cells)

    def points(self):
        out = []
        for f in self.faces:
            for p in f:
                if p not in out:
                    out.append(p)
        return out

    def vertices(self):
        """The boundary's corners: points of its faces' cycles on three or
        more distinct planes of them (not the points where a face is split
        within itself or along an edge)."""
        planes = {}
        for f in self.faces:
            n = normal_of(f)
            k = next(x for x in n if x != 0)
            key = (tuple(x/abs(k) for x in n), dot(n, f[0])/abs(k))
            for p in f:
                planes.setdefault(p, set()).add(key)
        return [p for p, ps in planes.items() if len(ps) >= 3]

    def edges(self):
        """The boundary's edges as segments: the sides of its faces' cycles
        on a plane of another face too (pieces of an edge where a face was
        split along it), each once."""
        keys = []
        for f in self.faces:
            n = normal_of(f)
            k = next(x for x in n if x != 0)
            keys.append((tuple(x/abs(k) for x in n), dot(n, f[0])/abs(k)))
        planes = sorted(set(keys))
        out = set()
        for f, own in zip(self.faces, keys):
            for p, q in zip(f, f[1:]+f[:1]):
                if any(key != own and dot(key[0], p) == key[1] and dot(key[0], q) == key[1] for key in planes):
                    out.add((p, q) if p <= q else (q, p))
        return sorted(out)

    def size(self):
        return max([F(1)]+[abs(x) for p in self.points() for x in p])


def prism_body(case):
    return Body(pr.Prism(case).cells)


def wedge_points(frame, dx, dy, dz, xmin, zmin, xmax, zmax):
    """`BRepPrimAPI_MakeWedge(gp_Ax2, DX, DY, DZ, XMIN, ZMIN, XMAX, ZMAX)`'s
    corners on the frame's stored axes (`o + u x + v y + w n`, its face at
    `v = 0` the rectangle `[0, DX] x [0, DZ]` in `(u, w)`, its face at `v =
    DY` the rectangle `[XMIN, XMAX] x [ZMIN, ZMAX]`)."""
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    at = lambda u, v, w: tuple(o[i]+F(u)*x[i]+F(v)*y[i]+F(w)*n[i] for i in range(3))
    base = [at(u, 0, w) for u in (0, dx) for w in (0, dz)]
    top = [at(u, dy, w) for u in (xmin, xmax) for w in (zmin, zmax)]
    return base+top


def wedge_closed(frame, dx, dy, dz, xmin, zmin, xmax, zmax):
    """The wedge's volume by the prismatoid formula (`h (A0 + 4 Am + A1) /
    6` of its rectangles' sections) times the stored axes' determinant."""
    _, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    a0, a1 = F(dx)*F(dz), (F(xmax)-F(xmin))*(F(zmax)-F(zmin))
    am = (F(dx)+F(xmax)-F(xmin))/2*(F(dz)+F(zmax)-F(zmin))/2
    return F(dy)*(a0+4*am+a1)/6*pr.det3(x, y, n)


def prismatoid(bottom, top, height_axis):
    """The volume of a polyhedron whose vertices lie in two parallel planes
    normal to a coordinate axis (`height_axis`), its bottom and top
    polygons' vertices corresponding (the middle section their mean):
    `h (A0 + 4 Am + A1) / 6` with the polygons' areas projected on the
    other two coordinates."""
    i, j = [(1, 2), (2, 0), (0, 1)][height_axis]
    area = lambda poly: abs(pr.signed_area2([(p[i], p[j]) for p in poly]))/2
    mid = [tuple((F(a[k])+F(b[k]))/2 for k in range(3)) for a, b in zip(bottom, top)]
    h = abs(F(top[0][height_axis])-F(bottom[0][height_axis]))
    return h*(area(bottom)+4*area(mid)+area(top))/6


# ------------------------------------------------------------------ the Boolean

def boolean(x, y, operation):
    """A Boolean's result as a body: S9b's cells."""
    common = [c for p in x.cells for q in y.cells for c in [pr.intersect(p, q)] if c]
    if operation == 'common':
        return Body(common), common
    cut = []
    for p in x.cells:
        pieces = [p]
        for q in y.cells:
            pieces = [d for piece in pieces for d in pr.difference(piece, q)]
        cut += pieces
    if operation == 'cut':
        return Body(cut), common
    return Body(cut+list(y.cells)), common


def classified_area(x, y, operation):
    """The result's area from the inputs' faces classified against each
    other (`polyhedral_reference.Pair.boundary_area`'s rule)."""
    f = OPS[operation]
    total = mp.mpf(0)
    for a, b, first in ((x, y, True), (y, x, False)):
        for face in a.faces:
            n = normal_of(face)
            for part in pr.split_face(face, b.cells):
                c = pr.centroid(part)
                front, back = pr.inside(b.cells, c, n), pr.inside(b.cells, c, scale(n, -1))
                if first:
                    counts = f(False, front) != f(True, back)
                else:
                    counts = front == back and f(front, False) != f(back, True)
                if counts:
                    total += pr.area(part)
    return total


def divergence_volume(faces):
    """A closed outward surface's volume by signed tetrahedra from the
    first face's first point."""
    r = faces[0][0]
    vol = F(0)
    for f in faces:
        for a, b in zip(f[1:], f[2:]):
            vol += pr.det3(sub(f[0], r), sub(a, r), sub(b, r))/6
    return vol


def result(x, y, operation):
    """`(solids, volume, area, centre, body, checks)`: the result's solid
    count, exact volume and centre (Fractions), area (mpmath), body, and the
    two-way and three-way deviations (relative to the size)."""
    body, common = boolean(x, y, operation)
    va, ma = x.measure()
    vb, mb = y.measure()
    vc, mc = measure(common)
    if operation == 'common':
        vol, mom = vc, mc
    elif operation == 'cut':
        vol, mom = va-vc, sub(ma, mc)
    else:
        vol, mom = va+vb-vc, sub(add(ma, mb), mc)
    cv, cm = body.measure()
    assert (cv, cm) == (vol, mom), 'the cells miss the inclusion-exclusion measures'
    size = max(x.size(), y.size())
    checks = {}
    if vol == 0:
        assert not body.faces, 'an empty result with a boundary'
        return 0, F(0), mp.mpf(0), None, body, checks
    checks['volume_three_ways'] = abs(divergence_volume(body.faces)-vol)/size**3
    a1 = sum((pr.area(f) for f in body.faces), mp.mpf(0))
    a2 = classified_area(x, y, operation)
    checks['area_two_ways'] = abs(a1-a2)/size**2
    return pr.solids(body.cells), vol, a2, scale(mom, 1/vol), body, checks


# ------------------------------------------------------------------ margins

def point_face_distance(p, face):
    """A point's distance from a closed convex planar cycle (mpmath)."""
    P = tuple(pr.M(c) for c in p)
    Q = [tuple(pr.M(c) for c in q) for q in face]
    n = tuple(pr.M(c) for c in normal_of(face))
    nn = mp.sqrt(sum(c*c for c in n))
    n = tuple(c/nn for c in n)
    h = sum(n[k]*(P[k]-Q[0][k]) for k in range(3))
    foot = tuple(P[k]-h*n[k] for k in range(3))
    m = len(Q)
    inside = True
    for i in range(m):
        a, b = Q[i], Q[(i+1) % m]
        e = tuple(b[k]-a[k] for k in range(3))
        w = tuple(foot[k]-a[k] for k in range(3))
        c = (e[1]*w[2]-e[2]*w[1], e[2]*w[0]-e[0]*w[2], e[0]*w[1]-e[1]*w[0])
        if sum(c[k]*n[k] for k in range(3)) < 0:
            inside = False
            break
    if inside:
        return abs(h)
    best = None
    for i in range(m):
        a, b = Q[i], Q[(i+1) % m]
        e = tuple(b[k]-a[k] for k in range(3))
        w = tuple(P[k]-a[k] for k in range(3))
        t = sum(e[k]*w[k] for k in range(3))/sum(e[k]*e[k] for k in range(3))
        t = min(max(t, mp.mpf(0)), mp.mpf(1))
        d = mp.sqrt(sum((P[k]-a[k]-t*e[k])**2 for k in range(3)))
        best = d if best is None else min(best, d)
    return best


def on_face(p, face):
    """Whether an exact point lies on a closed convex planar cycle, exactly."""
    n = normal_of(face)
    if dot(n, sub(p, face[0])) != 0:
        return False
    m = len(face)
    return all(dot(cross(sub(face[(i+1) % m], face[i]), sub(p, face[i])), n) >= 0 for i in range(m))


def segment_distance(e, f):
    """The least distance between two segments (mpmath)."""
    P = [tuple(pr.M(c) for c in p) for p in e]
    Q = [tuple(pr.M(c) for c in p) for p in f]
    d1 = tuple(P[1][k]-P[0][k] for k in range(3))
    d2 = tuple(Q[1][k]-Q[0][k] for k in range(3))
    r = tuple(P[0][k]-Q[0][k] for k in range(3))
    dt = lambda u, v: sum(u[k]*v[k] for k in range(3))
    a, e2, f2 = dt(d1, d1), dt(d2, d2), dt(d2, r)
    c, b = dt(d1, r), dt(d1, d2)
    denom = a*e2-b*b
    clamp = lambda t: min(max(t, mp.mpf(0)), mp.mpf(1))
    s = clamp((b*f2-c*e2)/denom) if denom > mp.mpf(10)**-60*a*e2 else mp.mpf(0)
    t = (b*s+f2)/e2
    if t < 0:
        t, s = mp.mpf(0), clamp(-c/a)
    elif t > 1:
        t, s = mp.mpf(1), clamp((b-c)/a)
    return mp.sqrt(sum((P[0][k]+s*d1[k]-Q[0][k]-t*d2[k])**2 for k in range(3)))


def margins(x, y):
    """`(clearance, contacts, sine)`: the least positive distance of either
    input's vertices from the other's faces and of their edges from each
    other relative to the size, the vertices exactly on a face of the other
    and the edges meeting, and the least sine between faces of different
    inputs whose planes each cross the other face."""
    size = pr.M(max(x.size(), y.size()))
    least, contacts = mp.mpf(1), 0
    for a, b in ((x, y), (y, x)):
        for p in a.vertices():
            on = False
            for g in b.faces:
                if on_face(p, g):
                    on = True
                    continue
                d = point_face_distance(p, g)/size
                if d > 0:
                    least = min(least, d)
            contacts += on
            if on and DEBUG:
                print('contact', [float(c) for c in p])
    for e in x.edges():
        for f in y.edges():
            d = segment_distance(e, f)/size
            if d > mp.mpf(10)**-30:
                least = min(least, d)
            else:
                contacts += 1
                if DEBUG:
                    print('edge contact', [[float(c) for c in p] for p in e], [[float(c) for c in p] for p in f])
    sine = mp.mpf(1)

    def crosses(face, n, a):
        s = [dot(n, sub(p, a)) for p in face]
        return min(s) < 0 < max(s)
    for f in x.faces:
        nf = normal_of(f)
        for g in y.faces:
            ng = normal_of(g)
            c = cross(nf, ng)
            if c == (0, 0, 0) or not crosses(f, ng, g[0]) or not crosses(g, nf, f[0]):
                continue
            s = mp.sqrt(pr.M(dot(c, c))/(pr.M(dot(nf, nf))*pr.M(dot(ng, ng))))
            sine = min(sine, s)
    return least, contacts, sine


def shares_plane(x, y):
    """Whether a face of one lies on a plane of a face of the other."""
    for f in x.faces:
        nf = normal_of(f)
        for g in y.faces:
            if cross(nf, normal_of(g)) == (0, 0, 0) and dot(nf, sub(g[0], f[0])) == 0:
                return True
    return False


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(inputs, expr, n, seed):
    """Volume and centre estimates of a Boolean expression over inputs
    (`expr`: an input's index, or `(op, e1, e2)`) from uniform points in
    the inputs' common box, membership by each input's cells' half-spaces
    in binary64: `(volume, centre, volume_sigma, centre_sigmas)`."""
    import random
    import math
    rnd = random.Random(seed)
    cells = []
    for body in inputs:
        cs = []
        for c in body.cells:
            cs.append([(tuple(float(v) for v in nn), tuple(float(v) for v in a)) for nn, a in pr.planes(c)])
        cells.append(cs)
    pts = [p for b in inputs for p in b.points()]
    lo = [min(float(p[k]) for p in pts) for k in range(3)]
    hi = [max(float(p[k]) for p in pts) for k in range(3)]
    box = (hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2])

    def member(k, X):
        for cell in cells[k]:
            if all(sum(nn[i]*(X[i]-a[i]) for i in range(3)) <= 0 for nn, a in cell):
                return True
        return False

    def holds(e, X):
        if isinstance(e, int):
            return member(e, X)
        op, e1, e2 = e
        return OPS[op](holds(e1, X), holds(e2, X))
    k, s, s2 = 0, [0.0]*3, [0.0]*3
    for _ in range(n):
        X = [lo[i]+(hi[i]-lo[i])*rnd.random() for i in range(3)]
        if holds(expr, X):
            k += 1
            for i in range(3):
                s[i] += X[i]
                s2[i] += X[i]*X[i]
    p = k/n
    sv = box*math.sqrt(max(p*(1-p), 0.0)/n)
    if k == 0:
        return 0.0, None, sv, None
    cm = [s[i]/k for i in range(3)]
    sc = [math.sqrt(max(s2[i]/k-cm[i]*cm[i], 0.0)/k) for i in range(3)]
    return box*p, cm, sv, sc


def number(x):
    return pr.number(x)


def rows(n, vol, area, centre):
    if n == 0:
        return ['empty']
    return [' '.join(['result', str(n), number(vol), number(area)]+[number(c) for c in centre])]
