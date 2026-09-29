#!/usr/bin/env python3
"""Independent reference for S9d.1 of REVIEW_NOTES.md: Booleans of a sphere
(whole, a cap or a zone between two latitudes, `Solid::sphere_with`)
against a polyhedral prism (a profile of lines, every face a plane) in any
relative position, every pair of faces meeting in a line (two planes) or a
circle (a plane and the sphere).

Each input is its construction's exact model, in rationals from the stored
binary64 data (`stored_axes`, the kernel's `Frame3::new`, not exactly
orthonormal). The prism: the points `o + u x + v y + w n` with `(u, v)` in
the profile (its stored points, `identity_reference.stored`) and `w`
between the offsets. The sphere: `|X - o|^2 <= R^2` about its frame's
stored origin, and for a cap or a zone the end planes `w = R sin(latitude)`
in the same affine coordinates (the heights as `Solid::sphere_with` stores
them, binary64 products; a pole has none): with `m = x * y`, `m . (X - o)`
between `h_lo (m . n)` and `h_hi (m . n)`. (In an exact frame whose origin
plus `h n` is exact the plane through `o + h n` normal to `n`, as the
kernel's cap face is built, is the same plane.) Nothing here uses the
kernel, a surface/surface intersection or an arrangement: a nonconvex
profile is cut into triangles (ear clipping, no new vertices, so
neighbouring pieces share whole faces), each prism piece a convex
polyhedron of rational half-spaces.

* **Slicing (volume, first moments, the sphere's face).** Both solids are
  sliced by the planes `d . X = s`, `d = m` (the zone's normal: its end
  planes are slices; any rational `d` is accepted, and a second direction
  checks the first). In orthonormal coordinates `(g, h)` of the plane the
  sphere's section is a disc (centre the centre's projection, radius
  `sqrt(R^2 - z^2)`, `z` the plane's distance from the centre) cut by the
  zone's end planes' lines when they are not slices, and each prism piece's
  section a convex polygon (its edges' crossings of the plane). Each
  slice's regions are bounded by segments and arcs of the one circle, and
  their areas and first moments are Green's theorem over them in closed
  form. The four operations are taken apart: the disc's boundary cut at
  every piece's edge lines and each piece's edges cut at the disc's circle
  and lines, every piece classified at its midpoint against the other
  region (`D` inside `P`, `P` inside `D`): common = `dD n P + dP n D`,
  fuse = `dD \\ P + dP \\ D`, `D - P` = `dD \\ P - dP n D`, `P - D` =
  `dP \\ D - dD n P` (pieces' shared faces cancel); the common a second
  way, each piece's polygon clipped by the zone's lines and then by the
  disc (its edges' chords and the circle's arcs inside it). The sphere's
  own face: by Archimedes the area element is `R dz dtheta`, so its part
  inside the prism is `R / |d|` times the integral of the angle of the
  section's circle inside the pieces (the arcs above, classified), and its
  part outside the rest. The integrands are analytic between breakpoints,
  all roots of exact polynomials of degree at most two in `s`: the plane
  tangent to the sphere; through a piece's vertex; through a zone plane's
  meeting with a piece's edge; through an edge's meetings with the sphere;
  tangent to the circle of a face's or a zone plane's section of the sphere
  (the line of the two planes at distance `R` from the centre); through the
  meetings of a zone plane, a face's plane and the sphere. Each interval is
  integrated by Gauss-Legendre after `s = a + (b - a)(1 - cos t)/2` (end
  point square roots analytic), doubling the rule (12 to 96 nodes) until
  successive estimates agree within 1e-33 of the case's size to the fourth,
  else halving the interval.
* **Planar faces.** Each face of the prism (a wall's parallelogram, a
  cap's convex pieces) in orthonormal coordinates of its plane against the
  sphere's section there (a disc, cut by the zone's lines), in closed form:
  inside by clipping each convex piece, outside by the boundaries cut and
  classified as above (so their sum checks the face's area), or on a zone's
  end plane (the planes equal in rationals) with the same or the opposite
  orientation. Each end disc of a cap or zone against the prism's section
  by its plane (every piece crossing it cut into a polygon, a piece's face
  in the plane on it with its orientation), the same way.
* **Solids** (combinatorics, exact where it can be). A piece meets the
  sphere's interior when the squared distance from the centre to the
  piece's polyhedron cut by the zone (every half-space moved inwards by
  1e-20 of the case's size, so a polyhedron without interior is empty) is
  below `R^2`, exactly (Fractions, the projection's active set enumerated).
  Common: the pieces meeting the sphere, joined across a shared face whose
  part inside the sphere has positive area (the same test on the face).
  Fuse: one solid when the inputs overlap or share a face of opposite
  orientation, else two. `S - K` (a convex prism): the sphere's boundary
  outside `K` is the union of the connected sets `dS n H_f` (the open
  outer half-spaces of `K`'s faces); by radial projection from a point
  inside both, its components are those of `S - K`, and two such sets meet
  when `S` meets `H_f n H_g` (the same distance test). `K - S`: likewise
  the components of `dK - S`, each face minus the sphere's convex section
  one component per maximal run of its boundary outside it (edges cut where
  they enter the sphere, runs continued through the vertices outside);
  pieces joined through their shared edges. Regions meeting along a curve
  or at a point are separate solids (the regularized Boolean).

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as the S9a to S9c
references. A result keeps (object `A`, tool `B`): fuse `A`'s outside and
shared same-orientation pieces and `B`'s outside, cut `A`'s outside and
shared opposite pieces and `B`'s inside, common `A`'s inside and shared
same pieces and `B`'s inside. Nothing is shared with another reference but
`stored` and `stored_axes`.
"""
from fractions import Fraction as F
import itertools
import math

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

from identity_reference import stored
from curve_surface_reference import stored_axes

mp.mp.dps = 40

OPS = ('fuse', 'cut', 'common')
HALF_PI = math.pi/2
TAU = None

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


def is_zero(v):
    return all(x == 0 for x in v)


def solve(A, b):
    """Gaussian elimination in Fractions; None when singular."""
    n = len(b)
    a = [list(A[i])+[b[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            return None
        a[col], a[piv] = a[piv], a[col]
        for r in range(n):
            if r != col and a[r][col] != 0:
                f = a[r][col]/a[col][col]
                a[r] = [x-f*y for x, y in zip(a[r], a[col])]
    return [a[i][n]/a[i][i] for i in range(n)]


def tau():
    global TAU
    if TAU is None or TAU[0] != mp.mp.prec:
        TAU = (mp.mp.prec, 2*mp.pi)
    return TAU[1]


# ------------------------------------------------------------------ exact quadratics in s

def compose(A, B, C, p0, p1):
    """`A t^2 + B t + C` at `t = p0 + p1 s`: ascending coefficients in s."""
    return [A*p0*p0+B*p0+C, 2*A*p0*p1+B*p1, A*p1*p1]


def roots2(c):
    """The real roots (mpf) of an exact polynomial of degree at most two."""
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    if len(c) <= 1:
        return []
    if len(c) == 2:
        return [-M(c[0])/M(c[1])]
    a, b, k = c[2], c[1], c[0]
    disc = b*b-4*a*k
    if disc < 0:
        return []
    if disc == 0:
        return [-M(b)/(2*M(a))]
    sq = mp.sqrt(M(disc))
    return [(-M(b)-sq)/(2*M(a)), (-M(b)+sq)/(2*M(a))]


# ------------------------------------------------------------------ exact distance to a polyhedron

def dist2(c, ineqs, eqs=()):
    """The squared distance (Fraction) from `c` to `{a . X >= b} n {a . X =
    b}`, or None when it is empty: the projection's active set enumerated
    (every set of at most three independent constraints, the equalities
    always in, multipliers of the inequalities nonnegative, the point
    feasible)."""
    best = None
    k = 3-len(eqs)
    for size in range(0, k+1):
        for subset in itertools.combinations(range(len(ineqs)), size):
            active = list(eqs)+[ineqs[i] for i in subset]
            if active:
                G = [[dot(a1, a2) for a2, _ in active] for a1, _ in active]
                lam = solve(G, [b-dot(a, c) for a, b in active])
                if lam is None:
                    continue
                if any(l < 0 for l in lam[len(eqs):]):
                    continue
                X = c
                for l, (a, _) in zip(lam, active):
                    X = add(X, scale(a, l))
            else:
                X = c
            if any(dot(a, X) != b for a, b in eqs) or any(dot(a, X) < b for a, b in ineqs):
                continue
            d = dot(sub(X, c), sub(X, c))
            if best is None or d < best:
                best = d
    return best


# ------------------------------------------------------------------ inputs

class Sphere:
    """A sphere, cap or zone on its exact model: centre `c`, radius `r`, and
    the end planes as half-spaces `a . X >= b` (the zone on their positive
    sides), each with its end (`low`, `high`)."""

    def __init__(self, case):
        o, x, y, n = (tuple(F(v) for v in w) for w in stored_axes(case.frame))
        radius, low, high = case.sphere
        self.case = case
        self.c, self.r = o, F(radius)
        self.r2 = self.r*self.r
        self.m = cross(x, y)
        mn = dot(self.m, n)
        assert mn > 0
        self.heights = (None if low == -HALF_PI else F(radius*math.sin(low)),
                        None if high == HALF_PI else F(radius*math.sin(high)))
        self.planes = []
        mo = dot(self.m, o)
        if self.heights[0] is not None:
            self.planes.append((self.m, mo+self.heights[0]*mn, 'low'))
        if self.heights[1] is not None:
            self.planes.append((scale(self.m, -1), -(mo+self.heights[1]*mn), 'high'))
        self.mn = mn
        self.cm, self.rm = Mv(self.c), M(self.r)

    def z_range(self):
        """The ends' signed distances from the centre along `m` (mpf)."""
        mm = mp.sqrt(M(dot(self.m, self.m)))
        lo = -self.rm if self.heights[0] is None else M(self.heights[0]*self.mn)/mm
        hi = self.rm if self.heights[1] is None else M(self.heights[1]*self.mn)/mm
        return lo, hi

    def closed(self):
        """Volume, first moments and area in closed form."""
        r = self.rm
        lo, hi = self.z_range()
        V = mp.pi*((r*r*hi-hi**3/3)-(r*r*lo-lo**3/3))
        Mz = mp.pi*((r*r*hi**2/2-hi**4/4)-(r*r*lo**2/2-lo**4/4))
        mhat = scale(Mv(self.m), 1/mp.sqrt(M(dot(self.m, self.m))))
        mom = tuple(self.cm[i]*V+mhat[i]*Mz for i in range(3))
        area = 2*mp.pi*r*(hi-lo)
        for z, end in ((lo, self.heights[0]), (hi, self.heights[1])):
            if end is not None:
                area += mp.pi*(r*r-z*z)
        return V, mom, area

    def ineqs(self):
        return [(a, b) for a, b, _ in self.planes]


def ear_clip(ring):
    """Triangles (index triples, counter-clockwise) of a simple polygon given
    counter-clockwise, no new vertices."""
    idx = list(range(len(ring)))
    out = []
    while len(idx) > 3:
        for k in range(len(idx)):
            i, j, l = idx[k-1], idx[k], idx[(k+1) % len(idx)]
            a, b, c = ring[i], ring[j], ring[l]
            if cross2(sub(b, a), sub(c, b)) <= 0:
                continue
            inside = False
            for q in idx:
                if q in (i, j, l):
                    continue
                p = ring[q]
                if cross2(sub(b, a), sub(p, a)) >= 0 and cross2(sub(c, b), sub(p, b)) >= 0 \
                        and cross2(sub(a, c), sub(p, c)) >= 0:
                    inside = True
                    break
            if not inside:
                out.append([i, j, l])
                idx.pop(k)
                break
        else:
            raise AssertionError('no ear')
    out.append(idx)
    return out


class Piece:
    """A convex prism piece: its vertices (global keys `(k, level)`), faces
    (a loop of keys, the half-space `a . X >= b` holding the piece, and a
    tag: `('cap', level)`, `('wall', k)` for the profile's edge `k`, or
    `('internal', (i, j))` for a diagonal shared with another piece) and
    edges (sorted key pairs)."""

    def __init__(self, prism, indices):
        self.keys = [(k, lv) for lv in (0, 1) for k in indices]
        P = prism.P
        centroid = scale(sum_vec([P[k] for k in self.keys]), F(1, len(self.keys)))
        self.centroid = centroid
        n = len(indices)
        N = len(prism.ring)
        faces = [([(k, 0) for k in reversed(indices)], ('cap', 0)),
                 ([(k, 1) for k in indices], ('cap', 1))]
        for t in range(n):
            i, j = indices[t], indices[(t+1) % n]
            tag = ('wall', i) if j == (i+1) % N else ('internal', (min(i, j), max(i, j)))
            faces.append(([(i, 0), (j, 0), (j, 1), (i, 1)], tag))
        self.faces = []
        for loop, tag in faces:
            p0, p1, p2 = P[loop[0]], P[loop[1]], P[loop[2]]
            a = cross(sub(p1, p0), sub(p2, p0))
            b = dot(a, p0)
            if dot(a, centroid) < b:
                a, b = scale(a, -1), -b
            assert dot(a, centroid) > b
            self.faces.append({'loop': loop, 'a': a, 'b': b, 'tag': tag})
        self.edges = set()
        for f in self.faces:
            L = f['loop']
            for t in range(len(L)):
                self.edges.add(tuple(sorted((L[t], L[(t+1) % len(L)]))))

    def ineqs(self):
        return [(f['a'], f['b']) for f in self.faces]


def sum_vec(vs):
    out = (0, 0, 0)
    for v in vs:
        out = add(out, v)
    return out


class Prism:
    """A prism of one polygon (lines only) on its exact model, its convex
    pieces and its faces."""

    def __init__(self, case):
        self.case = case
        o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(case.frame))
        self.o, self.x, self.y, self.n = o, x, y, n
        self.lo, self.hi = sorted((F(case.start), F(case.end)))
        assert len(case.boundaries) == 1, 'S9d.1 reference: one polygon'
        b = case.boundaries[0]
        pts, _ = stored(b, case.tolerance)
        if b.segments is not None:
            points, segs = pts
            assert all(s is None for s in segs), 'S9d.1: lines only'
        else:
            points = pts
        self.ring = [(F(u), F(v)) for u, v in points]
        N = len(self.ring)
        self.det = det3(x, y, n)
        assert self.det > 0
        levels = (self.lo, self.hi)
        self.P = {(k, lv): self.world(self.ring[k], levels[lv]) for k in range(N) for lv in (0, 1)}
        convex = all(cross2(sub(self.ring[k], self.ring[k-1]), sub(self.ring[(k+1) % N], self.ring[k])) > 0
                     for k in range(N))
        self.convex = convex
        self.tris = [list(range(N))] if convex else ear_clip(self.ring)
        self.pieces = [Piece(self, t) for t in self.tris]

    def world(self, uv, w):
        return tuple(self.o[i]+uv[0]*self.x[i]+uv[1]*self.y[i]+w*self.n[i] for i in range(3))

    def profile_area(self, ring=None):
        ring = ring or self.ring
        return sum((cross2(ring[k-1], ring[k]) for k in range(len(ring))), F(0))/2

    def closed(self):
        """Volume, first moments and area in closed form."""
        A = self.profile_area()
        N = len(self.ring)
        mu = sum(((self.ring[k-1][0]+self.ring[k][0])*cross2(self.ring[k-1], self.ring[k]) for k in range(N)),
                 F(0))/6
        mv = sum(((self.ring[k-1][1]+self.ring[k][1])*cross2(self.ring[k-1], self.ring[k]) for k in range(N)),
                 F(0))/6
        h = self.hi-self.lo
        V = A*h*self.det
        cen = self.world((mu/A, mv/A), (self.lo+self.hi)/2)
        mom = tuple(M(V)*M(ci) for ci in cen)
        xy = cross(self.x, self.y)
        area = 2*M(A)*mp.sqrt(M(dot(xy, xy)))
        for k in range(N):
            e = sub(self.ring[(k+1) % N], self.ring[k])
            ew = add(scale(self.x, e[0]), scale(self.y, e[1]))
            en = cross(ew, self.n)
            area += mp.sqrt(M(dot(en, en)))*M(h)
        return M(V), mom, area

    def faces(self):
        """The prism's faces for areas: (tag, [convex polygons of 3D
        points], a piece holding it, exact area)."""
        out = []
        N = len(self.ring)
        h = self.hi-self.lo
        for k in range(N):
            quad = [self.P[(k, 0)], self.P[((k+1) % N, 0)], self.P[((k+1) % N, 1)], self.P[(k, 1)]]
            piece = next(p for p in self.pieces if any(f['tag'] == ('wall', k) for f in p.faces))
            e = sub(self.ring[(k+1) % N], self.ring[k])
            en = cross(add(scale(self.x, e[0]), scale(self.y, e[1])), self.n)
            out.append((('wall', k), [quad], piece, mp.sqrt(M(dot(en, en)))*M(h)))
        xy = cross(self.x, self.y)
        for lv in (0, 1):
            polys = [[self.P[(k, lv)] for k in t] for t in self.tris]
            out.append((('cap', lv), polys, self.pieces[0], M(self.profile_area())*mp.sqrt(M(dot(xy, xy)))))
        return out


# ------------------------------------------------------------------ 2D: segments, arcs, regions

def green_seg(p, q):
    """Green's integrals over the segment from `p` to `q`: `(x dy - y dx)/2`,
    `x^2 dy / 2` and `-y^2 dx / 2` (the arcs' fields: pieces of one closed
    boundary must share them)."""
    cr = p[0]*q[1]-q[0]*p[1]
    dx, dy = q[0]-p[0], q[1]-p[1]
    return (cr/2, dy*(p[0]*p[0]+p[0]*q[0]+q[0]*q[0])/6, -dx*(p[1]*p[1]+p[1]*q[1]+q[1]*q[1])/6)


def green_arc(C, rho, a, b):
    """Green's integrals over the circle `C + rho (cos t, sin t)`, `t` from
    `a` to `b` (counter-clockwise when `b > a`)."""
    cu, cv, r = C[0], C[1], rho
    sa, sb, ca, cb = mp.sin(a), mp.sin(b), mp.cos(a), mp.cos(b)
    A = (r*cu*(sb-sa)-r*cv*(cb-ca)+r*r*(b-a))/2
    ic = sb-sa
    ic2 = (b/2+mp.sin(2*b)/4)-(a/2+mp.sin(2*a)/4)
    ic3 = (sb-sb**3/3)-(sa-sa**3/3)
    is1 = -cb+ca
    is2 = (b/2-mp.sin(2*b)/4)-(a/2-mp.sin(2*a)/4)
    is3 = (-cb+cb**3/3)-(-ca+ca**3/3)
    Mu = (cu*cu*r*ic+2*cu*r*r*ic2+r**3*ic3)/2
    Mv_ = (cv*cv*r*is1+2*cv*r*r*is2+r**3*is3)/2
    return (A, Mu, Mv_)


def vadd(a, b, k=1):
    return tuple(x+k*y for x, y in zip(a, b))


ZERO3 = None


def zero3():
    return (mp.mpf(0), mp.mpf(0), mp.mpf(0))


def line_of(p, q):
    """The line through `p`, `q` as `alpha x + beta y >= gamma` holding its
    left side."""
    al, be = -(q[1]-p[1]), q[0]-p[0]
    return (al, be, al*p[0]+be*p[1])


def circle_line_angles(C, rho, line):
    al, be, ga = line
    nn = mp.sqrt(al*al+be*be)
    if nn == 0:
        return []
    val = (ga-al*C[0]-be*C[1])/(rho*nn)
    if abs(val) >= 1:
        return []
    phi = mp.atan2(be, al)
    w = mp.acos(val)
    return [phi-w, phi+w]


def seg_circle_params(p, q, C, rho2):
    """Parameters `t` where `p + t (q - p)` meets the circle."""
    e = (q[0]-p[0], q[1]-p[1])
    w = (p[0]-C[0], p[1]-C[1])
    A = e[0]*e[0]+e[1]*e[1]
    B = 2*(w[0]*e[0]+w[1]*e[1])
    Cc = w[0]*w[0]+w[1]*w[1]-rho2
    disc = B*B-4*A*Cc
    if disc <= 0 or A == 0:
        return []
    sq = mp.sqrt(disc)
    return [(-B-sq)/(2*A), (-B+sq)/(2*A)]


class Disc2:
    """A disc `|p - C|^2 < rho2` cut by half-planes `alpha x + beta y >=
    gamma` (convex; empty when `rho2 <= 0`)."""

    def __init__(self, C, rho2, lines=()):
        self.C, self.rho2, self.lines = C, rho2, list(lines)
        self.rho = mp.sqrt(rho2) if rho2 > 0 else mp.mpf(0)

    def empty(self):
        return self.rho2 <= 0

    def contains(self, p):
        if self.rho2 <= 0:
            return False
        if (p[0]-self.C[0])**2+(p[1]-self.C[1])**2 >= self.rho2:
            return False
        return all(al*p[0]+be*p[1] > ga for al, be, ga in self.lines)

    def point(self, t):
        return (self.C[0]+self.rho*mp.cos(t), self.C[1]+self.rho*mp.sin(t))

    def in_lines(self, p):
        return all(al*p[0]+be*p[1] > ga for al, be, ga in self.lines)

    def boundary(self):
        """[('arc', a, b) | ('seg', p, q)], the region on the left."""
        if self.rho2 <= 0:
            return []
        angles = []
        for ln in self.lines:
            angles += circle_line_angles(self.C, self.rho, ln)
        T = tau()
        out = []
        if not angles:
            if self.in_lines(self.point(mp.mpf(0))):
                return [('arc', mp.mpf(0), T)]
            return []
        angles = sorted(a % T for a in angles)
        for i in range(len(angles)):
            a = angles[i]
            b = angles[i+1] if i+1 < len(angles) else angles[0]+T
            if b-a <= 0:
                continue
            if self.in_lines(self.point((a+b)/2)):
                out.append(('arc', a, b))
        for k, (al, be, ga) in enumerate(self.lines):
            nn2 = al*al+be*be
            foot = (al*ga/nn2, be*ga/nn2)
            t = (be, -al)
            # |foot + tau t - C|^2 = rho2
            w = (foot[0]-self.C[0], foot[1]-self.C[1])
            A = nn2
            B = 2*(w[0]*t[0]+w[1]*t[1])
            Cc = w[0]**2+w[1]**2-self.rho2
            disc = B*B-4*A*Cc
            if disc <= 0:
                continue
            sq = mp.sqrt(disc)
            lo, hi = (-B-sq)/(2*A), (-B+sq)/(2*A)
            for j, (a2, b2, g2) in enumerate(self.lines):
                if j == k:
                    continue
                # a2 (foot + tau t) >= g2
                k1 = a2*t[0]+b2*t[1]
                k0 = a2*foot[0]+b2*foot[1]-g2
                if k1 == 0:
                    if k0 < 0:
                        lo, hi = mp.mpf(1), mp.mpf(0)
                    continue
                tt = -k0/k1
                if k1 > 0:
                    lo = max(lo, tt)
                else:
                    hi = min(hi, tt)
            if hi > lo:
                out.append(('seg', (foot[0]+lo*t[0], foot[1]+lo*t[1]), (foot[0]+hi*t[0], foot[1]+hi*t[1])))
        return out

    def split_lines(self):
        return self.lines


class Polys2:
    """A union of convex polygons with disjoint interiors (counter-
    clockwise vertex lists)."""

    def __init__(self, polys):
        self.polys = [p for p in polys if len(p) >= 3]
        self.plines = [[line_of(p[i], p[(i+1) % len(p)]) for i in range(len(p))] for p in self.polys]

    def contains(self, q):
        return any(all(al*q[0]+be*q[1] > ga for al, be, ga in ls) for ls in self.plines)

    def boundary(self):
        return [('seg', p[i], p[(i+1) % len(p)]) for p in self.polys for i in range(len(p))]

    def split_lines(self):
        return [ln for ls in self.plines for ln in ls]


def split_classify(curve, circle, other, other_circle):
    """Pieces of `curve` (on `circle` when an arc) cut at `other`'s lines (and
    its circle), each with whether its midpoint is inside `other`:
    [(green triple, inside, angle)]."""
    out = []
    if curve[0] == 'seg':
        _, p, q = curve
        ts = [mp.mpf(0), mp.mpf(1)]
        for al, be, ga in other.split_lines():
            fp, fq = al*p[0]+be*p[1]-ga, al*q[0]+be*q[1]-ga
            if (fp < 0 < fq) or (fq < 0 < fp):
                ts.append(fp/(fp-fq))
        if other_circle is not None and other_circle[1] > 0:
            ts += [t for t in seg_circle_params(p, q, other_circle[0], other_circle[1]) if 0 < t < 1]
        ts.sort()
        for t0, t1 in zip(ts, ts[1:]):
            if t1-t0 <= 0:
                continue
            a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
            b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
            mid = ((a[0]+b[0])/2, (a[1]+b[1])/2)
            out.append((green_seg(a, b), other.contains(mid), mp.mpf(0)))
        return out
    _, a0, a1 = curve
    C, rho = circle
    T = tau()
    angs = [a0, a1]
    for ln in other.split_lines():
        for t in circle_line_angles(C, rho, ln):
            k = mp.floor((a0-t)/T)+1
            t = t+k*T
            while t < a1:
                if t > a0:
                    angs.append(t)
                t += T
    angs.sort()
    for t0, t1 in zip(angs, angs[1:]):
        if t1-t0 <= 0:
            continue
        m = (t0+t1)/2
        mid = (C[0]+rho*mp.cos(m), C[1]+rho*mp.sin(m))
        out.append((green_arc(C, rho, t0, t1), other.contains(mid), t1-t0))
    return out


def classify(D, P):
    """Green's triples of `dD n P`, `dD \\ P`, `dP n D`, `dP \\ D`, and the
    angles of `D`'s arcs inside `P` and in all."""
    din, dout, pin, pout = zero3(), zero3(), zero3(), zero3()
    th_in, th_all = mp.mpf(0), mp.mpf(0)
    circ = (D.C, D.rho)
    for curve in D.boundary():
        for g, inside, ang in split_classify(curve, circ, P, None):
            if inside:
                din = vadd(din, g)
                th_in += ang
            else:
                dout = vadd(dout, g)
            th_all += ang
    dcirc = (D.C, D.rho2) if not D.empty() else None
    for curve in P.boundary():
        for g, inside, _ in split_classify(curve, None, D, dcirc):
            if inside and not D.empty():
                pin = vadd(pin, g)
            else:
                pout = vadd(pout, g)
    return din, dout, pin, pout, th_in, th_all


def clip_poly(poly, line):
    """The part of a convex polygon where `alpha x + beta y >= gamma`."""
    al, be, ga = line
    out = []
    n = len(poly)
    vals = [al*p[0]+be*p[1]-ga for p in poly]
    for i in range(n):
        p, q, fp, fq = poly[i], poly[(i+1) % n], vals[i], vals[(i+1) % n]
        if fp >= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp/(fp-fq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return out if len(out) >= 3 else []


def convex_disc(poly, D):
    """Green's triple of a convex polygon (counter-clockwise) inside `D`: the
    polygon clipped by `D`'s lines, then its edges' chords inside the circle
    and the circle's arcs inside it."""
    if D.empty():
        return zero3()
    for ln in D.lines:
        poly = clip_poly(poly, ln)
        if not poly:
            return zero3()
    lines = [line_of(poly[i], poly[(i+1) % len(poly)]) for i in range(len(poly))]
    inside = lambda q: all(al*q[0]+be*q[1] > ga for al, be, ga in lines)
    total = zero3()
    angs = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i+1) % n]
        ts = seg_circle_params(p, q, D.C, D.rho2)
        pin = (p[0]-D.C[0])**2+(p[1]-D.C[1])**2 < D.rho2
        if len(ts) == 2:
            t0, t1 = max(ts[0], mp.mpf(0)), min(ts[1], mp.mpf(1))
        else:
            t0, t1 = (mp.mpf(0), mp.mpf(1)) if pin else (mp.mpf(1), mp.mpf(0))
        if t1 > t0:
            a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
            b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
            total = vadd(total, green_seg(a, b))
        for t in ts:
            if 0 < t < 1:
                X = (p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1]))
                angs.append(mp.atan2(X[1]-D.C[1], X[0]-D.C[0]))
    T = tau()
    if not angs:
        if inside(D.point(mp.mpf(0))):
            total = vadd(total, green_arc(D.C, D.rho, mp.mpf(0), T))
        return total
    angs = sorted(a % T for a in angs)
    for i in range(len(angs)):
        a = angs[i]
        b = angs[i+1] if i+1 < len(angs) else angs[0]+T
        if b > a and inside(D.point((a+b)/2)):
            total = vadd(total, green_arc(D.C, D.rho, a, b))
    return total


def hull2(points):
    """Counter-clockwise convex polygon of points in convex position (sorted
    by angle about their mean, duplicates dropped)."""
    if len(points) < 3:
        return []
    mx = sum(p[0] for p in points)/len(points)
    my = sum(p[1] for p in points)/len(points)
    pts = sorted(points, key=lambda p: mp.atan2(p[1]-my, p[0]-mx))
    out = []
    for p in pts:
        if not out or abs(p[0]-out[-1][0])+abs(p[1]-out[-1][1]) > mp.mpf(10)**-30:
            out.append(p)
    if len(out) > 2 and abs(out[0][0]-out[-1][0])+abs(out[0][1]-out[-1][1]) <= mp.mpf(10)**-30:
        out.pop()
    return out if len(out) >= 3 else []


def plane_basis(a):
    """Orthonormal (g, h) of the plane normal to `a` (mpf), `g x h` along `a`."""
    am = Mv(a)
    nn = mp.sqrt(dot(am, am))
    ah = scale(am, 1/nn)
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    e = min(axes, key=lambda v: abs(dot(ah, v)))
    g = cross(ah, e)
    g = scale(g, 1/mp.sqrt(dot(g, g)))
    h = cross(ah, g)
    return g, h


# ------------------------------------------------------------------ quadrature

_NODES = {}


def nodes(degree):
    key = (degree, mp.mp.prec)
    if key not in _NODES:
        _NODES[key] = GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)
    return _NODES[key]


def integrate(f, a, b, tol, depth=0):
    """The integral over [a, b] of the vector function `f` by Gauss-Legendre
    after `s = a + (b - a)(1 - cos t)/2`, doubling the rule (12, 24, 48, 96
    nodes) until two successive estimates agree within `tol`, else halving
    the interval: (estimate, largest difference)."""
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
            est = [jac*c for c in v] if est is None else [e+jac*c for e, c in zip(est, v)]
        if prev is not None:
            diff = max(abs(x-y) for x, y in zip(est, prev))
            if diff <= tol:
                return est, diff
        prev = est
    assert depth < 60, 'quadrature did not converge'
    m = (a+b)/2
    e1, d1 = integrate(f, a, m, tol*3/5, depth+1)
    e2, d2 = integrate(f, m, b, tol*3/5, depth+1)
    return [x+y for x, y in zip(e1, e2)], d1+d2


# ------------------------------------------------------------------ slicing

class Slicing:
    """Both solids sliced by `d . X = s`."""

    def __init__(self, S, P, d, size):
        self.S, self.P, self.d = S, P, tuple(F(v) for v in d)
        self.size = size
        self.dd = dot(self.d, self.d)
        self.dn = mp.sqrt(M(self.dd))
        self.g, self.h = plane_basis(self.d)
        self.C2 = (dot(S.cm, self.g), dot(S.cm, self.h))
        self.dc = dot(self.d, S.c)
        # Zone planes: lines in each slice, or bounds on s.
        self.zlines, self.sbounds = [], []
        for a, b, _ in S.planes:
            if is_zero(cross(a, self.d)):
                k = dot(a, self.d)/self.dd
                self.sbounds.append((k, b))
            else:
                am = Mv(a)
                self.zlines.append((dot(am, self.g), dot(am, self.h), b, dot(a, self.d)/self.dd))
        self.level = {key: dot(self.d, p) for key, p in P.P.items()}
        self.Pm = {k: Mv(v) for k, v in P.P.items()}
        self.quad_error = mp.mpf(0)
        self.breaks = None

    def disc(self, s):
        for k, b in self.sbounds:
            if M(k)*s <= M(b):
                return Disc2(self.C2, mp.mpf(-1))
        z = (s-M(self.dc))/self.dn
        rho2 = M(self.S.r2)-z*z
        lines = [(al, be, M(b)-s*M(k)) for al, be, b, k in self.zlines]
        return Disc2(self.C2, rho2, lines)

    def sections(self, s):
        out = []
        Pm = self.Pm
        for piece in self.P.pieces:
            pts = []
            for i, j in piece.edges:
                li, lj = M(self.level[i]), M(self.level[j])
                if (li < s < lj) or (lj < s < li):
                    t = (s-li)/(lj-li)
                    X = add(Pm[i], scale(sub(Pm[j], Pm[i]), t))
                    pts.append((dot(X, self.g), dot(X, self.h)))
            poly = hull2(pts)
            if poly:
                out.append(poly)
        return out

    def integrand(self, s):
        D = self.disc(s)
        polys = self.sections(s)
        P = Polys2(polys)
        din, dout, pin, pout, th_in, th_all = classify(D, P)
        conv = zero3()
        for poly in polys:
            conv = vadd(conv, convex_disc(poly, D))
        out = []
        for g in (din, dout, pin, pout, conv):
            out += [g[0], s*g[0], g[1], g[2]]
        return out+[th_in, th_all]

    def breakpoints(self):
        S, P, d = self.S, self.P, self.d
        polys = []
        polys.append([self.dc*self.dc-S.r2*self.dd, -2*self.dc, F(1)])
        levels = set(self.level.values())
        planes = [(a, b) for a, b, _ in S.planes]
        faces = [(f['a'], f['b']) for piece in P.pieces for f in piece.faces]
        edges = set(e for piece in P.pieces for e in piece.edges)
        for k, b in self.sbounds:
            levels.add(b/k)
        for i, j in edges:
            vi, vj = P.P[i], P.P[j]
            e = sub(vj, vi)
            de = dot(d, e)
            for a, b in planes:
                ae = dot(a, e)
                if ae != 0:
                    t = (b-dot(a, vi))/ae
                    if 0 <= t <= 1:
                        levels.add(dot(d, vi)+t*de)
            if de != 0:
                w = sub(vi, S.c)
                polys.append(compose(dot(e, e), 2*dot(w, e), dot(w, w)-S.r2, -dot(d, vi)/de, 1/de))
        for a, b in faces+planes:
            if is_zero(cross(a, d)):
                continue
            aa, ad = dot(a, a), dot(a, d)
            det = aa*self.dd-ad*ad
            r1 = b-dot(a, S.c)
            # aa (s - dc)^2 - 2 ad r1 (s - dc) + dd r1^2 - R^2 det = 0
            p = compose(aa, -2*ad*r1, self.dd*r1*r1-S.r2*det, -self.dc, F(1))
            polys.append(p)
        for (az, bz), (af, bf) in itertools.product(planes, faces):
            l = cross(az, af)
            if is_zero(l):
                continue
            X0 = solve([az, af, l], [bz, bf, dot(l, S.c)])
            dl = dot(d, l)
            if dl == 0:
                levels.add(dot(d, X0))
                continue
            w = sub(X0, S.c)
            polys.append(compose(dot(l, l), 2*dot(w, l), dot(w, w)-S.r2, -dot(d, X0)/dl, 1/dl))
        pts = [M(v) for v in levels]
        for p in polys:
            pts += roots2(p)
        lo_s = min(M(v) for v in self.level.values())
        hi_s = max(M(v) for v in self.level.values())
        slo, shi = M(self.dc)-M(S.r)*self.dn, M(self.dc)+M(S.r)*self.dn
        a, b = min(lo_s, slo), max(hi_s, shi)
        pts = sorted(p for p in pts if a <= p <= b)
        out = []
        gap = mp.mpf(10)**-30*self.size
        for p in [a]+pts+[b]:
            if not out or p-out[-1] > gap:
                out.append(p)
        self.raw_breaks = pts
        return out

    def measure(self):
        """Integrals over s of the integrand's components."""
        if self.breaks is None:
            self.breaks = self.breakpoints()
        tol = mp.mpf(10)**-33*self.size**4
        total = None
        for a, b in zip(self.breaks, self.breaks[1:]):
            est, diff = integrate(self.integrand, a, b, tol)
            if est is None:
                continue
            self.quad_error = max(self.quad_error, diff)
            total = est if total is None else [x+y for x, y in zip(total, est)]
        return total

    def results(self):
        """{name: (volume, moments)} for `S`, `P`, `common`, `fuse`, `S-P`,
        `P-S` and `common2` (by convex clipping), and the sphere's face's
        areas inside and outside the prism."""
        v = self.measure()
        G = {}
        for k, name in enumerate(('din', 'dout', 'pin', 'pout', 'conv')):
            G[name] = v[4*k:4*k+4]
        combos = {'S': [('din', 1), ('dout', 1)], 'P': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'S-P': [('dout', 1), ('pin', -1)], 'P-S': [('pout', 1), ('din', -1)],
                  'common2': [('conv', 1)]}
        dm = Mv(self.d)
        out = {}
        for name, parts in combos.items():
            A = sum(k*G[p][0] for p, k in parts)
            sA = sum(k*G[p][1] for p, k in parts)
            mx = sum(k*G[p][2] for p, k in parts)
            my = sum(k*G[p][3] for p, k in parts)
            vol = A/self.dn
            mom = tuple((dm[i]*sA/M(self.dd)+self.g[i]*mx+self.h[i]*my)/self.dn for i in range(3))
            out[name] = (vol, mom)
        th_in, th_all = v[20], v[21]
        k = M(self.S.r)/self.dn
        out['sphere_face'] = {'in': k*th_in, 'out': k*(th_all-th_in), 'total': k*th_all}
        return out


# ------------------------------------------------------------------ planar faces

def plane_disc(S, a, b):
    """The sphere's section by the plane `a . X = b` in the plane's
    orthonormal coordinates: (Disc2 or None, basis, coincident end plane's
    outward normal or None). None when the plane is outside the zone."""
    g, h = plane_basis(a)
    aa = dot(a, a)
    t = b-dot(a, S.c)
    rho2 = M(S.r2-t*t/aa)
    C = (dot(S.cm, g), dot(S.cm, h))
    lines, on = [], None
    for az, bz, _ in S.planes:
        if is_zero(cross(az, a)):
            k = dot(az, a)/aa
            level = k*b
            if level == bz:
                on = scale(az, -1)
            elif level < bz:
                return None, (g, h), None
            continue
        azm = Mv(az)
        X0 = scale(Mv(a), M(b/aa))
        lines.append((dot(azm, g), dot(azm, h), M(bz)-dot(azm, X0)))
    return Disc2(C, rho2, lines), (g, h), on


def to2(X, basis):
    Xm = Mv(X)
    return (dot(Xm, basis[0]), dot(Xm, basis[1]))


def ccw2(poly):
    A = sum(cross2(poly[i-1], poly[i]) for i in range(len(poly)))
    return poly if A > 0 else list(reversed(poly))


def region_classes(D, polys):
    """(inside by convex clipping, outside by the boundaries classified) of
    the union of convex polygons against a Disc2."""
    inside = zero3()
    for p in polys:
        inside = vadd(inside, convex_disc(p, D))
    din, dout, pin, pout, _, _ = classify(D, Polys2(polys))
    outside = vadd(pout, din, -1)
    return inside[0], outside[0]


def prism_face_classes(S, prism):
    """[(tag, {in, out, same, opp}, exact area)] for the prism's faces."""
    out = []
    for tag, polys, piece, area in prism.faces():
        p0, p1, p2 = polys[0][0], polys[0][1], polys[0][2]
        a = cross(sub(p1, p0), sub(p2, p0))
        b = dot(a, p0)
        if dot(a, piece.centroid) > b:
            a, b = scale(a, -1), -b
        # `a` outward now; the plane a . X = b.
        D, basis, on = plane_disc(S, a, b)
        p2d = [ccw2([to2(X, basis) for X in poly]) for poly in polys]
        cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        if D is None or D.empty():
            cls['out'] = sum(abs(sum(cross2(p[i-1], p[i]) for i in range(len(p)))/2) for p in p2d)
        else:
            inside, outside = region_classes(D, p2d)
            if on is not None:
                cls['same' if dot(on, a) > 0 else 'opp'] = inside
            else:
                cls['in'] = inside
            cls['out'] = outside
        out.append((tag, cls, area))
    return out


def end_disc_classes(S, prism):
    """[(end, {in, out, same, opp}, exact area)] for the zone's end discs."""
    out = []
    for az, bz, end in S.planes:
        D, basis, _ = plane_disc(S, az, bz)
        D = Disc2(D.C, D.rho2)
        area = mp.pi*D.rho2
        groups = {'in': [], 'same': [], 'opp': []}
        for piece in prism.pieces:
            vals = {k: dot(az, prism.P[k])-bz for k in piece.keys}
            pos = any(v > 0 for v in vals.values())
            neg = any(v < 0 for v in vals.values())
            zeros = [k for k, v in vals.items() if v == 0]
            if pos and neg:
                pts = [to2(prism.P[k], basis) for k in zeros]
                for i, j in piece.edges:
                    if (vals[i] < 0 < vals[j]) or (vals[j] < 0 < vals[i]):
                        t = vals[i]/(vals[i]-vals[j])
                        pts.append(to2(add(prism.P[i], scale(sub(prism.P[j], prism.P[i]), t)), basis))
                poly = hull2(pts)
                if poly:
                    groups['in'].append(poly)
            elif len(zeros) >= 3:
                poly = hull2([to2(prism.P[k], basis) for k in zeros])
                # A piece on the zone's side has its face's outward normal
                # along the disc's (-a); on the other side, opposite.
                groups['same' if pos else 'opp'].append(poly)
        cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        allp = groups['in']+groups['same']+groups['opp']
        for key, polys in groups.items():
            for p in polys:
                cls[key] += convex_disc(p, D)[0]
        if allp:
            din, dout, pin, pout, _, _ = classify(D, Polys2(allp))
            cls['out'] = dout[0]-pin[0]
        else:
            cls['out'] = area
        out.append((end, cls, area))
    return out


# ------------------------------------------------------------------ solids

class Solids:
    """Solid counts by convexity (see the module's note)."""

    def __init__(self, S, P, size):
        self.S, self.P = S, P
        self.delta = F(1, 10**20)*F(int(size)+1)

    def shrink(self, ineqs):
        out = []
        for a, b in ineqs:
            out.append((a, b+self.delta*max(abs(x) for x in a)))
        return out

    def meets(self, ineqs, eqs=()):
        """Whether the polyhedron (shrunk) cut by the zone meets the open ball."""
        d = dist2(self.S.c, self.shrink(ineqs+self.S.ineqs()), eqs)
        return d is not None and d < self.S.r2

    def overlap(self, piece):
        return self.meets(piece.ineqs())

    def common(self):
        pieces = [k for k, p in enumerate(self.P.pieces) if self.overlap(p)]
        parent = {k: k for k in pieces}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k
        for i, j in itertools.combinations(pieces, 2):
            pi, pj = self.P.pieces[i], self.P.pieces[j]
            shared = [f for f in pi.faces if f['tag'][0] == 'internal'
                      and any(g['tag'] == f['tag'] for g in pj.faces)]
            for f in shared:
                rest = [(g['a'], g['b']) for g in pi.faces if g is not f]
                if self.meets(rest, [(f['a'], f['b'])]):
                    parent[find(i)] = find(j)
        return len({find(k) for k in pieces})

    def fuse(self, opp_area):
        if any(self.overlap(p) for p in self.P.pieces) or opp_area > 0:
            return 1
        return 2

    def sphere_minus(self):
        """Components of `S - K`, `K` convex."""
        assert len(self.P.pieces) == 1, 'S9d.1 reference: S - P for a convex prism only'
        piece = self.P.pieces[0]
        if not self.overlap(piece):
            return 1
        outer = [(scale(f['a'], -1), -f['b']) for f in piece.faces]
        nodes = [k for k, h in enumerate(outer) if self.meets([h])]
        parent = {k: k for k in nodes}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k
        for i, j in itertools.combinations(nodes, 2):
            if self.meets([outer[i], outer[j]]):
                parent[find(i)] = find(j)
        return len({find(k) for k in nodes})

    def edge_outside(self, i, j):
        """The parts of the edge from key `i` to key `j` outside the sphere:
        (piece at i or None, piece at j or None), the same object when the
        whole edge is outside."""
        S = self.S
        vi, vj = Mv(self.P.P[i]), Mv(self.P.P[j])
        e = sub(vj, vi)
        w = sub(vi, S.cm)
        A, B, C = dot(e, e), 2*dot(w, e), dot(w, w)-M(S.r2)
        disc = B*B-4*A*C
        if disc <= 0:
            return ('whole', i, j), ('whole', i, j)
        sq = mp.sqrt(disc)
        t0, t1 = (-B-sq)/(2*A), (-B+sq)/(2*A)
        for a, b in S.ineqs():
            am = Mv(a)
            k1, k0 = dot(am, e), dot(am, vi)-M(b)
            if k1 == 0:
                if k0 < 0:
                    t0, t1 = mp.mpf(1), mp.mpf(0)
                continue
            tt = -k0/k1
            if k1 > 0:
                t0 = max(t0, tt)
            else:
                t1 = min(t1, tt)
        t0, t1 = max(t0, mp.mpf(0)), min(t1, mp.mpf(1))
        if t1-t0 <= mp.mpf(10)**-30:
            return ('whole', i, j), ('whole', i, j)
        return (('start', i, j) if t0 > 0 else None), (('end', i, j) if t1 < 1 else None)

    def prism_minus(self):
        """Components of `P - S`."""
        ends = {}
        for piece in self.P.pieces:
            for i, j in piece.edges:
                if (i, j) not in ends:
                    ends[(i, j)] = self.edge_outside(i, j)
        nodes = set(x for v in ends.values() for x in v if x is not None)
        parent = {k: k for k in nodes}

        def find(k):
            while parent[k] != k:
                k = parent[k]
            return k

        def at(edge_from, edge_to, key):
            e = tuple(sorted((edge_from, edge_to)))
            pa, pb = ends[e]
            return pa if key == e[0] else pb
        for piece in self.P.pieces:
            for f in piece.faces:
                L = f['loop']
                n = len(L)
                for t in range(n):
                    prev, cur, nxt = L[t-1], L[t], L[(t+1) % n]
                    a, b = at(prev, cur, cur), at(cur, nxt, cur)
                    if a is not None and b is not None:
                        parent[find(a)] = find(b)
        return len({find(k) for k in nodes})


# ------------------------------------------------------------------ the pair

def case_size(S, P):
    vals = [F(1)]
    vals += [abs(x) for p in P.P.values() for x in p]
    vals += [abs(x)+S.r for x in S.c]
    return M(max(vals))


class Pair:
    """A sphere and a prism, in either order (`obj`, `tool` Boolean cases)."""

    def __init__(self, obj, tool, d=None):
        self.sphere_first = obj.sphere is not None
        s_case, p_case = (obj, tool) if self.sphere_first else (tool, obj)
        assert s_case.sphere is not None and p_case.sphere is None, 'S9d.1: a sphere and a prism'
        self.S, self.P = Sphere(s_case), Prism(p_case)
        self.size = case_size(self.S, self.P)
        self.slicing = Slicing(self.S, self.P, d or self.S.m, self.size)
        self._res = self._faces = None

    def sliced(self):
        if self._res is None:
            self._res = self.slicing.results()
        return self._res

    def faces(self):
        """[(input 'S' or 'P', tag, classes, exact area)]."""
        if self._faces is None:
            out = [('P', tag, cls, area) for tag, cls, area in prism_face_classes(self.S, self.P)]
            out += [('S', ('end', end), cls, area) for end, cls, area in end_disc_classes(self.S, self.P)]
            sf = self.sliced()['sphere_face']
            lo, hi = self.S.z_range()
            out.append(('S', ('sphere',), {'in': sf['in'], 'out': sf['out'], 'same': mp.mpf(0),
                                           'opp': mp.mpf(0)}, 2*mp.pi*self.S.rm*(hi-lo)))
            self._faces = out
        return self._faces

    def classes(self, which):
        tot = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        for w, _, cls, _ in self.faces():
            if w == which:
                for k in tot:
                    tot[k] += cls[k]
        return tot

    def roles(self):
        return ('S', 'P') if self.sphere_first else ('P', 'S')

    def area(self, op):
        a, b = self.roles()
        ca, cb = self.classes(a), self.classes(b)
        keep = {'fuse': (('out', 'same'), ('out',)), 'cut': (('out', 'opp'), ('in',)),
                'common': (('in', 'same'), ('in',))}[op]
        return sum(ca[k] for k in keep[0])+sum(cb[k] for k in keep[1])

    def volume(self, op):
        r = self.sliced()
        key = {'fuse': 'fuse', 'common': 'common', 'cut': 'S-P' if self.sphere_first else 'P-S'}[op]
        return r[key]

    def solids(self, op):
        sol = Solids(self.S, self.P, self.size)
        if op == 'common':
            return sol.common()
        if op == 'fuse':
            return sol.fuse(self.classes('S')['opp'])
        return sol.sphere_minus() if self.sphere_first else sol.prism_minus()

    def result(self, op):
        vol, mom = self.volume(op)
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        n = self.solids(op)
        assert n > 0, 'a result of positive volume without solids'
        return n, vol, self.area(op), tuple(m/vol for m in mom)


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
