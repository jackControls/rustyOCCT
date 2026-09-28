"""Independent reference for tessellation (T-a of REVIEW_NOTES.md).

From an identity case alone (identity_reference.Case: a prism of a profile
with lines, arcs, circles and holes, a box, a cone, a sphere or a torus
builder's arguments, or a face body of a profile), never from a kernel or
OCCT result, this module derives the body's boundary and its exact
properties:

* the boundary as patches with exact point distances: a prism's caps
  (planar regions with arcs; the distance combines the height with the
  in-plane distance to the region), its walls (each profile piece times the
  height range: rectangles and cylinder patches, a product of sets in
  orthogonal subspaces), and a solid of revolution's meridian pieces in the
  half-plane of the point (the nearest point of a surface of revolution
  lies in the point's own meridian half-plane); a torus wedge adds its end
  discs and restricts its tube to the swept angle;
* exact area and volume (mpmath: profile areas with circular segments and
  arc lengths; primitive_reference.py for cones, spheres and tori), the
  Euler characteristic of the boundary (2 - 2 genus for a closed shell,
  1 - holes for a face) and a face's number of boundary loops;
* point membership in the solid (for the orientation of triangles);
* each patch's whole surface (its line or circle times the height, or
  revolved), with exact point distances.

`check` measures a mesh against them: every mesh edge in exactly two
triangles in opposite directions (a face: interior edges so, boundary edges
once, in the expected number of loops), the Euler characteristic, every
node on the boundary, 12 barycentric samples of every triangle within the
request and within the triangle's reported bound of its face's whole
surface (the patch nearest the face's triangles) and within twice the
request of the boundary (the slivers between pcurves and their chords lie on
the surface just outside the face), edge polylines within the request of the
boundary, each solid triangle's normal leaving the solid, and the enclosed
volume within `δ (A + A_mesh)` of the exact one. Distances of samples are
taken in binary64 (the bounds are at least 1e-4 of the size; rounding is
covered by `1e-12` of the size), exact properties in mpmath. Nothing here
reads a kernel result except the mesh under test.
"""
from dataclasses import dataclass
import math

import mpmath as mp

import identity_reference as ident
import primitive_reference as prim

mp.mp.dps = 40
HALF_PI = 1.5707963267948966
TWO_PI = 6.283185307179586
# Rounding allowance relative to the case's scale.
ROUNDING = 1e-12
# Displacement of a point before a membership test, relative to its size.
NUDGE = 1e-14
# Barycentric samples of a triangle: the order-4 lattice without vertices.
SAMPLES = [(i/4, j/4, (4-i-j)/4) for i in range(5) for j in range(5-i)
           if (i, j, 4-i-j) not in ((4, 0, 0), (0, 4, 0), (0, 0, 4))]


# ------------------------------------------------------------------ 2D paths

@dataclass
class Line:
    p: tuple
    q: tuple


@dataclass
class Arc:
    """From p to q about c with radius r, starting at angle a0 and turning by
    `sweep` (positive counter-clockwise); a full circle has sweep 2π and no
    chord."""
    c: tuple
    r: float
    a0: float
    sweep: float
    p: tuple
    q: tuple
    full: bool = False


def arc_between(p, q, c, r, ccw):
    """An arc of a stored path from p to q (identity_reference.arc_sweep)."""
    sweep = float(ident.arc_sweep(p, q, (c[0], c[1], r, ccw)))
    return Arc(tuple(map(float, c)), float(r), math.atan2(p[1]-c[1], p[0]-c[0]), sweep,
               tuple(map(float, p)), tuple(map(float, q)))


def circle(c, r):
    c = tuple(map(float, c))
    p = (c[0]+r, c[1])
    return Arc(c, float(r), 0.0, TWO_PI, p, p, full=True)


def segment_distance(x, p, q):
    dx, dy = q[0]-p[0], q[1]-p[1]
    ex, ey = x[0]-p[0], x[1]-p[1]
    ll = dx*dx+dy*dy
    t = 0.0 if ll == 0 else max(0.0, min(1.0, (ex*dx+ey*dy)/ll))
    return math.hypot(ex-t*dx, ey-t*dy)


def on_arc(a, x):
    """Whether x's direction from the centre lies within the arc's sweep."""
    if a.full:
        return True
    t = math.atan2(x[1]-a.c[1], x[0]-a.c[0])
    if a.sweep > 0:
        return (t-a.a0) % TWO_PI <= a.sweep
    return (a.a0-t) % TWO_PI <= -a.sweep


def piece_distance(piece, x):
    if isinstance(piece, Line):
        return segment_distance(x, piece.p, piece.q)
    if on_arc(piece, x):
        return abs(math.hypot(x[0]-piece.c[0], x[1]-piece.c[1])-piece.r)
    return min(math.dist(x, piece.p), math.dist(x, piece.q))


def cross2(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])


def arc_midpoint(a):
    m = a.a0+a.sweep/2
    return (a.c[0]+a.r*math.cos(m), a.c[1]+a.r*math.sin(m))


def in_segment(a, x):
    """Inside the region between an arc and its chord (the whole disc for a
    full circle)."""
    if math.hypot(x[0]-a.c[0], x[1]-a.c[1]) >= a.r:
        return False
    if a.full:
        return True
    return (cross2(a.p, a.q, x) > 0) == (cross2(a.p, a.q, arc_midpoint(a)) > 0)


class Path:
    """A closed 2D boundary: pieces joined end to end (or one full circle).
    A point is inside when it is inside the polygon of the pieces' ends
    exclusive-or inside an odd number of arcs' circular segments, which holds
    for a simple path whose segments meet nothing else. The test point is
    moved by `NUDGE` of its size first, so a point on a chord (inside the
    region, not on its boundary) is decided; a point within the nudge of
    the true boundary is within it of the boundary either way."""
    def __init__(self, pieces):
        self.pieces = pieces
        self.points = [pc.p for pc in pieces if not (isinstance(pc, Arc) and pc.full)]

    def inside(self, x):
        k = NUDGE*(1.0+abs(x[0])+abs(x[1]))
        x = (x[0]+0.7548776662466927*k, x[1]+0.5698402909980532*k)
        odd = False
        pts = self.points
        for i in range(len(pts)):
            a, b = pts[i], pts[(i+1) % len(pts)]
            if (a[1] > x[1]) != (b[1] > x[1]):
                t = a[0]+(x[1]-a[1])*(b[0]-a[0])/(b[1]-a[1])
                if t > x[0]:
                    odd = not odd
        for pc in self.pieces:
            if isinstance(pc, Arc) and in_segment(pc, x):
                odd = not odd
        return odd

    def distance(self, x):
        return min(piece_distance(pc, x) for pc in self.pieces)

    def bounds(self):
        xs, ys = [], []
        for pc in self.pieces:
            if isinstance(pc, Line):
                xs += [pc.p[0], pc.q[0]]
                ys += [pc.p[1], pc.q[1]]
            else:
                # The ends and the extreme points the sweep reaches.
                xs += [pc.p[0], pc.q[0]]
                ys += [pc.p[1], pc.q[1]]
                for dx, dy in ((1, 0), (0, 1), (-1, 0), (0, -1)):
                    x = (pc.c[0]+pc.r*dx, pc.c[1]+pc.r*dy)
                    if on_arc(pc, x):
                        xs.append(x[0])
                        ys.append(x[1])
        return (min(xs), min(ys)), (max(xs), max(ys))


def profile_paths(boundaries, tolerance):
    """The stored paths of a profile's boundaries (outer first) and each one's
    exact signed area and length (mpmath)."""
    out = []
    for b in boundaries:
        if b.circle is not None:
            cx, cy, r = b.circle
            R = mp.mpf(r)
            out.append((Path([circle((cx, cy), r)]), mp.pi*R*R, 2*mp.pi*R))
            continue
        pts, _ = ident.stored(b, tolerance)
        if b.segments is not None:
            pts, segs = pts
        else:
            segs = [None]*len(pts)
        pieces, length = [], mp.mpf(0)
        n = len(pts)
        for i in range(n):
            p, q = pts[i], pts[(i+1) % n]
            if segs[i] is None:
                pieces.append(Line(tuple(map(float, p)), tuple(map(float, q))))
                length += mp.sqrt((mp.mpf(q[0])-p[0])**2+(mp.mpf(q[1])-p[1])**2)
            else:
                cx, cy, r, ccw = segs[i]
                pieces.append(arc_between(p, q, (cx, cy), r, ccw))
                length += mp.mpf(r)*abs(ident.arc_sweep(p, q, segs[i]))
        area = abs(ident.path_area(pts, segs))/2
        out.append((Path(pieces), area, length))
    return out


# ------------------------------------------------------------------ vectors

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


# ------------------------------------------------------------------ bodies

def line_distance(x, p, q):
    """To the whole line through p and q."""
    dx, dy = q[0]-p[0], q[1]-p[1]
    return abs((x[0]-p[0])*dy-(x[1]-p[1])*dx)/math.hypot(dx, dy)


def piece_surface_distance(piece, x):
    """To the piece's whole line or circle."""
    if isinstance(piece, Line):
        return line_distance(x, piece.p, piece.q)
    return abs(math.hypot(x[0]-piece.c[0], x[1]-piece.c[1])-piece.r)


class Prism:
    """A profile between two heights on a frame (a box is a rectangle on the
    world frame); `sheet` for a face body at height 0. Its patches are the
    caps and one wall per profile piece, each as (distance to the patch,
    distance to its whole surface) of a point."""
    def __init__(self, frame_axes, paths, start, end, sheet=False):
        self.o, self.x, self.y, self.n = frame_axes
        self.paths = paths
        self.z0, self.z1 = min(start, end), max(start, end)
        self.sheet = sheet
        self.patches = []
        for c in ([self.z0] if sheet else [self.z0, self.z1]):
            self.patches.append((lambda p, c=c: self.cap(p, c),
                                 lambda p, c=c: abs(self.local(p)[1]-c)))
        if not sheet:
            for pth in paths:
                for pc in pth.pieces:
                    self.patches.append((lambda p, pc=pc: self.wall(p, pc),
                                         lambda p, pc=pc: piece_surface_distance(pc, self.local(p)[0])))

    def local(self, p):
        w = sub(p, self.o)
        return (dot(w, self.x), dot(w, self.y)), dot(w, self.n)

    def region_distance(self, q):
        outer, holes = self.paths[0], self.paths[1:]
        if outer.inside(q) and not any(h.inside(q) for h in holes):
            return 0.0
        return min(pth.distance(q) for pth in self.paths)

    def cap(self, p, c):
        q, z = self.local(p)
        return math.hypot(self.region_distance(q), z-c)

    def wall(self, p, piece):
        q, z = self.local(p)
        return math.hypot(piece_distance(piece, q), max(0.0, self.z0-z, z-self.z1))

    def distance(self, p):
        return min(trimmed(p) for trimmed, _ in self.patches)

    def inside(self, p):
        q, z = self.local(p)
        return self.z0 < z < self.z1 and self.region_distance(q) == 0.0 \
            and self.paths[0].inside(q)


class Revolved:
    """A meridian region in the half-plane (rho >= 0, z) revolved about the
    frame's normal through `angle` from its x axis. `path` bounds the region;
    its pieces on the axis are not boundary. A piece's whole surface is its
    whole line or circle revolved, the mirror image included."""
    def __init__(self, frame_axes, path, angle):
        self.o, self.x, self.y, self.n = frame_axes
        self.path = path
        self.angle = angle
        self.boundary = [pc for pc in path.pieces
                         if not (isinstance(pc, Line) and pc.p[0] == 0.0 and pc.q[0] == 0.0)]
        self.patches = []
        for pc in self.boundary:
            whole = (lambda p, pc=pc: min(piece_surface_distance(pc, (s*self.local(p)[0], self.local(p)[1]))
                                          for s in (1.0, -1.0)))
            if angle == TWO_PI:
                self.patches.append((lambda p, pc=pc: piece_distance(pc, self.local(p)[:2]), whole))
            else:
                self.patches.append((self.tube, whole))
        if angle != TWO_PI:
            for u in (0.0, angle):
                self.patches.append((lambda p, u=u: self.disc(p, u)[0],
                                     lambda p, u=u: abs(self.disc(p, u)[1])))

    def local(self, p):
        w = sub(p, self.o)
        a, b, z = dot(w, self.x), dot(w, self.y), dot(w, self.n)
        return math.hypot(a, b), z, math.atan2(b, a) % TWO_PI

    def disc(self, p, u):
        """(distance to a wedge's end disc, signed distance to its plane)."""
        (tube,) = self.boundary
        e = tuple(math.cos(u)*a+math.sin(u)*b for a, b in zip(self.x, self.y))
        t = tuple(-math.sin(u)*a+math.cos(u)*b for a, b in zip(self.x, self.y))
        w = sub(sub(p, self.o), tuple(tube.c[0]*k for k in e))
        a, b, off = dot(w, e), dot(w, self.n), dot(w, t)
        return math.hypot(max(0.0, math.hypot(a, b)-tube.r), off), off

    def tube(self, p):
        """To a wedge's tube: in the point's half-plane within the swept
        angle, else through the nearer end circle."""
        (tube,) = self.boundary
        rho, z, phi = self.local(p)
        if phi <= self.angle:
            return abs(math.hypot(rho-tube.c[0], z-tube.c[1])-tube.r)
        best = math.inf
        for u in (0.0, self.angle):
            e = tuple(math.cos(u)*a+math.sin(u)*b for a, b in zip(self.x, self.y))
            t = tuple(-math.sin(u)*a+math.cos(u)*b for a, b in zip(self.x, self.y))
            w = sub(sub(p, self.o), tuple(tube.c[0]*k for k in e))
            best = min(best, math.hypot(math.hypot(dot(w, e), dot(w, self.n))-tube.r, dot(w, t)))
        return best

    def distance(self, p):
        return min(trimmed(p) for trimmed, _ in self.patches)

    def inside(self, p):
        rho, z, phi = self.local(p)
        return phi < self.angle and self.path.inside((rho, z))


@dataclass
class Body:
    shape: object
    solid: bool
    euler: int
    loops: int          # a face body's boundary loops (0 for a solid)
    area: object        # mpmath
    volume: object      # mpmath (0 for a face)
    scale: float        # characteristic size (deflections are relative to it)
    extent: float       # size plus the largest coordinate: the rounding scale


def frame_axes(frame):
    o, x, y, n = ident.frame_axes(tuple(frame))
    return tuple(map(float, o)), x, y, n


def meridian_body(case):
    """A cone, sphere or torus case as a Revolved shape and its exact mass."""
    fr = frame_axes(case.frame)
    o, normal, xh = case.frame[0:3], case.frame[3:6], case.frame[6:9]
    if case.cone is not None:
        r1, r2, h = case.cone
        pts = [(0.0, 0.0), (r1, 0.0), (r2, h), (0.0, h)]
        pts = [p for k, p in enumerate(pts) if p != pts[k-1]]
        path = Path([Line(pts[k], pts[(k+1) % len(pts)]) for k in range(len(pts))])
        volume, area, _, _ = prim.mass(prim.Cone('', o, normal, xh, r1, r2, h))
        return Revolved(fr, path, TWO_PI), volume, area, 2, max(r1, r2, h)
    if case.sphere is not None:
        R, lo, hi = case.sphere
        pieces = []
        a, b = (R*math.cos(lo), R*math.sin(lo)), (R*math.cos(hi), R*math.sin(hi))
        a = (0.0, -R) if lo == -HALF_PI else a
        b = (0.0, R) if hi == HALF_PI else b
        if lo != -HALF_PI:
            pieces.append(Line((0.0, a[1]), a))
        pieces.append(Arc((0.0, 0.0), R, lo, hi-lo, a, b))
        if hi != HALF_PI:
            pieces.append(Line(b, (0.0, b[1])))
        pieces.append(Line(pieces[-1].q, pieces[0].p))
        volume, area, _, _ = prim.sphere_mass(prim.Sphere('', o, normal, xh, R, lo, hi))
        return Revolved(fr, Path(pieces), TWO_PI), volume, area, 2, R
    R, r, lo, hi, angle = case.torus
    t = prim.Torus('', o, normal, xh, R, r, lo, hi, angle)
    volume, area, _, _ = prim.torus_mass(t)
    if hi-lo == TWO_PI:
        path = Path([circle((R, 0.0), r)])
        euler = 0 if angle == TWO_PI else 2
    else:
        a = (R+r*math.cos(lo), r*math.sin(lo))
        b = (R+r*math.cos(hi), r*math.sin(hi))
        path = Path([Line((0.0, a[1]), a), Arc((R, 0.0), r, lo, hi-lo, a, b),
                     Line(b, (0.0, b[1])), Line((0.0, b[1]), (0.0, a[1]))])
        euler = 2
    return Revolved(fr, path, angle), volume, area, euler, R+r


def body(case):
    """The independent body of an identity case."""
    if case.box is not None:
        origin, size = case.box
        lo = tuple(o+min(s, 0.0) for o, s in zip(origin, size))
        w, d, h = (abs(s) for s in size)
        rect = ident.Boundary(points=[(0.0, 0.0), (w, 0.0), (w, d), (0.0, d)])
        (path, area, length), = profile_paths([rect], case.tolerance)
        shape = Prism((lo, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), [path], 0.0, h)
        scale = max(w, d, h)
        return Body(shape, True, 2, 0, 2*area+length*h, area*h, scale,
                    scale+max(abs(v) for v in lo))
    if case.cone is not None or case.sphere is not None or case.torus is not None:
        shape, volume, area, euler, scale = meridian_body(case)
        return Body(shape, True, euler, 0, area, volume, scale,
                    scale+max(abs(v) for v in case.frame[0:3]))
    paths = profile_paths(case.boundaries, case.tolerance)
    area = paths[0][1]-sum(a for _, a, _ in paths[1:])
    length = sum(l for _, _, l in paths)
    (x0, y0), (x1, y1) = paths[0][0].bounds()
    fr = frame_axes(case.frame)
    extent = max(abs(v) for v in case.frame[0:3])+max(abs(x0), abs(x1), abs(y0), abs(y1))
    holes = len(paths)-1
    if case.make == 'face':
        scale = max(x1-x0, y1-y0)
        shape = Prism(fr, [p for p, _, _ in paths], 0.0, 0.0, sheet=True)
        return Body(shape, False, 1-holes, 1+holes, area, mp.mpf(0), scale, scale+extent)
    h = abs(case.end-case.start)
    scale = max(x1-x0, y1-y0, h)
    shape = Prism(fr, [p for p, _, _ in paths], case.start, case.end)
    extent += max(abs(case.start), abs(case.end))
    return Body(shape, True, 2-2*holes, 0, 2*area+length*h, area*h, scale, scale+extent)


# ------------------------------------------------------------------ meshes

@dataclass
class Mesh:
    nodes: list
    triangles: list     # (a, b, c, face, deflection bound, angle bound)
    edges: list         # (edge, [node], deflection bound, angle bound)
    deflection: float
    angle: float


def parse_meshes(text):
    """{(case, setting): Mesh} from the probe's text: `mesh CASE SETTING
    deflection angle`, `v x y z` rows, `t a b c face bound angle` rows,
    `e edge bound angle nodes...` rows, `end`."""
    out, current, key = {}, None, None
    for line in text.splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == 'mesh':
            key = (w[1], w[2])
            current = Mesh([], [], [], float(w[3]), float(w[4]))
        elif w[0] == 'v':
            current.nodes.append(tuple(float(x) for x in w[1:4]))
        elif w[0] == 't':
            current.triangles.append((int(w[1]), int(w[2]), int(w[3]), int(w[4]), float(w[5]), float(w[6])))
        elif w[0] == 'e':
            current.edges.append((int(w[1]), [int(x) for x in w[4:]], float(w[2]), float(w[3])))
        elif w[0] == 'end':
            out[key] = current
        elif w[0] == 'error':
            out[(w[1], w[2])] = ' '.join(w[3:])
        else:
            raise ValueError('malformed mesh line '+line)
    return out


def boundary_loops(edges):
    """The number of cycles formed by directed boundary edges."""
    succ = {}
    for a, b in edges:
        succ.setdefault(a, []).append(b)
    seen, loops = set(), 0
    for a, b in edges:
        if (a, b) in seen:
            continue
        loops += 1
        cur = (a, b)
        while cur not in seen:
            seen.add(cur)
            nxt = succ.get(cur[1], [])
            nxt = [n for n in nxt if (cur[1], n) not in seen]
            if not nxt:
                break
            cur = (cur[1], nxt[0])
    return loops


def face_patches(shape, mesh):
    """{face: patch index}: the patch whose whole surface is nearest, summed
    over the face's triangle centroids (coincident surfaces are equivalent)."""
    sums = {}
    for a, b, c, face, _, _ in mesh.triangles:
        centre = tuple((x+y+z)/3 for x, y, z in zip(mesh.nodes[a], mesh.nodes[b], mesh.nodes[c]))
        row = sums.setdefault(face, [0.0]*len(shape.patches))
        for k, (_, whole) in enumerate(shape.patches):
            row[k] += whole(centre)
    return {face: min(range(len(row)), key=row.__getitem__) for face, row in sums.items()}


def check(case, deflection, angle, mesh):
    """(measurements, failures) of a mesh against the independent body.

    Per triangle, 12 samples: their distance to the whole surface of the
    triangle's face (its patch, chosen by `face_patches`) must be within the
    request and within the triangle's reported bound; their distance to the
    body's boundary within twice the request (the slivers between a face's
    pcurves and their chords lie on the surface outside the face)."""
    b = body(case)
    eps = ROUNDING*b.extent
    failures = []
    if mesh.deflection > deflection or mesh.angle > angle:
        failures.append('reported_over_request')
    directed = {}
    used = set()
    for t in mesh.triangles:
        a, bb, c = t[:3]
        if len({a, bb, c}) < 3:
            failures.append('degenerate_triangle')
            continue
        used.update((a, bb, c))
        for e in ((a, bb), (bb, c), (c, a)):
            directed[e] = directed.get(e, 0)+1
    undirected = {}
    for (a, bb), k in directed.items():
        key = (min(a, bb), max(a, bb))
        undirected[key] = undirected.get(key, 0)+k
    open_edges = [e for e, k in directed.items() if (e[1], e[0]) not in directed]
    if any(k > 1 for k in directed.values()):
        failures.append('misoriented')
    if any(k > 2 for k in undirected.values()):
        failures.append('non_manifold')
    loops = boundary_loops(open_edges)
    if b.solid and open_edges:
        failures.append('not_closed')
    if not b.solid and loops != b.loops:
        failures.append('boundary_loops')
    euler = len(used)-len(undirected)+len(mesh.triangles)
    if euler != b.euler:
        failures.append('euler')
    node_gap = max((b.shape.distance(mesh.nodes[i]) for i in used), default=0.0)
    if node_gap > eps:
        failures.append('node_off_boundary')
    patch = face_patches(b.shape, mesh)
    worst, off_face, unsound, area, volume, inward = 0.0, 0.0, 0, 0.0, 0.0, 0
    for a, bb, c, face, bound, _ in mesh.triangles:
        p, q, r = mesh.nodes[a], mesh.nodes[bb], mesh.nodes[c]
        _, whole = b.shape.patches[patch[face]]
        for l1, l2, l3 in SAMPLES:
            s = tuple(l1*x+l2*y+l3*z for x, y, z in zip(p, q, r))
            d = whole(s)
            worst = max(worst, d)
            if d > bound+eps:
                unsound += 1
            off_face = max(off_face, b.shape.distance(s))
        nrm = cross(sub(q, p), sub(r, p))
        area += norm(nrm)/2
        volume += dot(p, cross(q, r))/6
        if b.solid and norm(nrm) > 0:
            centre = tuple((x+y+z)/3 for x, y, z in zip(p, q, r))
            step = 2*deflection/norm(nrm)
            if b.shape.inside(tuple(x+step*v for x, v in zip(centre, nrm))):
                inward += 1
    if worst > deflection+eps:
        failures.append('deflection_exceeded')
    if off_face > 2*deflection+eps:
        failures.append('off_face')
    if unsound:
        failures.append('bound_unsound')
    if inward:
        failures.append('normal_inward')
    edge_worst = 0.0
    for _, nodes, bound, _ in mesh.edges:
        for i, j in zip(nodes, nodes[1:]):
            for t in (0.25, 0.5, 0.75):
                s = tuple((1-t)*x+t*y for x, y in zip(mesh.nodes[i], mesh.nodes[j]))
                d = b.shape.distance(s)
                edge_worst = max(edge_worst, d)
                if d > min(deflection, bound)+eps:
                    failures.append('edge_deflection')
    exact_area, exact_volume = float(b.area), float(b.volume)
    if b.solid and not abs(volume-exact_volume) <= deflection*(exact_area+area)+eps*exact_area:
        failures.append('volume')
    return {'nodes': len(used), 'triangles': len(mesh.triangles), 'euler': euler,
            'deflection': worst, 'off_face': off_face, 'edge_deflection': edge_worst,
            'node_gap': node_gap, 'area_error': (area-exact_area)/exact_area,
            'volume_error': (volume-exact_volume)/exact_volume if b.solid else 0.0,
            'inward': inward}, sorted(set(failures))
