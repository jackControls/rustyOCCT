"""Independent reference for tessellation (T-a and T-b of REVIEW_NOTES.md).

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

T-b adds spline bodies (`SplineCase`, from their builder's arguments): prisms
whose profiles have planar B-spline pieces (cut exactly into Bézier pieces,
distances by projection, membership by the parity of each piece's roots,
exact areas and lengths by quadrature), a face in a periodic spline ring, a
box under a spline graph and face bodies on spline surfaces with a
rectangular hole in their parameters (distances by projection onto the
surface over its domain, trimmed through the projection's parameters, and
areas by Gauss-Legendre quadrature of the exact patches). A projection can
only overstate a distance, so a mesh found within its bound is.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import math

import mpmath as mp

import identity_reference as ident
import primitive_reference as prim
import spline_cell_reference as spl

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
    if isinstance(piece, Spline):
        return project_curve(piece.pieces, x)[0]
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
        self.points = [pc.p for pc in pieces if not (isinstance(pc, (Arc, Spline)) and pc.full)]

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
            elif isinstance(pc, Spline) and spline_segment(pc, x):
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
            elif isinstance(pc, Spline):
                # The control hull (the curve lies in it).
                for b in pc.pieces:
                    for h in b.h:
                        xs.append(h[0]/h[-1])
                        ys.append(h[1]/h[-1])
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


def nearest(patches, p, within=None):
    """The distance from p to the nearest trimmed patch; with `within`, the
    first patch at most that far (enough to decide a check)."""
    best = math.inf
    for trimmed, _ in patches:
        best = min(best, trimmed(p))
        if within is not None and best <= within:
            break
    return best


def piece_surface_distance(piece, x):
    """To the piece's whole line or circle (a spline piece: its whole curve,
    the fixtures' splines being their whole domains)."""
    if isinstance(piece, Spline):
        return project_curve(piece.pieces, x)[0]
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

    def distance(self, p, within=None):
        return nearest(self.patches, p, within)

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

    def distance(self, p, within=None):
        return nearest(self.patches, p, within)

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
    """The independent body of an identity case (or of a spline case, T-b)."""
    if isinstance(case, SplineCase):
        return spline_body(case)
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
    # A spline body's distances are projections: a sample of each face's
    # triangles (every k-th, at most 24) decides as well.
    stride = {}
    if getattr(shape, 'slow', False):
        for t in mesh.triangles:
            stride[t[3]] = stride.get(t[3], 0)+1
        stride = {f: -(-n//24) for f, n in stride.items()}
    seen = {}
    for a, b, c, face, _, _ in mesh.triangles:
        seen[face] = seen.get(face, -1)+1
        if stride and seen[face] % stride[face]:
            continue
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
    # A spline body decides a check at the first patch close enough (its
    # distances are projections); the analytic bodies measure every patch.
    slow = getattr(b.shape, 'slow', False)
    node_gap = max((b.shape.distance(mesh.nodes[i], eps if slow else None) for i in used), default=0.0)
    if node_gap > eps:
        failures.append('node_off_boundary')
    patch = face_patches(b.shape, mesh)
    worst, off_face, unsound, area, volume, inward = 0.0, 0.0, 0, 0.0, 0.0, 0
    for a, bb, c, face, bound, _ in mesh.triangles:
        p, q, r = mesh.nodes[a], mesh.nodes[bb], mesh.nodes[c]
        own, whole = b.shape.patches[patch[face]]
        here = 0.0
        for l1, l2, l3 in SAMPLES:
            s = tuple(l1*x+l2*y+l3*z for x, y, z in zip(p, q, r))
            d = whole(s)
            here = max(here, d)
            if d > bound+eps:
                unsound += 1
            if slow:
                off = own(s)
                off_face = max(off_face, off if off <= 2*deflection+eps else b.shape.distance(s))
            else:
                off_face = max(off_face, b.shape.distance(s))
        worst = max(worst, here)
        nrm = cross(sub(q, p), sub(r, p))
        area += norm(nrm)/2
        volume += dot(p, cross(q, r))/6
        if b.solid and norm(nrm) > 0:
            # Twice the triangle's own deviation, at least the request's:
            # a centroid that far inside the solid is still left by an
            # outward normal.
            centre = tuple((x+y+z)/3 for x, y, z in zip(p, q, r))
            step = 2*max(deflection, here)/norm(nrm)
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
                d = b.shape.distance(s, min(deflection, bound)+eps if slow else None)
                edge_worst = max(edge_worst, d)
                if d > min(deflection, bound)+eps:
                    failures.append('edge_deflection')
    exact_area, exact_volume = float(b.area), float(b.volume)
    if b.solid and not (volume > 0 and abs(volume-exact_volume) <= deflection*(exact_area+area)+eps*exact_area):
        failures.append('volume')
    return {'nodes': len(used), 'triangles': len(mesh.triangles), 'euler': euler,
            'deflection': worst, 'off_face': off_face, 'edge_deflection': edge_worst,
            'node_gap': node_gap, 'area_error': (area-exact_area)/exact_area,
            'volume_error': (volume-exact_volume)/exact_volume if b.solid else 0.0,
            'inward': inward}, sorted(set(failures))


# ------------------------------------------------------------------ splines (T-b)
#
# A spline is cut exactly (Fractions, spline_cell_reference's blossoms) into
# its Bézier pieces between knots; distances are projections in binary64:
# the nearest of a few samples per piece, then Newton's iteration on
# (C - x) . C' within the piece (Gauss-Newton on a surface's patches, clamped
# to the domain). A projection's distance is the length of a difference, so
# it can only overstate the true distance; a start in the wrong basin shows
# as a larger distance, never a smaller one. Exact areas and lengths are
# mpmath quadratures of the exact pieces.

@dataclass
class SplineCase:
    """A T-b fixture body, from its builder's arguments alone.

    kind 'prism': `boundaries` (generate_brep_fixtures.prism's piece lists:
    ('line', p), ('arc', p, centre, ccw), ('circle', centre, r, ccw),
    ('spline', p, interior poles, Basis, weights)) between heights z0 and z1,
    then the rigid `motion` (unit axis, angle, shift) if any; `variant`
    'stadium_edge' or 'stadium_pcurve' replaces a stadium's vertical edge or
    its pcurve on the cylinder by a degree-1 spline (the same geometry).
    kind 'face': a face body on the plane z = 0 bounded by the periodic
    `ring` (a spline_cell_reference.BSpline2 over its whole period).
    kind 'dome': the box [0, w] x [0, d] x [0, h] under `surface`, a spline
    graph S(u, v) = (u, v, z(u, v)) over the same rectangle equal to h on
    its boundary. kind 'sheet': a face body on the whole domain of
    `surface`, less the parameter rectangle `hole` (u0, u1, v0, v1) if any."""
    name: str
    kind: str
    boundaries: list = None
    z0: float = 0.0
    z1: float = 1.0
    motion: tuple = None
    variant: str = None
    ring: object = None
    surface: object = None
    size: tuple = None
    hole: tuple = None


def _bern(ctrl, t):
    """De Casteljau's value of Bernstein controls (tuples) at t."""
    row = list(ctrl)
    s = 1.0-t
    while len(row) > 1:
        row = [tuple(s*a+t*b for a, b in zip(x, y)) for x, y in zip(row, row[1:])]
    return row[0]


def _diff(ctrl, scale):
    if len(ctrl) < 2:
        return [tuple(0.0 for _ in ctrl[0])]
    return [tuple(scale*(b-a) for a, b in zip(x, y)) for x, y in zip(ctrl, ctrl[1:])]


def _mpf(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, F) else mp.mpf(x)


class Bezier:
    """A rational Bézier piece on [0, 1] in binary64: homogeneous controls
    (w x, w y[, w z], w) and their first and second differences; `exact`
    keeps the Fractions."""
    def __init__(self, hctrl):
        self.exact = hctrl
        self.h = [tuple(float(c) for c in x) for x in hctrl]
        p = len(self.h)-1
        self.d1 = _diff(self.h, p)
        self.d2 = _diff(self.d1, p-1)

    def point(self, t):
        a = _bern(self.h, t)
        return tuple(x/a[-1] for x in a[:-1])

    def jet(self, t):
        """(C, C', C'') at t, from C = A/w: C' = (A' - w' C)/w,
        C'' = (A'' - 2 w' C' - w'' C)/w."""
        a, a1, a2 = _bern(self.h, t), _bern(self.d1, t), _bern(self.d2, t)
        w, w1, w2 = a[-1], a1[-1], a2[-1]
        c = tuple(x/w for x in a[:-1])
        c1 = tuple((x-w1*y)/w for x, y in zip(a1[:-1], c))
        c2 = tuple((x-2*w1*y-w2*z)/w for x, y, z in zip(a2[:-1], c1, c))
        return c, c1, c2

    def mp_jet(self, t):
        """(C, C') at an mpmath t from the exact controls."""
        h = [tuple(_mpf(c) for c in x) for x in self.exact]
        p = len(h)-1
        a = _bern(h, t)
        a1 = _bern([tuple(p*(q-o) for o, q in zip(x, y)) for x, y in zip(h, h[1:])], t) if p else \
            tuple(mp.mpf(0) for _ in h[0])
        w, w1 = a[-1], a1[-1]
        c = tuple(x/w for x in a[:-1])
        return c, tuple((x-w1*y)/w for x, y in zip(a1[:-1], c))


def curve_pieces(curve):
    """The Bézier pieces of a spline_cell_reference curve between its knots,
    in its span's direction (exact blossoms)."""
    cuts = [F(0)]+spl.knot_fractions(curve)+[F(1)]
    return [Bezier(spl.blossom_piece(curve, a, b)) for a, b in zip(cuts, cuts[1:])]


def _vsub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def _vdot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _newton_curve(b, x, t):
    """Newton's iteration on (C - x) . C' from t, within [0, 1]."""
    for _ in range(60):
        c, c1, c2 = b.jet(t)
        r = _vsub(c, x)
        f = _vdot(r, c1)
        fp = _vdot(c1, c1)+_vdot(r, c2)
        if fp <= 0:
            fp = _vdot(c1, c1)
        if fp == 0:
            break
        nt = min(1.0, max(0.0, t-f/fp))
        if abs(nt-t) <= 1e-16:
            t = nt
            break
        t = nt
    return math.dist(b.point(t), x), t


def project_curve(pieces, x, samples=8):
    """(distance, piece, t): the nearest point of the pieces to x."""
    best = (math.inf, 0, 0.0)
    for k, b in enumerate(pieces):
        start = min(((math.dist(b.point(i/samples), x), i/samples) for i in range(samples+1)))
        d, t = _newton_curve(b, x, start[1])
        d = min(d, start[0])
        if d < best[0]:
            best = (d, k, t)
    return best


def _bernstein_halves(c):
    row = list(c)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(a+b)/2 for a, b in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return left, right[::-1]


def _crossings(coeffs, lo=0.0, hi=1.0, depth=0):
    """Midpoints of the intervals, halved to 2^-50, where a Bernstein
    polynomial changes sign between its ends; none where its coefficients
    keep a strict sign. Their number has the parity of its roots' count."""
    if all(c > 0 for c in coeffs) or all(c < 0 for c in coeffs):
        return []
    if depth == 50:
        return [(lo+hi)/2] if (coeffs[0] > 0) != (coeffs[-1] > 0) else []
    left, right = _bernstein_halves(coeffs)
    mid = (lo+hi)/2
    return _crossings(left, lo, mid, depth+1)+_crossings(right, mid, hi, depth+1)


def spline_segment(piece, x):
    """Inside the region between a spline piece and its chord (the curve's
    own interior for a closed ring): the parity of the crossings of the ray
    from x along +u with the curve (roots of w (y - x_v) on each Bézier
    piece) and the chord back from q to p."""
    odd = False
    for b in piece.pieces:
        coeffs = [h[1]-x[1]*h[-1] for h in b.h]
        for t in _crossings(coeffs):
            if b.point(t)[0] > x[0]:
                odd = not odd
    if not piece.full:
        a, bb = piece.q, piece.p
        if (a[1] > x[1]) != (bb[1] > x[1]):
            t = a[0]+(x[1]-a[1])*(bb[0]-a[0])/(bb[1]-a[1])
            if t > x[0]:
                odd = not odd
    return odd


class Spline:
    """A planar B-spline piece of a path from p to q (the whole closed curve
    when `full`)."""
    def __init__(self, curve, full=False):
        self.pieces = curve_pieces(curve)
        self.p = self.pieces[0].point(0.0)
        self.q = self.pieces[-1].point(1.0)
        self.full = full

    def exact_terms(self):
        """(∫ (x y' - y x'), length) over the curve, mpmath quadrature of
        the exact pieces."""
        twice, length = mp.mpf(0), mp.mpf(0)
        for b in self.pieces:
            def cross(t, b=b):
                (x, y), (dx, dy) = b.mp_jet(t)
                return x*dy-y*dx

            def speed(t, b=b):
                _, (dx, dy) = b.mp_jet(t)
                return mp.sqrt(dx*dx+dy*dy)
            twice += mp.quad(cross, [0, 1])
            length += mp.quad(speed, [0, 1])
        return twice, length


def spline_paths(boundaries):
    """[(Path, area, length)] of generate_brep_fixtures piece lists, exact
    (mpmath): the signed areas by ∮ (x dy - y dx) / 2 over each piece."""
    out = []
    for boundary in boundaries:
        if boundary[0][0] == 'circle':
            _, c, r, _ = boundary[0]
            R = mp.mpf(r)
            out.append((Path([circle(c, r)]), mp.pi*R*R, 2*mp.pi*R))
            continue
        if boundary[0][0] == 'ring':
            piece = Spline(boundary[0][1], full=True)
            twice, length = piece.exact_terms()
            out.append((Path([piece]), abs(twice)/2, length))
            continue
        n = len(boundary)
        pieces, twice, length = [], mp.mpf(0), mp.mpf(0)
        for i, piece in enumerate(boundary):
            p, q = piece[1], boundary[(i+1) % n][1]
            if piece[0] == 'line':
                pieces.append(Line(tuple(map(float, p)), tuple(map(float, q))))
                twice += mp.mpf(p[0])*q[1]-mp.mpf(q[0])*p[1]
                length += mp.sqrt((mp.mpf(q[0])-p[0])**2+(mp.mpf(q[1])-p[1])**2)
            elif piece[0] == 'arc':
                _, _, c, ccw = piece
                r = math.hypot(p[0]-c[0], p[1]-c[1])
                phi = ident.arc_sweep(p, q, (c[0], c[1], r, ccw))
                pieces.append(Arc(tuple(map(float, c)), r, math.atan2(p[1]-c[1], p[0]-c[0]), float(phi),
                                  tuple(map(float, p)), tuple(map(float, q))))
                # x dy - y dx over the arc: c x (q - p) + r^2 phi.
                twice += (mp.mpf(c[0])*(mp.mpf(q[1])-p[1])-mp.mpf(c[1])*(mp.mpf(q[0])-p[0])
                          + mp.mpf(r)**2*phi)
                length += mp.mpf(r)*abs(phi)
            else:
                _, _, interior, basis, weights = piece
                poles = [tuple(p), *[tuple(x) for x in interior], tuple(q)]
                sp = Spline(spl.BSpline2(basis, poles, weights or [1.0]*len(poles)))
                pieces.append(sp)
                t, l = sp.exact_terms()
                twice += t
                length += l
        out.append((Path(pieces), abs(twice)/2, length))
    return out


class PatchSurface:
    """A nonperiodic B-spline surface as its exact Bézier patches, evaluated
    in binary64: S = A/w with S_u = (A_u - w_u S)/w and likewise S_v."""
    def __init__(self, s):
        self.s = s
        nv = spl._pole_count(s.v)
        nu = len(s.poles)//nv
        h = spl.homogeneous(s.poles, s.weights)
        self.ubreaks = _breaks(s.u)
        self.vbreaks = _breaks(s.v)
        p, q = s.u.degree, s.v.degree
        self.patches = []
        for ua, ub in zip(self.ubreaks, self.ubreaks[1:]):
            rows = []
            for m in range(nv):
                rows.append(_bezier(s.u, [h[i*nv+m] for i in range(nu)], ua, ub))
            row = []
            for va, vb in zip(self.vbreaks, self.vbreaks[1:]):
                grid = [_bezier(s.v, [rows[m][i] for m in range(nv)], va, vb) for i in range(p+1)]
                row.append(_Patch(grid, float(ua), float(ub), float(va), float(vb)))
            self.patches.append(row)
        self.domain = (float(self.ubreaks[0]), float(self.ubreaks[-1]),
                       float(self.vbreaks[0]), float(self.vbreaks[-1]))
        self.starts = []
        for row in self.patches:
            for pt in row:
                for i in range(5):
                    for j in range(5):
                        u = pt.u0+(pt.u1-pt.u0)*i/4
                        v = pt.v0+(pt.v1-pt.v0)*j/4
                        self.starts.append((self.jet(u, v)[0], (u, v)))
        self.memo = {}

    def patch(self, u, v):
        i = min(max(_locate(self.ubreaks, u), 0), len(self.patches)-1)
        j = min(max(_locate(self.vbreaks, v), 0), len(self.patches[0])-1)
        return self.patches[i][j]

    def jet(self, u, v):
        return self.patch(u, v).jet(u, v)

    def clamp(self, u, v):
        a, b, c, d = self.domain
        return min(max(u, a), b), min(max(v, c), d)

    def project(self, x, box=None):
        """(distance, (u, v)) of the nearest point of the surface over the
        parameter box (the whole domain by default) to x."""
        key = (x, box)
        if key in self.memo:
            return self.memo[key]
        a, b, c, d = box or self.domain
        starts = [st for st in self.starts if a <= st[1][0] <= b and c <= st[1][1] <= d] or \
            [(self.jet(a, c)[0], (a, c))]
        _, (u, v) = min(starts, key=lambda st: math.dist(st[0], x))
        for _ in range(80):
            S, Su, Sv = self.jet(u, v)
            r = _vsub(S, x)
            g = (-_vdot(r, Su), -_vdot(r, Sv))
            a11, a12, a22 = _vdot(Su, Su), _vdot(Su, Sv), _vdot(Sv, Sv)
            det = a11*a22-a12*a12
            if det <= 0:
                break
            du, dv = (g[0]*a22-g[1]*a12)/det, (a11*g[1]-a12*g[0])/det
            nu, nv = u+du, v+dv
            # A step leaving the box: fix that coordinate at its bound and
            # move the other alone.
            if not a <= nu <= b:
                nu = min(max(nu, a), b)
                nv = min(max(v+g[1]/a22, c), d) if a22 > 0 else v
            elif not c <= nv <= d:
                nv = min(max(nv, c), d)
                nu = min(max(u+g[0]/a11, a), b) if a11 > 0 else u
            done = abs(nu-u) <= 1e-16*(b-a) and abs(nv-v) <= 1e-16*(d-c)
            u, v = nu, nv
            if done:
                break
        out = (math.dist(self.jet(u, v)[0], x), (u, v))
        if len(self.memo) > 4096:
            self.memo.clear()
        self.memo[key] = out
        return out

    def distance(self, x):
        return self.project(x)[0]


def _breaks(basis):
    """The distinct knots of a nonperiodic basis in its domain (Fractions)."""
    _, _, (a, e) = basis.flat()
    return [a]+[F(k) for k in basis.knots if a < F(k) < e]+[e]


def _locate(breaks, x):
    k = 0
    while k+2 < len(breaks) and x >= breaks[k+1]:
        k += 1
    return k


def _bezier(basis, hctrl, a, b):
    """The homogeneous Bézier controls on [a, b] (inside one span) of the
    spline of `basis` with controls `hctrl`, by blossoming (Fractions)."""
    flat, n, _ = basis.flat()
    p = basis.degree
    k = max(j for j in range(p, len(flat)-p-1) if flat[j] <= a < flat[j+1])
    ctrl = {j: hctrl[j % n if basis.periodic else j] for j in range(k-p, k+1)}
    out = []
    for i in range(p+1):
        args = [a]*(p-i)+[b]*i
        d = dict(ctrl)
        for r in range(1, p+1):
            for j in range(k, k-p+r-1, -1):
                al = (args[r-1]-flat[j])/(flat[j+p-r+1]-flat[j])
                d[j] = tuple((1-al)*x+al*y for x, y in zip(d[j-1], d[j]))
        out.append(d[k])
    return out


class _Patch:
    """A rational Bézier patch over [u0, u1] x [v0, v1]: grid[i][j]."""
    def __init__(self, grid, u0, u1, v0, v1):
        self.exact = grid
        self.g = [[tuple(float(c) for c in x) for x in row] for row in grid]
        self.u0, self.u1, self.v0, self.v1 = u0, u1, v0, v1
        p, q = len(self.g)-1, len(self.g[0])-1
        self.gu = [[tuple(p/(u1-u0)*(b-a) for a, b in zip(x, y)) for x, y in zip(r0, r1)]
                   for r0, r1 in zip(self.g, self.g[1:])]
        self.gv = [_diff(row, q/(v1-v0)) for row in self.g]

    def jet(self, u, v):
        s = (u-self.u0)/(self.u1-self.u0)
        t = (v-self.v0)/(self.v1-self.v0)
        a = _bern([_bern(row, t) for row in self.g], s)
        au = _bern([_bern(row, t) for row in self.gu], s)
        av = _bern([_bern(row, t) for row in self.gv], s)
        w = a[3]
        S = tuple(x/w for x in a[:3])
        Su = tuple((x-au[3]*y)/w for x, y in zip(au[:3], S))
        Sv = tuple((x-av[3]*y)/w for x, y in zip(av[:3], S))
        return S, Su, Sv

    def mp_columns(self, u):
        """The controls along v of the patch and of its u-derivative at a
        local u in [0, 1], in mpmath from the exact grid."""
        if not hasattr(self, 'mp_grid'):
            self.mp_grid = [[tuple(_mpf(c) for c in x) for x in row] for row in self.exact]
        g = self.mp_grid
        p = len(g)-1
        cols = list(zip(*g))
        at = [_bern(list(col), u) for col in cols]
        du = [_bern([tuple(p*(b-a) for a, b in zip(x, y)) for x, y in zip(col, col[1:])], u) for col in cols]
        return at, du


def _mp_surface_jet(at, du, v):
    """(S_u, S_v) at local v from `_Patch.mp_columns`."""
    q = len(at)-1
    a = _bern(at, v)
    au = _bern(du, v)
    av = _bern([tuple(q*(b-x_) for x_, b in zip(x, y)) for x, y in zip(at, at[1:])], v)
    w = a[3]
    S = tuple(x/w for x in a[:3])
    return (tuple((x-au[3]*y)/w for x, y in zip(au[:3], S)),
            tuple((x-av[3]*y)/w for x, y in zip(av[:3], S)))


def _mp_area(surface, box):
    """∫∫ |S_u x S_v| over a parameter box inside the domain: Gauss-Legendre
    on each patch's part, 2 x 2 pieces of 24 x 24 nodes (mpmath)."""
    from cell_reference import gauss
    a, b, c, d = (F(x) for x in box)
    total = mp.mpf(0)
    for row in surface.patches:
        for pt in row:
            u0, u1 = max(a, F(pt.u0)), min(b, F(pt.u1))
            v0, v1 = max(c, F(pt.v0)), min(d, F(pt.v1))
            if u0 >= u1 or v0 >= v1:
                continue
            lu, lv = F(pt.u1)-F(pt.u0), F(pt.v1)-F(pt.v0)
            s0, s1 = _mpf((u0-F(pt.u0))/lu), _mpf((u1-F(pt.u0))/lu)
            t0, t1 = _mpf((v0-F(pt.v0))/lv), _mpf((v1-F(pt.v0))/lv)

            def inner(s, pt=pt, t0=t0, t1=t1):
                at, du = pt.mp_columns(s)

                def f(t):
                    return norm_mp(cross(*_mp_surface_jet(at, du, t)))
                return sum(gauss(f, t0+(t1-t0)*k/2, t0+(t1-t0)*(k+1)/2) for k in range(2))
            # The partials are in the patch's local parameters, so this is
            # the area of the part.
            total += sum(gauss(inner, s0+(s1-s0)*k/2, s0+(s1-s0)*(k+1)/2) for k in range(2))
    return total


def norm_mp(v):
    return mp.sqrt(sum(x*x for x in v))


class Dome:
    """The box [0, w] x [0, d] x [0, h] with its top replaced by a spline
    graph over the same rectangle, equal to h on its boundary: the bottom,
    four plane walls and the spline top as patches."""
    slow = True

    def __init__(self, size, surface):
        self.w, self.d, self.h = size
        self.surface = surface
        w, d, h = self.w, self.d, self.h
        rect = lambda a, b, x0, x1, y0, y1: math.hypot(max(0.0, x0-a, a-x1), max(0.0, y0-b, b-y1))
        self.patches = [
            (lambda p: math.hypot(rect(p[0], p[1], 0, w, 0, d), p[2]), lambda p: abs(p[2])),
            (lambda p: math.hypot(rect(p[0], p[2], 0, w, 0, h), p[1]), lambda p: abs(p[1])),
            (lambda p: math.hypot(rect(p[1], p[2], 0, d, 0, h), p[0]-w), lambda p: abs(p[0]-w)),
            (lambda p: math.hypot(rect(p[0], p[2], 0, w, 0, h), p[1]-d), lambda p: abs(p[1]-d)),
            (lambda p: math.hypot(rect(p[1], p[2], 0, d, 0, h), p[0]), lambda p: abs(p[0])),
            (surface.distance, surface.distance),
        ]

    def distance(self, p, within=None):
        return nearest(self.patches, p, within)

    def inside(self, p):
        x, y, z = p
        if not (0 < x < self.w and 0 < y < self.d and z > 0):
            return False
        # S(u, v) = (u, v, z(u, v)): the top's height at (x, y).
        return z < self.surface.jet(x, y)[0][2]


class SplineSheet:
    """A face body on a spline surface's whole domain, less a parameter
    rectangle: one patch, trimmed through the projection's parameters."""
    slow = True

    def __init__(self, surface, hole):
        self.surface = surface
        self.hole = hole
        self.patches = [(self.trimmed, surface.distance)]

    def trimmed(self, p):
        d, (u, v) = self.surface.project(p)
        if self.hole is None:
            return d
        u0, u1, v0, v1 = self.hole
        if not (u0 < u < u1 and v0 < v < v1):
            return d
        # The nearest point of the face lies on the hole's boundary.
        return min(self.surface.project(p, box)[0] for box in
                   ((u0, u0, v0, v1), (u1, u1, v0, v1), (u0, u1, v0, v0), (u0, u1, v1, v1)))

    def distance(self, p, within=None):
        return self.trimmed(p)

    def inside(self, p):
        return False


def moved_axes(motion):
    """(origin, x, y, n) of the world frame after a rigid motion (unit axis,
    angle, shift), as generate_tessellation_fixtures moves its models."""
    from brep_reference import cos_rn, sin_rn
    axis, angle, shift = motion
    k = [a/math.sqrt(sum(b*b for b in axis)) for a in axis]
    c, s = cos_rn(angle), sin_rn(angle)

    def rot(v):
        kv = sum(a*b for a, b in zip(k, v))
        kx = (k[1]*v[2]-k[2]*v[1], k[2]*v[0]-k[0]*v[2], k[0]*v[1]-k[1]*v[0])
        return tuple(v[i]*c+kx[i]*s+k[i]*kv*(1-c) for i in range(3))
    return (tuple(shift), rot((1.0, 0.0, 0.0)), rot((0.0, 1.0, 0.0)), rot((0.0, 0.0, 1.0)))


def _hull_extent(points):
    lo = [min(p[k] for p in points) for k in range(3)]
    hi = [max(p[k] for p in points) for k in range(3)]
    return max(h-l for h, l in zip(hi, lo)), max(max(abs(x) for x in lo), max(abs(x) for x in hi))


def spline_body(case):
    """The independent body of a T-b spline case."""
    if case.kind in ('prism', 'face'):
        boundaries = case.boundaries if case.kind == 'prism' else [[('ring', case.ring)]]
        paths = spline_paths(boundaries)
        area = paths[0][1]-sum(a for _, a, _ in paths[1:])
        length = sum(l for _, _, l in paths)
        (x0, y0), (x1, y1) = paths[0][0].bounds()
        holes = len(paths)-1
        fr = moved_axes(case.motion) if case.motion else ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                                          (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        extent = max(abs(v) for v in fr[0])+max(abs(x0), abs(x1), abs(y0), abs(y1))
        if case.kind == 'face':
            scale = max(x1-x0, y1-y0)
            shape = Prism(fr, [p for p, _, _ in paths], 0.0, 0.0, sheet=True)
            shape.slow = True
            return Body(shape, False, 1-holes, 1+holes, area, mp.mpf(0), scale, scale+extent)
        h = case.z1-case.z0
        scale = max(x1-x0, y1-y0, h)
        shape = Prism(fr, [p for p, _, _ in paths], case.z0, case.z1)
        shape.slow = True
        extent += max(abs(case.z0), abs(case.z1))
        return Body(shape, True, 2-2*holes, 0, 2*area+length*h, area*h, scale, scale+extent)
    surface = PatchSurface(case.surface)
    scale, extent = _hull_extent(case.surface.poles)
    if case.kind == 'dome':
        w, d, h = case.size
        # The top's volume above h: z - h integrates over the patches'
        # parameters (x = u, y = v) as the mean of its controls (nonrational
        # patches) times each patch's area.
        above = mp.mpf(0)
        for row in surface.patches:
            for pt in row:
                zs = [F(x[2])/F(x[3]) for r in pt.exact for x in r]
                assert all(F(x[3]) == F(pt.exact[0][0][3]) for r in pt.exact for x in r)
                above += _mpf((sum(zs)/len(zs)-F(h))*(F(pt.u1)-F(pt.u0))*(F(pt.v1)-F(pt.v0)))
        volume = mp.mpf(w)*d*h+above
        area = mp.mpf(w)*d+2*(mp.mpf(w)+d)*h+_mp_area(surface, surface.domain)
        return Body(Dome(case.size, surface), True, 2, 0, area, volume, max(scale, w, d), scale+extent)
    area = _mp_area(surface, surface.domain)
    if case.hole is not None:
        area -= _mp_area(surface, case.hole)
    holes = 0 if case.hole is None else 1
    return Body(SplineSheet(surface, case.hole), False, 1-holes, 1+holes, area, mp.mpf(0), scale,
                scale+extent)
