#!/usr/bin/env python3
"""Fixtures for S9a of REVIEW_NOTES.md: Booleans of two prisms in one frame.

`boolean-cases.txt` lists each case in the protocol of `identity-cases.txt`
joined by a `boolean OP ID` row (`identity_reference.encode_boolean_case`):
the object's prism (op 91), then the tool's (op 92), the Boolean's
operation id 93. Both frames have bitwise-equal axes (`XY`, or the tilted
`TILT` of the split fixtures) and the tool's origin is the object's moved by
`a x + b y + c n` with `(a, b, c)` binary64 and the world offset exact in
binary64 (`shifted`), so the tool is the object's frame's prism of its
profile translated by `(a, b)` over its offsets moved by `c`.

`boolean-expected.tsv` gives per case, from `boolean_reference.py`:

* `expect KIND STEP`: the kernel's expected outcome, declared with the case
  from the recorded decisions (`prisms`: every result solid a prism of the
  object's frame, one run of slabs; `stack`: a general body, `OutOfDomain`
  until S9a.2; `empty`: no solids; `degenerate`: the result touches itself
  at a point or along an edge, `Degenerate` until the kernel holds
  non-manifold bodies), and the sub-step the reference computes (`S9a.1`:
  empty or one run of slabs; `S9a.2`: a stack), checked against the kind.
* `result N volume area cx cy cz`: the totals over the result's N solids
  (maximal connected regions of the regularized Boolean; regions meeting at
  a point or along an edge are separate), their volume, surface area and
  common centre (world coordinates), or `empty` when the result has no
  volume. Totals, not per solid: every solid's measures sum to them, and the
  solid count is checked apart.
* `slab w0 w1 area perimeter` per run of slabs with one region (heights in
  the object's frame, bottom to top; runs of equal regions merged, empty
  slabs dropped): the region's area and boundary length.

Touching conventions (the regularized Boolean): prisms sharing a wall fuse
into one solid, their common is empty and a cut keeps the object whole; a
common or cut leaving only a face, an edge or a point is empty; caps at
equal heights merge (the M3 fuse). A fuse of prisms touching along a
vertical edge (a vertex or a tangency in the profiles) is two solids sharing
the edge, and a cut leaving a hole tangent to the outer boundary a solid
touching itself along an edge: both `degenerate`.

Before writing, `reference_checks` compares the reference with independent
computations: each profile's area and first moments by Green's theorem in
closed form against the slicing's atoms, and its perimeter against its
classified pieces (and the shared pieces seen from either side); the three
operations' volumes and moments against `fuse = A + B - common` and `cut = A
- common` with `A` and `B` from Green's theorem; every pair of polygons'
common atom against exact Fraction clipping (Sutherland-Hodgman of each
object boundary by each triangle of an ear-clipping of each tool boundary,
holes subtracted); every pair of single circles against the closed-form
lens (area, centre and arc lengths), and hand-computed results. Each case's
solid count is declared with it and must equal the reference's.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path

from identity_reference import Boundary, Case, encode_boolean_case, stored_points
from curve_surface_reference import stored_axes
import boolean_reference as ref

ROOT = Path(__file__).resolve().parents[1]
XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (1.0, -2.0, 0.5, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)
BOOLEAN_OPERATION = 93


def shifted(frame, a, b, c):
    """The frame with its origin moved by `a x + b y + c n` (stored axes),
    required exact in binary64."""
    o, x, y, n = stored_axes(frame)
    out = []
    for i in range(3):
        exact = F(o[i])+F(a)*F(x[i])+F(b)*F(y[i])+F(c)*F(n[i])
        assert F(float(exact)) == exact, 'the offset origin is not binary64'
        out.append(float(exact))
    return tuple(out)+tuple(frame[3:])


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def disc(cx, cy, r):
    return Boundary(circle=(cx, cy, r))


def path(points, segments):
    return Boundary(points=points, segments=segments)


def arc(cx, cy, r, ccw=True):
    return (cx, cy, r, ccw)


def prism(boundaries, frame=XY, start=0.0, end=5.0, op=91):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Boolean:
    def __init__(self, name, operation, obj, tool, offset=(0.0, 0.0, 0.0), solids=1, kind='prisms'):
        self.name, self.operation, self.solids, self.kind = name, operation, solids, kind
        self.obj = prism(obj[0], obj[1], *obj[2], op=91)
        self.obj.name = name
        frame = shifted(obj[1], *offset)
        self.tool = prism(tool[0], frame, *tool[1], op=92)
        self.tool.name = name
        self.offset = offset

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def cases():
    """Every class of the S9a decisions, in both frames (`TILT` cases end in
    `_tilted`): overlapping rectangles, a rectangle and a circle or a
    stadium, lenses, one inside the other, identical, disjoint, touching
    along an edge (whole and partly) and at a vertex, collinear overlapping
    edges, arcs of one circle, tangent circles (outside and inside), a
    vertex on the other's edge, holes (cut into a solid, a tool through a
    hole, a hole filled, a common with a hole), equal and touching cap
    heights, a tool missing, spanning and partly covering the object's
    heights (S9a.2), and origins offset along the axis."""
    H = (0.0, 5.0)
    sq = [square(0.0, 0.0, 10.0, 10.0)]
    holed = [square(0.0, 0.0, 10.0, 10.0), square(3.0, 3.0, 7.0, 7.0)]
    rounded = [path([(1.0, 0.0), (7.0, 0.0), (8.0, 1.0), (8.0, 5.0), (7.0, 6.0), (1.0, 6.0), (0.0, 5.0),
                     (0.0, 1.0)],
                    [None, arc(7.0, 1.0, 1.0), None, arc(7.0, 5.0, 1.0), None, arc(1.0, 5.0, 1.0), None,
                     arc(1.0, 1.0, 1.0)])]
    stadium = [path([(-3.0, -1.0), (3.0, -1.0), (3.0, 1.0), (-3.0, 1.0)],
                    [None, arc(3.0, 0.0, 1.0), None, arc(-3.0, 0.0, 1.0)])]
    upper_half = [path([(-4.0, 0.0), (4.0, 0.0)], [None, arc(0.0, 0.0, 4.0)])]
    right_half = [path([(0.0, 4.0), (0.0, -4.0)], [None, arc(0.0, 0.0, 4.0)])]
    wedge = [Boundary(points=[(5.0, 10.0), (12.0, 4.0), (12.0, 12.0)])]
    B = Boolean
    return [
        # Overlapping rectangles.
        B('rects_fuse', 'fuse', (sq, XY, H), (sq, H), (5.0, 5.0, 0.0)),
        B('rects_cut_tilted', 'cut', (sq, TILT, H), (sq, H), (5.0, 4.0, 0.0)),
        B('bar_cut_two_pieces', 'cut', (sq, XY, H), ([square(0.0, 0.0, 6.0, 14.0)], (-1.0, 6.0)),
          (2.0, -2.0, 0.0), solids=2),
        # A rectangle and a circle or a stadium.
        B('rect_disc_cut', 'cut', (sq, XY, H), ([disc(0.0, 0.0, 3.0)], (-1.0, 6.0)), (10.0, 5.0, 0.0)),
        B('rect_stadium_common_tilted', 'common', (sq, TILT, H), (stadium, (-2.0, 7.0)), (2.0, 8.0, 0.0)),
        # Two circles: lenses.
        B('lens_common', 'common', ([disc(0.0, 0.0, 4.0)], XY, H), ([disc(0.0, 0.0, 4.0)], H), (5.0, 0.0, 0.0)),
        B('lens_fuse_tilted', 'fuse', ([disc(0.0, 0.0, 4.0)], TILT, H), ([disc(0.0, 0.0, 4.0)], H),
          (0.0, 4.0, 0.0)),
        B('lens_cut', 'cut', ([disc(0.0, 0.0, 4.0)], XY, H), ([disc(0.0, 0.0, 3.0)], (-1.0, 6.0)),
          (3.0, 4.0, 0.0)),
        # One inside the other.
        B('disc_hole_cut', 'cut', (sq, XY, H), ([disc(5.0, 5.0, 2.0)], (-1.0, 6.0))),
        B('disc_inside_common', 'common', (sq, XY, H), ([disc(5.0, 5.0, 2.0)], (-1.0, 6.0))),
        # Identical profiles.
        B('identical_cut', 'cut', (rounded, XY, H), (rounded, H), solids=0, kind='empty'),
        B('identical_common_tilted', 'common', (rounded, TILT, H), (rounded, H)),
        # Disjoint.
        B('disjoint_fuse', 'fuse', (sq, XY, H), ([disc(0.0, 0.0, 2.0)], H), (20.0, 5.0, 0.0), solids=2),
        B('disjoint_common', 'common', (sq, XY, H), ([disc(0.0, 0.0, 2.0)], H), (20.0, 5.0, 0.0),
          solids=0, kind='empty'),
        # Touching along an edge: a shared wall.
        B('edge_touch_fuse', 'fuse', (sq, XY, H), (sq, H), (10.0, 0.0, 0.0)),
        B('edge_touch_common', 'common', (sq, XY, H), (sq, H), (10.0, 0.0, 0.0), solids=0, kind='empty'),
        B('edge_touch_cut_tilted', 'cut', (sq, TILT, H), (sq, H), (10.0, 0.0, 0.0)),
        B('edge_partial_fuse', 'fuse', (sq, XY, H), ([square(0.0, 0.0, 6.0, 6.0)], H), (10.0, 2.0, 0.0)),
        # Touching at a vertex.
        B('vertex_touch_fuse', 'fuse', (sq, XY, H), (sq, H), (10.0, 10.0, 0.0), solids=2, kind='degenerate'),
        # Collinear overlapping edges: coincident walls.
        B('collinear_fuse', 'fuse', (sq, XY, H), ([square(0.0, 0.0, 10.0, 6.0)], H), (5.0, 0.0, 0.0)),
        B('collinear_cut', 'cut', (sq, XY, H), ([square(0.0, 0.0, 10.0, 6.0)], H), (5.0, 0.0, 0.0)),
        B('collinear_common_tilted', 'common', (sq, TILT, H), ([square(0.0, 0.0, 10.0, 6.0)], H),
          (5.0, 0.0, 0.0)),
        # Arcs of one circle overlapping.
        B('same_circle_common', 'common', (upper_half, XY, H), (right_half, H)),
        B('same_circle_cut_tilted', 'cut', (upper_half, TILT, H), (right_half, H)),
        # Tangent circles, outside and inside.
        B('tangent_out_fuse', 'fuse', ([disc(0.0, 0.0, 2.0)], XY, H), ([disc(0.0, 0.0, 2.0)], H),
          (4.0, 0.0, 0.0), solids=2, kind='degenerate'),
        B('tangent_out_common', 'common', ([disc(0.0, 0.0, 2.0)], XY, H), ([disc(0.0, 0.0, 2.0)], H),
          (4.0, 0.0, 0.0), solids=0, kind='empty'),
        B('tangent_in_cut', 'cut', ([disc(0.0, 0.0, 4.0)], XY, H), ([disc(0.0, 0.0, 2.0)], H),
          (2.0, 0.0, 0.0), kind='degenerate'),
        B('tangent_in_common', 'common', ([disc(0.0, 0.0, 4.0)], XY, H), ([disc(0.0, 0.0, 2.0)], H),
          (2.0, 0.0, 0.0)),
        B('tangent_in_fuse_tilted', 'fuse', ([disc(0.0, 0.0, 4.0)], TILT, H), ([disc(0.0, 0.0, 2.0)], H),
          (2.0, 0.0, 0.0)),
        # A vertex on the other's edge.
        B('vertex_on_edge_cut', 'cut', (sq, XY, H), (wedge, H)),
        B('vertex_on_edge_fuse_tilted', 'fuse', (sq, TILT, H), (wedge, H)),
        # Holes.
        B('holed_pin_cut', 'cut', (holed, XY, H), ([disc(5.0, 5.0, 1.0)], (-1.0, 6.0))),
        B('holed_fill_fuse', 'fuse', (holed, XY, H), ([square(3.0, 3.0, 7.0, 7.0)], H)),
        B('holed_disc_common_tilted', 'common', (holed, TILT, H), ([disc(5.0, 5.0, 4.0)], H)),
        # Heights: the tool missing, spanning and partly covering the
        # object's; caps at equal heights; origins offset along the axis.
        B('tool_misses_heights_cut', 'cut', (sq, XY, H), ([disc(5.0, 5.0, 2.0)], (6.0, 8.0))),
        B('common_offset_heights', 'common', (sq, XY, H), (sq, H), (5.0, 5.0, 3.0)),
        B('stacked_fuse_tilted', 'fuse', (sq, TILT, H), (sq, H), (0.0, 0.0, 4.0)),
        B('cap_touch_fuse', 'fuse', (sq, XY, H), (sq, H), (0.0, 0.0, 5.0)),
        B('cap_touch_common', 'common', (sq, XY, H), (sq, H), (0.0, 0.0, 5.0), solids=0, kind='empty'),
        B('step_fuse', 'fuse', (sq, XY, H), (sq, (0.0, 8.0)), (5.0, 0.0, 0.0), kind='stack'),
        B('pocket_cut', 'cut', (sq, XY, H), ([disc(5.0, 5.0, 2.0)], (3.0, 6.0)), kind='stack'),
        B('pocket_cut_tilted', 'cut', (rounded, TILT, H), ([disc(0.0, 0.0, 1.5)], (2.0, 6.0)), (4.0, 2.0, 0.0),
          kind='stack'),
        B('slot_cut_two_solids', 'cut', (sq, XY, H), ([square(-1.0, -1.0, 11.0, 11.0)], (2.0, 3.0)),
          solids=2, kind='stack'),
        B('cap_touch_stack_fuse', 'fuse', (sq, XY, H), ([disc(5.0, 5.0, 2.0)], (5.0, 8.0)), kind='stack'),
        B('partial_heights_cut_tilted', 'cut', (sq, TILT, H), (sq, (2.0, 8.0)), (5.0, 4.0, 0.0), kind='stack'),
    ]


# ------------------------------------------------------------------ checks

def green(elements):
    """Area and first moments of a profile by Green's theorem in closed
    form: lines exactly, arcs by their angle integrals."""
    M, mp = ref.M, ref.mp
    area = mu = mv = M(0)
    for e in elements:
        if e.kind == 'L':
            (x0, y0), (x1, y1) = e.p, e.q
            area += M(x0*y1-x1*y0)/2
            mu += M((y1-y0)*(x0*x0+x0*x1+x1*x1))/6
            mv -= M((x1-x0)*(y0*y0+y0*y1+y1*y1))/6
            continue
        cx, cy, r = M(e.c[0]), M(e.c[1]), M(e.r)
        a, b = e.theta0, e.theta0+e.sweep
        S = lambda t: mp.sin(t)
        C = lambda t: mp.cos(t)
        d = lambda f: f(b)-f(a)
        # x dy - y dx = (r^2 + r cx cos + r cy sin) dt
        area += (r*r*(b-a)+r*cx*d(S)-r*cy*d(C))/2
        # x^2 dy / 2 = (cx + r cos)^2 r cos dt / 2
        c1, c2, c3 = d(S), (b-a)/2+d(lambda t: mp.sin(2*t))/4, d(lambda t: mp.sin(t)-mp.sin(t)**3/3)
        mu += r*(cx*cx*c1+2*cx*r*c2+r*r*c3)/2
        # -y^2 dx / 2 = (cy + r sin)^2 r sin dt / 2
        s1, s2, s3 = -d(C), (b-a)/2-d(lambda t: mp.sin(2*t))/4, d(lambda t: -mp.cos(t)+mp.cos(t)**3/3)
        mv += r*(cy*cy*s1+2*cy*r*s2+r*r*s3)/2
    return area, mu, mv


def _left(a, b, p):
    return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])


def triangulate(points):
    """Ear clipping of a simple counter-clockwise polygon in Fractions."""
    pts = list(points)
    out = []
    while len(pts) > 3:
        n = len(pts)
        for i in range(n):
            a, b, c = pts[i-1], pts[i], pts[(i+1) % n]
            if _left(a, b, c) <= 0:
                continue
            if any(_left(a, b, p) >= 0 and _left(b, c, p) >= 0 and _left(c, a, p) >= 0
                   for p in pts if p not in (a, b, c)):
                continue
            out.append((a, b, c))
            del pts[i]
            break
        else:
            raise ValueError('no ear')
    out.append(tuple(pts))
    return out


def clip(subject, triangle):
    """Sutherland-Hodgman of a polygon by a counter-clockwise triangle."""
    out = list(subject)
    for i in range(3):
        a, b = triangle[i], triangle[(i+1) % 3]
        pts, out = out, []
        for j, p in enumerate(pts):
            q = pts[(j+1) % len(pts)]
            fp, fq = _left(a, b, p), _left(a, b, q)
            if fp >= 0:
                out.append(p)
            if fp*fq < 0:
                t = fp/(fp-fq)
                out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
        if not out:
            return []
    return out


def shoelace(pts):
    area = mx = my = F(0)
    for i, p in enumerate(pts):
        q = pts[(i+1) % len(pts)]
        c = p[0]*q[1]-q[0]*p[1]
        area += c/2
        mx += (p[0]+q[0])*c/6
        my += (p[1]+q[1])*c/6
    return area, mx, my


def exact_common(case, pair):
    """The common atom of two polygon profiles in exact Fractions."""
    a, b, _ = pair.offset
    total = [F(0)]*3
    for i, bo in enumerate(case.obj.boundaries):
        subject = [(F(x), F(y)) for x, y in stored_points(bo, case.obj.tolerance)]
        for j, bt in enumerate(case.tool.boundaries):
            sign = (1 if i == 0 else -1)*(1 if j == 0 else -1)
            tool = [(F(x)+a, F(y)+b) for x, y in stored_points(bt, case.tool.tolerance)]
            for tri in triangulate(tool):
                piece = clip(subject, tri)
                if piece:
                    for k, v in enumerate(shoelace(piece)):
                        total[k] += sign*v
    return total


def exact_classes(X, Y):
    """The boundary classes of polygon X's edges against polygon Y in exact
    Fractions: each edge cut where Y's edges cross it or Y's vertices lie
    on it (rational parameters), each piece classified at its midpoint (on
    an edge of Y with the same or the opposite direction, else inside or
    outside by the half-open crossing rule). Returns {class: [squared
    length, parameter span]} summed as `sqrt(squared) * span` per edge."""
    out = {}
    for e in X:
        ts = {F(0), F(1)}
        dd = e.d[0]**2+e.d[1]**2
        for f in Y:
            den = _cross2(e.d, f.d)
            w = (f.p[0]-e.p[0], f.p[1]-e.p[1])
            if den != 0:
                t, s = _cross2(w, f.d)/den, _cross2(w, e.d)/den
                if 0 <= t <= 1 and 0 <= s <= 1:
                    ts.add(t)
            if _cross2(w, e.d) == 0:
                t = (w[0]*e.d[0]+w[1]*e.d[1])/dd
                if 0 <= t <= 1:
                    ts.add(t)
        ts = sorted(ts)
        for t0, t1 in zip(ts, ts[1:]):
            tm = (t0+t1)/2
            m = (e.p[0]+tm*e.d[0], e.p[1]+tm*e.d[1])
            cls = None
            for f in Y:
                w = (m[0]-f.p[0], m[1]-f.p[1])
                if _cross2(w, f.d) == 0 and 0 <= (w[0]*f.d[0]+w[1]*f.d[1]) <= f.d[0]**2+f.d[1]**2:
                    cls = 'on_same' if e.d[0]*f.d[0]+e.d[1]*f.d[1] > 0 else 'on_opposite'
            if cls is None:
                inside = False
                for f in Y:
                    (x0, y0), (x1, y1) = f.p, f.q
                    if (y0 > m[1]) != (y1 > m[1]) and m[0] < x0+(m[1]-y0)*(x1-x0)/(y1-y0):
                        inside = not inside
                cls = 'inside' if inside else 'outside'
            out.setdefault(cls, []).append((dd, t1-t0))
    return out


def _cross2(a, b):
    return a[0]*b[1]-a[1]*b[0]


def lens(c1, r1, c2, r2):
    """Area, centre and the arc lengths of each circle inside the other
    for two circles (closed forms)."""
    M, mp = ref.M, ref.mp
    d = mp.sqrt((c2[0]-c1[0])**2+(c2[1]-c1[1])**2)
    pi = mp.pi
    if d >= r1+r2:
        return M(0), None, M(0), M(0)
    if d <= abs(r1-r2):
        small, c = (r1, c1) if r1 <= r2 else (r2, c2)
        return pi*small**2, c, (2*pi*r1 if r1 <= r2 else M(0)), (2*pi*r2 if r2 < r1 else M(0))
    a1 = mp.acos((d*d+r1*r1-r2*r2)/(2*d*r1))
    a2 = mp.acos((d*d+r2*r2-r1*r1)/(2*d*r2))
    seg1 = r1*r1*(a1-mp.sin(a1)*mp.cos(a1))
    seg2 = r2*r2*(a2-mp.sin(a2)*mp.cos(a2))
    arm1 = M(2)/3*r1**3*mp.sin(a1)**3/seg1
    arm2 = M(2)/3*r2**3*mp.sin(a2)**3/seg2
    e = ((c2[0]-c1[0])/d, (c2[1]-c1[1])/d)
    p1 = tuple(c1[i]+arm1*e[i] for i in range(2))
    p2 = tuple(c2[i]-arm2*e[i] for i in range(2))
    area = seg1+seg2
    centre = tuple((seg1*p1[i]+seg2*p2[i])/area for i in range(2))
    return area, centre, 2*r1*a1, 2*r2*a2


def reference_checks():
    """Returns the largest relative deviations by check; raises beyond
    1e-30."""
    M = ref.M
    worst, covered = {}, {}

    def near(kind, got, want):
        dev = abs(M(got)-M(want))/max(1, abs(M(want)))
        worst[kind] = max(worst.get(kind, M(0)), dev)
    for case in cases():
        pair = ref.Pair(case.obj, case.tool)
        # Green's theorem against the slicing's atoms.
        for name, elements, atoms in (('A', pair.A, {'AB', 'A'}), ('B', pair.B, {'AB', 'B'})):
            for got, want in zip(pair.measure(atoms), green(elements)):
                near('green', got, want)
            whole = sum((e.length() for e in elements), M(0))
            parts = sum((pair.length(name, c) for c in ('outside', 'inside', 'on_same', 'on_opposite')), M(0))
            near('perimeter', parts, whole)
        for c in ('on_same', 'on_opposite'):
            near('perimeter', pair.length('B', c), pair.length('A', c))
        # fuse = A + B - common, cut = A - common (volumes and moments).
        gA, gB = green(pair.A), green(pair.B)
        (a0, a1), (b0, b1) = pair.heights['A'], pair.heights['B']
        hA, hB = M(a1-a0), M(b1-b0)
        mA = (gA[0]*hA, gA[1]*hA, gA[2]*hA, gA[0]*hA*M(a0+a1)/2)
        mB = (gB[0]*hB, gB[1]*hB, gB[2]*hB, gB[0]*hB*M(b0+b1)/2)

        def moments(operation):
            out = [M(0)]*4
            for w0, w1, atoms in pair.slabs(operation):
                h = M(w1-w0)
                area, su, sv = pair.measure(atoms)
                for k, v in enumerate((h*area, h*su, h*sv, h*area*M(w0+w1)/2)):
                    out[k] += v
            return out
        fuse, cut, common = moments('fuse'), moments('cut'), moments('common')
        for k in range(4):
            near('identities', fuse[k], mA[k]+mB[k]-common[k])
            near('identities', cut[k], mA[k]-common[k])
        # Polygons against exact clipping.
        if all(e.kind == 'L' for e in pair.A+pair.B):
            for got, want in zip(pair.measure({'AB'}), exact_common(case, pair)):
                near('exact_clipping', got, want)
            covered['exact_clipping'] = covered.get('exact_clipping', 0)+1
            for name, X, Y in (('A', pair.A, pair.B), ('B', pair.B, pair.A)):
                exact = exact_classes(X, Y)
                for cls in ('outside', 'inside', 'on_same', 'on_opposite'):
                    want = sum((ref.mp.sqrt(M(dd))*M(span) for dd, span in exact.get(cls, [])), M(0))
                    near('exact_classes', pair.length(name, cls), want)
        # Single circles against the lens.
        if len(pair.A) == 1 and len(pair.B) == 1 and pair.A[0].kind == pair.B[0].kind == 'A' \
                and pair.A[0].full and pair.B[0].full:
            e, f = pair.A[0], pair.B[0]
            area, centre, inA, inB = lens((M(e.c[0]), M(e.c[1])), M(e.r), (M(f.c[0]), M(f.c[1])), M(f.r))
            got = pair.measure({'AB'})
            near('lens', got[0], area)
            if centre is not None:
                near('lens', got[1], area*centre[0])
                near('lens', got[2], area*centre[1])
            near('lens', pair.length('A', 'inside')+pair.length('A', 'on_same'), inA)
            near('lens', pair.length('B', 'inside'), inB)
            covered['lens'] = covered.get('lens', 0)+1
        # The declared solids and kinds.
        solids = pair.result(case.operation)[0]
        assert solids == case.solids, (case.name, solids, case.solids)
        step = ref.step(pair, case.operation)
        assert (case.kind == 'empty') == (solids == 0), case.name
        assert case.kind != 'prisms' or step == 'S9a.1', case.name
        assert case.kind != 'stack' or step == 'S9a.2', case.name
    # Hand-computed results.
    by_name = {c.name: c for c in cases()}
    mp = ref.mp
    hand = {
        # [0,10]^2 u [5,15]^2, height 5: area 175, outline 60.
        'rects_fuse': (875, 2*175+60*5, (M(15)/2, M(15)/2, M(5)/2)),
        # [0,10]^2 \ [2,8] x R: two 2 x 10 bars.
        'bar_cut_two_pieces': (200, 2*(2*20+24*5), (5, 5, M(5)/2)),
        # Two unit-height-5 squares side by side: a 20 x 10 prism.
        'edge_touch_fuse': (1000, 2*200+60*5, (10, 5, M(5)/2)),
        # The square with a cylindrical hole of radius 2 through it.
        'disc_hole_cut': (5*(100-4*mp.pi), 2*(100-4*mp.pi)+5*(40+4*mp.pi), (5, 5, M(5)/2)),
        # A quarter disc of radius 4.
        'same_circle_common': (5*4*mp.pi, 2*4*mp.pi+5*(8+2*mp.pi), (16/(3*mp.pi), 16/(3*mp.pi), M(5)/2)),
        # Stacked: a 10 x 10 prism of height 10 (touching caps merged).
        'cap_touch_fuse': (1000, 200+400, (5, 5, 5)),
        # A step: [0,10]^2 x [0,5] u [5,15] x [0,10] x [0,8].
        'step_fuse': (5*150+3*100, 150+100+50+5*50+3*40, None),
    }
    covered['hand'] = len(hand)
    for name, (V, A, centre) in hand.items():
        case = by_name[name]
        solids, gotV, gotA, gotC, _ = ref.Pair(case.obj, case.tool).result(case.operation)
        near('hand', gotV, V)
        near('hand', gotA, A)
        if centre is not None:
            for g, w in zip(gotC, centre):
                near('hand', g, w)
    if max(worst.values()) > M(10)**-30:
        raise SystemExit(f'boolean reference checks failed: {worst}')
    return worst, covered


def generate():
    blocks = []
    out = ['# case\trow (S9a, boolean_reference.py: expect KIND STEP, then result N volume area cx cy cz '
           'or empty, then slab w0 w1 area perimeter per run of slabs)']
    for case in cases():
        blocks.append(case.encode())
        rows, pair = ref.rows(case.obj, case.operation, case.tool)
        out.append(f'{case.name}\texpect {case.kind} {ref.step(pair, case.operation)}')
        for row in rows:
            out.append(f'{case.name}\t{row}')
    return {'boolean-cases.txt': '\n'.join(blocks)+'\n', 'boolean-expected.tsv': '\n'.join(out)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    names = [c.name for c in cases()]
    assert len(names) == len(set(names)), 'duplicate case names'
    worst, covered = reference_checks()
    files = generate()
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    listed = cases()
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')))
    print('reference checks (largest relative deviation):',
          ', '.join(f'{k} {ref.mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases compared with exact clipping and classes:', covered['exact_clipping'],
          '- with closed-form lenses:', covered['lens'], '- with hand-computed results:', covered['hand'])


if __name__ == '__main__':
    main()
