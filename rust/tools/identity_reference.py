"""Independent reference for value identity and history (IDENTITY_AND_HISTORY.md).

This module implements, without Rust or OCCT, the version-1 derivation byte
encoding, the 128-bit FNV-1a digest that turns a derivation into an entity id,
the structural enumeration of every entity an extrusion creates, and the
relations each operation must report. Fixture generators and the native
comparison import it; nothing here reads a Rust result.

Encoding (all integers little-endian):
    b"RSID", version u8 = 1, operation u64, operation kind u8, entity kind u8,
    role u8, ordinal u32, parent count u32, then each parent:
        Label:   tag 1, u64 label
        Profile: tag 2, boundary u32, element u8 (0 boundary, 1 segment,
                 2 vertex), index u32
        Entity:  tag 3, 16 id bytes
The id is the FNV-1a-128 digest of those bytes, stored big-endian.

Case protocol (`encode_case`): a block of rows `case NAME TOL`, `op ID`, a
construction (`frame` and `offsets` or `make`, `box`, `cone`, `sphere` or
`torus`), `boundary` rows, `transform` rows and `end`; a split case (S8) adds
a `split` row before `end`, and a Boolean case (S9a, `encode_boolean_case`)
is two prisms' rows joined by a `boolean OP ID` row.
"""
from dataclasses import dataclass, field
from fractions import Fraction as F
import struct

FNV_OFFSET = 0x6c62272e07bb014262b821756295c58d
FNV_PRIME = 0x0000000001000000000000000000013B
MASK = (1 << 128)-1

KIND = {'extrude': 1, 'transform': 2, 'external': 3, 'composite': 4, 'height_split': 5,
        'stacked_fuse': 6, 'revolve': 7, 'make_face': 8, 'make_wire': 9, 'plane_split': 10,
        'fuse': 11, 'cut': 12, 'common': 13}
ENTITY = {'vertex': 1, 'edge': 2, 'face': 3, 'body': 4, 'region': 5}
DIMENSION = {'vertex': 0, 'edge': 1, 'face': 2, 'body': 3, 'region': 3}
ROLE = {'start_cap': 1, 'end_cap': 2, 'wall': 3, 'bottom_edge': 4, 'top_edge': 5,
        'vertical': 6, 'seam': 7, 'bottom_vertex': 8, 'top_vertex': 9,
        'seam_vertex': 10, 'body': 11, 'external': 12, 'region': 13, 'cut_face': 14,
        'cut_edge': 15, 'cut_vertex': 16, 'apex': 17, 'pole': 18, 'face': 19, 'edge': 20,
        'vertex': 21}
ELEMENT = {'boundary': 0, 'segment': 1, 'vertex': 2}
RELATION = {'unchanged': 1, 'modified': 2, 'generated': 3, 'split': 4, 'merged': 5,
            'deleted': 6}


def fnv128(data):
    h = FNV_OFFSET
    for byte in data:
        h ^= byte
        h = (h*FNV_PRIME) & MASK
    return h.to_bytes(16, 'big')


# Parents are tuples: ('label', n), ('profile', boundary, element, index),
# ('entity', id_bytes).
def encode_parent(p):
    if p[0] == 'label':
        return struct.pack('<BQ', 1, p[1])
    if p[0] == 'profile':
        return struct.pack('<BIBI', 2, p[1], ELEMENT[p[2]], p[3])
    assert p[0] == 'entity' and len(p[1]) == 16
    return struct.pack('<B', 3)+p[1]


def encode_parents(parents):
    return struct.pack('<I', len(parents))+b''.join(encode_parent(p) for p in parents)


@dataclass(frozen=True)
class Derivation:
    operation: int
    kind: str
    entity: str
    role: str
    ordinal: int
    parents: tuple

    def encode(self):
        return (b'RSID'+struct.pack('<BQBBBI', 1, self.operation, KIND[self.kind],
                                    ENTITY[self.entity], ROLE[self.role], self.ordinal)
                + encode_parents(self.parents))

    def id(self):
        return fnv128(self.encode())


def hexid(i):
    return i.hex()


def parent_text(p):
    if p[0] == 'label':
        return f'L{p[1]}'
    if p[0] == 'profile':
        return f'P{p[1]}.{p[2][0]}{p[3]}'
    return 'E'+p[1].hex()


# ------------------------------------------------------------------ profiles

@dataclass
class Boundary:
    """Caller input: a polygon (points), a circle, or a path (points with a
    segment per point, S5), with optional labels in the caller's input
    order."""
    points: list = None          # [(x, y)] for a polygon or a path
    circle: tuple = None         # (cx, cy, r)
    labels: tuple = None         # (boundary, [segment...], [vertex...])
    # For a path, segment j from point j to j+1: None for a line,
    # (cx, cy, r, ccw) for an arc about the profile normal, or a Spline
    # (S8b) whose first and last poles are points j and j+1 exactly.
    segments: list = None


@dataclass(frozen=True)
class Spline:
    """S8b: a nonrational planar B-spline path segment of `degree` with
    `poles` ((x, y) binary64 pairs, the first and last the segment's points
    exactly), distinct increasing `knots` and their `mults` (clamped: the
    ends degree + 1, interior ones 1 to degree)."""
    degree: int
    poles: tuple
    knots: tuple
    mults: tuple

    def check(self, p, q):
        n, d = len(self.poles), self.degree
        assert 1 <= d <= 7 and n >= d+1, 'spline degree and poles'
        assert tuple(self.poles[0]) == tuple(p) and tuple(self.poles[-1]) == tuple(q), \
            'spline poles must start and end at the segment points exactly'
        assert len(self.knots) == len(self.mults) >= 2
        assert all(a < b for a, b in zip(self.knots, self.knots[1:])), 'knots increase'
        assert self.mults[0] == self.mults[-1] == d+1 and all(1 <= m <= d for m in self.mults[1:-1])
        assert sum(self.mults) == n+d+1, 'knot count'

    def reversed(self):
        """The same curve backwards on the same domain: poles reversed,
        knots mirrored `k -> a + b - k` (exact in binary64, asserted),
        multiplicities reversed."""
        a, b = F(self.knots[0]), F(self.knots[-1])
        knots = []
        for k in reversed(self.knots):
            m = a+b-F(k)
            assert F(float(m)) == m, 'mirrored knot not exact in binary64'
            knots.append(float(m))
        return Spline(self.degree, tuple(reversed(self.poles)), tuple(knots), tuple(reversed(self.mults)))

    def flat_knots(self):
        return [F(k) for k, m in zip(self.knots, self.mults) for _ in range(m)]

    def blossom(self, span, args):
        """The polar form of the span's polynomial at `args` (degree many
        Fractions): de Boor's recurrence with one argument per level."""
        U, d = self.flat_knots(), self.degree
        P = [(F(x), F(y)) for x, y in self.poles]
        pts = {j: P[j] for j in range(span-d, span+1)}
        for r in range(1, d+1):
            t = args[r-1]
            for j in range(span, span-d+r-1, -1):
                a = (t-U[j])/(U[j+d+1-r]-U[j])
                pts[j] = tuple((1-a)*pts[j-1][c]+a*pts[j][c] for c in range(2))
        return pts[span]

    def pieces(self):
        """The Bezier control points of each nonempty knot span, exactly, by
        blossoming: point `i` of span `[s, e]` is the blossom at `s`
        repeated `degree - i` times and `e` repeated `i` times."""
        U, d = self.flat_knots(), self.degree
        out = []
        for span in range(d, len(U)-d-1):
            s, e = U[span], U[span+1]
            if s < e:
                out.append(tuple(self.blossom(span, [s]*(d-i)+[e]*i) for i in range(d+1)))
        return out

    def twice_area(self):
        """The exact integral of `x dy - y dx` along the spline."""
        total = F(0)
        for ctrl in self.pieces():
            total += bezier_cross_integral(ctrl)
        return total


def bernstein_product(n, i, m, j):
    """The integral over [0, 1] of B(n, i) B(m, j)."""
    from math import comb
    return F(comb(n, i)*comb(m, j), comb(n+m, i+j)*(n+m+1))


def bezier_cross_integral(ctrl):
    """The exact integral of `x dy - y dx` along a Bezier curve with
    Fraction control points."""
    n = len(ctrl)-1
    d = [tuple(n*(ctrl[j+1][c]-ctrl[j][c]) for c in range(2)) for j in range(n)]
    total = F(0)
    for i in range(n+1):
        for j in range(n):
            w = bernstein_product(n, i, n-1, j)
            total += w*(ctrl[i][0]*d[j][1]-ctrl[i][1]*d[j][0])
    return total


def arc_sweep(p, q, arc):
    """The signed sweep of an arc from p to q (mpmath): in (0, 2π) turning
    counter-clockwise, in (-2π, 0) clockwise."""
    import mpmath as mp
    cx, cy, _, ccw = arc
    a = mp.atan2(mp.mpf(p[1])-mp.mpf(cy), mp.mpf(p[0])-mp.mpf(cx))
    b = mp.atan2(mp.mpf(q[1])-mp.mpf(cy), mp.mpf(q[0])-mp.mpf(cx))
    turn = (b-a) % (2*mp.pi) if ccw else -((a-b) % (2*mp.pi))
    assert turn != 0
    return turn


def path_area(pts, segments):
    """Twice the signed area of a path: the polygon of its points plus each
    arc's circular segment, (r^2)(φ - sin φ) for the signed sweep φ, and
    each spline's exact `x dy - y dx` integral less its chord's."""
    import mpmath as mp
    n = len(pts)
    twice = mp.mpf(0)
    for i in range(n):
        p, q = pts[i], pts[(i+1) % n]
        twice += mp.mpf(p[0])*q[1]-mp.mpf(q[0])*p[1]
        if isinstance(segments[i], Spline):
            segments[i].check(p, q)
            extra = segments[i].twice_area()-(F(p[0])*F(q[1])-F(q[0])*F(p[1]))
            twice += mp.mpf(extra.numerator)/extra.denominator
        elif segments[i] is not None:
            phi = arc_sweep(p, q, segments[i])
            twice += mp.mpf(segments[i][2])**2*(phi-mp.sin(phi))
    return twice


def reversed_segment(seg):
    """A path segment traversed backwards: a line stays one, an arc turns
    the other way, a spline is reversed (Spline.reversed)."""
    if seg is None:
        return None
    if isinstance(seg, Spline):
        return seg.reversed()
    return (*seg[:3], not seg[3])


def stored(boundary, tolerance):
    """Counter-clockwise stored points and labels mapped into stored order.

    Boundary::polygon drops a closing point within tolerance of the first,
    then reverses points[1..] when the polygon is clockwise. Segment j of the
    stored polygon runs from stored point j to j+1. A path (S5) keeps every
    point, is oriented by its area with arcs' and splines' bulges, and
    reversed segments flip their arcs' directions and reverse their splines
    (S8b)."""
    if boundary.circle is not None:
        return None, boundary.labels
    if boundary.segments is not None:
        pts, segs = list(boundary.points), list(boundary.segments)
        n = len(pts)
        twice = path_area(pts, segs)
        import mpmath as mp
        assert abs(twice) > mp.mpf(2)**-40, 'path orientation too close to call'
        labels = boundary.labels
        if twice < 0:
            pts = [pts[0]]+pts[:0:-1]
            segs = [reversed_segment(segs[n-1-j]) for j in range(n)]
            if labels is not None:
                b, seg, vert = labels
                labels = (b, [seg[n-1-j] for j in range(n)], [vert[(n-j) % n] for j in range(n)])
        return (pts, segs), labels
    pts = list(boundary.points)
    if len(pts) > 1:
        dx, dy = F(pts[0][0])-F(pts[-1][0]), F(pts[0][1])-F(pts[-1][1])
        if dx*dx+dy*dy <= F(tolerance)**2:
            pts.pop()
    n = len(pts)
    area = sum(F(pts[i][0])*F(pts[(i+1) % n][1])-F(pts[(i+1) % n][0])*F(pts[i][1])
               for i in range(n))
    assert area != 0
    labels = boundary.labels
    if area < 0:
        pts = [pts[0]]+pts[:0:-1]
        if labels is not None:
            b, seg, vert = labels
            labels = (b, [seg[n-1-j] for j in range(n)], [vert[(n-j) % n] for j in range(n)])
    return pts, labels


def stored_points(boundary, tolerance):
    """The stored points of a polygon or path (None for a circle)."""
    pts, _ = stored(boundary, tolerance)
    return pts[0] if boundary.segments is not None else pts


@dataclass
class Case:
    name: str
    tolerance: float
    operation: int
    frame: tuple                 # origin, normal, x hint (9 floats)
    start: float
    end: float
    boundaries: list             # outer first, then holes
    transforms: list = field(default_factory=list)   # ('T', v3) | ('R', origin, axis, angle)
    box: tuple = None            # (origin3, size3) for Solid::box_at instead of a profile
    cone: tuple = None           # (r1, r2, height) for Solid::cone_with on the frame
    sphere: tuple = None         # (radius, low, high) for Solid::sphere_with on the frame
    torus: tuple = None          # (major, minor, low, high, angle) for Solid::torus_with
    make: str = None             # 'face' (Body::face_from_profile) or 'wire' (the first boundary's)


def number(x):
    return repr(float(x))


def segment_words(seg, lift, dim):
    """A path segment's words: L, A cx cy r ccw, or B p n poles k knots
    mults (S8b) with each pole mapped by `lift` to `dim` coordinates."""
    if seg is None:
        return ['L']
    if isinstance(seg, Spline):
        words = ['B', str(seg.degree), str(len(seg.poles))]
        for pole in seg.poles:
            v = lift(pole)
            assert len(v) == dim
            words += [number(x) for x in v]
        words += [str(len(seg.knots)), *map(number, seg.knots), *map(str, seg.mults)]
        return words
    return ['A', *map(number, lift(seg[:2])), number(seg[2]), '1' if seg[3] else '0']


def encode_case(c):
    out = [f'case {c.name} {number(c.tolerance)}', f'op {c.operation}']
    if c.box is not None:
        out.append('box '+' '.join(number(x) for x in (*c.box[0], *c.box[1])))
    elif c.cone is not None:
        out.append('frame '+' '.join(number(x) for x in c.frame))
        out.append('cone '+' '.join(number(x) for x in c.cone))
    elif c.sphere is not None:
        out.append('frame '+' '.join(number(x) for x in c.frame))
        out.append('sphere '+' '.join(number(x) for x in c.sphere))
    elif c.torus is not None:
        out.append('frame '+' '.join(number(x) for x in c.frame))
        out.append('torus '+' '.join(number(x) for x in c.torus))
    else:
        out.append('frame '+' '.join(number(x) for x in c.frame))
        out.append(f'make {c.make}' if c.make else f'offsets {number(c.start)} {number(c.end)}')
        for b in c.boundaries:
            if b.circle is not None:
                row = 'boundary C '+' '.join(number(x) for x in b.circle)
            elif b.segments is not None:
                # S5: each point, then its segment to the next: L, or
                # A cx cy r and 1 (counter-clockwise) or 0, or (S8b) a
                # nonrational spline B p n x0 y0 ... x(n-1) y(n-1) k u0 ...
                # u(k-1) m0 ... m(k-1): degree p, n poles (the first the
                # point, the last the next point, exactly), k distinct knots
                # and their multiplicities (clamped: the ends p + 1, interior
                # ones 1 to p).
                words = [f'boundary S {len(b.points)}']
                for p, seg in zip(b.points, b.segments):
                    words += [number(p[0]), number(p[1])]
                    words += segment_words(seg, lambda v: tuple(v), 2)
                row = ' '.join(words)
            else:
                row = f'boundary P {len(b.points)} '+' '.join(number(x) for p in b.points for x in p)
            if b.labels is not None:
                bl, seg, vert = b.labels
                row += f' labels {bl} '+' '.join(map(str, seg))+' | '+' '.join(map(str, vert))
            out.append(row)
    for t in c.transforms:
        if t[0] == 'T':
            out.append('transform T '+' '.join(number(x) for x in t[1]))
        else:
            out.append('transform R '+' '.join(number(x) for x in (*t[1], *t[2], t[3])))
    out.append('end')
    return '\n'.join(out)


BOOLEAN_OPERATIONS = ('fuse', 'cut', 'common')


def encode_boolean_case(obj, operation, tool, boolean_operation):
    """S9a: a Boolean of two prisms in the protocol above. The object's
    block (`encode_case(obj)`) without its `end`, then a row `boolean OP ID`
    (`OP` one of fuse, cut, common: `Solid::fuse(ID, tool)` and so on, the
    object the receiver, `ID` the Boolean's operation id), then the tool's
    rows (`encode_case(tool)` without its `case` row and its `end`: `op`,
    `frame`, `offsets`, `boundary`...), then `end`:

        case NAME TOL
        op 91
        frame ...
        offsets ...
        boundary ...
        boolean fuse 93
        op 92
        frame ...
        offsets ...
        boundary ...
        end

    A reader splits the block at the `boolean` row and parses each side as
    an identity case, the tool's with the object's `case` row (its name and
    tolerance). Blocks without a `boolean` row are unchanged: every existing
    reader of a case still reads them.

    S9d.1: either input may be a sphere, a cap or a zone (its `frame` and
    `sphere R LOW HIGH` rows, `Solid::sphere_with` on the frame, the
    latitudes in radians), read as any identity case's. S9d.3a: likewise a
    cone or frustum (its `frame` and `cone BOTTOM TOP HEIGHT` rows,
    `Solid::cone_with` on the frame)."""
    assert operation in BOOLEAN_OPERATIONS, operation
    assert obj.make is None and tool.make is None and tool.box is None \
        and tool.torus is None, 'S9a: two prisms (S9d.1: or a sphere; S9d.3a: or a cone)'
    assert not obj.transforms and not tool.transforms, 'S9a: prisms in place'
    first = encode_case(obj).rsplit('\nend', 1)[0]
    second = encode_case(tool).split('\n')[1:-1]
    return '\n'.join([first, f'boolean {operation} {boolean_operation}', *second, 'end'])


def native_boolean_case(obj, operation, tool):
    """The explicit OCCT rows of a Boolean (occt_boolean_oracle.cpp): the
    object's `native_case` rows without `end`, a `boolean OP` row, the
    tool's rows without their `case` row, and `end`."""
    assert operation in BOOLEAN_OPERATIONS, operation
    first = native_case(obj).rsplit('\nend', 1)[0]
    second = native_case(tool).split('\n')[1:-1]
    return '\n'.join([first, f'boolean {operation}', *second, 'end'])


def box_boundaries(size):
    """Solid::box_at builds a cuboid from a rectangle at the origin."""
    w, d = abs(size[0]), abs(size[1])
    return [Boundary(points=[(0.0, 0.0), (w, 0.0), (w, d), (0.0, d)])]


# ------------------------------------------------------------------ extrusion

@dataclass
class Entity:
    kind: str
    derivation: Derivation
    locator: tuple   # (boundary, element, index, side) or ('cap', side)

    @property
    def id(self):
        return self.derivation.id()


def extrude_entities(c):
    """Every vertex, edge and face of the extrusion, independently of the
    Rust builder. 'Bottom' and 'top' mean the start and end sides."""
    op, tol = c.operation, c.tolerance
    boundaries = box_boundaries(c.box[1]) if c.box is not None else c.boundaries
    if c.box is not None:
        op = 0
    ents = []

    def add(kind, role, parents, locator, ordinal=0):
        ents.append(Entity(kind, Derivation(op, 'extrude', kind, role, ordinal, tuple(parents)), locator))

    labels = [l for b in boundaries if b.labels for l in (b.labels[0], *b.labels[1], *b.labels[2])]
    assert len(labels) == len(set(labels)), f'{c.name}: duplicate labels'
    cap_parents = []
    for b, boundary in enumerate(boundaries):
        _, labels = stored(boundary, tol)
        cap_parents.append(('label', labels[0]) if labels else ('profile', b, 'boundary', 0))
    add('face', 'start_cap', cap_parents, ('cap', 'start'))
    add('face', 'end_cap', cap_parents, ('cap', 'end'))
    # The solid region (TOPOLOGY_MODEL.md, T1).
    add('region', 'region', cap_parents, ('region',))
    for b, boundary in enumerate(boundaries):
        pts, labels = stored(boundary, tol)
        if boundary.segments is not None:
            # A path's segments and points are a polygon's (S5).
            pts = pts[0]
        n = 1 if pts is None else len(pts)
        seg = (lambda j: ('label', labels[1][j])) if labels else (lambda j: ('profile', b, 'segment', j))
        vert = (lambda j: ('label', labels[2][j])) if labels else (lambda j: ('profile', b, 'vertex', j))
        if pts is None:
            # Seamless (T1): two ring edges and the wall; no seam, no vertices.
            add('edge', 'bottom_edge', [seg(0)], (b, 'segment', 0, 'start'))
            add('edge', 'top_edge', [seg(0)], (b, 'segment', 0, 'end'))
            add('face', 'wall', [seg(0)], (b, 'segment', 0, 'both'))
            continue
        for j in range(n):
            add('vertex', 'bottom_vertex', [vert(j)], (b, 'vertex', j, 'start'))
            add('vertex', 'top_vertex', [vert(j)], (b, 'vertex', j, 'end'))
            add('edge', 'bottom_edge', [seg(j)], (b, 'segment', j, 'start'))
            add('edge', 'top_edge', [seg(j)], (b, 'segment', j, 'end'))
            add('edge', 'vertical', [vert(j)], (b, 'vertex', j, 'both'))
            add('face', 'wall', [seg(j)], (b, 'segment', j, 'both'))
    ids = [e.id for e in ents]
    assert len(set(ids)) == len(ids), f'{c.name}: id collision'
    return ents


def cone_entities(c):
    """Every entity of Solid::cone_with (S3), independently of the Rust
    builder. The meridian is boundary 0, the polygon (0, 0), (r1, 0),
    (r2, h), (0, h) in (radius, height) with segment j from point j to j+1:
    the rim points 1 and 2 revolve into the bottom and top rings, or are
    apices when their radius is 0; segments 0 and 2 into the discs, segment 1
    into the wall, the boundary into the solid region. The axis (points 0
    and 3, segment 3) generates nothing."""
    op = c.operation
    r1, r2, _ = c.cone
    ents = []

    def add(kind, role, parents, locator):
        ents.append(Entity(kind, Derivation(op, 'revolve', kind, role, 0, tuple(parents)), locator))

    def meridian(element, index):
        return ('profile', 0, element, index)
    add('region', 'region', [meridian('boundary', 0)], ('region',))
    for r, side, cap, segment, edge, rim in ((r1, 'start', 'start_cap', 0, 'bottom_edge', 1),
                                             (r2, 'end', 'end_cap', 2, 'top_edge', 2)):
        if r == 0:
            add('vertex', 'apex', [meridian('vertex', rim)], ('apex', side))
        else:
            add('edge', edge, [meridian('vertex', rim)], ('ring', side))
            add('face', cap, [meridian('segment', segment)], ('cap', side))
    add('face', 'wall', [meridian('segment', 1)], ('wall',))
    ids = [e.id for e in ents]
    assert len(set(ids)) == len(ids), f'{c.name}: id collision'
    return ents


HALF_PI = 1.5707963267948966


def sphere_entities(c):
    """Every entity of Solid::sphere_with (S3), independently of the Rust
    builder. The meridian is boundary 0 as for the cone: points (0, z1),
    (r1, z1), (r2, z2), (0, z2), the arc from point 1 to point 2 being
    segment 1. An end at latitude +-pi/2 (the binary64 value) is a pole: its
    rim point is a vertex (role pole) when the other end is not a pole, and
    nothing when both are (the whole sphere has no loops); any other end
    gives a ring and a disc."""
    op = c.operation
    _, low, high = c.sphere
    ends = [(low, low == -HALF_PI, 'start', 'start_cap', 0, 'bottom_edge', 1),
            (high, high == HALF_PI, 'end', 'end_cap', 2, 'top_edge', 2)]
    poles = sum(1 for e in ends if e[1])
    ents = []

    def add(kind, role, parents, locator):
        ents.append(Entity(kind, Derivation(op, 'revolve', kind, role, 0, tuple(parents)), locator))

    def meridian(element, index):
        return ('profile', 0, element, index)
    add('region', 'region', [meridian('boundary', 0)], ('region',))
    for _, pole, side, cap, segment, edge, rim in ends:
        if pole:
            if poles == 1:
                add('vertex', 'pole', [meridian('vertex', rim)], ('pole', side))
        else:
            add('edge', edge, [meridian('vertex', rim)], ('ring', side))
            add('face', cap, [meridian('segment', segment)], ('cap', side))
    add('face', 'wall', [meridian('segment', 1)], ('wall',))
    ids = [e.id for e in ents]
    assert len(set(ids)) == len(ids), f'{c.name}: id collision'
    return ents


TWO_PI = 6.283185307179586


def torus_entities(c):
    """Every entity of Solid::torus_with (S3), independently of the Rust
    builder, on the cone's meridian convention (the tube's arc from latitude
    low to high is segment 1): a whole torus is its wall and region; a
    v-segment has rings from rim points 1 and 2 and discs from segments 0 and
    2; a wedge has discs from the boundary (the meridian's start and end
    copies) and the tube's circles from the arc (segment 1's copies)."""
    op = c.operation
    _, _, low, high, angle = c.torus
    closed, turn = high-low == TWO_PI, angle == TWO_PI
    ents = []

    def add(kind, role, parents, locator):
        ents.append(Entity(kind, Derivation(op, 'revolve', kind, role, 0, tuple(parents)), locator))

    def meridian(element, index):
        return ('profile', 0, element, index)
    add('region', 'region', [meridian('boundary', 0)], ('region',))
    if not closed:
        for side, cap, segment, edge, rim in (('start', 'start_cap', 0, 'bottom_edge', 1),
                                              ('end', 'end_cap', 2, 'top_edge', 2)):
            add('edge', edge, [meridian('vertex', rim)], ('ring', side))
            add('face', cap, [meridian('segment', segment)], ('cap', side))
    elif not turn:
        for side, cap, edge in (('start', 'start_cap', 'bottom_edge'), ('end', 'end_cap', 'top_edge')):
            add('edge', edge, [meridian('segment', 1)], ('ring', side))
            add('face', cap, [meridian('boundary', 0)], ('cap', side))
    add('face', 'wall', [meridian('segment', 1)], ('wall',))
    ids = [e.id for e in ents]
    assert len(set(ids)) == len(ids), f'{c.name}: id collision'
    return ents


def sheet_entities(c):
    """S6: a face body's face (from every boundary), and each boundary's
    edges (from their segments) and vertices (from theirs); a wire body the
    edges and vertices of its one boundary. A circle is one ring edge."""
    op, tol, kind = c.operation, c.tolerance, 'make_'+c.make
    boundaries = c.boundaries if c.make == 'face' else c.boundaries[:1]
    ents = []

    def add(entity, role, parents, locator):
        ents.append(Entity(entity, Derivation(op, kind, entity, role, 0, tuple(parents)), locator))

    face_parents = []
    for b, boundary in enumerate(boundaries):
        pts, labels = stored(boundary, tol)
        if boundary.segments is not None:
            pts = pts[0]
        face_parents.append(('label', labels[0]) if labels else ('profile', b, 'boundary', 0))
        seg = (lambda j: ('label', labels[1][j])) if labels else (lambda j: ('profile', b, 'segment', j))
        vert = (lambda j: ('label', labels[2][j])) if labels else (lambda j: ('profile', b, 'vertex', j))
        if pts is None:
            add('edge', 'edge', [seg(0)], (b, 'segment', 0))
            continue
        for j in range(len(pts)):
            add('vertex', 'vertex', [vert(j)], (b, 'vertex', j))
            add('edge', 'edge', [seg(j)], (b, 'segment', j))
    if c.make == 'face':
        add('face', 'face', face_parents, ('face',))
    ids = [e.id for e in ents]
    assert len(set(ids)) == len(ids), f'{c.name}: id collision'
    return ents


def entities(c):
    """Every entity of a case's construction."""
    if c.make is not None:
        return sheet_entities(c)
    if c.torus is not None:
        return torus_entities(c)
    if c.cone is not None:
        return cone_entities(c)
    if c.sphere is not None:
        return sphere_entities(c)
    return extrude_entities(c)


def body_id(operation, kind='extrude', parents=(), ordinal=0):
    return Derivation(operation, kind, 'body', 'body', ordinal, tuple(parents)).id()


def entity_text(e):
    """One fixture row: id kind role ordinal parents locator."""
    d = e.derivation
    loc = ' '.join(map(str, e.locator))
    return f'{hexid(e.id)} {e.kind} {d.role} {d.ordinal} {",".join(parent_text(p) for p in d.parents)} {loc}'


# ------------------------------------------------------------------ relations

def relation_sort_key(r):
    """Canonical order: sources, then relation kind, then targets (and role)."""
    kind = r[0]
    if kind in ('unchanged', 'modified', 'deleted'):
        sources, targets = [('entity', r[1])], [r[-1]] if kind != 'deleted' else []
    elif kind == 'generated':
        sources, targets = list(r[1]), [r[2]]
    elif kind == 'split':
        sources, targets = [('entity', r[1])], list(r[2])
    else:
        sources, targets = [('entity', f) for f in r[1]], [r[2]]
    key = encode_parents(sources)+bytes([RELATION[kind]])+struct.pack('<I', len(targets))+b''.join(targets)
    if kind == 'generated':
        key += bytes([ROLE[r[3]]])
    return key


def relation_text(r):
    kind = r[0]
    if kind == 'generated':
        return f'generated {",".join(parent_text(p) for p in r[1])} {r[2].hex()} {r[3]}'
    if kind == 'modified':
        return f'modified {r[1].hex()} {r[2].hex()}'
    if kind in ('unchanged', 'deleted'):
        return f'{kind} {r[1].hex()}'
    if kind == 'split':
        return f'split {r[1].hex()} {",".join(t.hex() for t in r[2])}'
    return f'merged {",".join(f.hex() for f in r[1])} {r[2].hex()}'


def extrude_history(c):
    """The construction's relations: a prism's or a cone's."""
    rels = [('generated', e.derivation.parents, e.id, e.derivation.role) for e in entities(c)]
    return sorted(rels, key=relation_sort_key)


def transform_history(c):
    rels = [('modified', e.id, e.id) for e in entities(c)]
    return sorted(rels, key=relation_sort_key)


# ------------------------------------------------------------------ native rows

def _unit(v):
    from brep_reference import hypot_rn
    n = hypot_rn(*v)
    return tuple(x/n for x in v)


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def frame_axes(frame):
    """Frame3::new: unit normal, y = normal x hint, x = y x normal."""
    origin, normal, hint = frame[0:3], _unit(frame[3:6]), _unit(frame[6:9])
    y = _unit(_cross(normal, hint))
    x = _unit(_cross(y, normal))
    return origin, x, y, normal


def rotation(origin, axis, angle):
    """RigidTransform::rotation as a 3x4 row-major matrix."""
    from brep_reference import cos_rn, ltr_sum, sin_rn
    a = _unit(axis)
    c, s = cos_rn(angle), sin_rn(angle)
    cols = []
    for e in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        cr = _cross(a, e)
        d = ltr_sum(p*q for p, q in zip(a, e))
        cols.append(tuple(e[i]*c+cr[i]*s+a[i]*(d*(1.0-c)) for i in range(3)))
    rotated = tuple(ltr_sum(cols[j][i]*origin[j] for j in range(3)) for i in range(3))
    t = tuple(origin[i]-rotated[i] for i in range(3))
    return [[cols[0][i], cols[1][i], cols[2][i], t[i]] for i in range(3)]


def transform_matrix(t):
    if t[0] == 'T':
        return [[1.0, 0.0, 0.0, t[1][0]], [0.0, 1.0, 0.0, t[1][1]], [0.0, 0.0, 1.0, t[1][2]]]
    return rotation(t[1], t[2], t[3])


def native_case(c):
    """Explicit OCCT construction rows for occt_history_oracle.cpp (and
    occt_split_oracle.cpp). A face or wire body (S6, `make`) has its frame's
    plane (the start offset is 0), its boundaries' `wire` rows (a wire body
    its first boundary's only) and a `make face` or `make wire` row in place
    of the `prism` vector (S8e). A sphere, cap or zone (S9d.1) is one row
    `sphere ox oy oz nx ny nz xx xy xz R LOW HIGH`, its frame's origin,
    normal and x axis and `Solid::sphere_with`'s radius and latitudes
    (radians), as `BRepPrimAPI_MakeSphere(gp_Ax2, R, LOW, HIGH)` takes them.
    A cone or frustum (S9d.3a) is one row `cone ox oy oz nx ny nz xx xy xz
    R1 R2 H`, `Solid::cone_with`'s bottom and top radii and height, as
    `BRepPrimAPI_MakeCone(gp_Ax2, R1, R2, H)` takes them."""
    if c.sphere is not None:
        assert not c.transforms and c.make is None, f'{c.name}: a sphere in place'
        o, x, _, n = frame_axes(c.frame)
        return '\n'.join([f'case {c.name}', 'sphere '+' '.join(number(v) for v in (*o, *n, *x, *c.sphere)),
                          'end'])
    if c.cone is not None:
        assert not c.transforms and c.make is None, f'{c.name}: a cone in place'
        o, x, _, n = frame_axes(c.frame)
        return '\n'.join([f'case {c.name}', 'cone '+' '.join(number(v) for v in (*o, *n, *x, *c.cone)), 'end'])
    if c.box is not None:
        (ox, oy, oz), size = c.box
        frame = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
        start, end, boundaries = 0.0, abs(size[2]), box_boundaries(size)
        corner = (ox+min(size[0], 0.0), oy+min(size[1], 0.0), oz+min(size[2], 0.0))
        transforms = [('T', corner)]+list(c.transforms)
    else:
        frame, start, end, boundaries, transforms = c.frame, c.start, c.end, c.boundaries, c.transforms
    if c.make is not None:
        assert start == 0.0 and end == 0.0 and c.make in ('face', 'wire'), f'{c.name}: a face or wire body'
        if c.make == 'wire':
            boundaries = boundaries[:1]
    o, x, y, n = frame_axes(frame)
    at = lambda p: tuple(o[i]+x[i]*p[0]+y[i]*p[1]+n[i]*start for i in range(3))
    rows = [f'case {c.name}', 'plane '+' '.join(number(v) for v in (*at((0.0, 0.0)), *n, *x))]
    for b in boundaries:
        pts, _ = stored(b, c.tolerance)
        if pts is None:
            cx, cy, r = b.circle
            rows.append('wire C '+' '.join(number(v) for v in (*at((cx, cy)), r)))
        elif b.segments is not None:
            # Stored counter-clockwise: each point, then its segment to the
            # next (L, or A with the 3D centre, radius and 1 if it turns
            # counter-clockwise about the plane's normal, or B p n with the
            # 3D poles, k, the knots and multiplicities: S8b).
            points, segments = pts
            words = [f'wire S {len(points)}']
            for p, seg in zip(points, segments):
                words += [number(v) for v in at(p)]
                words += segment_words(seg, at, 3)
            rows.append(' '.join(words))
        else:
            rows.append(f'wire P {len(pts)} '+' '.join(number(v) for p in pts for v in at(p)))
    if c.make is not None:
        rows.append(f'make {c.make}')
    else:
        rows.append('prism '+' '.join(number(n[i]*(end-start)) for i in range(3)))
    for t in transforms:
        rows.append('transform '+' '.join(number(v) for row in transform_matrix(t) for v in row))
    rows.append('end')
    return '\n'.join(rows)
