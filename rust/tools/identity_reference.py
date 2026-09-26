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
"""
from dataclasses import dataclass, field
from fractions import Fraction as F
import struct

FNV_OFFSET = 0x6c62272e07bb014262b821756295c58d
FNV_PRIME = 0x0000000001000000000000000000013B
MASK = (1 << 128)-1

KIND = {'extrude': 1, 'transform': 2, 'external': 3, 'composite': 4, 'height_split': 5,
        'stacked_fuse': 6}
ENTITY = {'vertex': 1, 'edge': 2, 'face': 3, 'body': 4, 'region': 5}
DIMENSION = {'vertex': 0, 'edge': 1, 'face': 2, 'body': 3, 'region': 3}
ROLE = {'start_cap': 1, 'end_cap': 2, 'wall': 3, 'bottom_edge': 4, 'top_edge': 5,
        'vertical': 6, 'seam': 7, 'bottom_vertex': 8, 'top_vertex': 9,
        'seam_vertex': 10, 'body': 11, 'external': 12, 'region': 13, 'cut_face': 14,
        'cut_edge': 15, 'cut_vertex': 16}
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
    """Caller input: a polygon (points) or a circle, with optional labels in
    the caller's input order."""
    points: list = None          # [(x, y)] for a polygon
    circle: tuple = None         # (cx, cy, r)
    labels: tuple = None         # (boundary, [segment...], [vertex...])


def stored(boundary, tolerance):
    """Counter-clockwise stored points and labels mapped into stored order.

    Boundary::polygon drops a closing point within tolerance of the first,
    then reverses points[1..] when the polygon is clockwise. Segment j of the
    stored polygon runs from stored point j to j+1."""
    if boundary.circle is not None:
        return None, boundary.labels
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


def number(x):
    return repr(float(x))


def encode_case(c):
    out = [f'case {c.name} {number(c.tolerance)}', f'op {c.operation}']
    if c.box is not None:
        out.append('box '+' '.join(number(x) for x in (*c.box[0], *c.box[1])))
    else:
        out.append('frame '+' '.join(number(x) for x in c.frame))
        out.append(f'offsets {number(c.start)} {number(c.end)}')
        for b in c.boundaries:
            if b.circle is not None:
                row = 'boundary C '+' '.join(number(x) for x in b.circle)
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
    rels = [('generated', e.derivation.parents, e.id, e.derivation.role) for e in extrude_entities(c)]
    return sorted(rels, key=relation_sort_key)


def transform_history(c):
    rels = [('modified', e.id, e.id) for e in extrude_entities(c)]
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
    from brep_reference import cos_rn, sin_rn
    a = _unit(axis)
    c, s = cos_rn(angle), sin_rn(angle)
    cols = []
    for e in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        cr = _cross(a, e)
        d = sum(p*q for p, q in zip(a, e))
        cols.append(tuple(e[i]*c+cr[i]*s+a[i]*(d*(1.0-c)) for i in range(3)))
    rotated = tuple(sum(cols[j][i]*origin[j] for j in range(3)) for i in range(3))
    t = tuple(origin[i]-rotated[i] for i in range(3))
    return [[cols[0][i], cols[1][i], cols[2][i], t[i]] for i in range(3)]


def transform_matrix(t):
    if t[0] == 'T':
        return [[1.0, 0.0, 0.0, t[1][0]], [0.0, 1.0, 0.0, t[1][1]], [0.0, 0.0, 1.0, t[1][2]]]
    return rotation(t[1], t[2], t[3])


def native_case(c):
    """Explicit OCCT construction rows for occt_history_oracle.cpp."""
    if c.box is not None:
        (ox, oy, oz), size = c.box
        frame = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
        start, end, boundaries = 0.0, abs(size[2]), box_boundaries(size)
        corner = (ox+min(size[0], 0.0), oy+min(size[1], 0.0), oz+min(size[2], 0.0))
        transforms = [('T', corner)]+list(c.transforms)
    else:
        frame, start, end, boundaries, transforms = c.frame, c.start, c.end, c.boundaries, c.transforms
    o, x, y, n = frame_axes(frame)
    at = lambda p: tuple(o[i]+x[i]*p[0]+y[i]*p[1]+n[i]*start for i in range(3))
    rows = [f'case {c.name}', 'plane '+' '.join(number(v) for v in (*at((0.0, 0.0)), *n, *x))]
    for b in boundaries:
        pts, _ = stored(b, c.tolerance)
        if pts is None:
            cx, cy, r = b.circle
            rows.append('wire C '+' '.join(number(v) for v in (*at((cx, cy)), r)))
        else:
            rows.append(f'wire P {len(pts)} '+' '.join(number(v) for p in pts for v in at(p)))
    rows.append('prism '+' '.join(number(n[i]*(end-start)) for i in range(3)))
    for t in transforms:
        rows.append('transform '+' '.join(number(v) for row in transform_matrix(t) for v in row))
    rows.append('end')
    return '\n'.join(rows)
