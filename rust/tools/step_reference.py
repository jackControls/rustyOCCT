#!/usr/bin/env python3
"""Independent reference for the STEP import track (STEP-a of REVIEW_NOTES.md).

Three parts, none of which shares code or data with the kernel:

* `parse`: a Part 21 reader written from ISO 10303-21 (header, DATA
  sections, simple and complex instances, every parameter kind), returning
  the header records and the instances by entity number.
* `bodies`: the solids and surface models of a parsed file with the counts
  OCCT's reader gives them (distinct vertices, edges, wires, faces, shells,
  solids; one degenerated edge per passage of a loop through a pole of a
  sphere or the apex of a cone, which `STEPControl_Reader` adds), and the
  length unit and uncertainty of their representation context.
* closed forms of the volume, area and centre of every fixture shape, from
  its construction parameters (`generate_step_fixtures.py`), in mpmath at
  60 digits.
"""
from fractions import Fraction
import re

import mpmath

mpmath.mp.dps = 60


class Ref(int):
    """`#n`."""


class Enum(str):
    """`.NAME.`"""


class Typed:
    """`NAME(value)`: a typed parameter such as `LENGTH_MEASURE(1.)`."""

    def __init__(self, name, value):
        self.name, self.value = name, value

    def __repr__(self):
        return f'{self.name}({self.value!r})'


NULL = None
DERIVED = '*'

TOKEN = re.compile(r"""
    (?P<space>\s+)
  | (?P<comment>/\*.*?\*/)
  | (?P<ref>\#[0-9]+)
  | (?P<real>[+-]?[0-9]+\.[0-9]*(?:[Ee][+-]?[0-9]+)?)
  | (?P<int>[+-]?[0-9]+)
  | (?P<string>'(?:[^']|'')*')
  | (?P<enum>\.[A-Z_][A-Z0-9_]*\.)
  | (?P<binary>"[0-9A-F]*")
  | (?P<keyword>!?[A-Z_][A-Z0-9_]*(?:-[0-9A-Z_]+)*)
  | (?P<punct>[()=;,$*])
""", re.X | re.S)


def tokens(text):
    at = 0
    while at < len(text):
        m = TOKEN.match(text, at)
        if not m:
            raise ValueError(f'unexpected character at offset {at}: {text[at:at+20]!r}')
        at = m.end()
        kind = m.lastgroup
        if kind in ('space', 'comment'):
            continue
        yield kind, m.group()


class Stream:
    def __init__(self, text):
        self.items = list(tokens(text))
        self.at = 0

    def peek(self):
        return self.items[self.at] if self.at < len(self.items) else (None, None)

    def take(self, value=None):
        kind, text = self.peek()
        if kind is None or (value is not None and text != value):
            raise ValueError(f'expected {value!r}, found {text!r}')
        self.at += 1
        return kind, text


def parameter(s):
    kind, text = s.peek()
    if text == '(':
        s.take('(')
        out = []
        if s.peek()[1] != ')':
            while True:
                out.append(parameter(s))
                if s.peek()[1] == ',':
                    s.take(',')
                    continue
                break
        s.take(')')
        return out
    s.take()
    if kind == 'ref':
        return Ref(int(text[1:]))
    if kind == 'real':
        return float(text)
    if kind == 'int':
        return int(text)
    if kind == 'string':
        return text[1:-1].replace("''", "'")
    if kind == 'enum':
        return Enum(text[1:-1])
    if kind == 'binary':
        return text
    if text == '$':
        return NULL
    if text == '*':
        return DERIVED
    if kind == 'keyword':
        s.take('(')
        value = parameter(s)
        s.take(')')
        return Typed(text, value)
    raise ValueError('unexpected parameter '+text)


def record(s):
    _, name = s.take()
    return name, parameter(s)


def parse(text):
    """(header records, {entity number: [(NAME, [parameters]), ...]})."""
    s = Stream(text)
    s.take('ISO-10303-21')
    s.take(';')
    s.take('HEADER')
    s.take(';')
    header = []
    while s.peek()[1] != 'ENDSEC':
        header.append(record(s))
        s.take(';')
    s.take('ENDSEC')
    s.take(';')
    data = {}
    while s.peek()[1] == 'DATA':
        s.take('DATA')
        if s.peek()[1] == '(':
            parameter(s)
        s.take(';')
        while s.peek()[1] != 'ENDSEC':
            _, name = s.take()
            number = int(name[1:])
            s.take('=')
            if s.peek()[1] == '(':
                s.take('(')
                parts = []
                while s.peek()[1] != ')':
                    parts.append(record(s))
                s.take(')')
            else:
                parts = [record(s)]
            s.take(';')
            if number in data:
                raise ValueError(f'#{number} defined twice')
            data[number] = parts
        s.take('ENDSEC')
        s.take(';')
    s.take('END-ISO-10303-21')
    s.take(';')
    return header, data


def part(data, number, name):
    """The parameters of the named part of instance `number`."""
    found = [p for n, p in data[number] if n == name]
    if len(found) != 1:
        raise ValueError(f'#{number} has no single {name}')
    return found[0]


def kinds(data, number):
    return {n for n, _ in data[number]}


# --- units ----------------------------------------------------------------

PREFIX = {None: Fraction(1), 'MILLI': Fraction(1, 1000), 'CENTI': Fraction(1, 100),
          'DECI': Fraction(1, 10), 'KILO': Fraction(1000), 'MICRO': Fraction(1, 10**6)}


def unit_factor(data, number, kind):
    """The factor from the unit to millimetres (length) or radians (plane
    angle), exact where the file's numbers are, as a Fraction or mpf."""
    names = kinds(data, number)
    if 'SI_UNIT' in names:
        prefix, base = part(data, number, 'SI_UNIT')
        prefix = None if prefix is None else str(prefix)
        if kind == 'length':
            assert str(base) == 'METRE'
            return PREFIX[prefix]*1000
        assert str(base) == 'RADIAN' and prefix is None
        return Fraction(1)
    if 'CONVERSION_BASED_UNIT' in names:
        _, measure = part(data, number, 'CONVERSION_BASED_UNIT')
        (name, (value, unit)), = [(n, p) for n, p in data[measure]
                                  if n.endswith('MEASURE_WITH_UNIT')]
        return Fraction(value.value)*unit_factor(data, unit, kind)
    raise ValueError(f'#{number} is not a unit')


def context(data, number):
    """(millimetres per length unit, radians per angle unit, uncertainty in
    millimetres or None) of a representation context."""
    units = part(data, number, 'GLOBAL_UNIT_ASSIGNED_CONTEXT')[0]
    length = angle = None
    for u in units:
        names = kinds(data, u)
        if 'LENGTH_UNIT' in names:
            length = unit_factor(data, u, 'length')
        elif 'PLANE_ANGLE_UNIT' in names:
            angle = unit_factor(data, u, 'angle')
    uncertainty = None
    if 'GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT' in kinds(data, number):
        for u in part(data, number, 'GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT')[0]:
            value, unit, _, _ = part(data, u, 'UNCERTAINTY_MEASURE_WITH_UNIT')
            if 'LENGTH_UNIT' in kinds(data, unit):
                uncertainty = Fraction(value.value)*unit_factor(data, unit, 'length')
    return length, angle, uncertainty


# --- bodies and OCCT's counts --------------------------------------------

BODIES = ('MANIFOLD_SOLID_BREP', 'BREP_WITH_VOIDS', 'SHELL_BASED_SURFACE_MODEL')
REPRESENTATIONS = ('ADVANCED_BREP_SHAPE_REPRESENTATION', 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION',
                   'SHAPE_REPRESENTATION')


def point(data, number):
    return [float(x) for x in part(data, number, 'CARTESIAN_POINT')[1]]


def surface(data, number):
    (name, params), = data[number]
    return name, params


def normalized(v):
    n = mpmath.sqrt(sum(mpmath.mpf(x)**2 for x in v))
    return [mpmath.mpf(x)/n for x in v]


def pole_passages(data, face):
    """Passages of the face's loops through a pole: a sphere's poles or a
    cone's apex, where OCCT's reader adds a degenerated edge."""
    _, bounds, geometry, _ = part(data, face, 'ADVANCED_FACE')
    name, params = surface(data, geometry)
    if name not in ('SPHERICAL_SURFACE', 'CONICAL_SURFACE'):
        return 0
    _, place, radius, *rest = params
    _, location, axis, _ = part(data, place, 'AXIS2_PLACEMENT_3D')
    o = [mpmath.mpf(x) for x in point(data, location)]
    a = normalized(part(data, axis, 'DIRECTION')[1])
    poles = []
    if name == 'SPHERICAL_SURFACE':
        poles = [[o[i]+s*radius*a[i] for i in range(3)] for s in (1, -1)]
    else:
        # The apex, for the semi-angle in any angle unit: where the radius
        # R + v tan(alpha) vanishes along the axis, found from the file's
        # own vertices (any vertex on the axis is the apex).
        poles = None
    count = 0
    for b in bounds:
        loop = next(p for n, p in data[b] if n in ('FACE_BOUND', 'FACE_OUTER_BOUND'))[1]
        if 'EDGE_LOOP' not in kinds(data, loop):
            continue
        for oe in part(data, loop, 'EDGE_LOOP')[1]:
            _, _, _, edge, orientation = part(data, oe, 'ORIENTED_EDGE')
            _, start, end, _, _ = part(data, edge, 'EDGE_CURVE')
            head = end if (str(orientation) == 'T') else start
            p = [mpmath.mpf(x) for x in point(data, part(data, head, 'VERTEX_POINT')[1])]
            if poles is None:
                rel = [p[i]-o[i] for i in range(3)]
                h = sum(rel[i]*a[i] for i in range(3))
                radial = mpmath.sqrt(max(mpmath.mpf(0), sum(r*r for r in rel)-h*h))
                if radial <= mpmath.mpf(10)**-9*(1+abs(radius)):
                    count += 1
            elif any(mpmath.sqrt(sum((p[i]-q[i])**2 for i in range(3))) <= mpmath.mpf(10)**-9*(1+radius)
                     for q in poles):
                count += 1
    return count


def shell_faces(data, shell):
    names = kinds(data, shell)
    if 'ORIENTED_CLOSED_SHELL' in names:
        return shell_faces(data, part(data, shell, 'ORIENTED_CLOSED_SHELL')[2])
    for kind in ('CLOSED_SHELL', 'OPEN_SHELL'):
        if kind in names:
            return part(data, shell, kind)[1]
    raise ValueError(f'#{shell} is not a shell')


def counts(data, shells, solid):
    """OCCT's counts (V E W F SH SO) of a body made of these shells."""
    vertices, edges, wires, faces, poles = set(), set(), 0, set(), 0
    for shell in shells:
        for face in shell_faces(data, shell):
            faces.add(face)
            _, bounds, _, _ = part(data, face, 'ADVANCED_FACE')
            poles += pole_passages(data, face)
            for b in bounds:
                wires += 1
                loop = next(p for n, p in data[b] if n in ('FACE_BOUND', 'FACE_OUTER_BOUND'))[1]
                if 'VERTEX_LOOP' in kinds(data, loop):
                    vertices.add(part(data, loop, 'VERTEX_LOOP')[1])
                    continue
                for oe in part(data, loop, 'EDGE_LOOP')[1]:
                    edge = part(data, oe, 'ORIENTED_EDGE')[3]
                    edges.add(edge)
                    _, start, end, _, _ = part(data, edge, 'EDGE_CURVE')
                    vertices.update([start, end])
    return [len(vertices), len(edges)+poles, wires, len(faces), len(shells), 1 if solid else 0]


def bodies(data):
    """[(entity, 'solid'|'sheet', counts, (mm per unit, rad per unit,
    uncertainty mm))] in entity-number order; a surface model gives one
    body per shell."""
    owner = {}
    for number in sorted(data):
        for name, params in data[number]:
            if name in REPRESENTATIONS:
                for item in params[1]:
                    owner.setdefault(item, params[2])
    out = []
    for number in sorted(data):
        names = kinds(data, number)
        units = context(data, owner[number]) if number in owner else None
        if 'MANIFOLD_SOLID_BREP' in names:
            out.append((number, 'solid', counts(data, [part(data, number, 'MANIFOLD_SOLID_BREP')[1]], True),
                        units))
        elif 'BREP_WITH_VOIDS' in names:
            _, outer, voids = part(data, number, 'BREP_WITH_VOIDS')
            out.append((number, 'solid', counts(data, [outer]+list(voids), True), units))
        elif 'SHELL_BASED_SURFACE_MODEL' in names:
            for shell in part(data, number, 'SHELL_BASED_SURFACE_MODEL')[1]:
                out.append((number, 'sheet', counts(data, [shell], False), units))
    return out


def schema(header):
    return next(p[0][0] for n, p in header if n == 'FILE_SCHEMA')


# --- closed forms -----------------------------------------------------------

PI = mpmath.pi


def mp(x):
    return mpmath.mpf(x) if not isinstance(x, Fraction) else mpmath.mpf(x.numerator)/x.denominator


def box(lo, hi):
    """(volume, area, moment vector) of an axis-aligned box."""
    d = [mp(hi[i])-mp(lo[i]) for i in range(3)]
    v = d[0]*d[1]*d[2]
    a = 2*(d[0]*d[1]+d[1]*d[2]+d[2]*d[0])
    return v, a, [v*(mp(lo[i])+d[i]/2) for i in range(3)]


def polygon_prism(points, z0, z1):
    """A prism of a counter-clockwise polygon between two heights."""
    n = len(points)
    area = cx = cy = perimeter = mpmath.mpf(0)
    for i in range(n):
        (x0, y0), (x1, y1) = [[mp(c) for c in points[i]], [mp(c) for c in points[(i+1) % n]]]
        cross = x0*y1-x1*y0
        area += cross/2
        cx += (x0+x1)*cross/6
        cy += (y0+y1)*cross/6
        perimeter += mpmath.sqrt((x1-x0)**2+(y1-y0)**2)
    h = mp(z1)-mp(z0)
    v = area*h
    return v, 2*area+perimeter*h, [cx*h, cy*h, v*(mp(z0)+h/2)]


def cylinder(o, axis, r, h):
    """A right circular cylinder from the base centre along a unit axis."""
    r, h = mp(r), mp(h)
    v = PI*r*r*h
    return v, 2*PI*r*r+2*PI*r*h, [v*(mp(o[i])+mp(axis[i])*h/2) for i in range(3)]


def hole(o, axis, r, h):
    """A cylindrical through-hole: removed with sign -1 by `combine`, its
    area term is the wall added and the two discs removed."""
    v, _, m = cylinder(o, axis, r, h)
    r, h = mp(r), mp(h)
    return v, 2*PI*r*h-2*PI*r*r, m


def frustum(o, r0, r1, h):
    """A frustum (or apex cone, r1 = 0) on the base at o along +z."""
    r0, r1, h = mp(r0), mp(r1), mp(h)
    v = PI*h*(r0*r0+r0*r1+r1*r1)/3
    slant = mpmath.sqrt(h*h+(r0-r1)**2)
    a = PI*(r0*r0+r1*r1)+PI*(r0+r1)*slant
    zc = h*(r0*r0+2*r0*r1+3*r1*r1)/(4*(r0*r0+r0*r1+r1*r1))
    return v, a, [v*mp(o[0]), v*mp(o[1]), v*(mp(o[2])+zc)]


def apex_cone_down(o, r, h):
    """A cone on the base at o (radius r) with its apex h above."""
    return frustum(o, r, 0, h)


def sphere(o, r):
    r = mp(r)
    v = 4*PI*r**3/3
    return v, 4*PI*r*r, [v*mp(c) for c in o]


def hemisphere(o, r):
    """The upper half of a ball about o."""
    r = mp(r)
    v = 2*PI*r**3/3
    return v, 3*PI*r*r, [v*mp(o[0]), v*mp(o[1]), v*(mp(o[2])+3*r/8)]


def torus(o, big, small):
    big, small = mp(big), mp(small)
    v = 2*PI**2*big*small**2
    return v, 4*PI**2*big*small, [v*mp(c) for c in o]


def elbow(o, big, small):
    """The quarter of a solid torus about +z between the +x and +y
    half-planes, with its two meridian discs."""
    big, small = mp(big), mp(small)
    v = PI**2*big*small**2/2
    a = PI**2*big*small+2*PI*small**2
    # The centroid of a disc swept by a quarter turn: along the bisector at
    # (2 sin(theta/2)/theta) (R^2 + r^2/4)/R.
    d = 2*mpmath.sin(PI/4)/(PI/2)*(big**2+small**2/4)/big
    c = d/mpmath.sqrt(2)
    return v, a, [v*(mp(o[0])+c), v*(mp(o[1])+c), v*mp(o[2])]


def half_cylinder(o, r, h):
    """The half y >= 0 of a cylinder about +z through o."""
    r, h = mp(r), mp(h)
    v = PI*r*r*h/2
    a = PI*r*r+PI*r*h+2*r*h
    return v, a, [v*mp(o[0]), v*(mp(o[1])+4*r/(3*PI)), v*(mp(o[2])+h/2)]


def tetrahedron(o, a):
    """The corner tetrahedron o, o + a x, o + a y, o + a z."""
    a = mp(a)
    v = a**3/6
    area = 3*a*a/2+mpmath.sqrt(3)*a*a/2
    return v, area, [v*(mp(c)+a/4) for c in o]


def open_box(lo, hi):
    """The five faces of a box without its top: area and centre of area."""
    d = [mp(hi[i])-mp(lo[i]) for i in range(3)]
    faces = [  # (area, centre)
        (d[0]*d[1], [mp(lo[0])+d[0]/2, mp(lo[1])+d[1]/2, mp(lo[2])]),
        (d[1]*d[2], [mp(lo[0]), mp(lo[1])+d[1]/2, mp(lo[2])+d[2]/2]),
        (d[1]*d[2], [mp(hi[0]), mp(lo[1])+d[1]/2, mp(lo[2])+d[2]/2]),
        (d[0]*d[2], [mp(lo[0])+d[0]/2, mp(lo[1]), mp(lo[2])+d[2]/2]),
        (d[0]*d[2], [mp(lo[0])+d[0]/2, mp(hi[1]), mp(lo[2])+d[2]/2]),
    ]
    a = sum(f[0] for f in faces)
    return a, [sum(f[0]*f[1][i] for f in faces)/a for i in range(3)]


def combine(*parts):
    """Signed sums (volume, area, moments) -> (volume, area, centre)."""
    v = sum(s*p[0] for s, p in parts)
    a = sum(abs(s)*p[1] for s, p in parts)
    m = [sum(s*p[2][i] for s, p in parts) for i in range(3)]
    return v, a, [x/v for x in m]


def rn(x):
    """Round an mpf once to binary64."""
    return float(mpmath.mpf(x))
