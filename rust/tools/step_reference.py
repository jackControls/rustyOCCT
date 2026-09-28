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
  60 digits; for STEP-b's B-spline shapes, exact rational arithmetic
  (`Fraction`) where the measure is a polynomial integral, and mpmath
  Gauss-Legendre quadrature per knot span at 40 digits (its error estimate
  below `1e-30`) for arc lengths and areas under square roots.
* `geometry_gaps` (STEP-b): its own evaluation of every curve and surface
  of a parsed file (lines, circles, ellipses, B-splines by de Boor's
  algorithm, rational ones in homogeneous coordinates, planes, cylinders,
  cones, spheres, tori) and the largest distance of each edge's vertices
  from its curve, of its curve from the analytic surfaces of its faces, and,
  on a B-spline surface, of its curve from the image of the file's pcurve at
  the same fraction of both ranges.
"""
from fractions import Fraction
import math
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
  | (?P<string>'(?:''|\\\\|\\S\\.|[^'])*')
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
    """(name, params) of a simple surface instance; a complex one (a
    rational B-spline surface) is named by its first record."""
    (name, params), *_ = data[number]
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


# --- STEP-b: exact polynomial pieces --------------------------------------

def binomial(n, k):
    return Fraction(math.comb(n, k))


def bezier_pieces(degree, knots, mults, ctrl):
    """A clamped polynomial B-spline (end multiplicities degree + 1, knots
    distinct) as Bézier control polygons per span, by exact knot insertion
    (Boehm) in Fractions: [(first, last, [control tuples])]."""
    t = [Fraction(k) for k, m in zip(knots, mults) for _ in range(m)]
    P = [tuple(Fraction(c) for c in q) for q in ctrl]
    if t.count(t[0]) != degree+1 or t.count(t[-1]) != degree+1:
        raise ValueError('not clamped')
    for k in [Fraction(x) for x in knots[1:-1]]:
        while t.count(k) < degree:
            s = max(i for i in range(len(t)-1) if t[i] <= k < t[i+1])
            Q = []
            for i in range(len(P)+1):
                if i <= s-degree:
                    Q.append(P[i])
                elif i > s:
                    Q.append(P[i-1])
                else:
                    a = (k-t[i])/(t[i+degree]-t[i])
                    Q.append(tuple((1-a)*P[i-1][j]+a*P[i][j] for j in range(len(P[i]))))
            P = Q
            t.insert(s+1, k)
    return [(Fraction(knots[j]), Fraction(knots[j+1]), P[j*degree:j*degree+degree+1])
            for j in range(len(knots)-1)]


class Poly:
    """A polynomial in one variable with Fraction coefficients, power basis."""

    def __init__(self, c):
        self.c = list(c) or [Fraction(0)]

    def __add__(self, o):
        n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0)+(o.c[i] if i < len(o.c) else 0)
                     for i in range(n)])

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly([x*o for x in self.c])
        out = [Fraction(0)]*(len(self.c)+len(o.c)-1)
        for i, a in enumerate(self.c):
            for j, b in enumerate(o.c):
                out[i+j] += a*b
        return Poly(out)

    def deriv(self):
        return Poly([i*self.c[i] for i in range(1, len(self.c))])

    def integral01(self):
        return sum(x/(i+1) for i, x in enumerate(self.c))

    def __call__(self, x):
        out = mpmath.mpf(0)
        for a in reversed(self.c):
            out = out*x+mp(a)
        return out


def power(values):
    """Bernstein coefficients on [0, 1] to a Poly."""
    n = len(values)-1
    return Poly([binomial(n, j)*sum((-1)**(j-i)*binomial(j, i)*Fraction(values[i]) for i in range(j+1))
                 for j in range(n+1)])


def bernstein(poly, n):
    """A Poly of degree at most n to its Bernstein coefficients on [0, 1]."""
    a = poly.c+[Fraction(0)]*(n+1-len(poly.c))
    return [sum(binomial(k, j)/binomial(n, j)*a[j] for j in range(k+1)) for k in range(n+1)]


def curve_polys(ctrl):
    """A Bézier piece's coordinates as Polys in its local parameter."""
    return [power([q[i] for q in ctrl]) for i in range(len(ctrl[0]))]


def patch_polys(grid):
    """A Bézier patch (grid[i][j], i along u) as {(a, b): coefficient} per
    coordinate: the sum of c u^a v^b."""
    m, n = len(grid)-1, len(grid[0])-1
    out = []
    for k in range(len(grid[0][0])):
        # Along v for every row, then along u for every power of v.
        rows = [power([grid[i][j][k] for j in range(n+1)]).c for i in range(m+1)]
        coeffs = {}
        for b in range(n+1):
            column = power([rows[i][b] if b < len(rows[i]) else 0 for i in range(m+1)]).c
            for a, x in enumerate(column):
                coeffs[(a, b)] = x
        out.append(coeffs)
    return out


def bivariate(coeffs, u, v, du=0, dv=0):
    """A patch coordinate, or its first partial in u or v, at (u, v)."""
    out = mpmath.mpf(0)
    for (a, b), x in coeffs.items():
        if a < du or b < dv or x == 0:
            continue
        out += mp(x)*(a if du else 1)*(b if dv else 1)*u**(a-du)*v**(b-dv)
    return out


def compose_patch(coeffs, u, v):
    """The Poly t -> coordinate(u(t), v(t)) for Polys u and v."""
    out = Poly([Fraction(0)])
    for (a, b), x in coeffs.items():
        term = Poly([x])
        for _ in range(a):
            term = term*u
        for _ in range(b):
            term = term*v
        out = out+term
    return out


# --- STEP-b: closed forms and quadratures ---------------------------------

QUAD_DPS = 40
QUAD_ERROR = mpmath.mpf(10)**-30


def quad(f, *intervals):
    """Gauss-Legendre quadrature at 40 digits, its error estimate bounded."""
    with mpmath.workdps(QUAD_DPS):
        value, error = mpmath.quad(f, *intervals, method='gauss-legendre', error=True)
    if error > QUAD_ERROR*(1+abs(value)):
        raise ValueError(f'quadrature error {error}')
    return value


def profile(degree, poles, knots, mults):
    """The counter-clockwise region bounded by a clamped planar B-spline and
    the segment from its end back to its start: (area, centroid x, centroid
    y, curve length, segment length); area and centroid exact."""
    pieces = bezier_pieces(degree, knots, mults, poles)
    a = mx = my = Fraction(0)

    def green(x, y):
        # Area and first moments by Green's theorem along a piece.
        dx, dy = x.deriv(), y.deriv()
        return ((x*dy+(y*dx)*-1).integral01()/2, (x*x*dy).integral01()/2, (y*y*dx).integral01()/-2)

    (x1, y1), (x0, y0) = poles[-1], poles[0]
    segment = [Poly([Fraction(x1), Fraction(x0)-Fraction(x1)]), Poly([Fraction(y1), Fraction(y0)-Fraction(y1)])]
    for x, y in [curve_polys(ctrl) for _, _, ctrl in pieces]+[segment]:
        da, dmx, dmy = green(x, y)
        a, mx, my = a+da, mx+dmx, my+dmy
    length = mpmath.mpf(0)
    for _, _, ctrl in pieces:
        dx, dy = [p.deriv() for p in curve_polys(ctrl)]
        length += quad(lambda s: mpmath.sqrt(dx(s)**2+dy(s)**2), [0, 1])
    chord = mpmath.sqrt(mp(Fraction(x1)-Fraction(x0))**2+mp(Fraction(y1)-Fraction(y0))**2)
    return a, mx/a, my/a, length, chord


def bspline_plate(prof, z):
    """The planar sheet of a profile at height z: area and centre."""
    a, cx, cy, _, _ = prof
    return mp(a), [mp(cx), mp(cy), mp(z)]


def bspline_prism(prof, z0, z1):
    """The prism of a profile between two heights."""
    a, cx, cy, length, chord = prof
    h = mp(z1)-mp(z0)
    return mp(a)*h, 2*mp(a)+(length+chord)*h, [mp(cx)*mp(a)*h, mp(cy)*mp(a)*h, (mp(z0)+h/2)*mp(a)*h]


def patch_measure(grid, g=None):
    """Area and moments of a Bézier patch over the unit square, or of its
    part v <= g(u) for a Poly g, by quadrature of |S_u x S_v|."""
    polys = patch_polys(grid)

    def normal(u, v):
        su = [bivariate(c, u, v, du=1) for c in polys]
        sv = [bivariate(c, u, v, dv=1) for c in polys]
        n = [su[1]*sv[2]-su[2]*sv[1], su[2]*sv[0]-su[0]*sv[2], su[0]*sv[1]-su[1]*sv[0]]
        return mpmath.sqrt(sum(x*x for x in n))

    def integrand(k):
        def f(u, w):
            v, jac = (w, 1) if g is None else (g(u)*w, g(u))
            weight = normal(u, v)*jac
            return weight if k is None else weight*bivariate(polys[k], u, v)
        return f

    return quad(integrand(None), [0, 1], [0, 1]), [quad(integrand(k), [0, 1], [0, 1]) for k in range(3)]


def spline_sheet(pieces):
    """Area and centre of a sheet of Bézier patches [(grid, g or None)]."""
    area, moments = mpmath.mpf(0), [mpmath.mpf(0)]*3
    for grid, g in pieces:
        a, m = patch_measure(grid, g)
        area += a
        moments = [moments[i]+m[i] for i in range(3)]
    return area, [x/area for x in moments]


def half_ellipse(c, y, a, b):
    """The half of an ellipse (centre c, semi-axes a and b) on the side of
    its unit second axis y: area and centre."""
    a, b = mp(a), mp(b)
    d = 4*b/(3*PI)
    return PI*a*b/2, [mp(c[i])+d*mp(y[i]) for i in range(3)]


def oblique_cylinder(r, h, slope):
    """The cylinder of radius r about +z through the origin between z = 0
    and the plane z = h + slope y."""
    r, h, s = mp(r), mp(h), mp(slope)
    v = PI*r*r*h
    area = PI*r*r+PI*r*r*mpmath.sqrt(1+s*s)+2*PI*r*h
    return v, area, [mpmath.mpf(0), v*s*r*r/(4*h), v*(h*h+s*s*r*r/4)/(2*h)]


# --- STEP-b: evaluation of a parsed file ------------------------------------

def mpv(v):
    return [mpmath.mpf(x) for x in v]


def vsub(a, b):
    return [a[i]-b[i] for i in range(len(a))]


def vdot(a, b):
    return sum(a[i]*b[i] for i in range(len(a)))


def vcross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def unit_mp(v):
    n = mpmath.sqrt(vdot(v, v))
    return [x/n for x in v]


def axes(data, number):
    """AXIS2_PLACEMENT_3D by ISO 10303-42's build_axes: (o, x, y, z)."""
    _, loc, axis, ref = part(data, number, 'AXIS2_PLACEMENT_3D')
    o = mpv(point(data, loc))
    z = unit_mp(mpv(part(data, axis, 'DIRECTION')[1])) if axis is not None else mpv([0, 0, 1])
    if ref is not None:
        v = mpv(part(data, ref, 'DIRECTION')[1])
    else:
        v = mpv([0, 1, 0]) if [abs(c) for c in z] == [1, 0, 0] else mpv([1, 0, 0])
    k = vdot(v, z)
    x = unit_mp([v[i]-k*z[i] for i in range(3)])
    return o, x, vcross(z, x), z


def full_knots(knots, mults):
    return [mpmath.mpf(k) for k, m in zip(knots, mults) for _ in range(int(m))]


def basis(degree, t_knots, n, t):
    """(first index, the degree + 1 nonzero basis values at t) by Cox-de Boor
    (The NURBS Book, A2.2); t at the domain's end takes the last span."""
    spans = [i for i in range(degree, n) if t_knots[i] < t_knots[i+1]]
    inside = [i for i in spans if t_knots[i] <= t]
    s = inside[-1] if inside else spans[0]
    N = [mpmath.mpf(1)]+[mpmath.mpf(0)]*degree
    left, right = [0]*(degree+1), [0]*(degree+1)
    for j in range(1, degree+1):
        left[j] = t-t_knots[s+1-j]
        right[j] = t_knots[s+j]-t
        saved = mpmath.mpf(0)
        for r in range(j):
            temp = N[r]/(right[r+1]+left[j-r])
            N[r] = saved+right[r+1]*temp
            saved = left[j-r]*temp
        N[j] = saved
    return s-degree, N


class BSplineCurve:
    def __init__(self, degree, poles, weights, knots, mults):
        self.p = int(degree)
        self.poles = [mpv(q) for q in poles]
        self.w = mpv(weights) if weights is not None else [mpmath.mpf(1)]*len(poles)
        self.t = full_knots(knots, mults)
        n = len(self.poles)
        if len(self.t) != n+self.p+1:
            raise ValueError('not a clamped or unclamped B-spline')
        self.domain = (self.t[self.p], self.t[n])
        self.periodic = False

    def point(self, t):
        i0, N = basis(self.p, self.t, len(self.poles), t)
        num, den = [mpmath.mpf(0)]*len(self.poles[0]), mpmath.mpf(0)
        for r, b in enumerate(N):
            w = self.w[i0+r]*b
            den += w
            num = [num[k]+w*self.poles[i0+r][k] for k in range(len(num))]
        return [x/den for x in num]


class BSplineSurface:
    def __init__(self, degrees, grid, weights, knots, mults):
        self.p = [int(d) for d in degrees]
        self.grid = [[mpv(q) for q in row] for row in grid]
        self.w = [mpv(row) for row in weights] if weights is not None else \
            [[mpmath.mpf(1)]*len(row) for row in grid]
        self.t = [full_knots(knots[k], mults[k]) for k in range(2)]
        self.n = [len(grid), len(grid[0])]
        self.domain = [(self.t[k][self.p[k]], self.t[k][self.n[k]]) for k in range(2)]

    def point(self, u, v):
        iu, Nu = basis(self.p[0], self.t[0], self.n[0], u)
        iv, Nv = basis(self.p[1], self.t[1], self.n[1], v)
        num, den = [mpmath.mpf(0)]*3, mpmath.mpf(0)
        for a, bu in enumerate(Nu):
            for b, bv in enumerate(Nv):
                w = self.w[iu+a][iv+b]*bu*bv
                den += w
                num = [num[k]+w*self.grid[iu+a][iv+b][k] for k in range(3)]
        return [x/den for x in num]


class Line:
    def __init__(self, o, d):
        self.o, self.d = o, unit_mp(d)
        self.domain, self.periodic = None, False

    def point(self, t):
        return [self.o[i]+t*self.d[i] for i in range(len(self.o))]

    def locate(self, p):
        return vdot(vsub(p, self.o), self.d)


class Conic:
    """A circle or ellipse: o + a1 cos t x + a2 sin t y."""

    def __init__(self, frame, a1, a2):
        self.o, self.x, self.y, _ = frame
        self.a1, self.a2 = mpmath.mpf(a1), mpmath.mpf(a2)
        self.domain, self.periodic = (mpmath.mpf(0), 2*PI), True

    def point(self, t):
        c, s = mpmath.cos(t), mpmath.sin(t)
        return [self.o[i]+self.a1*c*self.x[i]+self.a2*s*self.y[i] for i in range(3)]

    def locate(self, p):
        q = vsub(p, self.o)
        return mpmath.atan2(vdot(q, self.y)/self.a2, vdot(q, self.x)/self.a1)


def bspline_parts(data, number, kind):
    """(degree(s), control points, weights|None, knots, multiplicities) of a
    `B_SPLINE_{kind}_WITH_KNOTS`, simple or a complex rational instance."""
    names = kinds(data, number)
    if len(data[number]) == 1:
        params = part(data, number, f'B_SPLINE_{kind}_WITH_KNOTS')
        if kind == 'CURVE':
            _, degree, poles, _, _, _, mults, knots, _ = params
            return degree, poles, None, knots, mults
        _, du, dv, poles, _, _, _, _, um, vm, uk, vk, _ = params
        return (du, dv), poles, None, (uk, vk), (um, vm)
    head = part(data, number, f'B_SPLINE_{kind}')
    tail = part(data, number, f'B_SPLINE_{kind}_WITH_KNOTS')
    weights = part(data, number, f'RATIONAL_B_SPLINE_{kind}')[0] \
        if f'RATIONAL_B_SPLINE_{kind}' in names else None
    if kind == 'CURVE':
        degree, poles, _, _, _ = head
        mults, knots, _ = tail
        return degree, poles, weights, knots, mults
    du, dv, poles, _, _, _, _ = head
    um, vm, uk, vk, _ = tail
    return (du, dv), poles, weights, (uk, vk), (um, vm)


def curve_of(data, number):
    """(curve, pcurve entities, reversed) of an edge's geometry."""
    names = kinds(data, number)
    for kind in ('SURFACE_CURVE', 'SEAM_CURVE'):
        if kind in names:
            _, c3, associated, _ = part(data, number, kind)
            curve, _, reversed_ = curve_of(data, c3)
            return curve, [a for a in associated if 'PCURVE' in kinds(data, a)], reversed_
    if 'TRIMMED_CURVE' in names:
        _, basis_curve, _, _, sense, _ = part(data, number, 'TRIMMED_CURVE')
        curve, pcurves, reversed_ = curve_of(data, basis_curve)
        return curve, pcurves, reversed_ != (str(sense) == 'F')
    return plain_curve(data, number), [], False


def plain_curve(data, number):
    names = kinds(data, number)
    if 'LINE' in names:
        _, o, vec = part(data, number, 'LINE')
        return Line(mpv(point(data, o)), mpv(part(data, part(data, vec, 'VECTOR')[1], 'DIRECTION')[1]))
    if 'CIRCLE' in names:
        _, place, r = part(data, number, 'CIRCLE')
        return Conic(axes(data, place), r, r)
    if 'ELLIPSE' in names:
        _, place, a1, a2 = part(data, number, 'ELLIPSE')
        return Conic(axes(data, place), a1, a2)
    if 'B_SPLINE_CURVE_WITH_KNOTS' in names:
        degree, poles, weights, knots, mults = bspline_parts(data, number, 'CURVE')
        return BSplineCurve(degree, [point(data, q) for q in poles], weights, knots, mults)
    raise ValueError(f'#{number}: no curve the reference evaluates')


def analytic_gap(data, number, angle, p):
    """The distance of p from an elementary surface; None for a spline."""
    if 'B_SPLINE_SURFACE_WITH_KNOTS' in kinds(data, number):
        return None
    name, params = surface(data, number)
    o, x, y, z = axes(data, params[1])
    q = vsub(p, o)
    h = vdot(q, z)
    radial = mpmath.sqrt(max(mpmath.mpf(0), vdot(q, q)-h*h))
    if name == 'PLANE':
        return abs(h)
    if name == 'CYLINDRICAL_SURFACE':
        return abs(radial-params[2])
    if name == 'CONICAL_SURFACE':
        a = mpmath.mpf(params[3])*mp(angle)
        return abs(radial-(params[2]+h*mpmath.tan(a)))*mpmath.cos(a)
    if name == 'SPHERICAL_SURFACE':
        return abs(mpmath.sqrt(vdot(q, q))-params[2])
    if name == 'TOROIDAL_SURFACE':
        return abs(mpmath.sqrt((radial-params[2])**2+h*h)-params[3])
    raise ValueError(f'#{number}: no surface the reference evaluates')


def spline_surface(data, number):
    degrees, poles, weights, knots, mults = bspline_parts(data, number, 'SURFACE')
    grid = [[point(data, q) for q in row] for row in poles]
    return BSplineSurface(degrees, grid, weights, knots, mults)


def locate(f, lo, hi, target, samples=256):
    """The parameter in [lo, hi] where f comes nearest the target: the best
    sample, then golden-section search between its neighbours."""
    def d(t):
        return sum(x*x for x in vsub(f(t), target))
    ts = [lo+(hi-lo)*k/samples for k in range(samples+1)]
    k = min(range(len(ts)), key=lambda i: d(ts[i]))
    a, b = ts[max(k-1, 0)], ts[min(k+1, samples)]
    g = (mpmath.sqrt(5)-1)/2
    for _ in range(160):
        c, e = b-g*(b-a), a+g*(b-a)
        if d(c) < d(e):
            b = e
        else:
            a = c
    return (a+b)/2


def curve_range(curve, start, end, closed):
    """The parameters of an edge's vertices on its curve, running along it."""
    if isinstance(curve, BSplineCurve):
        lo, hi = curve.domain
        if closed:
            return lo, hi
        return locate(curve.point, lo, hi, start), locate(curve.point, lo, hi, end)
    s, e = curve.locate(start), curve.locate(end)
    if curve.periodic:
        e = s+2*PI if closed else s+(e-s) % (2*PI)
    return s, e


def pcurve_curve(data, pcurve):
    """(basis surface, 2D curve) of a PCURVE."""
    _, basis_surface, rep = part(data, pcurve, 'PCURVE')
    items = part(data, rep, 'DEFINITIONAL_REPRESENTATION')[1]
    return basis_surface, plain_curve(data, items[0])


def pcurve_range(surf, c2, start, end):
    """The parameters of the edge's vertices on the image of a pcurve."""
    def f(t):
        return surf.point(*c2.point(t))
    if isinstance(c2, BSplineCurve):
        lo, hi = c2.domain
    else:
        # A line's part inside the surface's domain.
        ends = []
        for k in range(2):
            a, b = surf.domain[k]
            if c2.d[k] != 0:
                ends.append(sorted([(a-c2.o[k])/c2.d[k], (b-c2.o[k])/c2.d[k]]))
        lo, hi = max(e[0] for e in ends), min(e[1] for e in ends)
    return locate(f, lo, hi, start), locate(f, lo, hi, end)


def geometry_gaps(data, angle=1):
    """The largest gap of any edge: of its vertices from its curve, of its
    curve from each face's elementary surface, and from the image of its
    pcurve on each face's spline surface at the same fraction of both
    ranges; with the number of spline uses checked through a pcurve."""
    with mpmath.workdps(30):
        faces_of = {}
        for number in sorted(data):
            if 'ADVANCED_FACE' not in kinds(data, number):
                continue
            _, bounds, geometry, _ = part(data, number, 'ADVANCED_FACE')
            for b in bounds:
                loop = next(p for n, p in data[b] if n in ('FACE_BOUND', 'FACE_OUTER_BOUND'))[1]
                if 'EDGE_LOOP' not in kinds(data, loop):
                    continue
                for oe in part(data, loop, 'EDGE_LOOP')[1]:
                    faces_of.setdefault(part(data, oe, 'ORIENTED_EDGE')[3], set()).add(geometry)
        worst, spline_uses = mpmath.mpf(0), 0
        fractions = [mpmath.mpf(k)/8 for k in range(9)]

        def gap(p, q):
            return mpmath.sqrt(sum(x*x for x in vsub(p, q)))
        for edge in sorted(faces_of):
            _, v1, v2, geometry, same = part(data, edge, 'EDGE_CURVE')
            curve, pcurves, reversed_ = curve_of(data, geometry)
            if (str(same) == 'T') == reversed_:
                v1, v2 = v2, v1
            start, end = [mpv(point(data, part(data, v, 'VERTEX_POINT')[1])) for v in (v1, v2)]
            s, e = curve_range(curve, start, end, v1 == v2)
            worst = max(worst, gap(curve.point(s), start), gap(curve.point(e), end))
            samples = [curve.point(s+f*(e-s)) for f in fractions]
            for surf_number in sorted(faces_of[edge]):
                if analytic_gap(data, surf_number, angle, samples[0]) is not None:
                    worst = max([worst]+[analytic_gap(data, surf_number, angle, p) for p in samples])
                    continue
                surf = spline_surface(data, surf_number)
                found = [c2 for basis_surface, c2 in (pcurve_curve(data, p) for p in pcurves)
                         if basis_surface == surf_number]
                if not found:
                    raise ValueError(f'edge #{edge} has no pcurve on #{surf_number}')
                a, b = pcurve_range(surf, found[0], start, end)
                spline_uses += 1
                for f, p in zip(fractions, samples):
                    worst = max(worst, gap(surf.point(*found[0].point(a+f*(b-a))), p))
        return worst, spline_uses
