#!/usr/bin/env python3
"""Fixtures for STEP-a and STEP-b of REVIEW_NOTES.md: small STEP files
authored here.

Each case is built from its construction parameters by the B-rep builders
below (vertices, edges on lines and circles, faces on planes, cylinders,
cones, spheres and tori, loops oriented as ISO 10303-42 requires: the face
on the left seen from the face's normal) and written by the Part 21 writer
below in the style of an AP214 file (product structure, an advanced B-rep
or manifold surface shape representation, a context with units and an
uncertainty); one file, `syntax.stp`, is written by hand to exercise the
syntax (comments, complex instances, forward references, entity numbers out
of order, `$` and `*`, escaped strings, unnormalised directions). Variants
exercise the orientation flags (`same_sense` of faces and edge curves,
bound orientations), void shells, surface models, several bodies, lengths
in metres and inches and angles in degrees.

STEP-b's cases (`STEP_B`) add ellipses, B-spline curves and surfaces
(rational ones as complex instances, as OCCT writes them) and edges as
`SURFACE_CURVE`s carrying `PCURVE`s on spline surfaces (and, on the prism,
on planes): a half-ellipse sheet, a cylinder cut by an oblique plane, a
sheet and a prism bounded by a two-span cubic, a two-span patch sheet, a
patch trimmed by a curve whose pcurve is a quadratic B-spline (its 3D curve
the exact composition, of degree 9), and a cylinder of rational
half-circles. Every file's geometry is checked by the reference
(`step_reference.geometry_gaps`): each edge meets its vertices and lies on
its faces' surfaces, through the file's pcurve on a spline surface at the
same fraction of both ranges, within `1e-12` of the case's size.

`rust/fixtures/step/NAME.stp` holds each file and `step-expected.tsv` each
body: its entity number, class, OCCT's counts (`step_reference.bodies`,
from the independent parser) and its volume, area and centre from the
closed forms of `step_reference.py` (a sheet has no volume: `-`), exact
rational integrals or its quadrature. `--check` regenerates everything and
fails on any difference. No Rust or native result supplies an expectation.
"""
import argparse
from fractions import Fraction
from pathlib import Path

import mpmath

import step_reference as ref

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'fixtures/step'
EXPECTED = ROOT/'fixtures/step-expected.tsv'


def real(x):
    """A Part 21 real: digits, a point, an optional exponent."""
    r = repr(float(x))
    if r in ('inf', '-inf', 'nan'):
        raise ValueError('not finite')
    mantissa, _, exponent = r.partition('e')
    if '.' not in mantissa:
        mantissa += '.'
    if mantissa.endswith('.0'):
        mantissa = mantissa[:-1]
    return mantissa+('E'+exponent if exponent else '')


def add(a, b):
    return tuple(a[i]+b[i] for i in range(3))


def scale(a, s):
    return tuple(a[i]*s for i in range(3))


def neg(a):
    return tuple(-x for x in a)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


# --- the B-rep model --------------------------------------------------------

class Vertex:
    def __init__(self, p):
        self.p = tuple(float(c) for c in p)


class Edge:
    """A line from `start` to `end`, or a circle (`frame`: centre, axis, x
    direction; radius) traversed counter-clockwise about its axis from
    `start` to `end` (the same vertex for a closed circle), or another
    `curve` (STEP-b): ('ELLIPSE', frame, a1, a2) or ('BSPLINE', degree,
    poles, weights | None, knots, multiplicities), run along its parameter.
    `pcurves`: [(face, 2D curve)], each ('LINE', point, direction) or a 2D
    ('BSPLINE', ...), written with the curve as a `SURFACE_CURVE`."""

    def __init__(self, start, end, circle=None, curve=None):
        self.start, self.end, self.circle, self.curve = start, end, circle, curve
        self.pcurves = []


class Face:
    """`surface`: ('PLANE', frame) | ('CYLINDRICAL_SURFACE', frame, r) |
    ('CONICAL_SURFACE', frame, r, semi_angle) | ('SPHERICAL_SURFACE', frame,
    r) | ('TOROIDAL_SURFACE', frame, R, r), a frame being (origin, axis, x),
    or ('B_SPLINE_SURFACE', (du, dv), rows of poles in u, weights | None,
    (u knots, v knots), (u multiplicities, v multiplicities));
    `same_sense`: whether the face's normal is the surface's; `bounds`: lists
    of (edge, forward) in order, the first the outer bound, each counter-
    clockwise about the face's normal (inner bounds clockwise)."""

    def __init__(self, surface, same_sense, bounds):
        self.surface, self.same_sense, self.bounds = surface, same_sense, bounds


XY = ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))


def frame_at(origin, axis=(0.0, 0.0, 1.0), x=(1.0, 0.0, 0.0)):
    return (tuple(float(c) for c in origin), axis, x)


def prism(points, z0, z1, holes=()):
    """Faces of the prism of a counter-clockwise polygon in the xy-plane,
    with circular through-holes (centre x, y, radius): each hole's wall is
    a cylinder whose face normal points into the hole (`same_sense` false)."""
    n = len(points)
    bottom = [Vertex((x, y, z0)) for x, y in points]
    top = [Vertex((x, y, z1)) for x, y in points]
    be = [Edge(bottom[i], bottom[(i+1) % n]) for i in range(n)]
    te = [Edge(top[i], top[(i+1) % n]) for i in range(n)]
    ve = [Edge(bottom[i], top[i]) for i in range(n)]
    bottom_bounds = [[(be[i], False) for i in reversed(range(n))]]
    top_bounds = [[(te[i], True) for i in range(n)]]
    faces = []
    walls = []
    for cx, cy, r in holes:
        a, b = Vertex((cx+r, cy, z0)), Vertex((cx+r, cy, z1))
        cb = Edge(a, a, (frame_at((cx, cy, z0)), r))
        ct = Edge(b, b, (frame_at((cx, cy, z1)), r))
        seam = Edge(a, b)
        bottom_bounds.append([(cb, True)])
        top_bounds.append([(ct, False)])
        walls.append(Face(('CYLINDRICAL_SURFACE', frame_at((cx, cy, z0)), r), False,
                          [[(cb, False), (seam, True), (ct, True), (seam, False)]]))
    faces.append(Face(('PLANE', frame_at((0.0, 0.0, z0), (0.0, 0.0, -1.0))), True, bottom_bounds))
    faces.append(Face(('PLANE', frame_at((0.0, 0.0, z1))), True, top_bounds))
    for i in range(n):
        (x0, y0), (x1, y1) = points[i], points[(i+1) % n]
        d = (x1-x0, y1-y0)
        length = max(abs(d[0]), abs(d[1]))
        # Axis-aligned or Pythagorean sides only: exact unit directions.
        normal = (d[1]/length, -d[0]/length, 0.0)
        faces.append(Face(('PLANE', frame_at((x0, y0, z0), normal, (d[0]/length, d[1]/length, 0.0))), True,
                          [[(be[i], True), (ve[(i+1) % n], True), (te[i], False), (ve[i], False)]]))
    return faces+walls


def box_faces(lo, hi):
    return prism([(lo[0], lo[1]), (hi[0], lo[1]), (hi[0], hi[1]), (lo[0], hi[1])], lo[2], hi[2])


def cylinder_faces(o, axis, x, r, h):
    y = cross(axis, x)
    a = Vertex(add(o, scale(x, r)))
    top = add(o, scale(axis, h))
    b = Vertex(add(top, scale(x, r)))
    cb = Edge(a, a, ((o, axis, x), r))
    ct = Edge(b, b, ((top, axis, x), r))
    seam = Edge(a, b)
    del y
    return [
        Face(('PLANE', (o, neg(axis), x)), True, [[(cb, False)]]),
        Face(('PLANE', (top, axis, x)), True, [[(ct, True)]]),
        Face(('CYLINDRICAL_SURFACE', (o, axis, x), r), True,
             [[(cb, True), (seam, True), (ct, False), (seam, False)]]),
    ]


def apex_cone_faces(o, r, h, semi):
    """A cone on the base at o (radius r, normal -z) with its apex h above:
    the conical surface about -z at the base, as the cone narrows upward."""
    a = Vertex((o[0]+r, o[1], o[2]))
    apex = Vertex((o[0], o[1], o[2]+h))
    cb = Edge(a, a, (frame_at(o), r))
    seam = Edge(a, apex)
    return [
        Face(('PLANE', frame_at(o, (0.0, 0.0, -1.0))), True, [[(cb, False)]]),
        Face(('CONICAL_SURFACE', frame_at(o, (0.0, 0.0, -1.0)), r, semi), True,
             [[(cb, True), (seam, True), (seam, False)]]),
    ]


def frustum_faces(o, r0, r1, h, semi):
    """A frustum widening upward: radius r0 at o, r1 at h above."""
    a = Vertex((o[0]+r0, o[1], o[2]))
    b = Vertex((o[0]+r1, o[1], o[2]+h))
    top = (o[0], o[1], o[2]+h)
    cb = Edge(a, a, (frame_at(o), r0))
    ct = Edge(b, b, (frame_at(top), r1))
    seam = Edge(a, b)
    return [
        Face(('PLANE', frame_at(o, (0.0, 0.0, -1.0))), True, [[(cb, False)]]),
        Face(('PLANE', frame_at(top)), True, [[(ct, True)]]),
        Face(('CONICAL_SURFACE', frame_at(o), r0, semi), True,
             [[(cb, True), (seam, True), (ct, False), (seam, False)]]),
    ]


# The meridian plane through +x: a circle about -y from x towards +z.
MERIDIAN_X = ((0.0, -1.0, 0.0), (1.0, 0.0, 0.0))


def sphere_faces(o, r):
    s, n = Vertex((o[0], o[1], o[2]-r)), Vertex((o[0], o[1], o[2]+r))
    seam = Edge(s, n, ((o,)+MERIDIAN_X, r))
    return [Face(('SPHERICAL_SURFACE', frame_at(o), r), True, [[(seam, True), (seam, False)]])]


def hemisphere_faces(o, r):
    e, n = Vertex((o[0]+r, o[1], o[2])), Vertex((o[0], o[1], o[2]+r))
    equator = Edge(e, e, (frame_at(o), r))
    seam = Edge(e, n, ((o,)+MERIDIAN_X, r))
    return [
        Face(('PLANE', frame_at(o, (0.0, 0.0, -1.0))), True, [[(equator, False)]]),
        Face(('SPHERICAL_SURFACE', frame_at(o), r), True, [[(equator, True), (seam, True), (seam, False)]]),
    ]


def torus_faces(o, big, small):
    v = Vertex((o[0]+big+small, o[1], o[2]))
    parallel = Edge(v, v, (frame_at(o), big+small))
    meridian = Edge(v, v, (((o[0]+big, o[1], o[2]),)+MERIDIAN_X, small))
    return [Face(('TOROIDAL_SURFACE', frame_at(o), big, small), True,
                 [[(parallel, True), (meridian, True), (parallel, False), (meridian, False)]])]


def elbow_faces(o, big, small):
    """A quarter of a solid torus about +z, from the +x to the +y half-plane."""
    v0 = Vertex((o[0]+big+small, o[1], o[2]))
    v1 = Vertex((o[0], o[1]+big+small, o[2]))
    m0 = Edge(v0, v0, (((o[0]+big, o[1], o[2]),)+MERIDIAN_X, small))
    m1 = Edge(v1, v1, (((o[0], o[1]+big, o[2]), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)), small))
    arc = Edge(v0, v1, (frame_at(o), big+small))
    return [
        Face(('TOROIDAL_SURFACE', frame_at(o), big, small), True,
             [[(arc, True), (m1, True), (arc, False), (m0, False)]]),
        Face(('PLANE', ((o[0]+big, o[1], o[2]),)+MERIDIAN_X), True, [[(m0, True)]]),
        Face(('PLANE', ((o[0], o[1]+big, o[2]), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0))), True, [[(m1, False)]]),
    ]


def half_cylinder_faces(o, r, h):
    a0, b0 = Vertex((o[0]+r, o[1], o[2])), Vertex((o[0]-r, o[1], o[2]))
    a1, b1 = Vertex((o[0]+r, o[1], o[2]+h)), Vertex((o[0]-r, o[1], o[2]+h))
    top = (o[0], o[1], o[2]+h)
    ab = Edge(a0, b0, (frame_at(o), r))
    at = Edge(a1, b1, (frame_at(top), r))
    lb, lt = Edge(b0, a0), Edge(b1, a1)
    va, vb = Edge(a0, a1), Edge(b0, b1)
    return [
        Face(('CYLINDRICAL_SURFACE', frame_at(o), r), True, [[(ab, True), (vb, True), (at, False), (va, False)]]),
        Face(('PLANE', frame_at(o, (0.0, 0.0, -1.0))), True, [[(ab, False), (lb, False)]]),
        Face(('PLANE', frame_at(top)), True, [[(lt, True), (at, True)]]),
        Face(('PLANE', frame_at(o, (0.0, -1.0, 0.0))), True, [[(lb, True), (va, True), (lt, False), (vb, False)]]),
    ]


# --- STEP-b's shapes ----------------------------------------------------------

def f3(p):
    return tuple(float(c) for c in p)


def half_ellipse_faces(c, x, a, b):
    """The half on the side of y = z x x of an ellipse in the plane z =
    c[2] sharing the plane's axes: an arc from c + a x to c - a x and the
    segment back."""
    z = (0, 0, 1)
    frame = (f3(c), f3(z), f3(x))
    v0 = Vertex([c[i]+a*x[i] for i in range(3)])
    v1 = Vertex([c[i]-a*x[i] for i in range(3)])
    arc = Edge(v0, v1, curve=('ELLIPSE', frame, a, b))
    return [Face(('PLANE', frame), True, [[(arc, True), (Edge(v1, v0), True)]])]


def oblique_cylinder_faces(r, h, normal, major):
    """The cylinder of radius r about +z through the origin between z = 0
    and the plane through (0, 0, h) of unit normal (0, -s, c): its section
    the ellipse of semi-axes r / c along `major` = (0, c, s) and r along
    -x, from the vertex (r, 0, h)."""
    a, b = Vertex((r, 0, 0)), Vertex((r, 0, h))
    top = (0.0, 0.0, float(h))
    circle = Edge(a, a, (frame_at((0, 0, 0)), r))
    ellipse = Edge(b, b, curve=('ELLIPSE', (top, f3(normal), f3(major)), Fraction(r)/major[1], r))
    seam = Edge(a, b)
    return [
        Face(('PLANE', frame_at((0, 0, 0), (0.0, 0.0, -1.0))), True, [[(circle, False)]]),
        Face(('PLANE', (top, f3(normal), (1.0, 0.0, 0.0))), True, [[(ellipse, True)]]),
        Face(('CYLINDRICAL_SURFACE', frame_at((0, 0, 0)), r), True,
             [[(circle, True), (seam, True), (ellipse, False), (seam, False)]]),
    ]


# The profile of the B-spline plate and prism: a two-span cubic from
# (10, 0) round to the origin above the x axis, closed by the segment.
PROFILE = (3, [(10, 0), (11, 5), (5, 10), (-1, 5), (0, 0)], (0, 1, 2), (4, 1, 4))


def profile_curve(z):
    degree, poles, knots, mults = PROFILE
    return ('BSPLINE', degree, [(x, y, z) for x, y in poles], None, knots, mults)


def bspline_plate_faces(z):
    v0, v1 = Vertex((0, 0, z)), Vertex((10, 0, z))
    segment, curve = Edge(v0, v1), Edge(v1, v0, curve=profile_curve(z))
    return [Face(('PLANE', frame_at((0, 0, z))), True, [[(segment, True), (curve, True)]])]


def bspline_prism_faces(h):
    """The prism of PROFILE from 0 to h: its side along the curve is the
    B-spline surface C(u) + v h z with the file's pcurves (lines), and its
    curves carry pcurves on the caps too."""
    degree, poles, knots, mults = PROFILE
    v0b, v1b, v0t, v1t = Vertex((0, 0, 0)), Vertex((10, 0, 0)), Vertex((0, 0, h)), Vertex((10, 0, h))
    lb, lt = Edge(v0b, v1b), Edge(v0t, v1t)
    cb, ct = Edge(v1b, v0b, curve=profile_curve(0)), Edge(v1t, v0t, curve=profile_curve(h))
    up0, up1 = Edge(v1b, v1t), Edge(v0b, v0t)
    side = ('B_SPLINE_SURFACE', (degree, 1), [[(x, y, 0), (x, y, h)] for x, y in poles], None,
            (knots, (0, 1)), (mults, (2, 2)))
    bottom = Face(('PLANE', frame_at((0, 0, 0), (0.0, 0.0, -1.0))), True, [[(cb, False), (lb, False)]])
    top = Face(('PLANE', frame_at((0, 0, h))), True, [[(lt, True), (ct, True)]])
    wall = Face(('PLANE', ((0.0, 0.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0))), True,
                [[(lb, True), (up0, True), (lt, False), (up1, False)]])
    spline = Face(side, True, [[(cb, True), (up1, True), (ct, False), (up0, False)]])
    # The bottom plane's coordinates are (x, -y).
    cb.pcurves = [(spline, ('LINE', (0, 0), (1, 0))),
                  (bottom, ('BSPLINE', degree, [(x, -y) for x, y in poles], None, knots, mults))]
    ct.pcurves = [(spline, ('LINE', (0, 1), (1, 0))),
                  (top, ('BSPLINE', degree, poles, None, knots, mults))]
    up0.pcurves = [(spline, ('LINE', (0, 0), (0, 1)))]
    up1.pcurves = [(spline, ('LINE', (2, 0), (0, 1)))]
    return [bottom, top, wall, spline]


# A two-span cubic by quadratic patch: x and y linear in u and v (the poles
# at the Greville abscissae, x = 6u, y = 8v), z varying.
PATCH = ((3, 2), [[(x, y, z) for y, z in zip((0, 4, 8), zs)]
                  for x, zs in zip((0, 2, 6, 10, 12), ([0, 1, 0], [1, 2, 1], [2, 4, 1], [1, 2, 2], [0, 1, 0]))],
         ((0, 1, 2), (0, 1)), ((4, 1, 4), (3, 3)))


def patch_faces():
    """The whole patch, bounded by its four boundary curves: lines as the
    pcurves along u, degree-1 B-splines along v."""
    (du, dv), rows, (uk, vk), (um, vm) = PATCH
    corner = {(i, j): Vertex(rows[-i][-j]) for i in (0, 1) for j in (0, 1)}
    column = lambda i: ('BSPLINE', dv, rows[i], None, vk, vm)
    row = lambda j: ('BSPLINE', du, [r[j] for r in rows], None, uk, um)
    bottom = Edge(corner[0, 0], corner[1, 0], curve=row(0))
    right = Edge(corner[1, 0], corner[1, 1], curve=column(-1))
    top = Edge(corner[0, 1], corner[1, 1], curve=row(-1))
    left = Edge(corner[0, 0], corner[0, 1], curve=column(0))
    face = Face(('B_SPLINE_SURFACE', (du, dv), rows, None, (uk, vk), (um, vm)), True,
                [[(bottom, True), (right, True), (top, False), (left, False)]])
    segment = lambda u: ('BSPLINE', 1, [(u, 0), (u, 1)], None, (0, 1), (2, 2))
    bottom.pcurves = [(face, ('LINE', (0, 0), (1, 0)))]
    top.pcurves = [(face, ('LINE', (0, 1), (1, 0)))]
    right.pcurves = [(face, segment(2))]
    left.pcurves = [(face, segment(0))]
    return [face]


def patch_pieces():
    """PATCH as Bézier patches (per u span) for the reference."""
    (du, dv), rows, (uk, vk), (um, vm) = PATCH
    columns = [ref.bezier_pieces(du, uk, um, [r[j] for r in rows]) for j in range(len(rows[0]))]
    return [([[columns[j][s][2][i] for j in range(len(rows[0]))] for i in range(du+1)], None)
            for s in range(len(uk)-1)]


# A single bicubic Bézier patch over [0, 9] x [0, 6] and the trim v <= g(u)
# = 1/2 + u (1 - u) / 2, a quadratic in (u, v).
TRIM_GRID = [[(3*i, 2*j, z) for j, z in enumerate(zs)]
             for i, zs in enumerate(([0, 1, 1, 0], [1, 2, 2, 1], [0, 2, 3, 1], [0, 1, 1, 0]))]
TRIM = [Fraction(1, 2), Fraction(3, 4), Fraction(1, 2)]


def trimmed_faces():
    polys = ref.patch_polys(TRIM_GRID)
    g = ref.power(TRIM)

    def at(u, v):
        return [ref.compose_patch(c, ref.Poly([Fraction(u)]), ref.Poly([Fraction(v)])).c[0] for c in polys]
    a, b, c, d = Vertex(at(0, 0)), Vertex(at(1, 0)), Vertex(at(1, TRIM[0])), Vertex(at(0, TRIM[0]))
    column = lambda i: ('BSPLINE', 3, TRIM_GRID[i], None, (0, 1), (4, 4))
    # The trim's 3D curve: S(t, g(t)), a polynomial of degree 9.
    lifted = [ref.bernstein(ref.compose_patch(k, ref.Poly([0, 1]), g), 9) for k in polys]
    trim = ('BSPLINE', 9, list(zip(*lifted)), None, (0, 1), (10, 10))
    bottom = Edge(a, b, curve=('BSPLINE', 3, [r[0] for r in TRIM_GRID], None, (0, 1), (4, 4)))
    right, top, left = Edge(b, c, curve=column(3)), Edge(d, c, curve=trim), Edge(a, d, curve=column(0))
    face = Face(('B_SPLINE_SURFACE', (3, 3), TRIM_GRID, None, ((0, 1), (0, 1)), ((4, 4), (4, 4))), True,
                [[(bottom, True), (right, True), (top, False), (left, False)]])
    bottom.pcurves = [(face, ('LINE', (0, 0), (1, 0)))]
    right.pcurves = [(face, ('LINE', (1, 0), (0, 1)))]
    top.pcurves = [(face, ('BSPLINE', 2, [(Fraction(k, 2), TRIM[k]) for k in range(3)], None, (0, 1), (3, 3)))]
    left.pcurves = [(face, ('LINE', (0, 0), (0, 1)))]
    return [face]


# A quarter-arc half circle as a rational quadratic: poles on the square
# about the circle, weights 1 and sqrt(2)/2 (rounded once).
HALF_WEIGHT = float(mpmath.sqrt(2)/2)


def rational_cylinder_faces(r, h):
    """The cylinder of radius r about +z from 0 to h made of two rational
    half-cylinders (y >= 0 and y <= 0) and two caps bounded by rational
    half circles, every curve and surface a complex instance."""
    halves = [[(r, 0), (r, r), (0, r), (-r, r), (-r, 0)], [(-r, 0), (-r, -r), (0, -r), (r, -r), (r, 0)]]
    weights = [1, HALF_WEIGHT, 1, HALF_WEIGHT, 1]
    arc = lambda poles, z: ('BSPLINE', 2, [(x, y, z) for x, y in poles], weights, (0, 1, 2), (3, 2, 3))
    p, q, pt, qt = Vertex((r, 0, 0)), Vertex((-r, 0, 0)), Vertex((r, 0, h)), Vertex((-r, 0, h))
    ab, bb = Edge(p, q, curve=arc(halves[0], 0)), Edge(q, p, curve=arc(halves[1], 0))
    at, bt = Edge(pt, qt, curve=arc(halves[0], h)), Edge(qt, pt, curve=arc(halves[1], h))
    up_p, up_q = Edge(p, pt), Edge(q, qt)
    wall = lambda poles: ('B_SPLINE_SURFACE', (2, 1), [[(x, y, 0), (x, y, h)] for x, y in poles],
                          [[w, w] for w in weights], ((0, 1, 2), (0, 1)), ((3, 2, 3), (2, 2)))
    fa = Face(wall(halves[0]), True, [[(ab, True), (up_q, True), (at, False), (up_p, False)]])
    fb = Face(wall(halves[1]), True, [[(bb, True), (up_p, True), (bt, False), (up_q, False)]])
    bottom = Face(('PLANE', frame_at((0, 0, 0), (0.0, 0.0, -1.0))), True, [[(ab, False), (bb, False)]])
    top = Face(('PLANE', frame_at((0, 0, h))), True, [[(at, True), (bt, True)]])
    line = lambda u, v, du, dv: ('LINE', (u, v), (du, dv))
    ab.pcurves, at.pcurves = [(fa, line(0, 0, 1, 0))], [(fa, line(0, 1, 1, 0))]
    bb.pcurves, bt.pcurves = [(fb, line(0, 0, 1, 0))], [(fb, line(0, 1, 1, 0))]
    up_q.pcurves = [(fa, line(2, 0, 0, 1)), (fb, line(0, 0, 0, 1))]
    up_p.pcurves = [(fa, line(0, 0, 0, 1)), (fb, line(2, 0, 0, 1))]
    return [fa, fb, bottom, top]


# --- the Part 21 writer ------------------------------------------------------

LENGTHS = {
    'mm': ['( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )'],
    'm': ['( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )'],
    'inch': ['( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )',
             'DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)',
             "( CONVERSION_BASED_UNIT('INCH',#{m}) LENGTH_UNIT() NAMED_UNIT(#{e}) )",
             'LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#{base})'],
}
SCHEMAS = {
    'ap214': ("'AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'", 'automotive_design'),
    'ap203': ("'CONFIG_CONTROL_DESIGN'", 'config_control_design'),
    'ap242': ("'AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'",
              'ap242_managed_model_based_3d_engineering'),
}


class Writer:
    def __init__(self, flip_faces=False, flip_edges=False, flip_bounds=False):
        self.lines = []
        self.flip_faces, self.flip_edges, self.flip_bounds = flip_faces, flip_edges, flip_bounds
        self.vertices, self.edges = {}, {}
        # STEP-b: surfaces written before their shell's faces (a pcurve names
        # its basis surface) and the parameter space context of pcurves.
        self.surfaces, self.context2d = {}, None

    def add(self, text):
        self.lines.append(text)
        return len(self.lines)

    def point(self, p):
        return self.add("CARTESIAN_POINT('',(%s))" % ','.join(real(c) for c in p))

    def direction(self, d):
        return self.add("DIRECTION('',(%s))" % ','.join(real(c) for c in d))

    def placement(self, frame):
        o, axis, x = frame
        return self.add("AXIS2_PLACEMENT_3D('',#%d,#%d,#%d)"
                        % (self.point(o), self.direction(axis), self.direction(x)))

    def vertex(self, v):
        if id(v) not in self.vertices:
            self.vertices[id(v)] = self.add("VERTEX_POINT('',#%d)" % self.point(v.p))
        return self.vertices[id(v)]

    def edge(self, e):
        """An edge curve; with `flip_edges` every other edge is written with
        its vertices swapped and `same_sense` false (the use flips)."""
        if id(e) in self.edges:
            return self.edges[id(e)]
        start, end = self.vertex(e.start), self.vertex(e.end)
        if e.curve is not None:
            curve = self.curve(e.curve)
        elif e.circle is None:
            d = tuple(e.end.p[i]-e.start.p[i] for i in range(3))
            vector = self.add("VECTOR('',#%d,1.)" % self.direction(d))
            curve = self.add("LINE('',#%d,#%d)" % (self.point(e.start.p), vector))
        else:
            frame, r = e.circle
            curve = self.add("CIRCLE('',#%d,%s)" % (self.placement(frame), real(r)))
        if e.pcurves:
            pcurves = [self.pcurve(face, c2) for face, c2 in e.pcurves]
            curve = self.add("SURFACE_CURVE('',#%d,(%s),.CURVE_3D.)"
                             % (curve, ','.join('#%d' % p for p in pcurves)))
        flip = self.flip_edges and len(self.edges) % 2 == 1
        if flip:
            number = self.add("EDGE_CURVE('',#%d,#%d,#%d,.F.)" % (end, start, curve))
        else:
            number = self.add("EDGE_CURVE('',#%d,#%d,#%d,.T.)" % (start, end, curve))
        self.edges[id(e)] = (number, flip)
        return self.edges[id(e)]

    def curve(self, c):
        """An ellipse, or a B-spline curve of 2 or 3 dimensions: simple, or
        rational as the complex instance OCCT writes."""
        if c[0] == 'ELLIPSE':
            _, frame, a1, a2 = c
            return self.add("ELLIPSE('',#%d,%s,%s)" % (self.placement(frame), real(a1), real(a2)))
        if c[0] == 'LINE':
            _, p, d = c
            vector = self.add("VECTOR('',#%d,1.)" % self.direction(d))
            return self.add("LINE('',#%d,#%d)" % (self.point(p), vector))
        _, degree, poles, weights, knots, mults = c
        points = '(%s)' % ','.join('#%d' % self.point(p) for p in poles)
        mults = '(%s)' % ','.join(str(m) for m in mults)
        knots = '(%s)' % ','.join(real(k) for k in knots)
        if weights is None:
            return self.add("B_SPLINE_CURVE_WITH_KNOTS('',%d,%s,.UNSPECIFIED.,.F.,.F.,%s,%s,.UNSPECIFIED.)"
                            % (degree, points, mults, knots))
        return self.add("( BOUNDED_CURVE() B_SPLINE_CURVE(%d,%s,.UNSPECIFIED.,.F.,.F.) "
                        "B_SPLINE_CURVE_WITH_KNOTS(%s,%s,.UNSPECIFIED.) CURVE() "
                        "GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((%s)) "
                        "REPRESENTATION_ITEM('') )"
                        % (degree, points, mults, knots, ','.join(real(w) for w in weights)))

    def pcurve(self, face, c2):
        if self.context2d is None:
            self.context2d = self.add("( GEOMETRIC_REPRESENTATION_CONTEXT(2) "
                                      "PARAMETRIC_REPRESENTATION_CONTEXT() "
                                      "REPRESENTATION_CONTEXT('2D SPACE','') )")
        curve = self.curve(c2)
        rep = self.add("DEFINITIONAL_REPRESENTATION('',(#%d),#%d)" % (curve, self.context2d))
        return self.add("PCURVE('',#%d,#%d)" % (self.surfaces[id(face)], rep))

    def spline_surface(self, s):
        _, (du, dv), rows, weights, (uk, vk), (um, vm) = s
        grid = '(%s)' % ','.join('(%s)' % ','.join('#%d' % self.point(p) for p in row) for row in rows)
        ints = lambda xs: '(%s)' % ','.join(str(x) for x in xs)
        reals = lambda xs: '(%s)' % ','.join(real(x) for x in xs)
        if weights is None:
            return self.add("B_SPLINE_SURFACE_WITH_KNOTS('',%d,%d,%s,.UNSPECIFIED.,.F.,.F.,.F.,%s,%s,%s,%s,"
                            ".UNSPECIFIED.)" % (du, dv, grid, ints(um), ints(vm), reals(uk), reals(vk)))
        return self.add("( BOUNDED_SURFACE() B_SPLINE_SURFACE(%d,%d,%s,.UNSPECIFIED.,.F.,.F.,.F.) "
                        "B_SPLINE_SURFACE_WITH_KNOTS(%s,%s,%s,%s,.UNSPECIFIED.) "
                        "GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((%s)) "
                        "REPRESENTATION_ITEM('') SURFACE() )"
                        % (du, dv, grid, ints(um), ints(vm), reals(uk), reals(vk),
                           ','.join(reals(row) for row in weights)))

    def surface(self, s, flip):
        if s[0] == 'B_SPLINE_SURFACE':
            return self.spline_surface(s)
        kind, frame, *radii = s
        o, axis, x = frame
        if flip:
            # The same plane with the opposite normal.
            frame = (o, neg(axis), x)
        text = ','.join(['#%d' % self.placement(frame)]+[real(r) for r in radii])
        return self.add("%s('',%s)" % (kind, text))

    def face(self, f, index):
        flip = self.flip_faces and index % 2 == 1 and f.surface[0] == 'PLANE'
        bounds = []
        for k, loop in enumerate(f.bounds):
            reverse = self.flip_bounds and (index+k) % 2 == 0
            uses = [(e, not fwd) for e, fwd in reversed(loop)] if reverse else loop
            oriented = []
            for e, fwd in uses:
                number, swapped = self.edge(e)
                oriented.append(self.add("ORIENTED_EDGE('',*,*,#%d,%s)"
                                         % (number, '.T.' if fwd != swapped else '.F.')))
            lp = self.add("EDGE_LOOP('',(%s))" % ','.join('#%d' % o for o in oriented))
            kind = 'FACE_OUTER_BOUND' if k == 0 else 'FACE_BOUND'
            bounds.append(self.add("%s('',#%d,%s)" % (kind, lp, '.F.' if reverse else '.T.')))
        surface = self.surfaces[id(f)] if id(f) in self.surfaces else self.surface(f.surface, flip)
        sense = f.same_sense != flip
        return self.add("ADVANCED_FACE('',(%s),#%d,%s)"
                        % (','.join('#%d' % b for b in bounds), surface, '.T.' if sense else '.F.'))

    def shell(self, faces, kind='CLOSED_SHELL'):
        # With pcurves, every face's surface first (none flips: STEP-b's
        # cases have no flipped variants).
        if any(e.pcurves for f in faces for loop in f.bounds for e, _ in loop):
            for f in faces:
                self.surfaces[id(f)] = self.surface(f.surface, False)
        numbers = [self.face(f, i) for i, f in enumerate(faces)]
        return self.add("%s('',(%s))" % (kind, ','.join('#%d' % n for n in numbers)))

    def document(self, name, bodies, length='mm', angle='rad', uncertainty=1e-7, schema='ap214'):
        """`bodies`: [('solid', faces) | ('voids', outer, [void faces]) |
        ('sheet', [faces per shell])]."""
        file_schema, protocol = SCHEMAS[schema]
        app = self.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
        self.add("APPLICATION_PROTOCOL_DEFINITION('international standard','%s',2000,#%d)" % (protocol, app))
        pcontext = self.add("PRODUCT_CONTEXT('',#%d,'mechanical')" % app)
        product = self.add("PRODUCT('%s','%s','',(#%d))" % (name, name, pcontext))
        formation = self.add("PRODUCT_DEFINITION_FORMATION('','',#%d)" % product)
        dcontext = self.add("PRODUCT_DEFINITION_CONTEXT('part definition',#%d,'design')" % app)
        definition = self.add("PRODUCT_DEFINITION('design','',#%d,#%d)" % (formation, dcontext))
        shape = self.add("PRODUCT_DEFINITION_SHAPE('','',#%d)" % definition)
        # Units and the context.
        unit_lines = LENGTHS[length]
        if length == 'inch':
            base = self.add(unit_lines[0])
            exponents = self.add(unit_lines[1])
            measure = self.add(unit_lines[3].format(base=base))
            length_unit = self.add(unit_lines[2].format(m=measure, e=exponents))
        else:
            length_unit = self.add(unit_lines[0])
        if angle == 'deg':
            radian = self.add('( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )')
            measure = self.add('PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925199433),#%d)'
                               % radian)
            exponents = self.add('DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)')
            angle_unit = self.add("( CONVERSION_BASED_UNIT('DEGREE',#%d) NAMED_UNIT(#%d) PLANE_ANGLE_UNIT() )"
                                  % (measure, exponents))
        else:
            angle_unit = self.add('( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )')
        solid_angle = self.add('( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )')
        accuracy = self.add("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(%s),#%d,'distance_accuracy_value',"
                            "'confusion accuracy')" % (real(uncertainty), length_unit))
        context = self.add("( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#%d)) "
                           "GLOBAL_UNIT_ASSIGNED_CONTEXT((#%d,#%d,#%d)) REPRESENTATION_CONTEXT('Context #1',"
                           "'3D Context with UNIT and UNCERTAINTY') )"
                           % (accuracy, length_unit, angle_unit, solid_angle))
        origin = self.placement(XY)
        items = []
        sheet = False
        for body in bodies:
            if body[0] == 'solid':
                items.append(self.add("MANIFOLD_SOLID_BREP('',#%d)" % self.shell(body[1])))
            elif body[0] == 'voids':
                outer = self.shell(body[1])
                voids = []
                for faces in body[2]:
                    inner = self.shell(faces)
                    voids.append(self.add("ORIENTED_CLOSED_SHELL('',*,#%d,.F.)" % inner))
                items.append(self.add("BREP_WITH_VOIDS('',#%d,(%s))"
                                      % (outer, ','.join('#%d' % v for v in voids))))
            else:
                sheet = True
                shells = [self.shell(faces, 'OPEN_SHELL') for faces in body[1]]
                items.append(self.add("SHELL_BASED_SURFACE_MODEL('',(%s))" % ','.join('#%d' % s for s in shells)))
        kind = 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION' if sheet else 'ADVANCED_BREP_SHAPE_REPRESENTATION'
        representation = self.add("%s('',(%s),#%d)" % (kind, ','.join('#%d' % i for i in [origin]+items), context))
        self.add("SHAPE_DEFINITION_REPRESENTATION(#%d,#%d)" % (shape, representation))
        head = ['ISO-10303-21;', 'HEADER;', "FILE_DESCRIPTION(('rustyOCCT STEP fixture'),'2;1');",
                "FILE_NAME('%s.stp','2026-09-28T00:00:00',('rustyOCCT'),('rustyOCCT'),"
                "'generate_step_fixtures.py','generate_step_fixtures.py','');" % name,
                'FILE_SCHEMA((%s));' % file_schema, 'ENDSEC;', 'DATA;']
        body = ['#%d=%s;' % (i+1, line) for i, line in enumerate(self.lines)]
        return '\n'.join(head+body+['ENDSEC;', 'END-ISO-10303-21;'])+'\n'


# The hand-written file: a corner tetrahedron of edge 10 with the syntax a
# conforming reader must take (see the module docstring).
SYNTAX = r"""ISO-10303-21;
HEADER;
/* hand-written for the STEP import track: a corner tetrahedron */
FILE_DESCRIPTION(('corner tetrahedron, it''s hand-written'),'2;1');
FILE_NAME('syntax.stp','2026-09-28T00:00:00',('rustyOCCT'),('\X2\00E9\X0\quipe'),
  'hand','hand','');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#900 = SHAPE_DEFINITION_REPRESENTATION ( #901 , #950 ) ;
#901=PRODUCT_DEFINITION_SHAPE('','',#902);
#902=PRODUCT_DEFINITION('design','',#903,#906);
#903=PRODUCT_DEFINITION_FORMATION('','',#904);
#904=PRODUCT('tet','tet','',(#905));
#905=PRODUCT_CONTEXT('',#907,'mechanical');
#906=PRODUCT_DEFINITION_CONTEXT('part definition',#907,'design');
#907=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#950=ADVANCED_BREP_SHAPE_REPRESENTATION('tet',(#1,#960),#970);
#960=AXIS2_PLACEMENT_3D('',#961,$,$);
#961=CARTESIAN_POINT('',(0.,0.,0.));
#970=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#973))
  GLOBAL_UNIT_ASSIGNED_CONTEXT((#971,#972,#974))REPRESENTATION_CONTEXT('',''));
#971=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#972=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));
#974=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());
#973=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#971,'distance_accuracy_value','');
#1=MANIFOLD_SOLID_BREP('',#2);
#2=CLOSED_SHELL('',(#10,#20,#30,#40));
/* the face on the slanted plane x + y + z = 10, normal unnormalised */
#40=ADVANCED_FACE('slanted',(#41),#45,.T.);
#41=FACE_OUTER_BOUND('',#42,.T.);
#42=EDGE_LOOP('',(#43,#44,#46));
#43=ORIENTED_EDGE('',*,*,#112,.T.);
#44=ORIENTED_EDGE('',*,*,#123,.T.);
#46=ORIENTED_EDGE('',*,*,#131,.T.);
#45=PLANE('',#47);
#47=AXIS2_PLACEMENT_3D('',#48,#49,#50);
#48=CARTESIAN_POINT('',(10.,0.,0.));
#49=DIRECTION('',(1.,1.,1.));
#50=DIRECTION('',(-1.,1.,0.));
/* x = 0, normal -x */
#10=ADVANCED_FACE('',(#11),#15,.T.);
#11=FACE_OUTER_BOUND('',#12,.T.);
#12=EDGE_LOOP('',(#13,#14,#16));
#13=ORIENTED_EDGE('',*,*,#103,.T.);
#14=ORIENTED_EDGE('',*,*,#123,.F.);
#16=ORIENTED_EDGE('',*,*,#102,.F.);
#15=PLANE('',#17);
#17=AXIS2_PLACEMENT_3D('',#18,#19,#53);
#53=DIRECTION('',(0.,0.,1.));
#18=CARTESIAN_POINT('',(0.,0.,0.));
#19=DIRECTION('',(-1.,0.,0.));
/* y = 0, normal -y, its surface's normal +y */
#20=ADVANCED_FACE('',(#21),#25,.F.);
#21=FACE_BOUND('',#22,.F.);
#22=EDGE_LOOP('',(#23,#24,#26));
#23=ORIENTED_EDGE('',*,*,#103,.T.);
#24=ORIENTED_EDGE('',*,*,#131,.T.);
#26=ORIENTED_EDGE('',*,*,#101,.F.);
#25=PLANE('',#27);
#27=AXIS2_PLACEMENT_3D('',#28,#29,#51);
#28=CARTESIAN_POINT('',(0.,0.,0.));
#29=DIRECTION('',(0.,2.,0.));
#51=DIRECTION('',(0.,0.,1.));
/* z = 0, normal -z */
#30=ADVANCED_FACE('',(#31),#35,.T.);
#31=FACE_OUTER_BOUND('',#32,.T.);
#32=EDGE_LOOP('',(#33,#34,#36));
#33=ORIENTED_EDGE('',*,*,#102,.T.);
#34=ORIENTED_EDGE('',*,*,#112,.F.);
#36=ORIENTED_EDGE('',*,*,#101,.F.);
#35=PLANE('',#37);
/* no reference direction: ISO 10303-42's default x */
#37=AXIS2_PLACEMENT_3D('',#38,#39,$);
#38=CARTESIAN_POINT('',(0.,0.,0.));
#39=DIRECTION('',(0.,0.,-1.));
/* edges: O-X, O-Y, O-Z, X-Y, Y-Z, Z-X; the last with same_sense false */
#101=EDGE_CURVE('',#201,#202,#301,.T.);
#102=EDGE_CURVE('',#201,#203,#302,.T.);
#103=EDGE_CURVE('',#201,#204,#303,.T.);
#112=EDGE_CURVE('',#202,#203,#312,.T.);
#123=EDGE_CURVE('',#203,#204,#323,.T.);
#131=EDGE_CURVE('',#204,#202,#331,.F.);
#201=VERTEX_POINT('',#211);
#202=VERTEX_POINT('',#212);
#203=VERTEX_POINT('',#213);
#204=VERTEX_POINT('',#214);
#211=CARTESIAN_POINT('',(0.,0.,0.));
#212=CARTESIAN_POINT('',(10.,0.,0.));
#213=CARTESIAN_POINT('',(0.,10.,0.));
#214=CARTESIAN_POINT('',(0.,0.,10.));
#301=LINE('',#211,#401);
#302=LINE('',#211,#402);
#303=LINE('',#211,#403);
#312=LINE('',#212,#412);
#323=LINE('',#213,#423);
#331=LINE('',#212,#431);
#401=VECTOR('',#501,10.);
#402=VECTOR('',#502,1.);
#403=VECTOR('',#503,1.);
#412=VECTOR('',#512,1.);
#423=VECTOR('',#523,1.);
#431=VECTOR('',#531,1.);
#501=DIRECTION('',(1.,0.,0.));
#502=DIRECTION('',(0.,1.,0.));
#503=DIRECTION('',(0.,0.,1.));
#512=DIRECTION('',(-1.,1.,0.));
#523=DIRECTION('',(0.,-1.,1.));
#531=DIRECTION('',(-1.,0.,1.));
ENDSEC;
END-ISO-10303-21;
"""

SEMI = mpmath.atan(mpmath.mpf(1)/2)  # the frusta's and the apex cone's half-angle


def cases():
    """[(name, text, [(closed form, ...)])]: one closed form per body, in
    the order the file's bodies are numbered."""
    out = []

    def case(name, bodies, forms, **kw):
        flips = {k: kw.pop(k) for k in ('flip_faces', 'flip_edges', 'flip_bounds') if k in kw}
        out.append((name, Writer(**flips).document(name, bodies, **kw), forms))

    lo, hi = (0.0, 0.0, 0.0), (10.0, 20.0, 30.0)
    box = ref.combine((1, ref.box(lo, hi)))
    case('box', [('solid', box_faces(lo, hi))], [box])
    case('box_ap203', [('solid', box_faces(lo, hi))], [box], schema='ap203')
    case('box_ap242', [('solid', box_faces(lo, hi))], [box], schema='ap242')
    case('box_flipped', [('solid', box_faces(lo, hi))], [box], flip_faces=True, flip_edges=True,
         flip_bounds=True)
    # Lengths in metres and inches: the same box scaled to millimetres by
    # one binary64 multiplication per coordinate.
    metres = [(0.0, 0.0, 0.0), (0.01, 0.02, 0.03)]
    case('box_metre', [('solid', box_faces(*metres))],
         [ref.combine((1, ref.box(*[[c*1000.0 for c in p] for p in metres])))],
         length='m', uncertainty=1e-10)
    inches = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
    # The file's numbers are binary64, and so is the factor (25.4 as read).
    case('box_inch', [('solid', box_faces(*inches))],
         [ref.combine((1, ref.box(*[[c*25.4 for c in p] for p in inches])))],
         length='inch', uncertainty=1e-7/25.4)
    ell = [(0.0, 0.0), (12.0, 0.0), (12.0, 4.0), (4.0, 4.0), (4.0, 9.0), (0.0, 9.0)]
    case('l_prism', [('solid', prism(ell, -2.0, 3.0))], [ref.combine((1, ref.polygon_prism(ell, -2.0, 3.0)))])
    rect = [(0.0, 0.0), (40.0, 0.0), (40.0, 30.0), (0.0, 30.0)]
    case('plate_hole', [('solid', prism(rect, 0.0, 5.0, holes=[(15.0, 12.0, 4.0)]))],
         [ref.combine((1, ref.polygon_prism(rect, 0.0, 5.0)),
                      (-1, ref.hole((15.0, 12.0, 0.0), (0, 0, 1), 4.0, 5.0)))])
    case('box_void', [('voids', box_faces((0.0, 0.0, 0.0), (20.0, 20.0, 20.0)),
                       [box_faces((5.0, 6.0, 7.0), (10.0, 12.0, 14.0))])],
         [ref.combine((1, ref.box((0, 0, 0), (20, 20, 20))), (-1, ref.box((5, 6, 7), (10, 12, 14))))])
    z, x = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
    case('cylinder', [('solid', cylinder_faces((0.0, 0.0, 0.0), z, x, 5.0, 12.0))],
         [ref.combine((1, ref.cylinder((0, 0, 0), (0, 0, 1), 5.0, 12.0)))])
    # A tilted, displaced axis: a Pythagorean direction, so every written
    # coordinate is exact.
    tilt, tx = (0.0, 0.6, 0.8), (1.0, 0.0, 0.0)
    case('cylinder_tilted', [('solid', cylinder_faces((1.0, 2.0, 3.0), tilt, tx, 4.0, 10.0))],
         [ref.combine((1, ref.cylinder((1, 2, 3), (0, Fraction('0.6'), Fraction('0.8')), 4.0, 10.0)))])
    case('half_cylinder', [('solid', half_cylinder_faces((0.0, 0.0, 0.0), 6.0, 8.0))],
         [ref.combine((1, ref.half_cylinder((0, 0, 0), 6.0, 8.0)))])
    case('cone', [('solid', apex_cone_faces((0.0, 0.0, 0.0), 10.0, 20.0, ref.rn(SEMI)))],
         [ref.combine((1, ref.apex_cone_down((0, 0, 0), 10.0, 20.0)))])
    case('frustum', [('solid', frustum_faces((0.0, 0.0, 0.0), 5.0, 10.0, 10.0, ref.rn(SEMI)))],
         [ref.combine((1, ref.frustum((0, 0, 0), 5.0, 10.0, 10.0)))])
    case('frustum_degree', [('solid', frustum_faces((0.0, 0.0, 0.0), 5.0, 10.0, 10.0,
                                                    ref.rn(SEMI*180/mpmath.pi)))],
         [ref.combine((1, ref.frustum((0, 0, 0), 5.0, 10.0, 10.0)))], angle='deg')
    case('sphere', [('solid', sphere_faces((1.0, -2.0, 3.0), 7.0))],
         [ref.combine((1, ref.sphere((1, -2, 3), 7.0)))])
    case('hemisphere', [('solid', hemisphere_faces((0.0, 0.0, 0.0), 5.0))],
         [ref.combine((1, ref.hemisphere((0, 0, 0), 5.0)))])
    case('torus', [('solid', torus_faces((0.0, 0.0, 0.0), 10.0, 3.0))],
         [ref.combine((1, ref.torus((0, 0, 0), 10.0, 3.0)))])
    case('elbow', [('solid', elbow_faces((0.0, 0.0, 0.0), 10.0, 3.0))],
         [ref.combine((1, ref.elbow((0, 0, 0), 10.0, 3.0)))])
    case('two_solids', [('solid', box_faces((0.0, 0.0, 0.0), (4.0, 4.0, 4.0))),
                        ('solid', cylinder_faces((10.0, 0.0, 0.0), z, x, 2.0, 4.0))],
         [ref.combine((1, ref.box((0, 0, 0), (4, 4, 4)))),
          ref.combine((1, ref.cylinder((10, 0, 0), (0, 0, 1), 2.0, 4.0)))])
    open_box = box_faces((0.0, 0.0, 0.0), (10.0, 10.0, 5.0))
    del open_box[1]  # the top
    case('open_box', [('sheet', [open_box])], [ref.open_box((0, 0, 0), (10, 10, 5))])
    out.append(('syntax', SYNTAX, [ref.combine((1, ref.tetrahedron((0, 0, 0), 10.0)))]))
    # STEP-b.
    c, x = (1, 2, 3), (Fraction('0.6'), Fraction('0.8'), 0)
    case('ellipse_sheet', [('sheet', [half_ellipse_faces(c, x, 5, 3)])],
         [ref.half_ellipse(c, (-x[1], x[0], 0), 5, 3)])
    normal, major = (0, Fraction('-0.6'), Fraction('0.8')), (0, Fraction('0.8'), Fraction('0.6'))
    case('cylinder_oblique', [('solid', oblique_cylinder_faces(5, 10, normal, major))],
         [ref.combine((1, ref.oblique_cylinder(5, 10, Fraction(3, 4))))])
    prof = ref.profile(*[PROFILE[i] for i in (0, 1, 2, 3)])
    case('bspline_plate', [('sheet', [bspline_plate_faces(0)])], [ref.bspline_plate(prof, 0)])
    case('bspline_prism', [('solid', bspline_prism_faces(4))], [ref.combine((1, ref.bspline_prism(prof, 0, 4)))])
    case('bspline_patch', [('sheet', [patch_faces()])], [ref.spline_sheet(patch_pieces())])
    case('bspline_trimmed', [('sheet', [trimmed_faces()])], [ref.spline_sheet([(TRIM_GRID, ref.power(TRIM))])])
    case('rational_cylinder', [('solid', rational_cylinder_faces(5, 8))],
         [ref.combine((1, ref.cylinder((0, 0, 0), (0, 0, 1), 5, 8)))])
    return out


# The cases of STEP-b, captured natively apart from STEP-a's.
STEP_B = ('ellipse_sheet', 'cylinder_oblique', 'bspline_plate', 'bspline_prism', 'bspline_patch',
          'bspline_trimmed', 'rational_cylinder')


# Per case, the largest geometry gap relative to its size and the spline
# uses checked through a pcurve (`rows`).
GAPS = {}


def rows(made=None):
    out = ['case\tentity\tclass\tV E W F SH SO\tvolume\tarea\tcx cy cz']
    for name, text, forms in made or cases():
        _, data = ref.parse(text)
        found = ref.bodies(data)
        if len(found) != len(forms):
            raise ValueError(f'{name}: {len(found)} bodies, {len(forms)} closed forms')
        # Every edge on its vertices and its faces (in the file's units).
        gap, spline_uses = ref.geometry_gaps(data, found[0][3][1])
        size = max([1.0]+[abs(c) for n in data if 'CARTESIAN_POINT' in ref.kinds(data, n)
                          for c in ref.point(data, n)])
        if gap > 1e-12*size:
            raise ValueError(f'{name}: geometry gap {mpmath.nstr(gap, 3)}')
        GAPS[name] = (gap/size, spline_uses)
        for (entity, kind, counts, _), form in zip(found, forms):
            if kind == 'sheet':
                area, centre = form
                volume = '-'
            else:
                v, area, centre = form
                volume = repr(ref.rn(v))
            out.append('\t'.join([name, str(entity), kind, ' '.join(map(str, counts)), volume,
                                  repr(ref.rn(area)), ' '.join(repr(ref.rn(c)) for c in centre)]))
    return '\n'.join(out)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    made = cases()
    files = {name+'.stp': text for name, text, _ in made}
    table = rows(made)
    worst = max(GAPS, key=lambda n: GAPS[n][0])
    print(f'largest geometry gap {mpmath.nstr(GAPS[worst][0], 3)} of the size ({worst}); '
          f'{sum(u for _, u in GAPS.values())} spline uses checked through their pcurves')
    if args.check:
        stale = [n for n, t in files.items() if not (OUT/n).exists() or (OUT/n).read_text() != t]
        stale += sorted(p.name for p in OUT.glob('*.stp') if p.name not in files)
        if EXPECTED.read_text() != table:
            stale.append(EXPECTED.name)
        if stale:
            raise SystemExit('stale STEP fixtures: '+', '.join(stale))
        print(f'{len(files)} STEP fixtures up to date')
        return
    OUT.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (OUT/name).write_text(text)
    EXPECTED.write_text(table)
    print(f'wrote {len(files)} STEP fixtures')


if __name__ == '__main__':
    main()
