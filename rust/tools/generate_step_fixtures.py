#!/usr/bin/env python3
"""Fixtures for STEP-a of REVIEW_NOTES.md: small STEP files authored here.

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

`rust/fixtures/step/NAME.stp` holds each file and `step-expected.tsv` each
body: its entity number, class, OCCT's counts (`step_reference.bodies`,
from the independent parser) and its volume, area and centre from the
closed forms of `step_reference.py` (a sheet has no volume: `-`). `--check`
regenerates everything and fails on any difference. No Rust or native
result supplies an expectation.
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
    `start` to `end` (the same vertex for a closed circle)."""

    def __init__(self, start, end, circle=None):
        self.start, self.end, self.circle = start, end, circle


class Face:
    """`surface`: ('PLANE', frame) | ('CYLINDRICAL_SURFACE', frame, r) |
    ('CONICAL_SURFACE', frame, r, semi_angle) | ('SPHERICAL_SURFACE', frame,
    r) | ('TOROIDAL_SURFACE', frame, R, r), a frame being (origin, axis, x);
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
        if e.circle is None:
            d = tuple(e.end.p[i]-e.start.p[i] for i in range(3))
            vector = self.add("VECTOR('',#%d,1.)" % self.direction(d))
            curve = self.add("LINE('',#%d,#%d)" % (self.point(e.start.p), vector))
        else:
            frame, r = e.circle
            curve = self.add("CIRCLE('',#%d,%s)" % (self.placement(frame), real(r)))
        flip = self.flip_edges and len(self.edges) % 2 == 1
        if flip:
            number = self.add("EDGE_CURVE('',#%d,#%d,#%d,.F.)" % (end, start, curve))
        else:
            number = self.add("EDGE_CURVE('',#%d,#%d,#%d,.T.)" % (start, end, curve))
        self.edges[id(e)] = (number, flip)
        return self.edges[id(e)]

    def surface(self, s, flip):
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
        surface = self.surface(f.surface, flip)
        sense = f.same_sense != flip
        return self.add("ADVANCED_FACE('',(%s),#%d,%s)"
                        % (','.join('#%d' % b for b in bounds), surface, '.T.' if sense else '.F.'))

    def shell(self, faces, kind='CLOSED_SHELL'):
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
    return out


def rows():
    out = ['case\tentity\tclass\tV E W F SH SO\tvolume\tarea\tcx cy cz']
    for name, text, forms in cases():
        _, data = ref.parse(text)
        found = ref.bodies(data)
        if len(found) != len(forms):
            raise ValueError(f'{name}: {len(found)} bodies, {len(forms)} closed forms')
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
    files = {name+'.stp': text for name, text, _ in cases()}
    table = rows()
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
