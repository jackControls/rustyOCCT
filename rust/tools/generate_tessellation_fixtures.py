#!/usr/bin/env python3
"""Fixtures for tessellation (T-a and T-b of REVIEW_NOTES.md).

`tessellation-cases.txt` holds identity-protocol case blocks
(identity_reference.encode_case: prisms of profiles with lines, arcs,
circles and holes, a box, cones, spheres, tori and face bodies), each with
`mesh SETTING DEFLECTION ANGLE` rows before its `end`: a coarse setting (a
hundredth of the case's scale, 0.5 rad, BRepMesh's default angle) and a fine
one (a thousandth, 0.3 rad). `tessellation-expected.tsv` gives what
tessellation_reference.py derives from each case alone: whether it is a
solid, the Euler characteristic of its boundary, a face's number of boundary
loops, exact area and volume, and the scale the settings are relative to.
No Rust or OCCT result supplies an expectation.

T-b: `tessellation-spline-cases.txt` holds twelve spline bodies in the B-rep
line protocol (cell_reference.encode of generate_brep_fixtures' builder,
each checked valid by the reference validator), with the same settings, and
`tessellation-spline-expected.tsv` what tessellation_reference.py derives
from each `SplineCase` alone.
"""
import argparse
import copy
import dataclasses
from fractions import Fraction as F
import math
from pathlib import Path

import mpmath as mp

from brep_reference import Arc3, Cylinder, Edge, Face, Frame, Line2, Line3, Model, Plane, Use, cos_rn, sin_rn
from cell_reference import declare, encode as encode_cell, to_cell, validate as validate_cell
import generate_brep_fixtures as gbf
from identity_reference import Boundary, Case, encode_case
import spline_cell_reference as spl
from spline_cell_reference import Basis, BSpline2, BSpline3, BSplineSurface, homogeneous as spl_homogeneous
import tessellation_reference as ref
from tessellation_reference import SplineCase

ROOT = Path(__file__).resolve().parents[1]
Z = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (1.0, -2.0, 0.5, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)
FAR = (10000.0, -20000.0, 5000.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
H = ref.HALF_PI
T = ref.TWO_PI


def rect(w, d, x=0.0, y=0.0):
    return Boundary(points=[(x, y), (x+w, y), (x+w, y+d), (x, y+d)])


def path(points_segments):
    """[(point, segment)]: segment None for a line, (cx, cy, r, ccw) for an arc."""
    return Boundary(points=[p for p, _ in points_segments], segments=[s for _, s in points_segments])


def stadium(x0, x1, r):
    return path([((x0, -r), None), ((x1, -r), (x1, 0.0, r, True)),
                 ((x1, r), None), ((x0, r), (x0, 0.0, r, True))])


def prism(name, boundaries, start, end, frame=Z):
    return Case(name, 1e-7, 1, frame, start, end, boundaries)


def cases():
    rounded = path([((1.0, 0.0), None), ((7.0, 0.0), (7.0, 1.0, 1.0, True)), ((8.0, 1.0), None),
                    ((8.0, 5.0), (7.0, 5.0, 1.0, True)), ((7.0, 6.0), None),
                    ((1.0, 6.0), (1.0, 5.0, 1.0, True)), ((0.0, 5.0), None),
                    ((0.0, 1.0), (1.0, 1.0, 1.0, True))])
    notch = path([((0.0, 0.0), None), ((10.0, 0.0), None), ((10.0, 6.0), None),
                  ((6.0, 6.0), (5.0, 6.0, 1.0, False)), ((4.0, 6.0), None), ((0.0, 6.0), None)])
    lens = path([((0.0, -3.0), (-4.0, 0.0, 5.0, True)), ((0.0, 3.0), (4.0, 0.0, 5.0, True))])
    scallop = path([((0.0, 0.0), (2.0, 1.5, 2.5, True)), ((4.0, 0.0), (2.5, 2.0, 2.5, True)),
                    ((4.0, 4.0), (2.0, 2.5, 2.5, True)), ((0.0, 4.0), (1.5, 2.0, 2.5, True))])
    square = rect(20.0, 20.0, -10.0, -10.0)
    mixed = [square,
             path([((-6.0, -6.0), None), ((-2.0, -6.0), (-2.0, -4.0, 2.0, True)), ((-2.0, -2.0), None)]),
             Boundary(circle=(5.0, 5.0, 1.5)),
             path([((4.0, -6.0), (5.0, -6.0, 1.0, True)), ((6.0, -6.0), (5.0, -6.0, 1.0, True))])]
    ell = Boundary(points=[(0.0, 0.0), (12.0, 0.0), (12.0, 4.0), (4.0, 4.0), (4.0, 10.0), (0.0, 10.0)])
    out = [
        Case('box', 1e-7, 1, Z, 0.0, 0.0, [], box=((1.0, 2.0, 3.0), (40.0, 20.0, 10.0))),
        prism('plate_hole', [rect(40.0, 20.0), Boundary(circle=(10.0, 10.0, 3.0))], 0.0, 5.0),
        prism('cylinder', [Boundary(circle=(0.0, 0.0, 5.0))], 0.0, 12.0),
        prism('tall_cylinder', [Boundary(circle=(0.0, 0.0, 1.0))], 0.0, 40.0),
        prism('thin_disc', [Boundary(circle=(0.0, 0.0, 10.0))], 0.0, 0.5),
        prism('stadium_slot', [stadium(0.0, 12.0, 4.0), rect(2.0, 2.0, 5.0, -1.0)], 0.0, 2.0),
        prism('rounded_rectangle', [rounded, Boundary(circle=(4.0, 3.0, 1.0))], 0.0, 3.0),
        prism('notch', [notch], 0.0, 4.0),
        prism('lens', [lens], 0.0, 2.0),
        prism('scallop', [scallop], 0.0, 1.0),
        prism('holes_mixed', mixed, 0.0, 3.0),
        prism('tilted_prism', [ell, Boundary(circle=(2.0, 7.0, 1.0))], -1.0, 4.0, TILT),
        prism('far_prism', [Boundary(circle=(0.0, 0.0, 2.0))], 0.0, 3.0, FAR),
        Case('cone_apex', 1e-7, 1, Z, 0.0, 0.0, [], cone=(2.0, 0.0, 3.0)),
        Case('frustum', 1e-7, 1, Z, 0.0, 0.0, [], cone=(2.0, 1.0, 3.0)),
        Case('cone_inverted', 1e-7, 1, Z, 0.0, 0.0, [], cone=(0.0, 1.5, 2.0)),
        Case('cone_tilted', 1e-7, 1, TILT, 0.0, 0.0, [], cone=(1.0, 2.5, 0.5)),
        Case('sphere', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, -H, H)),
        Case('hemisphere', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, 0.0, H)),
        Case('sphere_zone', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, -0.5, 0.8)),
        Case('sphere_far', 1e-7, 1, FAR, 0.0, 0.0, [], sphere=(2.0, -H, H)),
        Case('torus', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, 0.0, T, T)),
        Case('torus_segment', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, -1.0, 1.0, T)),
        Case('torus_inner_half', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, H, 4.71238898038469, T)),
        Case('torus_wedge', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, 0.0, T, 2.0)),
        Case('torus_tilted', 1e-7, 1, TILT, 0.0, 0.0, [], torus=(3.0, 1.0, 0.0, T, T)),
        Case('face_holes', 1e-7, 1, Z, 0.0, 0.0,
             [rect(20.0, 20.0, -10.0, -10.0), rect(4.0, 4.0, -7.0, -7.0),
              Boundary(circle=(5.0, 5.0, 1.5)), stadium(0.0, 3.0, 1.0)], make='face'),
        Case('face_lens_tilted', 1e-7, 1, TILT, 0.0, 0.0, [lens], make='face'),
    ]
    return out


def settings(case):
    """(label, deflection, angle) rows of a case, relative to its scale."""
    scale = ref.body(case).scale
    return [('coarse', scale/100, 0.5), ('fine', scale/1000, 0.3)]


def number(x):
    return repr(float(x))


def encode(case):
    block = encode_case(case)
    rows = [f'mesh {label} {number(d)} {number(a)}' for label, d, a in settings(case)]
    head, end = block.rsplit('\nend', 1)
    return head+'\n'+'\n'.join(rows)+'\nend'+end


def generate():
    blocks = []
    rows = ['# case\tkind\teuler\tloops\tarea\tvolume\tscale (tessellation_reference.py)']
    for c in cases():
        blocks.append(encode(c))
        b = ref.body(c)
        rows.append('\t'.join([c.name, 'solid' if b.solid else 'face', str(b.euler), str(b.loops),
                               mp.nstr(b.area, 17), mp.nstr(b.volume, 17), number(b.scale)]))
    return {'tessellation-cases.txt': '\n'.join(blocks)+'\n',
            'tessellation-expected.tsv': '\n'.join(rows)+'\n', **spline_generate()}


# ------------------------------------------------------------------ splines (T-b)

QUADRATIC = Basis(2, [0.0, 1.0], [3, 3])
CUBIC_KNOT = Basis(3, [0.0, 0.5, 1.0], [4, 1, 4])
WAVE = Basis(3, [0.0, 0.25, 0.5, 0.75, 1.0], [4, 1, 1, 1, 4])
# A translation far from the origin (a rotation's rounded frames leave the
# reference validator unable to bound a spline use on a moved plane).
MOTION = ((0.0, 0.0, 1.0), 0.0, (1000.0, -2000.0, 500.0))


def bulge(spline):
    return [('line', (0.0, 0.0)), spline, ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]


def ring_curve():
    """A periodic cubic of six poles, uniform: a smooth closed blob."""
    basis = Basis(3, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [1]*7, periodic=True)
    poles = [(3.0, 0.0), (2.0, 2.5), (-1.0, 2.75), (-3.0, 0.5), (-2.0, -2.5), (1.0, -2.75)]
    return BSpline2(basis, poles, [1.0]*6)


def dome_surface():
    """S(u, v) = (u, v, z) over [0, 3] x [0, 2]: the Greville abscissae as
    poles reproduce u and v exactly; z is 1 on the boundary rows."""
    gu, gv = [0.0, 0.5, 1.5, 2.5, 3.0], [0.0, 1.0, 2.0]
    lift = {(1, 1): 1.5, (2, 1): 1.75, (3, 1): 1.375}
    poles = [(x, y, lift.get((i, j), 1.0)) for i, x in enumerate(gu) for j, y in enumerate(gv)]
    return BSplineSurface(Basis(3, [0.0, 1.5, 3.0], [4, 1, 4]), Basis(2, [0.0, 2.0], [3, 3]), poles,
                          [1.0]*len(poles))


def wavy_surface():
    """A bicubic of 5 x 4 poles, a knot at u = 1: a warped hill."""
    xs, ys = [0.0, 1.0, 2.5, 4.0, 5.0], [0.0, 1.25, 2.5, 3.75]
    z = [[0.0, 0.25, 0.25, 0.0], [0.5, 1.25, 1.0, 0.25], [0.25, 1.5, 1.625, 0.5], [0.5, 0.875, 1.25, 0.25],
         [0.0, 0.25, 0.375, 0.125]]
    poles = [(x+0.25*(j == 1)-0.125*(j == 2), y+0.125*(i == 2), z[i][j])
             for i, x in enumerate(xs) for j, y in enumerate(ys)]
    return BSplineSurface(Basis(3, [0.0, 1.0, 2.0], [4, 1, 4]), Basis(3, [0.0, 1.0], [4, 4]), poles,
                          [1.0]*len(poles))


def rational_surface():
    """A rational biquadratic patch with weights far from equal."""
    poles = [(0.0, 0.0, 0.0), (0.0, 1.5, 0.75), (0.0, 3.0, 0.0),
             (1.5, 0.0, 1.0), (1.5, 1.5, 2.0), (1.5, 3.0, 1.0),
             (3.0, 0.0, 0.0), (3.0, 1.5, 0.625), (3.0, 3.0, 0.0)]
    weights = [1.0, 0.625, 1.0, 2.5, 1.5, 2.5, 1.0, 0.625, 1.0]
    return BSplineSurface(QUADRATIC, QUADRATIC, poles, weights)


def spline_cases():
    """The T-b bodies: arguments only (tessellation_reference.SplineCase)."""
    cubic = ('spline', (3.0, 0.0), [(3.5, 0.25), (4.0, 1.0), (3.5, 1.75)], CUBIC_KNOT, None)
    wave = [('line', (0.0, 0.0)), ('line', (8.0, 0.0)),
            ('spline', (8.0, 4.0), [(7.0, 5.0), (5.5, 3.25), (4.0, 4.75), (2.5, 3.25), (1.0, 5.0)], WAVE, None),
            ('line', (0.0, 4.0))]
    circle_hole = [('circle', (2.0, 2.0), 0.625, False)]
    lens_hole = [('spline', (5.0, 2.0), [(5.5, 2.75)], QUADRATIC, None),
                 ('spline', (6.0, 2.0), [(5.5, 1.25)], QUADRATIC, None)]
    corner = [('line', (0.0, 0.0)), ('line', (3.0, 0.0)),
              ('spline', (3.0, 1.5), [(3.0, 2.0)], QUADRATIC, [1.0, 0.7071067811865476, 1.0]),
              ('line', (2.5, 2.0)), ('line', (0.0, 2.0))]
    sharp = [('line', (0.0, 0.0)), ('spline', (4.0, 0.0), [(4.0, 4.0)], QUADRATIC, [1.0, 4.0, 1.0]),
             ('line', (0.0, 4.0))]
    return [
        SplineCase('spline_bulge', 'prism', [bulge(('spline', (3.0, 0.0), [(4.0, 1.0)], QUADRATIC, None))]),
        SplineCase('spline_cubic_bulge', 'prism', [bulge(cubic)]),
        SplineCase('spline_rounded_corner', 'prism', [corner], -0.5, 0.5),
        SplineCase('spline_wave_holes', 'prism', [wave, circle_hole, lens_hole], 0.0, 1.5),
        SplineCase('spline_sharp_rational', 'prism', [sharp], 0.0, 2.0),
        SplineCase('spline_bulge_far', 'prism', [bulge(cubic)], motion=MOTION),
        SplineCase('spline_stadium_edge', 'prism', [gbf.stadium_boundary(3.0, 1.0)], variant='stadium_edge'),
        SplineCase('spline_stadium_pcurve', 'prism', [gbf.stadium_boundary(3.0, 1.0)], variant='stadium_pcurve'),
        SplineCase('face_spline_ring', 'face', ring=ring_curve()),
        SplineCase('spline_dome', 'dome', surface=dome_surface(), size=(3.0, 2.0, 1.0)),
        SplineCase('sheet_spline_hole', 'sheet', surface=wavy_surface(), hole=(0.625, 1.375, 0.25, 0.75)),
        SplineCase('sheet_rational', 'sheet', surface=rational_surface()),
    ]


def moved(m, motion):
    """A rigid motion (unit axis, angle, shift) of every 3D datum, splines
    included (generate_brep_fixtures.transformed has lines and arcs only)."""
    axis, angle, shift = motion
    k = [a/math.sqrt(sum(b*b for b in axis)) for a in axis]
    c, s = cos_rn(angle), sin_rn(angle)

    def rot(v):
        kv = sum(a*b for a, b in zip(k, v))
        kx = (k[1]*v[2]-k[2]*v[1], k[2]*v[0]-k[0]*v[2], k[0]*v[1]-k[1]*v[0])
        return tuple(v[i]*c+kx[i]*s+k[i]*kv*(1-c) for i in range(3))
    point = lambda p: tuple(a+b for a, b in zip(rot(p), shift))
    frame = lambda f: Frame(point(f.origin), rot(f.normal), rot(f.x))
    out = copy.deepcopy(m)
    out.vertices = [point(v) for v in m.vertices]
    for e in out.edges:
        cv = e.curve
        if isinstance(cv, Line3):
            e.curve = Line3(point(cv.start), point(cv.end))
        elif isinstance(cv, BSpline3):
            e.curve = dataclasses.replace(cv, poles=[point(p) for p in cv.poles])
        else:
            e.curve = Arc3(frame(cv.frame), cv.radius, cv.start, cv.sweep)
    for f in out.faces:
        sf = f.surface
        if isinstance(sf, BSplineSurface):
            f.surface = dataclasses.replace(sf, poles=[point(p) for p in sf.poles])
        elif isinstance(sf, Plane):
            f.surface = Plane(frame(sf.frame))
        else:
            f.surface = Cylinder(frame(sf.frame), sf.radius)
    return out


def _value(basis, hctrl, t):
    """The homogeneous value of a spline at t, exactly (de Boor)."""
    t = F(t)
    flat, n, (_, e) = basis.flat()
    p = basis.degree
    k = max(j for j in range(p, len(flat)-p-1) if flat[j] <= t and (t < flat[j+1] or flat[j+1] == e))
    d = {j: hctrl[j] for j in range(k-p, k+1)}
    for r in range(1, p+1):
        for j in range(k, k-p+r-1, -1):
            al = (t-flat[j])/(flat[j+p-r+1]-flat[j])
            d[j] = tuple((1-al)*x+al*y for x, y in zip(d[j-1], d[j]))
    return d[k]


def iso_curve(s, axis, c, rng=None):
    """The iso-curve u = c (axis 0, a curve in v) or v = c (axis 1, in u) of
    a nonperiodic spline surface, exactly, its poles rounded to binary64."""
    h = spl_homogeneous(s.poles, s.weights)
    nv = s.v.flat()[1]
    nu = len(s.poles)//nv
    if axis == 0:
        ctrl = [_value(s.u, [h[i*nv+m] for i in range(nu)], c) for m in range(nv)]
        basis = s.v
    else:
        ctrl = [_value(s.v, [h[i*nv+m] for m in range(nv)], c) for i in range(nu)]
        basis = s.u
    poles = [tuple(float(x/q[-1]) for x in q[:-1]) for q in ctrl]
    return BSpline3(basis, poles, [float(q[-1]) for q in ctrl], rng)


def surface_point(s, u, v):
    h = spl_homogeneous(s.poles, s.weights)
    nv = s.v.flat()[1]
    nu = len(s.poles)//nv
    q = _value(s.u, [_value(s.v, [h[i*nv+m] for m in range(nv)], v) for i in range(nu)], u)
    return tuple(float(x/q[-1]) for x in q[:-1])


def sheet_model(name, s, hole):
    """A face body on the whole domain of s, less the parameter rectangle
    `hole`: every edge an iso-curve (over a range on the hole), every
    pcurve a line; the outer loop counter-clockwise in (u, v), the hole's
    clockwise."""
    (a, b), (c, d) = (float(x) for x in s.u.flat()[2]), (float(x) for x in s.v.flat()[2])
    m = Model(name, 1e-7, kind='sheet')
    face = Face(s, True)

    def rectangle(u0, u1, v0, v1, outer):
        base = len(m.vertices)
        m.vertices += [surface_point(s, u, v) for u, v in ((u0, v0), (u1, v0), (u1, v1), (u0, v1))]
        rng = lambda lo, hi, lo0, hi0: None if (lo, hi) == (lo0, hi0) else (lo, hi)
        e = len(m.edges)
        m.edges += [Edge(base, base+1, iso_curve(s, 1, v0, rng(u0, u1, a, b))),
                    Edge(base+1, base+2, iso_curve(s, 0, u1, rng(v0, v1, c, d))),
                    Edge(base+3, base+2, iso_curve(s, 1, v1, rng(u0, u1, a, b))),
                    Edge(base, base+3, iso_curve(s, 0, u0, rng(v0, v1, c, d)))]
        loop = [Use(e, True, Line2((u0, v0), (u1, v0))), Use(e+1, True, Line2((u1, v0), (u1, v1))),
                Use(e+2, False, Line2((u1, v1), (u0, v1))), Use(e+3, False, Line2((u0, v1), (u0, v0)))]
        if not outer:
            loop = [Use(u.edge, not u.forward, Line2(u.pcurve.end, u.pcurve.start)) for u in reversed(loop)]
        face.loops.append(loop)
    rectangle(a, b, c, d, True)
    if hole is not None:
        rectangle(*hole, False)
    m.faces.append(face)
    return m


def spline_model(case):
    """The B-rep model of a T-b case (generate_brep_fixtures' builder)."""
    if case.variant is not None:
        models = {mm.name: mm for mm in gbf.spline_models()}
        m = copy.deepcopy(models[case.name])
    elif case.kind == 'prism':
        m = gbf.prism(case.name, case.boundaries, case.z0, case.z1)
    elif case.kind == 'face':
        ring = case.ring
        start = tuple(float(x) for x in spl.end_point(ring, 0))
        m = Model(case.name, 1e-7, kind='sheet')
        m.vertices = [(start[0], start[1], 0.0)]
        m.edges = [Edge(0, 0, BSpline3(ring.basis, [(x, y, 0.0) for x, y in ring.poles], ring.weights))]
        m.faces = [Face(Plane(Frame((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))), True,
                        [[Use(0, True, ring)]])]
    elif case.kind == 'dome':
        w, d, h = case.size
        m = gbf.prism(case.name, [gbf.rectangle(0.0, 0.0, w, d)], 0.0, h)
        # The top cap's pcurves are its (x, y), which are the graph's (u, v).
        m.faces[1].surface = case.surface
    else:
        m = sheet_model(case.name, case.surface, case.hole)
    if case.motion is not None:
        m = moved(m, case.motion)
    return declare(to_cell(m))


def spline_generate():
    blocks = []
    rows = ['# case\tkind\teuler\tloops\tarea\tvolume\tscale (tessellation_reference.py)']
    for c in spline_cases():
        cell = spline_model(c)
        issues = validate_cell(cell)
        assert not issues, (c.name, issues)
        rows_ = [f'mesh {label} {number(d)} {number(a)}' for label, d, a in settings(c)]
        head, end = encode_cell(cell).rsplit('\nend', 1)
        blocks.append(head+'\n'+'\n'.join(rows_)+'\nend'+end)
        b = ref.body(c)
        rows.append('\t'.join([c.name, 'solid' if b.solid else 'face', str(b.euler), str(b.loops),
                               mp.nstr(b.area, 17), mp.nstr(b.volume, 17), number(b.scale)]))
    return {'tessellation-spline-cases.txt': '\n'.join(blocks)+'\n',
            'tessellation-spline-expected.tsv': '\n'.join(rows)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = generate()
    for name, contents in files.items():
        path_ = ROOT/'fixtures'/name
        if args.check:
            if path_.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path_.write_text(contents)
    print(len(cases()), 'cases,', len(spline_cases()), 'spline cases')


if __name__ == '__main__':
    main()
