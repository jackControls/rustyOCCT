#!/usr/bin/env python3
"""Generic B-rep validator fixtures from an independent builder and oracle.

Valid solids come from an independent prism builder (lines, convex and concave
arcs, full circles with seams, holes and inverted cavity shells), optionally
rigidly rotated and translated. Each mutation targets one contract check; the
independent oracle in brep_reference.py computes the complete expected issue
list, including cascades. No Rust or native OCCT result supplies an expectation.
"""
import argparse
import copy
import math
from pathlib import Path

import mpmath as mp

from fractions import Fraction as F

from cell_reference import (body_class, declare, encode as encode_cell, face_mass_terms, gap_bounds,
                            mass_properties, polynomial_face,
                            to_cell, validate as validate_cell)
import spline_cell_reference as spline
from spline_cell_reference import Basis, BSpline2, BSpline3, BSplineSurface

from brep_reference import (Arc2, Arc3, Cone, Cylinder, Edge, Face, Frame, Line2, Line3, Model, Sphere, Torus,
                            Plane, TAU, Use, atan2_rn, cos_rn, encode, hypot_rn, number, sin_rn,
                            validate)

ROOT = Path(__file__).resolve().parents[1]
Z = (0.0, 0.0, 1.0)
X = (1.0, 0.0, 0.0)


def angle_of(p, c):
    return atan2_rn(p[1]-c[1], p[0]-c[0])


def prism(name, boundaries, z0=0.0, z1=1.0, tolerance=1e-7):
    """Boundaries: lists of pieces ('line', p) or ('arc', p, center, ccw) or
    ('circle', center, radius, ccw), or ('spline', p, interior poles, basis,
    weights) from p to the next piece's point, clamped (S4): its wall is the
    ruled spline surface C(u) + v z over v in [0, h]. Outer boundaries run
    counter-clockwise about +z and holes clockwise; material is always on
    the left."""
    m = Model(name, tolerance)
    h = z1-z0
    bottom = Face(Plane(Frame((0.0, 0.0, z0), (0.0, 0.0, -1.0), X)), True)
    top = Face(Plane(Frame((0.0, 0.0, z1), Z, X)), True)
    m.faces += [bottom, top]
    uv_bottom = lambda p: (p[0], -p[1])
    for boundary in boundaries:
        bottom_loop, top_loop = [], []
        if boundary[0][0] == 'circle':
            _, c, r, ccw = boundary[0]
            sweep = TAU if ccw else -TAU
            p = (c[0]+r, c[1])
            vb, vt = len(m.vertices), len(m.vertices)+1
            m.vertices += [(p[0], p[1], z0), (p[0], p[1], z1)]
            eb, et, es = len(m.edges), len(m.edges)+1, len(m.edges)+2
            m.edges += [Edge(vb, vb, Arc3(Frame((c[0], c[1], z0), Z, X), r, 0.0, sweep)),
                        Edge(vt, vt, Arc3(Frame((c[0], c[1], z1), Z, X), r, 0.0, sweep)),
                        Edge(vb, vt, Line3((p[0], p[1], z0), (p[0], p[1], z1)))]
            bottom_loop.append(Use(eb, False, Arc2(uv_bottom(c), r, -sweep, sweep)))
            top_loop.append(Use(et, True, Arc2(c, r, 0.0, sweep)))
            wall = Face(Cylinder(Frame((c[0], c[1], z0), Z, X), r), ccw)
            wall.loops.append([
                Use(eb, True, Line2((0.0, 0.0), (sweep, 0.0))),
                Use(es, True, Line2((sweep, 0.0), (sweep, h))),
                Use(et, False, Line2((sweep, h), (0.0, h))),
                Use(es, False, Line2((0.0, h), (0.0, 0.0))),
            ])
            m.faces.append(wall)
        else:
            n = len(boundary)
            points = [piece[1] for piece in boundary]
            vb = [len(m.vertices)+i for i in range(n)]
            m.vertices += [(p[0], p[1], z0) for p in points]
            vt = [len(m.vertices)+i for i in range(n)]
            m.vertices += [(p[0], p[1], z1) for p in points]
            vertical = [len(m.edges)+i for i in range(n)]
            m.edges += [Edge(vb[i], vt[i], Line3((p[0], p[1], z0), (p[0], p[1], z1))) for i, p in enumerate(points)]
            for i, piece in enumerate(boundary):
                j = (i+1) % n
                p, q = points[i], points[j]
                eb, et = len(m.edges), len(m.edges)+1
                if piece[0] == 'spline':
                    _, _, interior, basis, weights = piece
                    poles = [p, *interior, q]
                    weights = weights or [1.0]*len(poles)
                    at = lambda z: [(x, y, z) for x, y in poles]
                    m.edges += [Edge(vb[i], vb[j], BSpline3(basis, at(z0), weights)),
                                Edge(vt[i], vt[j], BSpline3(basis, at(z1), weights))]
                    flipped = BSpline2(basis, [uv_bottom(x) for x in poles], weights)
                    bottom_loop.append(Use(eb, False, spline.reversed_curve(flipped)))
                    top_loop.append(Use(et, True, BSpline2(basis, poles, weights)))
                    rows = [x for pole in poles for x in ((pole[0], pole[1], z0), (pole[0], pole[1], z1))]
                    wall = Face(BSplineSurface(basis, Basis(1, [0.0, h], [2, 2]), rows,
                                               [w for w in weights for _ in range(2)]), True)
                    a, b = basis.knots[0], basis.knots[-1]
                    wall.loops.append([
                        Use(eb, True, Line2((a, 0.0), (b, 0.0))),
                        Use(vertical[j], True, Line2((b, 0.0), (b, h))),
                        Use(et, False, Line2((b, h), (a, h))),
                        Use(vertical[i], False, Line2((a, h), (a, 0.0))),
                    ])
                elif piece[0] == 'line':
                    m.edges += [Edge(vb[i], vb[j], Line3((p[0], p[1], z0), (q[0], q[1], z0))),
                                Edge(vt[i], vt[j], Line3((p[0], p[1], z1), (q[0], q[1], z1)))]
                    length = hypot_rn(q[0]-p[0], q[1]-p[1])
                    t = ((q[0]-p[0])/length, (q[1]-p[1])/length, 0.0)
                    normal = (t[1], -t[0], 0.0)  # t x z
                    wall = Face(Plane(Frame((p[0], p[1], z0), normal, t)), True)
                    bottom_loop.append(Use(eb, False, Line2(uv_bottom(q), uv_bottom(p))))
                    top_loop.append(Use(et, True, Line2(p, q)))
                    wall.loops.append([
                        Use(eb, True, Line2((0.0, 0.0), (length, 0.0))),
                        Use(vertical[j], True, Line2((length, 0.0), (length, h))),
                        Use(et, False, Line2((length, h), (0.0, h))),
                        Use(vertical[i], False, Line2((0.0, h), (0.0, 0.0))),
                    ])
                else:
                    _, _, c, ccw = piece
                    r = hypot_rn(p[0]-c[0], p[1]-c[1])
                    a = angle_of(p, c)
                    b = angle_of(q, c)
                    sweep = (b-a) % TAU if ccw else -((a-b) % TAU)
                    m.edges += [Edge(vb[i], vb[j], Arc3(Frame((c[0], c[1], z0), Z, X), r, a, sweep)),
                                Edge(vt[i], vt[j], Arc3(Frame((c[0], c[1], z1), Z, X), r, a, sweep))]
                    bottom_loop.append(Use(eb, False, Arc2(uv_bottom(c), r, -(a+sweep), sweep)))
                    top_loop.append(Use(et, True, Arc2(c, r, a, sweep)))
                    wall = Face(Cylinder(Frame((c[0], c[1], z0), Z, X), r), sweep > 0)
                    wall.loops.append([
                        Use(eb, True, Line2((a, 0.0), (a+sweep, 0.0))),
                        Use(vertical[j], True, Line2((a+sweep, 0.0), (a+sweep, h))),
                        Use(et, False, Line2((a+sweep, h), (a, h))),
                        Use(vertical[i], False, Line2((a, h), (a, 0.0))),
                    ])
                m.faces.append(wall)
        bottom.loops.append(list(reversed(bottom_loop)))
        top.loops.append(top_loop)
    m.shells.append(list(range(len(m.faces))))
    return m


def rectangle(x0, y0, x1, y1, hole=False):
    pts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    if hole:
        pts.reverse()
    return [('line', p) for p in pts]


def regular(n, r, c=(0.0, 0.0), hole=False):
    pts = [(c[0]+r*cos_rn(TAU*k/n), c[1]+r*sin_rn(TAU*k/n)) for k in range(n)]
    if hole:
        pts.reverse()
    return [('line', p) for p in pts]


def stadium_boundary(length, r):
    return stadium(length, r)


def stadium(length, r):
    return [('line', (0.0, -r)), ('arc', (length, -r), (length, 0.0), True),
            ('line', (length, r)), ('arc', (0.0, r), (0.0, 0.0), True)]


def notch(w, h, r):
    """A rectangle with a concave semicircular notch in its top side."""
    return [('line', (0.0, 0.0)), ('line', (w, 0.0)), ('line', (w, h)),
            ('arc', (w/2+r, h), (w/2, h), False), ('line', (w/2-r, h)), ('line', (0.0, h))]


def invert(m):
    """Reverse every face and loop: a shell bounding the complement."""
    for f in m.faces:
        f.forward = not f.forward
        f.loops = [[Use(u.edge, not u.forward, reverse_pcurve(u.pcurve)) for u in reversed(loop)] for loop in f.loops]
    return m


def reverse_pcurve(p):
    if isinstance(p, BSpline2):
        return spline.flip(p)
    if isinstance(p, Line2):
        return Line2(p.end, p.start)
    return Arc2(p.center, p.radius, p.start+p.sweep, -p.sweep)


def merge(name, *models):
    out = Model(name, models[0].tolerance)
    for m in models:
        dv, de, df = len(out.vertices), len(out.edges), len(out.faces)
        out.vertices += m.vertices
        out.edges += [Edge(e.start+dv, e.end+dv, e.curve) for e in m.edges]
        out.faces += [Face(f.surface, f.forward, [[Use(u.edge+de, u.forward, u.pcurve) for u in l] for l in f.loops]) for f in m.faces]
        out.shells += [[fi+df for fi in s] for s in m.shells]
    return out


def transformed(m, name, axis, angle, shift):
    """Rigid rotation about a unit axis, then translation, of every 3D datum."""
    k = [a/math.sqrt(sum(b*b for b in axis)) for a in axis]
    c, s = cos_rn(angle), sin_rn(angle)
    def rot(v):
        kv = sum(a*b for a, b in zip(k, v))
        kx = (k[1]*v[2]-k[2]*v[1], k[2]*v[0]-k[0]*v[2], k[0]*v[1]-k[1]*v[0])
        return tuple(v[i]*c+kx[i]*s+k[i]*kv*(1-c) for i in range(3))
    point = lambda p: tuple(a+b for a, b in zip(rot(p), shift))
    frame = lambda f: Frame(point(f.origin), rot(f.normal), rot(f.x))
    out = copy.deepcopy(m)
    out.name = name
    out.vertices = [point(v) for v in m.vertices]
    for e in out.edges:
        c3 = e.curve
        e.curve = Line3(point(c3.start), point(c3.end)) if isinstance(c3, Line3) else Arc3(frame(c3.frame), c3.radius, c3.start, c3.sweep)
    for f in out.faces:
        f.surface = Plane(frame(f.surface.frame)) if isinstance(f.surface, Plane) else Cylinder(frame(f.surface.frame), f.surface.radius)
    return out


def base_cases():
    box = prism('box', [rectangle(0, 0, 3, 2)])
    cylinder = prism('cylinder', [[('circle', (0.0, 0.0), 1.5, True)]], 0.0, 2.0)
    plate_square = prism('plate_square_hole', [rectangle(0, 0, 4, 3), rectangle(1, 1, 2, 2, hole=True)], 0.0, 0.5)
    plate_round = prism('plate_round_hole', [rectangle(0, 0, 4, 3), [('circle', (2.0, 1.5), 0.75, False)]], 0.0, 0.5)
    plate_two = prism('plate_two_holes', [rectangle(0, 0, 6, 3), [('circle', (1.5, 1.5), 0.5, False)],
                                          regular(6, 0.8, (4.2, 1.5), hole=True)], -0.25, 0.25)
    tube = prism('tube', [[('circle', (0.0, 0.0), 2.0, True)], [('circle', (0.0, 0.0), 1.0, False)]], 0.0, 3.0)
    triangle = prism('triangle', [regular(3, 1.0)], 0.0, 0.7)
    hexagon = prism('hexagon', [regular(6, 2.0, (1.0, -1.0))], -1.0, 1.0)
    stadium_prism = prism('stadium', [stadium(3.0, 1.0)], 0.0, 1.0)
    notched = prism('concave_notch', [notch(4.0, 2.0, 0.75)], 0.0, 1.5)
    half = prism('half_disc', [[('line', (-1.0, 0.0)), ('arc', (1.0, 0.0), (0.0, 0.0), True)]], 0.0, 0.4)
    cavity = merge('box_cavity', prism('o', [rectangle(0, 0, 4, 4)], 0.0, 4.0),
                   invert(prism('i', [rectangle(1, 1, 3, 3)], 1.0, 3.0)))
    two_cavities = merge('box_two_cavities', prism('o', [rectangle(0, 0, 6, 4)], 0.0, 4.0),
                         invert(prism('i', [rectangle(1, 1, 2, 3)], 1.0, 3.0)),
                         invert(prism('j', [[('circle', (4.5, 2.0), 0.8, True)]], 1.0, 3.0)))
    cylinder_cavity = merge('cylinder_box_cavity', prism('o', [[('circle', (0.0, 0.0), 3.0, True)]], 0.0, 3.0),
                            invert(prism('i', [rectangle(-1, -1, 1, 1)], 1.0, 2.0)))
    cases = [box, cylinder, plate_square, plate_round, plate_two, tube, triangle, hexagon,
             stadium_prism, notched, half, cavity, two_cavities, cylinder_cavity]
    cases.append(transformed(box, 'box_rotated', (1, 2, 3), 0.7, (10.0, -5.0, 2.5)))
    cases.append(transformed(cylinder, 'cylinder_rotated', (-2, 1, 0.5), 2.1, (0.0, 0.0, 0.0)))
    cases.append(transformed(stadium_prism, 'stadium_rotated', (0, 1, 1), -1.2, (3.0, 3.0, 3.0)))
    cases.append(transformed(cavity, 'box_cavity_far', (1, 1, 1), 0.3, (1.0e4, -2.0e4, 5.0e3)))
    tiny = prism('box_small', [rectangle(0, 0, 3e-3, 2e-3)], 0.0, 1e-3, tolerance=1e-9)
    cases.append(tiny)
    return {c.name: c for c in cases}


def mutated(base, name, change):
    m = copy.deepcopy(base)
    m.name = name
    change(m)
    return m


def shift_pcurve(delta):
    def change(m):
        u = m.faces[2].loops[0][0]
        p = u.pcurve
        u.pcurve = Line2(tuple(a+delta for a in p.start), tuple(a+delta for a in p.end)) if isinstance(p, Line2) else \
            Arc2((p.center[0]+delta, p.center[1]+delta), p.radius, p.start, p.sweep)
    return change


def mutations(bases):
    box, cyl, plate, stad, cav = bases['box'], bases['cylinder'], bases['plate_round_hole'], bases['stadium'], bases['box_cavity']
    out = []
    add = lambda base, name, change: out.append(mutated(base, name, change))

    def drop_face(m):
        m.shells[0].remove(3)
    add(box, 'box_face_outside_shell', drop_face)

    def delete_face(m):
        del m.faces[3]
        m.shells = [[f if f < 3 else f-1 for f in s if f != 3] for s in m.shells]
    add(box, 'box_missing_face', delete_face)
    add(cyl, 'cylinder_missing_cap', lambda m: (m.faces.pop(1), m.shells.__setitem__(0, [0, 1])))

    def flip_use(m):
        u = m.faces[2].loops[0][1]
        u.forward = not u.forward
    add(box, 'box_flipped_use', flip_use)
    add(box, 'box_flipped_face', lambda m: setattr(m.faces[4], 'forward', not m.faces[4].forward))

    add(box, 'box_reversed_face', lambda m: invert_face(m, 4))
    add(cyl, 'cylinder_reversed_wall', lambda m: invert_face(m, 2))
    add(box, 'box_moved_vertex', lambda m: m.vertices.__setitem__(0, (1e-3, 0.0, 0.0)))
    add(stad, 'stadium_moved_vertex', lambda m: m.vertices.__setitem__(1, (3.0, -1.0+2e-4, 0.0)))
    add(box, 'box_pcurve_shift', shift_pcurve(1e-3))
    add(box, 'box_pcurve_shift_below_tolerance', shift_pcurve(1e-9))
    add(cyl, 'cylinder_cap_pcurve_shift', lambda m: shift_first(m, 1, 1e-3))
    add(stad, 'stadium_arc_pcurve_shift', lambda m: shift_first(m, 3, 5e-4))

    def swap_seam(m):
        loop = m.faces[2].loops[0]
        loop[1].pcurve, loop[3].pcurve = reverse_pcurve(loop[3].pcurve), reverse_pcurve(loop[1].pcurve)
    add(cyl, 'cylinder_swapped_seam_pcurves', swap_seam)

    def wrong_radius(m):
        u = m.faces[1].loops[0][0]
        u.pcurve = Arc2(u.pcurve.center, u.pcurve.radius*1.001, u.pcurve.start, u.pcurve.sweep)
    add(cyl, 'cylinder_cap_wrong_pcurve_radius', wrong_radius)
    add(box, 'box_face_twice', lambda m: m.shells[0].append(2))

    def third_use(m):
        m.faces.append(copy.deepcopy(m.faces[2]))
        m.shells[0].append(len(m.faces)-1)
    add(box, 'box_non_manifold_edges', third_use)
    two = merge('two_boxes_one_shell', bases['box'], prism('b', [rectangle(5, 0, 6, 1)]))
    two.shells = [two.shells[0]+two.shells[1]]
    out.append(two)
    add(cav, 'cavity_not_inverted', lambda m: invert_shell(m, 1))
    add(cav, 'outer_inverted', lambda m: invert_shell(m, 0))
    out.append(merge('cavity_outside', prism('o', [rectangle(0, 0, 4, 4)], 0.0, 4.0),
                     invert(prism('i', [rectangle(5, 1, 7, 3)], 1.0, 3.0))))
    out.append(merge('nested_cavities', prism('o', [rectangle(0, 0, 8, 8)], 0.0, 8.0),
                     invert(prism('i', [rectangle(1, 1, 7, 7)], 1.0, 7.0)),
                     invert(prism('j', [rectangle(3, 3, 5, 5)], 3.0, 5.0))))
    out.append(prism('hole_outside_outline', [rectangle(0, 0, 4, 3), [('circle', (6.0, 1.5), 0.5, False)]], 0.0, 0.5))
    def cavity_uses_outer_edge(m):
        m.faces[m.shells[1][2]].loops[0][0].edge = 4
    add(cav, 'cavity_face_uses_outer_edge', cavity_uses_outer_edge)

    def arc_sweep(m):
        e = next(e for e in m.edges if isinstance(e.curve, Arc3))
        e.curve = Arc3(e.curve.frame, e.curve.radius, e.curve.start, e.curve.sweep*0.999)
    add(stad, 'stadium_arc_sweep_changed', arc_sweep)
    add(box, 'box_extra_vertex', lambda m: m.vertices.append((9.0, 9.0, 9.0)))
    add(box, 'box_extra_edge', lambda m: m.edges.append(Edge(0, 6, Line3(m.vertices[0], m.vertices[6]))))
    add(box, 'box_empty_shell', lambda m: m.shells.append([]))
    add(box, 'box_empty_loop', lambda m: m.faces[2].loops.append([]))
    add(box, 'box_bad_edge_reference', lambda m: setattr(m.faces[2].loops[0][0], 'edge', 99))
    add(box, 'box_bad_face_reference', lambda m: m.shells[0].append(42))

    def zero_line(m):
        e = m.edges[0]
        e.curve = Line3(m.vertices[e.start], m.vertices[e.start])
    add(box, 'box_zero_length_edge', zero_line)
    add(cyl, 'cylinder_zero_radius_surface', lambda m: setattr(m.faces[2].surface, 'radius', 0.0))
    add(plate, 'plate_hole_loop_wound_wrong', lambda m: flip_loop_geometry(m, 0, 1))
    tol_loose = copy.deepcopy(box)
    tol_loose.name = 'box_loose_tolerance_absorbs_shift'
    tol_loose.tolerance = 1e-2
    shift_pcurve(1e-3)(tol_loose)
    out.append(tol_loose)
    return out


def shift_first(m, face, delta):
    u = m.faces[face].loops[0][0]
    p = u.pcurve
    u.pcurve = Arc2((p.center[0]+delta, p.center[1]), p.radius, p.start, p.sweep) if isinstance(p, Arc2) else \
        Line2((p.start[0]+delta, p.start[1]), (p.end[0]+delta, p.end[1]))


def invert_face(m, fi):
    f = m.faces[fi]
    f.forward = not f.forward
    f.loops = [[Use(u.edge, not u.forward, reverse_pcurve(u.pcurve)) for u in reversed(l)] for l in f.loops]


def invert_shell(m, si):
    for fi in m.shells[si]:
        invert_face(m, fi)


def flip_loop_geometry(m, fi, li):
    """Reverse only the pcurve winding of one loop, keeping uses and edges."""
    loop = m.faces[fi].loops[li]
    # Mirror the loop's pcurves through its centroid in u: winding flips, 3D does not.
    for u in loop:
        p = u.pcurve
        if isinstance(p, Arc2):
            u.pcurve = Arc2(p.center, p.radius, math.pi-p.start, -p.sweep)


def cone_cell(name, origin, normal, hint, r1, r2, h, tolerance=1e-7):
    """A cone or frustum as the kernel's cone builder makes it, in the cell
    model: the lateral face on a Cone surface (v along the generatrix from
    the base), ring edges bounding discs at nonzero ends, and a pole (a
    vertex loop at the apex) at a zero-radius end."""
    from cell_reference import Cell, CEdge, CFace, Fin as CFin, Loop as CLoop, Region, Shell
    c = Cell(name, tolerance)
    frame = Frame(origin, normal, hint)
    from brep_reference import axes
    o, x, y, n = [tuple(float(v) for v in a) for a in axes(frame)]
    length = hypot_rn(h, r2-r1)
    angle = atan2_rn(r2-r1, h)
    top = tuple(o[i]+h*n[i] for i in range(3))
    lateral = CFace(Cone(frame, r1, angle), True, [], 0, 1)
    faces = [lateral]
    for r, centre, up, v in ((r1, o, False, 0.0), (r2, top, True, length)):
        if r == 0:
            c.vertices.append(centre)
            c.loops.append(CLoop([], 0, len(c.vertices)-1))
            lateral.loops.append(len(c.loops)-1)
            continue
        e = len(c.edges)
        c.edges.append(CEdge(None, None, Arc3(Frame(centre, normal, hint), r, 0.0, TAU)))
        disc = CFace(Plane(Frame(centre, normal if up else tuple(-v for v in n), hint)), True, [], 0, 1)
        k = len(c.fins)
        if up:
            c.fins.append(CFin(e, True, Arc2((0.0, 0.0), r, 0.0, TAU)))
            c.fins.append(CFin(e, False, Line2((TAU, v), (0.0, v))))
        else:
            c.fins.append(CFin(e, False, Arc2((0.0, 0.0), r, -TAU, TAU)))
            c.fins.append(CFin(e, True, Line2((0.0, v), (TAU, v))))
        c.edges[e].fins = [k, k+1]
        c.loops.append(CLoop([k], 0))
        disc.loops.append(len(c.loops)-1)
        c.loops.append(CLoop([k+1], -1 if up else 1))
        lateral.loops.append(len(c.loops)-1)
        faces.append(disc)
    c.faces = faces
    c.shells = [Shell(1, [(f, 'F') for f in range(len(faces))]), Shell(0, [(f, 'B') for f in range(len(faces))])]
    c.regions = [Region('void', [1]), Region('solid', [0])]
    return c


HALF_PI = 1.5707963267948966


def sphere_cell(name, origin, normal, hint, radius, low, high, tolerance=1e-7):
    """A sphere or zone as the kernel's sphere builder makes it, in the cell
    model: the lateral face on a Sphere surface, ring edges at latitudes that
    are not poles bounding discs, a pole (a vertex loop) when exactly one end
    is one, and no loops at all for the whole sphere."""
    from cell_reference import Cell, CEdge, CFace, Fin as CFin, Loop as CLoop, Region, Shell
    from brep_reference import axes
    c = Cell(name, tolerance)
    frame = Frame(origin, normal, hint)
    o, x, y, n = [tuple(float(v) for v in a) for a in axes(frame)]
    lateral = CFace(Sphere(frame, radius), True, [], 0, 1)
    faces = [lateral]
    ends = [(low, low == -HALF_PI, False), (high, high == HALF_PI, True)]
    poles = sum(1 for _, pole, _ in ends if pole)
    for a, pole, up in ends:
        if pole:
            if poles == 1:
                c.vertices.append(tuple(o[i]+(radius if up else -radius)*n[i] for i in range(3)))
                c.loops.append(CLoop([], 0, len(c.vertices)-1))
                lateral.loops.append(len(c.loops)-1)
            continue
        z, r = radius*sin_rn(a), radius*cos_rn(a)
        centre = tuple(o[i]+z*n[i] for i in range(3))
        e = len(c.edges)
        c.edges.append(CEdge(None, None, Arc3(Frame(centre, normal, hint), r, 0.0, TAU)))
        disc = CFace(Plane(Frame(centre, normal if up else tuple(-v for v in n), hint)), True, [], 0, 1)
        k = len(c.fins)
        if up:
            c.fins.append(CFin(e, True, Arc2((0.0, 0.0), r, 0.0, TAU)))
            c.fins.append(CFin(e, False, Line2((TAU, a), (0.0, a))))
        else:
            c.fins.append(CFin(e, False, Arc2((0.0, 0.0), r, -TAU, TAU)))
            c.fins.append(CFin(e, True, Line2((0.0, a), (TAU, a))))
        c.edges[e].fins = [k, k+1]
        c.loops.append(CLoop([k], 0))
        disc.loops.append(len(c.loops)-1)
        c.loops.append(CLoop([k+1], -1 if up else 1))
        lateral.loops.append(len(c.loops)-1)
        faces.append(disc)
    c.faces = faces
    c.shells = [Shell(1, [(f, 'F') for f in range(len(faces))]), Shell(0, [(f, 'B') for f in range(len(faces))])]
    c.regions = [Region('void', [1]), Region('solid', [0])]
    return c


def torus_cell(name, origin, normal, hint, major, minor, low, high, angle, tolerance=1e-7):
    """A torus, v-segment or wedge as the kernel's torus builder makes it, in
    the cell model: the whole torus without loops; a v-segment's latitude
    rings bounding discs, its wall reversed when the meridian arc bulges
    toward the axis; a wedge's tube circles at u = 0 and u = angle bounding
    discs, its wall wound in v."""
    from cell_reference import Cell, CEdge, CFace, Fin as CFin, Loop as CLoop, Region, Shell
    from brep_reference import axes
    c = Cell(name, tolerance)
    frame = Frame(origin, normal, hint)
    o, x, y, n = [tuple(float(v) for v in a) for a in axes(frame)]
    closed, turn = high-low == TAU, angle == TAU
    # The meridian region between the arc and the axis runs counterclockwise
    # when its signed area (the integral of rho dz along the arc) is positive.
    area = mp.quad(lambda v: (major+minor*mp.cos(v))*minor*mp.cos(v), [low, high])
    lateral = CFace(Torus(frame, major, minor), closed or area > 0, [], 0, 1)
    faces = [lateral]
    if not closed:
        heights = (minor*sin_rn(low), minor*sin_rn(high))
        for k, a in enumerate((low, high)):
            z, rho = heights[k], major+minor*cos_rn(a)
            lower = z < heights[1-k]
            centre = tuple(o[i]+z*n[i] for i in range(3))
            e = len(c.edges)
            c.edges.append(CEdge(None, None, Arc3(Frame(centre, normal, hint), rho, 0.0, TAU)))
            disc = CFace(Plane(Frame(centre, tuple(-v for v in n) if lower else n, hint)), True, [], 0, 1)
            kk = len(c.fins)
            if lower:
                c.fins.append(CFin(e, False, Arc2((0.0, 0.0), rho, -TAU, TAU)))
            else:
                c.fins.append(CFin(e, True, Arc2((0.0, 0.0), rho, 0.0, TAU)))
            # +u at the low end on a forward wall; the wall's traversal
            # runs the other way round when it is reversed.
            plus = (k == 0) == lateral.forward
            if plus:
                c.fins.append(CFin(e, True, Line2((0.0, a), (TAU, a))))
            else:
                c.fins.append(CFin(e, False, Line2((TAU, a), (0.0, a))))
            c.edges[e].fins = [kk, kk+1]
            c.loops.append(CLoop([kk], 0))
            disc.loops.append(len(c.loops)-1)
            c.loops.append(CLoop([kk+1], 1 if plus else -1))
            lateral.loops.append(len(c.loops)-1)
            faces.append(disc)
    elif not turn:
        for u, end in ((0.0, False), (angle, True)):
            e_u = tuple(x[i]*cos_rn(u)+y[i]*sin_rn(u) for i in range(3))
            tangent = tuple(-x[i]*sin_rn(u)+y[i]*cos_rn(u) for i in range(3))
            centre = tuple(o[i]+major*e_u[i] for i in range(3))
            ring_normal = tuple(-t for t in tangent)
            e = len(c.edges)
            c.edges.append(CEdge(None, None, Arc3(Frame(centre, ring_normal, e_u), minor, 0.0, TAU)))
            disc = CFace(Plane(Frame(centre, tangent if end else ring_normal, e_u)), True, [], 0, 1)
            kk = len(c.fins)
            if end:
                c.fins.append(CFin(e, False, Arc2((0.0, 0.0), minor, -TAU, TAU)))
                c.fins.append(CFin(e, True, Line2((u, 0.0), (u, TAU))))
            else:
                c.fins.append(CFin(e, True, Arc2((0.0, 0.0), minor, 0.0, TAU)))
                c.fins.append(CFin(e, False, Line2((u, TAU), (u, 0.0))))
            c.edges[e].fins = [kk, kk+1]
            c.loops.append(CLoop([kk], 0))
            disc.loops.append(len(c.loops)-1)
            c.loops.append(CLoop([kk+1], 0, winding_v=1 if end else -1))
            lateral.loops.append(len(c.loops)-1)
            faces.append(disc)
    c.faces = faces
    c.shells = [Shell(1, [(f, 'F') for f in range(len(faces))]), Shell(0, [(f, 'B') for f in range(len(faces))])]
    c.regions = [Region('void', [1]), Region('solid', [0])]
    return c


def cell_cases(bases):
    """Cell-model cases with no seamed form: the model's own failure modes
    (TOPOLOGY_MODEL.md) and seamless valid shapes. They have no OCCT rows."""
    from cell_reference import Loop as CLoop
    out = []

    def cell(base, name, change):
        c = to_cell(copy.deepcopy(bases[base]))
        c.name = name
        change(c)
        out.append(c)

    def seam_at_pi(c):
        # Rotate every circle parametrization by pi: edges, caps and the wall.
        for e in c.edges:
            if isinstance(e.curve, Arc3):
                e.curve = Arc3(e.curve.frame, e.curve.radius, e.curve.start+math.pi, e.curve.sweep)
        for fin in c.fins:
            p = fin.pcurve
            if isinstance(p, Arc2):
                fin.pcurve = Arc2(p.center, p.radius, p.start+math.pi, p.sweep)
            else:
                fin.pcurve = Line2((p.start[0]+math.pi, p.start[1]), (p.end[0]+math.pi, p.end[1]))
    cell('cylinder', 'cylinder_seamless_at_pi', seam_at_pi)

    def shift_fin_period(c):
        wall = next(f for f in c.faces if isinstance(f.surface, Cylinder))
        k = c.loops[wall.loops[0]].fins[0]
        p = c.fins[k].pcurve
        c.fins[k].pcurve = Line2((p.start[0]+TAU, p.start[1]), (p.end[0]+TAU, p.end[1]))
    cell('stadium', 'stadium_wall_fin_shifted_by_period', shift_fin_period)
    cell('cylinder', 'cylinder_winding_flipped', lambda c: setattr(c.loops[c.faces[2].loops[0]], 'winding', -c.loops[c.faces[2].loops[0]].winding))

    def swap_fins(c):
        c.edges[0].fins, c.edges[1].fins = c.edges[1].fins, c.edges[0].fins
    cell('box', 'box_fins_swapped_between_edges', swap_fins)
    cell('box', 'box_side_wrong_shell', lambda c: setattr(c.faces[2], 'front', c.faces[2].back))

    def vertex_loop(offset):
        def change(c):
            top = c.faces[1]
            z = c.vertices[c.edges[c.fins[c.loops[top.loops[0]].fins[0]].edge].start][2]
            c.vertices.append((1.0, 1.0, z+offset))
            c.loops.append(CLoop([], 0, len(c.vertices)-1))
            top.loops.append(len(c.loops)-1)
        return change
    cell('box', 'box_vertex_loop', vertex_loop(0.0))
    cell('box', 'box_vertex_loop_off_surface', vertex_loop(1e-3))
    cell('box', 'box_face_without_loops', lambda c: setattr(c.faces[3], 'loops', []))

    def one_vertex_ring(c):
        c.vertices.append((1.5, 0.0, 0.0))
        c.edges[0].start = len(c.vertices)-1
    cell('cylinder', 'cylinder_ring_edge_with_one_vertex', one_vertex_ring)
    cell('box', 'box_shell_not_in_its_region', lambda c: c.regions[1].shells.remove(0))

    def wall_period(span):
        # A ring fin's pcurve spans `span` in u while its circle sweeps 2 pi,
        # as in OCCT files that print 2 pi as 6.28318530717959.
        def change(c):
            wall = next(f for f in c.faces if isinstance(f.surface, Cylinder))
            k = c.loops[wall.loops[0]].fins[0]
            p = c.fins[k].pcurve
            sign = 1 if p.end[0] > p.start[0] else -1
            c.fins[k].pcurve = Line2(p.start, (p.start[0]+sign*span, p.end[1]))
        return change
    cell('cylinder', 'cylinder_ring_pcurve_printed_period', wall_period(6.28318530717959))
    cell('cylinder', 'cylinder_ring_pcurve_period_off', wall_period(TAU*(1+1e-6)))

    # Enclosures (M5): gaps just inside and just outside the resolution, and
    # declared bounds that are missing, out of range or unsound.
    def moved_vertex(fraction, bound=None):
        def change(c):
            x, y, z = c.vertices[0]
            c.vertices[0] = (x+fraction*c.tolerance, y, z)
            if bound is not None:
                c.enclosures[('v', 0)] = bound*c.tolerance
        return change
    cell('box', 'box_vertex_moved_inside_resolution', moved_vertex(0.5))
    cell('box', 'box_vertex_moved_outside_resolution', moved_vertex(1.5))
    cell('box', 'box_vertex_enclosure_unsound', moved_vertex(0.5, 0.25))

    def cap_pcurve(fraction, bound=None):
        def change(c):
            cap = next(f for f in c.faces if isinstance(f.surface, Plane))
            k = c.loops[cap.loops[0]].fins[0]
            p = c.fins[k].pcurve
            c.fins[k].pcurve = Arc2((p.center[0]+fraction*c.tolerance, p.center[1]), p.radius, p.start, p.sweep)
            if bound is not None:
                c.enclosures[('u', k)] = bound*c.tolerance
        return change
    cell('cylinder', 'cylinder_cap_pcurve_inside_resolution', cap_pcurve(0.5))
    cell('cylinder', 'cylinder_cap_fin_enclosure_unsound', cap_pcurve(0.5, 0.25))

    def side_pcurve(fraction, face_bound=None):
        # One fin of a box side moved in v: its own deviation and both of its
        # face's junction gaps become `fraction` of the tolerance.
        def change(c):
            fi = 2
            k = c.loops[c.faces[fi].loops[0]].fins[0]
            p = c.fins[k].pcurve
            d = fraction*c.tolerance
            c.fins[k].pcurve = Line2((p.start[0], p.start[1]+d), (p.end[0], p.end[1]+d))
            if face_bound is not None:
                c.enclosures[('f', fi)] = face_bound*c.tolerance
        return change
    cell('box', 'box_side_pcurve_inside_resolution', side_pcurve(0.5))
    cell('box', 'box_face_enclosure_unsound', side_pcurve(0.5, 0.1))
    cell('box', 'box_fin_enclosure_missing', lambda c: c.enclosures.__setitem__(('u', 0), None))
    cell('box', 'box_face_enclosure_exceeds_resolution',
         lambda c: c.enclosures.__setitem__(('f', 0), 2*c.tolerance))
    cell('box', 'box_vertex_enclosure_negative', lambda c: c.enclosures.__setitem__(('v', 1), -c.tolerance))

    # Cones (S3): poles at either end, frustums, rotated and far copies, and
    # the pole's own failure modes.
    z, x = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
    cones = {
        'cone_apex': cone_cell('cone_apex', (0.0, 0.0, 0.0), z, x, 2.0, 0.0, 3.0),
        'cone_apex_at_base': cone_cell('cone_apex_at_base', (0.0, 0.0, 0.0), z, x, 0.0, 1.5, 2.0),
        'cone_frustum': cone_cell('cone_frustum', (1.0, 2.0, 3.0), z, x, 2.0, 1.0, 3.0),
        'cone_widening': cone_cell('cone_widening', (0.0, 0.0, 0.0), z, x, 1.0, 2.5, 0.5),
        'cone_rotated': cone_cell('cone_rotated', (0.5, -1.0, 2.0), (0.3, -0.4, 0.8), (1.0, 0.2, 0.0), 1.25, 0.0, 2.0),
        'cone_far': cone_cell('cone_far', (10000.0, -20000.0, 5000.0), z, x, 2.0, 0.5, 3.0),
    }
    out.extend(cones.values())

    def cone_case(base, name, change):
        c = copy.deepcopy(cones[base])
        c.name = name
        change(c)
        out.append(c)

    def pole(c):
        return next(l for l in c.loops if l.vertex is not None)

    def move_pole(offset):
        def change(c):
            v = pole(c).vertex
            c.vertices[v] = tuple(a+b for a, b in zip(c.vertices[v], offset))
        return change
    # Along the generatrix (still on the surface), and radially off it.
    cone_case('cone_apex', 'cone_pole_off_apex', move_pole((2.0*1e-3/3.605551275463989, 0.0, -3.0*1e-3/3.605551275463989)))
    cone_case('cone_apex', 'cone_pole_off_surface', move_pole((1e-3, 0.0, 0.0)))
    cone_case('cone_apex', 'cone_pole_missing',
              lambda c: c.faces[0].loops.remove(next(l for l in c.faces[0].loops if c.loops[l].vertex is not None)))
    cone_case('cone_frustum', 'cone_winding_twice', lambda c: setattr(c.loops[c.faces[0].loops[0]], 'winding', 2))
    cone_case('cone_apex', 'cone_right_angle', lambda c: setattr(c.faces[0].surface, 'half_angle', math.pi/2))

    def shift_lateral(c):
        k = c.loops[c.faces[0].loops[0]].fins[0]
        p = c.fins[k].pcurve
        c.fins[k].pcurve = Line2((p.start[0], p.start[1]+1e-3), (p.end[0], p.end[1]+1e-3))
    cone_case('cone_frustum', 'cone_lateral_pcurve_shift', shift_lateral)

    # Spheres (S3): the whole sphere (no loops), hemispheres with a pole at
    # either end, a zone, rotated and far copies, and the failure modes.
    spheres = {
        'sphere_whole': sphere_cell('sphere_whole', (0.0, 0.0, 0.0), z, x, 2.0, -HALF_PI, HALF_PI),
        'sphere_upper': sphere_cell('sphere_upper', (0.0, 0.0, 0.0), z, x, 1.5, 0.0, HALF_PI),
        'sphere_lower': sphere_cell('sphere_lower', (1.0, 2.0, 3.0), z, x, 1.5, -HALF_PI, 0.0),
        'sphere_zone': sphere_cell('sphere_zone', (0.0, 0.0, 0.0), z, x, 3.0, -0.5, 0.7),
        'sphere_rotated': sphere_cell('sphere_rotated', (0.5, -1.0, 2.0), (0.3, -0.4, 0.8), (1.0, 0.2, 0.0),
                                      1.25, -HALF_PI, 0.3),
        'sphere_far': sphere_cell('sphere_far', (10000.0, -20000.0, 5000.0), z, x, 2.0, -HALF_PI, HALF_PI),
    }
    out.extend(spheres.values())

    def sphere_case(base, name, change):
        c = copy.deepcopy(spheres[base])
        c.name = name
        change(c)
        out.append(c)
    # The pole at the other pole: on the surface, off the band's pole.
    sphere_case('sphere_upper', 'sphere_pole_wrong_side', move_pole((0.0, 0.0, -3.0)))
    sphere_case('sphere_upper', 'sphere_pole_off_surface', move_pole((0.0, 0.0, 1e-3)))
    sphere_case('sphere_upper', 'sphere_pole_missing',
                lambda c: c.faces[0].loops.remove(next(l for l in c.faces[0].loops if c.loops[l].vertex is not None)))
    sphere_case('sphere_zone', 'sphere_winding_twice', lambda c: setattr(c.loops[c.faces[0].loops[0]], 'winding', 2))
    sphere_case('sphere_whole', 'sphere_zero_radius', lambda c: setattr(c.faces[0].surface, 'radius', 0.0))
    sphere_case('sphere_whole', 'sphere_whole_reversed', lambda c: setattr(c.faces[0], 'forward', False))
    sphere_case('sphere_zone', 'sphere_ring_pcurve_shift', shift_lateral)

    def immersed_vertex(c):
        # A vertex loop on a whole sphere closes no band: an immersed vertex.
        from cell_reference import Loop as CLoop
        c.vertices.append((0.0, 0.0, 2.0))
        c.loops.append(CLoop([], 0, len(c.vertices)-1))
        c.faces[0].loops.append(len(c.loops)-1)
    sphere_case('sphere_whole', 'sphere_whole_vertex_loop', immersed_vertex)
    # A v winding on a sphere is structurally wrong.
    sphere_case('sphere_zone', 'sphere_winding_in_v', lambda c: setattr(c.loops[c.faces[0].loops[0]], 'winding_v', 1))

    # Tori (S3): the whole torus, the outer and inner halves (the inner one's
    # wall reversed), a general segment, wedges, rotated and far copies, and
    # the failure modes of windings in v.
    tori = {
        'torus_whole': torus_cell('torus_whole', (0.0, 0.0, 0.0), z, x, 3.0, 1.0, 0.0, TAU, TAU),
        'torus_outer_half': torus_cell('torus_outer_half', (0.0, 0.0, 0.0), z, x, 3.0, 1.0, -HALF_PI, HALF_PI, TAU),
        'torus_inner_half': torus_cell('torus_inner_half', (0.0, 0.0, 0.0), z, x, 3.0, 1.0, HALF_PI, 3*HALF_PI,
                                       TAU),
        'torus_segment': torus_cell('torus_segment', (1.0, 2.0, 3.0), z, x, 4.0, 1.5, 0.3, 2.0, TAU),
        'torus_wedge': torus_cell('torus_wedge', (0.0, 0.0, 0.0), z, x, 3.0, 1.0, 0.0, TAU, HALF_PI),
        'torus_wide_wedge': torus_cell('torus_wide_wedge', (0.5, -1.0, 2.0), (0.3, -0.4, 0.8), (1.0, 0.2, 0.0),
                                       5.0, 2.0, 0.0, TAU, 5.0),
        'torus_far': torus_cell('torus_far', (10000.0, -20000.0, 5000.0), z, x, 3.0, 1.0, 0.0, TAU, TAU),
    }
    out.extend(tori.values())

    def torus_case(base, name, change):
        c = copy.deepcopy(tori[base])
        c.name = name
        change(c)
        out.append(c)
    wall_loop = lambda c, k: c.loops[c.faces[0].loops[k]]
    torus_case('torus_wedge', 'torus_winding_v_twice', lambda c: setattr(wall_loop(c, 1), 'winding_v', 2))
    torus_case('torus_wedge', 'torus_winding_v_unbalanced', lambda c: setattr(wall_loop(c, 0), 'winding_v', 1))

    def shift_meridian(c):
        k = wall_loop(c, 0).fins[0]
        p = c.fins[k].pcurve
        c.fins[k].pcurve = Line2((p.start[0]+1e-3, p.start[1]), (p.end[0]+1e-3, p.end[1]))
    torus_case('torus_wedge', 'torus_meridian_pcurve_shift', shift_meridian)
    torus_case('torus_whole', 'torus_whole_reversed', lambda c: setattr(c.faces[0], 'forward', False))
    torus_case('torus_whole', 'torus_spindle', lambda c: setattr(c.faces[0].surface, 'minor', 3.5))
    torus_case('torus_inner_half', 'torus_inner_half_wall_forward', lambda c: setattr(c.faces[0], 'forward', True))
    out.extend(spline_cases(bases))
    return out


def exact(values):
    """Fractions as binary64, which they must be exactly."""
    out = tuple(float(x) for x in values)
    assert all(F(x) == y for x, y in zip(out, values)), values
    return out


def along(a, b, fractions):
    """Exact points a + f (b - a)."""
    a, b = [F(x) for x in a], [F(x) for x in b]
    return [exact([x+f*(y-x) for x, y in zip(a, b)]) for f in fractions]


# A clamped quadratic with one knot of multiplicity 2 at 1/2: its Greville
# fractions reproduce a line at uniform speed, so it is C1 there.
QUADRATIC = Basis(2, [0.0, 0.5, 1.0], [3, 2, 3])
QUADRATIC_LINE = [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]
# The same with the first interior pole moved along the line: the speed jumps.
QUADRATIC_KINK = [F(0), F(3, 8), F(1, 2), F(3, 4), F(1)]


def spline_cases(bases):
    """R4 of REVIEW_NOTES.md: spline edges, pcurves and faces, C1 or not in
    their own parameterisation. Every other check of spline geometry is
    uncertified before the rest of S4, so none of these is valid."""
    out = []

    def cell(base, name, change):
        c = to_cell(copy.deepcopy(bases[base]))
        c.name = name
        change(c)
        out.append(c)

    def edge_spline(k, fractions, basis=QUADRATIC, weights=None):
        def change(c):
            e = c.edges[k]
            a, b = c.vertices[e.start], c.vertices[e.end]
            poles = along(a, b, fractions)
            e.curve = BSpline3(basis, poles, weights or [1.0]*len(poles))
        return change
    cell('box', 'spline_edge_c1', edge_spline(0, QUADRATIC_LINE))
    cell('box', 'spline_edge_not_c1', edge_spline(0, QUADRATIC_KINK))
    # Degree 1: a knot of multiplicity 1 is tested; the midpoint keeps the
    # speed, a point a quarter along does not.
    linear = Basis(1, [0.0, 0.5, 1.0], [2, 1, 2])
    cell('box', 'spline_edge_linear_c1', edge_spline(0, [F(0), F(1, 2), F(1)], linear))
    cell('box', 'spline_edge_linear_kink', edge_spline(0, [F(0), F(1, 4), F(1)], linear))
    # Unclamped degree 1 over [1, 2]: the knots 0 and 3 lie beyond the
    # domain and are not tested.
    cell('box', 'spline_edge_unclamped', edge_spline(0, [F(0), F(1)], Basis(1, [0.0, 1.0, 2.0, 3.0], [1, 1, 1, 1])))

    # Rational: C1 of the homogeneous curve. With weights 1, 2, 3 around the
    # knot's pole, w2 P2 = (w1 P1 + w3 P3) / 2 keeps it; P2 on the line
    # otherwise breaks it.
    def rational(c1):
        def change(c):
            e = c.edges[0]
            a, b = [F(x) for x in c.vertices[e.start]], [F(x) for x in c.vertices[e.end]]
            p1 = [x+(y-x)/4 for x, y in zip(a, b)]
            p3 = [x+3*(y-x)/4 for x, y in zip(a, b)]
            p2 = ([(x+3*y)/4 for x, y in zip(p1, p3)] if c1
                  else [(x+y)/2 for x, y in zip(a, b)])
            poles = [exact(a), exact(p1), exact(p2), exact(p3), exact(b)]
            e.curve = BSpline3(QUADRATIC, poles, [1.0, 1.0, 2.0, 3.0, 1.0])
        return change
    cell('box', 'spline_edge_rational_c1', rational(True))
    cell('box', 'spline_edge_rational_not_c1', rational(False))
    cell('box', 'spline_edge_degenerate', edge_spline(0, [F(0)]*5))

    # Ranges (S4e): a line spline three times the edge's length, over the
    # range [1, 2] that is the edge itself (valid); with a corner knot at
    # 0.5 outside the range, neither tested nor used (valid); and a cap's
    # pcurve likewise.
    long_basis = Basis(2, [0.0, 1.0, 2.0, 3.0], [3, 1, 1, 3])
    long_line = [F(-1), F(-1, 2), F(1, 2), F(3, 2), F(2)]
    corner_basis = Basis(2, [0.0, 0.5, 1.0, 2.0, 3.0], [3, 2, 1, 1, 3])
    corner_line = [F(-1), F(-3, 4), F(-1, 2), F(-1, 4), F(1, 2), F(3, 2), F(2)]

    def ranged_edge(fractions, basis, kink=False):
        def change(c):
            e = c.edges[0]
            poles = along(c.vertices[e.start], c.vertices[e.end], fractions)
            if kink:
                x, y, z = poles[2]
                poles[2] = (x, y+0.25, z)
            e.curve = BSpline3(basis, poles, [1.0]*len(poles), (1.0, 2.0))
        return change
    cell('box', 'spline_edge_range', ranged_edge(long_line, long_basis))
    cell('box', 'spline_edge_range_corner_outside', ranged_edge(corner_line, corner_basis, kink=True))

    def ranged_pcurve(c):
        k = c.loops[c.faces[1].loops[0]].fins[0]
        p = c.fins[k].pcurve
        c.fins[k].pcurve = BSpline2(long_basis, along(p.start, p.end, long_line), [1.0]*5, (1.0, 2.0))
    cell('box', 'spline_pcurve_range', ranged_pcurve)

    # A reversed span: the line spline from the fin's end to its start,
    # traversed backwards over a range (the importer's reversed pcurves).
    def flagged_pcurve(c):
        k = c.loops[c.faces[1].loops[0]].fins[0]
        p = c.fins[k].pcurve
        back = BSpline2(long_basis, along(p.end, p.start, long_line), [1.0]*5, (1.0, 2.0))
        c.fins[k].pcurve = spline.flip(back)
    cell('box', 'spline_pcurve_reversed', flagged_pcurve)

    # A ring edge: a periodic quadratic of period 3 with every knot of
    # multiplicity 2; its poles 0, 2 and 4 lie on the curve (at 0, 1, 2),
    # each the midpoint of its neighbours for C1. Moving pole 0 breaks only
    # the seam. A nonperiodic spline cannot be a ring edge.
    def ring(seam_c1=True, periodic=True):
        def change(c):
            k = next(i for i, e in enumerate(c.edges) if e.start is None)
            z = F(c.edges[k].curve.frame.origin[2])
            corners = [(F(2), F(0)), (F(-1), F(2)), (F(-1), F(-2))]
            mid = lambda p, q: ((p[0]+q[0])/2, (p[1]+q[1])/2)
            ring = []
            for i in range(3):
                ring += [mid(corners[i-1], corners[i]), corners[i]]
            if not seam_c1:
                ring[0] = (ring[0][0]+F(1, 8), ring[0][1])
            poles = [exact((x, y, z)) for x, y in ring]
            basis = Basis(2, [0.0, 1.0, 2.0, 3.0], [2, 2, 2, 2], periodic=True)
            if not periodic:
                poles.append(poles[0])
                basis = Basis(2, [0.0, 1.0, 2.0, 3.0], [3, 2, 2, 3])
            c.edges[k].curve = BSpline3(basis, poles, [1.0]*len(poles))
        return change
    cell('cylinder', 'spline_ring_edge_c1', ring())
    cell('cylinder', 'spline_ring_edge_seam_not_c1', ring(seam_c1=False))
    cell('cylinder', 'spline_ring_edge_nonperiodic', ring(periodic=False))

    # A periodic quadratic of three poles with its seam of multiplicity 2:
    # one removal would leave fewer poles than the basis allows, so the
    # kernel refines elsewhere first. Pole 0 is the seam's point, C1 as the
    # midpoint of the others.
    def small_ring(c1):
        def change(c):
            k = next(i for i, e in enumerate(c.edges) if e.start is None)
            z = F(c.edges[k].curve.frame.origin[2])
            p1, p2 = (F(2), F(0)), (F(-2), F(2))
            p0 = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2) if c1 else (F(1, 8), F(1))
            poles = [exact((x, y, z)) for x, y in (p0, p1, p2)]
            basis = Basis(2, [0.0, 1.0, 2.0], [2, 1, 2], periodic=True)
            c.edges[k].curve = BSpline3(basis, poles, [1.0]*3)
        return change
    cell('cylinder', 'spline_ring_edge_small_basis_c1', small_ring(True))
    cell('cylinder', 'spline_ring_edge_small_basis_not_c1', small_ring(False))

    # A cap's pcurve as a spline over the same segment.
    def pcurve_spline(fractions):
        def change(c):
            k = c.loops[c.faces[1].loops[0]].fins[0]
            p = c.fins[k].pcurve
            c.fins[k].pcurve = BSpline2(QUADRATIC, along(p.start, p.end, fractions), [1.0]*5)
        return change
    cell('box', 'spline_pcurve_c1', pcurve_spline(QUADRATIC_LINE))
    cell('box', 'spline_pcurve_not_c1', pcurve_spline(QUADRATIC_KINK))
    cell('box', 'spline_pcurve_degenerate', pcurve_spline([F(0)]*5))

    # The top cap's surface as a spline: biquadratic, a knot of
    # multiplicity 2 at 1.5 in u, its rows lines at uniform speed (C1, and
    # the plane itself, so the cap stays valid with its line pcurves across
    # the knot line) or one row kinked. A vertex loop on it is uncertified.
    def surface_spline(kinked=False, vertex_loop=False):
        def change(c):
            fractions = QUADRATIC_LINE
            rows = []
            for i, fu in enumerate(fractions):
                for fv in (F(0), F(1, 2), F(1)):
                    if kinked and i == 1 and fv == 0:
                        fu = QUADRATIC_KINK[1]
                    rows.append(exact((3*fu, 2*fv, F(1))))
            # Over the cap's own (x, y): u in [0, 3] with the knot 1.5, v in
            # [0, 2]; the Greville poles make S(u, v) = (u, v, 1) exactly.
            c.faces[1].surface = BSplineSurface(Basis(2, [0.0, 1.5, 3.0], [3, 2, 3]), Basis(2, [0.0, 2.0], [3, 3]),
                                                rows, [1.0]*len(rows))
            if vertex_loop:
                from cell_reference import Loop as CLoop
                c.vertices.append((1.0, 1.0, 1.0))
                c.loops.append(CLoop([], 0, len(c.vertices)-1))
                c.faces[1].loops.append(len(c.loops)-1)
        return change
    cell('box', 'spline_face_c1', surface_spline())
    cell('box', 'spline_face_not_c1', surface_spline(kinked=True))
    cell('box', 'spline_face_vertex_loop', surface_spline(vertex_loop=True))

    # Containment against a spline side (S4c): a hole inside the bulge's
    # spline region, and one just outside it.
    quadratic = Basis(2, [0.0, 1.0], [3, 3])
    bulge = [('line', (0.0, 0.0)), ('spline', (3.0, 0.0), [(4.0, 1.0)], quadratic, None),
             ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]
    for name, x0 in (('spline_bulge_hole_inside', 3.125), ('spline_bulge_hole_outside', 3.625)):
        m = prism(name, [bulge, rectangle(x0, 0.875, x0+0.25, 1.125, hole=True)])
        out.append(to_cell(m))
    # S4d: the rational rounded corner shrunk by 2^-10 and moved to about
    # (-7.6, 4.5, 0.06), exactly (a brep_validation fuzz input): enclosed
    # from absolute rational jets, its flux's sign was undecided.
    k, (tx, ty, tz) = 2.0**-10, (-7.5625, 4.5, 0.0625)
    at = lambda p: (tx+k*p[0], ty+k*p[1])
    far = [('line', at((0.0, 0.0))), ('line', at((3.0, 0.0))),
           ('spline', at((3.0, 1.5)), [at((3.0, 2.0))], quadratic, [1.0, 0.7071067811865476, 1.0]),
           ('line', at((2.5, 2.0))), ('line', at((0.0, 2.0)))]
    out.append(to_cell(prism('spline_rounded_corner_far', [far], tz-0.5*k, tz+0.5*k, tolerance=1e-7*k)))
    return out


def spline_models():
    """S4: seamed models with spline edges, pcurves and walls and their
    explicit OCCT rows, for the native spline capture and bridge
    (compare_brep.py --family spline). Kept apart from the prism models,
    whose native inputs are pinned by the pre-implementation captures."""
    quadratic = Basis(2, [0.0, 1.0], [3, 3])
    bulge = [('line', (0.0, 0.0)), ('spline', (3.0, 0.0), [(4.0, 1.0)], quadratic, None),
             ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]
    cubic = Basis(3, [0.0, 0.5, 1.0], [4, 1, 4])
    cubic_bulge = [('line', (0.0, 0.0)),
                   ('spline', (3.0, 0.0), [(3.5, 0.25), (4.0, 1.0), (3.5, 1.75)], cubic, None),
                   ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]
    # A knot of multiplicity 2 on a quadratic whose poles around it are not
    # collinear at equal spacing: C0 there, a corner the kernel rejects.
    corner = Basis(2, [0.0, 0.5, 1.0], [3, 2, 3])
    c0_bulge = [('line', (0.0, 0.0)),
                ('spline', (3.0, 0.0), [(3.5, 0.25), (3.75, 1.0), (3.5, 1.75)], corner, None),
                ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]
    # A rounded corner: the rational quadratic of a quarter circle.
    rounded = [('line', (0.0, 0.0)), ('line', (3.0, 0.0)),
               ('spline', (3.0, 1.5), [(3.0, 2.0)], quadratic, [1.0, 0.7071067811865476, 1.0]),
               ('line', (2.5, 2.0)), ('line', (0.0, 2.0))]
    models = [prism('spline_bulge', [bulge]), prism('spline_cubic_bulge', [cubic_bulge]),
              prism('spline_c0_bulge', [c0_bulge]), prism('spline_rounded_corner', [rounded], -0.5, 0.5)]
    base = {m.name: m for m in models}

    def mutated_model(name, source, change):
        m = copy.deepcopy(base[source])
        m.name = name
        change(m)
        models.append(m)

    def shift_top_pcurve(m):
        top = m.faces[1]
        u = next(u for u in top.loops[0] if isinstance(u.pcurve, BSpline2))
        u.pcurve = BSpline2(u.pcurve.basis, [(x+1e-3, y) for x, y in u.pcurve.poles], u.pcurve.weights)
    mutated_model('spline_bulge_pcurve_shift', 'spline_bulge', shift_top_pcurve)
    mutated_model('spline_bulge_vertex_moved', 'spline_bulge',
                  lambda m: m.vertices.__setitem__(1, (3.0+1e-3, 0.0, 0.0)))
    mutated_model('spline_bulge_wall_reversed', 'spline_bulge',
                  lambda m: setattr(next(f for f in m.faces if isinstance(f.surface, BSplineSurface)), 'forward', False))

    # A stadium's vertical edge and its pcurve on the cylinder as degree-1
    # splines: spline geometry on a periodic analytic surface.
    linear = Basis(1, [0.0, 1.0], [2, 2])

    def cylinder_use(m):
        wall = next(f for f in m.faces if isinstance(f.surface, Cylinder))
        return next(u for u in wall.loops[0] if m.edges[u.edge].start != m.edges[u.edge].end
                    and isinstance(m.edges[u.edge].curve, Line3)
                    and m.edges[u.edge].curve.start[:2] == m.edges[u.edge].curve.end[:2])

    def spline_pcurve(du=0.0):
        def change(m):
            u = cylinder_use(m)
            u.pcurve = BSpline2(linear, [(u.pcurve.start[0]+du, u.pcurve.start[1]),
                                         (u.pcurve.end[0]+du, u.pcurve.end[1])], [1.0, 1.0])
        return change

    def spline_edge(m):
        e = m.edges[cylinder_use(m).edge]
        e.curve = BSpline3(linear, [e.curve.start, e.curve.end], [1.0, 1.0])
    stadium = prism('spline_stadium', [stadium_boundary(3.0, 1.0)])
    base['spline_stadium'] = stadium
    mutated_model('spline_stadium_pcurve', 'spline_stadium', spline_pcurve())
    mutated_model('spline_stadium_pcurve_shift', 'spline_stadium', spline_pcurve(1e-3))
    mutated_model('spline_stadium_edge', 'spline_stadium', spline_edge)
    return models


def quadrature_models():
    """F8: spline models for routes of the certified quadrature the S4
    models do not exercise, captured natively before the kernel code
    (compare_brep.py --family spline): a spline pcurve along a cylinder's
    parallel, where the Green integral runs with du != 0, and spline walls
    with a knot in v, whose upper patch has a column below it."""
    base = {m.name: m for m in spline_models()}
    # A narrow stadium at tolerance 1e-6: the kernel measures an arc edge's
    # use by a spline pcurve along a parallel (M5) by second-order Taylor
    # bounds on at most 2^8 pieces, about r (pi/256)^3/4 for a half turn, 4.6e-7
    # at r = 1: above 1e-7, and above the bound the reference declares (it
    # refines to a quarter of the tolerance). Validation itself refines
    # further and certifies the use at 1e-7 either way.
    base['spline_stadium'] = prism('spline_stadium', [stadium_boundary(3.0, 0.25)], tolerance=1e-6)
    models = []
    linear = Basis(1, [0.0, 1.0], [2, 2])

    def parallel(m):
        # The first cylinder wall's top use, (a + sweep, h) to (a, h).
        wall = next(f for f in m.faces if isinstance(f.surface, Cylinder))
        u = wall.loops[0][2]
        assert isinstance(u.pcurve, Line2) and u.pcurve.start[1] == u.pcurve.end[1] != 0.0
        u.pcurve = BSpline2(linear, [u.pcurve.start, u.pcurve.end], [1.0, 1.0])

    def split_wall(m):
        wall = next(f for f in m.faces if isinstance(f.surface, BSplineSurface))
        s = wall.surface
        h = s.v.knots[-1]
        poles = [s.poles[2*i] for i in range(len(s.poles)//2)]
        tops = [s.poles[2*i+1] for i in range(len(s.poles)//2)]
        rows = [x for p, q in zip(poles, tops) for x in (p, (p[0], p[1], (p[2]+q[2])/2), q)]
        weights = [s.weights[2*i] for i in range(len(s.poles)//2)]
        wall.surface = BSplineSurface(s.u, Basis(1, [0.0, h/2, h], [2, 1, 2]), rows,
                                      [w for w in weights for _ in range(3)])

    for name, source, change in [('spline_stadium_parallel', 'spline_stadium', parallel),
                                 ('spline_bulge_split_wall', 'spline_bulge', split_wall),
                                 ('spline_rounded_corner_split_wall', 'spline_rounded_corner', split_wall)]:
        m = copy.deepcopy(base[source])
        m.name = name
        change(m)
        models.append(m)
    return models


def extract(m, name, face_ids, kind):
    """S6: the faces `face_ids` of a model as a sheet ('sheet': a free face,
    or an open shell of several) or a closed shell without a solid
    ('shell'), with the edges and vertices they use."""
    edges = sorted({u.edge for fi in face_ids for loop in m.faces[fi].loops for u in loop})
    emap = {e: i for i, e in enumerate(edges)}
    verts = sorted({v for e in edges for v in (m.edges[e].start, m.edges[e].end)})
    vmap = {v: i for i, v in enumerate(verts)}
    out = Model(name, m.tolerance, kind=kind)
    out.vertices = [m.vertices[v] for v in verts]
    out.edges = [Edge(vmap[m.edges[e].start], vmap[m.edges[e].end], m.edges[e].curve) for e in edges]
    for fi in face_ids:
        f = m.faces[fi]
        out.faces.append(Face(f.surface, f.forward,
                              [[Use(emap[u.edge], u.forward, u.pcurve) for u in loop] for loop in f.loops]))
    out.shells = [] if kind == 'sheet' and len(face_ids) == 1 else [list(range(len(face_ids)))]
    return out


def wire_of(m, name, edge_ids):
    """S6: the edges `edge_ids` of a model, in order, as a wire body."""
    verts = sorted({v for e in edge_ids for v in (m.edges[e].start, m.edges[e].end)})
    vmap = {v: i for i, v in enumerate(verts)}
    out = Model(name, m.tolerance, kind='wire')
    out.vertices = [m.vertices[v] for v in verts]
    out.edges = [Edge(vmap[m.edges[e].start], vmap[m.edges[e].end], m.edges[e].curve) for e in edge_ids]
    out.wire = list(range(len(edge_ids)))
    return out


def sheet_models():
    """S6: sheets (free faces and open shells), closed shells without a
    solid, wires and an acorn, cut from the neutral prisms, with their
    explicit OCCT rows (compare_brep.py --family sheet), and mutations."""
    box = prism('box', [rectangle(0, 0, 3, 2)])
    plate = prism('plate', [rectangle(0, 0, 4, 3), [('circle', (2.0, 1.5), 0.75, False)]], 0.0, 0.5)
    stadium_prism = prism('stadium', [stadium(3.0, 1.0)], 0.0, 1.0)
    half = prism('half_disc', [[('line', (-1.0, 0.0)), ('arc', (1.0, 0.0), (0.0, 0.0), True)]], 0.0, 0.4)
    cylinder = prism('cylinder', [[('circle', (0.0, 0.0), 1.5, True)]], 0.0, 2.0)
    quadratic = Basis(2, [0.0, 1.0], [3, 3])
    bulge = prism('bulge', [[('line', (0.0, 0.0)), ('spline', (3.0, 0.0), [(4.0, 1.0)], quadratic, None),
                             ('line', (3.0, 2.0)), ('line', (0.0, 2.0))]])
    models = [
        extract(box, 'sheet_square', [1], 'sheet'),
        extract(box, 'sheet_square_bottom', [0], 'sheet'),
        extract(plate, 'sheet_round_hole', [1], 'sheet'),
        extract(stadium_prism, 'sheet_stadium', [1], 'sheet'),
        extract(half, 'sheet_half_disc', [1], 'sheet'),
        # A partial cylinder: the stadium's first arc wall.
        extract(stadium_prism, 'sheet_cylinder_patch', [3], 'sheet'),
        # The spline wall of the bulge.
        extract(bulge, 'sheet_spline_wall', [3], 'sheet'),
        # An open box: every face but the top.
        extract(box, 'sheet_open_box', [0, 2, 3, 4, 5], 'sheet'),
        # Closed shells bounding a void.
        extract(box, 'shell_box', [0, 1, 2, 3, 4, 5], 'shell'),
        extract(cylinder, 'shell_cylinder', [0, 1, 2], 'shell'),
        # Wires: the box's bottom edges, a circle, a line-arc-line run.
        wire_of(box, 'wire_square', [4, 6, 8, 10]),
        wire_of(cylinder, 'wire_circle', [0]),
        wire_of(stadium_prism, 'wire_open', [4, 6, 8]),
        wire_of(box, 'wire_edge', [4]),
    ]
    acorn = Model('acorn_point', 1e-7, kind='acorn')
    acorn.vertices = [(1.0, 2.0, 3.0)]
    models.append(acorn)
    base = {m.name: m for m in models}

    def mutated(name, source, change):
        m = copy.deepcopy(base[source])
        m.name = name
        change(m)
        models.append(m)
    mutated('sheet_square_vertex_moved', 'sheet_square',
            lambda m: m.vertices.__setitem__(0, (m.vertices[0][0]+1e-3, m.vertices[0][1], m.vertices[0][2])))

    def shift(m):
        u = m.faces[0].loops[0][0]
        u.pcurve = Line2((u.pcurve.start[0]+1e-3, u.pcurve.start[1]), (u.pcurve.end[0]+1e-3, u.pcurve.end[1]))
    mutated('sheet_open_box_pcurve_shift', 'sheet_open_box', shift)
    # Two edges that share no vertex.
    mutated('wire_disconnected', 'wire_square', lambda m: setattr(m, 'wire', [0, 2]))
    return models


# F8: the reference's own quadrature error bound, relative to the scale of
# each property (volume, area, centroid, inertia).
REFINED = mp.mpf('1e-20')
# Closed forms of spline models' properties (occt-spline-properties/NOTES.md):
# the bulge's profile is the rectangle 3 x 2 and a parabolic segment of area
# 2/3 whose arc length is sqrt 2 + asinh 1. (The rounded corner's weight is
# the binary64 value of sqrt(2)/2, so its arc is not exactly circular and the
# closed forms of a quarter disc hold only to about 1e-17.)
CLOSED_FORMS = {
    'spline_bulge': {0: mp.mpf(20)/3, 1: mp.mpf(40)/3+8+mp.sqrt(2)+mp.asinh(1)},
    'spline_bulge_split_wall': {0: mp.mpf(20)/3, 1: mp.mpf(40)/3+8+mp.sqrt(2)+mp.asinh(1)},
}


def spline_mass(c):
    """The mass properties of a valid spline case (S4d) as a row: volume,
    area, centroid, inertia about it; None when not integrated. F8: the
    case is integrated again with every Gauss-Legendre interval halved on
    the faces whose integrands are not polynomials the 24 nodes integrate
    exactly, and both must agree within REFINED of each property's scale;
    closed forms, where known, must hold to the same bound."""
    terms = {}

    def first(c, f, ref):
        terms[id(f)] = face_mass_terms(c, f, ref)
        return terms[id(f)]

    def halved(c, f, ref):
        return terms[id(f)] if polynomial_face(c, f) else face_mass_terms(c, f, ref, split=2)
    rows = []
    for terms_of in (first, halved):
        props = mass_properties(c, terms_of)
        if props is None:
            return None
        rows.append([props['volume'], props['area'], *props['centroid'], *sum(props['inertia'], [])])
    values, again = rows
    for group in (range(0, 1), range(1, 2), range(2, 5), range(5, 14)):
        scale = max(abs(values[i]) for i in group) or 1
        for i in group:
            assert abs(values[i]-again[i]) <= REFINED*scale, (c.name, i, values[i], again[i])
            if i in CLOSED_FORMS.get(c.name, {}):
                assert abs(values[i]-CLOSED_FORMS[c.name][i]) <= REFINED*scale, (c.name, i, values[i])
    return values


def generate():
    bases = base_cases()
    models = list(bases.values())+mutations(bases)
    cells = [declare(c) for c in [to_cell(m) for m in models]+cell_cases(bases)
             + [to_cell(m) for m in spline_models()]+[to_cell(m) for m in sheet_models()]
             + [to_cell(m) for m in quadrature_models()]]
    names = [c.name for c in cells]
    assert len(names) == len(set(names)), 'duplicate case names'
    text = '\n'.join(encode_cell(c) for c in cells)+'\n'
    rows, lows, masses, classes = [], [], [], []
    for c in cells:
        issues = validate_cell(c)
        rows.append(c.name+'\t'+';'.join(f'{k}:{e}' for k, e in issues))
        if not issues:
            classes.append(c.name+'\t'+body_class(c))
        if not issues and c.name.startswith('spline_'):
            # S4d: the mass properties of every valid spline case.
            values = spline_mass(c)
            if values is not None:
                # Quadrature noise far below any enclosure's width reads 0.
                masses.append(c.name+'\t'+' '.join('0' if abs(x) < 1e-25 else mp.nstr(x, 20)
                                                     for x in values))
        if not issues:
            # Certain lower values of every gap of a valid case, rounded down:
            # a measured enclosure below one is unsound.
            for (kind, i), (low, _) in sorted(gap_bounds(c).items()):
                value = float(low)
                if value > low:
                    value = math.nextafter(value, -math.inf)
                lows.append(f'{c.name}\t{kind} {i}\t{number(value)}')
    return models, {'brep-cases.txt': text,
                    'brep-expected.tsv': '# name\tsorted issues kind:entity separated by ;\n'+'\n'.join(rows)+'\n',
                    'brep-enclosure-lows.tsv': '# valid case\tv|u|f index (vertex, fin arena index, face)'
                    '\tcertain lower value of its gap\n'+'\n'.join(lows)+'\n',
                    'brep-spline-mass.tsv': '# valid spline case\tvolume, area, centroid, inertia about it'
                    ' (row-major), by nested Gauss-Legendre quadrature\n'+'\n'.join(masses)+'\n',
                    'brep-classes.tsv': '# valid case\tits class (D9: solid, sheet, wire, acorn or general)\n'
                    + '\n'.join(classes)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    models, files = generate()
    rows = files['brep-expected.tsv'].splitlines()[1:]
    valid = sum(1 for row in rows if row.endswith('\t'))
    print(f'{len(rows)} cases ({len(models)} with OCCT rows), {valid} valid')
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)


if __name__ == '__main__':
    main()
