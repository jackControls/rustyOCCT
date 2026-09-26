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

from cell_reference import declare, encode as encode_cell, gap_bounds, to_cell, validate as validate_cell

from brep_reference import (Arc2, Arc3, Cylinder, Edge, Face, Frame, Line2, Line3, Model,
                            Plane, TAU, Use, atan2_rn, cos_rn, encode, hypot_rn, number, sin_rn,
                            validate)

ROOT = Path(__file__).resolve().parents[1]
Z = (0.0, 0.0, 1.0)
X = (1.0, 0.0, 0.0)


def angle_of(p, c):
    return atan2_rn(p[1]-c[1], p[0]-c[0])


def prism(name, boundaries, z0=0.0, z1=1.0, tolerance=1e-7):
    """Boundaries: lists of pieces ('line', p) or ('arc', p, center, ccw) or
    ('circle', center, radius, ccw). Outer boundaries run counter-clockwise
    about +z and holes clockwise; material is always on the left."""
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
                if piece[0] == 'line':
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
    return out


def generate():
    bases = base_cases()
    models = list(bases.values())+mutations(bases)
    cells = [declare(c) for c in [to_cell(m) for m in models]+cell_cases(bases)]
    names = [c.name for c in cells]
    assert len(names) == len(set(names)), 'duplicate case names'
    text = '\n'.join(encode_cell(c) for c in cells)+'\n'
    rows, lows = [], []
    for c in cells:
        issues = validate_cell(c)
        rows.append(c.name+'\t'+';'.join(f'{k}:{e}' for k, e in issues))
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
                    '\tcertain lower value of its gap\n'+'\n'.join(lows)+'\n'}


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
