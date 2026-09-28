#!/usr/bin/env python3
"""Fixtures for S8a of REVIEW_NOTES.md: prisms of line and arc profiles
split by a plane.

`split-cases.txt` lists each case in the protocol of `identity-cases.txt`
(a prism: `case`, `op`, `frame`, `offsets`, `boundary` rows) with a `split ox
oy oz nx ny nz` row (the plane through a point with a normal) before `end`;
`split-expected.tsv` gives each case's sides from `split_reference.py`
(`side below|above|whole volume area cx cy cz`: the totals of the pieces on
each side), and `split-frames.tsv` each
prism's stored frame axes (`stored_axes`, the kernel's `Frame3::new`; its tests
check them bit for bit). Every class has a case: planes crossing the caps and
the walls at angles, normal to the axis (M3's height split), parallel to it
across line and arc walls, through two vertical edges, through one vertex
only, tangent to an arc wall along a ruling, lying in a cap, missing, and
crossing a hole, a non-convex profile in two pieces of section, and a
circle's whole wall in an ellipse. No Rust result supplies an expectation.

S8c: `split-primitive-cases.txt` gives one cone, frustum, sphere or zone per
line (`cone NAME tol ox oy oz nx ny nz xx xy xz r1 r2 h split ...` or
`sphere NAME tol ox oy oz nx ny nz xx xy xz R a1 a2 split ...`, as the
kernel's builders and `BRepPrimAPI_MakeCone`/`MakeSphere` take them) and
`split-primitive-expected.tsv` its sides (`split_reference.revolved_rows`):
planes normal to the axis (through a frustum, an apex cone, a zone, a whole
sphere, at an apex, in a cap, missing, touching a sphere), a whole sphere
by oblique planes, planes containing the axis, and the conic and circle
sections that wait for procedural edges (S8d). S8d: `split-torus-cases.txt`
the same with `torus NAME tol ... major minor angle` lines (whole tori, the
angle a full turn) and `split-torus-expected.tsv` (`split_reference.
torus_rows`): planes normal to the axis (through the equator, above it,
touching, missing, in a tilted frame), containing it, parallel to it
through the tube, and oblique.
"""
import argparse
import math
from pathlib import Path
import struct

from identity_reference import Boundary, Case, encode_case
from curve_surface_reference import stored_axes
import split_reference as ref

ROOT = Path(__file__).resolve().parents[1]
XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (1.0, -2.0, 0.5, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)


def arc(cx, cy, r, ccw=True):
    return (cx, cy, r, ccw)


def path(points, segments):
    return Boundary(points=points, segments=segments)


def prism(name, boundaries, frame=XY, start=0.0, end=5.0):
    return Case(name, 1e-7, 81, frame, start, end, boundaries)


def cases():
    square = Boundary(points=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    hole = Boundary(points=[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)])
    ushape = Boundary(points=[(0.0, 0.0), (9.0, 0.0), (9.0, 6.0), (6.0, 6.0), (6.0, 2.0), (3.0, 2.0),
                              (3.0, 6.0), (0.0, 6.0)])
    stadium = path([(0.0, -1.0), (3.0, -1.0), (3.0, 1.0), (0.0, 1.0)],
                   [None, arc(3.0, 0.0, 1.0), None, arc(0.0, 0.0, 1.0)])
    rounded = path([(1.0, 0.0), (7.0, 0.0), (8.0, 1.0), (8.0, 5.0), (7.0, 6.0), (1.0, 6.0), (0.0, 5.0), (0.0, 1.0)],
                   [None, arc(7.0, 1.0, 1.0), None, arc(7.0, 5.0, 1.0), None, arc(1.0, 5.0, 1.0), None,
                    arc(1.0, 1.0, 1.0)])
    disc = Boundary(circle=(0.0, 0.0, 2.0))
    slot = path([(-2.0, -1.0), (2.0, -1.0), (2.0, 1.0), (-2.0, 1.0)],
                [None, arc(2.0, 0.0, 1.0), None, arc(-2.0, 0.0, 1.0)])
    frame10 = Boundary(points=[(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)])
    return [
        (prism('sq_height', [square]), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        (prism('sq_tilted', [square]), (5.0, 5.0, 2.5, 1.0, 0.0, 2.0)),
        (prism('sq_skew', [square]), (4.0, 6.0, 2.0, 1.0, 2.0, 3.0)),
        (prism('sq_vertical', [square]), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('sq_diagonal_edges', [square]), (5.0, 5.0, 0.0, 1.0, 1.0, 0.0)),
        (prism('sq_through_corner', [square]), (10.0, 10.0, 5.0, 1.0, 1.0, 1.0)),
        (prism('sq_in_cap', [square]), (3.0, 4.0, 5.0, 0.0, 0.0, 1.0)),
        (prism('sq_corner_cut', [square]), (0.0, 0.0, 7.0, 1.0, 1.0, 4.0)),
        (prism('sq_miss', [square]), (0.0, 0.0, 7.0, 0.25, 0.25, 4.0)),
        (prism('sq_hole_crossed', [square, hole]), (5.0, 5.0, 2.5, 1.0, 0.5, 1.5)),
        (prism('sq_hole_vertical', [square, hole]), (5.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('u_two_sections', [ushape]), (0.0, 4.0, 0.0, 0.0, 1.0, 0.0)),
        (prism('u_tilted', [ushape], start=0.0, end=3.0), (4.5, 3.0, 1.5, 0.25, 1.0, 0.5)),
        (prism('stadium_height', [stadium], end=2.0), (0.0, 0.0, 0.75, 0.0, 0.0, 1.0)),
        (prism('stadium_tilted', [stadium], end=2.0), (1.5, 0.0, 1.0, 1.0, 0.0, 2.0)),
        (prism('stadium_arc_only', [stadium], end=2.0), (3.5, 0.0, 1.0, 1.0, 0.0, 0.25)),
        (prism('stadium_vertical_arc', [stadium], end=2.0), (3.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('stadium_tangent', [stadium], end=2.0), (4.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('rounded_skew', [rounded], end=3.0), (4.0, 3.0, 1.5, 1.0, 1.0, 2.0)),
        (prism('disc_wall_ellipse', [disc], end=4.0), (0.0, 0.0, 2.0, 1.0, 0.0, 4.0)),
        (prism('disc_through_caps', [disc], end=4.0), (0.0, 0.0, 2.0, 1.0, 0.0, 1.0)),
        (prism('disc_vertical', [disc], end=4.0), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (prism('slot_hole_tilted', [frame10, slot], frame=TILT, start=0.0, end=3.0),
         (1.0, -2.0, 2.0, 0.25, 0.5, 1.0)),
        (prism('tilted_frame_in_cap', [square], frame=TILT, start=0.0, end=4.0),
         (1.0, -2.0, 0.5, 0.0, 3.0, 4.0)),
        (prism('tilted_frame_height', [square], frame=TILT, start=0.0, end=4.0),
         (1.0, -0.8, 2.1, 0.0, 3.0, 4.0)),
        (prism('reversed_offsets', [stadium], start=2.0, end=-1.0), (1.0, 0.0, 0.5, 1.0, 1.0, 2.0)),
    ]


HALF_PI = math.pi/2


def primitive_cases():
    """(kind, name, frame, params, plane)."""
    cone, sphere = 'cone', 'sphere'
    return [
        (cone, 'frustum_height', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 1.25, 0.0, 0.0, 1.0)),
        (cone, 'apex_height', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.5, 0.0, 0.0, 1.0)),
        (cone, 'apex_up_height', XY, (0.0, 1.5, 2.0), (0.0, 0.0, 0.5, 0.0, 0.0, -1.0)),
        (cone, 'widening_height', XY, (1.0, 2.5, 0.5), (0.0, 0.0, 0.375, 0.0, 0.0, 1.0)),
        (cone, 'apex_touch', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 3.0, 0.0, 0.0, 1.0)),
        (cone, 'frustum_in_cap', XY, (2.0, 1.0, 3.0), (0.5, 0.0, 0.0, 0.0, 0.0, 1.0)),
        (cone, 'frustum_miss', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 4.0, 0.0, 0.0, 1.0)),
        (cone, 'tilted_height', TILT, (2.0, 0.5, 3.0), (1.0, -0.8, 2.1, 0.0, 3.0, 4.0)),
        (cone, 'frustum_meridian', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (cone, 'apex_meridian', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, 1.0, 1.0, 0.0)),
        (cone, 'apex_oblique', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, 1.0, 0.0, 2.0)),
        (cone, 'frustum_parallel', XY, (2.0, 1.0, 3.0), (0.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (sphere, 'sphere_equator', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_height', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 0.75, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_oblique', XY, (2.0, -HALF_PI, HALF_PI), (0.3, 0.2, 0.1, 1.0, 2.0, 2.0)),
        (sphere, 'sphere_tilted', TILT, (1.5, -HALF_PI, HALF_PI), (1.25, -2.0, 0.5, 1.0, -1.0, 0.5)),
        (sphere, 'sphere_touch', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        (sphere, 'sphere_miss', XY, (2.0, -HALF_PI, HALF_PI), (0.0, 3.0, 0.0, 0.0, 1.0, 0.0)),
        (sphere, 'zone_height', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.25, 0.0, 0.0, 1.0)),
        (sphere, 'cap_height', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 1.5, 0.0, 0.0, 1.0)),
        (sphere, 'zone_meridian', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.0, 0.0, 1.0, 0.0)),
        (sphere, 'cap_meridian', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 0.0, 1.0, 1.0, 0.0)),
        (sphere, 'zone_oblique', XY, (2.0, -0.5, 1.0), (0.0, 0.0, 0.5, 1.0, 0.0, 1.0)),
    ]


def conic_cases():
    """(kind, name, frame, params, plane): cones, frusta, zones and caps by
    planes neither normal to their axes nor containing them (S8d.2), beside
    S8c's apex_oblique, frustum_parallel and zone_oblique: a closed ellipse,
    a parabola, hyperbolas, a tongue at the top, the rulings through a
    frustum's virtual apex, an ellipse touching both rims, a zone's side cap
    and circle arcs across both rims."""
    cone, sphere = 'cone', 'sphere'
    return [
        (cone, 'apex_ellipse', XY, (2.0, 0.0, 3.0), (0.0, 0.0, 1.0, -0.2, 0.0, 1.0)),
        (cone, 'apex_parabola', XY, (2.0, 0.0, 2.0), (0.5, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (cone, 'apex_hyperbola', XY, (2.0, 0.0, 3.0), (0.8, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (cone, 'frustum_top_tongue', XY, (2.0, 1.0, 3.0), (0.6, 0.0, 3.0, 1.0, 0.0, 1.0)),
        (cone, 'frustum_rulings', XY, (2.0, 1.0, 3.0), (0.0, 0.0, 6.0, 1.0, 0.0, 0.1)),
        (cone, 'frustum_touch_both', XY, (2.0, 1.0, 3.0), (2.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
        (cone, 'widening_oblique', XY, (1.0, 2.5, 0.5), (0.0, 0.0, 0.25, 0.3, 0.2, 1.0)),
        (cone, 'tilted_oblique', TILT, (2.0, 0.5, 3.0), (1.0, -0.8, 2.1, 1.0, 3.0, 4.0)),
        (sphere, 'zone_side_cap', XY, (2.0, -0.5, 1.0), (1.8, 0.0, 0.0, 1.0, 0.0, 0.0)),
        (sphere, 'cap_oblique', XY, (2.0, 0.25, HALF_PI), (0.0, 0.0, 1.2, 1.0, 0.0, 1.0)),
        (sphere, 'zone_both_rims', XY, (2.0, -0.5, 1.0), (0.3, 0.0, 0.0, 1.0, 0.0, 0.2)),
        (sphere, 'tilted_zone_oblique', TILT, (1.5, -0.7, 0.9), (1.25, -2.0, 0.5, 1.0, -1.0, 0.5)),
    ]


def torus_cases():
    """(name, frame, (major, minor), plane): whole tori (S8d)."""
    return [
        ('torus_equator', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
        ('torus_height', XY, (3.0, 1.0), (0.0, 0.0, 0.5, 0.0, 0.0, 1.0)),
        ('torus_high', XY, (3.0, 1.0), (0.0, 0.0, 0.875, 0.0, 0.0, -1.0)),
        ('torus_touch', XY, (3.0, 1.0), (0.0, 0.0, 1.0, 0.0, 0.0, 1.0)),
        ('torus_miss', XY, (3.0, 1.0), (0.0, 0.0, 2.0, 0.0, 0.0, 1.0)),
        ('torus_meridian', XY, (3.0, 1.0), (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_meridian_skew', XY, (3.0, 1.0), (0.0, 0.0, 0.5, 1.0, 1.0, 0.0)),
        ('torus_tilted_height', TILT, (2.5, 0.75), (1.0, -1.82, 0.74, 0.0, 3.0, 4.0)),
        ('torus_parallel_outer', XY, (3.0, 1.0), (3.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_parallel_inner', XY, (3.0, 1.0), (2.5, 0.0, 0.0, 1.0, 0.0, 0.0)),
        ('torus_oblique', XY, (3.0, 1.0), (0.0, 0.0, 0.25, 0.25, 0.0, 1.0)),
        ('torus_oblique_tube', XY, (3.0, 1.0), (3.0, 0.0, 0.0, 1.0, 0.0, 1.0)),
    ]


def primitive_line(kind, name, frame, params, plane):
    words = [kind, name, '1e-07'] + [repr(float(v)) for v in frame+params] + ['split'] + \
        [repr(float(v)) for v in plane]
    return ' '.join(words)


def encode(case, plane):
    text = encode_case(case)
    body, end = text.rsplit('\nend', 1)
    return body+'\nsplit '+' '.join(repr(float(v)) for v in plane)+'\nend'


def generate():
    blocks, out = [], ['# case\trow (S8a, split_reference.py)']
    frames = ['# case\taxis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case, plane in cases():
        blocks.append(encode(case, plane))
        for row in ref.rows(case, plane):
            out.append(f'{case.name}\t{ref.text(row)}')
        _, x, y, n = stored_axes(case.frame)
        for key, v in (('n', n), ('x', x), ('y', y)):
            frames.append(f'{case.name}\t{key}\t'+' '.join(struct.pack('>d', c).hex() for c in v))
    lines, prim = [], ['# case\trow (S8c, split_reference.revolved_rows)']
    for kind, name, frame, params, plane in primitive_cases():
        lines.append(primitive_line(kind, name, frame, params, plane))
        for row in ref.revolved_rows(kind, frame, params, plane):
            prim.append(f'{name}\t{ref.text(row)}')
    clines, conics = [], ['# case\trow (S8d.2, split_reference.revolved_rows)']
    for kind, name, frame, params, plane in conic_cases():
        clines.append(primitive_line(kind, name, frame, params, plane))
        for row in ref.revolved_rows(kind, frame, params, plane):
            conics.append(f'{name}\t{ref.text(row)}')
    tlines, tori = [], ['# case\trow (S8d, split_reference.torus_rows)']
    for name, frame, params, plane in torus_cases():
        tlines.append(primitive_line('torus', name, frame, params+(2*math.pi,), plane))
        for row in ref.torus_rows(frame, params, plane):
            tori.append(f'{name}\t{ref.text(row)}')
    return {'split-cases.txt': '\n'.join(blocks)+'\n', 'split-expected.tsv': '\n'.join(out)+'\n',
            'split-frames.tsv': '\n'.join(frames)+'\n',
            'split-primitive-cases.txt': '\n'.join(lines)+'\n',
            'split-primitive-expected.tsv': '\n'.join(prim)+'\n',
            'split-conic-cases.txt': '\n'.join(clines)+'\n',
            'split-conic-expected.tsv': '\n'.join(conics)+'\n',
            'split-torus-cases.txt': '\n'.join(tlines)+'\n',
            'split-torus-expected.tsv': '\n'.join(tori)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = generate()
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                raise SystemExit(f'{path} is stale')
        else:
            path.write_text(contents)
    print(len(cases()), 'cases,', len(primitive_cases()), 'primitive cases,', len(conic_cases()),
          'conic cases,', len(torus_cases()), 'tori')


if __name__ == '__main__':
    main()
