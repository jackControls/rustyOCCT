#!/usr/bin/env python3
"""Cone, frustum, sphere and zone fixtures from the independent primitive
reference (S3).

primitive-cases.txt lists one cone per line (`cone NAME ox oy oz nx ny nz
xx xy xz r1 r2 h`, the arguments of BRepPrimAPI_MakeCone with a gp_Ax2);
sphere-cases.txt one sphere per line (`sphere NAME ox oy oz nx ny nz xx xy
xz R a1 a2`, BRepPrimAPI_MakeSphere's, latitudes in radians); torus-cases.txt
one torus per line (`torus NAME ox oy oz nx ny nz xx xy xz R r a1 a2 angle`,
BRepPrimAPI_MakeTorus's, radians). primitive-expected.tsv, sphere-expected.tsv
and torus-expected.tsv hold, per case, OCCT's
distinct subshape counts, the volume, surface area, centre of mass and the
six entries of the symmetric inertia matrix about it, and every face, edge
and vertex as `primitive_reference.py` derives them. No kernel or OCCT
result supplies an expectation.
"""
import argparse
import math
from pathlib import Path

import mpmath as mp

from primitive_reference import (HALF_PI, TWO_PI, Cone, Sphere, Torus, counts, edges, faces, mass,
                                 sphere_counts, sphere_edges, sphere_faces, sphere_mass, sphere_vertices,
                                 torus_counts, torus_edges, torus_faces, torus_mass, torus_vertices,
                                 vertices)

ROOT = Path(__file__).resolve().parents[1]


def number(x):
    return repr(float(x))


def corpus():
    cases = []
    z, x = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
    cases.append(Cone('apex', (0.0, 0.0, 0.0), z, x, 2.0, 0.0, 3.0))
    cases.append(Cone('apex_at_base', (0.0, 0.0, 0.0), z, x, 0.0, 1.5, 2.0))
    cases.append(Cone('frustum', (1.0, 2.0, 3.0), z, x, 2.0, 1.0, 3.0))
    cases.append(Cone('frustum_widening', (0.0, 0.0, 0.0), z, x, 1.0, 2.5, 0.5))
    cases.append(Cone('nearly_cylinder', (0.0, 0.0, 0.0), z, x, 1.0, 1.0009765625, 4.0))
    cases.append(Cone('flat', (0.0, 0.0, 0.0), z, x, 5.0, 4.0, 0.125))
    cases.append(Cone('needle', (0.0, 0.0, 0.0), z, x, 0.0625, 0.0, 8.0))
    cases.append(Cone('millimetre', (0.0, 0.0, 0.0), z, x, 0.002, 0.001, 0.003))
    cases.append(Cone('large', (0.0, 0.0, 0.0), z, x, 400.0, 150.0, 900.0))
    cases.append(Cone('far', (10000.0, -20000.0, 5000.0), z, x, 2.0, 0.5, 3.0))
    # Rotated frames from a fixed sequence of directions and x hints.
    for k in range(12):
        a, b = 0.37*k+0.1, 0.61*k+0.3
        normal = (round(math.cos(a)*math.sin(b), 6), round(math.sin(a)*math.sin(b), 6), round(math.cos(b), 6))
        hint = (round(math.cos(1.3*k+0.2), 6), round(math.sin(1.3*k+0.2), 6), 0.25)
        if abs(normal[0]*hint[0]+normal[1]*hint[1]+normal[2]*hint[2]) > 0.9:
            hint = (0.0, 0.0, 1.0)
        r1 = [1.0, 0.0, 2.0, 0.75][k % 4]
        r2 = [0.0, 1.25, 0.5, 3.0][k % 4]
        origin = (0.5*k-3.0, 1.5-0.25*k, 0.125*k)
        cases.append(Cone(f'rotated_{k}', origin, normal, hint, r1, r2, 1.0+0.5*k))
    return cases


def spheres():
    cases = []
    z, x = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
    o = (0.0, 0.0, 0.0)
    cases.append(Sphere('full', o, z, x, 2.0, -HALF_PI, HALF_PI))
    cases.append(Sphere('upper_hemisphere', o, z, x, 1.5, 0.0, HALF_PI))
    cases.append(Sphere('lower_hemisphere', (1.0, 2.0, 3.0), z, x, 1.5, -HALF_PI, 0.0))
    cases.append(Sphere('zone', o, z, x, 3.0, -0.5, 0.7))
    cases.append(Sphere('zone_above', o, z, x, 1.0, 0.25, 0.5))
    cases.append(Sphere('thin_equator', o, z, x, 1.0, -0.0009765625, 0.0009765625))
    cases.append(Sphere('north_cap', o, z, x, 2.5, 1.25, HALF_PI))
    cases.append(Sphere('south_cap', o, z, x, 2.5, -HALF_PI, -1.3))
    cases.append(Sphere('millimetre', o, z, x, 0.002, -HALF_PI, HALF_PI))
    cases.append(Sphere('large', o, z, x, 400.0, -0.25, HALF_PI))
    cases.append(Sphere('far', (10000.0, -20000.0, 5000.0), z, x, 2.0, -HALF_PI, HALF_PI))
    # Rotated frames from a fixed sequence of directions and x hints.
    latitudes = [(-HALF_PI, HALF_PI), (0.0, HALF_PI), (-0.75, 0.5), (-HALF_PI, 0.3)]
    for k in range(11):
        a, b = 0.41*k+0.2, 0.53*k+0.4
        normal = (round(math.cos(a)*math.sin(b), 6), round(math.sin(a)*math.sin(b), 6), round(math.cos(b), 6))
        hint = (round(math.cos(1.1*k+0.3), 6), round(math.sin(1.1*k+0.3), 6), 0.25)
        if abs(normal[0]*hint[0]+normal[1]*hint[1]+normal[2]*hint[2]) > 0.9:
            hint = (0.0, 0.0, 1.0)
        a1, a2 = latitudes[k % 4]
        origin = (0.5*k-3.0, 1.5-0.25*k, 0.125*k)
        cases.append(Sphere(f'rotated_{k}', origin, normal, hint, 0.5+0.25*k, a1, a2))
    return cases


def tori():
    cases = []
    z, x = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
    o = (0.0, 0.0, 0.0)
    cases.append(Torus('whole', o, z, x, 3.0, 1.0, 0.0, TWO_PI, TWO_PI))
    cases.append(Torus('whole_thin', (1.0, 2.0, 3.0), z, x, 10.0, 0.5, 0.0, TWO_PI, TWO_PI))
    cases.append(Torus('whole_fat', o, z, x, 2.0, 1.75, 0.0, TWO_PI, TWO_PI))
    cases.append(Torus('millimetre', o, z, x, 0.003, 0.001, 0.0, TWO_PI, TWO_PI))
    cases.append(Torus('far', (10000.0, -20000.0, 5000.0), z, x, 3.0, 1.0, 0.0, TWO_PI, TWO_PI))
    cases.append(Torus('outer_half', o, z, x, 3.0, 1.0, -HALF_PI, HALF_PI, TWO_PI))
    cases.append(Torus('inner_half', o, z, x, 3.0, 1.0, HALF_PI, 3*HALF_PI, TWO_PI))
    cases.append(Torus('segment', o, z, x, 4.0, 1.5, 0.3, 2.0, TWO_PI))
    cases.append(Torus('segment_below', o, z, x, 4.0, 1.5, -1.0, 0.5, TWO_PI))
    cases.append(Torus('quarter_wedge', o, z, x, 3.0, 1.0, 0.0, TWO_PI, HALF_PI))
    cases.append(Torus('half_wedge', o, z, x, 3.0, 1.0, 0.0, TWO_PI, 2*HALF_PI))
    cases.append(Torus('wide_wedge', (1.0, 2.0, 3.0), z, x, 5.0, 2.0, 0.0, TWO_PI, 5.0))
    # Rotated frames from a fixed sequence of directions and x hints.
    shapes = [(0.0, TWO_PI, TWO_PI), (-HALF_PI, HALF_PI, TWO_PI), (0.0, TWO_PI, 1.25), (0.4, 2.5, TWO_PI),
              (0.0, TWO_PI, 4.0)]
    for k in range(10):
        a, b = 0.29*k+0.15, 0.47*k+0.35
        normal = (round(math.cos(a)*math.sin(b), 6), round(math.sin(a)*math.sin(b), 6), round(math.cos(b), 6))
        hint = (round(math.cos(0.9*k+0.1), 6), round(math.sin(0.9*k+0.1), 6), 0.25)
        if abs(normal[0]*hint[0]+normal[1]*hint[1]+normal[2]*hint[2]) > 0.9:
            hint = (0.0, 0.0, 1.0)
        a1, a2, angle = shapes[k % 5]
        origin = (0.5*k-3.0, 1.5-0.25*k, 0.125*k)
        cases.append(Torus(f'rotated_{k}', origin, normal, hint, 2.0+0.25*k, 0.5+0.1*k, a1, a2, angle))
    return cases


def encode(c):
    if isinstance(c, Torus):
        return ' '.join(['torus', c.name] + [number(v) for v in (*c.origin, *c.normal, *c.x, c.major, c.minor,
                                                                 c.a1, c.a2, c.angle)])
    if isinstance(c, Sphere):
        return ' '.join(['sphere', c.name] + [number(v) for v in (*c.origin, *c.normal, *c.x, c.radius, c.a1, c.a2)])
    return ' '.join(['cone', c.name] + [number(v) for v in (*c.origin, *c.normal, *c.x, c.r1, c.r2, c.height)])


def expected(c):
    if isinstance(c, Torus):
        return shape_rows(c, torus_counts(c), torus_mass(c), torus_faces(c), torus_edges(c), torus_vertices(c))
    if isinstance(c, Sphere):
        return shape_rows(c, sphere_counts(c), sphere_mass(c), sphere_faces(c), sphere_edges(c), sphere_vertices(c))
    return shape_rows(c, counts(c), mass(c), faces(c), edges(c), vertices(c))


def shape_rows(c, counts, mass, faces, edges, vertices):
    rows = [f'{c.name}\tcounts\t' + ' '.join(map(str, counts))]
    volume, area, centre, inertia = mass
    entries = [inertia[i][j] for i in range(3) for j in range(i, 3)]
    rows.append(f'{c.name}\tprops\t' + ' '.join(number(v) for v in [volume, area, *centre, *entries]))
    for kind, a, centre in sorted(faces, key=lambda f: (f[0], float(f[1]))):
        rows.append(f'{c.name}\tface\t{kind} {number(a)} ' + ' '.join(number(v) for v in centre))
    for degenerate, closed, length, mid in sorted(edges, key=lambda e: (e[0], e[1], float(e[2]))):
        rows.append(f'{c.name}\tedge\t{degenerate} {closed} {number(length)} ' + ' '.join(number(v) for v in mid))
    for p in sorted(vertices, key=lambda p: [float(v) for v in p]):
        rows.append(f'{c.name}\tvertex\t' + ' '.join(number(v) for v in p))
    return rows


def generate():
    out = {}
    for stem, cases in (('primitive', corpus()), ('sphere', spheres()), ('torus', tori())):
        out[f'{stem}-cases.txt'] = '\n'.join(encode(c) for c in cases)+'\n'
        out[f'{stem}-expected.tsv'] = ('# case\tkind\tvalues\n'
                                      + '\n'.join(row for c in cases for row in expected(c))+'\n')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for name, text in generate().items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != text:
                raise SystemExit(f'{name} is stale')
        else:
            path.write_text(text)
    print(f'{len(corpus())} cones, {len(spheres())} spheres, {len(tori())} tori')


if __name__ == '__main__':
    main()
