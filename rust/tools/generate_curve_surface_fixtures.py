#!/usr/bin/env python3
"""Fixtures for S7c.1 of REVIEW_NOTES.md: lines and circles against planes,
cylinders, cones, spheres and tori.

`curve-surface-cases.txt` lists each case as `case NAME`, a `curve line x0
y0 z0 x1 y1 z1` or `curve circle ox oy oz nx ny nz hx hy hz radius` row (the
frame as `Frame3::new` takes it), a `surface` row in the protocol of
`analytic-intersection-cases.txt`, and `end`; `curve-surface-expected.tsv`
gives each case's rows from `curve_surface_reference.py`, and
`curve-surface-frames.tsv` the stored unit normals (every normal has a zero
coordinate and a Pythagorean pair). Every class has an exact case: a
crossing, a tangency (lines touching a cylinder, a sphere, a torus's top and
inner equator; circles touching a plane, a sphere, a cylinder, a torus's
outer equator), containment (a line in a plane or along a cylinder, a circle
in a plane, on a sphere, on a cylinder, a torus's meridian and parallel),
and a miss. No Rust result supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

import analytic_intersection_reference as ana
import curve_surface_reference as ref
from generate_analytic_intersection_fixtures import X, Y, Z, cone, cylinder, encode, plane, sphere, torus
from identity_reference import frame_axes

ROOT = Path(__file__).resolve().parents[1]
TILT = (0.0, 3.0, 4.0)


def line(p0, p1):
    return ('line', (*p0, *p1))


def circle(o, n, r, x=None):
    from generate_analytic_intersection_fixtures import hint
    return ('circle', (*o, *n, *(x or hint(n)), r))


def cases():
    o = (0.0, 0.0, 0.0)
    thin = torus(o, Z, 2.0, 0.5)
    return [
        # Lines.
        ('l_plane_cross', line((0.0, 0.0, -1.0), (1.0, 1.0, 1.0)), plane(o, Z)),
        ('l_plane_parallel', line((0.0, 0.0, 1.0), (1.0, 2.0, 1.0)), plane(o, Z)),
        ('l_plane_contained', line((1.0, 0.0, 0.0), (0.0, 3.0, 0.0)), plane(o, Z)),
        ('l_cylinder_two', line((-2.0, 0.25, 0.5), (2.0, 0.25, 0.5)), cylinder(o, Z, 1.0)),
        ('l_cylinder_tangent', line((-2.0, 1.0, 0.5), (2.0, 1.0, 1.5)), cylinder(o, Z, 1.0)),
        ('l_cylinder_ruling', line((1.0, 0.0, 0.0), (1.0, 0.0, 1.0)), cylinder(o, Z, 1.0)),
        ('l_cylinder_miss', line((-2.0, 3.0, 0.0), (2.0, 3.0, 1.0)), cylinder(o, Z, 1.0)),
        ('l_sphere_two', line((-3.0, 0.5, 0.25), (3.0, 0.5, 0.25)), sphere(o, 2.0)),
        ('l_sphere_tangent', line((-3.0, 2.0, 0.0), (3.0, 2.0, 0.0)), sphere(o, 2.0)),
        ('l_sphere_miss', line((-3.0, 3.0, 0.0), (3.0, 3.0, 1.0)), sphere(o, 2.0)),
        ('l_cone_two', line((-3.0, 0.25, 1.0), (3.0, 0.25, 1.0)), cone(o, Z, 0.0, 0.5)),
        ('l_cone_nappes', line((0.25, 0.0, -2.0), (0.25, 0.5, 2.0)), cone(o, Z, 0.0, 0.5)),
        ('l_cone_narrow', line((-3.0, 0.0, 0.25), (3.0, 0.0, 0.5)), cone(o, Z, 0.0, 0.05)),
        ('l_cone_miss', line((-3.0, 1.0, 0.25), (3.0, 1.0, 0.5)), cone(o, Z, 0.0, 0.05)),
        ('l_torus_four', line((-3.0, 0.0, 0.0), (3.0, 0.0, 0.0)), thin),
        ('l_torus_top', line((-3.0, 0.0, 0.5), (3.0, 0.0, 0.5)), thin),
        ('l_torus_inner', line((1.5, 0.0, -1.0), (1.5, 0.0, 1.0)), thin),
        ('l_torus_axis', line((0.0, 0.0, -1.0), (0.0, 0.0, 1.0)), thin),
        ('l_torus_tilted', line((-3.0, 0.5, 0.25), (3.0, 1.0, -0.25)), torus((0.0, 0.5, 0.0), TILT, 2.0, 0.5)),
        # Circles.
        ('c_plane_cross', circle(o, Z, 1.0), plane((0.5, 0.0, 0.0), X)),
        ('c_plane_tangent', circle(o, Z, 1.0), plane((1.0, 0.0, 0.0), X)),
        ('c_plane_contained', circle(o, Z, 1.0), plane((5.0, 7.0, 0.0), Z)),
        ('c_plane_parallel', circle(o, Z, 1.0), plane((0.0, 0.0, 1.0), Z)),
        ('c_plane_tilted', circle((0.0, 0.5, 0.0), TILT, 1.5), plane((0.25, 0.0, 0.0), X)),
        ('c_sphere_two', circle(o, Z, 1.0), sphere((1.0, 0.0, 0.0), 1.0)),
        ('c_sphere_tangent', circle(o, Z, 1.0), sphere((2.0, 0.0, 0.0), 1.0)),
        ('c_sphere_contained', circle((0.0, 0.0, 3.0), Z, 4.0), sphere(o, 5.0)),
        ('c_cylinder_four', circle(o, Y, 1.0), cylinder(o, Z, 0.5)),
        ('c_cylinder_tangent', circle((1.5, 0.0, 0.0), Z, 1.0), cylinder(o, Z, 0.5)),
        ('c_cylinder_contained', circle((0.0, 0.0, 2.0), Z, 0.5), cylinder(o, Z, 0.5)),
        ('c_cone_two', circle((1.0, 0.0, 1.0), X, 1.0), cone(o, Z, 0.0, 0.5)),
        ('c_torus_meridian', circle((2.0, 0.0, 0.0), Y, 0.5), thin),
        ('c_torus_parallel', circle((0.0, 0.0, 0.5), Z, 2.0), thin),
        ('c_torus_cross', circle(o, X, 2.0), thin),
        ('c_torus_tangent', circle((0.5, 0.0, 0.0), Z, 3.0), thin),
        ('c_torus_tilted', circle((0.0, 0.5, 0.25), TILT, 1.75), thin),
    ]


def encode(name, curve, surface):
    kind, values = curve
    rows = [f'case {name}', f'curve {kind} '+' '.join(repr(float(v)) for v in values)]
    skind, svalues = surface
    rows.append(f'surface {skind} '+' '.join(repr(float(v)) for v in svalues))
    return '\n'.join(rows)


def rows(curve, surface_block):
    kind, values = curve
    _, (surface,) = ana.parse(surface_block)
    if kind == 'line':
        p0, p1 = values[:3], values[3:]
        from fractions import Fraction as F
        return ref.line_rows(tuple(F(v) for v in p0), tuple(F(v) for v in p1), surface)
    return ref.circle_rows(tuple(values[:9]), values[9], surface)


def generate():
    blocks, out = [], ['# case\trow (S7c.1, curve_surface_reference.py)']
    frames = ['# case\tsurface\tstored unit normal (the reference\'s Frame3::new, as hex bits)']
    for name, curve, surface in cases():
        block = encode(name, curve, surface)
        blocks.append(block+'\nend')
        surface_block = f'case {name}\n'+block.splitlines()[2]
        for row in rows(curve, surface_block):
            out.append(f'{name}\t{ref.text(row)}')
        framed = [surface[1]] if curve[0] == 'line' else [curve[1][:9], surface[1]]
        for k, values in enumerate(framed):
            _, _, _, n = frame_axes(tuple(values[:9]))
            frames.append(f'{name}\t{k}\t'+' '.join(struct.pack('>d', v).hex() for v in n))
    return {'curve-surface-cases.txt': '\n'.join(blocks)+'\n',
            'curve-surface-expected.tsv': '\n'.join(out)+'\n',
            'curve-surface-frames.tsv': '\n'.join(frames)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = generate()
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)
    print(len(cases()), 'cases')


if __name__ == '__main__':
    main()
