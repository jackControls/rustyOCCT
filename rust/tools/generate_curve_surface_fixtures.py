#!/usr/bin/env python3
"""Fixtures for S7c of REVIEW_NOTES.md: lines, circles (S7c.1), ellipses,
hyperbolas and splines (S7c.2) against planes, cylinders, cones, spheres and
tori.

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
and a miss. S7c.2 adds ellipses and hyperbolas (`curve ellipse|hyperbola ox
oy oz nx ny nz hx hy hz major minor`: crossings, tangencies to a plane, a
sphere, a cylinder and a torus, containment in a plane, the hyperbola's
other branch, misses) and rational B-splines against tori and cones (`curve
spline degree n x y z w ... k knot mult ...`: crossings, a line's two
tangencies to a torus's top, an exact rational quarter circle on a torus's
parallel and on a cone's reference circle, a partial overlap, a rational
apex, misses). No Rust result supplies an expectation.
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


def conic(kind, o, n, major, minor, x=None):
    from generate_analytic_intersection_fixtures import hint
    return (kind, (*o, *n, *(x or hint(n)), major, minor))


def spline(degree, poles, knots, mults):
    """A clamped rational B-spline: poles `(x, y, z, w)`."""
    return ('spline', (degree, tuple(poles), tuple(knots), tuple(mults)))


def cases():
    o = (0.0, 0.0, 0.0)
    thin = torus(o, Z, 2.0, 0.5)
    quarter = [(2.0, 0.0, 0.5, 1.0), (2.0, 2.0, 0.5, 1.0), (0.0, 2.0, 0.5, 2.0)]
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
        # Ellipses (S7c.2).
        ('e_plane_cross', conic('ellipse', o, Z, 2.0, 1.0), plane((1.0, 0.0, 0.0), X)),
        ('e_plane_tangent', conic('ellipse', o, Z, 2.0, 1.0), plane((2.0, 0.0, 0.0), X)),
        ('e_plane_contained', conic('ellipse', o, Z, 2.0, 1.0), plane((3.0, 5.0, 0.0), Z)),
        ('e_plane_parallel', conic('ellipse', o, Z, 2.0, 1.0), plane((0.0, 0.0, 1.0), Z)),
        ('e_plane_tilted', conic('ellipse', (0.0, 0.5, 0.0), TILT, 1.5, 1.0), plane((0.25, 0.0, 0.0), X)),
        ('e_sphere_two', conic('ellipse', o, Z, 2.0, 1.0), sphere((1.0, 0.0, 0.0), 1.5)),
        ('e_sphere_tangent', conic('ellipse', o, Z, 2.0, 1.0), sphere((3.0, 0.0, 0.0), 1.0)),
        ('e_cylinder_four', conic('ellipse', o, Y, 2.0, 1.0), cylinder(o, Z, 1.5)),
        ('e_cylinder_tangent', conic('ellipse', o, Z, 2.0, 1.0), cylinder((3.0, 0.0, 0.0), Z, 1.0)),
        ('e_cone_two', conic('ellipse', (1.0, 0.0, 1.0), X, 1.5, 1.0), cone(o, Z, 0.0, 0.5)),
        ('e_cone_miss', conic('ellipse', (0.0, 0.0, 3.0), Z, 0.5, 0.25), cone(o, Z, 0.0, 0.5)),
        ('e_torus_cross', conic('ellipse', o, X, 2.25, 0.25, Y), thin),
        ('e_torus_tangent', conic('ellipse', o, Z, 2.5, 1.0), thin),
        ('e_torus_tilted', conic('ellipse', (0.25, 0.0, 0.0), TILT, 2.25, 1.75), thin),
        # Hyperbolas (S7c.2): one branch.
        ('h_plane_cross', conic('hyperbola', o, Z, 1.0, 1.0), plane((2.0, 0.0, 0.0), X)),
        ('h_plane_tangent', conic('hyperbola', o, Z, 1.0, 1.0), plane((1.0, 0.0, 0.0), X)),
        ('h_plane_other', conic('hyperbola', o, Z, 1.0, 1.0), plane((-2.0, 0.0, 0.0), X)),
        ('h_plane_one', conic('hyperbola', o, Z, 1.0, 0.5), plane((0.0, 0.25, 0.0), Y)),
        ('h_plane_contained', conic('hyperbola', o, Z, 1.0, 1.0), plane((1.0, 2.0, 0.0), Z)),
        ('h_sphere_two', conic('hyperbola', o, Z, 1.0, 1.0), sphere(o, 2.0)),
        ('h_sphere_tangent', conic('hyperbola', o, Z, 1.0, 1.0), sphere(o, 1.0)),
        ('h_cylinder_two', conic('hyperbola', o, Z, 1.0, 1.0), cylinder(o, X, 1.0)),
        ('h_cone_two', conic('hyperbola', o, Z, 1.0, 1.0), cone((0.0, 0.0, -3.0), Z, 0.0, 0.5)),
        ('h_cone_miss', conic('hyperbola', o, Z, 1.0, 1.0), cone((0.0, 0.0, -1.0), Z, 0.0, 0.5)),
        ('h_torus_cross', conic('hyperbola', o, Z, 1.0, 1.0), thin),
        ('h_torus_tilted', conic('hyperbola', (0.25, 0.0, 0.0), TILT, 1.75, 0.75), thin),
        # Splines against tori and cones (S7c.2).
        ('s_torus_cross', spline(2, [(-3.0, 0.0, 0.25, 1.0), (0.0, 0.0, -0.5, 1.0), (3.0, 0.0, 0.25, 1.0),
                                     (3.0, 3.0, 0.0, 1.0)], [0.0, 1.0, 2.0], [3, 1, 3]), thin),
        ('s_torus_top', spline(2, [(-3.0, 0.0, 0.5, 1.0), (0.0, 0.0, 0.5, 1.0), (3.0, 0.0, 0.5, 1.0)],
                               [0.0, 1.0], [3, 3]), thin),
        ('s_torus_contained', spline(2, quarter, [0.0, 1.0], [3, 3]), thin),
        ('s_torus_partial', spline(2, quarter+[(-2.0, 2.0, -0.5, 1.0), (-2.0, 0.0, -1.0, 1.0)],
                                   [0.0, 1.0, 2.0], [3, 2, 3]), thin),
        ('s_torus_miss', spline(3, [(5.0, 5.0, 0.0, 1.0), (6.0, 5.0, 1.0, 2.0), (6.0, 6.0, 1.0, 1.0),
                                    (5.0, 6.0, 0.0, 1.0)], [0.0, 1.0], [4, 4]), thin),
        ('s_cone_cross', spline(3, [(-2.0, 0.25, 1.0, 1.0), (-0.5, 0.25, -1.0, 1.0), (0.5, 0.5, 2.0, 2.0),
                                    (2.0, 0.25, 0.5, 1.0)], [0.0, 1.0], [4, 4]), cone(o, Z, 1.0, 0.5)),
        ('s_cone_apex', spline(2, [(-1.0, 1.0, -1.0, 1.0), (0.0, -1.0, 0.0, 1.0), (1.0, 1.0, 1.0, 1.0)],
                               [0.0, 1.0], [3, 3]), cone(o, Z, 0.0, 0.5)),
        ('s_cone_contained', spline(2, [(1.0, 0.0, 0.0, 1.0), (1.0, 1.0, 0.0, 1.0), (0.0, 1.0, 0.0, 2.0)],
                                    [0.0, 1.0], [3, 3]), cone(o, Z, 1.0, 0.5)),
        ('s_cone_miss', spline(2, [(3.0, 3.0, 0.0, 1.0), (4.0, 3.0, 0.5, 1.0), (4.0, 4.0, 0.25, 1.0)],
                               [0.0, 1.0], [3, 3]), cone(o, Z, 0.5, 0.25)),
    ]


def encode(name, curve, surface):
    kind, values = curve
    if kind == 'spline':
        degree, poles, knots, mults = values
        words = [str(degree), str(len(poles))]+[repr(float(v)) for p in poles for v in p]
        words += [str(len(knots))]+[w for k, m in zip(knots, mults) for w in (repr(float(k)), str(m))]
        rows = [f'case {name}', 'curve spline '+' '.join(words)]
    else:
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
    if kind in ('ellipse', 'hyperbola'):
        return ref.conic_rows(tuple(values[:9]), values[9], values[10], kind == 'hyperbola', surface)
    if kind == 'spline':
        return ref.spline_rows(values, surface)
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
        framed = [surface[1]] if curve[0] in ('line', 'spline') else [curve[1][:9], surface[1]]
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
