#!/usr/bin/env python3
"""Fixtures for S7b.1, S7b.2 and S7b.3a of REVIEW_NOTES.md: two cylinders
with crossing axes, a cylinder with a sphere off its axis, a cylinder with a
cone (axes not coaxial), a sphere with a cone (the centre off the axis), and
a torus with a plane, a sphere, or a coaxial cylinder, cone or torus.

`procedural-intersection-cases.txt` uses the case protocol of
`analytic-intersection-cases.txt`; `procedural-intersection-expected.tsv`
gives each case's canonical rows from `procedural_intersection_reference.py`
and `procedural-intersection-frames.tsv` each surface's stored unit normal as
the reference computes `Frame3::new`'s (the kernel must store the same bits;
near-degenerate cases are sensitive to one unit in the last place).
Every class (empty, tangent point, loop, figure-eight, rings) has an exact
case, beside near-tangent and near-figure-eight ones and tilted frames. No
Rust result supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

import analytic_intersection_reference as ana
import procedural_intersection_reference as ref
from identity_reference import frame_axes
from generate_analytic_intersection_fixtures import X, Y, Z, cone, cylinder, encode, plane, sphere, torus

ROOT = Path(__file__).resolve().parents[1]
TILT = (0.0, 3.0, 4.0)


def cases():
    o = (0.0, 0.0, 0.0)
    near = 2.0**-20
    thin, fat = torus(o, Z, 2.0, 0.5), torus(o, Z, 2.0, 1.0)
    return [
        ('cc_loop_equal', cylinder(o, Z, 1.0), cylinder((0.0, 0.5, 0.0), X, 1.0)),
        ('cc_loop_unequal', cylinder(o, Z, 1.0), cylinder((0.0, 0.8, 0.0), X, 0.5)),
        ('cc_rings', cylinder(o, Z, 2.0), cylinder((0.0, 0.5, 0.0), X, 1.0)),
        ('cc_rings_crossing', cylinder(o, Z, 2.0), cylinder(o, X, 1.0)),
        ('cc_tangent_point', cylinder(o, Z, 1.0), cylinder((0.0, 3.0, 0.0), X, 2.0)),
        ('cc_figure_eight', cylinder(o, Z, 2.0), cylinder((0.0, 1.0, 0.0), X, 1.0)),
        ('cc_empty', cylinder(o, Z, 1.0), cylinder((0.0, 5.0, 0.0), X, 1.0)),
        ('cc_oblique_loop', cylinder(o, Z, 1.0), cylinder((0.0, 0.25, 0.0), (1.0, 0.0, 1.0), 1.0)),
        ('cc_tilted', cylinder((1.0, 2.0, 3.0), TILT, 1.5), cylinder((1.0, 2.5, 3.0), X, 1.0)),
        ('cc_near_tangent', cylinder(o, Z, 1.0), cylinder((0.0, 3.0-near, 0.0), X, 2.0)),
        ('cc_near_figure_inside', cylinder(o, Z, 2.0), cylinder((0.0, 1.0-near, 0.0), X, 1.0)),
        ('cc_near_figure_outside', cylinder(o, Z, 2.0), cylinder((0.0, 1.0+near, 0.0), X, 1.0)),
        ('cs_loop', cylinder(o, Z, 1.0), sphere((1.5, 0.0, 0.0), 1.0)),
        ('cs_viviani', cylinder((1.0, 0.0, 0.0), Z, 1.0), sphere(o, 2.0)),
        ('cs_rings', cylinder((0.5, 0.0, 0.0), Z, 1.0), sphere(o, 3.0)),
        ('cs_tangent_outside', cylinder((3.0, 0.0, 0.0), Z, 1.0), sphere(o, 2.0)),
        ('cs_tangent_inside', cylinder((1.0, 0.0, 0.0), Z, 3.0), sphere(o, 2.0)),
        ('cs_empty', cylinder((5.0, 0.0, 0.0), Z, 1.0), sphere(o, 2.0)),
        ('cs_tilted', cylinder((1.0, 2.0, 3.0), TILT, 1.0), sphere((1.75, 2.0, 3.0), 1.5)),
        ('cs_near_viviani', cylinder((1.0, 0.0, 0.0), Z, 1.0), sphere(o, 2.0-near)),
        # S7b.2: a cylinder and a cone.
        ('ck_cross', cylinder(o, Z, 1.0), cone((-3.0, 0.0, 0.0), X, 0.0, 0.3)),
        ('ck_parallel', cylinder((0.5, 0.0, 0.0), Z, 1.0), cone(o, Z, 0.0, 0.4)),
        ('ck_skew', cylinder((1.0, 0.0, 0.0), Z, 0.5), cone((0.0, 0.0, 1.0), Y, 0.2, 0.5)),
        ('ck_miss', cylinder((0.0, 5.0, 0.0), X, 1.0), cone(o, Z, 0.0, 0.1)),
        ('ck_rings', cylinder((0.0, 0.0, 5.0), X, 0.5), cone(o, Z, 0.0, 1.2)),
        ('ck_tilted', cylinder((1.0, 2.0, 3.0), TILT, 0.75), cone((1.0, 2.5, 3.0), X, 0.5, 0.35)),
        ('ck_two_loops', cylinder((0.0, 0.0, 2.0), X, 1.0), cone(o, Z, 0.0, 0.4)),
        # S7b.2: a sphere and a cone, the apex rational or not.
        ('ks_apex_inside', cone(o, Z, 0.0, 0.5), sphere((0.5, 0.0, 0.5), 2.0)),
        ('ks_one_loop', cone(o, Z, 0.0, 0.5), sphere((3.0, 0.0, 5.0), 1.0)),
        ('ks_two_loops', cone(o, Z, 0.0, 0.5), sphere((1.5, 0.0, 0.0), 1.4)),
        # The centre is 1.5 cos 0.5 = 1.31637... from each generatrix.
        ('ks_near_tangent', cone(o, Z, 0.0, 0.5), sphere((1.5, 0.0, 0.0), 1.3163738428355591+near)),
        ('ks_rings_one_nappe', cone(o, Z, 0.0, 0.3), sphere((0.2, 0.0, 5.0), 3.0)),
        ('ks_irrational_apex', cone(o, Z, 0.5, 0.4), sphere((1.0, 0.0, 2.0), 1.0)),
        ('ks_tilted', cone((1.0, 2.0, 3.0), TILT, 0.25, 0.45), sphere((2.0, 2.0, 3.0), 1.5)),
        ('ks_miss', cone(o, Z, 0.0, 0.2), sphere((5.0, 0.0, 0.0), 1.0)),
        # S7b.3a: a torus and a plane (major 2, minor 0.5 about z).
        ('tp_normal_two', plane((0.0, 0.0, 0.25), Z), thin),
        ('tp_normal_tangent', plane((0.0, 0.0, 0.5), Z), thin),
        ('tp_normal_miss', plane((0.0, 0.0, 1.0), Z), thin),
        ('tp_axis', plane(o, X), thin),
        ('tp_axis_offset', plane((0.0, 0.0, 3.0), (1.0, 1.0, 0.0)), thin),
        ('tp_side_loop', plane((2.25, 0.0, 0.0), X), thin),
        ('tp_side_node', plane((1.5, 0.0, 0.0), X), thin),
        ('tp_side_two_loops', plane((1.0, 0.0, 0.0), X), thin),
        ('tp_side_tangent', plane((2.5, 0.0, 0.0), X), thin),
        ('tp_side_miss', plane((3.0, 0.0, 0.0), X), thin),
        ('tp_near_side_tangent', plane((2.5-near, 0.0, 0.0), X), thin),
        ('tp_near_node_outside', plane((1.5+near, 0.0, 0.0), X), thin),
        ('tp_near_node_inside', plane((1.5-near, 0.0, 0.0), X), thin),
        ('tp_oblique_rings', plane(o, (0.0, 1.0, 8.0)), thin),
        ('tp_oblique_two_loops', plane(o, (0.0, 1.0, 1.0)), thin),
        ('tp_near_villarceau', plane(o, TILT), torus(o, Z, 2.0, 1.2)),
        ('tp_oblique_offset', plane((0.0, 1.0, 0.25), TILT), thin),
        ('tp_tilted', plane((1.0, 2.0, 3.0), X), torus((1.5, 2.0, 3.0), TILT, 2.0, 0.5)),
        # A torus and a sphere (major 2, minor 1 about z): every class exactly.
        ('ts_rings', sphere((0.25, 0.0, 0.0), 1.5), fat),
        ('ts_loop', sphere((0.25, 0.0, 0.0), 1.0), fat),
        ('ts_loop_across', sphere((0.5, 0.0, 0.0), 3.0), fat),
        ('ts_two_loops', sphere((1.25, 0.0, 0.0), 2.0), fat),
        ('ts_touching_zero', sphere((1.25, 0.0, 0.0), 1.75), fat),
        ('ts_touching_pi', sphere((1.25, 0.0, 0.0), 2.25), fat),
        ('ts_figure_zero', sphere((0.25, 0.0, 0.0), 2.75), fat),
        ('ts_figure_pi', sphere((0.25, 0.0, 0.0), 1.25), fat),
        ('ts_point_zero', sphere((0.25, 0.0, 0.0), 0.75), fat),
        ('ts_point_pi', sphere((0.25, 0.0, 0.0), 3.25), fat),
        ('ts_empty', sphere((0.25, 0.0, 0.0), 0.25), fat),
        ('ts_inside_tube', sphere((1.75, 0.0, 0.0), 0.25), fat),
        ('ts_villarceau_flat', sphere((1.0, 0.0, 0.0), 2.0), fat),
        ('ts_villarceau', sphere((2.0, 0.0, 3.0), 4.0), fat),
        ('ts_meridian', sphere((2.0, 0.0, 0.0), 1.0), fat),
        ('ts_two_meridians', sphere((2.0, 0.75, 0.0), 1.25), fat),
        ('ts_near_figure', sphere((0.25, 0.0, 0.0), 2.75+near), fat),
        ('ts_near_point', sphere((0.25, 0.0, 0.0), 0.75+near), fat),
        ('ts_near_villarceau', sphere((2.0, 0.0, 3.0), 4.0-near), fat),
        ('ts_coaxial_two', sphere(o, 2.0), fat),
        ('ts_coaxial_inner', sphere(o, 1.0), fat),
        ('ts_coaxial_outer', sphere(o, 3.0), fat),
        ('ts_coaxial_miss', sphere(o, 4.0), fat),
        ('ts_coaxial_above', sphere((0.0, 0.0, 1.0), 2.0), fat),
        ('ts_tilted', sphere((2.0, 2.5, 3.0), 1.5), torus((1.0, 2.0, 3.0), TILT, 2.0, 0.75)),
        # Coaxial pairs (the torus of major 2, minor 0.5 about z).
        ('tx_cylinder_two', cylinder(o, Z, 2.25), thin),
        ('tx_cylinder_outer', cylinder((0.0, 0.0, 5.0), Z, 2.5), thin),
        ('tx_cylinder_inner', cylinder(o, (0.0, 0.0, -1.0), 1.5), thin),
        ('tx_cylinder_miss', cylinder(o, Z, 3.0), thin),
        ('tx_cone_two', cone((0.0, 0.0, -4.0), Z, 0.0, 0.5), thin),
        ('tx_cone_four', cone(o, Z, 0.0, 1.373400766945016), thin),
        ('tx_cone_reversed', cone((0.0, 0.0, 4.0), (0.0, 0.0, -1.0), 0.0, 0.5), thin),
        ('tx_cone_irrational_apex', cone(o, Z, 1.75, 0.5), thin),
        ('tx_cone_miss', cone(o, Z, 0.0, 0.1), thin),
        ('tx_torus_two', torus((0.0, 0.0, 0.5), Z, 2.0, 0.5), thin),
        ('tx_torus_tangent', torus((0.0, 0.0, 1.0), Z, 2.0, 0.5), thin),
        ('tx_torus_nested_tangent', torus(o, Z, 2.5, 1.0), thin),
        ('tx_torus_miss', torus(o, Z, 4.0, 0.5), thin),
        ('tx_torus_concentric', torus(o, Z, 2.0, 0.25), thin),
        ('tx_torus_same', torus(o, (0.0, 0.0, -1.0), 2.0, 0.5), thin),
    ]


def generate():
    blocks, rows = [], ['# case\tcanonical row (S7b, procedural_intersection_reference.py)']
    frames = ['# case\tsurface\tstored unit normal (the reference\'s Frame3::new, as hex bits)']
    for name, a, b in cases():
        block = encode(name, a, b)
        blocks.append(block+'\nend')
        _, surfaces = ana.parse(block)
        found = ref.rows(*surfaces)
        assert found is not None, name
        for row in found:
            rows.append(f'{name}\t{ref.text(row)}')
        for k, surface in enumerate(surfaces):
            _, _, _, n = frame_axes(surface.frame)
            frames.append(f'{name}\t{k}\t'+' '.join(struct.pack('>d', v).hex() for v in n))
    return {'procedural-intersection-cases.txt': '\n'.join(blocks)+'\n',
            'procedural-intersection-expected.tsv': '\n'.join(rows)+'\n',
            'procedural-intersection-frames.tsv': '\n'.join(frames)+'\n'}


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
