#!/usr/bin/env python3
"""Fixtures for S7b.1 of REVIEW_NOTES.md: two cylinders with crossing axes
and a cylinder with a sphere off its axis.

`procedural-intersection-cases.txt` uses the case protocol of
`analytic-intersection-cases.txt`; `procedural-intersection-expected.tsv`
gives each case's canonical rows from `procedural_intersection_reference.py`.
Every class (empty, tangent point, loop, figure-eight, rings) has an exact
case, beside near-tangent and near-figure-eight ones and tilted frames. No
Rust result supplies an expectation.
"""
import argparse
from pathlib import Path

import analytic_intersection_reference as ana
import procedural_intersection_reference as ref
from generate_analytic_intersection_fixtures import X, Z, cylinder, encode, sphere

ROOT = Path(__file__).resolve().parents[1]
TILT = (0.0, 3.0, 4.0)


def cases():
    o = (0.0, 0.0, 0.0)
    near = 2.0**-20
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
    ]


def generate():
    blocks, rows = [], ['# case\tcanonical row (S7b.1, procedural_intersection_reference.py)']
    for name, a, b in cases():
        block = encode(name, a, b)
        blocks.append(block+'\nend')
        _, surfaces = ana.parse(block)
        found = ref.rows(*surfaces)
        assert found is not None, name
        for row in found:
            rows.append(f'{name}\t{ref.text(row)}')
    return {'procedural-intersection-cases.txt': '\n'.join(blocks)+'\n',
            'procedural-intersection-expected.tsv': '\n'.join(rows)+'\n'}


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
