#!/usr/bin/env python3
"""Fixtures for S7b.4 of REVIEW_NOTES.md: two cones (not coaxial), and a cone
whose rational apex lies on a sphere or a cylinder.

`ruled-curve-cases.txt` uses the case protocol of
`analytic-intersection-cases.txt`; `ruled-curve-expected.tsv` gives each
case's canonical rows from `ruled_curve_reference.py` (folds, components with
their winding numbers and crossings of infinity, rings' points, the apex),
and `ruled-curve-frames.tsv` each surface's stored unit normal as the
reference computes `Frame3::new`'s (every normal has a zero coordinate and a
Pythagorean pair, which the platform's `hypot` cannot round differently).
Exact cases: an apex on a sphere (isolated), on a cylinder (crossing, and
with parallel axes), parallel cones of equal half-angles (the circle at
infinity common), two cones with one apex; unbounded components where
rulings are parallel to generatrices. A cone and a sphere or a cylinder with
the apex off it are S7b.2's (`procedural-intersection-cases.txt`). No Rust
result supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

import analytic_intersection_reference as ana
import ruled_curve_reference as ref
import torus_curve_reference as tc
from generate_analytic_intersection_fixtures import X, Y, Z, cone, cylinder, encode, sphere
from identity_reference import frame_axes

ROOT = Path(__file__).resolve().parents[1]
TILT = (0.0, 3.0, 4.0)


def cases():
    o = (0.0, 0.0, 0.0)
    return [
        # Two cones.
        ('kk_cross', cone(o, Z, 0.0, 0.5), cone((3.0, 0.0, 0.0), X, 0.0, 0.4)),
        ('kk_parallel', cone(o, Z, 0.0, 0.5), cone((1.0, 0.0, 0.0), Z, 0.0, 0.3)),
        ('kk_unbounded', cone(o, Z, 0.0, 0.9), cone((2.0, 0.0, 0.0), X, 0.0, 0.8)),
        ('kk_tilted', cone((1.0, 2.0, 3.0), TILT, 0.25, 0.4), cone((3.0, 2.0, 3.0), X, 0.0, 0.3)),
        ('kk_irrational_apexes', cone(o, Z, 0.5, 0.5), cone((3.0, 0.0, 1.0), X, 0.25, 0.4)),
        ('kk_miss', cone(o, Z, 0.0, 0.3), cone((0.0, 2.0, 1.0), X, 0.25, 0.35)),
        ('kk_twins', cone(o, Z, 0.0, 0.5), cone((1.0, 0.0, 0.5), Z, 0.0, 0.5)),
        ('kk_common_apex', cone(o, Z, 0.0, 0.9), cone(o, X, 0.0, 0.8)),
        ('kk_common_apex_point', cone(o, Z, 0.0, 0.5), cone(o, X, 0.0, 0.9)),
        # A rational apex on a sphere or a cylinder.
        ('ka_sphere', cone(o, Z, 0.0, 0.5), sphere((0.75, 0.0, 1.0), 1.25)),
        ('ka_cylinder', cone(o, Z, 0.0, 0.5), cylinder((1.0, 0.0, 0.0), Y, 1.0)),
        ('ka_cylinder_parallel', cone(o, Z, 0.0, 0.2), cylinder((1.0, 0.0, 0.0), Z, 1.0)),
    ]


def generate():
    blocks, rows = [], ['# case\tcanonical row (S7b.4, ruled_curve_reference.py)']
    frames = ['# case\tsurface\tstored unit normal (the reference\'s Frame3::new, as hex bits)']
    for name, a, b in cases():
        block = encode(name, a, b)
        blocks.append(block+'\nend')
        _, surfaces = ana.parse(block)
        for row in ref.rows(*surfaces):
            rows.append(f'{name}\t{tc.text(row)}')
        for k, surface in enumerate(surfaces):
            _, _, _, n = frame_axes(surface.frame)
            frames.append(f'{name}\t{k}\t'+' '.join(struct.pack('>d', v).hex() for v in n))
    return {'ruled-curve-cases.txt': '\n'.join(blocks)+'\n',
            'ruled-curve-expected.tsv': '\n'.join(rows)+'\n',
            'ruled-curve-frames.tsv': '\n'.join(frames)+'\n'}


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
