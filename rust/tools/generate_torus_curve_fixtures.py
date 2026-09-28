#!/usr/bin/env python3
"""Fixtures for S7b.3b of REVIEW_NOTES.md: a torus with a cylinder, a cone or
another torus off its axis.

`torus-curve-cases.txt` uses the case protocol of
`analytic-intersection-cases.txt`; `torus-curve-expected.tsv` gives each
case's canonical rows from `torus_curve_reference.py` (folds, tangencies,
components and their winding numbers, rings' points), and
`torus-curve-frames.tsv` each surface's stored unit normal as the reference
computes `Frame3::new`'s (the kernel must store the same bits: every
normal has a zero coordinate and Pythagorean others, which the platform's
`hypot` cannot round differently). Every class
has an exact case (tangencies from outside, inside the hole, over the tube,
at a saddle where two branches cross) beside near ones. No Rust result
supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

import analytic_intersection_reference as ana
import torus_curve_reference as ref
from generate_analytic_intersection_fixtures import X, Y, Z, cone, cylinder, encode, torus
from identity_reference import frame_axes

ROOT = Path(__file__).resolve().parents[1]
TILT = (0.0, 3.0, 4.0)


def cases():
    o = (0.0, 0.0, 0.0)
    near = 2.0**-20
    thin = torus(o, Z, 2.0, 0.5)
    return [
        # A torus (major 2, minor 0.5 about z) and a cylinder.
        ('tc_through_tube', thin, cylinder(o, X, 0.3)),
        ('tc_hole_cut', thin, cylinder((1.0, 0.0, 0.0), Z, 1.0)),
        ('tc_rings', thin, cylinder((0.25, 0.0, 0.0), Z, 2.0)),
        ('tc_across_hole', thin, cylinder(o, Y, 1.75)),
        ('tc_enclosing', thin, cylinder(o, X, 2.25)),
        ('tc_vertical_through', thin, cylinder((2.0, 0.0, 0.0), Z, 0.25)),
        ('tc_skew', thin, cylinder((1.5, 0.0, 0.25), TILT, 0.4)),
        ('tc_tilted', torus((1.0, 2.0, 3.0), TILT, 2.0, 0.5), cylinder((1.0, 2.0, 3.0), X, 0.3)),
        ('tc_miss', thin, cylinder((0.0, 0.0, 5.0), X, 0.5)),
        # Exact tangencies: isolated points and a crossing at a saddle.
        ('tc_tangent_hole', thin, cylinder((0.5, 0.0, 0.0), Z, 1.0)),
        ('tc_tangent_outer', thin, cylinder((3.5, 0.0, 0.0), Z, 1.0)),
        ('tc_tangent_over', thin, cylinder((0.0, 0.0, 1.0), X, 0.5)),
        ('tc_saddle', thin, cylinder((3.0, 0.0, 0.0), Z, 1.5)),
        ('tc_near_outer_in', thin, cylinder((3.5, 0.0, 0.0), Z, 1.0+near)),
        ('tc_near_outer_out', thin, cylinder((3.5, 0.0, 0.0), Z, 1.0-near)),
        ('tc_near_over', thin, cylinder((0.0, 0.0, 1.0), X, 0.5+near)),
        ('tc_near_saddle_wide', thin, cylinder((3.0, 0.0, 0.0), Z, 1.5+near)),
        ('tc_near_saddle_narrow', thin, cylinder((3.0, 0.0, 0.0), Z, 1.5-near)),
        # A torus and a cone (both nappes; never exactly tangent).
        ('tk_skew', thin, cone((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0), 0.0, 0.2)),
        ('tk_parallel', thin, cone((0.5, 0.0, 3.0), (0.0, 0.0, -1.0), 0.0, 0.6)),
        ('tk_through_tube', thin, cone((2.0, 0.0, 2.0), (0.0, 0.0, -1.0), 0.0, 0.1)),
        ('tk_irrational_apex', thin, cone((2.0, 0.0, 0.0), Z, 0.2, 0.3)),
        ('tk_tilted', torus((1.0, 2.0, 3.0), TILT, 2.0, 0.5), cone((3.0, 2.0, 3.0), X, 0.25, 0.3)),
        ('tk_miss', thin, cone((0.0, 0.0, 5.0), X, 0.0, 0.1)),
        # S7b.3b.2: two tori off a common axis.
        ('tt_side', thin, torus((3.0, 0.0, 0.0), Z, 1.0, 0.25)),
        ('tt_link_cut', thin, torus((2.0, 0.0, 0.0), Y, 2.0, 1.625)),
        ('tt_ring_round_tube', thin, torus((2.0, 0.0, 0.0), X, 1.0, 0.625)),
        ('tt_parallel_offset', thin, torus((1.0, 0.0, 0.25), Z, 2.0, 0.5)),
        ('tt_through_hole', thin, torus(o, X, 1.25, 0.375)),
        ('tt_over', thin, torus((0.0, 0.0, 1.0), X, 2.0, 0.5)),
        ('tt_tilted', torus((1.0, 2.0, 3.0), TILT, 2.0, 0.5), torus((3.0, 2.0, 3.0), X, 1.0, 0.375)),
        ('tt_miss', thin, torus((0.0, 0.0, 5.0), X, 1.0, 0.25)),
        # Exact tangencies: outer equators, a saddle and an outer equator, two
        # saddles (crossing), the top of one tube and the bottom of another.
        ('tt_side_touch', thin, torus((5.0, 0.0, 0.0), Z, 2.0, 0.5)),
        ('tt_saddle_touch', thin, torus((3.0, 0.0, 0.0), Z, 1.125, 0.375)),
        ('tt_saddles_touch', thin, torus((3.0, 0.0, 0.0), Z, 2.25, 0.75)),
        ('tt_perpendicular_touch', thin, torus((2.0, 0.0, 2.0), X, 1.0, 0.5)),
        ('tt_side_touch_near', thin, torus((5.0-near, 0.0, 0.0), Z, 2.0, 0.5)),
        ('tt_saddles_near', thin, torus((3.0, 0.0, 0.0), Z, 2.25, 0.75+near)),
    ]


def generate():
    blocks, rows = [], ['# case\tcanonical row (S7b.3b, torus_curve_reference.py)']
    frames = ['# case\tsurface\tstored unit normal (the reference\'s Frame3::new, as hex bits)']
    for name, a, b in cases():
        block = encode(name, a, b)
        blocks.append(block+'\nend')
        _, surfaces = ana.parse(block)
        for row in ref.rows(*surfaces):
            rows.append(f'{name}\t{ref.text(row)}')
        for k, surface in enumerate(surfaces):
            _, _, _, n = frame_axes(surface.frame)
            frames.append(f'{name}\t{k}\t'+' '.join(struct.pack('>d', v).hex() for v in n))
    return {'torus-curve-cases.txt': '\n'.join(blocks)+'\n',
            'torus-curve-expected.tsv': '\n'.join(rows)+'\n',
            'torus-curve-frames.tsv': '\n'.join(frames)+'\n'}


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
