#!/usr/bin/env python3
"""Fixtures for S8a of REVIEW_NOTES.md: prisms of line and arc profiles
split by a plane.

`split-cases.txt` lists each case in the protocol of `identity-cases.txt`
(a prism: `case`, `op`, `frame`, `offsets`, `boundary` rows) with a `split ox
oy oz nx ny nz` row (the plane through a point with a normal) before `end`;
`split-expected.tsv` gives each case's sides from `split_reference.py`
(`side below|above|whole volume area cx cy cz`: the totals of the pieces on
each side), and `split-frames.tsv` each
prism's stored frame axes (the reference's `frame_axes`; the kernel's tests
check them bit for bit). Every class has a case: planes crossing the caps and
the walls at angles, normal to the axis (M3's height split), parallel to it
across line and arc walls, through two vertical edges, through one vertex
only, tangent to an arc wall along a ruling, lying in a cap, missing, and
crossing a hole, a non-convex profile in two pieces of section, and a
circle's whole wall in an ellipse. No Rust result supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

from identity_reference import Boundary, Case, encode_case, frame_axes
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


def encode(case, plane):
    text = encode_case(case)
    body, end = text.rsplit('\nend', 1)
    return body+'\nsplit '+' '.join(repr(float(v)) for v in plane)+'\nend'


def generate():
    blocks, out = [], ['# case\trow (S8a, split_reference.py)']
    frames = ['# case\taxis (n, x, y)\tstored unit vector (the reference\'s frame_axes, as hex bits)']
    for case, plane in cases():
        blocks.append(encode(case, plane))
        for row in ref.rows(case, plane):
            out.append(f'{case.name}\t{ref.text(row)}')
        _, x, y, n = frame_axes(case.frame)
        for key, v in (('n', n), ('x', x), ('y', y)):
            frames.append(f'{case.name}\t{key}\t'+' '.join(struct.pack('>d', c).hex() for c in v))
    return {'split-cases.txt': '\n'.join(blocks)+'\n', 'split-expected.tsv': '\n'.join(out)+'\n',
            'split-frames.tsv': '\n'.join(frames)+'\n'}


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
    print(len(cases()), 'cases')


if __name__ == '__main__':
    main()
