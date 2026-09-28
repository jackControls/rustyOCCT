#!/usr/bin/env python3
"""Fixtures for tessellation (T-a of REVIEW_NOTES.md).

`tessellation-cases.txt` holds identity-protocol case blocks
(identity_reference.encode_case: prisms of profiles with lines, arcs,
circles and holes, a box, cones, spheres, tori and face bodies), each with
`mesh SETTING DEFLECTION ANGLE` rows before its `end`: a coarse setting (a
hundredth of the case's scale, 0.5 rad, BRepMesh's default angle) and a fine
one (a thousandth, 0.3 rad). `tessellation-expected.tsv` gives what
tessellation_reference.py derives from each case alone: whether it is a
solid, the Euler characteristic of its boundary, a face's number of boundary
loops, exact area and volume, and the scale the settings are relative to.
No Rust or OCCT result supplies an expectation.
"""
import argparse
from pathlib import Path

import mpmath as mp

from identity_reference import Boundary, Case, encode_case
import tessellation_reference as ref

ROOT = Path(__file__).resolve().parents[1]
Z = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
TILT = (1.0, -2.0, 0.5, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0)
FAR = (10000.0, -20000.0, 5000.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
H = ref.HALF_PI
T = ref.TWO_PI


def rect(w, d, x=0.0, y=0.0):
    return Boundary(points=[(x, y), (x+w, y), (x+w, y+d), (x, y+d)])


def path(points_segments):
    """[(point, segment)]: segment None for a line, (cx, cy, r, ccw) for an arc."""
    return Boundary(points=[p for p, _ in points_segments], segments=[s for _, s in points_segments])


def stadium(x0, x1, r):
    return path([((x0, -r), None), ((x1, -r), (x1, 0.0, r, True)),
                 ((x1, r), None), ((x0, r), (x0, 0.0, r, True))])


def prism(name, boundaries, start, end, frame=Z):
    return Case(name, 1e-7, 1, frame, start, end, boundaries)


def cases():
    rounded = path([((1.0, 0.0), None), ((7.0, 0.0), (7.0, 1.0, 1.0, True)), ((8.0, 1.0), None),
                    ((8.0, 5.0), (7.0, 5.0, 1.0, True)), ((7.0, 6.0), None),
                    ((1.0, 6.0), (1.0, 5.0, 1.0, True)), ((0.0, 5.0), None),
                    ((0.0, 1.0), (1.0, 1.0, 1.0, True))])
    notch = path([((0.0, 0.0), None), ((10.0, 0.0), None), ((10.0, 6.0), None),
                  ((6.0, 6.0), (5.0, 6.0, 1.0, False)), ((4.0, 6.0), None), ((0.0, 6.0), None)])
    lens = path([((0.0, -3.0), (-4.0, 0.0, 5.0, True)), ((0.0, 3.0), (4.0, 0.0, 5.0, True))])
    scallop = path([((0.0, 0.0), (2.0, 1.5, 2.5, True)), ((4.0, 0.0), (2.5, 2.0, 2.5, True)),
                    ((4.0, 4.0), (2.0, 2.5, 2.5, True)), ((0.0, 4.0), (1.5, 2.0, 2.5, True))])
    square = rect(20.0, 20.0, -10.0, -10.0)
    mixed = [square,
             path([((-6.0, -6.0), None), ((-2.0, -6.0), (-2.0, -4.0, 2.0, True)), ((-2.0, -2.0), None)]),
             Boundary(circle=(5.0, 5.0, 1.5)),
             path([((4.0, -6.0), (5.0, -6.0, 1.0, True)), ((6.0, -6.0), (5.0, -6.0, 1.0, True))])]
    ell = Boundary(points=[(0.0, 0.0), (12.0, 0.0), (12.0, 4.0), (4.0, 4.0), (4.0, 10.0), (0.0, 10.0)])
    out = [
        Case('box', 1e-7, 1, Z, 0.0, 0.0, [], box=((1.0, 2.0, 3.0), (40.0, 20.0, 10.0))),
        prism('plate_hole', [rect(40.0, 20.0), Boundary(circle=(10.0, 10.0, 3.0))], 0.0, 5.0),
        prism('cylinder', [Boundary(circle=(0.0, 0.0, 5.0))], 0.0, 12.0),
        prism('tall_cylinder', [Boundary(circle=(0.0, 0.0, 1.0))], 0.0, 40.0),
        prism('thin_disc', [Boundary(circle=(0.0, 0.0, 10.0))], 0.0, 0.5),
        prism('stadium_slot', [stadium(0.0, 12.0, 4.0), rect(2.0, 2.0, 5.0, -1.0)], 0.0, 2.0),
        prism('rounded_rectangle', [rounded, Boundary(circle=(4.0, 3.0, 1.0))], 0.0, 3.0),
        prism('notch', [notch], 0.0, 4.0),
        prism('lens', [lens], 0.0, 2.0),
        prism('scallop', [scallop], 0.0, 1.0),
        prism('holes_mixed', mixed, 0.0, 3.0),
        prism('tilted_prism', [ell, Boundary(circle=(2.0, 7.0, 1.0))], -1.0, 4.0, TILT),
        prism('far_prism', [Boundary(circle=(0.0, 0.0, 2.0))], 0.0, 3.0, FAR),
        Case('cone_apex', 1e-7, 1, Z, 0.0, 0.0, [], cone=(2.0, 0.0, 3.0)),
        Case('frustum', 1e-7, 1, Z, 0.0, 0.0, [], cone=(2.0, 1.0, 3.0)),
        Case('cone_inverted', 1e-7, 1, Z, 0.0, 0.0, [], cone=(0.0, 1.5, 2.0)),
        Case('cone_tilted', 1e-7, 1, TILT, 0.0, 0.0, [], cone=(1.0, 2.5, 0.5)),
        Case('sphere', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, -H, H)),
        Case('hemisphere', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, 0.0, H)),
        Case('sphere_zone', 1e-7, 1, Z, 0.0, 0.0, [], sphere=(5.0, -0.5, 0.8)),
        Case('sphere_far', 1e-7, 1, FAR, 0.0, 0.0, [], sphere=(2.0, -H, H)),
        Case('torus', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, 0.0, T, T)),
        Case('torus_segment', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, -1.0, 1.0, T)),
        Case('torus_inner_half', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, H, 4.71238898038469, T)),
        Case('torus_wedge', 1e-7, 1, Z, 0.0, 0.0, [], torus=(6.0, 2.0, 0.0, T, 2.0)),
        Case('torus_tilted', 1e-7, 1, TILT, 0.0, 0.0, [], torus=(3.0, 1.0, 0.0, T, T)),
        Case('face_holes', 1e-7, 1, Z, 0.0, 0.0,
             [rect(20.0, 20.0, -10.0, -10.0), rect(4.0, 4.0, -7.0, -7.0),
              Boundary(circle=(5.0, 5.0, 1.5)), stadium(0.0, 3.0, 1.0)], make='face'),
        Case('face_lens_tilted', 1e-7, 1, TILT, 0.0, 0.0, [lens], make='face'),
    ]
    return out


def settings(case):
    """(label, deflection, angle) rows of a case, relative to its scale."""
    scale = ref.body(case).scale
    return [('coarse', scale/100, 0.5), ('fine', scale/1000, 0.3)]


def number(x):
    return repr(float(x))


def encode(case):
    block = encode_case(case)
    rows = [f'mesh {label} {number(d)} {number(a)}' for label, d, a in settings(case)]
    head, end = block.rsplit('\nend', 1)
    return head+'\n'+'\n'.join(rows)+'\nend'+end


def generate():
    blocks = []
    rows = ['# case\tkind\teuler\tloops\tarea\tvolume\tscale (tessellation_reference.py)']
    for c in cases():
        blocks.append(encode(c))
        b = ref.body(c)
        rows.append('\t'.join([c.name, 'solid' if b.solid else 'face', str(b.euler), str(b.loops),
                               mp.nstr(b.area, 17), mp.nstr(b.volume, 17), number(b.scale)]))
    return {'tessellation-cases.txt': '\n'.join(blocks)+'\n',
            'tessellation-expected.tsv': '\n'.join(rows)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = generate()
    for name, contents in files.items():
        path_ = ROOT/'fixtures'/name
        if args.check:
            if path_.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path_.write_text(contents)
    print(len(cases()), 'cases')


if __name__ == '__main__':
    main()
