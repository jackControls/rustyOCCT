#!/usr/bin/env python3
"""Fixtures for S7a of REVIEW_NOTES.md: analytic surface pairs whose
intersection is empty, the same surface, points, lines or a conic, and
pairs that are not (S7b).

`analytic-intersection-cases.txt` lists each case's two surfaces as
`surface KIND ox oy oz nx ny nz xx xy xz [radius [half-angle]]` (the frame
as `Frame3::new` takes it); `analytic-intersection-expected.tsv` gives each
case's canonical items from `analytic_intersection_reference.py`, and
`analytic-intersection-frames.tsv` each surface's stored unit normal as the
reference computes `Frame3::new`'s (the kernel must store the same bits:
inputs are chosen so that it does, e.g. `(0, 3, 4)` rather than
`(0, 0.6, 0.8)`, whose normalisation depends on the platform's hypot). Every
degeneracy class has an exact case (the predicate holds on the stored
binary64 data) and a near one. No Rust result supplies an expectation.
"""
import argparse
from pathlib import Path
import struct

import analytic_intersection_reference as ref
from identity_reference import frame_axes

ROOT = Path(__file__).resolve().parents[1]
Z, X, Y = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
TILT = (0.0, 3.0, 4.0)
TILT_ACROSS = (0.0, 4.0, -3.0)


def plane(o, n, x=None):
    return ('plane', (*o, *n, *(x or hint(n))))


def cylinder(o, n, r, x=None):
    return ('cylinder', (*o, *n, *(x or hint(n)), r))


def sphere(o, r):
    return ('sphere', (*o, *Z, *X, r))


def cone(o, n, r, a, x=None):
    return ('cone', (*o, *n, *(x or hint(n)), r, a))


def torus(o, n, major, minor, x=None):
    return ('torus', (*o, *n, *(x or hint(n)), major, minor))


def hint(n):
    """An x hint not parallel to the normal."""
    return Y if abs(n[0]) >= max(abs(n[1]), abs(n[2])) else X


def cases():
    o = (0.0, 0.0, 0.0)
    near = 2.0**-30
    a = 0.5
    out = [
        # Plane/plane.
        ('pp_line', plane(o, Z), plane(o, X)),
        ('pp_tilted', plane((1.0, 2.0, 3.0), TILT), plane(o, (1.0, 1.0, 0.0))),
        ('pp_parallel', plane(o, Z), plane((0.0, 0.0, 1.0), Z)),
        ('pp_same', plane(o, Z), plane((5.0, -3.0, 0.0), Z)),
        ('pp_same_opposite', plane(o, Z), plane((2.0, 7.0, 0.0), (0.0, 0.0, -1.0))),
        ('pp_near_parallel', plane(o, Z), plane((0.0, 0.0, 1.0), (near, 0.0, 1.0))),
        # Plane/sphere.
        ('ps_circle', plane((0.0, 0.0, 1.0), Z), sphere(o, 2.0)),
        ('ps_tangent', plane((0.0, 0.0, 2.0), Z), sphere(o, 2.0)),
        ('ps_miss', plane((0.0, 0.0, 3.0), Z), sphere(o, 2.0)),
        ('ps_near_tangent', plane((0.0, 0.0, 2.0-near), Z), sphere(o, 2.0)),
        ('ps_tilted', plane((1.0, 2.0, 3.0), TILT), sphere((1.0, 1.5, 3.25), 4.0)),
        ('ps_through_centre', plane((1.0, 1.0, 1.0), TILT), sphere((1.0, 1.0, 1.0), 3.0)),
        # Plane/cylinder.
        ('pc_perpendicular', plane((0.0, 0.0, 2.0), Z), cylinder(o, Z, 1.5)),
        ('pc_oblique', plane((0.0, 0.0, 2.0), TILT), cylinder(o, Z, 1.5)),
        ('pc_parallel_two', plane((0.0, 0.0, 0.5), Z), cylinder(o, X, 1.0)),
        ('pc_parallel_tangent', plane((0.0, 0.0, 1.0), Z), cylinder(o, X, 1.0)),
        ('pc_parallel_miss', plane((0.0, 0.0, 2.0), Z), cylinder(o, X, 1.0)),
        ('pc_axis_in_plane', plane(o, Z), cylinder((0.0, 3.0, 0.0), X, 1.0)),
        ('pc_parallel_tilted', plane((1.0, 2.0, 3.0), TILT), cylinder((1.0, 2.0, 3.0), TILT_ACROSS, 2.0)),
        ('pc_near_parallel', plane(o, Z), cylinder((0.0, 0.0, 0.5), (1.0, 0.0, near), 1.0)),
        ('pc_near_perpendicular', plane(o, Z), cylinder(o, (near, 0.0, 1.0), 1.0)),
        # Plane/cone (radius 1 at the origin, half-angle 0.5, unless apexed).
        ('pk_circle', plane((0.0, 0.0, 2.0), Z), cone(o, Z, 1.0, a)),
        ('pk_circle_other_nappe', plane((0.0, 0.0, -5.0), Z), cone(o, Z, 1.0, a)),
        ('pk_ellipse', plane((0.0, 0.0, 2.0), (0.0, 0.3, 1.0)), cone(o, Z, 1.0, a)),
        ('pk_hyperbola', plane((0.0, 0.0, 2.0), (0.0, 1.0, 0.1)), cone(o, Z, 1.0, a)),
        ('pk_axis_plane', plane(o, X), cone(o, Z, 1.0, a)),
        ('pk_apex_point', plane(o, Z), cone(o, Z, 0.0, a)),
        ('pk_apex_shallow', plane(o, (0.0, 0.2, 1.0)), cone(o, Z, 0.0, a)),
        ('pk_apex_steep', plane(o, (0.0, 1.0, 0.2)), cone(o, Z, 0.0, a)),
        ('pk_apex_circle', plane((0.0, 0.0, 1.0), Z), cone(o, Z, 0.0, a)),
        ('pk_near_parabola', plane((0.0, 0.0, 3.0), (0.8775825618903728, 0.0, -0.479425538604203)),
         cone(o, Z, 1.0, a)),
        ('pk_tilted', plane((1.0, 2.0, 3.0), TILT), cone((1.0, 2.0, 0.0), Z, 0.75, 0.3)),
        # Sphere/sphere.
        ('ss_circle', sphere(o, 2.0), sphere((2.0, 0.0, 0.0), 2.0)),
        ('ss_external_tangent', sphere(o, 1.0), sphere((3.0, 0.0, 0.0), 2.0)),
        ('ss_internal_tangent', sphere(o, 3.0), sphere((1.0, 0.0, 0.0), 2.0)),
        ('ss_disjoint', sphere(o, 1.0), sphere((5.0, 0.0, 0.0), 2.0)),
        ('ss_contained', sphere(o, 5.0), sphere((1.0, 0.0, 0.0), 1.0)),
        ('ss_same', sphere((1.0, 2.0, 3.0), 2.0), sphere((1.0, 2.0, 3.0), 2.0)),
        ('ss_concentric', sphere((1.0, 2.0, 3.0), 2.0), sphere((1.0, 2.0, 3.0), 1.0)),
        ('ss_tilted', sphere((1.0, 2.0, 3.0), 2.5), sphere((2.0, 3.5, 2.25), 1.75)),
        ('ss_near_tangent', sphere(o, 1.0), sphere((3.0-near, 0.0, 0.0), 2.0)),
        # Cylinder/cylinder.
        ('cc_parallel_two', cylinder(o, Z, 1.0), cylinder((1.0, 0.0, 0.0), Z, 1.0)),
        ('cc_parallel_external', cylinder(o, Z, 1.0), cylinder((3.0, 0.0, 0.0), Z, 2.0)),
        ('cc_parallel_internal', cylinder(o, Z, 3.0), cylinder((1.0, 0.0, 0.0), Z, 2.0)),
        ('cc_parallel_miss', cylinder(o, Z, 1.0), cylinder((5.0, 0.0, 0.0), Z, 1.0)),
        ('cc_coaxial_same', cylinder(o, Z, 1.0), cylinder((0.0, 0.0, 7.0), Z, 1.0)),
        ('cc_coaxial_other', cylinder(o, Z, 1.0), cylinder((0.0, 0.0, 7.0), Z, 2.0)),
        ('cc_parallel_tilted', cylinder((1.0, 2.0, 3.0), TILT, 1.0),
         cylinder((1.0, 2.0, 3.0+1.5), TILT, 1.25)),
        ('cc_equal_cross', cylinder(o, Z, 1.0), cylinder(o, X, 1.0)),
        ('cc_equal_cross_oblique', cylinder(o, Z, 1.0), cylinder((0.0, 0.0, 2.0), (1.0, 0.0, 1.0), 1.0)),
        ('cc_unequal_cross', cylinder(o, Z, 1.0), cylinder(o, X, 0.5)),
        ('cc_skew', cylinder(o, Z, 1.0), cylinder((0.0, 0.5, 0.0), X, 1.0)),
        ('cc_near_parallel', cylinder(o, Z, 1.0), cylinder((1.0, 0.0, 0.0), (near, 0.0, 1.0), 1.0)),
        # Coaxial pairs and spheres on an axis.
        ('cs_two_circles', cylinder(o, Z, 1.0), sphere((0.0, 0.0, 1.0), 2.0)),
        ('cs_tangent', cylinder(o, Z, 1.0), sphere((0.0, 0.0, 1.0), 1.0)),
        ('cs_miss', cylinder(o, Z, 1.0), sphere((0.0, 0.0, 1.0), 0.5)),
        ('cs_tilted', cylinder((1.0, 2.0, 3.0), TILT, 1.0), sphere((1.0, 2.0, 3.0), 2.0)),
        ('cs_off_axis', cylinder(o, Z, 1.0), sphere((0.5, 0.0, 0.0), 2.0)),
        ('ks_two_circles', cone(o, Z, 0.0, a), sphere((0.0, 0.0, 3.0), 2.0)),
        ('ks_through_apex', cone(o, Z, 0.0, a), sphere((0.0, 0.0, 2.0), 2.0)),
        ('ks_off_axis', cone(o, Z, 1.0, a), sphere((0.5, 0.0, 0.0), 2.0)),
        ('ck_two_circles', cylinder(o, Z, 2.0), cone(o, Z, 1.0, a)),
        ('ck_parallel_off_axis', cylinder((0.5, 0.0, 0.0), Z, 2.0), cone(o, Z, 1.0, a)),
        ('kk_two_circles', cone(o, Z, 1.0, a), cone(o, Z, 2.0, 0.3)),
        ('kk_same', cone(o, Z, 1.0, a), cone(o, Z, 1.0, a)),
        ('kk_other_origin', cone(o, Z, 1.0, a), cone((0.0, 0.0, 1.0), Z, 1.5, a)),
        ('kk_crossing', cone(o, Z, 1.0, a), cone(o, X, 1.0, a)),
    ]
    return out


def encode(name, *surfaces):
    rows = [f'case {name}']
    for kind, values in surfaces:
        rows.append(f'surface {kind} '+' '.join(repr(float(v)) for v in values))
    return '\n'.join(rows)


def generate():
    blocks, rows = [], ['# case\tcanonical item (S7a, analytic_intersection_reference.py)']
    frames = ['# case\tsurface\tstored unit normal (the reference\'s Frame3::new, as hex bits)']
    for name, a, b in cases():
        block = encode(name, a, b)
        blocks.append(block+'\nend')
        _, surfaces = ref.parse(block)
        for item in ref.canonical(ref.intersect(*surfaces)):
            rows.append(f'{name}\t{ref.text(item)}')
        for k, surface in enumerate(surfaces):
            _, _, _, n = frame_axes(surface.frame)
            frames.append(f'{name}\t{k}\t'+' '.join(struct.pack('>d', v).hex() for v in n))
    return {'analytic-intersection-cases.txt': '\n'.join(blocks)+'\n',
            'analytic-intersection-expected.tsv': '\n'.join(rows)+'\n',
            'analytic-intersection-frames.tsv': '\n'.join(frames)+'\n'}


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
