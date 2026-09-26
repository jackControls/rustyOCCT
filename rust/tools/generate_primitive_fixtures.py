#!/usr/bin/env python3
"""Cone and frustum fixtures from the independent primitive reference (S3).

primitive-cases.txt lists one cone per line (`cone NAME ox oy oz nx ny nz
xx xy xz r1 r2 h`, the arguments of BRepPrimAPI_MakeCone with a gp_Ax2).
primitive-expected.tsv holds, per case, OCCT's distinct subshape counts,
the volume, surface area, centre of mass and the six entries of the
symmetric inertia matrix about it, and every face, edge and vertex as
`primitive_reference.py` derives them. No kernel or OCCT result supplies an
expectation.
"""
import argparse
import math
from pathlib import Path

import mpmath as mp

from primitive_reference import Cone, counts, edges, faces, mass, vertices

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


def encode(c):
    return ' '.join(['cone', c.name] + [number(v) for v in (*c.origin, *c.normal, *c.x, c.r1, c.r2, c.height)])


def expected(c):
    rows = [f'{c.name}\tcounts\t' + ' '.join(map(str, counts(c)))]
    volume, area, centre, inertia = mass(c)
    entries = [inertia[i][j] for i in range(3) for j in range(i, 3)]
    rows.append(f'{c.name}\tprops\t' + ' '.join(number(v) for v in [volume, area, *centre, *entries]))
    for kind, a, centre in sorted(faces(c), key=lambda f: (f[0], float(f[1]))):
        rows.append(f'{c.name}\tface\t{kind} {number(a)} ' + ' '.join(number(v) for v in centre))
    for degenerate, closed, length, mid in sorted(edges(c), key=lambda e: (e[0], e[1], float(e[2]))):
        rows.append(f'{c.name}\tedge\t{degenerate} {closed} {number(length)} ' + ' '.join(number(v) for v in mid))
    for p in sorted(vertices(c), key=lambda p: [float(v) for v in p]):
        rows.append(f'{c.name}\tvertex\t' + ' '.join(number(v) for v in p))
    return rows


def generate():
    cases = corpus()
    return {'primitive-cases.txt': '\n'.join(encode(c) for c in cases)+'\n',
            'primitive-expected.tsv': '# case\tkind\tvalues\n'
            + '\n'.join(row for c in cases for row in expected(c))+'\n'}


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
    print(f'{len(corpus())} cones')


if __name__ == '__main__':
    main()
