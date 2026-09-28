#!/usr/bin/env python3
"""Fixtures for S7d.1 of REVIEW_NOTES.md: pairs of lines, circles, ellipses
and hyperbolas' branches.

`curve-curve-cases.txt` lists each case as `case NAME`, two `curve` rows in
the protocol of `curve-surface-cases.txt` (`line x0 y0 z0 x1 y1 z1`; `circle
ox oy oz nx ny nz hx hy hz radius`; `ellipse|hyperbola ... major minor`) and
`end`; `curve-curve-expected.tsv` gives each case's rows from
`curve_curve_reference.py`, and `curve-curve-frames.tsv` every framed
curve's stored normal and axes (`stored_axes`, the kernel's `Frame3::new`;
its tests check them bit for bit). Every class has an exact case: crossing
and parallel, skew and coincident lines; a line piercing a conic's plane on
the curve and off it, lying in it across the curve, touching it, missing
it, meeting the hyperbola's other branch; two conics in one plane crossing,
touching inside and outside, coincident, apart, on opposite branches; two
conics in crossing planes meeting at two points, touching with a common
tangent line, meeting once, missing. No Rust result supplies an
expectation.
"""
import argparse
from pathlib import Path
import struct

import curve_curve_reference as ref
from curve_surface_reference import stored_axes

ROOT = Path(__file__).resolve().parents[1]
X, Y, Z = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
TILT = (0.0, 3.0, 4.0)


def hint(n):
    return Y if abs(n[0]) >= max(abs(n[1]), abs(n[2])) else X


def line(p0, p1):
    return ('line', (*p0, *p1))


def circle(o, n, r, x=None):
    return ('circle', (*o, *n, *(x or hint(n)), r))


def ellipse(o, n, a, b, x=None):
    return ('ellipse', (*o, *n, *(x or hint(n)), a, b))


def hyperbola(o, n, a, b, x=None):
    return ('hyperbola', (*o, *n, *(x or hint(n)), a, b))


def cases():
    o = (0.0, 0.0, 0.0)
    unit = circle(o, Z, 1.0)
    e21 = ellipse(o, Z, 2.0, 1.0)
    h11 = hyperbola(o, Z, 1.0, 1.0)
    return [
        # Two lines.
        ('ll_cross', line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)), line((0.5, -1.0, 0.0), (0.5, 1.0, 0.0))),
        ('ll_parallel', line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)), line((-1.0, 1.0, 0.0), (3.0, 1.0, 0.0))),
        ('ll_skew', line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)), line((0.0, -1.0, 1.0), (0.0, 1.0, 1.0))),
        ('ll_coincident', line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)), line((3.0, 0.0, 0.0), (2.0, 0.0, 0.0))),
        # A line and a conic.
        ('lc_pierce', line((0.0, 1.0, -1.0), (0.0, 1.0, 1.0)), unit),
        ('lc_pierce_off', line((0.25, 0.5, -1.0), (0.25, 0.5, 1.0)), unit),
        ('lc_two', line((-2.0, 0.5, 0.0), (2.0, 0.5, 0.0)), unit),
        ('lc_tangent', line((-2.0, 1.0, 0.0), (2.0, 1.0, 0.0)), unit),
        ('lc_miss', line((-2.0, 2.0, 0.0), (2.0, 2.0, 0.0)), unit),
        ('lc_parallel', line((-2.0, 0.0, 1.0), (2.0, 0.0, 1.0)), unit),
        ('lc_tilted', line((-2.0, 0.25, 0.5), (2.0, 0.75, -0.25)), circle((0.0, 0.5, 0.0), TILT, 1.5)),
        ('le_tangent', line((2.0, -1.0, 0.0), (2.0, 1.0, 0.0)), e21),
        ('le_two', line((-3.0, 0.5, 0.0), (3.0, 0.25, 0.0)), e21),
        ('le_pierce', line((2.0, 0.0, -1.0), (2.0, 0.0, 1.0)), e21),
        ('lh_two', line((2.0, -3.0, 0.0), (2.0, 3.0, 0.0)), h11),
        ('lh_other', line((-2.0, -3.0, 0.0), (-2.0, 3.0, 0.0)), h11),
        ('lh_asymptote', line((0.0, 0.5, 0.0), (1.0, 1.5, 0.0)), h11),
        ('lh_tangent', line((1.0, -1.0, 0.0), (1.0, 1.0, 0.0)), h11),
        # Two conics in one plane.
        ('cc_two', unit, circle((1.0, 0.0, 0.0), Z, 1.0)),
        ('cc_outside', unit, circle((2.0, 0.0, 0.0), Z, 1.0)),
        ('cc_inside', unit, circle((0.5, 0.0, 0.0), Z, 0.5)),
        ('cc_coincident', unit, circle(o, (0.0, 0.0, -2.0), 1.0, Y)),
        ('cc_apart', unit, circle((3.0, 0.0, 0.0), Z, 1.0)),
        ('cc_concentric', unit, circle(o, Z, 0.5)),
        ('ce_four', circle(o, Z, 1.5), e21),
        ('ce_tangent', circle(o, Z, 2.0), e21),
        ('ee_four', e21, ellipse(o, Z, 2.0, 1.0, Y)),
        ('ee_coincident', e21, ellipse(o, Z, 2.0, 1.0, (-1.0, 0.0, 0.0))),
        ('eh_four', ellipse(o, Z, 2.0, 1.5), h11),
        ('hh_opposite', h11, hyperbola(o, Z, 1.0, 1.0, (-1.0, 0.0, 0.0))),
        ('hh_coincident', h11, hyperbola(o, (0.0, 0.0, -1.0), 1.0, 1.0)),
        ('ch_two', unit, hyperbola((-0.5, 0.0, 0.0), Z, 1.0, 1.0)),
        # Two conics in crossing planes.
        ('cc_across_two', unit, circle((0.5, 0.0, 0.5), X, 1.0)),
        ('cc_across_tangent', unit, circle((1.0, 0.0, 1.0), X, 1.0)),
        ('cc_across_one', circle(o, Z, 5.0), circle((3.0, 4.0, 1.0), X, 1.0)),
        ('cc_across_miss', unit, circle((0.5, 0.0, 2.0), X, 1.0)),
        ('cc_parallel_planes', unit, circle((0.0, 0.0, 1.0), Z, 1.0)),
        ('ce_across', ellipse(o, Z, 5.0, 2.5), circle((3.0, 2.0, 1.0), X, 1.0)),
        ('ce_across_miss', e21, circle((1.0, 0.0, 0.0), X, 1.0)),
        ('eh_across', e21, hyperbola((0.0, 0.0, 0.0), Y, 1.0, 0.5)),
        ('cc_tilted', circle((0.0, 0.5, 0.0), TILT, 1.5), circle((0.25, 0.0, 0.0), X, 1.25)),
    ]


def encode(name, a, b):
    rows = [f'case {name}']
    for kind, values in (a, b):
        rows.append(f'curve {kind} '+' '.join(repr(float(v)) for v in values))
    return '\n'.join(rows)


def generate():
    blocks, out = [], ['# case\trow (S7d.1, curve_curve_reference.py)']
    frames = ['# case\tcurve (0 the first, 1 the second)\taxis (n, x, y)\tstored unit vector (the reference\'s Frame3::new, as hex bits)']
    for name, a, b in cases():
        blocks.append(encode(name, a, b)+'\nend')
        for row in ref.rows(ref.Curve(*a), ref.Curve(*b)):
            out.append(f'{name}\t{ref.text(row)}')
        for k, (kind, values) in enumerate((a, b)):
            if kind == 'line':
                continue
            _, x, y, n = stored_axes(tuple(values[:9]))
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{name}\t{k}\t{key}\t'+' '.join(struct.pack('>d', c).hex() for c in v))
    return {'curve-curve-cases.txt': '\n'.join(blocks)+'\n',
            'curve-curve-expected.tsv': '\n'.join(out)+'\n',
            'curve-curve-frames.tsv': '\n'.join(frames)+'\n'}


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
