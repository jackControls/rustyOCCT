#!/usr/bin/env python3
"""Reproduce exact proximity expectations with independent analytic formulas."""
import argparse
from fractions import Fraction as Q
import math
from pathlib import Path
import random
import struct
from compare_proximity import cases
from proximity_reference import COUNTS, bits, distance2, enclosure, exact, parse, validate

OUTPUT = Path(__file__).resolve().parents[1]/'fixtures/proximity.tsv'


def inputs():
    for row in cases().splitlines():
        name, shapes, rest = parse(row)
        assert not rest
        yield name, shapes
    rng = random.Random(0xD157A3CE)
    for a in COUNTS:
        for b in COUNTS:
            for mode in range(6):
                def number():
                    if mode == 5:
                        # Independent exponents across coordinates, including
                        # subnormals and near-maximum finite values.
                        word = rng.getrandbits(64)
                        if word & 0x7ff0000000000000 == 0x7ff0000000000000:
                            word ^= 0x0010000000000000
                        return struct.unpack('>d', word.to_bytes(8, 'big'))[0]
                    exponent = [0, -500, 450, -1070, 950][mode]
                    return math.ldexp(float(rng.randrange(-32, 33)), exponent)
                shapes = [(k, [tuple(number() for _ in range(3)) for _ in range(COUNTS[k])]) for k in (a, b)]
                yield f'random_{a}{b}_{mode}', shapes
    tiny, maximum = math.ulp(0.), float.fromhex('0x1.fffffffffffffp+1023')
    point = lambda p: ('P', [p])
    origin = point((0., 0., 0.))
    tri = ('T', [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)])
    yield 'distance_overflow', [origin, point((maximum, maximum, maximum))]
    yield 'squared_overflow', [origin, point((maximum, 0., 0.))]
    yield 'subnormal_distance', [point((tiny, 0., 0.)), ('L', [(0., 0., 0.), (tiny, tiny, 0.)])]
    yield 'parameter_overflow', [point((1., 0., 0.)), ('L', [(0., 0., 0.), (tiny, 0., 0.)])]
    yield 'point_overflow', [('L', [(0., 0., 0.), (1., tiny, 0.)]), ('L', [(0., 1., 0.), (1., 1., 0.)])]
    yield 'huge_difference', [('S', [(-maximum, 0., 0.), (maximum, 0., 0.)]), point((0., tiny, 0.))]
    yield 'subnormal_triangle', [('T', [(0., 0., 0.), (tiny, 0., 0.), (0., tiny, 0.)]), point((tiny, tiny, tiny))]
    for delta in [-tiny, 0., tiny]:
        yield f'near_parallel_{bits(delta)}', [('L', [(0., 0., 0.), (1., delta, 0.)]), ('L', [(0., 1., 1.), (1., 1., 1.)])]
        yield f'triangle_edge_{bits(delta)}', [tri, point((.5, delta, 1.))]
    for k in ['L', 'S', 'T', 'F']:
        yield f'collapsed_{k}', [(k, [(0., 0., 0.)]*COUNTS[k]), origin]
    for k in ['F', 'T']:
        yield f'collinear_{k}', [(k, [(0., 0., 0.), (1., 1., 1.), (2., 2., 2.)]), origin]
    for k in COUNTS:
        for value in [math.nan, math.inf, -math.inf]:
            p = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)][:COUNTS[k]]
            p[-1] = (p[-1][0], p[-1][1], value)
            yield f'nonfinite_{k}_{bits(value)}', [(k, p), origin]


def generate():
    rows = ['# Input binary64 bits; analytic Fraction distance squared; minimal squared/distance bounds.',
            '# R = exact result; D = degenerate line/plane/triangle; N = nonfinite; U = enclosure overflow.']
    for name, shapes in inputs():
        status = next((s for s in map(validate, shapes) if s != 'R'), 'R')
        result = [status]
        if status == 'R':
            d = distance2(*map(exact, shapes))
            result += [str(d.numerator)+'/'+str(d.denominator), *enclosure(d), *enclosure(d, sqrt=True)]
        row = [name]
        for k, points in shapes:
            row += [k, *[bits(x) for p in points for x in p]]
        rows.append(' '.join([*row, *result]))
    return '\n'.join(rows)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    content = generate()
    if args.check:
        if OUTPUT.read_text() != content: raise SystemExit('proximity fixtures differ; investigate before updating')
    else:
        OUTPUT.write_text(content)
    print(f'{len(content.splitlines())-2} independent proximity fixtures verified')


if __name__ == '__main__': main()
