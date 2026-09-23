#!/usr/bin/env python3
"""Independent algebraic ordering by irreducible factor identity and VAS bounds.

Production compares signs in an isolating interval using Sturm-Tarski queries.
Here irreducible factors plus canonical real-root indices decide equality;
continued-fraction intervals separate distinct values. No Rust output is read.
"""
import argparse
from fractions import Fraction as F
from functools import cmp_to_key
from pathlib import Path
import random

from exact_polynomial_oracle import Root, multiply, poly
from generate_spatial_fixtures import bits, value

ROOT = Path(__file__).resolve().parents[1]


def compare(a, b):
    if a.identity == b.identity:
        return 0
    for _ in range(4096):
        if a.high < b.low:
            return -1
        if b.high < a.low:
            return 1
        a.refine()
        b.refine()
    raise ArithmeticError('distinct irreducible-root identities did not separate')


def roots(coefficients, domain):
    result = []
    for factor, multiplicity in poly(coefficients).factor_list()[1]:
        factor = factor.monic()
        for index, (interval, _) in enumerate(factor.intervals()):
            root = Root(factor, multiplicity, interval)
            root.identity = (tuple(factor.all_coeffs()), index)
            if domain is not None and (root.compare(domain[0]) < 0 or root.compare(domain[1]) > 0):
                continue
            result.append(root)
    return sorted(result, key=cmp_to_key(compare))


def generate():
    source = []
    for row in (ROOT/'fixtures/real-roots.tsv').read_text().splitlines():
        if row.startswith('#'):
            continue
        name, np, _, lower, upper, *words = row.split()
        p = [value(int(word, 16)) for word in words[:int(np)]]
        if not p or not any(p):
            continue
        domain = None if lower == '*' else (value(int(lower, 16)), value(int(upper, 16)))
        source.append((name, p, domain))
    pairs = []
    for i, a in enumerate(source):
        pairs.append((a, source[(i+1) % len(source)]))
    for degree in (9, 25):
        cluster = next(a for a in source if a[0] == f'cluster_{degree}')
        # Same roots with opposite leading coefficient and different square-free
        # scaling; the tight clusters share binary64 enclosures.
        pairs.append((cluster, (f'negative_cluster_{degree}', [-c for c in cluster[1]], cluster[2])))
    square = ('square', [-2., 0., 1.], None)
    for exponent in (-70, -500, -1074):
        pairs.append((square, (f'close_{exponent}', [-2., -2.**exponent, 1.], None)))
    pairs.append((('linear_huge', [-value(0x7fefffffffffffff), value(1)], None),
                  ('linear_half_huge', [-value(0x7fefffffffffffff), value(2)], None)))
    rng = random.Random(0x93_47_FC)
    for i in range(30):
        common = [-2, 0, 1] if i % 2 else [-3, 0, 0, 1]
        if i % 3 == 0:
            common = multiply(common, common)
        def side(label):
            extra = [F(rng.randrange(-4, 5)) for _ in range(1+i % 7)]+[F(1)]
            p = multiply(common, extra)
            assert all(F(float(c)) == c for c in p)
            return (f'shared_{i}_{label}', list(map(float, p)), None)
        pairs.append((side('left'), side('right')))
    rows = []
    for a, b in pairs:
        ar, br = roots(a[1], a[2]), roots(b[1], b[2])
        row = [a[0]+'__'+b[0]]
        for (_, p, domain), found in ((a, ar), (b, br)):
            row += [str(len(p)), *(map(bits, domain) if domain else ('*', '*')),
                    *map(bits, p), str(len(found))]
        row += [str(compare(x, y)) for x in ar for y in br]
        rows.append(' '.join(row))
    return '# name [n lower upper coefficients root_count]x2 row-major-exact-comparisons\n'+'\n'.join(rows)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = generate()
    path = ROOT/'fixtures/root-comparisons.tsv'
    if args.check:
        if path.read_text() != text:
            parser.error('root comparisons changed; investigate before updating')
    else:
        path.write_text(text)
    print(f'root-comparisons.tsv: {len(text.splitlines())-1} independent equation pairs')


if __name__ == '__main__':
    main()
