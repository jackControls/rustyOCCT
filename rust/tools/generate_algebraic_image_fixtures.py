#!/usr/bin/env python3
"""Independent rational-image equations using SymPy's exact resultants.

The production kernel uses quotient-ring powers and integer row elimination.
This oracle eliminates the parameter with a resultant, then removes repeated
factors. Every coefficient is checked, including relative numerator/denominator
scales; no production output or duplicate of its construction is used.
"""
import argparse
from pathlib import Path
import random

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
x, y = sp.symbols('x y')


def generate():
    rows = []

    def add(name, p, n, d):
        p, n, d = (sp.Poly(v, x, domain=sp.QQ) for v in (p, n, d))
        p = p.sqf_part()
        kept = p.exquo(sp.gcd(p, d))
        assert kept.degree() >= 1
        image = sp.Poly(sp.resultant(kept.as_expr(), n.as_expr()-y*d.as_expr(), x), y)
        image = image.sqf_part().monic()
        row = [name]
        for poly in (p, n, d, image):
            coefficients = list(reversed(poly.all_coeffs())) if poly else []
            row += [str(len(coefficients)), *map(str, coefficients)]
        rows.append(' '.join(row))

    add('identity', x*x-2, x, 1)
    add('constant_image', x*x-2, x*x, 1)
    add('different_rational_images', (x*x-2)*(x*x-3), x*x, 1)
    add('complex_poles', (x*x+1)*(x-2), x*(x*x+1), x*x+1)
    add('fractional_image', x**3-2, x+1, x-1)
    add('degree25_collapse', x**25-1, x**5, 1)
    add('zero_image', x*x-2, 0, 7)
    add('negative_denominator', x**4-3, x**3+x, -7)
    add('relative_fractional_scales', x**3-2, (x+1)/3, (x-1)/7)
    add('scaled_constant_image', x*x-2, sp.Rational(1, 2)*x*x, sp.Rational(1, 3))
    add('real_pole_factor', (x-1)*(x*x-2), x+2, x-1)
    add('mixed_pole_factors', (x-2)*(x+2)*(x*x+1)*(x*x-3), x**3+1, (x-2)*(x*x+1))
    add('linear_input', 3*x-2, x*x+7, 5*x-1)
    add('high_numerator', x**3-2, x**17+sp.Rational(2, 7)*x-1, x*x+sp.Rational(3, 11))
    rng = random.Random(7433)
    for i in range(80):
        degree = 2+i % 6

        def polynomial(d):
            return sum(rng.randrange(-5, 6)*x**j for j in range(d))+x**d

        add(f'random_{i}', polynomial(degree), polynomial(degree+1), polynomial(degree-1))
    return '# name [count ascending_rational_coefficients] for P,N,D,image\n'+'\n'.join(rows)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = generate()
    path = ROOT/'fixtures/algebraic-images.tsv'
    if args.check:
        if path.read_text() != text:
            parser.error('algebraic image fixtures changed; investigate before updating')
    else:
        path.write_text(text)
    print(f'algebraic-images.tsv: {len(text.splitlines())-1} independent whole equations')


if __name__ == '__main__':
    main()
