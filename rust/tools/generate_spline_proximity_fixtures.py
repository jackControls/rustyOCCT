#!/usr/bin/env python3
"""Complete closed-range spline minima from an independent exact oracle.

Cox basis coefficients, irreducible QQ factorization, VAS isolation/refinement,
and resultant images provide the expected answer. Interval separation can prove
a unique minimum without building a high-degree image polynomial. No Rust or
native OCCT answer supplies an expected parameter, distance, or completeness.
"""
import argparse
from fractions import Fraction as F
from functools import cmp_to_key
import json
from pathlib import Path

import sympy as sp

from exact_polynomial_oracle import X, Root, evaluate, poly, roots
from knot_editing_reference import polynomial, value

ROOT = Path(__file__).resolve().parents[1]
V = sp.Symbol('v')


def coefficients(p):
    return list(map(F, reversed(p.all_coeffs()))) if p else []


def rational_root(x):
    return Root(poly([-x, F(1)]), 1, (x, x))


def compare(a, b):
    if a.factor.monic() == b.factor.monic():
        lo, hi = max(a.low, b.low), min(a.high, b.high)
        if lo <= hi and a.factor.count_roots(lo, hi) == 1:
            return 0
    for _ in range(4096):
        if a.high < b.low:
            return -1
        if b.high < a.low:
            return 1
        a.refine()
        b.refine()
    raise ArithmeticError('distinct image roots failed to separate')


def image(root, n, w):
    p = root.factor.monic()
    nr, dr = n.rem(p), (w*w).rem(p)
    assert dr
    quotient, remainder = sp.div(nr, dr)
    if not remainder and quotient.degree() <= 0:
        return rational_root(F(quotient.nth(0)))
    resultant = sp.Poly(sp.resultant(p.as_expr(), nr.as_expr()-V*dr.as_expr(), X), V)
    factors = resultant.factor_list()[1]
    assert len(factors) == 1  # image of an irreducible field element
    im = sp.Poly.from_list(factors[0][0].monic().all_coeffs(), X, domain=sp.QQ)
    selected = []
    for interval, _ in im.intervals():
        lo, hi = map(F, interval)
        if root.sign_at(coefficients(n-lo*w*w)) >= 0 and root.sign_at(coefficients(n-hi*w*w)) <= 0:
            selected.append(Root(im, 1, interval))
    assert len(selected) == 1
    return selected[0]


def interval(p, a, b):
    lo = hi = F(0)
    for c in reversed(p):
        products = [lo*a, lo*b, hi*a, hi*b]
        lo, hi = min(products)+c, max(products)+c
    return lo, hi


def distance_range(candidate):
    cell, root = candidate['cell'], candidate['root']
    a, b = root.low, root.high
    low = high = F(0)
    for delta in cell['delta']:
        lo, hi = interval(delta, a, b)
        low += 0 if lo <= 0 <= hi else min(lo*lo, hi*hi)
        high += max(lo*lo, hi*hi)
    lo, hi = interval(coefficients(cell['w']), a, b)
    lo, hi = max(lo, cell['weights'][0]), min(hi, cell['weights'][1])
    assert 0 < lo <= hi
    return low/(hi*hi), high/(lo*lo)


def equations(case):
    degree, periodic = case['degree'], case['periodic']
    knots = tuple(map(F, case['knots']))
    controls = tuple((F(x)*F(w), F(y)*F(w), F(z)*F(w), F(w)) for x, y, z, w in case['controls'])
    curve = degree, periodic, knots, tuple(case['multiplicities']), controls
    query = list(map(F, case['query']))
    first, last = map(F, case['range'])
    if first == last:
        h = value(curve, first)
        # The singleton uses a constant exact Cartesian representation.
        h = [[h[i]/h[3]] for i in range(3)]+[[F(1)]]
        domains = [(first, first, h)]
        weights = (F(1), F(1))
    else:
        cuts = {first, last}
        if periodic:
            period = knots[-1]-knots[0]
            for turn in range(int((first-knots[0])//period)-1, int((last-knots[0])//period)+2):
                cuts.update(k+turn*period for k in knots if first < k+turn*period < last)
        else:
            cuts.update(k for k in knots if first < k < last)
        cuts = sorted(cuts)
        domains = [(a, b, polynomial(curve, a, b)) for a, b in zip(cuts, cuts[1:])]
        weights = min(c[3] for c in controls), max(c[3] for c in controls)
    result = []
    for lo, hi, h in domains:
        h = list(map(poly, h))
        w = h[3]
        delta = [h[i]-query[i]*w for i in range(3)]
        n = sum((p*p for p in delta), poly([]))
        f = n.diff()*w-2*n*w.diff()
        assert f.is_zero or f.degree() <= 3*degree-2
        result.append(dict(lo=lo, hi=hi, h=h, w=w, n=n, f=f,
                           delta=list(map(coefficients, delta)), weights=weights))
    return result


def minima(cells):
    candidates = []
    for cell in cells:
        def add(root, span=False):
            candidates.append(dict(cell=cell, root=root, span=span))
        if cell['lo'] == cell['hi']:
            add(rational_root(F(0)))
        elif cell['f'].is_zero:
            add(rational_root(F(0)), True)
        else:
            add(rational_root(F(0)))
            add(rational_root(F(1)))
            for r in roots(coefficients(cell['f']), F(0), F(1)):
                if r.compare(0) and r.compare(1):
                    add(r)
    zeros = [c for c in candidates if not c['cell']['n'].rem(c['root'].factor)]
    if zeros:
        winners = zeros
    else:
        # Discard a candidate only when its lower bound is strictly greater
        # than an attained candidate's upper bound. Equal minima cannot vanish.
        active = candidates.copy()
        for _ in range(64):
            bounds = [distance_range(c) for c in active]
            upper = min(b for _, b in bounds)
            active = [c for c, (a, _) in zip(active, bounds) if a <= upper]
            if len(active) == 1:
                break
            if all(c['root'].low == c['root'].high for c in active):
                break
            # Low-degree unresolved ties are cheap to prove with resultants.
            if _ >= 16 and max(c['root'].factor.degree() for c in active) <= 12:
                break
            for c in active:
                c['root'].refine()
        if len(active) == 1:
            winners = active
        else:
            for c in active:
                c['distance'] = image(c['root'], c['cell']['n'], c['cell']['w'])
            best = min((c['distance'] for c in active), key=cmp_to_key(compare))
            winners = [c for c in active if compare(c['distance'], best) == 0]
    spans = []
    for c in sorted((c for c in winners if c['span']), key=lambda c: c['cell']['lo']):
        lo, hi = c['cell']['lo'], c['cell']['hi']
        if spans and spans[-1][1] == lo:
            spans[-1][1] = hi
        else:
            spans.append([lo, hi])
    points = []
    for c in winners:
        if c['span']:
            continue
        cell, r = c['cell'], c['root']
        a, length = cell['lo'], cell['hi']-cell['lo']
        if r.low == r.high:
            parameter = a+length*r.low
            if any(lo <= parameter <= hi for lo, hi in spans):
                continue
            encoded = (F(0), F(1), [-parameter, F(1)], parameter, parameter)
        else:
            encoded = (a, a+length, coefficients(r.factor.monic()), r.low, r.high)
        if encoded not in points:
            points.append(encoded)
    points.sort(key=lambda item: item[0]+(item[1]-item[0])*item[3])
    return points, spans


def generate():
    cases = json.loads((ROOT/'fixtures/spline-proximity-inputs.json').read_text())['cases']
    rows, inputs = [], []
    for c in cases:
        words = [c['name'], c['degree'], int(c['periodic']), len(c['controls']), len(c['knots']), *c['range'], *c['query']]
        words.extend(x for p in c['controls'] for x in p)
        words.extend(x for pair in zip(c['knots'], c['multiplicities']) for x in pair)
        inputs.append(' '.join(format(x, '.17g') if isinstance(x, float) else str(x) for x in words))
        cells = equations(c)
        points, spans = minima(cells)
        row = [c['name'], len(cells)]

        def polynomial_words(coeffs):
            row.extend([len(coeffs), *coeffs])

        for cell in cells:
            row.extend([cell['lo'], cell['hi']])
            for h in cell['h']:
                polynomial_words(coefficients(h))
        row.extend([len(points), len(spans)])
        for a, b, p, lo, hi in points:
            row.extend([a, b, lo, hi])
            polynomial_words(p)
        for span in spans:
            row.extend(span)
        rows.append(' '.join(map(str, row)))
        print(f"{c['name']}: {len(points)} isolated minima, {len(spans)} minimum intervals", flush=True)
    return {'spline-proximity-inputs.txt': '\n'.join(inputs)+'\n',
            'spline-proximity.tsv': '# name cells [first last Hx Hy Hz W] points intervals [cell_first cell_last root_lo root_hi P] [first last]\n'+'\n'.join(rows)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for name, contents in generate().items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)


if __name__ == '__main__':
    main()
