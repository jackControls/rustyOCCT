#!/usr/bin/env python3
"""Exact rational intersection inputs and complete closed-form certificates.

Python Fraction manufacturing and Bernstein-to-power verification are independent
of Rust's de Boor polynomials and algebraic root isolation. No Rust output is read.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import math
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[1]


def atom(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'


def enclosure(x):
    try:
        value = float(x)
    except OverflowError:
        return 'X'
    if not math.isfinite(value):
        return 'X'
    lo = math.nextafter(value, -math.inf) if F(value) > x else value
    hi = math.nextafter(value, math.inf) if F(value) < x else value
    if not all(map(math.isfinite, (lo, hi))):
        return 'X'
    return ':'.join(struct.pack('>d', v).hex() for v in (lo, hi))


def multiply(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x*y
    return c


def bernstein_power(values):
    n = len(values)-1
    return [sum(values[i]*math.comb(n, i)*math.comb(n-i, k-i)*(-1)**(k-i)
                for i in range(k+1)) for k in range(n+1)]


def point(t, coordinates, order, first, last, crossing):
    values = [t, *coordinates]
    contact = 'B' if t in (first, last) else 'C' if crossing else 'T'
    return [*map(atom, values), *map(enclosure, values), contact,
            str(0 if t == first else order), str(0 if t == last else order)]


def row(name, family, degree, periodic, controls, knots, mults, first, last, points, overlaps):
    inputs = [name, str(family), str(degree), str(int(periodic)), str(len(controls)),
              str(len(knots)), atom(first), atom(last)]
    inputs += [atom(x) for control in controls for x in control]
    inputs += [x for k, m in zip(knots, mults) for x in (atom(k), str(m))]
    outputs = [str(len(points)), str(len(overlaps))]
    outputs += [x for p in points for x in p]
    outputs += [x for a, b in overlaps for x in (atom(a), atom(b), enclosure(a), enclosure(b))]
    return ' '.join(inputs) + ' | ' + ' '.join(outputs)


def factored(name, roots, family, a, width, lo, hi):
    n = len(roots)
    numerator = [F(1)]
    for root in roots:
        numerator = multiply(numerator, [-root, F(1)])
    values = [sum(numerator[j]*F(math.comb(i, j), math.comb(n, j)) for j in range(i+1))
              for i in range(n+1)]
    assert bernstein_power(values) == numerator
    controls = []
    for i, value in enumerate(values):
        t = F(i, n)
        w = 1+t
        controls.append([t, F(0), value, w] if family == 0 else
                        [value, w, F(0) if family == 1 else t, w])
    polynomials = [bernstein_power([c[i] for c in controls]) for i in range(4)]
    assert polynomials[3] == [F(1), F(1)] + [F(0)]*(n-1)
    if family:
        implicit = [F(0)]*(2*n+1)
        for c in range(4 if family == 1 else 2):
            for i, value in enumerate(multiply(polynomials[c], polynomials[c])):
                implicit[i] += (-1 if c == 3 else 1)*value
        if family == 2:
            for i, value in enumerate(multiply(polynomials[3], polynomials[3])):
                implicit[i] -= value
        assert implicit == multiply(numerator, numerator)
    first, last = a+width*lo, a+width*hi
    points = []
    for t, m in sorted(Counter(roots).items()):
        if not lo <= t <= hi:
            continue
        coordinates = [t/(1+t), F(0), F(0)] if family == 0 else \
                      [F(0), F(1), F(0) if family == 1 else t/(1+t)]
        order = m if family == 0 else 2*m
        points.append(point(a+width*t, coordinates, order, first, last, order % 2 == 1))
    return row(name, family, n, False, controls, [a, a+width], [n+1]*2,
               first, last, points, [])


def generate():
    rows = ['# rational input | exact contacts/overlaps and independently rounded bounds; X means unrepresentable']
    domains = [(F(1,3), F(5,7)), (F(2**2048), F(3)),
               (-F(2**2048), F(1,2**2048)), (-F(1,2**2048), F(1,2**2047))]
    candidates = [F(-1,5), F(0), F(1,7), F(1,3), F(2,3), F(1), F(6,5)]
    for n in [1,2,3,5,8,13,25]:
        roots = [candidates[(i+n) % len(candidates)] for i in range(n)]
        for family in range(3):
            for di, (a, width) in enumerate(domains):
                for ti, (lo, hi) in enumerate([(F(0),F(1)), (F(1,7),F(2,3)), (F(2,7),F(5,7))]):
                    rows.append(factored(f'd{n}_f{family}_d{di}_t{ti}', roots, family, a, width, lo, hi))
    for family in range(3):
        rows.append(factored(f'near_roots_{family}', [F(1,3),F(1,3)+F(1,2**2048)],
                             family, F(0), F(1), F(0), F(1)))
        rows.append(factored(f'mixed_range_{family}', [F(1,2**2048),F(1,2)],
                             family, F(0), F(2**2048), F(0), F(1)))
    for i, (a, width) in enumerate(domains):
        for family in range(3):
            controls = [[F(0),F(0),F(0),F(1)]]*3 if family == 0 else \
                       [[F(1),F(0),F(0),F(1)],[F(1),F(1),F(0),F(1)],[F(0),F(2),F(0),F(2)]]
            first, last = a+width/7, a+width*F(5,6)
            rows.append(row(f'overlap_{family}_{i}', family, 2, False, controls,
                            [a,a+width], [3,3], first, last, [], [(first,last)]))
        # A cyclic piecewise line has Z=0 exactly at a+k*3*width/2.
        first, last = a-width*F(5,2), a+width*F(13,2)
        roots = [a+width*F(3*k,2) for k in range(-3,6)]
        points = [point(t,[F(0)]*3,1,first,last,True) for t in roots if first <= t <= last]
        controls = [[F(0),F(0),F(0),F(1)],[F(1),F(0),F(1),F(1)],[-F(1),F(0),-F(1),F(1)]]
        rows.append(row(f'periodic_{i}',0,1,True,controls,[a+j*width for j in range(4)],
                        [1]*4,first,last,points,[]))
    return '\n'.join(rows)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = ROOT/'fixtures/exact-spline-intersections.tsv'
    text = generate()
    if args.check:
        if path.read_text() != text:
            parser.error('exact intersection fixtures changed; investigate before updating')
    else:
        path.write_text(text)
    print(f'{path.name}: {len(text.splitlines())-1} independent exact cases, {len(text)} bytes')


if __name__ == '__main__':
    main()
