#!/usr/bin/env python3
"""Complete spline/line and spline/segment preimages from an independent oracle.

Cox basis power coefficients give each span's homogeneous polynomials. A
common QQ gcd of the cross-product equations, irreducible factorization and
VAS isolation give isolated parameters. Contained spans are clipped by exact
rational SAMPLE POINTS between separated boundary roots; production instead
uses first-nonzero Taylor derivatives at each boundary. No Rust or native
OCCT answer supplies an expected parameter, interval or completeness claim.

The first 31 cases are the binary64 inputs captured from native OCCT before a
Rust port existed. The remaining cases use arbitrary rationals.
"""
import argparse
from fractions import Fraction as F
from functools import cmp_to_key
from pathlib import Path
import sys

import sympy as sp

from exact_polynomial_oracle import X, Root, evaluate, poly, roots
from generate_spline_fixtures import axis
from knot_editing_reference import polynomial, value

ROOT = Path(__file__).resolve().parents[1]
sys.setrecursionlimit(10000)


def coefficients(p):
    return list(map(F, reversed(p.all_coeffs()))) if not p.is_zero else []


def rational_root(x):
    return Root(poly([-x, F(1)]), 1, (x, x))


def power(expressions, domain=(-2, 2), degree=None, rounded=True):
    """Exact Bernstein controls of polynomial Cartesian components."""
    a, b = map(sp.Rational, domain)
    ps = [sp.Poly(sp.expand(sp.sympify(e).subs(X, a+(b-a)*X)), X) for e in expressions]
    degree = degree or max(p.degree() for p in ps)
    controls = []
    for i in range(degree+1):
        c = [F(str(sum(p.nth(j)*sp.binomial(i, j)/sp.binomial(degree, j) for j in range(i+1)))) for p in ps]
        controls.append([F(float(x)) if rounded else x for x in c]+[F(1)])
    return dict(degree=degree, periodic=False, controls=controls,
                knots=[F(a), F(b)], multiplicities=[degree+1]*2, range=[F(a), F(b)])


def make_cases():
    result = []

    def add(name, curve, endpoints, kind='S', domain=None):
        c = dict(curve)
        c.update(name=name, kind=kind, native=len(result) < 31,
                 endpoints=[[F(x) for x in p] for p in endpoints])
        if domain is not None:
            c['range'] = [F(x) for x in domain]
        result.append(c)

    # Native preflight inputs, in their original order.
    parabola = power([X, X**2, 0])
    add('parabola_secant', parabola, [[-2, 1, 0], [2, 1, 0]])
    add('parabola_secant_reversed', parabola, [[2, 1, 0], [-2, 1, 0]])
    add('parabola_one_hit', parabola, [[0, 1, 0], [2, 1, 0]])
    add('parabola_tangent', parabola, [[-1, 0, 0], [1, 0, 0]])
    add('parabola_skew', parabola, [[-2, 1, 1], [2, 1, 1]])
    add('parabola_outside_segment', parabola, [[2, 1, 0], [3, 1, 0]])
    add('parabola_line', parabola, [[2, 1, 0], [3, 1, 0]], kind='L')
    add('parabola_trimmed', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(0, 2))
    add('parabola_closed_ends', parabola, [[-1, 1, 0], [1, 1, 0]], domain=(-1, 1))
    for exponent in [-40, -20]:
        for sign in [-1, 1]:
            c = power([X, X**2+sign*sp.Rational(2)**exponent, 0], (-1, 1))
            add(f'near_tangent_{exponent}_{sign}', c, [[-1, 0, 0], [1, 0, 0]])
    # Degree-three Bernstein rounding makes the ideal tangent an exact miss.
    add('rounded_space_near_tangent', power([X, X**2, X**3]), [[-1, 0, 0], [1, 0, 0]])
    space = power([3*X, 3*X**2, X**3])
    add('space_tangent', space, [[-1, 0, 0], [1, 0, 0]])
    add('space_endpoint_pair', space, [[-3, 3, -1], [3, 3, 1]])
    add('space_skew', space, [[-3, 3, 0], [3, 3, 0]])
    backwards = power([X**2, 0, 0])
    add('collinear_disconnected_parameters', backwards, [[1, 0, 0], [2, 0, 0]])
    add('collinear_whole_line', backwards, [[1, 0, 0], [2, 0, 0]], kind='L')
    add('collinear_connected_parameters', backwards, [[0, 0, 0], [1, 0, 0]])
    add('collinear_touch_only', backwards, [[-1, 0, 0], [0, 0, 0]])
    circle = dict(degree=2, periodic=False,
                  controls=[[F(1), F(0), F(0), F(1)], [F(1), F(1), F(0), F(1)], [F(0), F(1), F(0), F(2)]],
                  knots=[F(0), F(1)], multiplicities=[3, 3], range=[F(0), F(1)])
    add('rational_circle_diagonal', circle, [[0, 0, 0], [1, 1, 0]])
    add('rational_circle_endpoints', circle, [[1, 0, 0], [0, 1, 0]])
    add('rational_circle_tangent', circle, [[1, -1, 0], [1, 1, 0]])
    polyline = dict(degree=1, periodic=False,
                    controls=[[F(-1), F(1), F(0), F(1)], [F(0), F(0), F(0), F(1)],
                              [F(1), F(0), F(0), F(1)], [F(2), F(1), F(0), F(1)]],
                    knots=[F(0), F(1), F(2), F(3)], multiplicities=[2, 1, 1, 2], range=[F(0), F(3)])
    add('polyline_partial_overlap', polyline, [[F(1, 2), 0, 0], [2, 0, 0]])
    add('polyline_knot_hit', polyline, [[0, -1, 0], [0, 1, 0]])
    periodic = dict(degree=1, periodic=True,
                    controls=[[F(0), F(0), F(0), F(1)], [F(1), F(0), F(0), F(1)], [F(0), F(1), F(0), F(1)]],
                    knots=[F(0), F(1), F(2), F(3)], multiplicities=[1]*4, range=[F(0), F(3)])
    add('periodic_one_turn', periodic, [[-1, 0, 0], [2, 0, 0]])
    add('periodic_multiple_turns', periodic, [[-1, 0, 0], [2, 0, 0]], domain=(-3, 6))
    for degree in [3, 8, 25]:
        add(f'chebyshev_{degree}', power([0, sp.chebyshevt(degree, X), 0], (-1, 1)), [[-1, 0, 0], [1, 0, 0]])
    assert len(result) == 31

    # Exact rational trims, including closed endpoint contacts.
    add('rational_trim_hit', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(F(-7, 5), F(1, 3)))
    add('rational_trim_endpoint', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(-1, F(1, 3)))
    add('rational_trim_miss', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(F(-1, 3), F(5, 7)))
    # Singletons and collapsed segments.
    add('singleton_hit', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(1, 1))
    add('singleton_miss', parabola, [[-2, 1, 0], [2, 1, 0]], domain=(F(1, 2), F(1, 2)))
    add('singleton_knot_point', polyline, [[0, 0, 0], [0, 0, 0]], domain=(1, 1))
    add('collapsed_on_curve', power([X, X**2, 0], rounded=False), [[F(1, 3), F(1, 9), 0]]*2)
    add('collapsed_off_curve', parabola, [[0, 1, 0]]*2)
    add('collapsed_periodic_aliases', periodic, [[0, 0, 0]]*2, domain=(-3, 6))
    add('collapsed_backtracking_pair', backwards, [[1, 0, 0]]*2)
    constant = dict(degree=2, periodic=False,
                    controls=[[F(1, 3), F(2), F(-4), F(w)] for w in (1, 3, 2)],
                    knots=[F(1, 7), F(10, 7)], multiplicities=[3, 3], range=[F(1, 7), F(10, 7)])
    add('constant_collapsed_segment', constant, [[F(1, 3), 2, -4]]*2)
    add('constant_on_line', constant, [[0, 2, -4], [1, 2, -4]], kind='L')
    add('constant_inside_segment', constant, [[0, 2, -4], [1, 2, -4]])
    add('constant_outside_segment', constant, [[0, 2, -4], [F(1, 4), 2, -4]])
    # Oblique 3D and irrational contacts.
    twisted = power([3*X, 3*X**2, X**3], rounded=False)
    add('twisted_cubic_oblique_secant', twisted, [[-3, 3, -1], [6, 12, 8]], kind='L')
    add('twisted_cubic_oblique_segment', twisted, [[-3, 3, -1], [F(3, 2), F(3, 4), F(1, 8)]])
    add('circle_irrational_slope', circle, [[0, 0, 0], [1, 2, 0]], kind='L')
    add('circle_irrational_segment', circle, [[0, 0, 0], [F(2, 3), F(4, 3), 0]])
    # Exact parameters outside binary64 range and sub-float separation.
    huge = F(2)**1100
    shifted = dict(parabola)
    shifted.update(knots=[huge-2, huge+2], range=[huge-2, huge+2])
    add('huge_parameter_domain', shifted, [[-2, 1, 0], [2, 1, 0]])
    tiny = F(1, 2**1100)
    small = dict(parabola)
    small.update(controls=[[x*tiny, y*tiny, z*tiny, w] for x, y, z, w in parabola['controls']])
    add('tiny_coordinates', small, [[-2*tiny, tiny, 0], [2*tiny, tiny, 0]])
    epsilon = F(1, 2**180)
    gap = dict(degree=2, periodic=False,
               controls=[[x, F(0), F(0), F(1)] for x in (1+epsilon, F(-1), 1-epsilon)],
               knots=[F(0), F(2)], multiplicities=[3, 3], range=[F(0), F(2)])
    add('sub_float_excluded_gap', gap, [[0, 0, 0], [2, 0, 0]])
    add('sub_float_point_pair', gap, [[0, 0, 0]]*2)
    # Contact multiplicity six at a segment endpoint.
    add('high_order_endpoint_touch', power([X, X**6, 0], (-1, 1), rounded=False), [[0, 0, 0], [1, 0, 0]])
    # Unclamped quadratic whose active domain is [0,2].
    unclamped = dict(degree=2, periodic=False,
                     controls=[[F(-3), F(0), F(0), F(1)], [F(-1), F(0), F(0), F(2)],
                               [F(1), F(0), F(0), F(3)], [F(3), F(0), F(0), F(4)]],
                     knots=[F(k) for k in range(-2, 5)], multiplicities=[1]*7, range=[F(0), F(2)])
    add('unclamped_on_line', unclamped, [[0, 0, 0], [1, 0, 0]], kind='L')
    add('unclamped_clipped', unclamped, [[F(-1, 3), 0, 0], [F(1, 2), 0, 0]])
    # Periodic quadratic crossing seams with rational trims.
    square = dict(degree=2, periodic=True,
                  controls=[[F(1), F(0), F(0), F(1)], [F(0), F(1), F(0), F(1)],
                            [F(-1), F(0), F(0), F(2)], [F(0), F(-1), F(0), F(1)]],
                  knots=[F(k) for k in range(5)], multiplicities=[1]*5, range=[F(0), F(4)])
    add('periodic_quadratic_turns', square, [[-2, F(1, 5), 0], [2, F(1, 5), 0]], domain=(F(-9, 2), F(17, 3)))
    add('periodic_quadratic_skew', square, [[0, 0, 1], [0, 1, 1]], domain=(-4, 8))
    # Overlaps continuing across interior knots and a rational-weight retrace.
    straight = dict(degree=1, periodic=False,
                    controls=[[F(x), F(0), F(0), F(w)] for x, w in ((0, 1), (1, 3), (3, 2), (2, 1), (5, 1))],
                    knots=[F(k) for k in range(5)], multiplicities=[2, 1, 1, 1, 2], range=[F(0), F(4)])
    add('overlap_across_knots', straight, [[F(1, 2), 0, 0], [F(5, 2), 0, 0]])
    weighted = dict(degree=2, periodic=False,
                    controls=[[F(4), F(0), F(0), F(1)], [F(-4), F(0), F(0), F(3)], [F(4), F(0), F(0), F(2)]],
                    knots=[F(-2), F(2)], multiplicities=[3, 3], range=[F(-2), F(2)])
    add('weighted_retrace_segment', weighted, [[F(1, 3), 0, 0], [F(3, 2), 0, 0]])
    return result


def homogeneous(case):
    return tuple(tuple([x*w, y*w, z*w, w]) for x, y, z, w in case['controls'])


def curve_tuple(case):
    return (case['degree'], case['periodic'], tuple(case['knots']),
            tuple(case['multiplicities']), homogeneous(case))


def cells(case):
    """Common span partition of the closed query, in original parameter units."""
    curve = curve_tuple(case)
    first, last = case['range']
    if first == last:
        h = value(curve, first)
        return [(first, first, [[h[i]/h[3]] for i in range(3)]+[[F(1)]])]
    flat, a, b = axis(case['degree'], case['knots'], case['multiplicities'], case['periodic'])
    cuts = {first, last}
    if case['periodic']:
        period = b-a
        for turn in range(int((first-a)//period)-1, int((last-a)//period)+2):
            cuts.update(k+turn*period for k in case['knots'] if first < k+turn*period < last)
    else:
        assert a <= first and last <= b
        cuts.update(k for k in flat if first < k < last)
    cuts = sorted(cuts)
    return [(lo, hi, polynomial(curve, lo, hi)) for lo, hi in zip(cuts, cuts[1:])]


def separated(items):
    def compare(a, b):
        for _ in range(4096):
            if a.high < b.low:
                return -1
            if b.high < a.low:
                return 1
            a.refine()
            b.refine()
        raise ArithmeticError('boundary roots failed to separate')
    items = sorted(items, key=cmp_to_key(compare))
    for a, b in zip(items, items[1:]):
        assert compare(a, b) < 0
    return items


def sign(x):
    return (x > 0)-(x < 0)


def local_sets(h, a, d, kind):
    """Local points and closed intervals of one cell, variable x in [0,1]."""
    w = poly(h[3])
    delta = [poly(h[i])-a[i]*w for i in range(3)]
    norm = sum(v*v for v in d)
    if norm:
        equations = [delta[(i+1) % 3]*d[(i+2) % 3]-delta[(i+2) % 3]*d[(i+1) % 3] for i in range(3)]
    else:
        equations = delta
    common = equations[0].gcd(equations[1]).gcd(equations[2])
    along = sum((delta[i]*d[i] for i in range(3)), poly([]))
    end = along-norm*w
    # Non-strict inequalities: along >= 0 and along - norm*W <= 0.
    constraints = [(coefficients(along), 1), (coefficients(end), -1)] if kind == 'S' and norm else []

    def allowed_root(r):
        return all(s*r.sign_at(p) >= 0 for p, s in constraints)

    def allowed_value(x):
        return all(s*sign(evaluate(p, x)) >= 0 for p, s in constraints)

    if not common.is_zero:
        found = roots(coefficients(common), F(0), F(1)) or []
        return [r for r in found if allowed_root(r)], [], (along, norm, w)
    boundaries = [rational_root(F(0)), rational_root(F(1))]
    for p, _ in constraints:
        if p:
            boundaries += [r for r in roots(p, F(0), F(1)) if r.compare(0) > 0 and r.compare(1) < 0]
    boundaries = separated(boundaries)
    accepted = [allowed_root(r) for r in boundaries]
    cell = [allowed_value((x.high+y.low)/2) for x, y in zip(boundaries, boundaries[1:])]
    intervals, points = [], []
    i = 0
    while i < len(boundaries):
        if i < len(cell) and cell[i]:
            j = i
            while j < len(cell) and cell[j]:
                j += 1
            assert accepted[i] and accepted[j]
            intervals.append([boundaries[i], boundaries[j]])
            i = j+1
            continue
        if accepted[i] and not (i > 0 and cell[i-1]):
            points.append(boundaries[i])
        i += 1
    return points, intervals, (along, norm, w)


def linear_bounds(root, along, norm, w):
    """Rigorous bracket for s=along/(norm*W) at the root; exact when rational."""
    if not norm:
        return F(0), F(0)
    if root.low == root.high:
        x = root.low
        s = evaluate(coefficients(along), x)/(norm*evaluate(coefficients(w), x))
        return s, s
    width = max(abs(root.low), abs(root.high), F(1))
    for _ in range(256):
        if root.high-root.low < width/F(2)**80:
            break
        root.refine()

    def interval(p):
        lo = hi = F(0)
        for c in reversed(coefficients(p)):
            products = [lo*root.low, lo*root.high, hi*root.low, hi*root.high]
            lo, hi = min(products)+c, max(products)+c
        return lo, hi
    (al, ah), (wl, wh) = interval(along), interval(w)
    assert wl > 0
    candidates = [x/(norm*y) for x in (al, ah) for y in (wl, wh)]
    return min(candidates), max(candidates)


def encode_root(first, last, root, extra):
    # Global parameter = first + (last-first)*x for the unique root x of P in
    # [lo, hi]. Rational values use the unit map for a canonical identity.
    if root.low == root.high:
        u = first+(last-first)*root.low
        return [F(0), F(1), 2, -u, F(1), u, u, *extra], u
    p = coefficients(root.factor.monic())
    return [first, last, len(p), *p, root.low, root.high, *extra], None


def preimage(case):
    a, b = case['endpoints']
    d = [y-x for x, y in zip(a, b)]
    if case['kind'] == 'L':
        assert any(d)
    points, intervals = [], []
    for lo, hi, h in cells(case):
        local_points, local_intervals, (along, norm, w) = local_sets(h, a, d, case['kind'])
        if lo == hi:
            # Singleton: the constant Cartesian cell has no finite roots.
            assert not local_intervals or len(local_intervals) == 1
            if local_intervals or local_points:
                local_points = [rational_root(F(0))]
            local_intervals = []

        def encoded(r):
            return encode_root(lo, hi, r, linear_bounds(r, along, norm, w))
        for r in local_points:
            points.append(encoded(r))
        for x, y in local_intervals:
            intervals.append([encoded(x), encoded(y)])
    # Merge only by rational global identity. Across cells, only shared cuts
    # can coincide; irrational endpoints lie strictly inside one cell.
    merged = []
    for interval in intervals:
        if merged and merged[-1][1][1] is not None and merged[-1][1][1] == interval[0][1]:
            merged[-1][1] = interval[1]
        else:
            merged.append(interval)
    unique = []
    for p in points:
        if unique and p[1] is not None and unique[-1][1] == p[1]:
            continue
        if p[1] is not None and any(p[1] in (x[1], y[1]) for x, y in merged):
            continue
        unique.append(p)
    return [p[0] for p in unique], [[x[0], y[0]] for x, y in merged]


def words(values):
    return ' '.join(str(v) for v in values)


def generate():
    inputs, rows = [], []
    counts = [0, 0]
    for case in make_cases():
        values = [*case['range'], *case['endpoints'][0], *case['endpoints'][1],
                  *(x for c in case['controls'] for x in c),
                  *(x for pair in zip(case['knots'], case['multiplicities']) for x in pair)]
        # Binary64 entry points are also checked when every input is exact.
        binary = all(abs(v) <= F(2)**1023 and F(float(v)) == v for v in values)
        assert binary or not case['native'], case['name']
        inputs.append(words([case['name'], case['kind'], int(binary), case['degree'],
                             int(case['periodic']), len(case['controls']), len(case['knots']), *values]))
        points, intervals = preimage(case)
        row = [case['name'], len(points), len(intervals)]
        for p in points:
            row += p
        for x, y in intervals:
            row += x+y
        rows.append(words(row))
        counts[0] += len(points)
        counts[1] += len(intervals)
        print(f"{case['name']}: {len(points)} points, {len(intervals)} intervals", flush=True)
    print(f'total: {len(inputs)} cases, {counts[0]} points, {counts[1]} intervals')
    return {
        'spline-linear-inputs.txt': '# name kind binary degree periodic poles knots first last A B [x y z w] [knot mult]\n'+'\n'.join(inputs)+'\n',
        'spline-linear.tsv': '# name points intervals; each parameter: first last n P lo hi s_lo s_hi\n'+'\n'.join(rows)+'\n',
    }


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
