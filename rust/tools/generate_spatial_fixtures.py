#!/usr/bin/env python3
"""Independent Fraction oracles for exact spatial predicates and constructions.

Uses rational matrix elimination and barycentric coordinates. Rust uses integer
determinant expansion and edge half-planes. No expectations come from Rust/OCCT.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import random
import struct

ROOT = Path(__file__).resolve().parents[1]


def bits(x):
    return f"{struct.unpack('>Q', struct.pack('>d', x))[0]:016x}"


def value(bits):
    return struct.unpack('>d', struct.pack('>Q', bits))[0]


def determinant(rows):
    a = [list(row) for row in rows]
    result = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        result *= a[i][i]
        for j in range(i+1, len(a)):
            ratio = a[j][i] / a[i][i]
            a[j] = [x-ratio*y for x, y in zip(a[j], a[i])]
    return result


def orientation(points):
    return determinant([[F(1), *p] for p in points])


def sign(x):
    return (x > 0) - (x < 0)


def sphere(points):
    o = orientation(points[:4])
    if not o:
        return 'D'
    d = determinant([[F(1), *p, sum(x*x for x in p)] for p in points])
    return str(-sign(d)*sign(o))  # inside = 1


def interval(x):
    # CPython converts Fraction to nearest binary64 using integer division.
    # Expand on the appropriate side to the smallest closed finite enclosure.
    maximum = F.from_float(value(0x7fefffffffffffff))
    if abs(x) > maximum:
        raise OverflowError
    nearest = float(x)
    represented = F.from_float(nearest)
    lo = math.nextafter(nearest, -math.inf) if represented > x else nearest
    hi = math.nextafter(nearest, math.inf) if represented < x else nearest
    return [bits(lo), bits(hi)]


def point_result(p, q, t):
    return interval(t) + [b for x, y in zip(p, q) for b in interval(x+t*(y-x))]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def sub(a, b):
    return [x-y for x, y in zip(a, b)]


def intersection(kind, points):
    a, b, c, p, q = points
    normal = cross(sub(b, a), sub(c, a))
    if not any(normal):
        return ['D']
    dp, dq = [orientation([a, b, c, x]) for x in [p, q]]
    if kind == 'l':
        if p == q:
            return ['D']
        if dp == dq:
            return ['C' if dp == 0 else 'N']
        return ['P', *point_result(p, q, dp/(dp-dq))]
    if p == q and kind == 'p':
        return ['P', *point_result(p, q, F(0))] if not dp else ['N']
    if dp and dq and sign(dp) == sign(dq):
        return ['N']
    if kind == 'p':
        if not dp and not dq:
            return ['C']
        return ['P', *point_result(p, q, dp/(dp-dq))]
    # Triangle: solve for two barycentric coordinates in a nondegenerate
    # projected triangle. Clip their affine interpolation and their complement.
    drop = next(i for i, x in enumerate(normal) if x)
    i, j = [i for i in range(3) if i != drop]
    u, v = sub(b, a), sub(c, a)
    det = u[i]*v[j]-u[j]*v[i]

    def barycentric(point):
        w = sub(point, a)
        beta = (w[i]*v[j]-w[j]*v[i])/det
        gamma = (u[i]*w[j]-u[j]*w[i])/det
        return [1-beta-gamma, beta, gamma]

    bp, bq = barycentric(p), barycentric(q)
    if dp or dq:
        t = dp/(dp-dq)
        if any(x+t*(y-x) < 0 for x, y in zip(bp, bq)):
            return ['N']
        return ['P', *point_result(p, q, t)]
    low, high = F(0), F(1)
    for x, y in zip(bp, bq):
        if x < 0 and y < 0:
            return ['N']
        if y > x:
            low = max(low, -x/(y-x))
        elif y < x:
            high = min(high, -x/(y-x))
    if low > high:
        return ['N']
    if p == q:
        low = high = F(0)
    if low == high:
        return ['P', *point_result(p, q, low)]
    return ['S', *point_result(p, q, low), *point_result(p, q, high)]


def generate():
    rng = random.Random(0x0CC7_3D)
    tiny = value(1)
    maximum = value(0x7fefffffffffffff)
    unit = [(0.,0.,0.), (1.,0.,0.), (0.,1.,0.), (0.,0.,1.), (0.5,0.5,0.5)]
    cases = [('unit', unit)]
    for scale in [tiny, 2**-1022, 2**-200, 1., 2**200, maximum]:
        for query in [(0.,0.,0.), (1.,1.,1.), (0.5,0.5,0.5), (2.,0.,0.)]:
            p = [[x*scale for x in row] for row in unit[:4]+[query]]
            if all(math.isfinite(x) for row in p for x in row):
                cases.append((f'scaled_{len(cases)}', p))
    for index in range(400):
        def finite_random():
            while True:
                x = value(rng.getrandbits(64))
                if math.isfinite(x):
                    return x
        cases.append((f'bits_{index}', [[finite_random() for _ in range(3)] for _ in range(5)]))
    for index in range(400):
        scale = math.ldexp(1., rng.randint(-1000, 970))
        a = [float(rng.randint(-100,100))*scale for _ in range(3)]
        b = [float(rng.randint(-100,100))*scale for _ in range(3)]
        c = [float(rng.randint(-100,100))*scale for _ in range(3)]
        d = [x+(y-x)*0.5 for x,y in zip(a,b)]
        if index % 3:
            d[index % 3] = math.nextafter(d[index % 3], math.inf)
        q = a if index % 4 == 0 else [x+(y-x)*0.5 for x,y in zip(b,c)]
        cases.append((f'near_plane_{index}', [a,b,c,d,q]))
    predicate_rows = []
    for label, points in cases:
        rational = [[F.from_float(x) for x in p] for p in points]
        words = [bits(x) for p in points for x in p]
        predicate_rows.append(' '.join(['o', label, *words[:12], str(sign(orientation(rational[:4])))]))
        predicate_rows.append(' '.join(['s', label, *words, sphere(rational)]))

    plane = unit[:3]
    intersection_cases = []
    for label, p, q in [
        ('cross',(.25,.25,-1.),(.25,.25,2.)),
        ('round_third',(.1,.2,-1.),(.3,.6,2.)),
        ('edge',(.5,.5,-1.),(.5,.5,1.)),
        ('vertex',(0.,0.,-1.),(0.,0.,1.)),
        ('miss',(2.,2.,-1.),(2.,2.,1.)),
        ('parallel',(0.,0.,1.),(1.,1.,1.)),
        ('contained',(-1.,.25,0.),(2.,.25,0.)),
        ('edge_overlap',(-1.,0.,0.),(2.,0.,0.)),
        ('touch',(-1.,1.,0.),(1.,-1.,0.)),
        ('point_inside',(.25,.25,0.),(.25,.25,0.)),
        ('point_outside',(2.,2.,0.),(2.,2.,0.)),
        ('point_above',(.25,.25,1.),(.25,.25,1.)),
        ('shallow',(0.,0.,-tiny),(1.,1.,tiny)),
        ('overflow',(-maximum,0.,-1.),(maximum,0.,1.)),
        ('unrepresentable_coordinate',(0.,0.,1.),(maximum,0.,.5)),
        ('unrepresentable_parameter',(0.,0.,1.),(1.,0.,1.+2**-52)),
    ]:
        intersection_cases.append((label, plane+[p,q]))
    # A far-away plane and a line whose exact parameter exceeds f64::MAX.
    intersection_cases.append(('parameter_overflow',
        [(0.,0.,maximum),(1.,0.,maximum),(0.,1.,maximum),(0.,0.,0.),(0.,0.,tiny)]))
    for index in range(240):
        scale = math.ldexp(1., rng.randint(-1070, 1000))
        points = [[rng.randint(-8,8)*scale for _ in range(3)] for _ in range(5)]
        if index % 3 == 0:  # coplanar with many edge/vertex overlaps
            points = [[p[0],p[1],0.] for p in points]
        if index % 11 == 0:
            points[2] = points[1]
        intersection_cases.append((f'integer_{index}',points))
    for index in range(64):
        intersection_cases.append((f'bits_{index}',cases[24+index][1]))
    intersection_rows = []
    for label, points in intersection_cases:
        rational = [[F.from_float(x) for x in p] for p in points]
        words = [bits(x) for p in points for x in p]
        for kind in ['l','p','t']:
            try:
                expected = intersection(kind,rational)
            except OverflowError:
                expected = ['U']
            intersection_rows.append(' '.join([kind,label,*words,*expected]))
    return {
        'predicates3d.tsv': '# Exact Fraction matrix oracle: kind label input_bits expected\n'+'\n'.join(predicate_rows)+'\n',
        'intersections.tsv': '# Exact Fraction barycentric oracle: kind label input_bits status [t_lo t_hi x_lo x_hi y_lo y_hi z_lo z_hi]...\n'+'\n'.join(intersection_rows)+'\n',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for name, expected in generate().items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != expected:
                parser.error(f'{name} differs; investigate before changing expectations')
        else:
            path.write_text(expected)
        print(f'{name}: {len(expected.splitlines())-1} exact oracle cases')


if __name__ == '__main__':
    main()
