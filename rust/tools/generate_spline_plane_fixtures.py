#!/usr/bin/env python3
"""Exact spline/plane oracle using basis polynomials and independent real roots."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import random
from compare_spline_plane import cases, encode
from exact_polynomial_oracle import add, multiply, roots
from generate_spline_fixtures import axis
from generate_curved_fixtures import enclosure
from generate_spatial_fixtures import bits, value

ROOT = Path(__file__).resolve().parents[1]


def parse(row):
    w = row.split()
    degree, np, nk = map(int, w[2:5])
    plane = [list(map(lambda x: F(float(x)), w[5+3*i:8+3*i])) for i in range(3)]
    data = list(map(lambda x: F(float(x)), w[14:14+4*np]))
    tail = w[14+4*np:]
    if len(tail) != 2*nk: raise ValueError('invalid knot count')
    return (degree, plane, [data[i:i+3] for i in range(0, len(data), 4)], data[3::4],
            list(map(float, tail[::2])), list(map(int, tail[1::2])), w[1] == 'P')


def homogeneous(degree, poles, weights, flat, span):
    start, length = flat[span], flat[span+1]-flat[span]
    @lru_cache(None)
    def basis(i, p):
        if not p: return [F(1)] if i == span else []
        result = []
        left, right = flat[i+p]-flat[i], flat[i+p+1]-flat[i+1]
        if left:
            result = multiply([(start-flat[i])/left, length/left], basis(i, p-1))
        if right:
            result = add(result, multiply([(flat[i+p+1]-start)/right, -length/right], basis(i+1, p-1)))
        return result
    h = [[] for _ in range(4)]
    for i in range(len(flat)-degree-1):
        j = i % len(poles)
        for c in range(4):
            scale = weights[j]*(poles[j][c] if c < 3 else 1)
            h[c] = add(h[c], [scale*x for x in basis(i, degree)])
    return h


def expected(row):
    degree, plane, poles, weights, knots, mults, periodic = parse(row)
    flat, start, end = axis(degree, knots, mults, periodic)
    u, v = [[p[c]-plane[0][c] for c in range(3)] for p in plane[1:]]
    normal = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
    candidates, overlaps = [], []
    for span in range(degree, len(flat)-degree-1):
        a, b = flat[span:span+2]
        if a == b or a < start or b > end: continue
        h = homogeneous(degree, poles, weights, flat, span)
        f = []
        for c in range(3):
            f = add(f, [normal[c]*x for x in add(h[c], [-plane[0][c]*w for w in h[3]])])
        local = roots(f, F(0), F(1))
        if local is None:
            if overlaps and overlaps[-1][1] == a: overlaps[-1][1] = b
            else: overlaps.append([a,b])
            continue
        # Determine side signs by exact rational probes BETWEEN adjacent roots,
        # independently of production's derivative/multiplicity parity rule.
        for left, right in zip(local, local[1:]):
            while left.high >= right.low: left.refine(); right.refine()
        for i, root in enumerate(local):
            from exact_polynomial_oracle import evaluate, sign
            at_start, at_end = root.compare(0) == 0, root.compare(1) == 0
            knot = a if at_start else b if at_end else None
            while not at_start and root.low <= 0: root.refine()
            while not at_end and root.high >= 1: root.refine()
            left_probe = ((local[i-1].high if i else F(0)) + root.low)/2
            right_probe = (root.high + (local[i+1].low if i+1 < len(local) else F(1)))/2
            orders = [0 if at_start else root.multiplicity, 0 if at_end else root.multiplicity]
            signs = [None if at_start else sign(evaluate(f,left_probe)), None if at_end else sign(evaluate(f,right_probe))]
            assert all(x in [-1,1] for x in signs if x is not None)
            if knot is not None and candidates and candidates[-1]['knot'] == knot:
                candidates[-1]['orders'][1] = orders[1]; candidates[-1]['signs'][1] = signs[1]
                continue
            bounds = enclosure(lambda x: root.compare((x-a)/(b-a)))
            for c in range(3):
                bounds += enclosure(lambda x: root.sign_at(add(h[c], [-x*w for w in h[3]])))
            candidates.append({'knot':knot, 'bounds':bounds, 'orders':orders, 'signs':signs})
    points = []
    for item in candidates:
        if item['knot'] is not None and any(a <= item['knot'] <= b for a,b in overlaps): continue
        left,right = item['signs']
        contact = 'B' if None in [left,right] else 'T' if left == right else 'C'
        points.append([*item['bounds'], contact, *map(str,item['orders'])])
    return [str(len(points)), str(len(overlaps)), *[x for p in points for x in p], *[bits(float(x)) for pair in overlaps for x in pair]]


def inputs():
    result = cases().splitlines()
    plane = [(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)]
    tiny, maximum = value(1), value(0x7fefffffffffffff)
    for name, d, poles, weights, knots, mults in [
        ('huge_linear',1,[(-maximum,0.,-maximum),(maximum,0.,maximum)],[1.,1.],[-maximum,maximum],[2,2]),
        ('tiny_span',2,[(0.,0.,1.),(1.,2.,-1.),(2.,0.,1.)],[1.]*3,[0.,tiny],[3,3]),
        ('extreme_weights',2,[(0.,0.,-1.),(1.,2.,-1.),(2.,0.,1.)],[tiny,1.,maximum],[0.,1.],[3,3]),
        ('unequal_orders',2,[(0.,0.,-1.),(1.,0.,0.),(2.,0.,0.),(3.,0.,1.),(4.,0.,2.)],[1.]*5,[0.,1.,2.],[3,2,3]),
        ('overlap_chain',1,[(float(i),0.,z) for i,z in enumerate([-1.,0.,0.,0.,1.,0.,0.])],[1.]*7,list(map(float,range(7))),[2,1,1,1,1,1,2]),
    ]: result.append(encode(name,'S',d,plane,poles,weights,knots,mults))
    rng = random.Random(0x1A73_235)
    for i in range(48):
        d = 1+i%6
        periodic = i%3 == 0
        mults = [1,d,1,1] if periodic else [d+1,d,d+1]
        knots = [0.,.25,1.5,3.] if periodic else [-1.,.5,2.]
        np = sum(mults[:-1]) if periodic else sum(mults)-d-1
        poles = [[float(rng.randint(-8,8)) for _ in range(3)] for _ in range(np)]
        weights = [float(rng.randint(1,5)) for _ in range(np)]
        result.append(encode(f'random_{i}','P' if periodic else 'S',d,plane,poles,weights,knots,mults))
    return result


def generate():
    return '# input text | point_count overlap_count [t/x/y/z lo hi bits contact left_order right_order] [overlap endpoints bits]\n'+'\n'.join(row+' | '+' '.join(expected(row)) for row in inputs())+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    text=generate(); path=ROOT/'fixtures/spline-plane.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('spline-plane fixtures changed; investigate before updating')
    else: path.write_text(text)
    print(f'spline-plane.tsv: {len(text.splitlines())-1} independent exact cases')


if __name__=='__main__': main()
