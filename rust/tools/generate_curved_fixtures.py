#!/usr/bin/env python3
"""Exact independent oracles for degree-two roots and curved intersections.

Compare roots by polynomial evaluation and vertex position using Fraction.
No floating square root or Rust observation determines an expected result.
Cylinder equations use axial projection; production uses cross products.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import random
from generate_spatial_fixtures import bits, value, sign

ROOT = Path(__file__).resolve().parents[1]


class Root:
    def __init__(self, a=0, b=0, c=0, upper=False, rational=None, multiplicity=1):
        self.a,self.b,self.c,self.upper = a,b,c,upper
        self.rational,self.multiplicity = rational,multiplicity

    def compare(self, x):
        if self.rational is not None:
            return sign(self.rational-x)
        vertex = -self.b/(2*self.a)
        evaluated = self.a*x*x+self.b*x+self.c
        if self.upper:
            return 1 if x < vertex else -sign(evaluated)
        return -1 if x > vertex else sign(evaluated)


def roots(a,b,c):
    if not a:
        if b: return [Root(rational=-c/b)]
        return None if not c else []
    if a < 0: a,b,c = -a,-b,-c
    discriminant = b*b-4*a*c
    if discriminant < 0: return []
    if not discriminant: return [Root(rational=-b/(2*a),multiplicity=2)]
    return [Root(a,b,c),Root(a,b,c,upper=True)]


def enclosure(compare):
    negative = compare(F(0)) < 0
    def cmp_magnitude(bits):
        x = F.from_float(value(bits))
        return -compare(-x) if negative else compare(x)
    low,high = 0,0x7fefffffffffffff
    if cmp_magnitude(high) > 0: raise OverflowError
    while low < high:
        middle = (low+high+1)//2
        if cmp_magnitude(middle) < 0: high = middle-1
        else: low = middle
    upper = low if cmp_magnitude(low) == 0 else low+1
    return [bits(-value(upper)),bits(-value(low))] if negative else [bits(value(low)),bits(value(upper))]


def root_record(root):
    try:
        return [str(root.multiplicity),*enclosure(root.compare)]
    except OverflowError:
        return [str(root.multiplicity),'U']


def dot(a,b): return sum(x*y for x,y in zip(a,b))
def sub(a,b): return [x-y for x,y in zip(a,b)]


def hit_record(root,p,d,kind):
    values = enclosure(root.compare)
    for offset,delta in zip(p,d):
        if delta:
            compare = lambda x: sign(delta)*root.compare((x-offset)/delta)
        else:
            compare = lambda x: sign(offset-x)
        values += enclosure(compare)
    return [kind,*values]


def intersect(kind, values):
    center,axis,radius,p,q = values[:3],values[3:6],values[6],values[7:10],values[10:13]
    if radius <= 0 or (kind.lower() != 's' and not any(axis)): return ['D']
    segment = kind.isupper()
    point = p == q
    if point and not segment: return ['D']
    w,d = sub(p,center),sub(q,p)
    if kind.lower() == 'y':
        length2 = dot(axis,axis)
        wd,dd = dot(w,axis),dot(d,axis)
        a = dot(d,d)-dd*dd/length2
        b = 2*(dot(w,d)-wd*dd/length2)
        c = dot(w,w)-wd*wd/length2-radius*radius
    else:
        a,b,c = dot(d,d),2*dot(w,d),dot(w,w)-radius*radius
    if kind.lower() == 'c':
        initial,direction = dot(w,axis),dot(d,axis)
        if not direction:
            if initial: return ['N']
            candidates = roots(a,b,c)
        else:
            t = -initial/direction
            candidates = [Root(rational=t)] if a*t*t+b*t+c == 0 else []
    else:
        candidates = roots(a,b,c)
    if candidates is None:
        return ['1',*hit_record(Root(rational=F(0)),p,d,'P')] if point else ['A']
    if segment:
        candidates = [root for root in candidates if root.compare(F(0)) >= 0 and root.compare(F(1)) <= 0]
    if not candidates: return ['N']
    try:
        records = [hit_record(root,p,d,'T' if root.multiplicity == 2 else 'K') for root in candidates]
        return [str(len(records)),*[x for record in records for x in record]]
    except OverflowError:
        return ['U']


def generate():
    tiny,maximum = value(1),value(0x7fefffffffffffff)
    rng = random.Random(0x0CC7_C0DE)
    base = [(0.,0.,0.),(0.,0.,1.),(0.,2.,-1.),(1.,-2.,1.),(1.,0.,-2.),
            (1.,0.,1.),(1.,-1e16,1.),(1.,2.,-3.),(1.,0.,0.),(tiny,1.,1.),
            (tiny,maximum,-tiny),(maximum,tiny,-maximum),
            (1.,-2.,math.nextafter(1.,0.)),(1.,-2.,math.nextafter(1.,math.inf)),
            (0.,tiny,-maximum),(1.,-1.,tiny)]
    polynomials = [(f'special_{i}',p) for i,p in enumerate(base)]
    for power in [-1070,-1000,-500,-200,0,200,500,970]:
        for i,p in enumerate(base[:9]):
            for factor in [-1.,1.]:
                polynomials.append((f'scaled_{power}_{i}_{factor}',[math.ldexp(x*factor,power) for x in p]))
    def random_finite():
        while True:
            x = value(rng.getrandbits(64))
            if math.isfinite(x): return x
    for i in range(192):
        polynomials.append((f'bits_{i}',[random_finite() for _ in range(3)]))
    for i in range(96):
        a = rng.randint(1,1000)
        t = rng.randint(-1000,1000)
        c = float(a*t*t)
        if i % 3: c = math.nextafter(c,math.inf if i % 3 == 1 else -math.inf)
        polynomials.append((f'near_double_{i}',[float(a),float(-2*a*t),c]))
    rows = []
    for label,p in polynomials:
        expected = roots(*map(F.from_float,p))
        fields = ['A'] if expected is None else [str(len(expected)),*[x for root in expected for x in root_record(root)]]
        rows.append(' '.join([label,*map(bits,p),*fields]))

    curves = []
    center,axis = [0.,0.,0.],[0.,0.,1.]
    for label,p,q in [
        ('diameter',[-2.,0.,0.],[2.,0.,0.]),('irrational',[-2.,.5,0.],[2.,.5,0.]),
        ('tangent',[-2.,1.,0.],[2.,1.,0.]),('outside',[-2.,2.,0.],[2.,2.,0.]),
        ('near_inside',[-2.,math.nextafter(1.,0.),0.],[2.,math.nextafter(1.,0.),0.]),
        ('near_outside',[-2.,math.nextafter(1.,math.inf),0.],[2.,math.nextafter(1.,math.inf),0.]),
        ('interior',[0.,0.,0.],[.5,0.,0.]),('endpoint',[-1.,0.,0.],[0.,0.,0.]),
        ('both_endpoints',[-1.,0.,0.],[1.,0.,0.]),('clipped',[0.,0.,0.],[2.,0.,0.]),
        ('generator',[1.,0.,-2.],[1.,0.,2.]),('axis',[0.,0.,-2.],[0.,0.,2.]),
        ('skew',[1.,0.,-1.],[1.,0.,1.]),('skew_miss',[2.,0.,-1.],[2.,0.,1.]),
        ('point_on',[1.,0.,0.],[1.,0.,0.]),('point_in',[0.,0.,0.],[0.,0.,0.]),
        ('point_out',[2.,0.,0.],[2.,0.,0.]),
        ('collapsed_parameters',[1e100,0.,0.],[value(int(bits(1e100),16)+3),0.,0.]),
    ]:
        curves.append((label,center+axis+[1.]+p+q))
    for power in [-1070,-1000,-500,-200,0,200,500,970]:
        scale = math.ldexp(1.,power)
        curves += [(f'scale_{power}',[0.,0.,0.]+axis+[scale]+[-2*scale,.5*scale,0.,2*scale,.5*scale,0.]),
                   (f'axis_scale_{power}',center+[0.,0.,scale]+[1.]+[-2.,.5,0.,2.,.5,0.])]
    curves += [('overflow_points',[maximum,0.,0.]+axis+[maximum]+[0.,0.,0.,maximum,0.,0.]),
               ('tiny_direction',[maximum,0.,0.]+axis+[1.]+[0.,0.,0.,tiny,0.,0.]),
               ('zero_radius',center+axis+[0.]+[-2.,0.,0.,2.,0.,0.]),
               ('zero_axis',center+center+[1.]+[-2.,0.,0.,2.,0.,0.])]
    for i in range(96):
        scale = math.ldexp(1.,rng.randint(-1000,990))
        values = [rng.randint(-8,8)*scale for _ in range(13)]
        values[6] = rng.randint(1,8)*scale
        if i % 3 == 0:
            values[2]=values[9]=values[12]=0.
            values[3:6] = [0.,0.,scale]
        curves.append((f'integer_{i}',values))
    for i in range(48):
        values = [random_finite() for _ in range(13)]
        values[6] = abs(values[6]) or tiny
        curves.append((f'bits_{i}',values))
    curved_rows = []
    for label,values in curves:
        rational = list(map(F.from_float,values))
        for kind in ['s','S','y','Y','c','C']:
            curved_rows.append(' '.join([kind,label,*map(bits,values),*intersect(kind,rational)]))
    return {'quadratic.tsv':'# label a_bits b_bits c_bits count/A [multiplicity lo hi / U]...\n'+'\n'.join(rows)+'\n',
            'curved.tsv':'# kind label center axis radius p q (binary64 bits) status [kind t_lo t_hi x_lo x_hi y_lo y_hi z_lo z_hi]...\n'+'\n'.join(curved_rows)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    for name,expected in generate().items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != expected: parser.error(f'{name} changed; investigate before updating the baseline')
        else: path.write_text(expected)
        print(f'{name}: {len(expected.splitlines())-1} independent exact cases')


if __name__ == '__main__': main()
