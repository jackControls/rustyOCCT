#!/usr/bin/env python3
"""Independent Fraction Cox-de Boor basis/derivative oracle, not pole de Boor.

The basis derivative identity and closed quotient formulas provide independent
checks on the differentiated homogeneous interpolation used by production.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path
import random
from compare_splines import cases
from generate_spatial_fixtures import bits, value, sign
from generate_curved_fixtures import enclosure

ROOT=Path(__file__).resolve().parents[1]


def basis_jet(degree, flat, u, span, order):
    @lru_cache(None)
    def basis(i,p,r):
        if r>p: return F(0)
        if p==0: return F(i==span)
        left,right=flat[i+p]-flat[i],flat[i+p+1]-flat[i+1]
        if r:
            return ((p*basis(i,p-1,r-1)/left) if left else 0)-((p*basis(i+1,p-1,r-1)/right) if right else 0)
        return (((u-flat[i])*basis(i,p-1,0)/left) if left else 0)+(((flat[i+p+1]-u)*basis(i+1,p-1,0)/right) if right else 0)
    return [[basis(i,degree,r) for r in range(order+1)] for i in range(len(flat)-degree-1)]


def jet(degree, poles, weights, flat, u, span, order):
    basis=basis_jet(degree,flat,u,span,order)
    h=[[sum(basis[i][r]*w*(p[c] if c<3 else 1) for i,(p,w) in enumerate(zip(poles,weights))) for c in range(4)] for r in range(order+1)]
    w=h[0][3]
    result=[[x/w for x in h[0][:3]]]
    if order>=1:
        result.append([(h[1][c]*w-h[0][c]*h[1][3])/(w*w) for c in range(3)])
    if order>=2:
        result.append([(h[2][c]*w*w-h[0][c]*h[2][3]*w-2*h[1][c]*w*h[1][3]+2*h[0][c]*h[1][3]*h[1][3])/(w*w*w) for c in range(3)])
    return result


def axis(degree,knots,mults,periodic):
    flat=[F(k) for k,m in zip(knots,mults) for _ in range(m)]
    if not periodic:
        return flat,flat[degree],flat[len(flat)-degree-1]
    start,end=F(knots[0]),F(knots[-1])
    cycle=flat[:-mults[-1]]
    # Infinite periodic knot function indexed by Euclidean quotient/remainder.
    # This differs from production's two slices of extended knot data.
    offset=degree+1-mults[0]
    flat=[cycle[i%len(cycle)]+(i//len(cycle))*(end-start) for i in range(-offset,len(cycle)+mults[0]+offset)]
    return flat,start,end


def parameters(flat,start,end,u,side,periodic):
    if not math.isfinite(u): raise ValueError('nonfinite')
    u=F(u)
    if periodic: u=start+(u-start)%(end-start)
    elif u<start or u>end or (u==start and side=='L') or (u==end and side=='R'): raise IndexError('outside domain')
    seam=periodic and u==start
    def at(x,left): return (x,max(i for i,k in enumerate(flat) if (k<x if left else k<=x)))
    selected=at(end,True) if seam and side=='L' else at(u,side=='L' or u==end)
    choices=[selected]
    # Always compare all actual jets at automatic interior knots/seams. The
    # production multiplicity shortcut is intentionally absent from this oracle.
    if side=='A' and (seam or (start<u<end and u in flat)):
        choices.append(at(end if seam else u,True))
    return choices


def expected(degree,poles,weights,knots,mults,u,side,order,periodic=False):
    if not math.isfinite(u): return ['F']
    flat,start,end=axis(degree,knots,mults,periodic)
    try: choices=parameters(flat,start,end,u,side,periodic)
    except IndexError: return ['O']
    rational_poles=[[F(x) for x in poles[i%len(poles)]] for i in range(len(flat)-degree-1)]
    rational_weights=[F(weights[i%len(weights)]) for i in range(len(rational_poles))]
    u,span=choices[0]
    v=jet(degree,rational_poles,rational_weights,flat,u,span,order)
    for u,span in choices[1:]:
        if v!=jet(degree,rational_poles,rational_weights,flat,u,span,order): return ['D']
    try:
        return ['P',*[x for row in v for val in row for x in enclosure(lambda x: sign(val-x))]]
    except OverflowError: return ['U']


def generate():
    inputs=[]
    for row in cases().splitlines():
        w=row.split()
        name,kind,degree,np,nk,u,side,order=w[:8]
        degree,np,nk,order=map(int,[degree,np,nk,order])
        data=list(map(float,w[8:8+4*np]))
        tail=w[8+4*np:]
        inputs.append((name,kind,degree,[data[i:i+3] for i in range(0,len(data),4)],data[3::4],list(map(float,tail[::2])),list(map(int,tail[1::2])),float(u),side,order))
    tiny,maximum=value(1),value(0x7fefffffffffffff)
    # C0 and C1 joins, including an actually smooth curve with high multiplicity.
    for name,degree,poles,mults in [
        ('corner',1,[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.)],[2,1,2]),
        ('c1',2,[(0.,0.,0.),(1.,0.,0.),(2.,1.,0.),(4.,3.,1.)],[3,1,3]),
        ('smooth_repeated',2,[(1.,2.,3.)]*5,[3,2,3]),
    ]:
        for side in ['A','L','R']:
            for order in range(3):
                inputs.append((f'{name}_{side}_{order}','S',degree,poles,[1.]*len(poles),[0.,1.,2.],mults,1.,side,order))
    for name,degree,poles,weights,knots,mults,parameters in [
        ('tiny_span',1,[(0.,0.,0.),(1.,0.,0.)],[1.,1.],[0.,tiny],[2,2],[0.,tiny]),
        ('huge_span',1,[(-maximum,0.,0.),(maximum,0.,0.)],[1.,1.],[-maximum,maximum],[2,2],[-maximum,0.,maximum]),
        ('cancellation',1,[(-maximum,0.,0.),(maximum,0.,0.)],[1.,1.],[0.,1.],[2,2],[0.,.5,1.]),
        ('homogeneous_overflow',2,[(maximum,1.,0.)]*3,[maximum,tiny,maximum],[0.,1.],[3,3],[0.,.5,1.]),
        ('unequal_tiny_weights',2,[(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)],[tiny,2*tiny,3*tiny],[0.,1.],[3,3],[0.,.5,1.]),
        ('neighbor_knots',2,[(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(3.,1.,0.)],[1.]*4,[1.,math.nextafter(1.,math.inf),2.],[3,1,3],[1.,math.nextafter(1.,math.inf),1.5,2.]),
    ]:
        for i,u in enumerate(parameters):
            for order in range(3):
                inputs.append((f'{name}_{i}_{order}','S',degree,poles,weights,knots,mults,u,'A',order))
    rng=random.Random(0xB5_911E)
    for i in range(90):
        degree=1+i%5
        scale=math.ldexp(1.,rng.randint(-1000,970))
        mults=[degree+1,1+i%degree,degree+1]
        np=sum(mults)-degree-1
        poles=[[rng.randint(-8,8)*scale for _ in range(3)] for _ in range(np)]
        weights=[math.ldexp(float(rng.randint(1,8)),rng.randint(-500,500)) for _ in range(np)]
        knots=[0.,scale,3*scale]
        u=rng.choice([0.,.25,.5,1.,2.,3.])*scale
        inputs.append((f'wide_{i}','S',degree,poles,weights,knots,mults,u,rng.choice(['A','L','R']),2))
    for degree,seam in [(1,1),(2,1),(2,2),(3,1),(3,3)]:
        mults=[seam,degree,1,seam]
        np=sum(mults[:-1])
        for constant in [False,True]:
            poles=[(2.,3.,4.) if constant else (float(i%3),float(i*i),float(2-i)) for i in range(np)]
            for u in [0.,3.,-6.,.5]:
                for side in ['A','L','R']:
                    for order in range(3):
                        inputs.append((f'periodic_continuity_{degree}_{seam}_{constant}_{u}_{side}_{order}','P',degree,poles,[1.+i%2 for i in range(np)],[0.,.5,2.,3.],mults,u,side,order))
    for name,knots,parameters_ in [
        ('subnormal_period',[0.,tiny,3*tiny],[maximum,-maximum,1.,-1.,tiny]),
        ('unrepresentable_period',[-maximum,0.,maximum],[-maximum,0.,maximum]),
        ('nonrepresentable_remainder',[.1,.4,1.3],[maximum,-maximum,3.,-3.]),
    ]:
        for u in parameters_:
            for side in ['A','L','R']:
                for order in [0,2]:
                    inputs.append((f'{name}_{u}_{side}_{order}','P',1,[(0.,0.,0.),(1.,2.,3.)],[1.,2.],knots,[1,1,1],u,side,order))
    rows=[]
    for name,kind,degree,poles,weights,knots,mults,u,side,order in inputs:
        fields=[name,kind,str(degree),str(len(poles)),str(len(knots)),bits(u),side,str(order)]
        fields += [bits(x) for p,w in zip(poles,weights) for x in [*p,w]]
        fields += [x for k,m in zip(knots,mults) for x in [bits(k),str(m)]]
        rows.append(' '.join(fields+expected(degree,poles,weights,knots,mults,u,side,order,kind=='P')))
    return '# label kind degree np nk u_bits side order [x y z w bits] [knot_bits mult] status [lo hi bits]...\n'+'\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    text=generate()
    path=ROOT/'fixtures/splines.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('spline fixture changed; investigate before updating')
    else: path.write_text(text)
    print(f'splines.tsv: {len(text.splitlines())-1} independent exact cases')


if __name__=='__main__': main()
