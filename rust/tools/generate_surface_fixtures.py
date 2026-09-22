#!/usr/bin/env python3
"""Independent exact tensor basis + closed quotient surface oracle."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import random
from compare_surfaces import cases, encode
from generate_spline_fixtures import axis, parameters, basis_jet
from generate_spatial_fixtures import bits, value, sign
from generate_curved_fixtures import enclosure

PARTIALS=[(0,0),(1,0),(0,1),(2,0),(0,2),(1,1)]


def parse(row):
    words=row.split()
    name,kind=words[:2]
    du,dv,nu,nv,ku,kv,pu,pv=map(int,words[2:10])
    u,v=map(float,words[10:12])
    su,sv=words[12:14]
    order=int(words[14])
    data=list(map(float,words[15:15+4*nu*nv]))
    tail=words[15+4*nu*nv:]
    assert len(tail)==2*(ku+kv)
    uk=list(map(float,tail[:2*ku:2])); um=list(map(int,tail[1:2*ku:2]))
    vk=list(map(float,tail[2*ku::2])); vm=list(map(int,tail[2*ku+1::2]))
    return name,kind,du,dv,[data[i:i+3] for i in range(0,len(data),4)],data[3::4],uk,um,vk,vm,pu,pv,u,v,su,sv,order


def jet(du,dv,poles,weights,nu,nv,uf,vf,uat,vat,order):
    ub=basis_jet(du,uf,*uat,order)
    vb=basis_jet(dv,vf,*vat,order)
    count=[1,3,6][order]
    h=[]
    for a,b in PARTIALS[:count]:
        h.append([sum(ub[i][a]*vb[j][b]*weights[(i%nu)*nv+j%nv]*(poles[(i%nu)*nv+j%nv][c] if c<3 else 1)
                      for i in range(len(ub)) for j in range(len(vb))) for c in range(4)])
    w=h[0][3]
    result=[[x/w for x in h[0][:3]]]
    if order>=1:
        result += [[(h[r][c]*w-h[0][c]*h[r][3])/(w*w) for c in range(3)] for r in [1,2]]
    if order>=2:
        result += [[(h[r][c]*w*w-h[0][c]*h[r][3]*w-2*h[first][c]*w*h[first][3]+2*h[0][c]*h[first][3]**2)/(w**3) for c in range(3)]
                   for r,first in [(3,1),(4,2)]]
        result.append([(h[5][c]*w*w-h[1][c]*w*h[2][3]-h[2][c]*w*h[1][3]-h[0][c]*w*h[5][3]+2*h[0][c]*h[1][3]*h[2][3])/(w**3) for c in range(3)])
    return result


def expected(shape):
    _,_,du,dv,poles,weights,uk,um,vk,vm,pu,pv,u,v,su,sv,order=shape
    if not all(map(math.isfinite,[u,v])): return ['F']
    uf,us,ue=axis(du,uk,um,pu)
    vf,vs,ve=axis(dv,vk,vm,pv)
    try:
        uq=parameters(uf,us,ue,u,su,pu)
        vq=parameters(vf,vs,ve,v,sv,pv)
    except IndexError: return ['O']
    nu=sum(um[:-1]) if pu else sum(um)-du-1
    nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
    poles=[list(map(F,p)) for p in poles]
    weights=list(map(F,weights))
    values=jet(du,dv,poles,weights,nu,nv,uf,vf,uq[0],vq[0],order)
    for i,a in enumerate(uq):
        for j,b in enumerate(vq):
            if (i or j) and values!=jet(du,dv,poles,weights,nu,nv,uf,vf,a,b,order): return ['D']
    try: return ['P',*[x for row in values for val in row for x in enclosure(lambda x: sign(val-x))]]
    except OverflowError: return ['U']


def generate():
    inputs=cases().splitlines()
    # Intersecting C0 knots/seams: all four side combinations and automatic.
    for pu,pv in [(False,False),(True,False),(False,True),(True,True)]:
        um=[1,1,1] if pu else [2,1,2]
        vm=[1,1,1] if pv else [2,1,2]
        nu=sum(um[:-1]) if pu else 3
        nv=sum(vm[:-1]) if pv else 3
        for constant in [False,True]:
            poles=[(1.,2.,3.) if constant else (float(i%2),float(j%2),float(i*j)) for i in range(nu) for j in range(nv)]
            weights=[1.+i%2 for i in range(nu*nv)]
            for u,v in [(0.,0.),(1.,1.),(3.,3.)]:
                for su,sv in [('A','A'),('L','A'),('A','R'),('L','L'),('L','R'),('R','L'),('R','R')]:
                    for order in [0,1,2]:
                        inputs.append(encode(f'join_{pu}_{pv}_{constant}_{u}_{su}{sv}_{order}','S',1,1,poles,weights,[0.,1.,3.],um,[0.,1.,3.],vm,pu,pv,u,v,su,sv,order))
    tiny,maximum=value(1),value(0x7fefffffffffffff)
    for name,knots in [('tiny',[0.,tiny]),('huge',[-maximum,maximum])]:
        for order in [0,1,2]:
            inputs.append(encode(f'{name}_{order}','S',1,1,[(0.,0.,0.),(1.,2.,3.),(1.,3.,4.),(2.,5.,6.)],[1.,2.,3.,4.],knots,[2,2],knots,[2,2],False,False,knots[0],knots[0],'A','A',order))
    for order in [0,1,2]:
        inputs.append(encode(f'homogeneous_overflow_{order}','B',1,1,[(maximum,2.,-maximum)]*4,[tiny,maximum,tiny,maximum],[0.,1.],[2,2],[0.,1.],[2,2],False,False,.5,.5,'A','A',order))
    for parameter in [maximum,-maximum,1.,-1.]:
        for order in [0,2]:
            inputs.append(encode(f'periodic_tiny_{parameter}_{order}','S',1,1,[(0.,0.,0.),(1.,2.,3.),(2.,3.,4.),(3.,4.,5.)],[1.,2.,3.,4.],[0.,tiny,3*tiny],[1]*3,[0.,tiny,3*tiny],[1]*3,True,True,parameter,-parameter,'R','L',order))
    rng=random.Random(0x5FACE)
    for i in range(60):
        du,dv=1+i%3,1+(i//3)%3
        pu,pv=bool(i%2),bool(i%4>=2)
        scale=math.ldexp(1.,rng.randint(-950,950))
        uk,vk=[0.,scale,3*scale],[0.,scale,3*scale]
        um=[1,du,1] if pu else [du+1,1,du+1]
        vm=[1,dv,1] if pv else [dv+1,1,dv+1]
        nu=sum(um[:-1]) if pu else sum(um)-du-1
        nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
        poles=[[rng.randint(-8,8)*scale for _ in range(3)] for _ in range(nu*nv)]
        weights=[math.ldexp(float(rng.randint(1,8)),rng.randint(-100,100)) for _ in poles]
        inputs.append(encode(f'wide_{i}','S',du,dv,poles,weights,uk,um,vk,vm,pu,pv,rng.choice([0.,.5,1.,3.])*scale,rng.choice([0.,.5,1.,3.])*scale,rng.choice(['L','R','A']),rng.choice(['L','R','A']),2))
    rows=[]
    for row in inputs:
        shape=parse(row)
        words=row.split()
        nu,nv,ku,kv=map(int,words[4:8])
        words[10:12]=[bits(float(x)) for x in words[10:12]]
        words[15:15+4*nu*nv]=[bits(float(x)) for x in words[15:15+4*nu*nv]]
        for i in range(15+4*nu*nv,len(words),2): words[i]=bits(float(words[i]))
        rows.append(' '.join(words+expected(shape)))
    return '# label kind du dv nu nv ku kv pu pv u_bits v_bits su sv order [xyzw bits] [uknot mult] [vknot mult] status [lo hi bits]...\n'+'\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    text=generate()
    path=Path(__file__).resolve().parents[1]/'fixtures/surfaces.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('surface fixture changed; investigate before updating')
    else: path.write_text(text)
    print(f'surfaces.tsv: {len(text.splitlines())-1} independent exact cases')


if __name__=='__main__': main()
