#!/usr/bin/env python3
"""Independent exact homogeneous controls, including unsupported native extremes."""
import argparse
import math
from pathlib import Path
import random
from bezier_editing_reference import expected, encode as encode_result
from compare_bezier_editing import cases, encode

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    rows=cases().splitlines()
    tiny=float.fromhex('0x0.0000000000001p-1022'); maximum=float.fromhex('0x1.fffffffffffffp+1023')
    for name,kind,d,poles,weights,knots,mults,first,last in [
        ('huge_span','S',1,[(-maximum,0.,0.),(maximum,0.,0.)],[1.,1.],[-maximum,maximum],[2,2],-maximum,maximum),
        ('tiny_span','S',2,[(0.,0.,0.),(1.,1.,0.),(2.,0.,1.)],[1.,2.,1.],[0.,tiny],[3,3],0.,tiny),
        ('extreme_homogeneous','S',2,[(maximum,0.,0.),(0.,-maximum,1.),(-maximum,0.,maximum)],[tiny,maximum,1.],[0.,1.],[3,3],0.,1.),
        ('adjacent_knots','S',2,[(0.,0.,0.),(1.,1.,0.),(2.,0.,1.),(3.,-1.,0.)],[1.,2.,1.,3.],[1.,math.nextafter(1.,math.inf),2.],[3,1,3],1.,2.),
        ('subulp_arcs','P',1,[(0.,0.,0.),(1.,2.,3.),(-1.,0.,1.)],[1.,2.,3.],[0.,.125,.375,.5],[1]*4,2.**53,2.**53+2.),
        ('subnormal_turns','P',1,[(0.,0.,0.),(1.,2.,3.)],[1.,2.],[0.,tiny,4*tiny],[1]*3,tiny,9*tiny),
        ('negative_turns','P',2,[(0.,0.,0.),(1.,2.,3.),(2.,0.,1.),(3.,1.,0.)],[1.,2.,3.,2.],[.1,.4,1.3],[2,2,2],-4.,3.),
    ]:
        for op in range(6): rows.append(encode(f'{name}_op{op}',kind,d,poles,weights,knots,mults,first,last,op,d+3))
    rng=random.Random(0xBE21_E2)
    for i in range(48):
        d=1+i%8; periodic=i%3==0; m=1+rng.randrange(d)
        mults=[m,d,1,m] if periodic else [d+1,m,d+1]
        knots=[-1.,.25,2.,3.] if periodic else [-1.,.25,3.]
        n=sum(mults[:-1]) if periodic else sum(mults)-d-1
        scale=math.ldexp(1.,rng.randint(-900,900))
        poles=[[scale*rng.randint(-8,8) for _ in range(3)] for _ in range(n)]
        weights=[math.ldexp(float(rng.randint(1,8)),rng.randint(-20,20)) for _ in range(n)]
        first,last=(-2.,5.) if periodic else (-.5,2.75)
        rows.append(encode(f'random_{i}','P' if periodic else 'S',d,poles,weights,knots,mults,first,last,i%6,d+2))
    return rows


def generate():
    return '# input | exact output: label R count [degree first last (wx wy wz w)*]\n'+'\n'.join(row+' | '+encode_result(row.split()[0],expected(row)) for row in inputs())+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    text=generate(); path=ROOT/'fixtures/bezier-editing.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('Bezier editing fixtures changed; investigate before updating')
    else: path.write_text(text)
    print(f'bezier-editing.tsv: {len(text.splitlines())-1} independent exact cases, {len(text)} bytes')


if __name__=='__main__': main()
