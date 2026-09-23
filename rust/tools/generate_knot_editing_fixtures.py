#!/usr/bin/env python3
"""Complete knot editing expectations from independent coefficient equations."""
import argparse
import math
from pathlib import Path
import random
from compare_knot_editing import cases, encode
from knot_editing_reference import expected, encode as encode_result

ROOT = Path(__file__).resolve().parents[1]


def inputs():
    rows = cases().splitlines()
    tiny = float.fromhex('0x0.0000000000001p-1022')
    maximum = float.fromhex('0x1.fffffffffffffp+1023')
    for name,d,poles,weights,knots,cut in [
        ('huge_domain',1,[(-maximum,0.,0.),(maximum,0.,0.)],[1.,1.],[-maximum,maximum],0.),
        ('tiny_domain',2,[(0.,0.,0.),(1.,2.,3.),(4.,0.,1.)],[1.,3.,2.],[0.,4*tiny],tiny),
        ('huge_homogeneous',2,[(maximum,0.,0.),(0.,-maximum,1.),(-maximum,0.,maximum)],[tiny,maximum,1.],[0.,1.],.25),
        ('adjacent_cut',2,[(0.,0.,0.),(1.,2.,3.),(4.,0.,1.)],[1.,3.,2.],[1.,2.],math.nextafter(1.,2.)),
    ]:
        shape = name,'S',d,poles,weights,knots,[d+1,d+1]
        rows.append(encode(shape,'insert',[('I',cut,d)]))
        rows.append(encode(shape,'roundtrip',[('I',cut,d),('R',cut,0)]))
    rng = random.Random(0xB5_91_1E)
    for i in range(48):
        d = 1+i%8; periodic = i%3==0; m = 1+rng.randrange(d)
        mults = [m,d,1,m] if periodic else [d+1,m,d+1]
        knots = [-1.,.25,2.,3.] if periodic else [-1.,.25,3.]
        n = sum(mults)-(mults[0] if periodic else d+1)
        scale = math.ldexp(1.,rng.randint(-900,900))
        poles = [[scale*rng.randint(-8,8) for _ in range(3)] for _ in range(n)]
        weights = [math.ldexp(float(rng.randint(1,8)),rng.randint(-20,20)) for _ in range(n)]
        shape = f'random_{i}','P' if periodic else 'S',d,poles,weights,knots,mults
        rows.append(encode(shape,'sequence',[('I',-.5,d),('I',1.,1),('R',-.5,0),('R',1.,0)]))
    return rows


def generate():
    rows = ['# input | exact output: flags, degree, periodic, complete homogeneous controls and knots']
    for i,row in enumerate(inputs()):
        flags,curve = expected(row)
        rows.append(row+' | '+encode_result(row.split()[0],flags,curve))
        if (i+1)%100==0: print(f'{i+1} independent knot editing cases',flush=True)
    return '\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    text=generate(); path=ROOT/'fixtures/knot-editing.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('knot editing fixtures changed; investigate before updating')
    else: path.write_text(text)
    print(f'knot-editing.tsv: {len(text.splitlines())-1} independent exact cases, {len(text)} bytes')


if __name__=='__main__': main()
