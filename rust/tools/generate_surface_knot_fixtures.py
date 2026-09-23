#!/usr/bin/env python3
"""Generate full surface knot controls with independent coefficient equations."""
import argparse
import math
from pathlib import Path
from compare_surface_knots import cases, encode
from surface_knot_reference import expected, encode as result
ROOT=Path(__file__).resolve().parents[1]


def inputs():
    yield from cases().splitlines()
    tiny=float.fromhex('0x0.0000000000001p-1022')
    huge=float.fromhex('0x1.fffffffffffffp+1023')
    for name,lo,hi,cut,scale,weights in [
        ('tiny_domain',0.,4*tiny,tiny,1.,[1.,2.,3.]),
        ('huge_domain',-huge,huge,0.,1.,[1.,2.,3.]),
        ('adjacent_cut',1.,2.,math.nextafter(1.,2.),1.,[1.,2.,3.]),
        ('huge_homogeneous',0.,1.,.25,huge,[tiny,huge,1.]),
    ]:
        poles=[((i-1)*scale,(j-1)*scale,((i+j)%3-1)*scale) for i in range(3) for j in range(3)]
        shape=(name,'S',2,2,poles,[weights[(i+j)%3] for i in range(3) for j in range(3)],[lo,hi],[3,3],[lo,hi],[3,3],False,False)
        yield encode(shape,'both',[('I','U',cut,2),('I','V',cut,2)])
        yield encode(shape,'roundtrip',[('I','V',cut,2),('I','U',cut,2),('R','V',cut,0),('R','U',cut,0)])


def generate():
    rows=['# input | exact output: flags, axes, full homogeneous grid with a shared denominator']
    for i,row in enumerate(inputs()):
        flags,surface=expected(row)
        rows.append(row+' | '+result(row.split()[0],flags,surface,compact=True))
        if (i+1)%100==0: print(f'{i+1} independent surface knot cases',flush=True)
    return '\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    text=generate(); path=ROOT/'fixtures/surface-knots.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('independent surface knot fixtures changed')
    else: path.write_text(text)
    print(f'{len(text.splitlines())-1} exact surface knot fixtures, {len(text)} bytes',flush=True)


if __name__=='__main__': main()
