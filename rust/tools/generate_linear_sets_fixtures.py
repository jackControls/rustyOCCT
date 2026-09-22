#!/usr/bin/env python3
"""Reproduce complete linear sets with the independent boundary oracle."""
import argparse
import math
from pathlib import Path
from compare_linear_sets import cases
from generate_proximity_fixtures import inputs as proximity_inputs
from linear_sets_reference import intersection, encode
from proximity_reference import bits, exact, parse, validate

OUTPUT=Path(__file__).resolve().parents[1]/'fixtures/linear_sets.tsv'

def inputs():
    for row in cases().splitlines():
        name,shapes,rest=parse(row)
        assert not rest
        yield name,shapes
    native_names={x.split()[0] for x in cases().splitlines()}
    for name,shapes in proximity_inputs():
        # Full exponent, invalid, near-degenerate and overflow queries retained.
        if name.startswith('random_') or name not in native_names:
            yield 'adversarial_'+name,shapes
    # A genuine five-vertex overlap, at normal, huge and subnormal scales.
    for exponent in [0,-1000,990,-1073]:
        for frame in range(3):
            def transform(p):
                p=p if frame==0 else (p[2],-p[0],p[1]) if frame==1 else (2*p[0]+p[2],p[0]+2*p[1]-2*p[2],p[1]+4*p[2])
                return tuple(math.ldexp(x,exponent) for x in p)
            a=[(0.,0.,0.),(6.,0.,0.),(0.,6.,0.)]
            b=[(-1.,2.,0.),(2.,-1.,0.),(2.,2.,0.)]
            yield f'five_vertices_{exponent}_{frame}', [('T',list(map(transform,a))),('T',list(map(transform,b)))]

def generate():
    rows=['# Binary64 input bits; independent exact boundary intersections and gift-wrapping hull.',
          '# R followed by E/P/S/G/L/F set; N nonfinite; D degenerate. Plane: normal, scalar offset.']
    for name,shapes in inputs():
        status=next((s for s in map(validate,shapes) if s!='R'),'R')
        result=[status]
        if status=='R': result.append(encode(intersection(*map(exact,shapes))))
        row=[name]
        for k,points in shapes: row += [k,*[bits(x) for p in points for x in p]]
        rows.append(' '.join([*row,*result]))
    return '\n'.join(rows)+'\n'

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    content=generate()
    if args.check:
        if OUTPUT.read_text()!=content: raise SystemExit('linear set fixtures differ; investigate before updating')
    else: OUTPUT.write_text(content)
    print(f'{len(content.splitlines())-2} independent complete intersection fixtures verified')
