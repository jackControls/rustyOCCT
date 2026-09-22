#!/usr/bin/env python3
"""Independent exact basis/continued-fraction oracle for spline quadrics."""
import argparse
import random
from fractions import Fraction as F
from pathlib import Path
from compare_spline_quadric import cases
from compare_spline_plane import encode
from generate_spline_plane_fixtures import expected as spline_expected
from exact_polynomial_oracle import add,multiply


def expected(row):
    primitive=row.split()[1].split(':')[0]
    def implicit(h,shape):
        center,axis,r=shape; radius=r[0]
        delta=[add(h[c],[-center[c]*w for w in h[3]]) for c in range(3)]
        if primitive=='cylinder':
            # Independent cross-product form, versus production dot projection.
            delta=[add([axis[b]*x for x in delta[a]],[-axis[a]*x for x in delta[b]]) for a,b in [(1,2),(2,0),(0,1)]]
            radius2=radius*radius*sum(x*x for x in axis)
        elif primitive=='sphere': radius2=radius*radius
        else: raise ValueError('invalid quadric')
        result=[-radius2*x for x in multiply(h[3],h[3])]
        for p in delta: result=add(result,multiply(p,p))
        return result
    return spline_expected(row,implicit)


def inputs():
    rows=cases().splitlines()
    # Full-exponent inputs remain exact even if squared intermediates overflow
    # or underflow. These are not claims about OCCT's tolerance-based domain.
    for primitive in ['sphere','cylinder']:
        for name,radius,axis in [('tiny',float.fromhex('0x0.0000000000001p-1022'),[0.,0.,1.]),('huge',2.**1022,[0.,0.,1.]),('axis_tiny',1.,[0.,0.,2.**-1074]),('axis_huge',1.,[0.,0.,2.**1023])]:
            poles=[(-2*radius,0.,0.),(2*radius,0.,0.)]
            rows.append(encode(f'{primitive}_{name}',primitive+':B',1,[[0.,0.,0.],axis,[radius,0.,0.]],poles,[1.,1.],[0.,1.],[2,2]))
        # An order-50 tangency at a query endpoint: x=t^25, y=1.
        rows.append(encode(primitive+'_order50',primitive+':B',25,[[0.,0.,0.],[0.,0.,1.],[1.,0.,0.]],[(0.,1.,0.)]*25+[(1.,1.,0.)],[1.]*26,[0.,1.],[26,26]))
        rows.append(encode(primitive+'_order50_interior',primitive+':B',25,[[0.,0.,0.],[0.,0.,1.],[1.,0.,0.]],[((-1.)**i,1.,0.) for i in range(26)],[1.]*26,[0.,1.],[26,26]))
        rows.append(encode(primitive+'_subulp',primitive+':PT',1,[[0.,0.,0.],[0.,0.,1.],[1.,0.,0.]],[(-2.,0.,0.),(2.,0.,0.)],[1.,1.],[0.,.25,.5],[1]*3)+f' {2.**53:.17g} {2.**53+2:.17g}')
        rng=random.Random(73891)
        for d in [2,3,4,5,8,25]:
            poles=[tuple(float(rng.randint(-3,3)) for _ in range(3)) for _ in range(d+1)]
            poles[0]=(-3.,0.,0.);poles[-1]=(3.,0.,0.)
            weights=[float(rng.randint(1,3)) for _ in poles]
            rows.append(encode(f'{primitive}_dense_{d}',primitive+':B',d,[[0.,0.,0.],[1.,2.,3.],[1.,0.,0.]],poles,weights,[0.,1.],[d+1]*2))
    return rows


def generate():
    return '# input text | point_count overlap_count [t/x/y/z lo hi bits contact left_order right_order] [overlap endpoints bits]\n'+'\n'.join(row+' | '+' '.join(expected(row)) for row in inputs())+'\n'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args()
    text=generate();path=Path(__file__).resolve().parents[1]/'fixtures/spline-quadric.tsv'
    if args.check:
        if path.read_text()!=text:p.error('spline-quadric fixtures changed; investigate before updating')
    else:path.write_text(text)
    print(f'spline-quadric.tsv: {len(text.splitlines())-1} independent exact cases')


if __name__=='__main__':main()
