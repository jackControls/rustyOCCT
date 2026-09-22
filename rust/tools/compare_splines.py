#!/usr/bin/env python3
"""Compare spline positions and first/second derivatives with native OCCT."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import struct
import sys
from compare_occt import ROOT, run


def cases():
    shapes=[]
    # Original Geom_BSplineCurve_Test.cxx SetUp poles and clamped cubic knots.
    shapes.append(('upstream_cubic','S',3,[(0.,0.,0.),(1.,1.,0.),(2.,1.,0.),(3.,0.,0.)],[1.]*4,[0.,1.],[4,4]))
    shapes.append(('quarter_conic','B',2,[(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)],[1.,math.sqrt(.5),1.],[0.,1.],[3,3]))
    for degree in [1,2,3,5,8,25]:
        for rational in [False,True]:
            for spans in [1,3]:
                knots=[float(i) for i in range(spans+1)]
                mults=[degree+1]+[degree if i==1 else 1 for i in range(1,spans)]+[degree+1]
                count=sum(mults)-degree-1
                poles=[(float(i-3),float((i*i)%11-5),float((i*7)%13-6)) for i in range(count)]
                weights=[float(1+i%3)/2 if rational else 1. for i in range(count)]
                kind='B' if spans==1 else 'S'
                shapes.append((f'd{degree}_r{int(rational)}_s{spans}',kind,degree,poles,weights,knots,mults))
    # Non-clamped, nonuniform domain [U[p], U[n+1]], not the extreme knots.
    shapes.append(('unclamped','S',2,[(float(i),float(i%2),float(-i)) for i in range(5)],[1.,2.,1.,.5,1.],[-3.,-2.,0.,.5,2.,3.,4.,7.],[1]*8))
    for degree in [1,2,3,5,8,25]:
        for seam in sorted({1,degree}):
            knots=[0.,.5,2.,3.]
            mults=[seam,degree,1,seam]
            count=sum(mults[:-1])
            poles=[(float(i-3),float(i*i%11-5),float(i*7%13-6)) for i in range(count)]
            weights=[float(1+i%3)/2 for i in range(count)]
            shapes.append((f'periodic_d{degree}_m{seam}','P',degree,poles,weights,knots,mults))
    rows=[]
    for name,kind,degree,poles,weights,knots,mults in shapes:
        flat=[k for k,m in zip(knots,mults) for _ in range(m)]
        start,end=(knots[0],knots[-1]) if kind=='P' else (flat[degree],flat[len(poles)])
        queries=[(start,'R'),(end,'L')]+[(start+(end-start)*x,'R') for x in [.125,.375,.625,.875]]
        queries += [(k,side) for k in knots if start<k<end for side in ['L','R']]
        if kind=='P':
            queries += [(u,side) for u in [start,end,-6.,9.,-5.5,9.5] for side in ['L','R']]
        for i,(u,side) in enumerate(queries):
            values=[name+'_'+str(i),kind,str(degree),str(len(poles)),str(len(knots)),format(u,'.17g'),side,'2']
            values += [format(x,'.17g') for p,w in zip(poles,weights) for x in [*p,w]]
            values += [x for k,m in zip(knots,mults) for x in [format(k,'.17g'),str(m)]]
            rows.append(' '.join(values))
    return '\n'.join(rows)+'\n'


def observations(text):
    rows={}
    for row in text.splitlines():
        label,*values=row.split()
        if label in rows or len(values)!=9: raise ValueError('malformed spline observation')
        values=list(map(float,values))
        if not all(math.isfinite(x) for x in values): raise ValueError('nonfinite spline observation')
        rows[label]=values
    return rows


def reviewed_divergence(oracle, label, component, native, rust, input_row):
    """A review never exempts Rust from the independent exact reference.

    Match the version, input hash and exact observed native bits, then recompute
    the mathematical value from basis functions. Any changed disagreement needs
    a new investigation; there is no broader epsilon or per-label bypass.
    """
    reviews=json.loads((ROOT/'rust/fixtures/occt-spline-divergences.json').read_text())
    for review in reviews:
        if (review['oracle'],review['case'],review['component'])!=(oracle,label,component): continue
        if review['input_sha256']!=hashlib.sha256(input_row.encode()).hexdigest(): continue
        if review['native_bits']!=struct.pack('>d',native).hex(): continue
        from generate_spline_fixtures import jet
        from generate_curved_fixtures import enclosure
        from generate_spatial_fixtures import value, sign
        words=input_row.split()
        degree,np,nk=map(int,words[2:5])
        data=list(map(lambda x:F(float(x)),words[8:8+4*np]))
        tail=words[8+4*np:]
        if len(tail)!=2*nk: return None
        flat=[F(float(k)) for k,m in zip(tail[::2],tail[1::2]) for _ in range(int(m))]
        u=F(float(words[5]))
        left=words[6]=='L' or u==flat[np]
        span=max(i for i,k in enumerate(flat) if (k<u if left else k<=u))
        exact=jet(degree,[data[i:i+3] for i in range(0,len(data),4)],data[3::4],flat,u,span,2)[component//3][component%3]
        if exact!=F(int(review['exact_numerator']),int(review['exact_denominator'])): return None
        lower,upper=[value(int(x,16)) for x in enclosure(lambda x:sign(exact-x))]
        if rust not in [lower,upper]: return None
        return {'id':review['id'],'case':label,'component':component,'native':native,'rust':rust,
                'exact_lower':lower,'exact_upper':upper,'input_sha256':review['input_sha256'],
                'reason':review['reason'],'evidence':review['evidence']}
    return None


def compare_observations(oracle, inputs, expected, actual):
    if actual.keys()!=expected.keys() or expected.keys()!=inputs.keys(): raise AssertionError('spline oracle case sets differ')
    largest=0.
    reviewed=[]
    unexpected=[]
    for label,values in expected.items():
        from review_spline_jets import reviewed_jet
        fractions=[abs(a-e)/(1e-10+2e-12*abs(e)) for a,e in zip(actual[label],values)]
        largest=max(largest,*fractions)
        if max(fractions)>1.:
            review=reviewed_jet('curves',oracle,label,values,actual[label],inputs[label])
            if review:
                reviewed.append({**review,'components':[i for i,f in enumerate(fractions) if f>1.]})
                continue
        for component,(a,e) in enumerate(zip(actual[label],values)):
            fraction=abs(a-e)/(1e-10+2e-12*abs(e))
            largest=max(largest,fraction)
            if fraction>1.:
                review=reviewed_divergence(oracle,label,component,e,a,inputs[label])
                if review: reviewed.append(review)
                else: unexpected.append({'case':label,'component':component,'native':e,'rust':a,'error_fraction':fraction})
    divergent={x['case'] for x in reviewed+unexpected}
    return {'matched_cases':len(expected)-len(divergent),
            'reviewed_divergence_cases':len({x['case'] for x in reviewed}),
            'unexpected_mismatch_cases':len({x['case'] for x in unexpected}),
            'largest_error_fraction_of_budget':largest,
            'reviewed_divergences':reviewed,'unexpected_mismatches':unexpected}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--native-only',action='store_true')
    parser.add_argument('--strict-native',action='store_true',help='Fail on reviewed native numerical divergences as well as unexpected mismatches')
    args=parser.parse_args()
    prefix=args.occt_root.resolve()
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'Geom_BSplineCurve.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKG3d.*'))),None)
    if not include or not lib: parser.error('OCCT modeling-data SDK not found')
    output=ROOT/'target/spline-oracle'
    output.mkdir(parents=True,exist_ok=True)
    executable=output/'occt_spline_oracle'
    run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_spline_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKG3d','-lTKMath','-lTKernel'],cwd=ROOT)
    data=cases()
    native=run([str(executable)],input=data,cwd=ROOT)
    (output/'inputs.txt').write_text(data)
    (output/'occt.tsv').write_text(native.stdout)
    (output/'native-version.txt').write_text(native.stderr)
    expected=observations(native.stdout)
    if args.native_only:
        print(f'{native.stderr.strip()}: captured {len(expected)} independent spline observations')
        return
    rust=run(['cargo','run','--quiet','--locked','--example','spline_oracle'],input=data,cwd=ROOT)
    (output/'rust.tsv').write_text(rust.stdout)
    actual=observations(rust.stdout)
    comparison=compare_observations(native.stderr.strip(),{row.split()[0]:row for row in data.splitlines()},expected,actual)
    report={'oracle':native.stderr.strip(),'cases':len(expected),**comparison,
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'domain':'positive-weight Bezier/B-spline positions and first/second derivatives; degree 1..25; clamped/unclamped/periodic; explicit one-sided repeated knots and seams',
            'deliberate_differences':'Exact knots/weights, no tolerance snapping or extrapolation. Full-exponent and discontinuity decisions use independent exact oracles.'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if comparison['unexpected_mismatches'] or (args.strict_native and comparison['reviewed_divergences']):
        raise SystemExit(1)


if __name__=='__main__':
    import subprocess
    try: main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr)
        raise SystemExit(error.returncode) from error
