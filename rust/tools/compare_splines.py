#!/usr/bin/env python3
"""Compare spline positions and first/second derivatives with native OCCT."""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
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
    rows=[]
    for name,kind,degree,poles,weights,knots,mults in shapes:
        flat=[k for k,m in zip(knots,mults) for _ in range(m)]
        start,end=flat[degree],flat[len(poles)]
        queries=[(start,'R'),(end,'L')]+[(start+(end-start)*x,'R') for x in [.125,.375,.625,.875]]
        queries += [(k,side) for k in knots if start<k<end for side in ['L','R']]
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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--native-only',action='store_true')
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
    if actual.keys()!=expected.keys() or len(expected)!=len(data.splitlines()): raise AssertionError('spline oracle case sets differ')
    largest=0.
    for label,values in expected.items():
        for a,e in zip(actual[label],values):
            fraction=abs(a-e)/(1e-10+2e-12*abs(e))
            largest=max(largest,fraction)
            if fraction>1.: raise AssertionError(f'{label}: {a} != {e}')
    report={'oracle':native.stderr.strip(),'cases':len(expected),'largest_error_fraction_of_budget':largest,
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'domain':'positive-weight, nonperiodic Bezier/B-spline positions and first/second derivatives; degree 1..25; clamped/unclamped; explicit one-sided repeated knots',
            'deliberate_differences':'Exact knots/weights, no tolerance snapping or extrapolation. Full-exponent and discontinuity decisions use independent exact oracles.'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    import subprocess
    try: main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr)
        raise SystemExit(error.returncode) from error
