#!/usr/bin/env python3
"""Compare rational tensor-product surface jets with independent native OCCT."""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import sys
from compare_occt import ROOT, run


def encode(name,kind,du,dv,poles,weights,uk,um,vk,vm,pu,pv,u,v,su,sv,order):
    nu=sum(um[:-1]) if pu else sum(um)-du-1
    nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
    assert len(poles)==len(weights)==nu*nv
    header=[name,kind,du,dv,nu,nv,len(uk),len(vk),int(pu),int(pv),format(u,'.17g'),format(v,'.17g'),su,sv,order]
    values=[format(x,'.17g') for p,w in zip(poles,weights) for x in [*p,w]]
    values += [x for knots,mults in [(uk,um),(vk,vm)] for k,m in zip(knots,mults) for x in [format(k,'.17g'),str(m)]]
    return ' '.join(map(str,header+values))


def cases():
    rows=[]
    for du,dv in [(1,1),(2,3),(3,2),(5,4),(8,2),(2,8),(25,1),(1,25),(25,25)]:
        for pu,pv in [(False,False),(True,False),(False,True),(True,True)]:
            if du==dv==25 and (pu or pv): continue
            def direction(d,periodic):
                return ([0.,1.,3.],[1,d,1]) if periodic else ([0.,1.,3.],[d+1,min(d,2),d+1])
            uk,um=direction(du,pu)
            vk,vm=direction(dv,pv)
            nu=sum(um[:-1]) if pu else sum(um)-du-1
            nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
            poles=[(float(i-3),float(j-2),float((i*7+j*3+i*j)%13-6)) for i in range(nu) for j in range(nv)]
            weights=[.5*(1+(i+2*j)%3) for i in range(nu) for j in range(nv)]
            queries=[(0.,0.,'R','R'),(3.,3.,'L','L'),(.5,2.,'R','R'),(1.,1.,'L','R'),(1.,1.,'R','L')]
            if pu or pv:
                queries += [(-6. if pu else 0.,9. if pv else 3.,'L' if pu else 'R','R' if pv else 'L')]
            for q,(u,v,su,sv) in enumerate(queries):
                rows.append(encode(f'd{du}_{dv}_p{int(pu)}{int(pv)}_{q}','S',du,dv,poles,weights,uk,um,vk,vm,pu,pv,u,v,su,sv,2))
    # Bezier is a separate OCCT evaluator; a bilinear saddle has nonzero mixed derivative.
    for du,dv in [(1,1),(2,2),(3,5)]:
        poles=[(float(i),float(j),float(i*j)) for i in range(du+1) for j in range(dv+1)]
        for rational in [False,True]:
            weights=[float(1+i%3) if rational else 1. for i in range(len(poles))]
            for q,(u,v) in enumerate([(0.,0.),(.25,.75),(1.,1.)]):
                rows.append(encode(f'bezier_{du}_{dv}_r{int(rational)}_{q}','B',du,dv,poles,weights,[0.,1.],[du+1]*2,[0.,1.],[dv+1]*2,False,False,u,v,'A','A',2))
    # Unclamped directions with different domains and nonuniform knots.
    uk,um=[-3.,-2.,0.,.5,2.,3.,4.,7.],[1]*8
    vk,vm=[-7.,-5.,-2.,0.,1.,3.,4.,8.,9.],[1]*9
    poles=[(float(i),float(j),float(i-j)) for i in range(5) for j in range(5)]
    for q,(u,v) in enumerate([(0.,0.),(.125,.375),(3.,3.)]):
        rows.append(encode(f'unclamped_{q}','S',2,3,poles,[1.+i%3 for i in range(25)],uk,um,vk,vm,False,False,u,v,'A','A',2))
    return '\n'.join(rows)+'\n'


def observations(text):
    rows={}
    for row in text.splitlines():
        label,*values=row.split()
        if label in rows or len(values)!=18: raise ValueError('malformed surface observation')
        values=list(map(float,values))
        if not all(math.isfinite(x) for x in values): raise ValueError('nonfinite surface observation')
        rows[label]=values
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--native-only',action='store_true')
    parser.add_argument('--strict-native',action='store_true',help='Fail on reviewed numerical divergences too')
    args=parser.parse_args()
    prefix=args.occt_root.resolve()
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'Geom_BSplineSurface.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKG3d.*'))),None)
    if not include or not lib: parser.error('OCCT modeling-data SDK not found')
    output=ROOT/'target/surface-oracle'
    output.mkdir(parents=True,exist_ok=True)
    executable=output/'occt_surface_oracle'
    run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_surface_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKG3d','-lTKMath','-lTKernel'],cwd=ROOT)
    data=cases()
    native=run([str(executable)],input=data,cwd=ROOT)
    (output/'inputs.txt').write_text(data)
    (output/'occt.tsv').write_text(native.stdout)
    (output/'native-version.txt').write_text(native.stderr)
    expected=observations(native.stdout)
    if args.native_only:
        print(f'{native.stderr.strip()}: captured {len(expected)} independent surface observations')
        return
    rust=run(['cargo','run','--quiet','--locked','--example','surface_oracle'],input=data,cwd=ROOT)
    (output/'rust.tsv').write_text(rust.stdout)
    actual=observations(rust.stdout)
    if actual.keys()!=expected.keys() or len(expected)!=len(data.splitlines()): raise AssertionError('surface oracle case sets differ')
    largest=0.
    mismatches=[]
    reviewed=[]
    inputs={row.split()[0]:row for row in data.splitlines()}
    for label,values in expected.items():
        from review_spline_jets import reviewed_jet
        fractions=[abs(a-e)/(1e-10+2e-12*abs(e)) for a,e in zip(actual[label],values)]
        largest=max(largest,*fractions)
        if max(fractions)>1.:
            review=reviewed_jet('surfaces',native.stderr.strip(),label,values,actual[label],inputs[label])
            if review:
                reviewed.append({**review,'components':[i for i,f in enumerate(fractions) if f>1.]})
                continue
        for c,(a,e) in enumerate(zip(actual[label],values)):
            fraction=abs(a-e)/(1e-10+2e-12*abs(e))
            largest=max(largest,fraction)
            if fraction>1.: mismatches.append({'case':label,'component':c,'native':e,'rust':a,'error_fraction':fraction})
    report={'oracle':native.stderr.strip(),'cases':len(expected),'matched_cases':len(expected)-len({x['case'] for x in mismatches+reviewed}),
            'reviewed_divergence_cases':len(reviewed),'unexpected_mismatch_cases':len({x['case'] for x in mismatches}),
            'reviewed_divergences':reviewed,
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},'largest_error_fraction_of_budget':largest,
            'domain':'positive-weight Bezier and B-spline surface positions, Du, Dv, Duu, Dvv, Duv; degree 1..25; clamped, unclamped, and independently periodic U/V',
            'deliberate_differences':'Exact knots/weights and parameter wrapping, no tolerance snapping or nonperiodic extrapolation; automatic jets must agree across all requested knot sides.',
            'unexpected_mismatches':mismatches}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if mismatches or (args.strict_native and reviewed): raise SystemExit(1)


if __name__=='__main__':
    import subprocess
    try: main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr)
        raise SystemExit(error.returncode) from error
