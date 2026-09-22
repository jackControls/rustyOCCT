#!/usr/bin/env python3
"""Capture native OCCT Bézier extraction and editing before Rust execution."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from compare_occt import ROOT, run
from compare_splines import cases as spline_cases


def shapes():
    seen=set()
    for row in spline_cases().splitlines():
        w=row.split(); name=w[0].rsplit('_',1)[0]
        if name in seen: continue
        seen.add(name)
        kind=w[1]; degree,np,nk=map(int,w[2:5]); data=list(map(float,w[8:8+4*np])); tail=w[8+4*np:]
        yield name,kind,degree,[tuple(data[i:i+3]) for i in range(0,len(data),4)],data[3::4],list(map(float,tail[::2])),list(map(int,tail[1::2]))
    # Exact source inputs: Geom_BezierCurve_Test RationalIncrease/RationalReverse.
    for label,weights in [('upstream_rational_increase',[1.,2.,1.]),('upstream_rational_reverse',[1.,3.,2.])]:
        yield label,'B',2,[(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)],weights,[0.,1.],[3,3]


def encode(name,kind,degree,poles,weights,knots,mults,first,last,op,elevation):
    w=[name,kind,str(degree),str(len(poles)),str(len(knots)),format(first,'.17g'),format(last,'.17g'),str(op),str(elevation)]
    w += [format(x,'.17g') for p,weight in zip(poles,weights) for x in [*p,weight]]
    w += [x for knot,m in zip(knots,mults) for x in [format(knot,'.17g'),str(m)]]
    return ' '.join(w)


def cases():
    rows=[]
    for name,kind,degree,poles,weights,knots,mults in shapes():
        flat=[k for k,m in zip(knots,mults) for _ in range(m)]
        start,end=(knots[0],knots[-1]) if kind=='P' else (flat[degree],flat[len(poles)])
        windows=[(start,end),(start+(end-start)/8,end-(end-start)/8)]
        if kind=='P': windows.append((start-(end-start)/4,start+(end-start)/2))
        for i,(first,last) in enumerate(windows):
            for op in range(6):
                rows.append(encode(f'{name}_w{i}_op{op}',kind,degree,poles,weights,knots,mults,first,last,op,min(25,degree+3)))
    return '\n'.join(rows)+'\n'


def capture(prefix,output):
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'GeomConvert_BSplineCurveToBezierCurve.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKGeomBase.*'))),None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True,exist_ok=True); executable=output/'occt_bezier_editing_oracle'
    command=shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_bezier_editing_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKGeomBase','-lTKG3d','-lTKMath','-lTKernel']
    try: build=run(command,cwd=ROOT)
    except subprocess.CalledProcessError as e:
        (output/'native-build.log').write_text(e.stdout+e.stderr);raise
    (output/'native-build.log').write_text(build.stdout+build.stderr)
    data=cases(); native=run([str(executable)],input=data,cwd=ROOT)
    for file,text in [('inputs.txt',data),('occt.tsv',native.stdout),('native-version.txt',native.stderr)]: (output/file).write_text(text)
    return data,native


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    p.add_argument('--output',type=Path,default=ROOT/'target/bezier-editing-oracle')
    p.add_argument('--capture-only',action='store_true')
    p.add_argument('--reuse-capture',action='store_true')
    p.add_argument('--strict-native',action='store_true');args=p.parse_args()
    if args.reuse_capture:
        data=(args.output/'inputs.txt').read_text()
        if data!=cases(): raise ValueError('captured corpus changed')
        native=(args.output/'occt.tsv').read_text()
        oracle=(args.output/'native-version.txt').read_text().splitlines()[0]
    else:
        data,result=capture(args.occt_root,args.output)
        native,oracle=result.stdout,result.stderr.splitlines()[0]
    if args.capture_only:
        print(f'{oracle}: {len(data.splitlines())} independent extraction/editing inputs captured')
        return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','bezier_editing_oracle'],input=data,cwd=ROOT).stdout
    (args.output/'rust.tsv').write_text(rust)
    path=ROOT/'rust/fixtures/occt-bezier-editing-divergences.json'
    reviews=json.loads(path.read_text())['reviews'] if path.exists() and not args.strict_native else []
    report=compare(data,native,rust,oracle,reviews)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v if k not in ['failures','reviewed_differences'] else len(v) for k,v in report.items()},indent=2))
    if report['failures']: raise SystemExit(1)


def fingerprint(curves):
    import struct
    def bits(x): return struct.pack('>d',x).hex()
    return [[d,bits(a),bits(b),[[bits(x) for x in c] for c in controls]] for d,a,b,controls in curves]


def compare(data,native,rust,oracle,reviews=()):
    from bezier_editing_reference import decode,encode,expected
    inputs={r.split()[0]:r for r in data.splitlines()}
    if len(inputs)!=len(data.splitlines()): raise ValueError('duplicate input')
    actual,observed=decode(rust),decode(native,native=True)
    if actual.keys()!=inputs.keys() or observed.keys()!=inputs.keys(): raise ValueError('case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,'matched_cases':0,
            'arcs':0,'largest_error_fraction_of_budget':0.,'reviewed_differences':[],'failures':[],
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871'}
    for name,row in inputs.items():
        exact=expected(row)
        if actual[name]!=exact:
            report['failures'].append({'case':name,'reason':'Rust differs from independent exact homogeneous polynomial'})
            continue
        report['independently_verified']+=1; report['arcs']+=len(exact)
        native_curves=observed[name]; differences=[]
        if len(exact)!=len(native_curves): differences.append('arc count')
        for i,((degree,a,b,controls),(nd,na,nb,nc)) in enumerate(zip(exact,native_curves)):
            if degree!=nd: differences.append(f'arc {i} degree'); continue
            expected_values=[float(a),float(b),*[float(x) for c in controls for x in [c[0]/c[3],c[1]/c[3],c[2]/c[3],c[3]/controls[0][3]]]]
            native_values=[na,nb,*[x for c in nc for x in [*c[:3],c[3]/nc[0][3]]]]
            errors=[abs(x-y)/(1e-10+2e-12*abs(y)) for x,y in zip(native_values,expected_values)]
            report['largest_error_fraction_of_budget']=max(report['largest_error_fraction_of_budget'],*errors)
            differences.extend(f'arc {i} field {j}' for j,e in enumerate(errors) if e>1.)
        if not differences:
            report['matched_cases']+=1; continue
        evidence={'case':name,'oracle':oracle,'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                  'native':fingerprint(native_curves),'exact_sha256':hashlib.sha256(encode(name,exact).encode()).hexdigest(),
                  'differences':differences}
        review=next((r for r in reviews if all(r.get(k)==v for k,v in evidence.items())),None)
        if review: report['reviewed_differences'].append(review)
        else: report['failures'].append(evidence)
    return report

if __name__=='__main__':
    try: main()
    except subprocess.CalledProcessError as e:
        print(e.stderr,file=sys.stderr);raise SystemExit(e.returncode) from e
