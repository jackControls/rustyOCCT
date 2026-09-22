#!/usr/bin/env python3
"""Capture independent OCCT surface extraction/editing before Rust execution."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import struct
import subprocess
import sys
from compare_occt import ROOT, run
from compare_surfaces import cases as surface_cases


def shapes():
    seen=set()
    for row in surface_cases().splitlines():
        w=row.split(); name=w[0].rsplit('_',1)[0]
        if name in seen: continue
        seen.add(name)
        du,dv,nu,nv,nuk,nvk,pu,pv=map(int,w[2:10])
        data=list(map(float,w[15:15+4*nu*nv])); tail=w[15+4*nu*nv:]
        yield (name,w[1],du,dv,[tuple(data[i:i+3]) for i in range(0,len(data),4)],data[3::4],
               list(map(float,tail[:2*nuk:2])),list(map(int,tail[1:2*nuk:2])),
               list(map(float,tail[2*nuk::2])),list(map(int,tail[2*nuk+1::2])),bool(pu),bool(pv))
    # Complete the degree-25 x degree-25 periodic combinations.
    for pu,pv in [(True,False),(False,True),(True,True)]:
        um=[1,25,1] if pu else [26,2,26]; vm=[1,25,1] if pv else [26,2,26]
        nu=sum(um[:-1]) if pu else sum(um)-26
        nv=sum(vm[:-1]) if pv else sum(vm)-26
        poles=[(float(i-3),float(j-2),float((i*7+j*3+i*j)%13-6)) for i in range(nu) for j in range(nv)]
        weights=[.5*(1+(i+2*j)%3) for i in range(nu) for j in range(nv)]
        yield f'd25_25_p{int(pu)}{int(pv)}','S',25,25,poles,weights,[0.,1.,3.],um,[0.,1.,3.],vm,pu,pv
    # Original Geom_BezierSurface_Test SetUp, RationalSegment, RationalIncrease,
    # RationalSurface_UIso and VIso_Rational input families.
    for rational in [False,True]:
        poles=[(float(i),float(j),(i+j)*.1) for i in range(1,4) for j in range(1,4)]
        weights=[1.+.3*(i+j-2) if rational else 1. for i in range(1,4) for j in range(1,4)]
        yield f'upstream_setup_r{int(rational)}','B',2,2,poles,weights,[0.,1.],[3,3],[0.,1.],[3,3],False,False
    for z in [0.,1.]:
        yield f'upstream_bilinear_z{int(z)}','B',1,1,[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(1.,1.,z)],[1.,2.,2.,1.],[0.,1.],[2,2],[0.,1.],[2,2],False,False
    # Weights varying in only one axis exercise isocurve rational flags.
    for axis in [0,1]:
        poles=[(float(i),float(j),float(i*j)) for i in range(3) for j in range(4)]
        weights=[float(1+(i if axis==0 else j)) for i in range(3) for j in range(4)]
        yield f'one_axis_weights_{axis}','B',2,3,poles,weights,[0.,1.],[3,3],[0.,1.],[4,4],False,False


def encode(shape,window,op,elevation=None,name=None):
    label,kind,du,dv,poles,weights,uk,um,vk,vm,pu,pv=shape
    nu=sum(um[:-1]) if pu else sum(um)-du-1
    nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
    assert len(poles)==len(weights)==nu*nv
    elevation=elevation or (min(25,du+2),min(25,dv+2))
    w=[name or label,kind,du,dv,nu,nv,len(uk),len(vk),int(pu),int(pv),
       *[format(x,'.17g') for x in window],op,*elevation]
    w += [format(x,'.17g') for p,weight in zip(poles,weights) for x in [*p,weight]]
    w += [x for knots,mults in [(uk,um),(vk,vm)] for k,m in zip(knots,mults) for x in [format(k,'.17g'),str(m)]]
    return ' '.join(map(str,w))


def cases():
    rows=[]
    for shape in shapes():
        name,_,du,dv,poles,_,uk,um,vk,vm,pu,pv=shape
        def domain(d,k,m,periodic):
            f=[x for x,n in zip(k,m) for _ in range(n)]
            return (k[0],k[-1]) if periodic else (f[d],f[len(f)-d-1])
        a,b=domain(du,uk,um,pu); c,d=domain(dv,vk,vm,pv)
        windows=[(a,b,c,d),(a+(b-a)/8,b-(b-a)/8,c+(d-c)/4,d-(d-c)/8)]
        if pu or pv:
            windows.append((a-(b-a)/4,a+(b-a)/2,c-(d-c)/4,c+(d-c)/2) if pu and pv else
                           (a-(b-a)/4,a+(b-a)/2,c,d) if pu else (a,b,c-(d-c)/4,c+(d-c)/2))
        for i,window in enumerate(windows):
            for op in range(11) if i==0 else [0,10]:
                rows.append(encode(shape,window,op,name=f'{name}_w{i}_op{op}'))
    return '\n'.join(rows)+'\n'


def capture(prefix,output):
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'GeomConvert_BSplineSurfaceToBezierSurface.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKGeomBase.*'))),None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True,exist_ok=True); executable=output/'occt_surface_editing_oracle'
    command=shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_surface_editing_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKGeomBase','-lTKG3d','-lTKMath','-lTKernel']
    try: build=run(command,cwd=ROOT)
    except subprocess.CalledProcessError as e:
        (output/'native-build.log').write_text(e.stdout+e.stderr); raise
    (output/'native-build.log').write_text(build.stdout+build.stderr)
    data=cases(); native=run([str(executable)],input=data,cwd=ROOT)
    for file,text in [('inputs.txt',data),('occt.tsv',native.stdout),('native-version.txt',native.stderr)]: (output/file).write_text(text)
    return data,native


def fingerprint(items):
    if items is None: return None
    bits=[[*item[:3],*[struct.pack('>d',x).hex() for x in item[3]],
             [[struct.pack('>d',x).hex() for x in p] for p in item[4]]] for item in items]
    return hashlib.sha256(json.dumps(bits,separators=(',',':')).encode()).hexdigest()


def field_description(item,indices):
    ranges=[]
    for i in indices:
        if ranges and i==ranges[-1][1]+1: ranges[-1][1]=i
        else: ranges.append([i,i])
    return f'item {item} fields '+','.join(str(a) if a==b else f'{a}-{b}' for a,b in ranges)


def compare(data,native,rust,oracle,reviews=()):
    from surface_editing_reference import decode,encode_result,expected
    inputs={r.split()[0]:r for r in data.splitlines()}
    if len(inputs)!=len(data.splitlines()): raise ValueError('duplicate input')
    actual,observed=decode(rust),decode(native,native=True)
    if actual.keys()!=inputs.keys() or observed.keys()!=inputs.keys(): raise ValueError('case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,'matched_cases':0,
            'patches':0,'isocurves':0,'largest_error_fraction_of_budget':0.,
            'reviewed_differences':[],'failures':[],
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871'}
    for name,row in inputs.items():
        exact=expected(row)
        if actual[name]!=exact:
            report['failures'].append({'case':name,'reason':'Rust differs from independent exact tensor polynomial'})
            continue
        report['independently_verified']+=1
        report['patches']+=sum(x[0]=='P' for x in exact); report['isocurves']+=sum(x[0]=='C' for x in exact)
        observed_items=observed[name]; differences=[]
        if observed_items is None: differences.append('native exception')
        else:
            if len(exact)!=len(observed_items): differences.append('item count')
            for i,(e,n) in enumerate(zip(exact,observed_items)):
                if e[:3]!=n[:3]: differences.append(f'item {i} kind/degrees'); continue
                if any(not math.isfinite(x) for x in [*n[3],*[x for p in n[4] for x in p]]) or any(p[3]<=0 for p in n[4]):
                    differences.append(f'item {i} invalid native values'); continue
                ev=[*map(float,e[3]),*[float(x) for p in e[4] for x in [p[0]/p[3],p[1]/p[3],p[2]/p[3],p[3]/e[4][0][3]]]]
                nv=[*n[3],*[x for p in n[4] for x in [*p[:3],p[3]/n[4][0][3]]]]
                errors=[abs(x-y)/(1e-10+2e-12*abs(y)) for x,y in zip(nv,ev)]
                report['largest_error_fraction_of_budget']=max(report['largest_error_fraction_of_budget'],*errors)
                outside=[j for j,error in enumerate(errors) if error>1.]
                if outside: differences.append(field_description(i,outside))
        if not differences:
            report['matched_cases']+=1; continue
        evidence={'case':name,'oracle':oracle,'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                  'native':fingerprint(observed_items),'exact_sha256':hashlib.sha256(encode_result(name,exact).encode()).hexdigest(),
                  'differences':differences}
        review=next((r for r in reviews if all(r.get(k)==v for k,v in evidence.items())),None)
        if review: report['reviewed_differences'].append(review)
        else: report['failures'].append(evidence)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    p.add_argument('--output',type=Path,default=ROOT/'target/surface-editing-oracle')
    p.add_argument('--capture-only',action='store_true'); p.add_argument('--reuse-capture',action='store_true')
    p.add_argument('--strict-native',action='store_true'); args=p.parse_args()
    if args.reuse_capture:
        data=(args.output/'inputs.txt').read_text()
        if data!=cases(): raise ValueError('captured corpus changed')
        native=(args.output/'occt.tsv').read_text(); oracle=(args.output/'native-version.txt').read_text().splitlines()[0]
    else:
        data,result=capture(args.occt_root,args.output); native,oracle=result.stdout,result.stderr.splitlines()[0]
    if args.capture_only:
        print(f'{oracle}: {len(data.splitlines())} independent surface-editing inputs captured'); return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','surface_editing_oracle'],input=data,cwd=ROOT).stdout
    (args.output/'rust.tsv').write_text(rust)
    path=ROOT/'rust/fixtures/occt-surface-editing-divergences.json'
    reviews=json.loads(path.read_text())['reviews'] if path.exists() and not args.strict_native else []
    report=compare(data,native,rust,oracle,reviews)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v if k not in ['failures','reviewed_differences'] else len(v) for k,v in report.items()},indent=2))
    if report['failures']: raise SystemExit(1)


if __name__=='__main__':
    try: main()
    except subprocess.CalledProcessError as e:
        print(e.stderr,file=sys.stderr); raise SystemExit(e.returncode) from e
