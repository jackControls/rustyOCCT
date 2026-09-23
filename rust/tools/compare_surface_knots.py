#!/usr/bin/env python3
"""Compare complete surface knot grids with OCCT and independent equations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from compare_occt import ROOT, run
from compare_surface_editing import shapes as existing_shapes


def shapes():
    yield from existing_shapes()
    # Every supported degree appears in each edited axis, with unclamped raw
    # support and independently periodic transverse directions.
    for degree in range(1,26):
        for edited_axis in range(2):
            degrees=[degree,2] if edited_axis==0 else [2,degree]
            knots=[list(map(float,range(-2,5))),[0.,1.,3.]]
            mults=[[degree,degree,1,degree,1,degree,degree],[1,2,1]]
            periodic=[False,True]
            if edited_axis: knots.reverse(); mults.reverse(); periodic.reverse()
            counts=[sum(m)-(m[0] if p else d+1) for d,m,p in zip(degrees,mults,periodic)]
            poles=[(float(i-2),float(j-1),float((3*i+7*j+i*j)%13-5)) for i in range(counts[0]) for j in range(counts[1])]
            weights=[.5*(1+(i+2*j)%3) for i in range(counts[0]) for j in range(counts[1])]
            yield (f'raw_support_d{degree}_axis{edited_axis}','S',*degrees,poles,weights,*[x for k,m in zip(knots,mults) for x in [k,m]],*periodic)
            # A nonconstant surface whose edited direction is redundant. Origin
            # removal must keep the other parameter and its geometry unchanged.
            knots=[list(map(float,range(degree+4))),[0.,1.]]
            mults=[[1]*(degree+4),[3,3]]; periodic=[True,False]
            if edited_axis: knots.reverse(); mults.reverse(); periodic.reverse()
            counts=[sum(m)-(m[0] if p else d+1) for d,m,p in zip(degrees,mults,periodic)]
            poles=[]; weights=[]
            for i in range(counts[0]):
                for j in range(counts[1]):
                    k=j if edited_axis==0 else i
                    poles.append((float(k),float(k*k),float(2*k-3))); weights.append(float(1+k))
            yield (f'redundant_d{degree}_axis{edited_axis}','S',*degrees,poles,weights,*[x for k,m in zip(knots,mults) for x in [k,m]],*periodic)


def encode(shape,suffix,operations):
    name,_,du,dv,poles,weights,uk,um,vk,vm,pu,pv=shape
    nu=sum(um)-(um[0] if pu else du+1); nv=sum(vm)-(vm[0] if pv else dv+1)
    assert len(poles)==len(weights)==nu*nv<=4096
    words=[name+'_'+suffix,du,dv,nu,nv,len(uk),len(vk),int(pu),int(pv),len(operations)]
    words += [format(x,'.17g') for point,w in zip(poles,weights) for x in [*point,w]]
    words += [x for knots,mults in [(uk,um),(vk,vm)] for k,m in zip(knots,mults) for x in [format(k,'.17g'),m]]
    words += [x for op,axis,u,m in operations for x in [op,axis,format(u,'.17g'),m]]
    return ' '.join(map(str,words))


def cases():
    rows=[]
    for shape in shapes():
        name,_,du,dv,_,_,uk,um,vk,vm,pu,pv=shape
        axes=[]
        for degree,knots,mults,periodic in [(du,uk,um,pu),(dv,vk,vm,pv)]:
            flat=[k for k,m in zip(knots,mults) for _ in range(m)]
            a,b=(knots[0],knots[-1]) if periodic else (flat[degree],flat[-degree-1])
            cut=(3*a+next(k for k in knots if k>a))/4
            axes.append((degree,knots,mults,periodic,a,b,cut))
        u,v=axes[0][-1],axes[1][-1]
        rows.append(encode(shape,'identity',[]))
        rows.append(encode(shape,'both',[('I','U',u,du),('I','V',v,dv)]))
        rows.append(encode(shape,'sequence',[('I','V',v,1),('I','U',u,1),('I','V',v,dv),('I','U',u,du),('I','V',v,1)]))
        rows.append(encode(shape,'both_roundtrip',[('I','U',u,du),('I','V',v,dv),('R','U',u,0),('R','V',v,0)]))
        for axis,(degree,knots,mults,periodic,a,b,cut) in zip('UV',axes):
            rows.append(encode(shape,axis+'_insert',[('I',axis,cut,degree)]))
            rows.append(encode(shape,axis+'_roundtrip',[('I',axis,cut,1),('R',axis,cut,0)]))
            interior=next((k for k in knots if a<k<b),None)
            count=sum(mults)-(mults[0] if periodic else degree+1)
            if interior is not None and count>degree+1:
                rows.append(encode(shape,axis+'_remove_existing',[('R',axis,interior,mults[knots.index(interior)]-1)]))
            if periodic:
                rows.append(encode(shape,axis+'_seam',[('I',axis,b,degree),('I',axis,a,degree)]))
                if mults[0]<degree:
                    rows.append(encode(shape,axis+'_seam_roundtrip',[('I',axis,b,degree),('R',axis,a,mults[0])]))
                if name.startswith('redundant_'):
                    rows.append(encode(shape,axis+'_remove_origin',[('R',axis,a,0)]))
                    rows.append(encode(shape,axis+'_remove_end_alias',[('R',axis,b,0)]))
            elif knots[0]<a:
                rows.append(encode(shape,axis+'_domain_ends',[('I',axis,a,degree),('I',axis,b,degree)]))
        if name=='upstream_setup_r0':
            for axis in 'UV':
                rows.append(encode(shape,'bspline_gtest_'+axis+'_insert',[('I',axis,.5,1)]))
                rows.append(encode(shape,'bspline_gtest_'+axis+'_remove',[('I',axis,.5,1),('R',axis,.5,0)]))
    return '\n'.join(rows)+'\n'


def capture(prefix,output):
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'Geom_BSplineSurface.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKG3d.*'))),None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True,exist_ok=True); executable=output/'occt_surface_knot_oracle'
    command=shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_surface_knot_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKG3d','-lTKMath','-lTKernel']
    try: build=run(command,cwd=ROOT)
    except subprocess.CalledProcessError as e:
        (output/'native-build.log').write_text(e.stdout+e.stderr); raise
    (output/'native-build.log').write_text(build.stdout+build.stderr)
    data=cases(); native=run([str(executable)],input=data,cwd=ROOT)
    for name,text in [('inputs.txt',data),('occt.tsv',native.stdout),('native-version.txt',native.stderr)]: (output/name).write_text(text)
    record={'head':run(['git','rev-parse','HEAD'],cwd=ROOT).stdout.strip(),
            'implementation_present':(ROOT/'rust/kernel/src/surface/knot_editing.rs').exists(),
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'cases':len(data.splitlines()),'input_sha256':hashlib.sha256(data.encode()).hexdigest(),
            'native_sha256':hashlib.sha256(native.stdout.encode()).hexdigest(),'oracle':native.stderr.splitlines()[0]}
    (output/'capture.json').write_text(json.dumps(record,indent=2)+'\n')
    return data,native.stdout,record


def fingerprint(observation):
    import struct
    flags,(axes,controls),bounds=observation
    bits=lambda x:struct.pack('>d',x).hex()
    return [list(map(int,flags)),[[p,periodic,list(map(bits,k)),list(m)] for p,periodic,k,m in axes],[[bits(x) for x in b] for b in bounds],[[bits(x) for x in c] for c in controls]]


def compare(data,native,rust,exact,oracle,reviews=()):
    from surface_knot_reference import decode,encode as result
    inputs={r.split()[0]:r for r in data.splitlines()}
    if len(inputs)!=len(data.splitlines()): raise ValueError('duplicate input')
    actual,observed,expected=decode(rust),decode(native,True),decode(exact)
    if actual.keys()!=inputs.keys() or observed.keys()!=inputs.keys() or expected.keys()!=inputs.keys(): raise ValueError('case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,'matched_cases':0,
            'largest_error_fraction_of_budget':0.,'reviewed_differences':[],'failures':[],
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'removal_contract':'Exact positive homogeneous equality across the complete grid; native tolerance 1e-9.'}
    for name,row in inputs.items():
        if actual[name]!=expected[name]:
            report['failures'].append({'case':name,'reason':'Rust differs from independent coefficient equations'}); continue
        report['independently_verified']+=1
        flags,(axes,controls),bounds=expected[name]; nf,(na,nc),nb=observed[name]
        differences=[]
        if flags!=nf: differences.append('operation success flags')
        if axes!=na: differences.append('axes, knots or multiplicities')
        if bounds!=nb: differences.append('parameter domains')
        if len(controls)!=len(nc): differences.append('control grid count')
        if not differences:
            values=[float(x) for c in controls for x in [c[0]/c[3],c[1]/c[3],c[2]/c[3],c[3]]]
            errors=[abs(a-b)/(1e-10+2e-12*abs(a)) for a,b in zip(values,[x for c in nc for x in c])]
            report['largest_error_fraction_of_budget']=max(report['largest_error_fraction_of_budget'],*errors)
            differences.extend(f'control field {i}' for i,e in enumerate(errors) if e>1.)
        if not differences:
            report['matched_cases']+=1; continue
        canonical=result(name,flags,(axes,controls),compact=True)
        evidence={'case':name,'oracle':oracle,'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                  'native':fingerprint(observed[name]),'exact_sha256':hashlib.sha256(canonical.encode()).hexdigest(),
                  'differences':differences}
        review=next((r for r in reviews if all(r.get(k)==v for k,v in evidence.items())),None)
        if review: report['reviewed_differences'].append(review)
        else: report['failures'].append(evidence)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--output',type=Path,default=ROOT/'target/surface-knot-oracle')
    parser.add_argument('--capture-only',action='store_true')
    parser.add_argument('--reuse-capture',action='store_true')
    parser.add_argument('--strict-native',action='store_true'); args=parser.parse_args()
    if args.reuse_capture:
        data=(args.output/'inputs.txt').read_text()
        if data!=cases(): raise ValueError('captured corpus changed')
        native=(args.output/'occt.tsv').read_text(); oracle=(args.output/'native-version.txt').read_text().splitlines()[0]
    else:
        data,native,record=capture(args.occt_root,args.output); oracle=record['oracle']
    if args.capture_only:
        print(f'{oracle}: {len(data.splitlines())} complete surface knot observations captured'); return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','surface_knot_oracle'],input=data,cwd=ROOT).stdout
    (args.output/'rust.tsv').write_text(rust)
    path=ROOT/'rust/fixtures/occt-surface-knot-divergences.json'
    reviews=json.loads(path.read_text())['reviews'] if path.exists() and not args.strict_native else []
    fixtures={}
    for row in (ROOT/'rust/fixtures/surface-knots.tsv').read_text().splitlines():
        if row.startswith('#'): continue
        source,result=row.split(' | ')
        if source in fixtures: raise ValueError('duplicate fixture')
        fixtures[source]=result
    exact='\n'.join(fixtures[row] for row in data.splitlines())+'\n'
    report=compare(data,native,rust,exact,oracle,reviews)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:len(v) if k in ['failures','reviewed_differences'] else v for k,v in report.items()},indent=2))
    if report['failures']: raise SystemExit(1)


if __name__=='__main__':
    try: main()
    except subprocess.CalledProcessError as e:
        print(e.stderr,file=sys.stderr); raise SystemExit(e.returncode) from e
