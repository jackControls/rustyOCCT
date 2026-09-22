#!/usr/bin/env python3
"""Native OCCT observations for complete linear-primitive distance queries."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sys
from compare_occt import ROOT, run
from proximity_reference import bits, certify, distance2, enclosure, exact, parse


def encode(name,a,b):
    return ' '.join([name,*[str(x) for kind,points in [a,b]
                            for x in [kind,*[format(v,'.17g') for p in points for v in p]]]])


def cases():
    left={'P':[(.5,1.25,-.5)], 'L':[(-2.,-1.,0.),(3.,2.,1.)],
          'S':[(-2.,-1.,0.),(3.,2.,1.)],
          'F':[(0.,0.,0.),(2.,0.,0.),(0.,2.,0.)],
          'T':[(0.,0.,0.),(2.,0.,0.),(0.,2.,0.)]}
    right={'P':[(1.,.5,3.)], 'L':[(-1.,3.,2.),(2.,-2.,-1.)],
           'S':[(-1.,3.,2.),(2.,-2.,-1.)],
           'F':[(.5,.5,1.),(2.,.5,-1.),(.5,3.,2.)],
           'T':[(.5,.5,1.),(2.,.5,-1.),(.5,3.,2.)]}
    base=[(a+b,(a,left[a]),(b,right[b])) for a in left for b in right]
    tri=('T',[(0.,0.,0.),(4.,0.,0.),(0.,4.,0.)])
    base += [
        ('parallel_lines',('L',[(0.,0.,0.),(2.,1.,0.)]),('L',[(0.,0.,3.),(4.,2.,3.)])),
        ('coincident_lines',('L',[(0.,0.,0.),(2.,1.,0.)]),('L',[(4.,2.,0.),(-2.,-1.,0.)])),
        ('parallel_planes',('F',tri[1]),('F',[(0.,0.,3.),(2.,1.,3.),(-1.,2.,3.)])),
        ('coincident_planes',('F',tri[1]),('F',[(1.,1.,0.),(2.,1.,0.),(1.,2.,0.)])),
        ('parallel_line_plane',('L',[(0.,0.,3.),(2.,1.,3.)]),('F',tri[1])),
        ('parallel_segments',('S',[(0.,0.,0.),(4.,0.,0.)]),('S',[(1.,1.,0.),(3.,1.,0.)])),
        ('overlap_segments',('S',[(0.,0.,0.),(4.,0.,0.)]),('S',[(3.,0.,0.),(1.,0.,0.)])),
        ('touch_segments',('S',[(0.,0.,0.),(1.,0.,0.)]),('S',[(1.,0.,0.),(2.,1.,0.)])),
        ('collapsed_segment',('S',[(1.,1.,3.)]*2),tri),
        ('parallel_segment_face',('S',[(-1.,1.,2.),(5.,1.,2.)]),tri),
        ('coplanar_segment_face',('S',[(-1.,1.,0.),(5.,1.,0.)]),tri),
        ('pierce_face',('S',[(1.,1.,-2.),(1.,1.,2.)]),tri),
        ('touch_face_vertex',('S',[(-2.,-2.,1.),(0.,0.,0.)]),tri),
        ('triangle_contained',('T',[(.5,.5,0.),(1.,.5,0.),(.5,1.,0.)]),tri),
        ('triangle_parallel',('T',[(.5,.5,2.),(1.,.5,2.),(.5,1.,2.)]),tri),
        ('triangle_vertex_touch',('T',[(4.,0.,0.),(6.,0.,0.),(4.,2.,0.)]),tri),
    ]
    rows=[]
    for scale in [.125,1.,32.]:
        for frame in range(3):
            def transform(p):
                p=p if frame==0 else (p[2],-p[0],p[1]) if frame==1 else (2*p[0]+p[2],p[0]+2*p[1]-2*p[2],p[1]+4*p[2])
                return [scale*(x+c) for x,c in zip(p,(2.,-3.,5.))]
            for name,a,b in base:
                rows.append(encode(f'{name}_{scale}_{frame}',(a[0],list(map(transform,a[1]))),(b[0],list(map(transform,b[1])))))
    # Original BRepExtrema_DistShapeShape_Test.cxx BUC60870 input.
    rows.append(encode('BUC60870',('S',[(0.,0.,0.),(0.,1.,0.)]),('P',[(0.,.3,1.)])))
    return '\n'.join(rows)+'\n'


def capture(prefix,output):
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'BRepExtrema_DistShapeShape.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKTopAlgo.*'))),None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True,exist_ok=True)
    executable=output/'occt_proximity_oracle'
    libraries=['TKTopAlgo','TKBRep','TKGeomAlgo','TKGeomBase','TKG3d','TKG2d','TKMath','TKernel']
    run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_proximity_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable)]+['-l'+x for x in libraries],cwd=ROOT)
    data=cases()
    # Capture original observations before invoking any Rust implementation.
    native=run([str(executable)],input=data,cwd=ROOT)
    affine=run([str(executable),'--affine'],input=data,cwd=ROOT)
    for name,text in [('inputs.txt',data),('occt.tsv',native.stdout),('occt-affine.tsv',affine.stdout),('occt-affine-diagnostics.txt',affine.stderr),('native-version.txt',native.stderr)]:
        (output/name).write_text(text)
    return data,native


def parse_rust(text):
    results={}
    for row in text.splitlines():
        w=iter(row.split()); name=next(w)
        if name in results or next(w)!='R': raise ValueError('invalid Rust result')
        d=Q(next(w)); bounds=[]
        for _ in range(2):
            status=next(w)
            if status not in ['U','E']: raise ValueError('invalid enclosure')
            bounds.append([status]+([next(w),next(w)] if status=='E' else []))
        points=[]; parameters=[]
        for _ in range(2):
            n=int(next(w))
            if not 0<=n<=2: raise ValueError('invalid parameter count')
            parameters.append([Q(next(w)) for _ in range(n)])
            points.append(tuple(Q(next(w)) for _ in range(3)))
        if list(w): raise ValueError('extra Rust result fields')
        results[name]=(d,bounds,points,parameters)
    return results


def parse_native(text,affine=False):
    results={}
    for row in text.splitlines():
        name,status,*fields=row.split()
        if name in results or status not in (['P','N','E'] if affine else ['P','F','E']):
            raise ValueError('invalid native status')
        if status=='P':
            distance=float(fields[0])
            if not math.isfinite(distance) or distance<0: raise ValueError('invalid native distance')
            count=0 if affine else int(fields[1])
            coordinates=list(map(float,fields[1 if affine else 2:]))
            if len(coordinates)!=6*count or not all(map(math.isfinite,coordinates)) or (not affine and count<1):
                raise ValueError('malformed native witnesses')
            results[name]=(status,distance,coordinates)
        else:
            if fields: raise ValueError('extra native result fields')
            results[name]=(status,None,[])
    return results


def fingerprint(observation):
    status,distance,_=observation
    return {'status':status,'distance_bits':bits(distance) if distance is not None else None}


def reviewed(oracle,api,row,observation,squared_distance,reviews):
    # Called only after independent exact-value AND witness certification.
    for r in reviews:
        if (r['oracle'],r['api'],r['case'])!=(oracle,api,row.split()[0]): continue
        if r['input_sha256']!=hashlib.sha256(row.encode()).hexdigest(): continue
        if r['native']!=fingerprint(observation): continue
        if Q(r['exact_squared_distance'])!=squared_distance: continue
        return r
    return None


def compare(data,native,affine,rust,oracle,reviews):
    inputs={row.split()[0]:row for row in data.splitlines()}
    if len(inputs)!=len(data.splitlines()): raise ValueError('duplicate input')
    observed={'brep':parse_native(native),'affine':parse_native(affine,True)}
    actual=parse_rust(rust)
    if any(result.keys()!=inputs.keys() for result in [actual,*observed.values()]):
        raise ValueError('oracle case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,
            'brep':{'matches':0,'nonresults':0,'witness_pairs':0},
            'affine':{'matches':0,'not_applicable':0},
            'reviewed_differences':[],'failures':[],'largest_matching_error_fraction_of_budget':0.,
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'occt_source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'domain':'25 ordered point/line/segment/plane/triangle minimum distances; one certified Rust witness; all native witnesses retained'}
    for name,row in inputs.items():
        _,shapes,_=parse(row); shapes=list(map(exact,shapes))
        d,bounds,points,parameters=actual[name]
        try:
            assert d==distance2(*shapes),'independent minimum distance'
            assert bounds==[enclosure(d),enclosure(d,True)],'minimal enclosures'
            certify(shapes,points,parameters,d)
        except AssertionError as error:
            report['failures'].append({'case':name,'reason':'Rust certificate: '+str(error)})
            continue
        report['independently_verified']+=1
        # This corpus is well scaled; float conversion is only for comparing
        # the native numerical observation, never for certifying Rust.
        exact_distance=math.sqrt(float(d))
        for api,results in observed.items():
            status,value,coordinates=observation=results[name]
            if api=='affine' and status=='N':
                if all(k in 'PLF' for k,_ in shapes):
                    report['failures'].append({'case':name,'api':api,'reason':'missing affine observation'})
                else: report[api]['not_applicable']+=1
                continue
            if api=='brep' and status=='F': report[api]['nonresults']+=1
            fraction=abs(value-exact_distance)/(1e-10+2e-12*abs(value)) if status=='P' else math.inf
            witness_failure=None
            if api=='brep' and status=='P':
                for offset in range(0,len(coordinates),6):
                    p,q=[coordinates[offset+i:offset+i+3] for i in (0,3)]
                    budget=1e-10+2e-12*max(map(abs,[*p,*q]))
                    for shape,point in zip(shapes,[p,q]):
                        if distance2(exact(('P',[point])),shape)>Q.from_float(budget)**2:
                            witness_failure='native witness outside its operand'
                    if abs(math.dist(p,q)-value)>1e-10+2e-12*abs(value):
                        witness_failure='native witness distance differs from reported minimum'
                    report['brep']['witness_pairs']+=1
            if witness_failure:
                report['failures'].append({'case':name,'api':api,'reason':witness_failure})
            elif fraction<=1.:
                report[api]['matches']+=1
                report['largest_matching_error_fraction_of_budget']=max(report['largest_matching_error_fraction_of_budget'],fraction)
            else:
                review=reviewed(oracle,api,row,observation,d,reviews)
                evidence={'case':name,'api':api,'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                          'native':fingerprint(observation),'exact_squared_distance':str(d)}
                if review:
                    report['reviewed_differences'].append({**evidence,'id':review['id'],'reason':review['reason']})
                else:
                    report['failures'].append({**evidence,'reason':'unreviewed native difference'})
    report['passed']=not report['failures']
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    p.add_argument('--capture-only',action='store_true')
    p.add_argument('--reuse-capture',action='store_true',help='use the existing local native capture')
    p.add_argument('--strict-native',action='store_true',help='fail also on reviewed native nonresults or numerical differences')
    args=p.parse_args()
    output=ROOT/'target/proximity-oracle'
    if not args.reuse_capture: capture(args.occt_root.resolve(),output)
    data=(output/'inputs.txt').read_text()
    if data!=cases(): raise ValueError('saved capture has different inputs')
    if args.capture_only: return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','proximity_oracle'],input=data,cwd=ROOT)
    (output/'rust.tsv').write_text(rust.stdout)
    review_path=ROOT/'rust/fixtures/occt-proximity-divergences.json'
    report=compare(data,(output/'occt.tsv').read_text(),(output/'occt-affine.tsv').read_text(),rust.stdout,
                   (output/'native-version.txt').read_text().strip(),json.loads(review_path.read_text()) if review_path.exists() else [])
    report['passed'] &= not (args.strict_native and report['reviewed_differences'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key not in ['reviewed_differences','failures']},indent=2))
    print(f"{len(report['reviewed_differences'])} reviewed differences; {len(report['failures'])} failures")
    if not report['passed']: raise SystemExit(1)


if __name__=='__main__': main()
