#!/usr/bin/env python3
"""Native COMMON/SECTION observations for complete linear intersections."""
import argparse
import json
import hashlib
import math
from fractions import Fraction as Q
import os
from pathlib import Path
import shlex
import sys
from compare_occt import ROOT, run
from compare_proximity import encode
from proximity_reference import exact, parse, bits, norm2, sub, dot, cross, on_line, point_triangle
from linear_sets_reference import intersection, decode, encode as encode_exact


def cases():
    tri=('T',[(0.,0.,0.),(4.,0.,0.),(0.,4.,0.)])
    left={'P':[(1.,1.,0.)], 'L':[(-1.,1.,0.),(5.,1.,0.)],
          'S':[(-1.,1.,0.),(5.,1.,0.)], 'F':tri[1], 'T':tri[1]}
    right={'P':[(1.,1.,0.)], 'L':[(1.,-1.,0.),(1.,5.,0.)],
           'S':[(1.,-1.,0.),(1.,5.,0.)],
           'F':[(1.,0.,0.),(1.,3.,0.),(1.,0.,3.)],
           'T':[(1.,0.,-2.),(1.,3.,-2.),(1.,1.,2.)]}
    base=[(a+b,(a,left[a]),(b,right[b])) for a in left for b in right]
    base += [
      ('parallel_lines',('L',[(0.,0.,0.),(1.,0.,0.)]),('L',[(0.,1.,0.),(1.,1.,0.)])),
      ('skew_lines',('L',[(0.,0.,0.),(1.,0.,0.)]),('L',[(0.,0.,1.),(0.,1.,1.)])),
      ('coincident_lines',('L',[(0.,0.,0.),(1.,0.,0.)]),('L',[(4.,0.,0.),(-2.,0.,0.)])),
      ('parallel_planes',('F',tri[1]),('F',[(0.,0.,1.),(1.,0.,1.),(0.,1.,1.)])),
      ('coincident_planes',('F',tri[1]),('F',[(1.,1.,0.),(2.,1.,0.),(1.,2.,0.)])),
      ('contained_line',('L',[(-1.,1.,0.),(1.,1.,0.)]),('F',tri[1])),
      ('parallel_line_plane',('L',[(-1.,1.,1.),(1.,1.,1.)]),('F',tri[1])),
      ('overlap_segments',('S',[(0.,0.,0.),(4.,0.,0.)]),('S',[(3.,0.,0.),(1.,0.,0.)])),
      ('touch_segments',('S',[(0.,0.,0.),(1.,0.,0.)]),('S',[(1.,0.,0.),(2.,1.,0.)])),
      ('disjoint_collinear_segments',('S',[(0.,0.,0.),(1.,0.,0.)]),('S',[(2.,0.,0.),(3.,0.,0.)])),
      ('collapsed_inside',('S',[(1.,1.,0.)]*2),tri),
      ('collapsed_outside',('S',[(1.,1.,1.)]*2),tri),
      ('triangle_contained',('T',[(.5,.5,0.),(1.,.5,0.),(.5,1.,0.)]),tri),
      ('triangle_identical',tri,tri),
      ('triangle_parallel',('T',[(0.,0.,1.),(4.,0.,1.),(0.,4.,1.)]),tri),
      ('triangle_vertex_touch',('T',[(4.,0.,0.),(6.,0.,0.),(4.,2.,0.)]),tri),
      ('triangle_edge_touch',('T',[(0.,0.,0.),(4.,0.,0.),(0.,-4.,0.)]),tri),
      ('triangle_hexagon',('T',[(0.,0.,0.),(6.,0.,0.),(3.,6.,0.)]),('T',[(0.,4.,0.),(6.,4.,0.),(3.,-2.,0.)])),
      ('triangle_quadrilateral',('T',[(0.,0.,0.),(6.,0.,0.),(0.,6.,0.)]),('T',[(1.,-1.,0.),(7.,3.,0.),(1.,5.,0.)])),
      ('vertex_pierce',('S',[(0.,0.,-1.),(0.,0.,1.)]),tri),
      ('triangle_edge_line',('L',[(-2.,0.,0.),(6.,0.,0.)]),tri),
      ('triangle_miss_line',('L',[(0.,5.,0.),(1.,5.,0.)]),tri),
      ('triangle_vertex_plane',('F',[(0.,0.,0.),(0.,0.,1.),(1.,-1.,0.)]),tri),
    ]
    rows=[]
    for scale in [.125,1.,32.]:
      for frame in range(3):
        def transform(p):
          p=p if frame==0 else (p[2],-p[0],p[1]) if frame==1 else (2*p[0]+p[2],p[0]+2*p[1]-2*p[2],p[1]+4*p[2])
          return [scale*(x+c) for x,c in zip(p,(2.,-3.,5.))]
        for name,a,b in base:
          rows.append(encode(f'{name}_{scale}_{frame}',(a[0],list(map(transform,a[1]))),(b[0],list(map(transform,b[1])))))
    # First plane/plane intersection from tests/lowalgos/intss/buc60815.
    # The later swept-surface commands are outside this API and are not replayed.
    rows.append(encode('BUC60815_initial_planes',('F',[(1.,0.,0.),(2.,0.,0.),(1.,1.,0.)]),('F',[(0.,0.,0.),(0.,1.,0.),(0.,0.,1.)])))
    return '\n'.join(rows)+'\n'


def capture(prefix,output):
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'BRepAlgoAPI_Common.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKBO.*'))),None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True,exist_ok=True)
    executable=output/'occt_linear_sets_oracle'
    libraries=['TKBO','TKBool','TKTopAlgo','TKBRep','TKGeomAlgo','TKGeomBase','TKG3d','TKG2d','TKMath','TKernel']
    run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_linear_sets_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable)]+['-l'+x for x in libraries],cwd=ROOT)
    data=cases()
    native=run([str(executable)],input=data,cwd=ROOT)
    for name,text in [('inputs.txt',data),('occt.tsv',native.stdout),('native-version.txt',native.stderr)]:
        (output/name).write_text(text)
    return data,native


def parse_native(text):
    results={}
    for row in text.splitlines():
        w=iter(row.split()); name=next(w); op=next(w); status=next(w)
        key=(name,op)
        if key in results or op not in ['C','S'] or status not in ['R','F','E']: raise ValueError('invalid native row')
        pieces=[]
        if status=='R':
            count=int(next(w))
            if not 0<=count<=32: raise ValueError('invalid native component count')
            for _ in range(count):
                kind=next(w)
                if kind not in ['P','S','G','L','F']: raise ValueError('invalid native component')
                n=int(next(w)) if kind=='G' else 1 if kind=='P' else 2
                if kind=='G' and not 3<=n<=12: raise ValueError('invalid native face vertex count')
                values=[tuple(float(next(w)) for _ in range(3)) for _ in range(n)]
                if not all(math.isfinite(x) for v in values for x in v): raise ValueError('nonfinite native coordinate')
                if kind in ['L','F'] and norm2(values[1])==0: raise ValueError('zero native direction/normal')
                if kind=='S' and values[0]==values[1]: raise ValueError('collapsed native edge')
                pieces.append((kind,values))
        if list(w): raise ValueError('extra native fields')
        results[key]=(status,pieces)
    return results


def parse_rust(text):
    result={}
    for row in text.splitlines():
        w=iter(row.split()); name=next(w)
        if name in result or next(w)!='R': raise ValueError('invalid Rust row')
        result[name]=decode(w)
        if list(w): raise ValueError('extra Rust fields')
    return result


def distance2(point,result):
    kind,p=result
    if kind=='E': return None
    if kind=='P': return norm2(sub(point,p[0]))
    if kind=='S': return norm2(sub(point,on_line(point,p,True)))
    if kind=='G': return min(point_triangle(point,[p[0],p[i],p[i+1]]) for i in range(1,len(p)-1))
    if kind=='L': return norm2(cross(sub(point,p[0]),p[1]))/norm2(p[1])
    if kind=='F': return (dot(point,p[0])-p[1][0])**2/norm2(p[0])
    raise ValueError('kind')


def subset(piece,expected,budget):
    kind,values=piece
    p=[tuple(Q(x) for x in v) for v in values]
    ek,ep=expected
    eps=Q(budget)**2
    def near(point):
        d=distance2(point,expected)
        return d is not None and d<=eps
    if kind in ['P','S','G']:
        if kind=='G' and ek not in ['G','F']: return False
        return all(near(x) for x in p)
    if kind=='L':
        if ek=='L':
            return near(p[0]) and norm2(cross(p[1],ep[1]))<=Q(2e-12)**2*norm2(p[1])*norm2(ep[1])
        return ek=='F' and near(p[0]) and dot(p[1],ep[0])**2<=Q(2e-12)**2*norm2(p[1])*norm2(ep[0])
    return ek=='F' and near(p[0]) and norm2(cross(p[1],ep[0]))<=Q(2e-12)**2*norm2(p[1])*norm2(ep[0])


def complete(pieces,expected,budget):
    """Verify coverage without turning disconnected points into a filled shape."""
    kind,values=expected
    if kind=='E': return not pieces
    if kind in ['P','L','F']:
        return any(p[0]==kind and subset(p,expected,budget) for p in pieces)
    if kind=='G':
        # Native convex face must have all and only the exact extreme vertices.
        # Requiring one complete face avoids treating just its boundary as area.
        return any(k=='G' and len(p)==len(values)
                   and all(any(norm2(sub(v,tuple(Q(x) for x in q)))<=Q(budget)**2 for q in p) for v in values)
                   and all(any(norm2(sub(v,tuple(Q(x) for x in q)))<=Q(budget)**2 for v in values) for q in p)
                   for k,p in pieces)
    d=sub(values[1],values[0]); length2=norm2(d)
    intervals=[]
    for k,p in pieces:
        if k=='S':
            t=sorted(dot(sub(tuple(Q(x) for x in q),values[0]),d)/length2 for q in p)
            intervals.append(t)
    intervals.sort()
    # Parameter budget is derived from the spatial budget, without a sqrt seed.
    tolerance=Q(budget)/max(abs(x) for x in d)
    end=Q(0)
    for lo,hi in intervals:
        if lo>end+tolerance: return False
        end=max(end,hi)
    return bool(intervals) and intervals[0][0]<=tolerance and end>=1-tolerance


def fingerprint(observations):
    return {op:{'status':status,'pieces':[[kind,[[bits(x) for x in v] for v in values]] for kind,values in pieces]}
            for op,(status,pieces) in observations.items()}


def compare(data,native,rust,oracle,reviews=()):
    rows=data.splitlines()
    inputs={row.split()[0]:row for row in rows}
    if len(rows)!=len(inputs) or not inputs: raise ValueError('invalid case set')
    actual=parse_rust(rust); observed=parse_native(native)
    if actual.keys()!=inputs.keys() or set(observed)!={(name,op) for name in inputs for op in ['C','S']}:
        raise ValueError('oracle case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,
            'common_matches':0,'common_lower_dimension_omissions':0,
            'section_subset_matches':0,'combined_complete_matches':0,
            'native_components':0,'reviewed_differences':[],'failures':[],
            'comparison_budget':{'absolute':1e-10,'relative_to_case_scale':2e-12,'angular':2e-12},
            'occt_source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871'}
    dims={'E':-1,'P':0,'S':1,'L':1,'G':2,'F':2,'T':2}
    for name,row in inputs.items():
        _,shapes,rest=parse(row)
        if rest: raise ValueError('extra input')
        expected=intersection(*map(exact,shapes))
        if actual[name]!=expected:
            report['failures'].append({'case':name,'reason':'Rust differs from independently computed complete exact set'})
            continue
        report['independently_verified']+=1
        scale=max([1.,*[abs(x) for _,p in shapes for v in p for x in v],*[abs(float(x)) for v in expected[1] for x in v]])
        budget=1e-10+2e-12*scale
        native_case={op:observed[name,op] for op in ['C','S']}
        common,section=native_case['C'][1],native_case['S'][1]
        report['native_components']+=len(common)+len(section)
        min_dimension=min(0 if k=='S' and p[0]==p[1] else dims[k] for k,p in shapes)
        common_expected=('E',[]) if dims[expected[0]]<min_dimension else expected
        problems=[]
        if any(status!='R' for status,_ in native_case.values()): problems.append('native operation did not return a result')
        common_ok=native_case['C'][0]=='R' and all(subset(p,common_expected,budget) for p in common) and complete(common,common_expected,budget)
        section_ok=native_case['S'][0]=='R' and all(subset(p,expected,budget) for p in section)
        combined_ok=all(status=='R' for status,_ in native_case.values()) and all(subset(p,expected,budget) for p in common+section) and complete(common+section,expected,budget)
        if not common_ok: problems.append('COMMON differs from its dimension-filtered exact set')
        if not section_ok: problems.append('SECTION includes geometry outside the exact intersection')
        if not combined_ok: problems.append('COMMON plus SECTION differs from the complete exact set')
        if common_ok: report['common_matches']+=1
        if common_ok and expected[0]!='E' and common_expected[0]=='E': report['common_lower_dimension_omissions']+=1
        if section_ok: report['section_subset_matches']+=1
        if combined_ok: report['combined_complete_matches']+=1
        if problems:
            review=next((r for r in reviews if r['case']==name and r['oracle']==oracle
                and r['input_sha256']==hashlib.sha256(row.encode()).hexdigest()
                and r['native']==fingerprint(native_case)
                and r['exact_result']==encode_exact(expected)
                and r['differences']==problems),None)
            if review: report['reviewed_differences'].append(review)
            else: report['failures'].append({'case':name,'reason':problems,'native':fingerprint(native_case),'exact_result':encode_exact(expected)})
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    p.add_argument('--output',type=Path,default=ROOT/'target/linear-sets-oracle')
    p.add_argument('--capture-only',action='store_true')
    p.add_argument('--reuse-capture',action='store_true')
    p.add_argument('--strict-native',action='store_true')
    args=p.parse_args()
    if args.reuse_capture:
        data=(args.output/'inputs.txt').read_text()
        if data!=cases(): raise ValueError('captured corpus differs')
        native=(args.output/'occt.tsv').read_text()
        oracle=(args.output/'native-version.txt').read_text().splitlines()[0]
    else:
        data,result=capture(args.occt_root,args.output)
        native=result.stdout; oracle=result.stderr.splitlines()[0]
    if args.capture_only:
        print(json.dumps({'cases':len(data.splitlines()),'observations':len(native.splitlines()),'oracle':oracle,'capture_only':True}))
        return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','linear_sets_oracle'],input=data,cwd=ROOT).stdout
    (args.output/'rust.tsv').write_text(rust)
    path=ROOT/'rust/fixtures/occt-linear-sets-divergences.json'
    reviews=json.loads(path.read_text())['reviews'] if path.exists() else []
    report=compare(data,native,rust,oracle,[] if args.strict_native else reviews)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v if k not in ['failures','reviewed_differences'] else len(v) for k,v in report.items()},indent=2))
    if report['failures']: raise SystemExit(1)

if __name__=='__main__': main()
