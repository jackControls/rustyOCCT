#!/usr/bin/env python3
"""Native OCCT spline/plane observations, with explicit supported-domain reporting."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sys
import struct
from compare_occt import ROOT, run


def encode(name,kind,degree,plane,poles,weights,knots,mults):
    header=[name,kind,str(degree),str(len(poles)),str(len(knots))]
    values=[format(x,'.17g') for p in plane for x in p]
    values += [format(x,'.17g') for p,w in zip(poles,weights) for x in [*p,w]]
    values += [x for k,m in zip(knots,mults) for x in [format(k,'.17g'),str(m)]]
    return ' '.join(header+values)


def bernstein_from_roots(roots):
    p=[F(1)]
    for r in roots:
        q=[F(0)]*(len(p)+1)
        for i,a in enumerate(p): q[i]-=a*r; q[i+1]+=a
        p=q
    n=len(p)-1
    b=[sum(p[k]*F(math.comb(i,k),math.comb(n,k)) for k in range(i+1)) for i in range(n+1)]
    scale=math.lcm(*(x.denominator for x in b))
    return [float(x*scale) for x in b]


def cases():
    shapes=[]
    for i,roots in enumerate([[F(1,2)],[F(1,2)]*2,[F(1,4),F(3,4)],
                              [F(1,4),F(1,2),F(3,4)], [F(1,4)]*2+[F(3,4)]*2,
                              [F(0),F(1)], [F(-1),F(2)], [F(1,2)]*7]):
        z=bernstein_from_roots(roots);d=len(z)-1
        shapes.append((f'factor_{i}','B',d,[(float(j),float(j%3),v) for j,v in enumerate(z)],[1.]*len(z),[0.,1.],[d+1]*2))
    for d in [2,3,5,8,25]:
        shapes.append((f'rational_{d}','B',d,[(float(i),float(i%3),float(i)-d/2) for i in range(d+1)],[1.+i%3 for i in range(d+1)],[0.,1.],[d+1]*2))
    for name,d,poles,knots,mults in [
        ('contained',3,[(float(i),float(i%2),0.) for i in range(4)],[0.,1.],[4,4]),
        ('disjoint',3,[(float(i),float(i%2),1.) for i in range(4)],[0.,1.],[4,4]),
        ('join',2,[(float(i),float(i%2),z) for i,z in enumerate([-1.,0.,1.,0.,-1.])],[0.,1.,2.],[3,2,3]),
        ('join_tangent',2,[(float(i),float(i%2),z) for i,z in enumerate([1.,1.,0.,1.,1.])],[0.,1.,2.],[3,2,3]),
        ('partial_overlap',1,[(float(i),0.,z) for i,z in enumerate([-1.,0.,0.,1.])],[0.,1.,2.,3.],[2,1,1,2]),
        ('unclamped',2,[(float(i),float(i%2),z) for i,z in enumerate([-1.,1.,-2.,2.,-1.])],[-2.,-1.,0.,1.,2.,3.,4.,5.],[1]*8),
    ]: shapes.append((name,'S',d,poles,[1.]*len(poles),knots,mults))
    for d in [1,2,3]:
        knots=[0.,1.,2.,3.];mults=[1,d,1,1];n=d+2
        poles=[(float(i%2),float(i),float(2*(i%2)-1)) for i in range(n)]
        shapes.append((f'periodic_{d}','P',d,poles,[1.]*n,knots,mults))
    rows=[]
    for scale in [.125,1.,8.]:
        for frame in ['axis','oblique']:
            a=[2.,-3.,5.]
            u,v,n=([1.,0.,0.],[0.,1.,0.],[0.,0.,1.]) if frame=='axis' else ([2.,1.,0.],[0.,2.,1.],[1.,-2.,4.])
            def at(p): return [a[c]+scale*sum(p[i]*[u,v,n][i][c] for i in range(3)) for c in range(3)]
            plane=[at([0,0,0]),at([1,0,0]),at([0,1,0])]
            for name,kind,d,poles,weights,knots,mults in shapes:
                rows.append(encode(f'{name}_{frame}_{scale}',kind,d,plane,list(map(at,poles)),weights,knots,mults))
    return '\n'.join(rows)+'\n'


def observations(text):
    rows={}
    for row in text.splitlines():
        w=row.split();name=w[0]
        if name in rows: raise ValueError('duplicate case')
        if w[1:]==['FAIL']: rows[name]=None;continue
        np,ns=map(int,w[1:3]);values=list(map(float,w[3:]))
        if min(np,ns)<0 or len(values)!=4*np+2*ns or not all(map(math.isfinite,values)): raise ValueError('invalid intersection observation')
        rows[name]={'points':[values[i:i+4] for i in range(0,4*np,4)],'overlaps':[values[i:i+2] for i in range(4*np,len(values),2)]}
    return rows


def certificates(text, hexadecimal=False):
    rows={}
    number=(lambda s: struct.unpack('>d',bytes.fromhex(s))[0]) if hexadecimal else float
    for row in text.splitlines():
        w=row.split();name=w[0];np,ns=map(int,w[1:3])
        if name in rows or min(np,ns)<0 or len(w)!=3+11*np+2*ns: raise ValueError('invalid certificate shape')
        points=[]
        for i in range(np):
            at=3+11*i;bounds=list(map(number,w[at:at+8]));contact=w[at+8];orders=list(map(int,w[at+9:at+11]))
            if contact not in ['C','T','B'] or any(not 0<=x<=25 for x in orders): raise ValueError('invalid contact')
            if not all(map(math.isfinite,bounds)) or any(a>b for a,b in zip(bounds[::2],bounds[1::2])): raise ValueError('invalid bounds')
            points.append({'bounds':bounds,'contact':contact,'orders':orders})
        tail=list(map(number,w[3+11*np:]))
        if not all(map(math.isfinite,tail)) or any(a>=b for a,b in zip(tail[::2],tail[1::2])): raise ValueError('invalid overlap')
        rows[name]={'points':points,'overlaps':[tail[i:i+2] for i in range(0,len(tail),2)]}
    return rows


def native_fingerprint(native):
    if native is None: return None
    return {'points':[[struct.pack('>d',x).hex() for x in p] for p in native['points']],
            'overlaps':[[struct.pack('>d',x).hex() for x in p] for p in native['overlaps']]}


def differences(native, certificate):
    if native is None: return ['native operation failed']
    reasons=[]
    if len(native['points'])!=len(certificate['points']): reasons.append('isolated point count')
    else:
        for p,c in zip(native['points'],certificate['points']):
            for x,a,b in zip(p,c['bounds'][::2],c['bounds'][1::2]):
                # Distance to the exact enclosure. Bounds are adjacent floats;
                # this budget only describes the native observation comparison.
                if x<a-(1e-10+2e-12*abs(a)) or x>b+(1e-10+2e-12*abs(b)):
                    reasons.append('point parameter or position');break
    if len(native['overlaps'])!=len(certificate['overlaps']): reasons.append('overlap interval count')
    elif any(abs(a-b)>1e-10+2e-12*abs(b) for n,c in zip(native['overlaps'],certificate['overlaps']) for a,b in zip(n,c)):
        reasons.append('overlap interval endpoints')
    return sorted(set(reasons))


@lru_cache(maxsize=256)
def exact_certificate(row):
    from generate_spline_plane_fixtures import expected
    tokens=expected(row)
    name=row.split()[0]
    return certificates(name+' '+' '.join(tokens),hexadecimal=True)[name],hashlib.sha256(' '.join(tokens).encode()).hexdigest()


def reviewed(oracle,row,native,rust,reviews):
    # Even a fully pinned native discrepancy never exempts Rust from the
    # independently calculated complete result, including bounds and contacts.
    exact,exact_hash=exact_certificate(row)
    if rust!=exact: return None
    for review in reviews:
        if (review['oracle'],review['case'])!=(oracle,row.split()[0]): continue
        if review['input_sha256']!=hashlib.sha256(row.encode()).hexdigest(): continue
        if review['native']!=native_fingerprint(native): continue
        if review['exact_certificate_sha256']!=exact_hash: continue
        return review
    return None


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--native-only',action='store_true')
    parser.add_argument('--strict-native',action='store_true',help='fail even independently reviewed native differences')
    args=parser.parse_args();prefix=args.occt_root.resolve()
    include=next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'GeomAPI_IntCS.hxx').exists()),None)
    lib=next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKGeomAlgo.*'))),None)
    if not include or not lib: parser.error('OCCT modeling-algorithms SDK not found')
    output=ROOT/'target/spline-plane-oracle';output.mkdir(parents=True,exist_ok=True)
    executable=output/'occt_spline_plane_oracle'
    run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_spline_plane_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKGeomAlgo','-lTKGeomBase','-lTKG3d','-lTKMath','-lTKernel'],cwd=ROOT)
    data=cases();native=run([str(executable)],input=data,cwd=ROOT)
    for name,text in [('inputs.txt',data),('occt.tsv',native.stdout),('native-version.txt',native.stderr)]: (output/name).write_text(text)
    expected=observations(native.stdout)
    if len(expected)!=len(data.splitlines()): raise AssertionError('incomplete native observations')
    if args.native_only:
        print(native.stderr.strip(),len(expected),'native cases;',sum(x is None for x in expected.values()),'native failures')
        return
    run(['cargo','build','--locked','--release','--example','spline_plane_oracle'],cwd=ROOT)
    rust_process=run([str(ROOT/'target/release/examples/spline_plane_oracle')],input=data,cwd=ROOT)
    (output/'rust.tsv').write_text(rust_process.stdout)
    actual=certificates(rust_process.stdout)
    if actual.keys()!=expected.keys(): raise AssertionError('incomplete Rust observations')
    review_path=ROOT/'rust/fixtures/occt-spline-plane-divergences.json'
    reviews=json.loads(review_path.read_text()) if review_path.exists() else []
    report={'oracle':native.stderr.strip(),'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'cases':len(expected),'independently_verified':0,'native_matches':0,'reviewed_differences':[], 'failures':[]}
    for row in data.splitlines():
        name=row.split()[0];exact,exact_hash=exact_certificate(row)
        if actual[name]!=exact:
            report['failures'].append({'case':name,'reason':'Rust differs from the independent exact certificate'});continue
        report['independently_verified']+=1
        diff=differences(expected[name],actual[name])
        if not diff: report['native_matches']+=1;continue
        review=reviewed(report['oracle'],row,expected[name],actual[name],reviews)
        if review:
            report['reviewed_differences'].append({'id':review['id'],'case':name,'reason':review['reason'],'observed_differences':diff})
        else:
            report['failures'].append({'case':name,'reason':'unreviewed native difference','observed_differences':diff,
                'input_sha256':hashlib.sha256(row.encode()).hexdigest(),'native':native_fingerprint(expected[name]),'exact_certificate_sha256':exact_hash})
    report['passed']=not report['failures'] and not (args.strict_native and report['reviewed_differences'])
    (output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f"{report['oracle']}: {report['native_matches']}/{report['cases']} native matches; {len(report['reviewed_differences'])} reviewed differences; {len(report['failures'])} failures")
    if not report['passed']: raise SystemExit(1)


if __name__=='__main__':
    import subprocess
    try: main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr);raise SystemExit(error.returncode) from error
