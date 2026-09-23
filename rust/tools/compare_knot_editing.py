#!/usr/bin/env python3
"""Capture OCCT knot editing before Rust execution; compare full representations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from compare_occt import ROOT, run
from compare_bezier_editing import shapes as existing_shapes


def shapes():
    yield from existing_shapes()
    # Unclamped ends, including controls inactive on the fundamental domain.
    for d in range(1, 26):
        knots = [float(i) for i in range(-2, 5)]
        mults = [d, d, 1, d, 1, d, d]
        n = sum(mults)-d-1
        poles = [(float(i-2), float(i*i%13-6), float(i*7%11-5)) for i in range(n)]
        yield f'unclamped_nonuniform_d{d}', 'S', d, poles, [float(1+i%3)/2 for i in range(n)], knots, mults
        # Enough poles remain after deleting the old seam. All knots are
        # algebraically redundant, making origin-shift removal an exact case.
        knots = [float(i) for i in range(d+4)]
        yield f'constant_periodic_d{d}', 'P', d, [(2.,-3.,5.)]*(d+3), [1.]*(d+3), knots, [1.]*(d+4)


def encode(shape, suffix, operations):
    name, kind, degree, poles, weights, knots, mults = shape
    words = [name+'_'+suffix, kind, str(degree), str(len(poles)), str(len(knots)), str(len(operations))]
    words += [format(v, '.17g') for p, w in zip(poles, weights) for v in [*p, w]]
    words += [s for k, m in zip(knots, mults) for s in [format(k, '.17g'), str(int(m))]]
    words += [s for op, u, m in operations for s in [op, format(u, '.17g'), str(m)]]
    return ' '.join(words)


def cases():
    rows = []
    for shape in shapes():
        name, kind, degree, poles, _, knots, mults = shape
        flat = [k for k, m in zip(knots, mults) for _ in range(int(m))]
        a, b = (knots[0], knots[-1]) if kind=='P' else (flat[degree], flat[len(poles)])
        first_span = next(k for k in knots if k>a)
        u, v = (3*a+first_span)/4, (a+3*first_span)/4
        rows.append(encode(shape, 'identity', []))
        rows.append(encode(shape, 'insert', [('I',u,degree)]))
        rows.append(encode(shape, 'sequence', [('I',v,1),('I',u,degree),('I',v,degree),('I',u,1)]))
        rows.append(encode(shape, 'roundtrip', [('I',u,1),('R',u,0)]))
        rows.append(encode(shape, 'full_roundtrip', [('I',u,degree),('R',u,0)]))
        interior = next((k for k in knots if a<k<b), None)
        if interior is not None and len(poles)>degree+1:
            m=int(mults[knots.index(interior)])
            rows.append(encode(shape, 'remove_existing', [('R',interior,m-1)]))
        if kind=='P':
            rows.append(encode(shape, 'seam', [('I',a,degree),('I',b,degree)]))
            if degree>mults[0]:
                rows.append(encode(shape, 'seam_roundtrip', [('I',b,degree),('R',a,int(mults[0]))]))
            if name.startswith('constant_'):
                rows.append(encode(shape, 'remove_origin', [('R',a,0)]))
                rows.append(encode(shape, 'remove_end_alias', [('R',b,0)]))
        elif knots[0]<a:
            rows.append(encode(shape, 'domain_ends', [('I',a,degree),('I',b,degree)]))
    # A nonconstant removable seam, reconstructed independently by the complete
    # Cox coefficient system from a uniform quadratic periodic curve. All
    # controls are dyadic, so native inputs preserve the redundant seam exactly.
    redundant=('nonconstant_redundant_seam','P',2,
               [(1.5,.75,0.),(1.75,1.5,.5),(1.,3.,2.),(-1.,1.,1.),(0.,0.,0.)],
               [1.]*5,[.5,1.,2.,3.,4.,4.5],[1]*6)
    rows.append(encode(redundant,'remove_origin',[('R',.5,0)]))
    rows.append(encode(redundant,'remove_end_alias',[('R',4.5,0)]))
    return '\n'.join(rows)+'\n'


def capture(prefix, output):
    include = next((p for p in [prefix/'include/opencascade', prefix/'inc', prefix/'include'] if (p/'Geom_BSplineCurve.hxx').exists()), None)
    lib = next((p for p in [prefix/'lib', prefix/'lib64', prefix/'lib/x86_64-linux-gnu', prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKG3d.*'))), None)
    if not include or not lib: raise ValueError('OCCT modeling SDK not found')
    output.mkdir(parents=True, exist_ok=True)
    executable = output/'occt_knot_editing_oracle'
    command = shlex.split(os.environ.get('CXX', 'c++'))+['-std=c++17', '-O2', str(ROOT/'rust/tools/occt_knot_editing_oracle.cpp'), '-I'+str(include), '-L'+str(lib), '-Wl,-rpath,'+str(lib), '-o', str(executable), '-lTKG3d', '-lTKMath', '-lTKernel']
    build = run(command, cwd=ROOT)
    (output/'native-build.log').write_text(build.stdout+build.stderr)
    data = cases()
    native = run([str(executable)], input=data, cwd=ROOT)
    for file, text in [('inputs.txt',data), ('occt.tsv',native.stdout), ('native-version.txt',native.stderr)]:
        (output/file).write_text(text)
    record = {'head':run(['git','rev-parse','HEAD'],cwd=ROOT).stdout.strip(),
              'implementation_present':(ROOT/'rust/kernel/src/curve/knot_editing.rs').exists(),
              'cases':len(data.splitlines()), 'input_sha256':hashlib.sha256(data.encode()).hexdigest(),
              'native_sha256':hashlib.sha256(native.stdout.encode()).hexdigest(),
              'oracle':native.stderr.splitlines()[0]}
    (output/'capture.json').write_text(json.dumps(record,indent=2)+'\n')
    return data, native.stdout, record['oracle']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform=='darwin' else '/usr')))
    parser.add_argument('--output', type=Path, default=ROOT/'target/knot-editing-oracle')
    parser.add_argument('--capture-only', action='store_true')
    parser.add_argument('--reuse-capture', action='store_true')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    if args.reuse_capture:
        data=(args.output/'inputs.txt').read_text()
        if data!=cases(): raise ValueError('captured corpus changed')
        native=(args.output/'occt.tsv').read_text()
        oracle=(args.output/'native-version.txt').read_text().splitlines()[0]
    else: data, native, oracle = capture(args.occt_root, args.output)
    if args.capture_only:
        print(f'{oracle}: {len(data.splitlines())} complete knot editing observations captured')
        return
    rust=run(['cargo','run','--quiet','--release','--locked','--example','knot_editing_oracle'],input=data,cwd=ROOT).stdout
    (args.output/'rust.tsv').write_text(rust)
    path=ROOT/'rust/fixtures/occt-knot-editing-divergences.json'
    reviews=json.loads(path.read_text())['reviews'] if path.exists() and not args.strict_native else []
    # Ordinary Rust tests use these independently solved coefficient fixtures;
    # CI regenerates all of them before running this comparison.
    fixtures={}
    for row in (ROOT/'rust/fixtures/knot-editing.tsv').read_text().splitlines():
        if row.startswith('#'): continue
        source,result=row.split(' | ')
        if source in fixtures: raise ValueError('duplicate fixture')
        fixtures[source]=result
    exact='\n'.join(fixtures[row] for row in data.splitlines())+'\n'
    report=compare(data,native,rust,exact,oracle,reviews)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:len(v) if k in ['failures','reviewed_differences'] else v for k,v in report.items()},indent=2))
    if report['failures']: raise SystemExit(1)


def fingerprint(observation):
    import struct
    flags,(degree,periodic,knots,mults,controls),domain=observation
    bits=lambda x:struct.pack('>d',x).hex()
    return [list(map(int,flags)),degree,periodic,list(map(bits,domain)),list(map(bits,knots)),list(mults),[[bits(x) for x in c] for c in controls]]


def compare(data,native,rust,exact,oracle,reviews=()):
    from knot_editing_reference import decode
    inputs={row.split()[0]:row for row in data.splitlines()}
    if len(inputs)!=len(data.splitlines()): raise ValueError('duplicate input')
    actual,observed,expected=decode(rust),decode(native,True),decode(exact)
    exact_rows={row.split()[0]:row for row in exact.splitlines()}
    if actual.keys()!=inputs.keys() or observed.keys()!=inputs.keys() or expected.keys()!=inputs.keys():
        raise ValueError('case sets differ')
    report={'oracle':oracle,'cases':len(inputs),'independently_verified':0,'matched_cases':0,
            'largest_error_fraction_of_budget':0.,'reviewed_differences':[],'failures':[],
            'comparison_budget':{'absolute':1e-10,'relative':2e-12},
            'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
            'removal_contract':'Exact positive homogeneous equality; native RemoveKnot uses tolerance 1e-9.'}
    for name,row in inputs.items():
        if actual[name]!=expected[name]:
            report['failures'].append({'case':name,'reason':'Rust differs from independent coefficient equations'})
            continue
        report['independently_verified']+=1
        flags,curve,domain=expected[name]; nflags,ncurve,ndomain=observed[name]
        differences=[]
        if flags!=nflags: differences.append('operation success flags')
        if curve[:2]!=ncurve[:2]: differences.append('degree or periodicity')
        if curve[2:4]!=ncurve[2:4]: differences.append('knots or multiplicities')
        if domain!=ndomain: differences.append('parameter domain')
        if len(curve[4])!=len(ncurve[4]): differences.append('control count')
        if not differences:
            values=[float(x) for c in curve[4] for x in [c[0]/c[3],c[1]/c[3],c[2]/c[3],c[3]]]
            native_values=[x for c in ncurve[4] for x in c]
            errors=[abs(a-b)/(1e-10+2e-12*abs(a)) for a,b in zip(values,native_values)]
            report['largest_error_fraction_of_budget']=max(report['largest_error_fraction_of_budget'],*errors)
            differences.extend(f'control field {i}' for i,e in enumerate(errors) if e>1.)
        if not differences:
            report['matched_cases']+=1
            continue
        evidence={'case':name,'oracle':oracle,'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                  'native':fingerprint(observed[name]),'exact_sha256':hashlib.sha256(exact_rows[name].encode()).hexdigest(),
                  'differences':differences}
        review=next((r for r in reviews if all(r.get(k)==v for k,v in evidence.items())),None)
        if review: report['reviewed_differences'].append(review)
        else: report['failures'].append(evidence)
    return report


if __name__=='__main__':
    try: main()
    except subprocess.CalledProcessError as e:
        print(e.stderr,file=sys.stderr)
        raise SystemExit(e.returncode) from e
