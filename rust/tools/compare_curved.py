#!/usr/bin/env python3
"""Compare quadratic roots and curved intersections with native OCCT.

The corpus covers the well-conditioned common domain and exact tangencies.
Rust's exact degree/discriminant decisions intentionally differ from OCCT's
threshold policies at extremely small coefficients and near-double roots.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import sys
from compare_occt import ROOT, run


def fixtures():
    rows = []
    polynomials = [(0.,0.,0.),(0.,0.,1.),(0.,2.,-1.),(1.,-2.,1.),
                   (1.,0.,-2.),(1.,0.,1.),(1.,-1e8,1.),(1.,2.,-3.),(1.,0.,0.)]
    for i,values in enumerate(polynomials):
        for factor in [1.,-1.,2**-20,2**20]:
            rows.append(f'q quadratic_{i}_{factor} '+' '.join(format(factor*x,'.17g') for x in values))
    for scale in [0.125,1.,16.]:
        for kind in ['s','y','c']:
            for axis in range(1 if kind == 'c' else 3):
                center = [3.*scale,-4.*scale,8.*scale]
                def world(p):
                    return [center[i]+scale*p[(i+axis)%3] for i in range(3)]
                direction = [[0.,0.,1.][(i+axis)%3] for i in range(3)]
                cases = [((-2.,0.,0.),(2.,0.,0.)),((-2.,1.,0.),(2.,1.,0.)),
                         ((-2.,2.,0.),(2.,2.,0.)),((-2.,0.,0.),(2.,0.5,0.))]
                if kind != 'c':
                    cases += [((0.,0.,-2.),(0.,0.,2.)),((1.,0.,-2.),(1.,0.,2.)),
                              ((-2.,0.,0.),(2.,0.5,1.))]
                for i,(p,q) in enumerate(cases):
                    values = center+direction+[scale]+world(p)+world(q)
                    rows.append(f'{kind} {kind}_{scale}_{axis}_{i} '+' '.join(format(x,'.17g') for x in values))
    return '\n'.join(rows)+'\n'


def observations(text):
    rows = {}
    for row in text.splitlines():
        label,status,*values = row.split()
        if label in rows: raise ValueError('duplicate oracle case')
        numbers = list(map(float,values))
        if not all(math.isfinite(x) for x in numbers): raise ValueError(f'{label}: non-finite oracle')
        rows[label] = (status,numbers)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform == 'darwin' else '/usr')))
    parser.add_argument('--native-only',action='store_true',help='Capture native observations before implementing the Rust comparison')
    args = parser.parse_args()
    prefix = args.occt_root.resolve()
    include = next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'IntAna_IntConicQuad.hxx').exists()),None)
    lib = next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKGeomBase.*'))),None)
    if not include or not lib: parser.error('OCCT modeling-data SDK not found')
    output = ROOT/'target/curved-oracle'
    output.mkdir(parents=True,exist_ok=True)
    executable = output/'occt_curved_oracle'
    command = shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_curved_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKGeomBase','-lTKG3d','-lTKG2d','-lTKMath','-lTKernel']
    run(command,cwd=ROOT)
    data = fixtures()
    native = run([str(executable)],input=data,cwd=ROOT)
    (output/'inputs.txt').write_text(data)
    (output/'occt.tsv').write_text(native.stdout)
    (output/'native-version.txt').write_text(native.stderr)
    expected = observations(native.stdout)
    if args.native_only:
        print(f'{native.stderr.strip()}: captured {len(expected)} independent native observations')
        return
    rust = run(['cargo','run','--quiet','--locked','--example','curved_oracle'],input=data,cwd=ROOT)
    (output/'rust.tsv').write_text(rust.stdout)
    actual = observations(rust.stdout)
    if actual.keys() != expected.keys() or len(expected) != 174: raise AssertionError('oracle case sets differ')
    largest = 0.
    for label,(status,values) in expected.items():
        if actual[label][0] != status or len(actual[label][1]) != len(values):
            raise AssertionError(f'{label}: root count/type differs: {actual[label]} != {expected[label]}')
        for a,e in zip(actual[label][1],values):
            budget = 1e-10+2e-12*abs(e)
            largest = max(largest,abs(a-e)/budget)
            if abs(a-e)>budget: raise AssertionError(f'{label}: {a} != {e}, budget {budget}')
    report = {'oracle':native.stderr.strip(),'cases':len(expected),'largest_error_fraction_of_budget':largest,
              'source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
              'domain':'quadratic roots; line/sphere and line/cylinder; coplanar XY line/circle; well-conditioned cases and exact tangencies',
              'comparison_budget':{'absolute':1e-10,'relative':2e-12},
              'deliberate_differences':'Exact polynomial degree/discriminant; no threshold merging; minimal bounds. General near-degenerate OCCT parity is not asserted.'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    import subprocess
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr)
        raise SystemExit(error.returncode) from error
