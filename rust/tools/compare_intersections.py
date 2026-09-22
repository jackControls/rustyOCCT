#!/usr/bin/env python3
"""Live line/plane parity on an explicitly well-conditioned OCCT domain.

Near-parallel decisions differ deliberately: Rust uses exact predicates, OCCT
uses angular tolerance. Extreme inputs are tested by independent exact oracles.
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
    for scale in [0.125,1.,32.]:
        for axis in range(3):
            for frame in ['axis','oblique']:
                a = [8.*scale,-3.*scale,5.*scale]
                u,v,n = ([2.,0.,0.],[0.,2.,0.],[0.,0.,1.]) if frame == 'axis' else ([2.,1.,0.],[0.,2.,1.],[1.,-2.,4.])
                def at(x,y,z):
                    return [a[(i+axis)%3]+scale*(x*u[(i+axis)%3]+y*v[(i+axis)%3]+z*n[(i+axis)%3]) for i in range(3)]
                for case,(p,q) in enumerate([
                    (at(.25,.5,-1.),at(.75,.25,2.)),
                    (at(0.,0.,1.),at(2.,1.,1.)),
                    (at(0.,0.,0.),at(2.,1.,0.)),
                    (at(2.,-3.,-4.),at(-2.,3.,-2.)),
                ]):
                    points = [at(0.,0.,0.),at(1.,0.,0.),at(0.,1.,0.),p,q]
                    rows.append(f'{frame}_{scale}_{axis}_{case} '+' '.join(format(x,'.17g') for p in points for x in p))
    return '\n'.join(rows)+'\n'


def parse(text):
    rows = {}
    for line in text.splitlines():
        label,status,*values = line.split()
        if label in rows or status not in ['P','C','N']: raise ValueError('malformed oracle output')
        if len(values) != (4 if status == 'P' else 0): raise ValueError('malformed oracle point')
        values = list(map(float,values))
        if not all(math.isfinite(x) for x in values): raise ValueError('nonfinite oracle output')
        rows[label] = (status,values)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root',type=Path,default=Path(os.environ.get('OCCT_ROOT','/opt/homebrew/opt/opencascade' if sys.platform == 'darwin' else '/usr')))
    args = parser.parse_args()
    prefix = args.occt_root.resolve()
    include = next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include'] if (p/'Standard_Version.hxx').exists()),None)
    lib = next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu'] if any(p.glob('libTKGeomBase.*'))),None)
    if not include or not lib: parser.error('OCCT modeling SDK not found')
    output = ROOT/'target/intersection-oracle'
    output.mkdir(parents=True,exist_ok=True)
    executable = output/'occt_intersection_oracle'
    command = shlex.split(os.environ.get('CXX','c++'))+['-std=c++17','-O2',str(ROOT/'rust/tools/occt_intersection_oracle.cpp'),'-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),'-lTKGeomBase','-lTKG3d','-lTKMath','-lTKernel']
    run(command,cwd=ROOT)
    data = fixtures()
    # Independent native observations are captured first.
    native = run([str(executable)],input=data,cwd=ROOT)
    rust = run(['cargo','run','--quiet','--locked','--example','intersection_oracle'],input=data,cwd=ROOT)
    expected,actual = parse(native.stdout),parse(rust.stdout)
    if expected.keys() != actual.keys() or len(expected) != 72: raise AssertionError('case sets differ')
    largest = 0.
    for label,(status,values) in expected.items():
        if actual[label][0] != status: raise AssertionError(f'{label}: intersection type differs')
        for a,e in zip(actual[label][1],values):
            budget = 1e-10+2e-12*abs(e)
            fraction = abs(a-e)/budget
            largest = max(largest,fraction)
            if fraction > 1.: raise AssertionError(f'{label}: {a} != {e}; budget {budget}')
    report = {'oracle':native.stderr.strip(),'cases':len(expected),'largest_error_fraction_of_budget':largest,
        'occt_source_reference':'3d097a0328e71b826377d4814ab05ec3c3d23871',
        'domain':'well-conditioned line/plane intersection, parallel and contained; affine endpoint parameter',
        'comparison_budget':{'absolute':1e-10,'relative':2e-12},
        'occt_tolerances':{'angular':1e-12,'linear':1e-12},
        'deliberate_difference':'Rust exact parallelism; OCCT angular tolerance. No parity claim for near-parallel or extreme inputs.'}
    for name,text in [('inputs.txt',data),('occt.tsv',native.stdout),('rust.tsv',rust.stdout)]: (output/name).write_text(text)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    import subprocess
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stderr,file=sys.stderr)
        raise SystemExit(error.returncode) from error
