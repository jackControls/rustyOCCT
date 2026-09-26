#!/usr/bin/env python3
"""Source-pinned BRepPrimAPI observations beside the primitive reference (S3).

`occt_primitive_oracle.cpp` builds every cone of `primitive-cases.txt` with
BRepPrimAPI_MakeCone and reports its BRepCheck verdict, distinct subshape
counts, BRepGProp mass properties and every face, edge and vertex. Each must
equal what `primitive_reference.py` derives from the specification alone:
counts and verdict exactly, numbers within `BOUND` relative to the case's
size, and faces, edges and vertices as unordered sets matched by geometry.
The kernel's cone builder is compared with the same expectations by
`rust/kernel/tests/cones.rs`, and its history with MakeRevol's by
`compare_revolve_history.py`.

`--capture` records the native observations in
`fixtures/occt-primitive-preimplementation` before any kernel cone code
exists; every later run must reproduce them.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from build_pinned_occt import SOURCE, digest
from compare_brep import run, write
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
from compare_occt import ROOT
import generate_primitive_fixtures

CAPTURE = ROOT/'rust/fixtures/occt-primitive-preimplementation'
SOURCE_FILE = ROOT/'rust/tools/occt_primitive_oracle.cpp'
TOOLKITS = ['TKPrim', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath',
            'TKernel']
TIMEOUT = 600
# Relative to the case's size (its largest coordinate, radius or height,
# at least 1) for points, and to the value itself for lengths, areas,
# volumes and inertia (scaled by size^5).
BOUND = 1e-9


def build(prefix, output):
    include, lib = prefix/'include/opencascade', prefix/'lib'
    executable = output/'oracle'
    command = shlex.split(os.environ.get('CXX', 'c++'))+[
        '-std=c++17', '-O2', str(SOURCE_FILE), '-I'+str(include), '-L'+str(lib),
        '-Wl,-rpath,'+str(lib), '-o', str(executable)]+['-l'+name for name in TOOLKITS]
    built = subprocess.run(command, text=True, capture_output=True)
    (output/'native-build.log').write_text(built.stdout+built.stderr)
    built.check_returncode()
    env = os.environ.copy()
    for key in ['LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']:
        env[key] = str(lib)+(os.pathsep+env[key] if env.get(key) else '')
    if sys.platform == 'darwin':
        loaded_text = subprocess.run([str(executable)], input='', text=True, capture_output=True,
                                     timeout=TIMEOUT, check=True,
                                     env=dict(env, DYLD_PRINT_LIBRARIES='1')).stderr
    else:
        loaded_text = subprocess.run(['ldd', str(executable)], text=True, capture_output=True,
                                     check=True, env=env).stdout
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    for toolkit in ['TKPrim', 'TKTopAlgo']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded library: '+toolkit)
    return executable, env, loaded, command


def parse_observations(text):
    """{case: {'verdict', 'counts', 'props', 'faces', 'edges', 'vertices'}}."""
    out, current = {}, None
    for line in text.splitlines():
        w = line.split()
        if w[0] == 'case':
            current = out.setdefault(w[1], {'faces': [], 'edges': [], 'vertices': []})
            if w[2] == 'exception':
                current['verdict'] = 'exception'
                continue
            current['verdict'] = w[2]
            current['counts'] = tuple(int(x) for x in w[3:9])
        elif w[0] == 'props':
            current['props'] = [float(x) for x in w[1:]]
        elif w[0] == 'face':
            current['faces'].append((w[2], float(w[3]), [float(x) for x in w[4:7]]))
        elif w[0] == 'edge':
            current['edges'].append((w[2], w[3], float(w[4]), [float(x) for x in w[5:8]]))
        elif w[0] == 'vertex':
            current['vertices'].append([float(x) for x in w[2:5]])
        elif w[0] != 'end':
            raise ValueError('malformed observation line '+line)
    return out


def parse_expected(text):
    out = {}
    for line in text.splitlines():
        if line.startswith('#'):
            continue
        name, kind, values = line.split('\t')
        case = out.setdefault(name, {'verdict': 'valid', 'faces': [], 'edges': [], 'vertices': []})
        w = values.split()
        if kind == 'counts':
            case['counts'] = tuple(int(x) for x in w)
        elif kind == 'props':
            case['props'] = [float(x) for x in w]
        elif kind == 'face':
            case['faces'].append((w[0], float(w[1]), [float(x) for x in w[2:5]]))
        elif kind == 'edge':
            case['edges'].append((w[0], w[1], float(w[2]), [float(x) for x in w[3:6]]))
        elif kind == 'vertex':
            case['vertices'].append([float(x) for x in w])
    return out


def sizes():
    return {c.name: max([abs(v) for v in (*c.origin, c.r1, c.r2, c.height)]+[1.0])
            for c in generate_primitive_fixtures.corpus()}


def near(a, b, scale):
    return all(abs(x-y) <= BOUND*scale for x, y in zip(a, b))


def differences(observed, expected, size):
    """Sorted difference classes of one case; empty when it matches."""
    out = set()
    if observed.get('verdict') != 'valid':
        return ['verdict']
    if observed['counts'] != expected['counts']:
        out.add('counts')
    p, q = observed['props'], expected['props']
    scales = [max(abs(q[0]), size**3), max(abs(q[1]), size**2)]+[size]*3+[size**5]*6
    if not all(abs(a-b) <= BOUND*s for a, b, s in zip(p, q, scales)):
        out.add('mass_properties')

    def matched(items, wanted, same):
        left = list(items)
        for w in wanted:
            k = next((i for i, x in enumerate(left) if same(x, w)), None)
            if k is None:
                return False
            left.pop(k)
        return not left
    if not matched(observed['faces'], expected['faces'],
                   lambda a, b: a[0] == b[0] and abs(a[1]-b[1]) <= BOUND*max(b[1], size**2)
                   and near(a[2], b[2], size)):
        out.add('faces')
    if not matched(observed['edges'], expected['edges'],
                   lambda a, b: a[:2] == b[:2] and abs(a[2]-b[2]) <= BOUND*max(b[2], size)
                   and near(a[3], b[3], size)):
        out.add('edges')
    if not matched(observed['vertices'], expected['vertices'], lambda a, b: near(a, b, size)):
        out.add('vertices')
    return sorted(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/primitive-oracle')
    parser.add_argument('--capture', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    verify_sdk(args.occt_root.resolve(), args.sdk_manifest)
    files = generate_primitive_fixtures.generate()
    for name, text in files.items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    cases = files['primitive-cases.txt']
    executable, env, loaded, command = build(args.occt_root.resolve(), output)
    record = run(executable, cases.rstrip('\n'), env)
    if record['exit_code'] != 0:
        raise SystemExit('native run failed: '+json.dumps(record)[:2000])
    oracle = next(iter(record['stderr'].splitlines()), None)
    if args.capture:
        status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                                capture_output=True, check=True).stdout.splitlines()
        if any(line[3:].startswith('rust/kernel') for line in status):
            raise SystemExit('capture requires a kernel tree without changes')
        CAPTURE.mkdir(parents=True, exist_ok=True)
        (CAPTURE/'inputs.txt').write_text(cases)
        (CAPTURE/'native.txt').write_text(record['stdout'])
        (CAPTURE/'oracle.cpp').write_text(SOURCE_FILE.read_text())
        revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                                  check=True).stdout.strip()
        write(CAPTURE/'capture.json', {
            'source_reference': SOURCE, 'oracle': oracle, 'rust_revision': revision,
            'platform': sys.platform, 'rust_cone_implementation_exists': False,
            'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(args.sdk_manifest),
            'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
            'observations_sha256': digest(CAPTURE/'native.txt')})
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('captured evidence changed: '+name)
    if metadata['rust_cone_implementation_exists'] or (CAPTURE/'inputs.txt').read_text() != cases:
        raise ValueError('the capture is not of these cases before implementation')
    observed = parse_observations(record['stdout'])
    captured = parse_observations((CAPTURE/'native.txt').read_text())
    expected = parse_expected(files['primitive-expected.tsv'])
    size = sizes()
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': len(expected), 'matches': [],
              'failures': []}
    for name in expected:
        found = differences(observed[name], expected[name], size[name])
        # The capture must reproduce within the same bound.
        if differences(observed[name], dict(captured[name], verdict='valid'), size[name]) != [] \
                and captured[name].get('verdict') == 'valid':
            found = sorted(set(found) | {'capture_not_reproduced'})
        if found:
            report['failures'].append({'case': name, 'differences': found, 'native': observed[name]})
        else:
            report['matches'].append(name)
    write(output/'report.json', report)
    write(output/'capture.json', {'source_reference': SOURCE, 'oracle': oracle, 'loaded_libraries': loaded,
                                  'build_command': command, 'source_sha256': digest(SOURCE_FILE)})
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        print(json.dumps(report['failures'][:3], indent=1)[:4000])
        raise SystemExit(1)


if __name__ == '__main__':
    main()
