#!/usr/bin/env python3
"""Source-pinned IntTools_EdgeEdge observations beside certified Rust preimages.

Rust rows are first certified against the independent exact generator. Native
queries then use documented structural adapters: a proven finite pole-hull
cover for unbounded lines, and one fundamental-period window per requested
periodic turn with its integer offset. Contract and tolerance differences need
a fingerprinted review; timeouts and other process failures cannot be reviewed.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from build_pinned_occt import SOURCE, digest
from compare_occt import ROOT
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
from generate_spline_linear_fixtures import generate
from spline_linear_reference import (ABSOLUTE, RELATIVE, compare_native, decode_native,
                                     merge_native, native_cases, native_inputs, native_line,
                                     verify_rust)

ORIGINAL = ROOT/'rust/fixtures/occt-spline-linear-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-spline-linear-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_spline_linear_oracle.cpp'
TOOLKITS = ['TKBO', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath', 'TKernel']
# The degree-25 whole-range query takes about 20 seconds locally. A generous
# per-process deadline keeps slower CI runners from turning it into a failure.
TIMEOUT = 120


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def review_for(evidence, reviews):
    return next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                 and all(r.get(k) == v for k, v in evidence.items())), None)


def original_capture():
    """The pre-implementation observations and their exact inputs are unchanged."""
    metadata = json.loads((ORIGINAL/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE or metadata['rust_worktree_dirty']
            or metadata['rust_spline_linear_implementation_exists']):
        raise ValueError('original capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('cases_sha256', 'cases.json'),
                      ('probe_source_sha256', 'oracle.cpp'), ('observations_sha256', 'native.json')]:
        if metadata[key] != digest(ORIGINAL/name):
            raise ValueError('original native evidence changed: '+name)
    if (ORIGINAL/'inputs.txt').read_text() != native_inputs():
        raise ValueError('current native corpus differs from the pre-implementation input bytes')


def queries(case):
    """(row, integer period offset) pairs covering the requested closed range."""
    first, last = case['range']
    if not case['periodic'] or (case['knots'][0] <= first and last <= case['knots'][-1]):
        return [(case, 0)]
    start, period = case['knots'][0], case['knots'][-1]-case['knots'][0]
    result = []
    for turn in range(math.floor((first-start)/period), math.ceil((last-start)/period)):
        lo, hi = max(first, start+turn*period), min(last, start+(turn+1)*period)
        if lo < hi:
            window = dict(case, range=[lo-turn*period, hi-turn*period])
            result.append((window, turn*period))
    return result


def run(executable, rows, env):
    data = ''.join(native_line(case)+'\n' for case, _ in rows)
    started = time.monotonic()
    try:
        result = subprocess.run([str(executable)], input=data, text=True, capture_output=True,
                                timeout=TIMEOUT, env=env)
        return {'input_sha256': sha(data), 'windows': len(rows),
                'offsets': [float(o) for _, o in rows], 'exit_code': result.returncode,
                'stdout': result.stdout, 'stderr': result.stderr, 'seconds': time.monotonic()-started}
    except subprocess.TimeoutExpired as error:
        decode = lambda v: v.decode(errors='replace') if isinstance(v, bytes) else v or ''
        return {'input_sha256': sha(data), 'windows': len(rows),
                'offsets': [float(o) for _, o in rows], 'exit_code': None, 'timeout_seconds': TIMEOUT,
                'stdout': decode(error.stdout), 'stderr': decode(error.stderr),
                'seconds': time.monotonic()-started}


def native_parts(case, record):
    rows = [r for r in record['stdout'].splitlines() if r]
    windows = queries(case)
    if len(rows) != len(windows):
        raise ValueError('native window count changed')
    return [decode_native(row, case, offset) for row, (_, offset) in zip(rows, windows)]


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
    link = ['otool', '-L', str(executable)] if sys.platform == 'darwin' else ['ldd', str(executable)]
    linked = subprocess.run(link, text=True, capture_output=True, check=True, env=env)
    loaded_text = linked.stdout+linked.stderr
    if sys.platform == 'darwin':
        first = native_line(native_cases()[0])+'\n'
        loaded_text = subprocess.run([str(executable)], input=first, text=True, capture_output=True,
                                     timeout=TIMEOUT, check=True,
                                     env=dict(env, DYLD_PRINT_LIBRARIES='1')).stderr
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    for toolkit in ['TKBO', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded edge/edge library: '+toolkit)
    return executable, env, loaded, command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/spline-linear-oracle')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    original_capture()
    for name, text in generate().items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    data = native_inputs()
    cases = native_cases()
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'spline_linear_probe'],
                   cwd=ROOT, check=True)
    rust = subprocess.run([str(ROOT/'target/release/examples/spline_linear_probe')], input=data,
                          text=True, capture_output=True, timeout=600, check=True).stdout.splitlines()
    (output/'rust.txt').write_text('\n'.join(rust)+'\n')
    if len(rust) != len(cases):
        raise ValueError('Rust row count changed')
    certified = {c['name']: verify_rust(c, line) for c, line in zip(cases, rust)}
    executable, env, loaded, command = build(prefix, output)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'rust_cases_independently_certified': len(certified),
              'isolated_points': sum(len(c['points']) for c in certified.values()),
              'overlap_intervals': sum(len(c['intervals']) for c in certified.values()),
              'native_comparison_budget': {'absolute': ABSOLUTE, 'relative': RELATIVE},
              'native_timeout_seconds': TIMEOUT, 'native_seconds': {}, 'matches': [], 'reviewed_differences': [], 'failures': []}
    observations = {}
    oracle = None
    for case in cases:
        name = case['name']
        record = run(executable, queries(case), env)
        observations[name] = record
        write(output/'native.json', observations)
        report['native_seconds'][name] = round(record['seconds'], 3)
        if record['exit_code'] != 0:
            reason = 'native timeout' if record['exit_code'] is None else 'native process failed'
            report['failures'].append({'case': name, 'reason': reason, 'record': record})
            continue
        oracle = oracle or next(iter(record['stderr'].splitlines()), None)
        try:
            differences = compare_native(certified[name], merge_native(native_parts(case, record)))
        except (ValueError, IndexError) as error:
            report['failures'].append({'case': name, 'reason': str(error)})
            continue
        if not differences:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'input_sha256': sha(native_line(case)),
                    'native_stdout_sha256': sha(record['stdout']),
                    'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(review or evidence)
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'input_sha256': sha(data), 'source_sha256': digest(SOURCE_FILE),
                'probe_sha256': digest(executable), 'observations_sha256': digest(output/'native.json'),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
