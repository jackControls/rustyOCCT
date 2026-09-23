#!/usr/bin/env python3
"""Source-pinned native observations plus independently certified Rust minima.

The legacy adapter includes TrimmedSquareDistances endpoints. ExtremaPC uses
PerformWithEndpoints; its other variant is also retained in the raw capture.
Original GTest inputs are adapted, not original assertions executed. Native
stationary candidates are checked for complete minimum-witness coverage; other
stationary points are permitted. Reviewed contract/accuracy differences remain
separate from matches, and no process failure can receive a review exemption.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from build_pinned_occt import SOURCE, digest
from compare_occt import ROOT
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
from generate_spline_proximity_fixtures import generate
from spline_proximity_reference import ABSOLUTE, RELATIVE, verify_rust, decode_native, compare_native

SOURCES = {'legacy': 'occt_spline_projection_oracle.cpp', 'extremapc': 'occt_spline_extremapc_oracle.cpp'}


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def review_for(evidence, reviews):
    return next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                 and all(r.get(k) == v for k, v in evidence.items())), None)


def original_capture():
    base = ROOT/'rust/fixtures/occt-spline-proximity-preimplementation'
    for suffix in ['29', 'degree73']:
        metadata = json.loads((base/f'capture-{suffix}.json').read_text())
        if metadata['source_reference'] != SOURCE or metadata['rust_worktree_dirty']:
            raise ValueError('original capture was not from the recorded clean reference')
        for label, name in [('input', 'inputs'), ('cases', 'cases')]:
            extension = 'txt' if label == 'input' else 'json'
            if metadata[label+'_sha256'] != digest(base/f'{name}-{suffix}.{extension}'):
                raise ValueError('original native inputs changed')
        families = metadata['families'] if suffix == '29' else [dict(metadata, family='legacy')]
        for record in families:
            family = record['family']
            if (record['probe_source_sha256'] != digest(ROOT/'rust/tools'/SOURCES[family])
                    or record['observations_sha256'] != digest(base/f'{family}-{suffix}.json')):
                raise ValueError('original native source or observations changed')
    original = (base/'inputs-29.txt').read_text()+(base/'inputs-degree73.txt').read_text()
    if original != (ROOT/'rust/fixtures/spline-proximity-inputs.txt').read_text():
        raise ValueError('current corpus differs from original pre-implementation input bytes')


def observe(executable, line, seconds, env):
    started = time.monotonic()
    try:
        result = subprocess.run([str(executable)], input=line+'\n', text=True,
                                capture_output=True, timeout=seconds, env=env)
        return {'case': line.split()[0], 'exit_code': result.returncode,
                'stdout': result.stdout, 'stderr': result.stderr, 'seconds': time.monotonic()-started}
    except subprocess.TimeoutExpired as error:
        def decode(value):
            return value.decode(errors='replace') if isinstance(value, bytes) else value or ''
        return {'case': line.split()[0], 'exit_code': None, 'timeout_seconds': seconds,
                'stdout': decode(error.stdout), 'stderr': decode(error.stderr), 'seconds': time.monotonic()-started}


def capture(family, data, prefix, output, sdk_manifest):
    output.mkdir(parents=True, exist_ok=True)
    include, lib = prefix/'include/opencascade', prefix/'lib'
    source = ROOT/'rust/tools'/SOURCES[family]
    executable = output/'oracle'
    command = shlex.split(os.environ.get('CXX', 'c++'))+[
        '-std=c++17', '-O2', str(source), '-I'+str(include), '-L'+str(lib),
        '-Wl,-rpath,'+str(lib), '-o', str(executable), '-lTKGeomAlgo', '-lTKGeomBase',
        '-lTKG3d', '-lTKMath', '-lTKernel']
    built = subprocess.run(command, text=True, capture_output=True)
    (output/'native-build.log').write_text(built.stdout+built.stderr)
    built.check_returncode()
    env = os.environ.copy()
    for key in ['LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']:
        env[key] = str(lib)+(os.pathsep+env[key] if env.get(key) else '')
    link = ['otool', '-L', str(executable)] if sys.platform == 'darwin' else ['ldd', str(executable)]
    linked = subprocess.run(link, text=True, capture_output=True, check=True, env=env)
    (output/'linked-libraries.txt').write_text(linked.stdout+linked.stderr)
    loaded_text = linked.stdout+linked.stderr
    if sys.platform == 'darwin':
        traced = subprocess.run([str(executable)], input=data.splitlines()[0]+'\n', text=True,
                                capture_output=True, timeout=20, check=True,
                                env=dict(env, DYLD_PRINT_LIBRARIES='1'))
        loaded_text = traced.stderr
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    for toolkit in ['TKGeomAlgo', 'TKGeomBase']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded extrema library: '+toolkit)
    records = []
    for line in data.splitlines():
        records.append(observe(executable, line, 20, env))
        write(output/'native.json', records)
    versions = {r['stderr'].splitlines()[0] for r in records if r['stderr'].splitlines()}
    if len(versions) != 1 or not next(iter(versions)).startswith('OCCT '):
        raise ValueError('missing or inconsistent native version')
    metadata = {'source_reference': SOURCE, 'oracle': versions.pop(),
                'sdk_manifest_sha256': digest(sdk_manifest), 'input_sha256': sha(data),
                'source_sha256': digest(source), 'probe_sha256': digest(executable),
                'observations_sha256': digest(output/'native.json'), 'loaded_libraries': loaded,
                'build_command': command}
    write(output/'capture.json', metadata)
    return records, metadata


def reuse(family, data, output, sdk_manifest):
    metadata = json.loads((output/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE or metadata['input_sha256'] != sha(data)
            or metadata['sdk_manifest_sha256'] != digest(sdk_manifest)
            or metadata['source_sha256'] != digest(ROOT/'rust/tools'/SOURCES[family])
            or metadata['observations_sha256'] != digest(output/'native.json')
            or not metadata['loaded_libraries']):
        raise ValueError('native capture changed')
    for p in metadata['loaded_libraries']:
        if digest(Path(p['path'])) != p['sha256']:
            raise ValueError('captured native library changed')
    return json.loads((output/'native.json').read_text()), metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/spline-proximity-oracle')
    parser.add_argument('--reuse-capture', action='store_true')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    original_capture()
    generated = generate()
    for name, text in generated.items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    data = generated['spline-proximity-inputs.txt']
    equations = {row.split()[0]: row for row in generated['spline-proximity.tsv'].splitlines() if not row.startswith('#')}
    source = {line.split()[0]: line for line in data.splitlines()}
    cases = json.loads((ROOT/'rust/fixtures/spline-proximity-inputs.json').read_text())['cases']
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'spline_proximity_probe'], cwd=ROOT, check=True)
    rust = {}
    for name, line in source.items():
        rust[name] = observe(ROOT/'target/release/examples/spline_proximity_probe', line, 120, os.environ.copy())
        write(output/'rust.json', rust)
    certified = {}
    for case in cases:
        record = rust[case['name']]
        if record['exit_code'] != 0:
            raise ValueError('Rust process failed: '+case['name'])
        certified[case['name']] = verify_rust(case, record['stdout'])
    reviews_path = ROOT/'rust/fixtures/occt-spline-proximity-divergences.json'
    reviews = json.loads(reviews_path.read_text())['reviews'] if reviews_path.exists() and not args.strict_native else []
    report = {'source_reference': SOURCE, 'rust_cases_independently_certified': len(certified),
              'isolated_minima': sum(len(c['points']) for c in certified.values()),
              'minimum_intervals': sum(len(c['intervals']) for c in certified.values()),
              'native_comparison_budget': {'absolute': ABSOLUTE, 'relative': RELATIVE},
              'matches': [], 'reviewed_differences': [], 'failures': []}
    for family in SOURCES:
        folder = output/family
        records, metadata = (reuse(family, data, folder, args.sdk_manifest) if args.reuse_capture
                              else capture(family, data, prefix, folder, args.sdk_manifest))
        if [r['case'] for r in records] != list(source):
            raise ValueError('native case set changed')
        for case, record in zip(cases, records):
            name = case['name']
            if record['exit_code'] != 0:
                report['failures'].append({'family': family, 'case': name, 'reason': 'native process failed', 'record': record})
                continue
            try:
                label, native = decode_native(family, record['stdout'], *case['range'])
                if label != name:
                    raise ValueError('wrong native identity')
                differences = compare_native(certified[name], native)
            except (ValueError, StopIteration, IndexError) as error:
                report['failures'].append({'family': family, 'case': name, 'reason': str(error)})
                continue
            if not differences:
                report['matches'].append({'family': family, 'case': name})
                continue
            evidence = {'family': family, 'case': name, 'source_reference': SOURCE,
                        'oracle': metadata['oracle'], 'input_sha256': sha(source[name]),
                        'exact_sha256': sha(equations[name]), 'native_stdout_sha256': sha(record['stdout']),
                        'native_stderr_sha256': sha(record['stderr']), 'differences': differences}
            review = review_for(evidence, reviews)
            report['reviewed_differences' if review else 'failures'].append(review or evidence)
        write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
