#!/usr/bin/env python3
"""Source-pinned BRepCheck_Analyzer observations beside the generic B-rep validator.

Rust issue lists are first certified equal to the independent reference
validator. Each representable case is then built natively from explicit OCCT
rows and analyzed. A case matches when the verdicts agree and every Rust issue
class with a BRepCheck counterpart has a corresponding native status. Native
seams and the vertices only they use are structure the seamless model does not
have: statuses on them are set aside by rule (brep_reference.structure_only).
Every case valid on both sides must also report, through the kernel's count
synthesis, the distinct subshape counts native OCCT gives its seamed encoding.
Contract differences need a fingerprinted review of the status row; timeouts,
crashes and malformed output cannot be reviewed.
"""
import argparse
import functools
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
import brep_reference as reference
import generate_brep_fixtures

ORIGINAL = ROOT/'rust/fixtures/occt-brep-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-brep-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_brep_check_oracle.cpp'
TOOLKITS = ['TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath', 'TKernel']
# Every preflight case finished within 0.02 s; the deadline only bounds hangs.
TIMEOUT = 120


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def review_for(evidence, reviews):
    return next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                 and all(r.get(k) == v for k, v in evidence.items())), None)


@functools.lru_cache(maxsize=None)
def generate():
    """The independent corpus: generated once, as it also runs the reference."""
    return generate_brep_fixtures.generate()


def native_rows():
    """(model, explicit OCCT rows, inexact pcurve uses) for representable cases."""
    models, _ = generate()
    return [(m, *reference.native(m)) for m in models if reference.representable(m)]


def original_capture():
    """The pre-implementation observations and their exact inputs are unchanged."""
    metadata = json.loads((ORIGINAL/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE
            or metadata['rust_brep_validator_implementation_exists']
            or any(line[3:].startswith('rust/kernel') for line in metadata['rust_worktree_uncommitted'])):
        raise ValueError('original capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.json'), ('expected_sha256', 'brep-expected.tsv')]:
        if metadata[key] != digest(ORIGINAL/name):
            raise ValueError('original native evidence changed: '+name)
    same_inputs((ORIGINAL/'inputs.txt').read_text(), '\n'.join(text for _, text, _ in native_rows())+'\n')


# The capture was generated with macOS libm trigonometry. Generation now uses
# correctly rounded trigonometry so that every host produces the same bytes;
# that moved a few frame and vertex components of two regular polygons by at
# most 1.2e-16. Structure must be identical and numbers within this bound.
NUMBER_DRIFT = 2.0**-50


def same_inputs(captured, current):
    old, new = captured.split(), current.split()
    if len(old) != len(new):
        raise ValueError('current native corpus differs structurally from the pre-implementation inputs')
    for a, b in zip(old, new):
        if a == b:
            continue
        try:
            x, y = float(a), float(b)
        except ValueError:
            raise ValueError(f'current native corpus changed token {a!r} to {b!r}') from None
        if not abs(x-y) <= NUMBER_DRIFT:
            raise ValueError(f'current native corpus moved {a} to {b}')


def run(executable, text, env):
    started = time.monotonic()
    try:
        result = subprocess.run([str(executable)], input=text+'\n', text=True, capture_output=True,
                                timeout=TIMEOUT, env=env)
        return {'exit_code': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr,
                'seconds': time.monotonic()-started}
    except subprocess.TimeoutExpired as error:
        decode = lambda v: v.decode(errors='replace') if isinstance(v, bytes) else v or ''
        return {'exit_code': None, 'timeout_seconds': TIMEOUT, 'stdout': decode(error.stdout),
                'stderr': decode(error.stderr), 'seconds': time.monotonic()-started}


def build(prefix, output, first):
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
        loaded_text = subprocess.run([str(executable)], input=first+'\n', text=True, capture_output=True,
                                     timeout=TIMEOUT, check=True,
                                     env=dict(env, DYLD_PRINT_LIBRARIES='1')).stderr
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    for toolkit in ['TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded BRepCheck library: '+toolkit)
    return executable, env, loaded, command


def rust_issues():
    """Rust issue lists, certified equal to the independent reference validator."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'brep_validation_probe'],
                   cwd=ROOT, check=True)
    cases = (ROOT/'rust/fixtures/brep-cases.txt').read_text()
    rust = subprocess.run([str(ROOT/'target/release/examples/brep_validation_probe')], input=cases,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    expected = (ROOT/'rust/fixtures/brep-expected.tsv').read_text().splitlines()[1:]
    if rust.splitlines() != expected:
        raise ValueError('Rust issue lists differ from the independent reference validator')
    rows = (line.split('\t', 1) for line in expected)
    issues = {name: [i for i in issues.split(';') if i] for name, issues in rows}
    counts = subprocess.run([str(ROOT/'target/release/examples/brep_validation_probe'), 'counts'],
                            input=cases, text=True, capture_output=True, timeout=600, check=True).stdout
    counts = dict(line.split('\t') for line in counts.splitlines())
    if {n for n, c in counts.items() if c != '-'} != {n for n, i in issues.items() if not i}:
        raise ValueError('Rust counts are not reported for exactly the valid cases')
    return issues, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/brep-oracle')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    original_capture()
    models, files = generate()
    for name, text in files.items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    issues, counts = rust_issues()
    cases = native_rows()
    executable, env, loaded, command = build(prefix, output, cases[0][1])
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'rust_cases_independently_certified': len(issues),
              'not_constructible_natively': [m.name for m in models if not reference.representable(m)],
              'native_timeout_seconds': TIMEOUT, 'native_seconds': {},
              'structure_only_statuses': {}, 'counts_verified': 0,
              'matches': [], 'reviewed_differences': [], 'failures': []}
    observations = {}
    oracle = None
    for m, text, inexact in cases:
        record = run(executable, text, env)
        observations[m.name] = dict(record, inexact_pcurve_encodings=inexact)
        write(output/'native.json', observations)
        report['native_seconds'][m.name] = round(record['seconds'], 3)
        if record['exit_code'] != 0:
            reason = 'native timeout' if record['exit_code'] is None else 'native process failed'
            report['failures'].append({'case': m.name, 'reason': reason, 'record': record})
            continue
        oracle = oracle or next(iter(record['stderr'].splitlines()), None)
        rows = record['stdout'].splitlines()
        try:
            if len(rows) != 2 or not rows[1].startswith(f'{m.name} N '):
                raise ValueError(f'malformed native output for {m.name}')
            native = reference.decode_native(rows[0].strip(), m.name)
            native_counts = rows[1].split(maxsplit=2)[2]
        except ValueError as error:
            report['failures'].append({'case': m.name, 'reason': str(error)})
            continue
        structure = reference.structure_only(m)
        set_aside = sorted(label for label in native[1] if label in structure)
        if set_aside:
            report['structure_only_statuses'][m.name] = set_aside
        differences = reference.compare_native(issues[m.name], native, inexact, structure)
        # Count synthesis: a valid seamless body reports the counts OCCT
        # gives for its seamed encoding.
        if not issues[m.name] and native[0]:
            report['counts_verified'] += 1
            if counts[m.name] != native_counts:
                differences = sorted(set(differences) | {'synthesized_counts'})
        if not differences:
            report['matches'].append(m.name)
            continue
        # Reviews fingerprint the status row, which the count row follows.
        evidence = {'case': m.name, 'source_reference': SOURCE, 'oracle': oracle,
                    'input_sha256': sha(text), 'native_stdout_sha256': sha(rows[0]+'\n'),
                    'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, rust=issues[m.name], native=record['stdout'].strip()))
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'input_sha256': sha('\n'.join(text for _, text, _ in cases)+'\n'),
                'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                'observations_sha256': digest(output/'native.json'),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()
                      if k != 'native_seconds'}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
