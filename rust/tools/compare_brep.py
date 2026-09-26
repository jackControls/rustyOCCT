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
# M5: BRep_Tool::Tolerance and OCCT's own measured deviations, captured
# before any Rust enclosure existed.
ENCLOSURES = ROOT/'rust/fixtures/occt-enclosure-preimplementation'
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


def decode_tolerances(row, name):
    """{label: (tolerance or None, measured or None)} from a `T` row: vertices
    `vI:tol:gap`, edges `eI:tol`, uses `uF.W.K:deviation`, faces `fI:tol`."""
    words = row.split()
    if words[:2] != [name, 'T']:
        raise ValueError('malformed tolerance row for '+name)
    out = {}
    for word in words[2:]:
        label, *values = word.split(':')
        numbers = [float(v) for v in values]
        if label[0] == 'v' and len(numbers) == 2:
            out[label] = (numbers[0], numbers[1])
        elif label[0] in 'ef' and len(numbers) == 1:
            out[label] = (numbers[0], None)
        elif label[0] == 'u' and len(numbers) == 1:
            out[label] = (None, numbers[0])
        else:
            raise ValueError(f'malformed tolerance entry {word!r} for {name}')
    return out


def case_size(m):
    """The largest coordinate or radius of a model (at least 1)."""
    return max([abs(x) for v in m.vertices for x in v]
               + [getattr(e.curve, 'radius', 0.0) for e in m.edges] + [1.0])


def same_tolerances(captured, current, allowance=1e-15):
    """Tolerances exactly; measured deviations within `allowance` (rounding
    at the case's scale: OCCT evaluates with the platform's libm) or 1e-9
    relative."""
    if captured.keys() != current.keys():
        return False
    for label, (tol, gap) in captured.items():
        now_tol, now_gap = current[label]
        if tol != now_tol:
            return False
        if gap is not None and not (abs(gap-now_gap) <= max(allowance, 1e-9*abs(gap)) or gap != gap and now_gap != now_gap):
            return False
    return True


def enclosure_capture(tolerances):
    """The M5 pre-implementation observations are unchanged and reproduce."""
    metadata = json.loads((ENCLOSURES/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE or metadata['rust_enclosure_implementation_exists']
            or any(line[3:].startswith('rust/kernel') for line in metadata['rust_worktree_uncommitted'])):
        raise ValueError('enclosure capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(ENCLOSURES/name):
            raise ValueError('enclosure evidence changed: '+name)
    same_inputs((ENCLOSURES/'inputs.txt').read_text(), '\n'.join(text for _, text, _ in native_rows())+'\n')
    captured = {}
    for row in (ENCLOSURES/'native.txt').read_text().splitlines():
        name = row.split()[0]
        captured[name] = decode_tolerances(row, name)
    sizes = {m.name: case_size(m) for m, _, _ in native_rows()}
    for name, observed in tolerances.items():
        if name not in captured or not same_tolerances(captured[name], observed, sizes[name]*2.0**-46):
            raise ValueError('native tolerance observations of '+name+' differ from the capture')


def capture_enclosures(executable, env, oracle_source, sdk_manifest):
    """Record the tolerance rows of every native case (run once, before M5)."""
    cases = native_rows()
    rows = []
    for m, text, _ in cases:
        record = run(executable, text, env)
        lines = record['stdout'].splitlines()
        if record['exit_code'] != 0 or len(lines) != 3:
            raise ValueError('native capture failed for '+m.name)
        decode_tolerances(lines[2], m.name)
        rows.append(lines[2])
    ENCLOSURES.mkdir(parents=True, exist_ok=True)
    (ENCLOSURES/'inputs.txt').write_text('\n'.join(text for _, text, _ in cases)+'\n')
    (ENCLOSURES/'native.txt').write_text('\n'.join(rows)+'\n')
    (ENCLOSURES/'oracle.cpp').write_text(oracle_source.read_text())
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    write(ENCLOSURES/'capture.json', {
        'source_reference': SOURCE, 'rust_revision': revision, 'platform': sys.platform,
        'rust_enclosure_implementation_exists': False, 'rust_worktree_uncommitted': status,
        'sdk_manifest_sha256': digest(sdk_manifest), 'input_sha256': digest(ENCLOSURES/'inputs.txt'),
        'probe_source_sha256': digest(ENCLOSURES/'oracle.cpp'),
        'observations_sha256': digest(ENCLOSURES/'native.txt')})


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


def rust_enclosures():
    """{case: (vertex, fin, face)} largest measured enclosures of each valid case."""
    cases = (ROOT/'rust/fixtures/brep-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/brep_validation_probe'), 'enclosures'],
                          input=cases, text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        name, row = line.split('\t')
        if row != '-':
            out[name] = tuple(float(x) for x in row.split())
    return out


def enclosure_differences(m, measured, observed):
    """T6 of IDENTITY_AND_HISTORY.md on a case valid on both sides: the
    kernel's measured vertex and fin enclosures are not below OCCT's own
    measurements of the same gaps (seam structure excluded), up to 2^-46 of
    the case's size for the frames each side rounds differently; and no
    enclosure exceeds the tolerance OCCT stores for these exactly
    representable constructions."""
    structure = reference.structure_only(m)
    seams = {int(label[1:]) for label in structure if label[0] == 'e'}
    allowance = case_size(m)*2.0**-46
    native_vertex = max([gap for label, (_, gap) in observed.items()
                         if label[0] == 'v' and label not in structure] + [0.0])
    native_use = 0.0
    for label, (_, deviation) in observed.items():
        if label[0] != 'u':
            continue
        f, w, k = (int(x) for x in label[1:].split('.'))
        if m.faces[f].loops[w][k].edge not in seams:
            native_use = max(native_use, deviation)
    stored = min(tol for label, (tol, _) in observed.items() if tol is not None)
    vertex, fin, face = measured
    out = []
    if native_vertex > vertex+allowance or native_use > fin+allowance:
        out.append('enclosure_below_native_measurement')
    if max(vertex, fin, face) > stored:
        out.append('enclosure_above_occt_tolerance')
    return out, {'native_vertex': native_vertex, 'native_use': native_use, 'occt_tolerance': stored,
                 'rust': measured, 'allowance': allowance}


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
    parser.add_argument('--capture-enclosures', action='store_true',
                        help='record the M5 tolerance observations (before implementation only)')
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
    enclosures = rust_enclosures()
    cases = native_rows()
    executable, env, loaded, command = build(prefix, output, cases[0][1])
    if args.capture_enclosures:
        capture_enclosures(executable, env, SOURCE_FILE, args.sdk_manifest)
        print('captured', len(cases), 'tolerance rows')
        return
    tolerances = {}
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'rust_cases_independently_certified': len(issues),
              'not_constructible_natively': [m.name for m in models if not reference.representable(m)],
              'native_timeout_seconds': TIMEOUT, 'native_seconds': {},
              'structure_only_statuses': {}, 'counts_verified': 0, 'enclosures_compared': 0,
              'enclosure_observations': {},
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
            if (len(rows) != 3 or not rows[1].startswith(f'{m.name} N ')
                    or not rows[2].startswith(f'{m.name} T ')):
                raise ValueError(f'malformed native output for {m.name}')
            tolerances[m.name] = decode_tolerances(rows[2], m.name)
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
            found, detail = enclosure_differences(m, enclosures[m.name], tolerances[m.name])
            report['enclosures_compared'] += 1
            report['enclosure_observations'][m.name] = detail
            if found:
                report['failures'].append({'case': m.name, 'reason': ' '.join(found), 'detail': detail})
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
    # The M5 observations must reproduce whatever the Rust side does.
    enclosure_capture(tolerances)
    report['tolerance_rows_reproduced'] = len(tolerances)
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, (list, dict)) and k != 'structure_only_statuses' else v
                      for k, v in report.items() if k != 'native_seconds'}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
