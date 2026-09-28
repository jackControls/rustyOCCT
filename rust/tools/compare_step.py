#!/usr/bin/env python3
"""Source-pinned STEPControl_Reader observations beside the STEP reference
(the STEP import track of REVIEW_NOTES.md).

The independent reference (step_reference.py through
generate_step_fixtures.py) gives every body of every fixture in
`fixtures/step` its class, OCCT's counts and its volume, area and centre in
closed form (`step-expected.tsv`); the native probe (occt_step_oracle.cpp)
reads each file with `STEPControl_Reader` at its defaults. The native bodies
must be the reference's (matched by centre), valid, with its counts, their
volumes and areas within 1e-9 relative (BRepGProp's accuracy) and centres
within 1e-9 of the case's size. `--capture` records the native observations
before the kernel's importer exists; later runs must reproduce them on every
platform (counts and verdicts exactly, measures within 1e-9). Differences
need a fingerprinted review. The kernel's bodies (`step_probe`) must be the reference's, valid,
with enclosures containing the reference's measures up to 1e-12 relative
(the file's decimal data are binary64: a surface and its edges agree only
to rounding), and OCCT's up to 1e-9, and their synthesized OCCT counts must
equal the native ones (a difference needs a review).
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from build_pinned_occt import SOURCE, digest
from compare_brep import review_for, run, sha, write
from compare_brep_io import build
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
import generate_step_fixtures as fixtures

SOURCE_FILE = ROOT/'rust/tools/occt_step_oracle.cpp'
REVIEWS = ROOT/'rust/fixtures/occt-step-divergences.json'
CAPTURES = {
    'step_a': (ROOT/'rust/fixtures/occt-step-preimplementation', ROOT/'rust/kernel/src/step/import.rs'),
}
TOOLKITS = ['TKDESTEP', 'TKXSBase', 'TKDE', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d',
            'TKG2d', 'TKMath', 'TKernel']
# BRepGProp on analytic faces (occt-step-preimplementation/NOTES.md).
BOUND = 1e-9
# The kernel against the closed forms: the file's decimal data are binary64,
# so a surface and its edges agree only to rounding.
KERNEL_BOUND = 1e-12


def native_input():
    return ''.join(f'case {name} {len(text.encode())}\n{text}' for name, text, _ in fixtures.cases())


def parse_native(stdout):
    """{case: (status, [(kind, counts, valid, volume|None, area, centre, tolerance)])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'B':
            current[1].append((w[1], [int(x) for x in w[2:8]], w[8] == '1',
                               None if w[9] == '-' else float(w[9]), float(w[10]),
                               [float(x) for x in w[11:14]], float(w[14])))
        else:
            current = (w[1], [])
            out[w[0]] = current
    return out


def expected_rows():
    """{case: [(entity, kind, counts, volume|None, area, centre)]}."""
    out = {}
    for line in fixtures.EXPECTED.read_text().splitlines()[1:]:
        name, entity, kind, counts, volume, area, centre = line.split('\t')
        out.setdefault(name, []).append((int(entity), kind, [int(x) for x in counts.split()],
                                         None if volume == '-' else float(volume), float(area),
                                         [float(x) for x in centre.split()]))
    return out


def case_scale(rows):
    return max([1.0]+[abs(x) for r in rows for x in r[5]]+[r[4]**0.5 for r in rows])


def matched(rows, bodies):
    """Pairs (row, body) by class and nearest centre, or None when the
    classes' numbers differ."""
    free = list(bodies)
    out = []
    for r in rows:
        same = [b for b in free if b[0] == r[1]]
        if not same:
            return None
        b = min(same, key=lambda b: math.dist(b[5], r[5]))
        free.remove(b)
        out.append((r, b))
    return None if free else out


def differences(native, rows):
    status, bodies = native
    if status != 'done':
        return ['not_done']
    pairs = matched(rows, bodies)
    if pairs is None:
        return ['bodies']
    out = []
    scale = case_scale(rows)
    for r, (kind, counts, valid, volume, area, centre, _) in pairs:
        if not valid:
            out.append('invalid')
        if counts != r[2]:
            out.append('counts')
        if (r[3] is not None and abs(volume-r[3]) > BOUND*abs(r[3])) or abs(area-r[4]) > BOUND*abs(r[4]):
            out.append('measure')
        if math.dist(centre, r[5]) > BOUND*scale:
            out.append('centre')
    return sorted(set(out))


def rust_rows():
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'step_probe'],
                   cwd=ROOT, check=True)
    paths = '\n'.join(str(fixtures.OUT/(name+'.stp')) for name, _, _ in fixtures.cases())+'\n'
    rows = subprocess.run([str(ROOT/'target/release/examples/step_probe')], input=paths,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        out.setdefault(w[0], []).append(w[1:])
    return out


def inside(interval, value, bound, size):
    lo, hi = interval
    slack = bound*max(abs(size), 1.0)
    return lo-slack <= value <= hi+slack


def rust_differences(rust, rows, native):
    """Every reference body imported, valid, its enclosures containing the
    reference (a failure otherwise) and OCCT's measures (a difference
    otherwise); its synthesized OCCT counts equal the native ones."""
    bodies = {}
    for w in rust:
        if w[0] == 'error':
            return ['rust_error']
        bodies[int(w[0])] = w[1:]
    out = []
    pairs = matched(rows, native[1]) if native[0] == 'done' else None
    scale = case_scale(rows)
    for r in rows:
        w = bodies.get(r[0])
        if w is None or w[0] != r[1]:
            return ['rust_bodies']
        if w[1] != 'ok':
            return ['rust_rejected']
        counts = [int(x) for x in w[2:8]]
        numbers = [float(x) for x in w[8:]]
        volume, area, centre = (numbers[0:2], numbers[2:4], [numbers[4+2*i:6+2*i] for i in range(3)]) \
            if r[1] == 'solid' else (None, numbers[0:2], [numbers[2+2*i:4+2*i] for i in range(3)])
        if (volume is not None and not inside(volume, r[3], KERNEL_BOUND, r[3])) \
                or not inside(area, r[4], KERNEL_BOUND, r[4]) \
                or not all(inside(centre[i], r[5][i], KERNEL_BOUND, scale) for i in range(3)):
            out.append('rust_outside_reference')
        if pairs is not None:
            body = next(b for q, b in pairs if q is r)
            if counts != body[1]:
                out.append('rust_counts')
            if (volume is not None and not inside(volume, body[3], BOUND, body[3])) \
                    or not inside(area, body[4], BOUND, body[4]) \
                    or not all(inside(centre[i], body[5][i], BOUND, scale) for i in range(3)):
                out.append('rust_outside_native')
    return sorted(set(out))


def capture(executable, env, key, sdk_manifest):
    CAPTURE, kernel_file = CAPTURES[key]
    text = native_input()
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native STEP run failed: '+json.dumps(record)[:2000])
    CAPTURE.mkdir(parents=True, exist_ok=True)
    (CAPTURE/'inputs.txt').write_text(text)
    (CAPTURE/'native.txt').write_text(record['stdout'])
    (CAPTURE/'oracle.cpp').write_text(SOURCE_FILE.read_text())
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    write(CAPTURE/'capture.json', {
        'source_reference': SOURCE, 'oracle': next(iter(record['stderr'].splitlines()), None),
        'platform': sys.platform, 'rust_revision': revision,
        'rust_step_import_exists': kernel_file.exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    for CAPTURE, _ in CAPTURES.values():
        metadata = json.loads((CAPTURE/'capture.json').read_text())
        if metadata['source_reference'] != SOURCE or metadata['rust_step_import_exists']:
            raise ValueError('STEP capture was not a clean pre-implementation reference')
        for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                          ('observations_sha256', 'native.txt')]:
            if metadata[key] != digest(CAPTURE/name):
                raise ValueError('STEP evidence changed: '+name)
        if (CAPTURE/'inputs.txt').read_text() != native_input():
            raise ValueError('the native inputs differ from the captured ones')
        # Reading these files is not near any degeneracy: every platform must
        # reproduce the capture's counts and verdicts exactly and its measures
        # within 1e-9 (VALIDATION.md's allowance table), no platform record.
        was = parse_native((CAPTURE/'native.txt').read_text())
        if set(was) != set(observed):
            raise ValueError('native cases differ from the capture')
        for name, (status, bodies) in was.items():
            now = observed[name]
            flat = lambda bs: [x for b in bs for x in [b[3] or 0.0, b[4], *b[5]]]
            if now[0] != status or len(now[1]) != len(bodies) \
                    or [b[:3] for b in now[1]] != [b[:3] for b in bodies] \
                    or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now[1]), flat(bodies))):
                raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/step-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', choices=sorted(CAPTURES))
    parser.add_argument('--native-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'step-oracle',
                                             ('TKDESTEP', 'TKXSBase'), TOOLKITS)
    if args.capture:
        capture(executable, env, args.capture, args.sdk_manifest)
        print('captured', args.capture)
        return
    record = run(executable, native_input(), env)
    if record['exit_code'] != 0:
        raise SystemExit('native STEP run failed: '+json.dumps(record)[:2000])
    (output/'native-observed.txt').write_text(record['stdout'])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    expected = expected_rows()
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'bodies': 0,
              'rust_within_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = None if args.native_only else rust_rows()
    for name, _, _ in fixtures.cases():
        report['cases'] += 1
        rows = expected[name]
        report['bodies'] += len(rows)
        native = observed[name]
        found = differences(native, rows)
        if rust is not None:
            wrong = rust_differences(rust.get(name, [['error']]), rows, native)
            if any(w in ('rust_error', 'rust_bodies', 'rust_rejected', 'rust_outside_reference') for w in wrong):
                report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust.get(name)})
                continue
            report['rust_within_reference'] += 1
            found = sorted(set(found+wrong))
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps(native)), 'differences': found}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, bodies=native[1], rust=None if rust is None else rust.get(name)))
    write(output/'capture.json', {'source_reference': SOURCE, 'oracle': oracle,
                                  'sdk_manifest_sha256': digest(args.sdk_manifest),
                                  'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                                  'loaded_libraries': loaded, 'build_command': command})
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
