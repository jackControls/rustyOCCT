#!/usr/bin/env python3
"""Source-pinned IntTools_EdgeEdge observations beside the S7d.1 reference:
pairs of lines, circles, ellipses and hyperbolas' branches.

The independent reference (curve_curve_reference.py) gives every case of
curve-curve-cases.txt its rows (empty, coincident, or points with both
curves' parameters and contacts); the native probe
(occt_curve_curve_oracle.cpp) runs IntTools_EdgeEdge on edges of the kernel's
stored curves (lines over [-10, 10], circles and ellipses whole, hyperbolas
over [-3, 3]). Every reference point must be matched by a native vertex
within 1e-6 of the case's size, and no native part may lie elsewhere;
coincidence must be a native edge part. `--capture` records the native
observations before any kernel code for these pairs exists; later runs must
reproduce them (on another platform, its reviewed record). Differences need
a fingerprinted review. The kernel's results (`curve_curve_probe`) must have
the reference's rows, every reference number inside the kernel's enclosures.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from build_pinned_occt import SOURCE, digest
from compare_brep import review_for, run, sha, write
from compare_brep_io import TOOLKITS, build
from compare_curve_surface import platform_record
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
from curve_surface_reference import stored_axes
import generate_curve_curve_fixtures as fixtures

SOURCE_FILE = ROOT/'rust/tools/occt_curve_curve_oracle.cpp'
REVIEWS = ROOT/'rust/fixtures/occt-curve-curve-divergences.json'
KERNEL_FILE = ROOT/'rust/kernel/src/intersection/curve_curve.rs'
CAPTURE = ROOT/'rust/fixtures/occt-curve-curve-preimplementation'
BOUND = 1e-6


def curve_row(kind, values):
    if kind == 'line':
        p0, p1 = values[:3], values[3:6]
        d = [b-a for a, b in zip(p0, p1)]
        return ' '.join(['curve line', *map(repr, (*p0, *d))])
    o, x, _, n = stored_axes(tuple(values[:9]))
    return ' '.join(['curve', kind, *map(repr, (*o, *n, *x)), *(repr(float(v)) for v in values[9:])])


def native_input():
    rows = []
    for name, a, b in fixtures.cases():
        rows += [f'case {name}', curve_row(*a), curve_row(*b), 'end']
    return '\n'.join(rows)+'\n'


def parse_native(stdout):
    """{case: (status, [(point, t1, t2)], [(a, b)])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'P':
            current[1].append(([float(x) for x in w[1:4]], float(w[4]), float(w[5])))
        elif w[0] == 'S':
            current[2].append((float(w[1]), float(w[2])))
        else:
            current = (w[1], [], [])
            out[w[0]] = current
    return out


def expected_rows():
    out = {}
    for line in (ROOT/'rust/fixtures/curve-curve-expected.tsv').read_text().splitlines()[1:]:
        name, row = line.split('\t')
        out.setdefault(name, []).append(row.split())
    return out


def case_scale(a, b):
    return max([1.0]+[abs(x) for x in list(a[1][:3])+list(b[1][:3])])


def differences(a, b, native, rows):
    status, points, segments = native
    if status != 'done':
        return ['not_done']
    tol = BOUND*case_scale(a, b)
    if rows[0][0] == 'empty':
        return [] if not points and not segments else ['spurious']
    if rows[0][0] == 'coincident':
        return [] if segments else ['missed_coincidence']
    out, used = [], set()
    for r in rows:
        p = [float(x) for x in r[3:6]]
        hit = [k for k, (q, _, _) in enumerate(points) if math.dist(p, q) <= tol]
        if not hit:
            out.append('missed_tangent_point' if r[6] == 'tangent' else 'missed_point')
        used.update(hit)
    if len(used) != len(points) or segments:
        out.append('spurious')
    return sorted(set(out))


def rust_rows():
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'curve_curve_probe'],
                   cwd=ROOT, check=True)
    text = (ROOT/'rust/fixtures/curve-curve-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/curve_curve_probe')], input=text,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        out.setdefault(w[0], []).append(w[1:])
    return out


def rust_differences(rust, want, turns):
    """Kinds, contacts and counts equal, every reference number inside the
    kernel's enclosure (up to 1e-25 relative; a circle's or an ellipse's
    angle modulo a turn)."""
    if [g[0] for g in rust] != [w[0] for w in want]:
        return ['rust_class']
    for got, row in zip(rust, want):
        if row[0] in ('empty', 'coincident'):
            continue
        if got[-1] != row[-1] or len(got) != 12:
            return ['rust_contact']
        for k, x in enumerate(float(v) for v in row[1:6]):
            lo, hi = float(got[1+2*k]), float(got[2+2*k])
            slack = 1e-25*abs(x)
            turn = 2*math.pi if k < 2 and turns[k] else 0.0
            if not any(lo-slack <= x+s <= hi+slack for s in {0.0, turn, -turn}):
                return ['rust_outside_reference']
    return []


def capture(executable, env, sdk_manifest):
    text = native_input()
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native curve/curve run failed: '+json.dumps(record)[:2000])
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
        'rust_curve_curve_exists': KERNEL_FILE.exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_curve_curve_exists']:
        raise ValueError('curve/curve capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('curve/curve evidence changed: '+name)
    if (CAPTURE/'inputs.txt').read_text() != native_input():
        raise ValueError('the native inputs differ from the captured ones')
    was = parse_native(platform_record(CAPTURE, metadata))
    if set(was) != set(observed):
        raise ValueError('native cases differ from the capture')
    for name, (status, points, segments) in was.items():
        now = observed[name]
        flat = lambda x: [v for p, s, t in x[1] for v in p+[s, t]]+[v for s in x[2] for v in s]
        if now[0] != status or len(now[1]) != len(points) or len(now[2]) != len(segments) \
                or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now), flat((status, points, segments)))):
            raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/curve-curve-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', action='store_true')
    parser.add_argument('--native-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'curve-curve-oracle', ('TKBO',),
                                             ['TKBO']+TOOLKITS)
    if args.capture:
        capture(executable, env, args.sdk_manifest)
        print('captured')
        return
    record = run(executable, native_input(), env)
    if record['exit_code'] != 0:
        raise SystemExit('native curve/curve run failed: '+json.dumps(record)[:2000])
    (output/'native-observed.txt').write_text(record['stdout'])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    expected = expected_rows()
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'native_points': 0,
              'rust_within_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = None if args.native_only else rust_rows()
    for name, a, b in fixtures.cases():
        report['cases'] += 1
        rows = expected[name]
        if rust is not None:
            turns = [c[0] in ('circle', 'ellipse') for c in (a, b)]
            wrong = rust_differences(rust[name], rows, turns)
            if wrong:
                report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
                continue
            report['rust_within_reference'] += 1
        native = observed[name]
        report['native_points'] += len(native[1])
        found = differences(a, b, native, rows)
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps(native)), 'differences': found}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, points=len(native[1]), segments=len(native[2])))
    write(output/'capture.json', {'source_reference': SOURCE, 'oracle': oracle,
                                  'sdk_manifest_sha256': digest(args.sdk_manifest),
                                  'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                                  'loaded_libraries': loaded, 'build_command': command})
    (output/'native-observed.txt').write_text(record['stdout'])
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
