#!/usr/bin/env python3
"""Source-pinned IntAna_QuadQuadGeo observations beside the S7a reference.

The independent reference (analytic_intersection_reference.py) classifies
every case of analytic-intersection-cases.txt and gives its canonical items;
the native probe (occt_analytic_intersection_oracle.cpp) runs
IntAna_QuadQuadGeo on the kernel's stored frames. Native results are put in
the same canonical form (a line's point nearest the origin, unit directions
with their first nonzero coordinate positive, a hyperbola's two branches
one item) and compared with the reference: the same items, numbers within
1e-9 of each item's magnitude. `--capture` records the native observations
before any kernel intersection code exists; later runs must reproduce them.
Differences need a fingerprinted review (IntAna snaps near-degenerate
configurations with its angular and linear tolerances).
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
import analytic_intersection_reference as ref
import generate_analytic_intersection_fixtures as fixtures
from identity_reference import frame_axes

SOURCE_FILE = ROOT/'rust/tools/occt_analytic_intersection_oracle.cpp'
CAPTURE = ROOT/'rust/fixtures/occt-analytic-intersection-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-analytic-intersection-divergences.json'
BOUND = 1e-9


def native_input():
    """The cases with the kernel's stored frame axes (origin, normal, x)."""
    rows = []
    for name, a, b in fixtures.cases():
        rows.append(f'case {name}')
        for kind, values in (a, b):
            o, x, _, n = frame_axes(tuple(values[:9]))
            rest = [repr(float(v)) for v in values[9:]]
            rows.append(' '.join(['surface', kind, *map(repr, (*o, *n, *x)), *rest]))
        rows.append('end')
    return '\n'.join(rows)+'\n'


def canonical_unit(v):
    lead = next((x for x in v if abs(x) > 1e-300), 1.0)
    s = 1.0 if lead > 0 else -1.0
    norm = math.sqrt(sum(x*x for x in v))
    return tuple(s*x/norm for x in v)


def nearest(p, d):
    k = sum(a*b for a, b in zip(p, d))/sum(x*x for x in d)
    return tuple(a-k*b for a, b in zip(p, d))


def parse_native(stdout):
    """{case: (type, [canonical items])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] in 'PLCEHB' and len(w[0]) == 1:
            v = [float(x) for x in w[1:]]
            if w[0] == 'P':
                item = ('point', tuple(v[:3]))
            elif w[0] == 'L':
                item = ('line', nearest(v[:3], v[3:6]), canonical_unit(v[3:6]))
            elif w[0] == 'C':
                item = ('circle', tuple(v[:3]), canonical_unit(v[3:6]), v[6])
            elif w[0] == 'E':
                item = ('ellipse', tuple(v[:3]), canonical_unit(v[3:6]), canonical_unit(v[6:9]), v[9], v[10])
            elif w[0] == 'H':
                item = ('hyperbola', tuple(v[:3]), canonical_unit(v[3:6]), canonical_unit(v[6:9]), v[9], v[10])
            else:
                item = ('parabola', tuple(v[:3]), canonical_unit(v[3:6]), canonical_unit(v[6:9]), v[9])
            if item not in current[1]:
                current[1].append(item)
        else:
            current = (w[1], [])
            out[w[0]] = current
    return out


def flat(item):
    values = []
    for part in item[1:]:
        values += [float(x) for x in part] if isinstance(part, tuple) else [float(part)]
    return values


def reference_items(name, a, b):
    _, surfaces = ref.parse(fixtures.encode(name, a, b))
    return ref.canonical(ref.intersect(*surfaces))


def differences(expected, native_type, native_items):
    """What separates the native result from the reference's."""
    words = [e for e in expected if isinstance(e, str)]
    if words:
        want = {'empty': 'empty', 'same': 'same', 'not_conic': 'no_geometric_solution'}[words[0]]
        return [] if native_type == want else ['classification']
    if native_type in ('empty', 'same', 'no_geometric_solution', 'not_done', 'failure'):
        return ['classification']
    got = ref.canonical(native_items)
    if [g[0] for g in got] != [e[0] for e in expected]:
        return ['items']
    out = []
    for g, e in zip(got, expected):
        ev, gv = flat(e), flat(g)
        scale = max([1.0]+[abs(x) for x in ev])
        if any(abs(x-y) > BOUND*scale for x, y in zip(ev, gv)):
            out.append('parameters')
    return sorted(set(out))


def capture(executable, env, text, sdk_manifest):
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native analytic intersection run failed: '+json.dumps(record)[:2000])
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
        'rust_analytic_intersection_exists': (ROOT/'rust/kernel/src/intersection/analytic.rs').exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    """The pre-implementation observations are unchanged and reproduce."""
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_analytic_intersection_exists']:
        raise ValueError('analytic intersection capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('analytic intersection evidence changed: '+name)
    if (CAPTURE/'inputs.txt').read_text() != native_input():
        raise ValueError('the native inputs differ from the captured ones')
    was = parse_native((CAPTURE/'native.txt').read_text())
    if set(was) != set(observed):
        raise ValueError('native cases differ from the capture')
    for name, (kind, items) in was.items():
        now = observed[name]
        if now[0] != kind or len(now[1]) != len(items) or any(
                any(abs(x-y) > 1e-12*max(1.0, abs(x)) for x, y in zip(flat(a), flat(b)))
                for a, b in zip(ref.canonical(now[1]), ref.canonical(items))):
            raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/analytic-intersection-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', action='store_true',
                        help='record the native observations (before implementation only)')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    for name, text in fixtures.generate().items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'analytic-oracle')
    text = native_input()
    if args.capture:
        capture(executable, env, text, args.sdk_manifest)
        print('captured', len(fixtures.cases()), 'analytic intersection cases')
        return
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native analytic intersection run failed: '+json.dumps(record)[:2000])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'matches': [],
              'reviewed_differences': [], 'failures': []}
    for name, a, b in fixtures.cases():
        report['cases'] += 1
        expected = reference_items(name, a, b)
        kind, items = observed[name]
        found = differences(expected, kind, items)
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps([kind, items])), 'differences': found}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, native=[kind, items], reference=[ref.text(e) for e in expected]))
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
