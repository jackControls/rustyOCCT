#!/usr/bin/env python3
"""Source-pinned GeomInt_IntSS observations beside the S7b.3b reference: a
torus and a cylinder, a cone or another torus off its axis.

The independent reference (torus_curve_reference.py) gives every case of
torus-curve-cases.txt its canonical rows (folds, tangencies, components with
their winding numbers, rings' points); the native probe
(occt_procedural_intersection_oracle.cpp, planes to tori) runs GeomInt_IntSS
on the kernel's stored frames and samples its lines. Every native sample
must lie on both surfaces within 1e-6 of the case's size, and every
component of the reference curve must carry native samples; an isolated
tangency must be matched by a native point or sample at it, an empty
reference by no line. `--capture` records the native observations before
any kernel code for these pairs exists; later runs must reproduce them.
Differences need a fingerprinted review. The kernel's certified curves
(`torus_curve_probe`) must have the reference's components and tangency
types, and contain its folds, tangencies and rings' points.
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
from compare_procedural_intersections import parse_native, platform_record
import analytic_intersection_reference as ana
import generate_torus_curve_fixtures as fixtures
import torus_curve_reference as ref
from identity_reference import frame_axes

SOURCE_FILE = ROOT/'rust/tools/occt_procedural_intersection_oracle.cpp'
REVIEWS = ROOT/'rust/fixtures/occt-torus-curve-divergences.json'
KERNEL_FILE = ROOT/'rust/kernel/src/intersection/torus_curves.rs'
BOUND = 1e-6
# Each sub-step's native observations were captured before its kernel code:
# (directory, case-name prefixes, whether its kernel code exists).
CAPTURES = {
    's7b3b1': (ROOT/'rust/fixtures/occt-torus-curve-preimplementation', ('tc_', 'tk_'),
               lambda: KERNEL_FILE.exists()),
    's7b3b2': (ROOT/'rust/fixtures/occt-torus-pair-preimplementation', ('tt_',),
               lambda: KERNEL_FILE.exists() and 'torus_pair' in KERNEL_FILE.read_text()),
}


def native_input(prefixes=None):
    rows = []
    for name, a, b in fixtures.cases():
        if prefixes and not name.startswith(prefixes):
            continue
        rows.append(f'case {name}')
        for kind, values in (a, b):
            o, x, _, n = frame_axes(tuple(values[:9]))
            rows.append(' '.join(['surface', kind, *map(repr, (*o, *n, *x)),
                                  *(repr(float(v)) for v in values[9:])]))
        rows.append('end')
    return '\n'.join(rows)+'\n'


def surfaces_of(name, a, b):
    _, surfaces = ana.parse(fixtures.encode(name, a, b))
    return surfaces


def surface_distance(s, p):
    """The distance from a point to a torus, cylinder or cone of a case (a
    cone: its nearer nappe)."""
    o, a = s.axes()
    o = [float(x) for x in o]
    a = [float(x) for x in a]
    la = math.sqrt(sum(x*x for x in a))
    a = [x/la for x in a]
    rel = [x-y for x, y in zip(p, o)]
    h = sum(x*y for x, y in zip(rel, a))
    rho = math.sqrt(max(0.0, sum(x*x for x in rel)-h*h))
    if s.kind == 'torus':
        return abs(math.hypot(rho-s.radius, h)-s.minor)
    if s.kind == 'cylinder':
        return abs(rho-s.radius)
    ca, sa = math.cos(s.angle), math.sin(s.angle)
    shift = s.radius*ca+h*sa
    return min(abs(rho*ca-shift), abs(rho*ca+shift))


def differences(surfaces, native, rows):
    """What separates the native result from the reference's."""
    status, lines, points = native
    if status != 'done':
        return ['not_done']
    samples = [p for line in lines for p in line]
    scale = max([1.0]+[abs(x) for s in surfaces for x in s.frame[:3]]+[s.radius for s in surfaces])
    tol = BOUND*scale
    if rows == ['empty']:
        return [] if not samples and not points else ['spurious']
    out = []
    if any(surface_distance(s, p) > tol for s in surfaces for p in samples+points):
        out.append('off_curve')
    comps, isolated = ref.curve_samples(*surfaces)
    for q in isolated:
        if not any(math.dist(p, q) <= 1e-4*scale for p in samples+points):
            out.append('missed_tangent_point')
            break
    covered = set()
    for p in samples:
        best = min(((math.dist(p, q), k) for k, comp in enumerate(comps) for q in comp), default=None)
        if best is not None and best[0] <= 0.1*scale:
            covered.add(best[1])
    if covered != set(range(len(comps))):
        out.append('missed_component')
    return out


def rust_rows():
    """{case: [(kind, [values])]} from the kernel's probe: enclosures as
    (lo, hi) pairs, integers and words as they are."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'torus_curve_probe'],
                   cwd=ROOT, check=True)
    text = (ROOT/'rust/fixtures/torus-curve-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/torus_curve_probe')], input=text,
                          text=True, capture_output=True, timeout=1200, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        out.setdefault(w[0], []).append(w[1:])
    return out


def rust_differences(rust, expected):
    """The kernel against the reference, row by row: kinds, words and counts
    equal, every reference number inside the kernel's enclosure (up to
    1e-25 relative, the printing; angles modulo a turn)."""
    want = [e.split() for e in expected]
    if len(rust) != len(want) or any(g[0] != w[0] for g, w in zip(rust, want)):
        return ['rust_class']
    for got, row in zip(rust, want):
        kind = row[0]
        if kind in ('empty',):
            continue
        if kind in ('component', 'cluster'):
            if got != row:
                return ['rust_components']
            continue
        numbers = [float(x) for x in row[1:] if x not in ('crossing', 'isolated')]
        words = [x for x in row[1:] if x in ('crossing', 'isolated')]
        values = [x for x in got[1:] if x not in ('crossing', 'isolated')]
        if words != [x for x in got[1:] if x in ('crossing', 'isolated')] or len(values) != 2*len(numbers):
            return ['rust_values']
        for k, x in enumerate(numbers):
            lo, hi = float(values[2*k]), float(values[2*k+1])
            slack = 1e-25*abs(x)
            turn = 2*math.pi if k < 2 else 0.0
            if not any(lo-slack <= x+s <= hi+slack for s in {0.0, turn, -turn}):
                return ['rust_outside_reference']
    return []


def capture(executable, env, key, sdk_manifest):
    CAPTURE, prefixes, exists = CAPTURES[key]
    text = native_input(prefixes)
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native torus curve run failed: '+json.dumps(record)[:2000])
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
        'rust_torus_curve_exists': exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    """Every capture reproduced (on another platform, its reviewed record:
    compare_procedural_intersections.platform_record)."""
    for CAPTURE, prefixes, _ in CAPTURES.values():
        metadata = json.loads((CAPTURE/'capture.json').read_text())
        if metadata['source_reference'] != SOURCE or metadata['rust_torus_curve_exists']:
            raise ValueError('torus curve capture was not a clean pre-implementation reference')
        for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                          ('observations_sha256', 'native.txt')]:
            if metadata[key] != digest(CAPTURE/name):
                raise ValueError('torus curve evidence changed: '+name)
        if (CAPTURE/'inputs.txt').read_text() != native_input(prefixes):
            raise ValueError('the native inputs differ from the captured ones')
        was = parse_native(platform_record(CAPTURE, metadata))
        if set(was) != {n for n in observed if n.startswith(prefixes)}:
            raise ValueError('native cases differ from the capture')
        for name, (status, lines, points) in was.items():
            now = observed[name]
            flat = lambda x: [v for p in [q for line in x[1] for q in line]+x[2] for v in p]
            if now[0] != status or [len(l) for l in now[1]] != [len(l) for l in lines] \
                    or len(now[2]) != len(points) \
                    or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now), flat((status, lines, points)))):
                raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/torus-curve-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', choices=sorted(CAPTURES),
                        help='record one sub-step\'s native observations (before its implementation only)')
    parser.add_argument('--native-only', action='store_true',
                        help='compare the native observations with the reference only (no kernel)')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    expected = {}
    for line in (ROOT/'rust/fixtures/torus-curve-expected.tsv').read_text().splitlines()[1:]:
        name, row = line.split('\t')
        expected.setdefault(name, []).append(row)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'torus-curve-oracle', ('TKGeomAlgo',))
    if args.capture:
        capture(executable, env, args.capture, args.sdk_manifest)
        print('captured', args.capture)
        return
    record = run(executable, native_input(), env)
    if record['exit_code'] != 0:
        raise SystemExit('native torus curve run failed: '+json.dumps(record)[:2000])
    (output/'native-observed.txt').write_text(record['stdout'])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'native_samples': 0,
              'rust_within_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = None if args.native_only else rust_rows()
    for name, a, b in fixtures.cases():
        report['cases'] += 1
        surfaces = surfaces_of(name, a, b)
        # The fixture's rows (the generator's --check keeps them the
        # reference's).
        rows = expected[name]
        if rust is not None:
            wrong = rust_differences(rust[name], rows)
            if wrong:
                report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
                continue
            report['rust_within_reference'] += 1
        native = observed[name]
        report['native_samples'] += sum(len(l) for l in native[1])+len(native[2])
        found = differences(surfaces, native, rows)
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps(native)), 'differences': found}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, lines=len(native[1]), points=len(native[2])))
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
