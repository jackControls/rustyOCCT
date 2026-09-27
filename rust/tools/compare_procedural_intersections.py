#!/usr/bin/env python3
"""Source-pinned GeomInt_IntSS (GeomAPI_IntSS's engine) observations beside the S7b.1 reference.

The independent reference (procedural_intersection_reference.py) classifies
every case of procedural-intersection-cases.txt (two cylinders with crossing
axes, a cylinder and a sphere off its axis) and gives its canonical rows;
the native probe (occt_procedural_intersection_oracle.cpp) runs
GeomInt_IntSS on the kernel's stored frames and samples its lines. Every
native sample must lie within 1e-6 of the case's size from the reference
curve, and every component of the reference curve (a loop, each ring) must
carry native samples; an empty or single-point reference must be matched by
no line or by a point at it. `--capture` records the native observations
before any kernel code for these pairs exists; later runs must reproduce
them. Differences need a fingerprinted review. The kernel's certified curves
(`procedural_intersection_probe`) must have the reference's class and
components, contain its loop ranges, ends and ring and node points, and put
the middle of each loop within 1e-12 of the reference's.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

from build_pinned_occt import SOURCE, digest
from compare_brep import review_for, run, sha, write
from compare_brep_io import build
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
import analytic_intersection_reference as ana
import generate_procedural_intersection_fixtures as fixtures
import procedural_intersection_reference as ref
from identity_reference import frame_axes

SOURCE_FILE = ROOT/'rust/tools/occt_procedural_intersection_oracle.cpp'
CAPTURE = ROOT/'rust/fixtures/occt-procedural-intersection-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-procedural-intersection-divergences.json'
KERNEL_FILE = ROOT/'rust/kernel/src/intersection/procedural.rs'
BOUND = 1e-6


def native_input():
    rows = []
    for name, a, b in fixtures.cases():
        rows.append(f'case {name}')
        for kind, values in (a, b):
            o, x, _, n = frame_axes(tuple(values[:9]))
            rows.append(' '.join(['surface', kind, *map(repr, (*o, *n, *x)), repr(float(values[9]))]))
        rows.append('end')
    return '\n'.join(rows)+'\n'


def parse_native(stdout):
    """{case: (status, [[line points]], [points])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'L':
            current[1].append([])
        elif w[0] == 'p':
            current[1][-1].append([float(x) for x in w[1:4]])
        elif w[0] == 'P':
            current[2].append([float(x) for x in w[1:4]])
        else:
            current = (w[1], [], [])
            out[w[0]] = current
    return out


def surfaces_of(name, a, b):
    _, surfaces = ana.parse(fixtures.encode(name, a, b))
    return surfaces


def components(cls, cv, rows):
    """The reference curve's components as point functions of (u, branch)."""
    if cls == 'loop':
        return [('loop', None)]
    if cls == 'rings':
        return [('ring', 1), ('ring', -1)]
    return [('figure_eight', None)]


def distance(cv, p, branch=None):
    """Distance from a point to the reference curve (or one branch): the
    point's angle on the ruled cylinder, then the nearer branch there."""
    pm = [mp.mpf(x) for x in p]
    rel = [x-y for x, y in zip(pm, cv.o)]
    u = mp.atan2(sum(a*b for a, b in zip(rel, cv.y)), sum(a*b for a, b in zip(rel, cv.x)))
    best = None
    for sign in ([branch] if branch else [1, -1]):
        q = cv.point(u, sign)
        d = mp.sqrt(sum((a-b)**2 for a, b in zip(pm, q)))
        best = d if best is None else min(best, d)
    return float(best)


def differences(name, surfaces, native):
    """What separates the native result from the reference's."""
    found = ref.curve(*surfaces)
    cls, cv = found
    rows = ref.rows(*surfaces)
    status, lines, points = native
    if status != 'done':
        return ['not_done']
    samples = [p for line in lines for p in line]
    scale = max([1.0]+[abs(x) for s in surfaces for x in s.frame[:3]]+[s.radius for s in surfaces])
    tol = BOUND*scale
    if cls == 'empty':
        return [] if not samples and not points else ['spurious']
    if cls == 'point':
        target = [float(x) for x in rows[0][1]]
        near = [p for p in samples+points if max(abs(a-b) for a, b in zip(p, target)) <= 1e-4*scale]
        if not near:
            return ['missed_tangent_point']
        return [] if len(near) == len(samples+points) else ['spurious']
    out = []
    if any(distance(cv, p) > tol for p in samples+points):
        out.append('off_curve')
    comps = components(cls, cv, rows)
    covered = set()
    for p in samples:
        if cls == 'rings':
            d = {k: distance(cv, p, branch) for k, (_, branch) in enumerate(comps)}
            k = min(d, key=d.get)
            if d[k] <= tol:
                covered.add(k)
        elif distance(cv, p) <= tol:
            covered.add(0)
    if covered != set(range(len(comps))):
        out.append('missed_component')
    return out


def rust_rows():
    """{case: (kind, [(lo, hi)])} from the kernel's probe."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example',
                    'procedural_intersection_probe'], cwd=ROOT, check=True)
    text = (ROOT/'rust/fixtures/procedural-intersection-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/procedural_intersection_probe')], input=text,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        v = [float(x) for x in w[2:]] if w[1] not in ('error',) else []
        out[w[0]] = (w[1], list(zip(v[::2], v[1::2])))
    return out


def rust_differences(rust, expected):
    """The kernel against the reference: its class, every exact reference
    number inside the kernel's enclosure (up to 1e-25 relative, the
    printing), a loop's middle within 1e-12."""
    kind, got = rust
    want = expected[0] if isinstance(expected[0], str) else expected[0][0]
    if kind != want:
        return ['rust_class']
    if kind == 'empty':
        return []
    # The reference's numbers as it prints them (exact zeros as zero).
    numbers = []
    for part in expected[0][1:]:
        numbers += [float(ref.number(x)) for x in (part if isinstance(part, list) else [part])]
    if len(numbers) != len(got):
        return ['rust_values']
    inside = lambda k: got[k][0]-1e-25*abs(numbers[k]) <= numbers[k] <= got[k][1]+1e-25*abs(numbers[k])
    exact = range(len(numbers)) if kind != 'loop' else range(8)
    if not all(inside(k) for k in exact):
        return ['rust_outside_reference']
    if kind == 'loop':
        scale = max([1.0]+[abs(x) for x in numbers])
        if any(abs((got[k][0]+got[k][1])/2-numbers[k]) > 1e-12*scale for k in range(8, 14)):
            return ['rust_loop_middle']
    return []


def capture(executable, env, text, sdk_manifest):
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native procedural intersection run failed: '+json.dumps(record)[:2000])
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
        'rust_procedural_intersection_exists': KERNEL_FILE.exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_procedural_intersection_exists']:
        raise ValueError('procedural intersection capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('procedural intersection evidence changed: '+name)
    if (CAPTURE/'inputs.txt').read_text() != native_input():
        raise ValueError('the native inputs differ from the captured ones')
    was = parse_native((CAPTURE/'native.txt').read_text())
    if set(was) != set(observed):
        raise ValueError('native cases differ from the capture')
    for name, (status, lines, points) in was.items():
        now = observed[name]
        flat = lambda x: [v for p in [q for line in x[1] for q in line]+x[2] for v in p]
        if now[0] != status or [len(l) for l in now[1]] != [len(l) for l in lines] or len(now[2]) != len(points) \
                or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now), flat((status, lines, points)))):
            raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/procedural-intersection-oracle')
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
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'procedural-oracle')
    text = native_input()
    if args.capture:
        capture(executable, env, text, args.sdk_manifest)
        print('captured', len(fixtures.cases()), 'procedural intersection cases')
        return
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native procedural intersection run failed: '+json.dumps(record)[:2000])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'native_samples': 0,
              'rust_within_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = rust_rows()
    for name, a, b in fixtures.cases():
        report['cases'] += 1
        wrong = rust_differences(rust[name], ref.rows(*surfaces_of(name, a, b)))
        if wrong:
            report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
            continue
        report['rust_within_reference'] += 1
        native = observed[name]
        report['native_samples'] += sum(len(l) for l in native[1])+len(native[2])
        found = differences(name, surfaces_of(name, a, b), native)
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
