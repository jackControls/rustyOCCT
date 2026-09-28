#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Splitter observations beside the S8 reference:
prisms of line and arc profiles split by a plane (S8a), cones, frusta,
spheres and zones (S8c: `split-primitive-cases.txt`, built natively by
BRepPrimAPI_MakeCone and MakeSphere on the kernel's stored frame axes), whole
tori (S8d) and cones and zones cut in conics (S8d.2:
`split-conic-cases.txt`, whose capture came after the kernel code and is
recorded as such).

The independent reference (split_reference.py) gives every case of
split-cases.txt the totals of each side (volume, area, centre); the native
probe (occt_split_oracle.cpp) builds each prism from
identity_reference.native_case and splits it with a large planar face. On
each side the native pieces' volumes and areas must sum to the reference's
within 2e-8 relative (BRepGProp's accuracy on elliptic faces) and their
common centre lie within 2e-8 of the case's size; the sides present must
agree. `--capture` records the native
observations, before any kernel code for the split exists except under a
post-implementation key; later runs must
reproduce them (on another platform, its reviewed record). Differences need a
fingerprinted review. The kernel's pieces (`split_probe`) must have the
reference's sides, their volumes', areas' and centres' sums inside the
kernel's enclosures, and the native piece counts per side and face, edge and
vertex counts (a difference needs a review). A case the kernel reports
`unsupported` (a plane or solid of a later sub-step of S8) is listed apart.
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
from identity_reference import native_case
import generate_split_fixtures as fixtures

SOURCE_FILE = ROOT/'rust/tools/occt_split_oracle.cpp'
REVIEWS = ROOT/'rust/fixtures/occt-split-divergences.json'
CAPTURES = {
    's8a': (ROOT/'rust/fixtures/occt-split-preimplementation', ROOT/'rust/kernel/src/solid/split.rs'),
    's8c': (ROOT/'rust/fixtures/occt-split-primitive-preimplementation',
            ROOT/'rust/kernel/src/solid/split/revolved.rs'),
    's8d': (ROOT/'rust/fixtures/occt-split-torus-preimplementation',
            ROOT/'rust/kernel/src/solid/split/torus.rs'),
    's8d2': (ROOT/'rust/fixtures/occt-split-conic-postimplementation',
             ROOT/'rust/kernel/src/solid/split/conic.rs'),
    's8d3': (ROOT/'rust/fixtures/occt-split-spiric-preimplementation',
             ROOT/'rust/kernel/src/solid/split/spiric.rs'),
}
# Captures taken after their kernel code (S8d.2's conic configurations beyond
# the three S8c captured before it): recorded as such, never as references
# made before the code.
POST_IMPLEMENTATION = {'s8d2'}
# OCCT's BRepGProp on the split pieces' elliptic faces (a plane across an arc
# or a circle's wall) errs by up to 8.8e-9 relative in the pre-implementation
# capture (disc_through_caps); planar pieces agree to 1e-15.
BOUND = 2e-8


def native_input(key='s8a'):
    blocks = []
    if key in ('s8d', 's8d3'):
        listed = fixtures.torus_cases() if key == 's8d' else fixtures.spiric_cases()
        for name, frame, params, plane in listed:
            o, x, _, n = stored_axes(frame)
            axis = ' '.join(repr(float(v)) for v in list(o)+list(n)+list(x)+list(params)+[2*math.pi])
            blocks.append(f'case {name}\ntorus {axis}\nsplit '
                          + ' '.join(repr(float(v)) for v in plane)+'\nend')
        return '\n'.join(blocks)+'\n'
    if key in ('s8c', 's8d2'):
        listed = fixtures.primitive_cases() if key == 's8c' else fixtures.conic_cases()
        for kind, name, frame, params, plane in listed:
            o, x, _, n = stored_axes(frame)
            axis = ' '.join(repr(float(v)) for v in list(o)+list(n)+list(x)+list(params))
            blocks.append(f'case {name}\n{kind} {axis}\nsplit '
                          + ' '.join(repr(float(v)) for v in plane)+'\nend')
        return '\n'.join(blocks)+'\n'
    for case, plane in fixtures.cases():
        text = native_case(case)
        body, _ = text.rsplit('\nend', 1)
        blocks.append(body+'\nsplit '+' '.join(repr(float(v)) for v in plane)+'\nend')
    return '\n'.join(blocks)+'\n'


def parse_native(stdout):
    """{case: (status, [(side, volume, area, centre, counts, valid)])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'S':
            current[1].append((w[1], float(w[2]), float(w[3]), [float(x) for x in w[4:7]],
                               [int(x) for x in w[7:10]], w[10] == '1'))
        else:
            current = (w[1], [])
            out[w[0]] = current
    return out


def expected_rows():
    out = {}
    text = (ROOT/'rust/fixtures/split-expected.tsv').read_text().splitlines()[1:] + \
        (ROOT/'rust/fixtures/split-primitive-expected.tsv').read_text().splitlines()[1:] + \
        (ROOT/'rust/fixtures/split-conic-expected.tsv').read_text().splitlines()[1:] + \
        (ROOT/'rust/fixtures/split-torus-expected.tsv').read_text().splitlines()[1:] + \
        (ROOT/'rust/fixtures/split-spiric-expected.tsv').read_text().splitlines()[1:]
    for line in text:
        name, row = line.split('\t')
        w = row.split()
        out.setdefault(name, []).append((w[1], float(w[2]), float(w[3]), [float(x) for x in w[4:7]]))
    return out


def totals(pieces):
    """{side: (volume, area, centre)} of pieces (side, volume, area, centre, ...)."""
    out = {}
    for p in pieces:
        v, a, c = out.get(p[0], (0.0, 0.0, [0.0, 0.0, 0.0]))
        out[p[0]] = (v+p[1], a+p[2], [c[i]+p[1]*p[3][i] for i in range(3)])
    return {k: (v, a, [x/v for x in c]) for k, (v, a, c) in out.items()}


def case_scale(case, plane):
    if isinstance(case, tuple):
        _, _, frame, params, _ = case
        return max([1.0]+[abs(x) for x in list(frame[:3])+list(plane[:3])+list(params)])
    return max([1.0]+[abs(x) for x in list(case.frame[:3])+list(plane[:3])]+[abs(case.end-case.start)])


def differences(case, plane, native, rows):
    status, pieces = native
    if status != 'done':
        return ['not_done']
    out = []
    if not all(p[5] for p in pieces):
        out.append('invalid_piece')
    want = {r[0]: r[1:] for r in rows}
    if 'whole' in want:
        # One piece, on whichever side its centre falls.
        got = totals([('whole',)+tuple(p[1:]) for p in pieces])
        if len(pieces) != 1:
            out.append('split_where_whole')
    else:
        got = totals(pieces)
    if set(got) != set(want):
        return sorted(set(out+['sides']))
    scale = case_scale(case, plane)
    for side, (v, a, c) in want.items():
        gv, ga, gc = got[side]
        if abs(gv-v) > BOUND*abs(v) or abs(ga-a) > BOUND*abs(a):
            out.append('measure')
        if math.dist(gc, c) > BOUND*scale:
            out.append('centre')
    return sorted(set(out))


def rust_rows():
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'split_probe'],
                   cwd=ROOT, check=True)
    text = (ROOT/'rust/fixtures/split-cases.txt').read_text() + \
        (ROOT/'rust/fixtures/split-primitive-cases.txt').read_text() + \
        (ROOT/'rust/fixtures/split-conic-cases.txt').read_text() + \
        (ROOT/'rust/fixtures/split-torus-cases.txt').read_text() + \
        (ROOT/'rust/fixtures/split-spiric-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/split_probe')], input=text,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        out.setdefault(w[0], []).append(w[1:])
    return out


def rust_differences(rust, rows, native):
    """Each side's sums of the kernel's enclosures contain the reference's
    totals (1e-20 relative slack for the reference's quadrature); the piece
    counts per side and each piece's face, edge and vertex counts equal the
    native ones."""
    out = []
    if rust and rust[0][0] == 'limit':
        return ['rust_limit']
    pieces = [(w[1], [float(x) for x in w[2:4]], [float(x) for x in w[4:6]],
               [[float(x) for x in w[6+2*i:8+2*i]] for i in range(3)], [int(x) for x in w[12:15]])
              for w in rust if w[0] == 'piece']
    want = {r[0]: r[1:] for r in rows}
    if 'whole' in want:
        if len(pieces) != 1:
            return ['rust_pieces']
        pieces = [('whole',)+p[1:] for p in pieces]
    sides = {}
    for side, vol, area, centre, counts in pieces:
        s = sides.setdefault(side, [[0.0, 0.0], [0.0, 0.0], [[0.0, 0.0] for _ in range(3)], []])
        s[0] = [s[0][0]+vol[0], s[0][1]+vol[1]]
        s[1] = [s[1][0]+area[0], s[1][1]+area[1]]
        # Moments: the centre times the volume, bounded from both ends.
        for i in range(3):
            prods = [a*b for a in vol for b in centre[i]]
            s[2][i] = [s[2][i][0]+min(prods), s[2][i][1]+max(prods)]
        s[3].append(counts)
    if set(sides) != set(want):
        return ['rust_sides']
    for side, (v, a, c) in want.items():
        (vlo, vhi), (alo, ahi), moments, _ = sides[side]
        slack = lambda x: 1e-20*abs(x)
        if not (vlo-slack(v) <= v <= vhi+slack(v)) or not (alo-slack(a) <= a <= ahi+slack(a)):
            out.append('rust_measure_outside_reference')
        for i in range(3):
            m = v*c[i]
            if not (moments[i][0]-1e-12*max(1.0, abs(m)) <= m <= moments[i][1]+1e-12*max(1.0, abs(m))):
                out.append('rust_centre_outside_reference')
    if native is not None and native[0] == 'done':
        mine = sorted((side, tuple(counts)) for side, _, _, _, counts in pieces)
        theirs = sorted(('whole' if 'whole' in want else p[0], tuple(p[4])) for p in native[1])
        if [m[0] for m in mine] != [t[0] for t in theirs]:
            out.append('rust_piece_count')
        elif mine != theirs:
            out.append('rust_entity_counts')
    return sorted(set(out))


def capture(executable, env, key, sdk_manifest):
    CAPTURE, kernel_file = CAPTURES[key]
    text = native_input(key)
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native split run failed: '+json.dumps(record)[:2000])
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
        'rust_split_exists': kernel_file.exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed_by_key):
    for key, (CAPTURE, _) in CAPTURES.items():
        observed = observed_by_key[key]
        metadata = json.loads((CAPTURE/'capture.json').read_text())
        if metadata['source_reference'] != SOURCE \
                or metadata['rust_split_exists'] != (key in POST_IMPLEMENTATION):
            raise ValueError('split capture was not the recorded pre- or post-implementation one')
        for field, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                            ('observations_sha256', 'native.txt')]:
            if metadata[field] != digest(CAPTURE/name):
                raise ValueError('split evidence changed: '+name)
        if (CAPTURE/'inputs.txt').read_text() != native_input(key):
            raise ValueError('the native inputs differ from the captured ones')
        was = parse_native(platform_record(CAPTURE, metadata))
        if set(was) != set(observed):
            raise ValueError('native cases differ from the capture')
        for name, (status, pieces) in was.items():
            now = observed[name]
            flat = lambda ps: [x for p in ps for x in [p[1], p[2], *p[3]]]
            if now[0] != status or len(now[1]) != len(pieces) \
                    or [p[4] for p in now[1]] != [p[4] for p in pieces] \
                    or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now[1]), flat(pieces))):
                raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/split-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', choices=sorted(CAPTURES))
    parser.add_argument('--native-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'split-oracle', ('TKBO', 'TKPrim'),
                                             ['TKBO', 'TKPrim']+TOOLKITS)
    if args.capture:
        capture(executable, env, args.capture, args.sdk_manifest)
        print('captured', args.capture)
        return
    by_key, stdout = {}, []
    for key in CAPTURES:
        record = run(executable, native_input(key), env)
        if record['exit_code'] != 0:
            raise SystemExit('native split run failed: '+json.dumps(record)[:2000])
        stdout.append(record['stdout'])
        by_key[key] = parse_native(record['stdout'])
    (output/'native-observed.txt').write_text(''.join(stdout))
    captured(by_key)
    observed = {k: v for d in by_key.values() for k, v in d.items()}
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    expected = expected_rows()
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'rust_within_reference': 0,
              'rust_unsupported': [], 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = None if args.native_only else rust_rows()
    everything = [(case, plane, case.name) for case, plane in fixtures.cases()] + \
        [(c, c[4], c[1]) for c in fixtures.primitive_cases()] + \
        [(c, c[4], c[1]) for c in fixtures.conic_cases()] + \
        [(('torus', c[0], c[1], c[2], c[3]), c[3], c[0]) for c in fixtures.torus_cases()] + \
        [(('torus', c[0], c[1], c[2], c[3]), c[3], c[0]) for c in fixtures.spiric_cases()]
    for case, plane, name in everything:
        report['cases'] += 1
        rows = expected[name]
        native = observed[name]
        found = differences(case, plane, native, rows)
        if rust is not None and rust[name] == [['unsupported']]:
            # A plane or solid of a later sub-step (S8a.2 on): listed, and
            # the native comparison still made.
            report['rust_unsupported'].append(name)
        elif rust is not None:
            wrong = rust_differences(rust[name], rows, native)
            if any(w in ('rust_limit', 'rust_pieces', 'rust_sides', 'rust_measure_outside_reference',
                         'rust_centre_outside_reference') for w in wrong):
                report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
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
            review or dict(evidence, pieces=native[1]))
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
