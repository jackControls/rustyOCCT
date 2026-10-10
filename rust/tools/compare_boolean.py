#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9a reference: Booleans of two prisms in one frame.

The independent reference (boolean_reference.py) gives every case of
boolean-cases.txt its result's solid count and totals (volume, area,
centre), or `empty`; the native probe (occt_boolean_oracle.cpp) builds both
prisms from identity_reference.native_boolean_case and runs the operation.
The native result must be valid (the whole result and each solid under
BRepCheck_Analyzer), have the reference's number of solids, and its solids'
volumes and areas must sum to the reference's within 2e-8 relative (split's
allowance, BRepGProp's accuracy on curved faces) and their common centre lie
within 2e-8 of the case's size. `--capture` records the native observations,
before any kernel code for the Booleans exists (`rust/kernel/src/solid/
boolean.rs` absent); later runs must reproduce them (on another platform,
its reviewed record `platform-<name>/`). Differences need a fingerprinted
review in occt-boolean-divergences.json.

`--splines` runs S9a.2's spline set instead (`boolean-spline-cases.txt`,
`generate_boolean_fixtures.spline_cases`, the oracle building each B-spline
segment as a Geom_BSplineCurve edge) against `boolean-spline-expected.tsv`
and its capture `occt-boolean-spline-preimplementation` (reviews in
occt-boolean-spline-divergences.json), recorded while the
kernel refuses spline profiles (`OutOfDomain("a Boolean of spline profiles
(S9a.2)")` in `rust/kernel/src/profile/boolean.rs`: the capture's
`rust_spline_boolean_exists` is false). While that refusal stands, the probe
must report `unsupported` on every spline case (a probe failure or any
other row is a failure) and the cases are listed under `rust_unsupported`;
once it is gone, the kernel's rows are compared as S9a's.

The kernel's side runs only when `rust/kernel/examples/boolean_probe.rs`
exists; until then every case is listed under `rust_unsupported` and the
report's `rust_probe_exists` is false. The probe reads the cases on stdin and
prints per case `NAME ROW` lines: `limit`, `unsupported` (`OutOfDomain`: a
stack before S9a.2), `refused` (a documented `Degenerate`), `empty` (no
solids), or per result solid `solid vol_lo vol_hi area_lo area_hi cx_lo
cx_hi cy_lo cy_hi cz_lo cz_hi faces edges vertices` (the enclosures of its
volume, area and centre, and `Topology::occt_counts`). The solids' sums of
enclosures must contain the reference's totals (1e-20 relative slack), their
number must be the reference's, `refused` is accepted only where the case
expects `degenerate`, and each solid's counts are compared with the native
solid's after ShapeUpgrade_UnifySameDomain (a difference needs a review).
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from build_pinned_occt import SOURCE, digest
import compare_brep
from compare_brep import review_for, run, sha, write

# The native Boolean runs of the larger sets (cones, tori, spline walls) take
# close to compare_brep's 120 s on a loaded host; the probe already has 600.
compare_brep.TIMEOUT = max(compare_brep.TIMEOUT, 600)
from compare_brep_io import TOOLKITS, build
from compare_curve_surface import platform_record
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
from identity_reference import Spline, native_boolean_case
import generate_boolean_fixtures as fixtures

SOURCE_FILE = ROOT/'rust/tools/occt_boolean_oracle.cpp'
KERNEL_FILE = ROOT/'rust/kernel/src/solid/boolean.rs'
PROFILE_BOOLEAN = ROOT/'rust/kernel/src/profile/boolean.rs'
SPLINE_REFUSAL = 'OutOfDomain("a Boolean of spline profiles (S9a.2)")'
PROBE = ROOT/'rust/kernel/examples/boolean_probe.rs'
# Split's allowance: BRepGProp's error on curved faces.
BOUND = 2e-8


class Set:
    """A fixture set: S9a's (the default) or S9a.2's spline profiles. A
    wrapper running another set (`compare_polyhedral.py`,
    `compare_curved_boolean.py`) replaces `make_set` with a subclass giving
    its cases, expected rows, capture, reviews, the capture's flag that the
    kernel code did not exist (`exists_key`, `exists`) and whether the kernel
    must still refuse every case (`before_code`)."""

    def __init__(self, splines):
        self.splines = splines
        self.cases = fixtures.spline_cases if splines else fixtures.cases
        self.expected = ROOT/('rust/fixtures/boolean-spline-expected.tsv' if splines
                              else 'rust/fixtures/boolean-expected.tsv')
        self.capture = ROOT/('rust/fixtures/occt-boolean-spline-preimplementation' if splines
                             else 'rust/fixtures/occt-boolean-preimplementation')
        self.output = ROOT/('target/boolean-spline-oracle' if splines else 'target/boolean-oracle')
        self.reviews = ROOT/('rust/fixtures/occt-boolean-spline-divergences.json' if splines
                             else 'rust/fixtures/occt-boolean-divergences.json')
        self.exists_key = 'rust_spline_boolean_exists' if splines else 'rust_boolean_exists'

    def exists(self):
        """Whether the kernel code the capture came before exists now."""
        return rust_spline_boolean_exists() if self.splines else KERNEL_FILE.exists()

    def before_code(self):
        """True while the kernel must refuse every case (`unsupported`)."""
        return self.splines and not rust_spline_boolean_exists()

    def must_support(self, case):
        """Whether the kernel may no longer report `case` `unsupported` (a
        later sub-step's case once its code exists, S9f.2b.2)."""
        return False

    def refused_before_code(self, case, rows):
        """Whether the kernel's rows for `case` before its code are the
        refusal the capture came before (`unsupported`; S9e.4b.3c.2 also
        takes an S9 refusal of a declared degenerate case raised first)."""
        return rows == [['unsupported']]


def make_set(splines):
    return Set(splines)


SET = Set(False)


def rust_spline_boolean_exists():
    """False while the kernel refuses every Boolean of spline profiles."""
    return SPLINE_REFUSAL not in PROFILE_BOOLEAN.read_text()


def native_input():
    return '\n'.join(native_boolean_case(c.obj, c.operation, c.tool) for c in SET.cases())+'\n'


def parse_native(stdout):
    """{case: (status, valid, warnings, [(volume, area, centre, counts, unified counts, valid)])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'S':
            current[3].append((float(w[1]), float(w[2]), [float(x) for x in w[3:6]],
                               [int(x) for x in w[6:9]], [int(x) for x in w[9:12]], w[12] == '1'))
        else:
            current = (w[1], w[3] == '1', w[4] == '1', [])
            out[w[0]] = current
    return out


def expected_rows():
    """{case: {'expect': (kind, step), 'result': (N, volume, area, centre) or None, 'slabs': [...]}}."""
    out = {}
    for line in SET.expected.read_text().splitlines()[1:]:
        name, row = line.split('\t')
        w = row.split()
        e = out.setdefault(name, {'result': None, 'slabs': []})
        if w[0] == 'expect':
            e['expect'] = (w[1], w[2])
        elif w[0] == 'result':
            e['result'] = (int(w[1]), float(w[2]), float(w[3]), [float(x) for x in w[4:7]])
        elif w[0] == 'slab':
            e['slabs'].append([float(x) for x in w[1:]])
    return out


def case_scale(case):
    values = [1.0, abs(case.obj.start), abs(case.obj.end), abs(case.tool.start), abs(case.tool.end)]
    values += [abs(x) for x in case.obj.frame[:3]+case.tool.frame[:3]]
    # S9d.1: a sphere's reach from the origin.
    values += [max(abs(x) for x in c.frame[:3])+c.sphere[0] for c in (case.obj, case.tool)
               if c.sphere is not None]
    # S9d.3a: a cone's reach from the origin.
    values += [max(abs(x) for x in c.frame[:3])+max(c.cone[0], c.cone[1])+c.cone[2] for c in (case.obj, case.tool)
               if c.cone is not None]
    # S9d.4a: a torus's reach from the origin.
    values += [max(abs(x) for x in c.frame[:3])+c.torus[0]+c.torus[1] for c in (case.obj, case.tool)
               if c.torus is not None]
    for b in case.obj.boundaries+case.tool.boundaries:
        if b.circle is not None:
            values += [abs(b.circle[0])+b.circle[2], abs(b.circle[1])+b.circle[2]]
        else:
            values += [abs(x) for p in b.points for x in p]
            values += [abs(x) for s in b.segments or [] if isinstance(s, Spline) for p in s.poles for x in p]
    return max(values)


def totals(solids):
    V = sum(s[0] for s in solids)
    A = sum(s[1] for s in solids)
    return V, A, [sum(s[0]*s[2][i] for s in solids)/V for i in range(3)]


def differences(case, native, expected):
    status, valid, _, solids = native
    if status != 'done':
        return [status]
    out = []
    if not valid:
        out.append('invalid_result')
    if not all(s[5] for s in solids):
        out.append('invalid_solid')
    want = expected['result']
    if want is None:
        return sorted(set(out+(['solids'] if solids else [])))
    count, V, A, c = want
    if len(solids) != count:
        out.append('solids')
    if not solids:
        return sorted(set(out+['solids']))
    gv, ga, gc = totals(solids)
    if abs(gv-V) > BOUND*abs(V) or abs(ga-A) > BOUND*abs(A):
        out.append('measure')
    if math.dist(gc, c) > BOUND*case_scale(case):
        out.append('centre')
    return sorted(set(out))


def rust_rows():
    """{case: [row words]} from the probe, one case at a time (a probe that
    fails on a case or prints an unknown row reports it `unsupported`), and
    the cases it failed on."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'boolean_probe'],
                   cwd=ROOT, check=True)
    out, failed = {}, []
    for case in SET.cases():
        run_ = subprocess.run([str(ROOT/'target/release/examples/boolean_probe')], input=case.encode()+'\n',
                              text=True, capture_output=True, timeout=600)
        rows = [line.split() for line in run_.stdout.splitlines()]
        if run_.returncode != 0 or not rows or any(
                w[0] != case.name or w[1] not in ('solid', 'empty', 'limit', 'unsupported', 'refused') for w in rows):
            failed.append(case.name)
            out[case.name] = [['unsupported']]
        else:
            out[case.name] = [w[1:] for w in rows]
    return out, failed


def rust_differences(rust, expected, native):
    """The kernel's solids against the reference's totals and the native
    unified counts."""
    kind = expected['expect'][0]
    if rust[0][0] == 'limit':
        return ['rust_limit']
    if rust[0][0] == 'refused':
        return ['rust_refused'] if kind == 'degenerate' else ['rust_refused_unexpectedly']
    want = expected['result']
    solids = [([float(x) for x in w[1:3]], [float(x) for x in w[3:5]],
               [[float(x) for x in w[5+2*i:7+2*i]] for i in range(3)], [int(x) for x in w[11:14]])
              for w in rust if w[0] == 'solid']
    if want is None:
        return [] if rust[0][0] == 'empty' and not solids else ['rust_solids']
    count, V, A, c = want
    if len(solids) != count:
        return ['rust_solids']
    vlo = sum(s[0][0] for s in solids)
    vhi = sum(s[0][1] for s in solids)
    alo = sum(s[1][0] for s in solids)
    ahi = sum(s[1][1] for s in solids)
    out = []
    slack = lambda x: 1e-20*abs(x)
    if not (vlo-slack(V) <= V <= vhi+slack(V)) or not (alo-slack(A) <= A <= ahi+slack(A)):
        out.append('rust_measure_outside_reference')
    for i in range(3):
        lo = sum(min(a*b for a in s[0] for b in s[2][i]) for s in solids)
        hi = sum(max(a*b for a in s[0] for b in s[2][i]) for s in solids)
        m = V*c[i]
        allow = 1e-12*max(1.0, abs(m))
        if not (lo-allow <= m <= hi+allow):
            out.append('rust_centre_outside_reference')
    if native is not None and native[0] == 'done':
        if sorted(tuple(s[3]) for s in solids) != sorted(tuple(s[4]) for s in native[3]):
            out.append('rust_entity_counts')
    return sorted(set(out))


def capture(executable, env, sdk_manifest):
    text = native_input()
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native Boolean run failed: '+json.dumps(record)[:2000])
    CAPTURE = SET.capture
    CAPTURE.mkdir(parents=True, exist_ok=True)
    (CAPTURE/'inputs.txt').write_text(text)
    (CAPTURE/'native.txt').write_text(record['stdout'])
    (CAPTURE/'oracle.cpp').write_text(SOURCE_FILE.read_text())
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    exists = {SET.exists_key: SET.exists()}
    write(CAPTURE/'capture.json', {
        'source_reference': SOURCE, 'oracle': next(iter(record['stderr'].splitlines()), None),
        'platform': sys.platform, 'rust_revision': revision, **exists,
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    CAPTURE = SET.capture
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    exists = metadata[SET.exists_key]
    if metadata['source_reference'] != SOURCE or exists:
        raise ValueError('Boolean capture was not the recorded pre-implementation one')
    for field, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                        ('observations_sha256', 'native.txt')]:
        if metadata[field] != digest(CAPTURE/name):
            raise ValueError('Boolean evidence changed: '+name)
    if (CAPTURE/'inputs.txt').read_text() != native_input():
        raise ValueError('the native inputs differ from the captured ones')
    was = parse_native(platform_record(CAPTURE, metadata))
    if set(was) != set(observed):
        raise ValueError('native cases differ from the capture')
    flat = lambda ss: [x for s in ss for x in [s[0], s[1], *s[2]]]
    for name, (status, valid, warnings, solids) in was.items():
        now = observed[name]
        if now[:3] != (status, valid, warnings) or len(now[3]) != len(solids) \
                or [s[3:] for s in now[3]] != [s[3:] for s in solids] \
                or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now[3]), flat(solids))):
            raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', action='store_true')
    parser.add_argument('--native-only', action='store_true')
    parser.add_argument('--splines', action='store_true', help="S9a.2's spline set")
    args = parser.parse_args()
    global SET
    SET = make_set(args.splines)
    output = (args.output or SET.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'boolean-oracle',
                                             ('TKBO', 'TKPrim', 'TKShHealing'),
                                             ['TKBO', 'TKPrim', 'TKShHealing']+TOOLKITS)
    if args.capture:
        capture(executable, env, args.sdk_manifest)
        print('captured')
        return
    record = run(executable, native_input(), env)
    if record['exit_code'] != 0:
        raise SystemExit('native Boolean run failed: '+json.dumps(record)[:2000])
    (output/'native-observed.txt').write_text(record['stdout'])
    observed = parse_native(record['stdout'])
    captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not SET.reviews.exists() else \
        json.loads(SET.reviews.read_text())['reviews']
    expected = expected_rows()
    report = {'source_reference': SOURCE, 'oracle': oracle, 'cases': 0, 'rust_probe_exists': PROBE.exists(),
              'rust_within_reference': 0, 'rust_unsupported': [], 'rust_probe_failed': [], 'rust_refused': [],
              'matches': [], 'reviewed_differences': [], 'failures': []}
    # A set before its kernel code (S9a.2's splines, S9c's arcs in frames
    # with different axes): every case unsupported.
    pre_splines = SET.before_code()
    if SET.exists_key != 'rust_boolean_exists':
        report[SET.exists_key] = SET.exists()
    rust = None
    if not args.native_only and PROBE.exists():
        rust, report['rust_probe_failed'] = rust_rows()
    for case in SET.cases():
        name = case.name
        report['cases'] += 1
        native = observed[name]
        found = differences(case, native, expected[name])
        if pre_splines and rust is not None and (not SET.refused_before_code(case, rust[name])
                                                 or name in report['rust_probe_failed']):
            reason = SET.exists_key[:-len('_exists')]+'_before_its_code'
            report['failures'].append({'case': name, 'reason': reason, 'rust': rust[name]})
            continue
        if rust is not None and rust[name] == [['unsupported']] and SET.must_support(case):
            report['failures'].append({'case': name, 'reason': 'rust_unsupported_after_its_code',
                                       'rust': rust[name]})
            continue
        if rust is None or rust[name] == [['unsupported']]:
            # No probe yet, or a case of a later sub-step: listed, and the
            # native comparison still made.
            report['rust_unsupported'].append(name)
        else:
            wrong = rust_differences(rust[name], expected[name], native)
            if wrong == ['rust_refused']:
                report['rust_refused'].append(name)
            elif any(w != 'rust_entity_counts' for w in wrong):
                report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
                continue
            else:
                report['rust_within_reference'] += 1
                found = sorted(set(found+wrong))
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps(native)), 'differences': found}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, solids=native[3]))
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
