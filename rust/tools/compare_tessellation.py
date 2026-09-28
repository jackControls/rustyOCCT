#!/usr/bin/env python3
"""Source-pinned BRepMesh_IncrementalMesh observations beside the tessellation reference (T-a).

The bodies of tessellation-cases.txt are written as `.brep` text by the
kernel's existing writer (brep_io_probe `prisms` and `bodies`), read
natively and meshed by BRepMesh_IncrementalMesh at each case's settings
(occt_tessellation_oracle.cpp): node and triangle counts, OCCT's own
deflections, the deflection measured by sampling against each face's
surface, watertightness after joining nodes through OCCT's edge polygons,
orientation, area and volume. `--capture` records these native rows before
any kernel tessellation code exists; later runs must reproduce them, and
must write the same `.brep` texts. The native meshes are also checked by the
independent reference (tessellation_reference.py): where OCCT exceeds the
requested deflection, leaves gaps, misorients a triangle or understates its
own deflection, the difference needs a fingerprinted review. The kernel's
meshes (`tessellation_probe`) must pass every check of the reference on
every case and setting; there is no review for them.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from build_pinned_occt import SOURCE, digest
from compare_brep import NUMBER_DRIFT, review_for, run, sha, write
from compare_brep_io import build
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
import generate_tessellation_fixtures as fixtures
import tessellation_reference as ref

SOURCE_FILE = ROOT/'rust/tools/occt_tessellation_oracle.cpp'
CAPTURE = ROOT/'rust/fixtures/occt-tessellation-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-tessellation-divergences.json'
KERNEL_FILE = ROOT/'rust/kernel/src/tessellation.rs'
CASES = ROOT/'rust/fixtures/tessellation-cases.txt'
TOOLKITS = ['TKMesh', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath',
            'TKernel']
# Native rows on every run: statuses and counts exact, reals within this
# relative bound (or 1e-15 absolute), as listed in VALIDATION.md's allowance
# table: BRepMesh and the probe's projections use platform trigonometry.
REAL_BOUND = 1e-9
FIELDS = ('status', 'faces', 'unmeshed', 'nodes', 'welded', 'triangles', 'degenerate', 'free',
          'nonmanifold', 'misoriented', 'weld_gap', 'face_deflection', 'edge_deflection',
          'face_distance', 'edge_distance', 'area', 'volume')


def split_cases(text):
    """[(name, identity block without mesh rows, [(setting, deflection, angle)])]."""
    out = []
    for block in text.split('\nend'):
        if not block.strip():
            continue
        rows = block.strip().splitlines()
        settings = [(w[1], float(w[2]), float(w[3])) for w in (r.split() for r in rows) if w[0] == 'mesh']
        body = '\n'.join(r for r in rows if not r.startswith('mesh '))
        out.append((rows[0].split()[1], body+'\nend\n', settings))
    return out


def case_objects():
    return {c.name: c for c in fixtures.cases()}


def brep_texts():
    """{case: .brep text} from the kernel's writer."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'brep_io_probe'],
                   cwd=ROOT, check=True)
    probe = str(ROOT/'target/release/examples/brep_io_probe')
    cases = split_cases(CASES.read_text())
    objects = case_objects()
    solids = ''.join(b for n, b, _ in cases if objects[n].make is None)
    faces = ''.join(b for n, b, _ in cases if objects[n].make is not None)
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        for mode, text in (('prisms', solids), ('bodies', faces)):
            subprocess.run([probe, mode, tmp], input=text, text=True, capture_output=True, check=True,
                           timeout=600)
        for name, _, _ in cases:
            out[name] = (Path(tmp)/f'{name}.brep').read_text()
    return out


def native_input(texts):
    rows = []
    cases = split_cases(CASES.read_text())
    for name, _, _ in cases:
        rows.append(f'body {name}')
        rows.append(texts[name].rstrip('\n'))
        rows.append('endbody')
    for name, _, settings in cases:
        for setting, deflection, angle in settings:
            rows.append(f'mesh {name} {setting} {deflection!r} {angle!r}')
    return '\n'.join(rows)+'\n'


def parse_native(stdout):
    """({(case, setting): row dict}, dump text)."""
    rows, dump = {}, []
    for line in stdout.splitlines():
        w = line.split()
        if w and w[0] == 'S':
            values = w[3:]
            row = {'status': values[0]}
            for key, value in zip(FIELDS[1:], values[1:]):
                row[key] = float(value) if '.' in value or 'e' in value else int(value)
            rows[(w[1], w[2])] = row
        else:
            dump.append(line)
    return rows, '\n'.join(dump)+'\n'


COUNTS = ('faces', 'unmeshed', 'nodes', 'welded', 'triangles', 'degenerate', 'free', 'nonmanifold',
          'misoriented')


def same_row(was, now):
    if set(was) != set(now) or was['status'] != now['status']:
        return False
    for key, a in was.items():
        b = now[key]
        if key == 'status':
            continue
        if key in COUNTS:
            if a != b:
                return False
        elif abs(a-b) > REAL_BOUND*max(abs(a), abs(b)) and abs(a-b) > 1e-15:
            return False
    return True


def capture(executable, env, sdk_manifest, text):
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native tessellation run failed: '+json.dumps(record)[:2000])
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
        'rust_tessellation_exists': KERNEL_FILE.exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def drifted_inputs(captured, current):
    """The same tokens, numbers within `NUMBER_DRIFT` relative to their size
    (at least 1): the kernel writes them with the platform's libm, whose
    rounding moves a value by an ulp or two of its own magnitude (one Linux
    coordinate -9.999999999999998 for macOS's -10)."""
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
        if not abs(x-y) <= NUMBER_DRIFT*max(1.0, abs(x)):
            raise ValueError(f'current native corpus moved {a} to {b}')


def captured(observed, text):
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_tessellation_exists']:
        raise ValueError('tessellation capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('tessellation evidence changed: '+name)
    if (CAPTURE/'oracle.cpp').read_text() != SOURCE_FILE.read_text():
        raise ValueError('the native probe differs from the captured one')
    # The kernel writes the .brep texts with the platform's trigonometry
    # (frames, rotations): structure exact, numbers within compare_brep's
    # NUMBER_DRIFT (VALIDATION.md's allowances).
    if (CAPTURE/'inputs.txt').read_text() != text:
        try:
            drifted_inputs((CAPTURE/'inputs.txt').read_text(), text)
        except ValueError as e:
            raise ValueError('the kernel\'s .brep texts or the settings differ from the captured ones: '
                             + str(e)) from None
    was, _ = parse_native((CAPTURE/'native.txt').read_text())
    if set(was) != set(observed):
        raise ValueError('native cases differ from the capture')
    for key, row in was.items():
        if not same_row(row, observed[key]):
            raise ValueError(f'native observation of {key} differs from the capture: {observed[key]}')


def rust_meshes():
    """{(case, setting): Mesh | error text} from the kernel's probe."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'tessellation_probe'],
                   cwd=ROOT, check=True)
    out = subprocess.run([str(ROOT/'target/release/examples/tessellation_probe')], input=CASES.read_text(),
                         text=True, capture_output=True, timeout=1800, check=True).stdout
    return ref.parse_meshes(out)


def native_differences(case, deflection, angle, row, mesh):
    """What separates OCCT's mesh from the request and the reference: the
    reference's checks of the welded native mesh (a face body's free edges
    are its boundary), and OCCT's degenerate triangles."""
    if row['status'] != 'done' or row['unmeshed']:
        return ['not_meshed'], {}
    measured, failures = ref.check(case, deflection, angle, mesh)
    names = [('not_closed', 'occt_not_watertight'), ('boundary_loops', 'occt_not_watertight'),
             ('non_manifold', 'occt_not_watertight'), ('misoriented', 'occt_misoriented'),
             ('deflection_exceeded', 'occt_deflection_exceeded'),
             ('edge_deflection', 'occt_deflection_exceeded'), ('off_face', 'occt_off_face'),
             ('bound_unsound', 'occt_deflection_understated'), ('normal_inward', 'occt_normal_inward'),
             ('volume', 'occt_volume'), ('euler', 'occt_euler'),
             ('node_off_boundary', 'occt_node_off_boundary')]
    found = sorted({new for old, new in names if old in failures})
    if row['degenerate']:
        found.append('occt_degenerate_triangles')
    return found, measured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/tessellation-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', action='store_true',
                        help='record the native observations (before the kernel code exists only)')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    for name, text in fixtures.generate().items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'tessellation-oracle',
                                             ('TKMesh', 'TKBRep'), TOOLKITS)
    text = native_input(brep_texts())
    if args.capture:
        if KERNEL_FILE.exists():
            raise SystemExit('the kernel\'s tessellation exists: a capture now would not precede it')
        capture(executable, env, args.sdk_manifest, text)
        print('captured', len(parse_native((CAPTURE/'native.txt').read_text())[0]), 'native meshes')
        return
    record = run(executable, text, env)
    dumped = subprocess.run([str(executable), 'dump'], input=text, text=True, capture_output=True,
                            timeout=1800, env=env)
    if record['exit_code'] != 0 or dumped.returncode != 0:
        raise SystemExit('native tessellation run failed: '+json.dumps(record)[:2000])
    observed, _ = parse_native(record['stdout'])
    captured(observed, text)
    native_meshes = ref.parse_meshes(parse_native(dumped.stdout)[1])
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'meshes': 0,
              'rust_passes_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': [],
              'counts': {}}
    rust = rust_meshes()
    objects = case_objects()
    for name, _, settings in split_cases(CASES.read_text()):
        case = objects[name]
        for setting, deflection, angle in settings:
            key = f'{name}/{setting}'
            report['meshes'] += 1
            mesh = rust.get((name, setting))
            if not isinstance(mesh, ref.Mesh):
                report['failures'].append({'case': key, 'reason': 'rust_error', 'rust': mesh})
                continue
            measured, failures = ref.check(case, deflection, angle, mesh)
            if failures:
                report['failures'].append({'case': key, 'reason': ' '.join(failures), 'rust': measured})
                continue
            report['rust_passes_reference'] += 1
            row = observed[(name, setting)]
            found, native = native_differences(case, deflection, angle, row, native_meshes[(name, setting)])
            report['counts'][key] = {
                'rust': [measured['nodes'], measured['triangles'], mesh.deflection],
                'occt': [row['welded'], row['triangles'], native.get('deflection')],
                'rust_area_error': measured['area_error'], 'occt_area_error': native.get('area_error'),
                'rust_volume_error': measured['volume_error'],
                'occt_volume_error': native.get('volume_error')}
            if not found:
                report['matches'].append(key)
                continue
            # The fingerprint is this run's row: a platform whose row
            # differs in the last bits gets its own review.
            evidence = {'case': key, 'source_reference': SOURCE, 'oracle': oracle,
                        'native_sha256': sha(json.dumps(row, sort_keys=True)), 'differences': found}
            review = review_for(evidence, reviews)
            report['reviewed_differences' if review else 'failures'].append(
                review or dict(evidence, native=native, row=row))
    write(output/'capture.json', {'source_reference': SOURCE, 'oracle': oracle,
                                  'sdk_manifest_sha256': digest(args.sdk_manifest),
                                  'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                                  'loaded_libraries': loaded, 'build_command': command})
    write(output/'report.json', report)
    summary = {k: len(v) if isinstance(v, (list, dict)) else v for k, v in report.items()}
    print(json.dumps(summary, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
