#!/usr/bin/env python3
"""Source-pinned BRepTools::Read observations beside the kernel's .brep interop (T2).

The independent reader (brep_io_reference.py) is first certified against its
committed fixture and against native OCCT: for every data/occ file, the
solids OCCT reaches from the root and their distinct subshape counts are
exactly the reader's. The kernel then writes

* every identity prism: OCCT must read it back as one BRepCheck-valid solid
  whose counts equal the kernel's synthesized counts and whose volume,
  surface area and centroid equal the kernel's exact mass properties;
* every data/occ solid it imports: OCCT must read the written text as one
  valid solid with the original's counts (as the reader states them) and the
  original's volume, surface area and centroid;
* (S6) every free shell, face, wire, edge and vertex of data/occ it imports,
  with OCCT's counts of the original and a certified area or length and
  centre containing OCCT's (rust/fixtures/occt-free-shape-capture), and every
  face and wire body of identity-sheet-cases.txt: OCCT must read each written
  file as one valid free shape with those counts and properties.

A free face whose certified enclosure excludes OCCT's default measure (its
fixed Gauss rule errs by up to 7e-6 relative on the corpus's spline faces)
is a reviewed difference, never a match: a second probe
(occt_free_shape_adaptive_oracle.cpp, the same SDK) integrates the file's
faces adaptively, and the review, fingerprinted by the original file and
both native rows, must record an independent high-precision measure that
the enclosure contains (to binary64 rounding).

The native observations of the unmodified upstream files were captured
after the kernel's reader existed (rust/fixtures/occt-brep-io-capture/NOTES.md);
they depend on no Rust output and must stay unchanged. Differences need a
fingerprinted review; timeouts, crashes and malformed output cannot be
reviewed.

The public dataset's version-3 texts (`CASCADE Topology V3, (c) Open
Cascade`, HANDOFF's item 3) were read natively by both probes, unchanged,
while the reader still refused their header (rust/fixtures/occt-brep-v3-capture,
`--capture-v3 --v3-dataset DIR`); the dataset is local-only (U1). With
`--v3-dataset DIR` they are compared too (`v3_*` in the report): the native
reading as captured, the independent reader's counts native's, every file
read by the kernel, every solid and free shape it imports with native's
counts and a free shape's certified measure enclosing native's, and each
written text read natively as one valid shape of the original's counts and
properties. A shape the kernel imports but cannot measure or write is
listed (`v3_not_measured`, `v3_unwritable`).
"""
import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from build_pinned_occt import SOURCE, digest
from compare_brep import review_for, run, sha, write
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
from compare_occt import ROOT
import brep_io_reference
import generate_brep_io_fixtures

CAPTURE = ROOT/'rust/fixtures/occt-brep-io-capture'
REVIEWS = ROOT/'rust/fixtures/occt-brep-io-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_brep_io_oracle.cpp'
CORPUS = ROOT/'data/occ'
TOOLKITS = ['TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath', 'TKernel']
# S6: the corpus's free shapes, read natively before any kernel code imports
# them (a separate probe: the T2 capture pins SOURCE_FILE's text).
FREE_SOURCE = ROOT/'rust/tools/occt_free_shape_oracle.cpp'
FREE_CAPTURE = ROOT/'rust/fixtures/occt-free-shape-capture'
# The same probe with adaptive surface integration (BRepGProp Eps 1e-12), run
# only on files holding a face whose enclosure excludes the default measure.
ADAPTIVE_SOURCE = ROOT/'rust/tools/occt_free_shape_adaptive_oracle.cpp'
# The public dataset's version-3 texts (HANDOFF's item 3): read natively
# with both probes above, unchanged, before the reader took their header.
# The dataset is local-only (U1): the capture holds names, digests and
# native observations, never a file's text.
V3_CAPTURE = ROOT/'rust/fixtures/occt-brep-v3-capture'
V3_LINE = b'CASCADE Topology V3, (c) Open Cascade'
# Both native runs together took under a second; the deadline only bounds hangs.
TIMEOUT = 600
# The worst observed difference was 2.5e-14: volume and area relative to
# themselves, centroids relative to the cube root of the volume.
PROPERTY_BOUND = 1e-11
# S6: OCCT integrates spline faces with a relative error near 1e-9
# (occt-sheet-preimplementation/NOTES.md); the kernel's enclosures must
# contain its areas, lengths and centres up to that.
MEASURE_BOUND = 1e-9


def build(prefix, output, source=SOURCE_FILE, name='oracle', required=('TKTopAlgo', 'TKBRep'), toolkits=TOOLKITS):
    include, lib = prefix/'include/opencascade', prefix/'lib'
    executable = output/name
    command = shlex.split(os.environ.get('CXX', 'c++'))+[
        '-std=c++17', '-O2', str(source), '-I'+str(include), '-L'+str(lib),
        '-Wl,-rpath,'+str(lib), '-o', str(executable)]+['-l'+name for name in toolkits]
    built = subprocess.run(command, text=True, capture_output=True)
    (output/'native-build.log').write_text(built.stdout+built.stderr)
    built.check_returncode()
    env = os.environ.copy()
    for key in ['LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']:
        env[key] = str(lib)+(os.pathsep+env[key] if env.get(key) else '')
    if sys.platform == 'darwin':
        loaded_text = subprocess.run([str(executable)], input='', text=True, capture_output=True,
                                     timeout=TIMEOUT, check=True,
                                     env=dict(env, DYLD_PRINT_LIBRARIES='1')).stderr
    else:
        loaded_text = subprocess.run(['ldd', str(executable)], text=True, capture_output=True,
                                     check=True, env=env).stdout
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    # The toolkits the probe uses (a linker with --as-needed drops the rest).
    for toolkit in required:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded library: '+toolkit)
    return executable, env, loaded, command


def parse(stdout, paths):
    """{path: [(verdict, counts, [volume, area, cx, cy, cz])]} in input order."""
    out, current = {}, None
    for line in stdout.splitlines():
        words = line.split()
        if words[0] == 'F':
            if len(words) != 3 or not words[2].isdigit():
                raise ValueError('native could not read '+line)
            current = out.setdefault(words[1], [])
        elif words[0] == 'S' and current is not None and len(words) == 13:
            current.append((words[1], tuple(int(w) for w in words[2:8]), [float(w) for w in words[8:]]))
        else:
            raise ValueError('malformed native line '+line)
    if list(out) != [str(p) for p in paths]:
        raise ValueError('native output does not cover every input in order')
    return out


def parse_free(stdout, paths):
    """{path: [(type, verdict, counts, [mass, cx, cy, cz])]} in input order."""
    out, current = {}, None
    for line in stdout.splitlines():
        words = line.split()
        if words[0] == 'F':
            if len(words) != 3 or not words[2].isdigit():
                raise ValueError('native could not read '+line)
            current = out.setdefault(words[1], [])
        elif words[0] == 'X' and current is not None and len(words) == 13:
            current.append((words[1], words[2], tuple(int(w) for w in words[3:9]),
                            [float(w) for w in words[9:]]))
        else:
            raise ValueError('malformed native line '+line)
    if list(out) != [str(p) for p in paths]:
        raise ValueError('native output does not cover every input in order')
    return out


def capture_free(prefix, output, sdk_manifest):
    """Record the corpus's free shapes natively (S6), certify the
    independent reader's enumeration of them (types and counts) against it,
    before any kernel code imports them."""
    executable, env, _, _ = build(prefix, output, FREE_SOURCE, 'free-oracle')
    corpus = sorted(CORPUS.glob('*.brep'))
    record = run(executable, '\n'.join(str(p) for p in corpus), env)
    if record['exit_code'] != 0:
        raise SystemExit('native free-shape run failed: '+json.dumps(record)[:2000])
    observed = {Path(k).name: v for k, v in parse_free(record['stdout'], corpus).items()}
    for path in corpus:
        ours = [(k, c) for _, k, _, c in brep_io_reference.free_shapes(path.read_text(errors='replace'))]
        theirs = [(k, c) for k, _, c, _ in observed[path.name]]
        if ours != theirs:
            raise ValueError('the independent reader disagrees with OCCT on the free shapes of '+path.name)
    lines = []
    for name, shapes in observed.items():
        lines.append(f'F {name} {len(shapes)}')
        lines += [' '.join(['X', k, v, *map(str, c), *map(repr, m)]) for k, v, c, m in shapes]
    FREE_CAPTURE.mkdir(parents=True, exist_ok=True)
    (FREE_CAPTURE/'native.txt').write_text('\n'.join(lines)+'\n')
    (FREE_CAPTURE/'oracle.cpp').write_text(FREE_SOURCE.read_text())
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    write(FREE_CAPTURE/'capture.json', {
        'source_reference': SOURCE, 'oracle': next(iter(record['stderr'].splitlines()), None),
        'platform': sys.platform, 'rust_revision': revision,
        'rust_free_shape_import_exists': False, 'rust_worktree_uncommitted': status,
        'sdk_manifest_sha256': digest(sdk_manifest),
        'probe_source_sha256': digest(FREE_CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(FREE_CAPTURE/'native.txt'),
        'inputs_sha256': {p.name: digest(p) for p in corpus}})
    return sum(len(v) for v in observed.values())


def free_capture(observed):
    """The S6 native observations of the corpus's free shapes are unchanged:
    types, verdicts and counts exactly, properties within MEASURE_BOUND."""
    metadata = json.loads((FREE_CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_free_shape_import_exists']:
        raise ValueError('free-shape capture is not a clean pre-implementation reference')
    for key, name in [('probe_source_sha256', 'oracle.cpp'), ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(FREE_CAPTURE/name):
            raise ValueError('captured free-shape evidence changed: '+name)
    if (FREE_CAPTURE/'oracle.cpp').read_text() != FREE_SOURCE.read_text():
        raise ValueError('the free-shape oracle differs from the captured one')
    captured, current = {}, None
    for line in (FREE_CAPTURE/'native.txt').read_text().splitlines():
        w = line.split()
        if w[0] == 'F':
            current = captured.setdefault(w[1], [])
        else:
            current.append((w[1], w[2], tuple(int(x) for x in w[3:9]), [float(x) for x in w[9:]]))
    if set(captured) != set(observed):
        raise ValueError('free-shape capture does not cover the corpus')
    for name, shapes in captured.items():
        now = observed[name]
        if [x[:3] for x in now] != [x[:3] for x in shapes] or any(
                measure_apart(a[3], b[3]) for a, b in zip(now, shapes)):
            raise ValueError('native free-shape observation of '+name+' differs from the capture')
    return captured


def v3_files(dataset):
    """The dataset's `.brep` texts whose header is OCCT's version-3 line, as
    TopTools_ShapeSet::Read matches it (a whole line, its trailing CR
    removed), in relative path order."""
    found = []
    for path in sorted(dataset.rglob('*.brep')):
        with path.open('rb') as f:
            head = f.read(512).split(b'\n')[:3]
        if any(line.rstrip(b'\r') == V3_LINE for line in head):
            found.append(path)
    return found


def probe_executable():
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'brep_io_probe'],
                   cwd=ROOT, check=True)
    return str(ROOT/'target/release/examples/brep_io_probe')


def rust_reading(paths):
    """{file name: the kernel's reading} (`brep_io_probe read`): `refused
    ERROR`, or `read VERSION solids OK/ALL free OK/ALL unsupported ...`."""
    rows = subprocess.run([probe_executable(), 'read', '-'], input='\n'.join(map(str, paths))+'\n',
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        name, verdict = line.split(' ', 1)
        out[name] = verdict
    if list(out) != [p.name for p in paths]:
        raise ValueError('the kernel probe does not cover every version-3 file in order')
    return out


def native_lines(observed, row):
    lines = []
    for name, shapes in observed.items():
        lines.append(f'F {name} {len(shapes)}')
        lines += [row(s) for s in shapes]
    return '\n'.join(lines)+'\n'


def capture_v3(prefix, output, sdk_manifest, dataset):
    """Record the native reading of the dataset's version-3 texts by both
    probes, unchanged, while the reader still refuses their header (each
    refusal recorded), before any kernel code reads them."""
    paths = v3_files(dataset)
    names = [p.name for p in paths]
    if not paths or len(set(names)) != len(names):
        raise SystemExit('no version-3 texts, or two with one name, under '+str(dataset))
    reading = rust_reading(paths)
    if any(not v.startswith('refused ') for v in reading.values()):
        raise SystemExit('the kernel reads a version-3 text already: not a pre-implementation capture')
    executable, env, _, _ = build(prefix, output)
    free_executable, free_env, _, _ = build(prefix, output, FREE_SOURCE, 'free-oracle')
    records = [run(executable, '\n'.join(map(str, paths)), env),
               run(free_executable, '\n'.join(map(str, paths)), free_env)]
    if any(r['exit_code'] != 0 for r in records):
        raise SystemExit('native version-3 run failed: '+json.dumps(records)[:2000])
    solids = {Path(k).name: v for k, v in parse(records[0]['stdout'], paths).items()}
    free = {Path(k).name: v for k, v in parse_free(records[1]['stdout'], paths).items()}
    V3_CAPTURE.mkdir(parents=True, exist_ok=True)
    (V3_CAPTURE/'native-solids.txt').write_text(native_lines(
        solids, lambda s: ' '.join(['S', s[0], *map(str, s[1]), *map(repr, s[2])])))
    (V3_CAPTURE/'native-free.txt').write_text(native_lines(
        free, lambda s: ' '.join(['X', s[0], s[1], *map(str, s[2]), *map(repr, s[3])])))
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    write(V3_CAPTURE/'capture.json', {
        'source_reference': SOURCE, 'oracle': next(iter(records[0]['stderr'].splitlines()), None),
        'platform': sys.platform, 'rust_revision': revision,
        'rust_v3_header_read': False, 'rust_refusals': reading, 'rust_worktree_uncommitted': status,
        'sdk_manifest_sha256': digest(sdk_manifest),
        'solid_probe_source_sha256': digest(SOURCE_FILE), 'free_probe_source_sha256': digest(FREE_SOURCE),
        'solid_observations_sha256': digest(V3_CAPTURE/'native-solids.txt'),
        'free_observations_sha256': digest(V3_CAPTURE/'native-free.txt'),
        'inputs_sha256': {str(p.relative_to(dataset)): digest(p) for p in paths}})
    return len(paths), sum(map(len, solids.values())), sum(map(len, free.values()))


def v3_capture():
    """The committed version-3 capture, checked as a clean pre-implementation
    record of the probes as pinned: ({name: solid rows}, {name: free rows},
    metadata)."""
    metadata = json.loads((V3_CAPTURE/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE or metadata['rust_v3_header_read']
            or any(not v.startswith('refused ') for v in metadata['rust_refusals'].values())):
        raise ValueError('version-3 capture is not a clean pre-implementation reference')
    for key, name in [('solid_observations_sha256', 'native-solids.txt'),
                      ('free_observations_sha256', 'native-free.txt')]:
        if metadata[key] != digest(V3_CAPTURE/name):
            raise ValueError('captured version-3 evidence changed: '+name)
    if (metadata['solid_probe_source_sha256'] != digest(SOURCE_FILE)
            or metadata['free_probe_source_sha256'] != digest(FREE_SOURCE)):
        raise ValueError('a probe differs from the one the version-3 capture ran')
    solids, free, current = {}, {}, None
    for line in (V3_CAPTURE/'native-solids.txt').read_text().splitlines():
        w = line.split()
        if w[0] == 'F':
            current = solids.setdefault(w[1], [])
        else:
            current.append((w[1], tuple(int(x) for x in w[2:8]), [float(x) for x in w[8:]]))
    for line in (V3_CAPTURE/'native-free.txt').read_text().splitlines():
        w = line.split()
        if w[0] == 'F':
            current = free.setdefault(w[1], [])
        else:
            current.append((w[1], w[2], tuple(int(x) for x in w[3:9]), [float(x) for x in w[9:]]))
    names = [Path(k).name for k in metadata['inputs_sha256']]
    if list(solids) != names or list(free) != names or list(metadata['rust_refusals']) != names:
        raise ValueError('version-3 capture does not cover its inputs in order')
    return solids, free, metadata


def compare_v3(dataset, output, solid_probe, free_probe, report):
    """The dataset's version-3 texts (HANDOFF's item 3): the native reading
    as captured, the independent reader's counts native's, every file read
    by the kernel, every solid and free shape it imports with native's
    counts (a free shape's certified measure enclosing native's), and its
    written text read natively valid with the same counts and properties."""
    captured_solids, captured_free, metadata = v3_capture()
    paths = v3_files(dataset)
    if [str(p.relative_to(dataset)) for p in paths] != list(metadata['inputs_sha256']):
        raise SystemExit('the version-3 texts differ from the captured ones')
    for p in paths:
        if digest(p) != metadata['inputs_sha256'][str(p.relative_to(dataset))]:
            raise SystemExit('version-3 text changed: '+str(p.relative_to(dataset)))
    records = [run(solid_probe[0], '\n'.join(map(str, paths)), solid_probe[1]),
               run(free_probe[0], '\n'.join(map(str, paths)), free_probe[1])]
    if any(r['exit_code'] != 0 for r in records):
        raise SystemExit('native version-3 run failed: '+json.dumps(records)[:2000])
    solids = {Path(k).name: v for k, v in parse(records[0]['stdout'], paths).items()}
    free = {Path(k).name: v for k, v in parse_free(records[1]['stdout'], paths).items()}
    for name in solids:
        if ([s[:2] for s in solids[name]] != [s[:2] for s in captured_solids[name]]
                or any(close(a[2], b[2]) > PROPERTY_BOUND for a, b in zip(solids[name], captured_solids[name]))
                or [s[:3] for s in free[name]] != [s[:3] for s in captured_free[name]]
                or any(measure_apart(a[3], b[3]) for a, b in zip(free[name], captured_free[name]))):
            raise ValueError('native observation of version-3 '+name+' differs from the capture')
    report.update(v3_files=len(paths), v3_reading={}, v3_solids_imported=0, v3_free_imported=0,
                  v3_matches=[], v3_unwritable=[], v3_not_measured=[])
    fail = lambda case, reason, **more: report['failures'].append(dict(case='v3 '+case, reason=reason, **more))
    # The independent reader against native OCCT: records in the root's order.
    native_solid, native_free = {}, {}
    for p in paths:
        text = p.read_text()
        _, rows = brep_io_reference.summary(text)
        free_rows = brep_io_reference.free_shapes(text)
        if [c for _, _, c in rows] != [c for _, c, _ in solids[p.name]]:
            fail(p.name, 'reader solid counts differ from native')
        if [(k, c) for _, k, _, c in free_rows] != [(k, c) for k, _, c, _ in free[p.name]]:
            fail(p.name, 'reader free shapes differ from native')
        native_solid.update({(p.name, r): (c, n[2]) for (r, _, c), n in zip(rows, solids[p.name])})
        native_free.update({(p.name, r): (k, c, n[3]) for (r, k, _, c), n in zip(free_rows, free[p.name])})
    # The kernel reads every file.
    for name, verdict in rust_reading(paths).items():
        report['v3_reading'][name] = verdict
        if not verdict.startswith('read V3 '):
            fail(name, 'the kernel does not read it', rust=verdict)
    written = output/'v3-written'
    written.mkdir(parents=True, exist_ok=True)
    for f in written.glob('*.brep'):
        f.unlink()
    probe = probe_executable()
    listing = '\n'.join(map(str, paths))+'\n'
    targets = []
    for line in subprocess.run([probe, 'corpus', str(written)], input=listing, text=True,
                               capture_output=True, timeout=600, check=True).stdout.splitlines():
        w = line.split()
        key, counts = (w[0], int(w[1])), tuple(int(x) for x in w[2:8])
        report['v3_solids_imported'] += 1
        if key not in native_solid or native_solid[key][0] != counts:
            fail(f'{key[0]}-{key[1]}', 'imported solid counts differ from native', rust=counts,
                 native=native_solid.get(key, [None])[0])
            continue
        if w[8:] == ['unwritable']:
            report['v3_unwritable'].append(f'{key[0]}-{key[1]}')
            continue
        targets.append(('solid', f'{key[0]}-{key[1]}', written/f'{key[0]}-{key[1]}.brep', native_solid[key]))
    for line in subprocess.run([probe, 'free', str(written)], input=listing, text=True,
                               capture_output=True, timeout=600, check=True).stdout.splitlines():
        w = line.split()
        unwritable = w[-1] == 'unwritable'
        w = w[:-1] if unwritable else w
        case = f'{w[0]}-free-{w[1]}'
        key, kind, counts = (w[0], int(w[1])), w[2], tuple(int(x) for x in w[3:9])
        enclosure = None if w[9:] == ['-'] else list(zip(map(float, w[9::2]), map(float, w[10::2])))
        report['v3_free_imported'] += 1
        native = native_free.get(key)
        if native is None or native[:2] != (kind, counts):
            fail(case, 'imported free shape differs from native', rust=[kind, counts], native=native and native[:2])
            continue
        # A shape the kernel cannot measure or write is listed: neither is
        # the header's.
        if enclosure is None:
            report['v3_not_measured'].append(case)
        elif measure_outside(enclosure, native[2]):
            fail(case, 'measure outside the enclosure: '+','.join(measure_outside(enclosure, native[2])),
                 rust=enclosure, native=native[2])
            continue
        if unwritable:
            report['v3_unwritable'].append(case)
            continue
        targets.append(('free', case, written/f'{case}.brep', (counts, native[2])))
    # Each written text read natively: one valid shape of the original's
    # counts and properties.
    for kind, probe_pair, rows_of in [('solid', solid_probe, lambda out, ps: parse(out, ps)),
                                      ('free', free_probe, lambda out, ps: parse_free(out, ps))]:
        chosen = [t for t in targets if t[0] == kind]
        if not chosen:
            continue
        ps = [t[2] for t in chosen]
        record = run(probe_pair[0], '\n'.join(map(str, ps)), probe_pair[1])
        if record['exit_code'] != 0:
            raise SystemExit('native run on written version-3 shapes failed: '+json.dumps(record)[:2000])
        rows = rows_of(record['stdout'], ps)
        for _, case, path, (counts, props) in chosen:
            got = rows[str(path)]
            differences = []
            if len(got) != 1:
                differences.append('shape_count')
            else:
                verdict, got_counts, got_props = (got[0] if kind == 'solid' else got[0][1:])
                differences += ['brepcheck_invalid'] if verdict != 'valid' else []
                differences += ['counts'] if got_counts != counts else []
                apart = (close(got_props, props) > PROPERTY_BOUND if kind == 'solid'
                         else measure_apart(got_props, props))
                differences += ['properties'] if apart else []
            if differences:
                fail(case, 'written text differs natively: '+','.join(differences), native=got)
            else:
                report['v3_matches'].append(case)


def measure_apart(a, b):
    """Whether two [area or length, cx, cy, cz] rows differ by more than
    MEASURE_BOUND of the larger row's largest magnitude."""
    scale = max([1.0]+[abs(x) for x in a+b])
    return any(abs(x-y) > MEASURE_BOUND*scale for x, y in zip(a, b))


def measure_outside(enclosure, native):
    """The labels of a certified [lo, hi] measure row that do not contain
    OCCT's value, up to MEASURE_BOUND of its largest magnitude (OCCT's own
    integration error)."""
    allowance = MEASURE_BOUND*max([1.0]+[abs(x) for x in native])
    return [label for label, (lo, hi), x in zip(['measure', 'cx', 'cy', 'cz'], enclosure, native)
            if not lo-allowance <= x <= hi+allowance]


def measure_differences(enclosure, default, adaptive):
    """The labels of a free face whose enclosure excludes OCCT's default
    measure: `default_integration` when the adaptive measure lies within it
    (up to MEASURE_BOUND), else what the adaptive measure leaves outside."""
    outside = measure_outside(enclosure, adaptive)
    return ['default_integration'] if not outside else ['adaptive_outside:'+','.join(outside)]


def independent_outside(enclosure, values):
    """The labels of a certified [lo, hi] row that do not contain a review's
    independent measure (decimal strings, exactly), up to binary64 rounding
    of the row's largest magnitude: the decimals round a high-precision
    value, and a plane's centre coordinate is enclosed with zero width."""
    slack = Fraction(2.0**-52)*max([1]+[abs(Fraction(x)) for x in values])
    return [label for label, (lo, hi), x in zip(['measure', 'cx', 'cy', 'cz'], enclosure, values)
            if not Fraction(lo)-slack <= Fraction(x) <= Fraction(hi)+slack]


def measure_fingerprint(default, adaptive):
    """Both native rows to ten significant digits: the reviewed
    discrepancies are 1e-9 or more, and a platform's last-bit drift (the
    S6 capture allows MEASURE_BOUND) keeps the fingerprint."""
    return sha(json.dumps([[f'{x:.10g}' for x in row] for row in (default, adaptive)]))


def review_measure(case, text, oracle, enclosure, default, adaptive, reviews):
    """(review, None) for a reviewed free face whose recorded independent
    measure the enclosure contains, else (None, failure)."""
    evidence = {'case': case, 'source_reference': SOURCE, 'oracle': oracle, 'input_sha256': sha(text),
                'native_sha256': measure_fingerprint(default, adaptive),
                'differences': measure_differences(enclosure, default, adaptive)}
    review = review_for(evidence, reviews)
    independent = review.get('independent_measure') if review else None
    if review and independent and len(independent) == 4 and not independent_outside(enclosure, independent):
        return review, None
    reason = ('unreviewed measure difference' if not review else
              'the review records no independent measure' if not independent or len(independent) != 4 else
              'independent measure outside the enclosure: '+','.join(independent_outside(enclosure, independent)))
    return None, dict(evidence, reason=reason, rust=enclosure, native=default, adaptive=adaptive)


def close(a, b):
    """Relative agreement of [volume, area, cx, cy, cz]."""
    scale = max(1.0, abs(b[0]))**(1/3)
    return max([abs(x-y)/max(1.0, abs(y)) for x, y in zip(a[:2], b[:2])]
               + [abs(x-y)/scale for x, y in zip(a[2:], b[2:])])


def rust_outputs(output):
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'brep_io_probe'],
                   cwd=ROOT, check=True)
    probe = str(ROOT/'target/release/examples/brep_io_probe')
    prisms, written = output/'prisms', output/'written'
    for d in (prisms, written):
        d.mkdir(parents=True, exist_ok=True)
        for f in d.glob('*.brep'):
            f.unlink()
    # The identity prisms and, since S5, the arc prisms.
    cases = ((ROOT/'rust/fixtures/identity-cases.txt').read_text()
             + (ROOT/'rust/fixtures/identity-arc-cases.txt').read_text())
    prism_rows = subprocess.run([probe, 'prisms', str(prisms)], input=cases, text=True,
                                capture_output=True, timeout=600, check=True).stdout
    corpus = '\n'.join(str(p) for p in sorted(CORPUS.glob('*.brep')))+'\n'
    corpus_rows = subprocess.run([probe, 'corpus', str(written)], input=corpus, text=True,
                                 capture_output=True, timeout=600, check=True).stdout
    prism = {}
    for line in prism_rows.splitlines():
        w = line.split()
        prism[w[0]] = (tuple(int(x) for x in w[1:7]), [float(x) for x in w[7:]])
    imported = {}
    for line in corpus_rows.splitlines():
        w = line.split()
        imported[(w[0], int(w[1]))] = tuple(int(x) for x in w[2:8])
    # S6: free shapes and face and wire bodies, with certified measures.
    measure = lambda w: None if w == ['-'] else list(zip(map(float, w[::2]), map(float, w[1::2])))
    free_rows = subprocess.run([probe, 'free', str(written)], input=corpus, text=True,
                               capture_output=True, timeout=600, check=True).stdout
    free = {}
    for line in free_rows.splitlines():
        w = line.split()
        free[(w[0], int(w[1]))] = (w[2], tuple(int(x) for x in w[3:9]), measure(w[9:]))
    body_rows = subprocess.run([probe, 'bodies', str(prisms)],
                               input=(ROOT/'rust/fixtures/identity-sheet-cases.txt').read_text(), text=True,
                               capture_output=True, timeout=600, check=True).stdout
    bodies = {}
    for line in body_rows.splitlines():
        w = line.split()
        bodies[w[0]] = (tuple(int(x) for x in w[1:7]), measure(w[7:]))
    return prism, prisms, imported, written, free, bodies


def reference_rows():
    """The committed independent fixture, regenerated: {file: [(record, verdict, counts)]},
    and (S6) {file: [(record, kind, verdict, counts)]} of its free shapes."""
    text = generate_brep_io_fixtures.rows()
    if text != (ROOT/'rust/fixtures/brep-io-expected.tsv').read_text():
        raise ValueError('independent fixture regeneration changed: brep-io-expected.tsv')
    files, free = {}, {}
    for line in text.splitlines():
        w = line.split('\t')
        if w[0] == 'file':
            files[w[1]], free[w[1]] = [], []
        elif w[0] == 'free':
            free[w[1]].append((int(w[2]), w[3], w[4], tuple(int(x) for x in w[5:11])))
        else:
            files[w[1]].append((int(w[2]), w[3], tuple(int(x) for x in w[4:10])))
    return files, free


def original_capture(observed):
    """The native observations of the unmodified corpus are unchanged."""
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE:
        raise ValueError('capture is not of the pinned source')
    for key, name in [('probe_source_sha256', 'oracle.cpp'), ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('captured native evidence changed: '+name)
    if (CAPTURE/'oracle.cpp').read_text() != SOURCE_FILE.read_text():
        raise ValueError('the oracle differs from the captured one')
    for path in sorted(CORPUS.glob('*.brep')):
        if metadata['inputs_sha256'][path.name] != digest(path):
            raise ValueError('corpus file changed: '+path.name)
    captured = {}
    for line in (CAPTURE/'native.txt').read_text().splitlines():
        w = line.split()
        if w[0] == 'F':
            current = captured.setdefault(w[1], [])
        else:
            current.append((w[1], tuple(int(x) for x in w[2:8]), [float(x) for x in w[8:]]))
    if set(captured) != set(observed):
        raise ValueError('capture does not cover the corpus')
    for name, solids in captured.items():
        now = observed[name]
        if [s[:2] for s in now] != [s[:2] for s in solids] or any(
                close(a[2], b[2]) > PROPERTY_BOUND for a, b in zip(now, solids)):
            raise ValueError('native observation of '+name+' differs from the capture')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/brep-io-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', action='store_true',
                        help='write the corpus observations to the capture directory instead of checking them')
    parser.add_argument('--capture-free', action='store_true',
                        help='record the corpus free shapes (S6, before implementation only)')
    parser.add_argument('--v3-dataset', type=Path,
                        help="the public dataset's directory: compare its version-3 texts (local only, U1)")
    parser.add_argument('--capture-v3', action='store_true',
                        help='record the version-3 texts natively (with --v3-dataset, before implementation only)')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    if args.capture_free:
        print('captured', capture_free(prefix, output, args.sdk_manifest), 'free shapes')
        return
    if args.capture_v3:
        if not args.v3_dataset:
            parser.error('--capture-v3 needs --v3-dataset')
        files, solids, free = capture_v3(prefix, output, args.sdk_manifest, args.v3_dataset.resolve())
        print('captured', files, 'version-3 texts:', solids, 'solids,', free, 'free shapes')
        return
    reference, reference_free = reference_rows()
    executable, env, loaded, command = build(prefix, output)
    free_executable, free_env, _, _ = build(prefix, output, FREE_SOURCE, 'free-oracle')
    corpus = sorted(CORPUS.glob('*.brep'))
    record = run(executable, '\n'.join(str(p) for p in corpus), env)
    if record['exit_code'] != 0:
        raise SystemExit('native corpus run failed: '+json.dumps(record)[:2000])
    oracle = next(iter(record['stderr'].splitlines()), None)
    observed = {Path(k).name: v for k, v in parse(record['stdout'], corpus).items()}
    if args.capture:
        CAPTURE.mkdir(parents=True, exist_ok=True)
        lines = []
        for name, solids in observed.items():
            lines.append(f'F {name} {len(solids)}')
            lines += [' '.join(['S', v, *map(str, c), *map(repr, p)]) for v, c, p in solids]
        (CAPTURE/'native.txt').write_text('\n'.join(lines)+'\n')
        (CAPTURE/'oracle.cpp').write_text(SOURCE_FILE.read_text())
        write(CAPTURE/'capture.json', {
            'source_reference': SOURCE, 'oracle': oracle, 'platform': sys.platform,
            'sdk_manifest_sha256': digest(args.sdk_manifest),
            'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
            'observations_sha256': digest(CAPTURE/'native.txt'),
            'inputs_sha256': {p.name: digest(p) for p in corpus}})
    original_capture(observed)
    corpus_free = run(free_executable, '\n'.join(str(p) for p in corpus), free_env)
    if corpus_free['exit_code'] != 0:
        raise SystemExit('native free-shape run failed: '+json.dumps(corpus_free)[:2000])
    observed_free = {Path(k).name: v for k, v in parse_free(corpus_free['stdout'], corpus).items()}
    free_capture(observed_free)
    prism, prisms, imported, written, free, bodies = rust_outputs(output)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'corpus_files': len(corpus),
              'reader_solids_certified': 0, 'prisms_verified': 0, 'written_solids_verified': 0,
              'reader_free_certified': 0, 'free_imported': 0, 'free_rejected_by_validator': [],
              'free_measures_contained': 0, 'free_measures_reviewed': 0, 'written_free_verified': 0,
              'bodies_verified': 0, 'worst_property_difference': 0.0, 'matches': [],
              'reviewed_differences': [], 'failures': []}

    # The independent reader against native OCCT, file by file.
    for name, solids in reference.items():
        native = observed[name]
        if [c for _, _, c in solids] != [c for _, c, _ in native]:
            report['failures'].append({'case': name, 'reason': 'reader counts differ from native',
                                       'reader': [c for _, _, c in solids],
                                       'native': [c for _, c, _ in native]})
        else:
            report['reader_solids_certified'] += len(solids)
    representable = {(name, r): (c, native[2]) for name, solids in reference.items()
                     for (r, verdict, c), native in zip(solids, observed[name]) if verdict == 'representable'}
    # Representable structure the kernel's validator does not certify yet
    # (S4), pinned as in occt_brep.rs.
    uncertified = {('Ball.brep', 108)}
    if set(representable) - uncertified != set(imported):
        report['failures'].append({'case': 'corpus', 'reason': 'imported solids differ from the representable ones',
                                   'imported': sorted(map(list, imported)),
                                   'representable': sorted(map(list, representable))})

    # S6: the reader's free shapes against native OCCT, then the kernel's
    # imports against both.
    representable_free = {}
    for name, shapes in reference_free.items():
        native = observed_free[name]
        if [(k, c) for _, k, _, c in shapes] != [(k, c) for k, _, c, _ in native]:
            report['failures'].append({'case': name, 'reason': 'reader free shapes differ from native'})
            continue
        report['reader_free_certified'] += len(shapes)
        for (r, k, verdict, c), (_, _, _, props) in zip(shapes, native):
            if verdict == 'representable':
                representable_free[(name, r)] = (k, c, props)
    if not set(free) <= set(representable_free):
        report['failures'].append({'case': 'corpus', 'reason': 'imported free shapes the reader does not represent',
                                   'extra': sorted(map(list, set(free)-set(representable_free)))})
    # Representable but rejected by the validator: pinned with their issue
    # kinds in occt_brep.rs (the hammer's faces).
    report['free_rejected_by_validator'] = sorted(f'{n}-{r}' for n, r in set(representable_free)-set(free))
    measure_cases = []
    for key, (kind, counts, enclosure) in sorted(free.items()):
        if key not in representable_free:
            continue
        k, c, props = representable_free[key]
        found = [] if (kind, counts) == (k, c) else ['kind_or_counts']
        found += ['not measured'] if enclosure is None else []
        report['free_imported'] += 1
        if found:
            report['failures'].append({'case': f'{key[0]}-free-{key[1]}', 'reason': ' '.join(found),
                                       'rust': [kind, counts, enclosure], 'native': [k, c, props]})
        elif measure_outside(enclosure, props):
            measure_cases.append((key, enclosure, props))
        else:
            report['free_measures_contained'] += 1
    # Enclosures that exclude OCCT's default measure: the adaptive probe on
    # those files, and a review with an independent measure for each.
    adaptive_oracle = None
    if measure_cases:
        adaptive_executable, adaptive_env, _, _ = build(prefix, output, ADAPTIVE_SOURCE, 'adaptive-oracle')
        files = sorted({n for (n, _), _, _ in measure_cases})
        record = run(adaptive_executable, '\n'.join(str(CORPUS/n) for n in files), adaptive_env)
        write(output/'native-adaptive.json', record)
        if record['exit_code'] != 0:
            raise SystemExit('native adaptive run failed: '+json.dumps(record)[:2000])
        adaptive_oracle = next(iter(record['stderr'].splitlines()), None)
        rows = parse_free(record['stdout'], [CORPUS/n for n in files])
        adaptive = {}
        for n in files:
            shapes = reference_free[n]
            native = rows[str(CORPUS/n)]
            if [(k, c) for _, k, _, c in shapes] != [(k, c) for k, _, c, _ in native]:
                raise SystemExit('the adaptive probe disagrees with the reader on the free shapes of '+n)
            adaptive.update({(n, r): props for (r, _, _, _), (_, _, _, props) in zip(shapes, native)})
        for key, enclosure, props in measure_cases:
            review, failure = review_measure(f'{key[0]}-free-{key[1]}', (CORPUS/key[0]).read_text(),
                                             adaptive_oracle, enclosure, props, adaptive[key], reviews)
            report['reviewed_differences' if review else 'failures'].append(review or failure)
            report['free_measures_reviewed'] += review is not None
    free_targets = [(f'written {n}-free-{r}', written/f'{n}-free-{r}.brep', representable_free[(n, r)][1],
                     representable_free[(n, r)][2], None)
                    for (n, r) in free if (n, r) in representable_free]
    free_targets += [(f'body {name}', prisms/f'{name}.brep', counts, None, enclosure)
                     for name, (counts, enclosure) in bodies.items()]
    paths = [p for _, p, _, _, _ in free_targets]
    written_free = run(free_executable, '\n'.join(map(str, paths)), free_env)
    if written_free['exit_code'] != 0:
        raise SystemExit('native run on written free shapes failed: '+json.dumps(written_free)[:2000])
    native_free = parse_free(written_free['stdout'], paths)
    for case, path, counts, props, enclosure in free_targets:
        shapes = native_free[str(path)]
        differences = []
        if len(shapes) != 1:
            differences.append('shape_count')
        else:
            _, verdict, got, got_props = shapes[0]
            if verdict != 'valid':
                differences.append('brepcheck_invalid')
            if got != counts:
                differences.append('counts')
            if props is not None and measure_apart(got_props, props):
                differences.append('properties')
            if enclosure is not None and measure_outside(enclosure, got_props):
                differences.append('measure_outside')
        if not differences:
            report['matches'].append(case)
            report['bodies_verified' if case.startswith('body') else 'written_free_verified'] += 1
            continue
        text = path.read_text()
        evidence = {'case': case, 'source_reference': SOURCE, 'oracle': oracle, 'input_sha256': sha(text),
                    'native_stdout_sha256': sha(json.dumps(shapes)), 'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, native=shapes, expected=[counts, props, enclosure]))

    # Written files: prisms against exact properties, corpus solids against
    # their originals.
    targets = [(f'prism {name}', prisms/f'{name}.brep', counts, props) for name, (counts, props) in prism.items()]
    targets += [(f'written {n}-{r}', written/f'{n}-{r}.brep', counts, representable[(n, r)][1])
                for (n, r), counts in imported.items() if (n, r) in representable]
    for (n, r), counts in imported.items():
        if (n, r) in representable and counts != representable[(n, r)][0]:
            report['failures'].append({'case': f'{n}-{r}', 'reason': 'imported counts differ from the original'})
    paths = [p for _, p, _, _ in targets]
    record = run(executable, '\n'.join(map(str, paths)), env)
    write(output/'native-written.json', record)
    if record['exit_code'] != 0:
        raise SystemExit('native run on written files failed: '+json.dumps(record)[:2000])
    native = parse(record['stdout'], paths)
    for case, path, counts, props in targets:
        solids = native[str(path)]
        differences = []
        if len(solids) != 1:
            differences.append('solid_count')
        else:
            verdict, got, got_props = solids[0]
            if verdict != 'valid':
                differences.append('brepcheck_invalid')
            if got != counts:
                differences.append('counts')
            d = close(got_props, props)
            report['worst_property_difference'] = max(report['worst_property_difference'], d)
            if not d <= PROPERTY_BOUND:
                differences.append('mass_properties')
        if not differences:
            report['matches'].append(case)
            report['prisms_verified' if case.startswith('prism') else 'written_solids_verified'] += 1
            continue
        text = path.read_text()
        evidence = {'case': case, 'source_reference': SOURCE, 'oracle': oracle, 'input_sha256': sha(text),
                    'native_stdout_sha256': sha(json.dumps(solids)), 'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, native=solids, expected=[counts, props]))
    if args.v3_dataset:
        compare_v3(args.v3_dataset.resolve(), output, (executable, env), (free_executable, free_env), report)
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                'adaptive_oracle': adaptive_oracle, 'adaptive_source_sha256': digest(ADAPTIVE_SOURCE),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) and k != 'free_rejected_by_validator' else v
                      for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
