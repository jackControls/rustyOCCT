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

The native observations of the unmodified upstream files were captured
after the kernel's reader existed (rust/fixtures/occt-brep-io-capture/NOTES.md);
they depend on no Rust output and must stay unchanged. Differences need a
fingerprinted review; timeouts, crashes and malformed output cannot be
reviewed.
"""
import argparse
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
# Both native runs together took under a second; the deadline only bounds hangs.
TIMEOUT = 600
# The worst observed difference was 2.5e-14: volume and area relative to
# themselves, centroids relative to the cube root of the volume.
PROPERTY_BOUND = 1e-11
# S6: OCCT integrates spline faces with a relative error near 1e-9
# (occt-sheet-preimplementation/NOTES.md); the kernel's enclosures must
# contain its areas, lengths and centres up to that.
MEASURE_BOUND = 1e-9


def build(prefix, output, source=SOURCE_FILE, name='oracle', required=('TKTopAlgo', 'TKBRep')):
    include, lib = prefix/'include/opencascade', prefix/'lib'
    executable = output/name
    command = shlex.split(os.environ.get('CXX', 'c++'))+[
        '-std=c++17', '-O2', str(source), '-I'+str(include), '-L'+str(lib),
        '-Wl,-rpath,'+str(lib), '-o', str(executable)]+['-l'+name for name in TOOLKITS]
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
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    if args.capture_free:
        print('captured', capture_free(prefix, output, args.sdk_manifest), 'free shapes')
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
              'free_measures_contained': 0, 'written_free_verified': 0, 'bodies_verified': 0,
              'worst_property_difference': 0.0, 'matches': [], 'reviewed_differences': [], 'failures': []}

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
    uncertified = {('Ball.brep', 108), ('Motor-c.brep', 378)}
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
    for key, (kind, counts, enclosure) in sorted(free.items()):
        if key not in representable_free:
            continue
        k, c, props = representable_free[key]
        found = [] if (kind, counts) == (k, c) else ['kind_or_counts']
        found += ['not measured'] if enclosure is None else measure_outside(enclosure, props)
        report['free_imported'] += 1
        if found:
            report['failures'].append({'case': f'{key[0]}-free-{key[1]}', 'reason': ' '.join(found),
                                       'rust': [kind, counts, enclosure], 'native': [k, c, props]})
        else:
            report['free_measures_contained'] += 1
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
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) and k != 'free_rejected_by_validator' else v
                      for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
