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
  original's volume, surface area and centroid.

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
import generate_brep_io_fixtures

CAPTURE = ROOT/'rust/fixtures/occt-brep-io-capture'
REVIEWS = ROOT/'rust/fixtures/occt-brep-io-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_brep_io_oracle.cpp'
CORPUS = ROOT/'data/occ'
TOOLKITS = ['TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath', 'TKernel']
# Both native runs together took under a second; the deadline only bounds hangs.
TIMEOUT = 600
# The worst observed difference was 2.5e-14: volume and area relative to
# themselves, centroids relative to the cube root of the volume.
PROPERTY_BOUND = 1e-11


def build(prefix, output):
    include, lib = prefix/'include/opencascade', prefix/'lib'
    executable = output/'oracle'
    command = shlex.split(os.environ.get('CXX', 'c++'))+[
        '-std=c++17', '-O2', str(SOURCE_FILE), '-I'+str(include), '-L'+str(lib),
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
    for toolkit in ['TKTopAlgo', 'TKBRep']:
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
    cases = (ROOT/'rust/fixtures/identity-cases.txt').read_text()
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
    return prism, prisms, imported, written


def reference_rows():
    """The committed independent fixture, regenerated: {file: [(record, verdict, counts)]}."""
    text = generate_brep_io_fixtures.rows()
    if text != (ROOT/'rust/fixtures/brep-io-expected.tsv').read_text():
        raise ValueError('independent fixture regeneration changed: brep-io-expected.tsv')
    files = {}
    for line in text.splitlines():
        w = line.split('\t')
        if w[0] == 'file':
            files[w[1]] = []
        else:
            files[w[1]].append((int(w[2]), w[3], tuple(int(x) for x in w[4:10])))
    return files


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
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    reference = reference_rows()
    executable, env, loaded, command = build(prefix, output)
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
    prism, prisms, imported, written = rust_outputs(output)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'corpus_files': len(corpus),
              'reader_solids_certified': 0, 'prisms_verified': 0, 'written_solids_verified': 0,
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
    if set(representable) != set(imported):
        report['failures'].append({'case': 'corpus', 'reason': 'imported solids differ from the representable ones',
                                   'imported': sorted(map(list, imported)),
                                   'representable': sorted(map(list, representable))})

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
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
