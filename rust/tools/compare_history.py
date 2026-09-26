#!/usr/bin/env python3
"""Source-pinned OCCT prism/transform history beside the Rust histories.

Rust histories are first certified equal to the independent enumeration
(fixture regeneration plus the history tests). For each case captured before
implementation, the native probe's BRepTools_History and MakePrism queries
are mapped to Rust relations by documented role correspondences:

    profile vertex: Generated -> Vertical, FirstShape -> BottomVertex,
                    LastShape  -> TopVertex
    profile edge:   Generated -> Wall, FirstShape -> BottomEdge,
                    LastShape -> TopEdge
    profile face:   FirstShape -> StartCap, LastShape -> EndCap,
                    Generated -> the solid region (signed as the body)

OCCT's seams and the vertices only they use (S rows: edges closed on a face,
vertices used only by those and by edges closed on them) are structure the
seamless cell model does not have. Queries, outputs and transform pairs that
reach them are classified structure-only by that rule, and the rule is
verified on every case by count synthesis: the kernel's synthesized OCCT
counts of the Rust body equal the native N row.

A case matches when every other native query has exactly one Rust relation
with an equal geometric signature, every Rust relation is reached, every
native output is covered, entity counts agree, and each transform step's
before/after pairs (the region's against the native bodies) correspond one to
one. Differences need a fingerprinted review; timeouts, crashes and malformed
output are never reviewable.
"""
import argparse
import functools
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from build_pinned_occt import SOURCE, digest
from compare_brep import same_inputs
from compare_occt import ROOT
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
import generate_history_fixtures
import generate_identity_fixtures
from identity_reference import encode_case, native_case

ORIGINAL = ROOT/'rust/fixtures/occt-history-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-history-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_history_oracle.cpp'
TOOLKITS = ['TKPrim', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d', 'TKG2d', 'TKMath', 'TKernel']
# Every preflight case finished within 0.03 s; the deadline only bounds hangs.
TIMEOUT = 120
RELATIVE = 1e-9

NATIVE_ROLE = {
    ('vertex', 'gen'): {('vertical', None)},
    ('vertex', 'first'): {('bottom_vertex', None)},
    ('vertex', 'last'): {('top_vertex', None)},
    ('edge', 'gen'): {('wall', None)},
    ('edge', 'first'): {('bottom_edge', None)},
    ('edge', 'last'): {('top_edge', None)},
    ('face', 'first'): {('start_cap', None)},
    ('face', 'last'): {('end_cap', None)},
}


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def review_for(evidence, reviews):
    return next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                 and all(r.get(k) == v for k, v in evidence.items())), None)


@functools.lru_cache(maxsize=None)
def chosen():
    """The captured cases, from the independent generator, in capture order."""
    names = json.loads((ORIGINAL/'capture.json').read_text())['cases']
    cases, _ = generate_identity_fixtures.generate()
    by_name = {c.name: c for c in cases}
    return tuple(by_name[n] for n in names)


def original_capture():
    """The pre-implementation observations and their inputs are unchanged."""
    metadata = json.loads((ORIGINAL/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE
            or metadata['rust_identity_or_history_implementation_exists']
            or any(line[3:].startswith('rust/kernel') for line in metadata['rust_worktree_uncommitted'])):
        raise ValueError('original capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.json')]:
        if metadata[key] != digest(ORIGINAL/name):
            raise ValueError('original native evidence changed: '+name)
    same_inputs((ORIGINAL/'inputs.txt').read_text(), '\n'.join(native_case(c) for c in chosen())+'\n')


# ------------------------------------------------------------------ parsing

def parse_native(stdout, name):
    """Queries, output coverage, body rows and transform pairs."""
    lines = stdout.splitlines()
    if not lines or not lines[0].startswith(f'R {name} ') or lines[-1] != 'end':
        raise ValueError(f'malformed native output for {name}')
    out = {'valid': lines[0].split()[2] == '1', 'queries': [], 'outputs': [], 'bodies': [],
           'transforms': {}, 'counts': None, 'structure': {}}
    i = 1
    while i < len(lines)-1:
        w = lines[i].split()
        if w[0] == 'Q':
            if w[4] == 'deleted':
                out['queries'].append((' '.join(w[1:4]), 'deleted', None))
            else:
                out['queries'].append((' '.join(w[1:4]), w[4], ' '.join(w[5:])))
        elif w[0] == 'O':
            out['outputs'].append((' '.join(w[1:-1]), w[-1]))
        elif w[0] == 'B':
            out['bodies'].append((' '.join(w[1:-1]), w[-1] == '1'))
        elif w[0] == 'N':
            out['counts'] = ' '.join(w[1:])
        elif w[0] == 'S':
            out['structure'].setdefault(w[1], set()).add(' '.join(w[2:]))
        elif w[0] == 'T':
            k, arrow = int(w[1]), w.index('->')
            before, count = ' '.join(w[2:arrow]), int(w[arrow+1])
            images = [lines[i+1+j].strip() for j in range(count)]
            out['transforms'].setdefault(k, []).append((before, images))
            i += count
        else:
            raise ValueError(f'unknown native row {lines[i]!r}')
        i += 1
    return out


def parse_rust(text):
    cases, current = {}, None
    for line in text.splitlines():
        w = line.split()
        if w[0] == 'R':
            current = cases[w[1]] = {'relations': [], 'bodies': [], 'transforms': {}, 'counts': None}
        elif w[0] == 'C':
            current['counts'] = ' '.join(w[1:])
        elif w[0] == 'G':
            bar = w.index('|')
            current['relations'].append((w[1], int(w[2]), ' '.join(w[3:bar]), ' '.join(w[bar+1:])))
        elif w[0] == 'B':
            current['bodies'].append(' '.join(w[1:]))
        elif w[0] == 'T':
            arrow = w.index('->')
            current['transforms'].setdefault(int(w[1]), []).append((' '.join(w[2:arrow]), ' '.join(w[arrow+1:])))
    return cases


# ------------------------------------------------------------------ comparison

def numbers(sig):
    return [float(x) for x in sig.split()[1:] if x not in ('line', 'circle', 'plane', 'cylinder', 'other', 'S')]


def kind(sig):
    w = sig.split()
    return (w[0], w[1]) if w[0] in ('E', 'F') else (w[0],)


def same(a, b, scale):
    """Equal geometric signatures: kinds match, numbers agree within the
    relative budget; an edge may run either way."""
    if kind(a) != kind(b):
        return False
    x, y = numbers(a), numbers(b)
    if len(x) != len(y):
        return False
    close = lambda p, q: all(abs(s-t) <= RELATIVE*max(scale, abs(s), abs(t)) for s, t in zip(p, q))
    if a.startswith('E '):
        return close(x, y) or close(x[6:9]+x[3:6]+x[0:3], y)
    return close(x, y)


def scale_of(rust):
    values = [abs(v) for _, _, _, sig in rust['relations'] for v in numbers(sig)]
    return max(values+[1.0])


def compare(native, rust):
    """Sorted difference classes; empty means the case matches."""
    diffs = set()
    scale = scale_of(rust)
    if not native['valid'] or not all(ok for _, ok in native['bodies']):
        diffs.add('native_invalid')
    if any(state != 'direct' for _, state in native['outputs']):
        diffs.add('native_history_incomplete')
    structure = native['structure'].get('-', set())
    if native['counts'] is None or native['counts'] != rust['counts']:
        diffs.add('synthesized_counts')
    reached = set()
    for label, query, sig in native['queries']:
        if query == 'deleted':
            diffs.add('native_deleted')
            continue
        if sig in structure:
            continue
        element = label.split()[1] if not label.startswith('face') else 'face'
        if (element, query) == ('face', 'gen'):
            regions = [k for k, (role, _, loc, _) in enumerate(rust['relations'])
                       if role == 'region' and loc == label]
            if len(regions) != 1 or not same(sig, rust['relations'][regions[0]][3], scale):
                diffs.add('region_signature')
            else:
                reached.add(regions[0])
            if not same(sig, rust['bodies'][0], scale):
                diffs.add('body_signature')
            continue
        wanted = NATIVE_ROLE.get((element, query))
        if wanted is None:
            diffs.add(f'unmapped_{element}_{query}')
            continue
        matches = [k for k, (role, ordinal, loc, s) in enumerate(rust['relations'])
                   if loc == label and ((role, None) in wanted or (role, ordinal) in wanted)]
        if len(matches) != 1:
            diffs.add(f'missing_rust_{element}_{query}')
            continue
        reached.add(matches[0])
        if not same(sig, rust['relations'][matches[0]][3], scale):
            diffs.add(f'signature_{element}_{query}')
    if len(reached) != len(rust['relations']):
        diffs.add('unreached_rust_relation')
    for t in ('V', 'E', 'F'):
        count = sum(1 for sig, _ in native['outputs'] if sig.startswith(t+' ') and sig not in structure)
        if count != sum(1 for *_, s in rust['relations'] if s.startswith(t+' ')):
            diffs.add(f'count_{t}')
    if sorted(native['transforms']) != sorted(rust['transforms']):
        diffs.add('transform_steps')
    for k, pairs in native['transforms'].items():
        skipped = native['structure'].get(str(k), set())
        pairs = [(before, images) for before, images in pairs if before not in skipped]
        mine = [(b, a) for b, a in rust['transforms'].get(k, []) if not b.startswith('S ')]
        # The region moves with the body.
        for b, a in rust['transforms'].get(k, []):
            if b.startswith('S ') and not (k+1 < len(native['bodies'])
                                           and same(b, native['bodies'][k][0], scale)
                                           and same(a, native['bodies'][k+1][0], scale)):
                diffs.add(f'transform{k}_region')
        if any(len(images) != 1 for _, images in pairs) or len(pairs) != len(mine):
            diffs.add(f'transform{k}_arity')
            continue
        for before, (after,) in pairs:
            hit = next((j for j, (b, a) in enumerate(mine) if same(before, b, scale) and same(after, a, scale)), None)
            if hit is None:
                diffs.add(f'transform{k}_pair')
            else:
                mine.pop(hit)
    for k, (sig, _) in enumerate(native['bodies'][1:]):
        if k+1 >= len(rust['bodies']) or not same(sig, rust['bodies'][k+1], scale):
            diffs.add('transformed_body')
    return sorted(diffs)


# ------------------------------------------------------------------ running

def run(executable, text, env):
    started = time.monotonic()
    try:
        result = subprocess.run([str(executable)], input=text+'\n', text=True, capture_output=True,
                                timeout=TIMEOUT, env=env)
        return {'exit_code': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr,
                'seconds': time.monotonic()-started}
    except subprocess.TimeoutExpired as error:
        decode = lambda v: v.decode(errors='replace') if isinstance(v, bytes) else v or ''
        return {'exit_code': None, 'timeout_seconds': TIMEOUT, 'stdout': decode(error.stdout),
                'stderr': decode(error.stderr), 'seconds': time.monotonic()-started}


def build(prefix, output, first):
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
    link = ['otool', '-L', str(executable)] if sys.platform == 'darwin' else ['ldd', str(executable)]
    linked = subprocess.run(link, text=True, capture_output=True, check=True, env=env)
    loaded_text = linked.stdout+linked.stderr
    if sys.platform == 'darwin':
        loaded_text = subprocess.run([str(executable)], input=first+'\n', text=True, capture_output=True,
                                     timeout=TIMEOUT, check=True,
                                     env=dict(env, DYLD_PRINT_LIBRARIES='1')).stderr
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    for toolkit in ['TKPrim', 'TKTopAlgo', 'TKBRep']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded history library: '+toolkit)
    return executable, env, loaded, command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/history-oracle')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    original_capture()
    for module in (generate_identity_fixtures, generate_history_fixtures):
        for name, text in module.generate()[1].items():
            if text != (ROOT/'rust/fixtures'/name).read_text():
                raise ValueError('independent fixture regeneration changed: '+name)
    cases = chosen()
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'history_probe'],
                   cwd=ROOT, check=True)
    probe = subprocess.run([str(ROOT/'target/release/examples/history_probe')],
                           input='\n'.join(encode_case(c) for c in cases)+'\n',
                           text=True, capture_output=True, timeout=600, check=True).stdout
    (output/'rust.txt').write_text(probe)
    rust = parse_rust(probe)
    if list(rust) != [c.name for c in cases]:
        raise ValueError('Rust probe cases changed')
    rows = [native_case(c) for c in cases]
    executable, env, loaded, command = build(prefix, output, rows[0])
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'cases': len(cases), 'native_timeout_seconds': TIMEOUT,
              'native_comparison_budget': {'relative': RELATIVE}, 'native_seconds': {},
              'structure_only_entities': {}, 'matches': [], 'reviewed_differences': [], 'failures': []}
    observations, oracle = {}, None
    for c, text in zip(cases, rows):
        record = run(executable, text, env)
        observations[c.name] = record
        write(output/'native.json', observations)
        report['native_seconds'][c.name] = round(record['seconds'], 3)
        if record['exit_code'] != 0:
            reason = 'native timeout' if record['exit_code'] is None else 'native process failed'
            report['failures'].append({'case': c.name, 'reason': reason, 'record': record})
            continue
        oracle = oracle or next(iter(record['stderr'].splitlines()), None)
        try:
            parsed = parse_native(record['stdout'], c.name)
            differences = compare(parsed, rust[c.name])
            if parsed['structure'].get('-'):
                report['structure_only_entities'][c.name] = len(parsed['structure']['-'])
        except (ValueError, IndexError) as error:
            report['failures'].append({'case': c.name, 'reason': str(error)})
            continue
        if not differences:
            report['matches'].append(c.name)
            continue
        evidence = {'case': c.name, 'source_reference': SOURCE, 'oracle': oracle,
                    'input_sha256': sha(text), 'native_stdout_sha256': sha(record['stdout']),
                    'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(review or evidence)
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'input_sha256': sha('\n'.join(rows)+'\n'), 'source_sha256': digest(SOURCE_FILE),
                'probe_sha256': digest(executable), 'observations_sha256': digest(output/'native.json'),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, (list, dict)) else v for k, v in report.items()
                      if k != 'native_seconds'}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
