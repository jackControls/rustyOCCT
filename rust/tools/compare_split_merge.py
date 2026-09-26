#!/usr/bin/env python3
"""Source-pinned OCCT height split and stacked fuse beside the Rust histories.

Rust histories are first certified equal to the independent enumeration
(fixture regeneration plus the split_merge tests). For each scenario captured
before implementation, the native probe's BRepTools_History of
BRepAlgoAPI_Splitter and of BRepAlgoAPI_Fuse followed by
ShapeUpgrade_UnifySameDomain is mapped to the Rust relations of the same stage
(`split`, `fuse`, and `composed` for a split whose pieces are fused back) by
geometric signature, under these correspondences:

    Rust Unchanged    native: no record and kept, or Modified into one shape
                      with the same signature (unification rebuilds shapes)
    Rust Split        native Modified into the same pieces (contract 3: OCCT
                      reports a split parent as Modified; by design)
    Rust Merged       native Modified into the merged shape
    Rust Modified     (a composed history) native Modified into that shape
    Rust Deleted      native IsRemoved without images
    Rust region       native solid: Split as above; Merged or Modified when
                      native removes the argument solid and a result body has
                      the target's signature
    Rust Generated    cut edges and vertices: a native Generated image of the
                      parent wall or vertical edge; cut faces: a result face no
                      native query reaches (OCCT's cut face is the tool's image)

OCCT's seams and the vertices only they use (S rows) are structure-only, as in
compare_history.py: their queries, images and result shapes are set aside, and
the per-body counts (native distinct subshapes against the kernel's
synthesized OCCT counts) verify the rule on every case.

In the composed stage native Modified images are compared whatever IsRemoved
says: BRepTools_History::Merge keeps the removal of intermediate pieces.
Every native query must correspond to a Rust relation and every Rust
relation to a native query; every result shape to a Rust output entity and
back (OCCT shares a split's cut face between both pieces, so this is by
signature); result bodies, validity and each body's distinct subshape counts
(the kernel's synthesized OCCT counts) must agree. Differences need a
fingerprinted review; timeouts, crashes and malformed output are never
reviewable.
"""
import argparse
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
from compare_history import same, numbers
from compare_occt import ROOT
from compare_degree_elevation import verify_sdk, verify_loaded_libraries
import generate_split_merge_fixtures as generator

ORIGINAL = ROOT/'rust/fixtures/occt-split-merge-preimplementation'
REVIEWS = ROOT/'rust/fixtures/occt-split-merge-divergences.json'
SOURCE_FILE = ROOT/'rust/tools/occt_split_merge_oracle.cpp'
TOOLKITS = ['TKShHealing', 'TKBO', 'TKPrim', 'TKTopAlgo', 'TKBRep', 'TKGeomAlgo', 'TKGeomBase', 'TKG3d',
            'TKG2d', 'TKMath', 'TKernel']
# Every preflight case finished within 0.05 s; the deadline only bounds hangs.
TIMEOUT = 120


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def review_for(evidence, reviews):
    return next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                 and all(r.get(k) == v for k, v in evidence.items())), None)


def native_scenarios():
    scenarios, files = generator.generate()
    return generator.native_rows(scenarios), files


def original_capture():
    """The pre-implementation observations and their inputs are unchanged."""
    metadata = json.loads((ORIGINAL/'capture.json').read_text())
    if (metadata['source_reference'] != SOURCE
            or metadata['rust_split_or_fuse_implementation_exists']
            or any(line[3:].startswith('rust/kernel') for line in metadata['rust_worktree_uncommitted'])):
        raise ValueError('original capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.json')]:
        if metadata[key] != digest(ORIGINAL/name):
            raise ValueError('original native evidence changed: '+name)
    rows, _ = native_scenarios()
    if metadata['cases'] != [s.name for s, _ in rows]:
        raise ValueError('native scenarios changed')
    same_inputs((ORIGINAL/'inputs.txt').read_text(), '\n'.join(text for _, text in rows)+'\n')


# ------------------------------------------------------------------ parsing

def _stage():
    return {'queries': [], 'bodies': [], 'counts': {}, 'outputs': [], 'generated': [], 'structure': set()}


def parse_native(stdout, name):
    lines = stdout.splitlines()
    if not lines or lines[0] != f'R {name}' or lines[-1] != 'end':
        raise ValueError(f'malformed native output for {name}')
    out = {'inputs_valid': True, 'stages': {}}
    i = 1
    while i < len(lines)-1:
        w = lines[i].split()
        if w[0] == 'I':
            out['inputs_valid'] &= w[-1] == '1'
        elif w[0] == 'Q':
            st = out['stages'].setdefault(w[1], _stage())
            bar = w.index('|')
            mod, gen, removed = int(w[6]), int(w[8]), w[10] == '1'
            images = [lines[i+1+j].split(None, 1) for j in range(mod+gen)]
            if any(len(x) != 2 or x[0] not in 'MG' for x in images):
                raise ValueError(f'malformed native images for {name}')
            st['queries'].append({'tag': w[3], 'sig': ' '.join(w[bar+1:]), 'removed': removed,
                                  'modified': [s for t, s in images if t == 'M'],
                                  'generated': [s for t, s in images if t == 'G']})
            i += mod+gen
        elif w[0] == 'B':
            out['stages'].setdefault(w[1], _stage())['bodies'].append((' '.join(w[3:-1]), w[-1] == '1'))
        elif w[0] == 'N':
            out['stages'].setdefault(w[1], _stage())['counts'][int(w[2])] = ' '.join(w[3:])
        elif w[0] == 'O':
            out['stages'].setdefault(w[1], _stage())['outputs'].append((' '.join(w[2:-1]), w[-1]))
        elif w[0] == 'S':
            out['stages'].setdefault(w[1], _stage())['structure'].add(' '.join(w[2:]))
        else:
            raise ValueError(f'unknown native row {lines[i]!r}')
        i += 1
    return out


def parse_rust(text):
    cases, current = {}, None
    for line in text.splitlines():
        w = line.split()
        if w[0] == 'R':
            current = cases[w[1]] = {}
        elif w[0] in ('Q', 'G'):
            st = current.setdefault(w[1], _stage())
            bar, arrow = w.index('|'), w.index('->')
            source = ' '.join(w[bar+1:arrow])
            targets = [t.strip() for t in ' '.join(w[arrow+1:]).split(' ; ') if t.strip()]
            if w[0] == 'Q':
                st['queries'].append({'kind': w[2], 'sig': source, 'targets': targets})
            else:
                st['generated'].append({'role': w[2], 'parents': source.split(' ; '), 'target': targets[0]})
        elif w[0] == 'B':
            current.setdefault(w[1], _stage())['bodies'].append(' '.join(w[3:]))
        elif w[0] == 'C':
            current.setdefault(w[1], _stage())['counts'][int(w[2])] = ' '.join(w[3:])
        elif w[0] == 'O':
            current.setdefault(w[1], _stage())['outputs'].append(' '.join(w[2:]))
    return cases


# ------------------------------------------------------------------ comparison

def _scale(rust):
    values = [abs(v) for st in rust.values() for s in st['outputs'] for v in numbers(s)]
    return max(values+[1.0])


def _multiset_match(xs, ys, scale):
    """Every x pairs with a distinct y of equal signature."""
    left = list(ys)
    for x in xs:
        hit = next((k for k, y in enumerate(left) if same(x, y, scale)), None)
        if hit is None:
            return False
        left.pop(hit)
    return not left


def _fits(r, q, stage, scale):
    kind, targets = r['kind'], r['targets']
    solid = r['sig'].startswith('S ')
    modified = q['modified'] if stage == 'composed' or not q['removed'] else None
    if kind == 'unchanged':
        return (not q['removed'] and (not q['modified'] or
                (len(q['modified']) == 1 and same(q['modified'][0], r['sig'], scale))))
    if kind == 'deleted':
        return q['removed'] and not q['modified']
    if kind == 'split':
        return modified is not None and _multiset_match(modified, targets, scale)
    if kind in ('merged', 'modified'):
        if solid and q['removed'] and not q['modified']:
            return 'body'
        return bool(modified) and len(modified) == 1 and same(modified[0], targets[0], scale)
    return False


def compare(native, rust):
    diffs = set()
    scale = _scale(rust)
    if not native['inputs_valid']:
        diffs.add('native_invalid_input')
    if sorted(native['stages']) != sorted(rust):
        return sorted(diffs | {'stages'})
    for stage, n in native['stages'].items():
        r = rust[stage]
        # Seams and the vertices only they use are structure the seamless
        # model does not have; counts verify the rule below.
        skip = n['structure']
        n = dict(n, queries=[dict(q, generated=[x for x in q['generated'] if x not in skip])
                             for q in n['queries'] if q['sig'] not in skip],
                 outputs=[(s, state) for s, state in n['outputs'] if s not in skip])
        # Result bodies, validity and counts.
        if len(n['bodies']) != len(r['bodies']):
            diffs.add(f'{stage}_body_count')
        if not all(ok for _, ok in n['bodies']):
            diffs.add(f'{stage}_native_invalid')
        for k, sig in enumerate(r['bodies']):
            hits = [j for j, (s, _) in enumerate(n['bodies']) if same(s, sig, scale)]
            if len(hits) != 1:
                diffs.add(f'{stage}_body_signature')
            elif n['counts'].get(hits[0]) != r['counts'].get(k):
                diffs.add(f'{stage}_counts')
        # Relations of every input entity.
        used = set()
        for rel in r['queries']:
            candidates = [j for j, q in enumerate(n['queries'])
                          if q['tag'][0] == ('s' if rel['sig'].startswith('S ') else rel['sig'][0].lower())
                          and same(q['sig'], rel['sig'], scale)]
            fit = next((j for j in candidates if j not in used and _fits(rel, n['queries'][j], stage, scale)), None)
            if fit is None:
                diffs.add(f"{stage}_{rel['kind']}_{rel['sig'][0]}")
                continue
            used.add(fit)
            if _fits(rel, n['queries'][fit], stage, scale) == 'body' and not any(
                    same(b, rel['targets'][0], scale) for b, _ in n['bodies']):
                diffs.add(f'{stage}_region_body')
        if len(used) != len(n['queries']):
            diffs.add(f'{stage}_unmatched_native_query')
        # Generated entities.
        unreached = [s for s, state in n['outputs'] if state == 'none']
        for g in r['generated']:
            if g['role'] == 'CutFace':
                ok = any(same(s, g['target'], scale) for s in unreached)
            else:
                parents = [q for q in n['queries'] if same(q['sig'], g['parents'][0], scale)]
                ok = any(same(x, g['target'], scale) for q in parents for x in q['generated'])
            if not ok:
                diffs.add(f"{stage}_generated_{g['role']}")
        for q in n['queries']:
            for x in q['generated']:
                if not any(same(x, g['target'], scale) for g in r['generated']):
                    diffs.add(f'{stage}_native_generated_unmatched')
        # Every result shape corresponds to an output entity and back.
        rust_outputs = [s for s in r['outputs'] if not s.startswith('S ')]
        for s, _ in n['outputs']:
            if not any(same(s, x, scale) for x in rust_outputs):
                diffs.add(f'{stage}_native_output_unmatched')
        for x in rust_outputs:
            if not any(same(s, x, scale) for s, _ in n['outputs']):
                diffs.add(f'{stage}_rust_output_unmatched')
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
    for toolkit in ['TKBO', 'TKShHealing', 'TKPrim', 'TKBRep']:
        if not any(Path(p['path']).name.startswith('lib'+toolkit+'.') for p in loaded):
            raise ValueError('missing loaded split/fuse library: '+toolkit)
    return executable, env, loaded, command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/split-merge-oracle')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    original_capture()
    rows, files = native_scenarios()
    for name, text in files.items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'split_merge_probe'],
                   cwd=ROOT, check=True)
    probe = subprocess.run([str(ROOT/'target/release/examples/split_merge_probe')],
                           input='\n'.join(generator.encode(generator.native_scenario(s)) for s, _ in rows)+'\n',
                           text=True, capture_output=True, timeout=600, check=True).stdout
    (output/'rust.txt').write_text(probe)
    rust = parse_rust(probe)
    if list(rust) != [s.name for s, _ in rows]:
        raise ValueError('Rust probe scenarios changed')
    executable, env, loaded, command = build(prefix, output, rows[0][1])
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'cases': len(rows), 'native_timeout_seconds': TIMEOUT,
              'native_seconds': {}, 'matches': [], 'reviewed_differences': [], 'failures': []}
    observations, oracle = {}, None
    for s, text in rows:
        record = run(executable, text, env)
        observations[s.name] = record
        write(output/'native.json', observations)
        report['native_seconds'][s.name] = round(record['seconds'], 3)
        if record['exit_code'] != 0:
            reason = 'native timeout' if record['exit_code'] is None else 'native process failed'
            report['failures'].append({'case': s.name, 'reason': reason, 'record': record})
            continue
        oracle = oracle or next(iter(record['stderr'].splitlines()), None)
        try:
            differences = compare(parse_native(record['stdout'], s.name), rust[s.name])
        except (ValueError, IndexError, KeyError) as error:
            report['failures'].append({'case': s.name, 'reason': f'{type(error).__name__}: {error}'})
            continue
        if not differences:
            report['matches'].append(s.name)
            continue
        evidence = {'case': s.name, 'source_reference': SOURCE, 'oracle': oracle,
                    'input_sha256': sha(text), 'native_stdout_sha256': sha(record['stdout']),
                    'differences': differences}
        review = review_for(evidence, reviews)
        report['reviewed_differences' if review else 'failures'].append(review or evidence)
    metadata = {'source_reference': SOURCE, 'oracle': oracle, 'sdk_manifest_sha256': digest(args.sdk_manifest),
                'input_sha256': sha('\n'.join(t for _, t in rows)+'\n'), 'source_sha256': digest(SOURCE_FILE),
                'probe_sha256': digest(executable), 'observations_sha256': digest(output/'native.json'),
                'loaded_libraries': loaded, 'build_command': command}
    write(output/'capture.json', metadata)
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, (list, dict)) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
