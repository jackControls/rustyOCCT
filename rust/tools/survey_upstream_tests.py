#!/usr/bin/env python3
"""Run unregistered original cases on both backends, to decide what to register.

A survey never changes the coverage contract: it runs each listed case in a
fresh process on the Rust adapter and native DRAW, with the fetched public
dataset when present, and reports every status, every unsupported construct
by name and every divergence. `--restore-only` selects the cases that only
restore data files and run check commands (REVIEW_NOTES.md S2). `--register`
then records in `fixtures/upstream-draw.json` every surveyed case that both
backends evaluated, with the SHA-256 of the case and its whole upstream
context and the observed statuses as expectations.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import shutil
import sys

from run_upstream_tests import (DATASET, MANIFEST, ROOT, EVALUATED, build_worker, case_files,
                                classify_missing, dataset_inventory, digest, run_case,
                                source_files)

CHECKS = {'checkshape', 'checkprops', 'checknbshapes', 'checkmaxtol', 'checkview', 'checkfreebounds',
          'checkreal', 'checkgravitycenter', 'checklength', 'checkarea', 'checktrinfo', 'checkcolor',
          'checkdump', 'checkpoint'}
TCL = {'puts', 'set', 'if', 'else', 'elseif', 'foreach', 'for', 'incr', 'expr', 'catch', 'regexp',
       'string', 'llength', 'lindex', 'list', 'lappend', 'append', 'proc', 'return', 'info', 'dict',
       'unset', 'global', 'upvar', 'cpulimit', 'pload', 'whatis', 'nbshapes', 'vprops', 'sprops',
       'lprops', 'dump', 'isdraw', 'dchrono', 'chrono'}


def viewer_commands():
    text = (ROOT / 'rust/fixtures/draw-viewer-commands.txt').read_text()
    return {line.strip() for line in text.splitlines() if line.strip() and not line.startswith('#')}


def restore_only():
    """Cases whose commands are `restore`, check procedures, viewer commands
    and plain Tcl, with at least one check."""
    allowed = {'restore'} | CHECKS | TCL | viewer_commands()
    out = []
    for path in case_files():
        # tests/<group>/<grid>/<case>; deeper files are data of other grids.
        if len(path.relative_to(ROOT / 'tests').parts) != 3:
            continue
        text = path.read_text(errors='replace')
        if 'restore' not in text:
            continue
        words = set()
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith('#') or stripped.startswith('}'):
                continue
            first = re.match(r'([A-Za-z_][\w:]*)', stripped)
            if first:
                words.add(first.group(1))
            words.update(m.group(1) for m in re.finditer(r'\[\s*([A-Za-z_]\w*)', stripped))
        words.discard('locate_data_file')
        if 'restore' in words and words <= allowed and words & CHECKS:
            out.append(str(path.relative_to(ROOT)))
    return out


def context(case):
    """The case, its begin/end scripts and parse rules, as recorded paths."""
    sources = source_files(case)
    group, grid = Path(case).parent.relative_to('tests').parts
    rules = [p for p in [ROOT / 'tests/parse.rules', ROOT / 'tests' / group / 'parse.rules',
                         ROOT / 'tests' / group / grid / 'parse.rules'] if p.is_file()]
    return sources, [str(p.relative_to(ROOT)) for p in sources + rules]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-only', action='store_true')
    parser.add_argument('--case', action='append', default=[])
    parser.add_argument('--draw-exe', default=shutil.which('DRAWEXE') or shutil.which('occt-draw')
                        or '/opt/homebrew/opt/opencascade/bin/DRAWEXE')
    parser.add_argument('--tclsh', default=shutil.which('tclsh') or 'tclsh')
    parser.add_argument('--timeout', type=float, default=60.0)
    parser.add_argument('--output', type=Path, default=ROOT / 'target/upstream-survey')
    parser.add_argument('--register', action='store_true')
    args = parser.parse_args()
    cases = (restore_only() if args.restore_only else []) + args.case
    inventory = dataset_inventory()
    data_dirs = [ROOT / 'data'] + ([DATASET] if inventory is not None else [])
    worker = build_worker()
    rows = []
    for k, case in enumerate(cases):
        sources, _ = context(case)
        names = Path(case).relative_to('tests').parts
        row = {'case': case}
        for backend in ['rust', 'occt']:
            result = run_case(backend, sources, args.output / backend / case, worker, args.draw_exe,
                              args.tclsh, args.timeout, data_dirs, names)
            result = classify_missing(result, inventory)
            row[backend] = result['status']
            row[backend + '_unsupported'] = result.get('unsupported', '')
            row[backend + '_error'] = result.get('error', '')
        rows.append(row)
        print(f"{k + 1:4}/{len(cases)} {row['rust']:15} {row['occt']:15} {case}", flush=True)
    evaluated = [r for r in rows if r['rust'] in EVALUATED and r['occt'] in EVALUATED]
    constructs = Counter()
    for r in rows:
        for line in r['rust_unsupported'].splitlines():
            m = re.search(r'unsupported constructs (\S+)', line)
            for item in (m.group(1).split(',') if m else []):
                constructs[item.split('=')[0]] += 1
    summary = {
        'dataset': 'fetched' if inventory is not None else 'not fetched',
        'cases': len(rows),
        'pairs': dict(Counter(f"{r['rust']} / {r['occt']}" for r in rows).most_common()),
        'evaluated_on_both_backends': len(evaluated),
        'rust_unsupported_restore_by_construct': dict(constructs.most_common()),
        'rust_unsupported_other': dict(Counter(
            line.split(':')[0] if not line.startswith('restore') else 'restore'
            for r in rows for line in r['rust_unsupported'].splitlines()
            if 'unsupported constructs' not in line).most_common(20)),
        'native_evaluated_rust_failed': [r['case'] for r in rows
                                         if r['occt'] in EVALUATED and r['rust'] in {'failed', 'timeout'}],
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'survey.json').write_text(json.dumps({'summary': summary, 'cases': rows}, indent=1) + '\n')
    print(json.dumps(summary, indent=1))
    if args.register:
        manifest = json.loads(MANIFEST.read_text())
        known = {c['path'] for c in manifest['cases']}
        for r in evaluated:
            if r['case'] in known:
                continue
            _, paths = context(r['case'])
            for path in paths:
                manifest['sources'].setdefault(path, digest(ROOT / path))
            manifest['cases'].append({
                'path': r['case'], 'data': True, 'expected_rust': r['rust'], 'expected_occt': r['occt'],
                'purpose': 'Original restore-and-check case on the public dataset (REVIEW_NOTES.md S2).'})
        manifest['sources'] = dict(sorted(manifest['sources'].items()))
        MANIFEST.write_text(json.dumps(manifest, indent=2) + '\n')
        print(f'registered {len(evaluated)} evaluated cases')
    return 0


if __name__ == '__main__':
    sys.exit(main())
