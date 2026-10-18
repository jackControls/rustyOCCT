#!/usr/bin/env python3
"""The kernel workflow's fixture-generator groups, run under an input audit.

Each group's commands are the `--check` runs of its independent references.
`run` executes them with generator_audit/sitecustomize.py loaded in every
Python process, which records each file opened for reading, each directory
listed and each path whose existence was asked; after a pass it writes a
marker: those of the inputs inside the checkout (outside target/ and .git/),
each with its contents' hash, existence or listing, plus the workflow file,
this file, the audit hook and the oracle's pinned requirements, the runner
image and the Python version, and the run it passed in. `check` recomputes
every one of those inputs from the checkout; only if all are identical, and
the group's commands too, may a push skip the group (R15: it passed on
identical inputs in the named run, so it counts as green). Any difference,
a missing marker, or a group that started a program other than Python,
runs it. The schedule and manual runs always run every group.
"""
import argparse
import ast
import functools
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT/'rust/tools'
HOOK = TOOLS/'generator_audit'
VENV_PYTHON = 'target/math-oracle-venv/bin/python'
# Always part of a group's inputs, whatever it opened.
FIXED_INPUTS = ['.github/workflows/rust-kernel.yml', 'rust/tools/generator_groups.py',
                'rust/tools/generator_audit/sitecustomize.py', 'rust/tools/math-oracle-requirements.txt']
# Generated during the job (the venv, the turned consistency oracle's
# output) or version control's own: not inputs.
EXCLUDED = ['target', '.git']
MARKER_VERSION = 1

# (interpreter, script, arguments...): 'python3' the runner's, 'venv' the
# independent polynomial oracle's (target/math-oracle-venv).
GROUPS = {
    'algebra': [
        ('python3', 'rust/tools/generate_predicate_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_spatial_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_proximity_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_linear_sets_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_curved_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_spline_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_bezier_editing_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_knot_editing_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_surface_knot_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_degree_elevation_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_exact_spline_intersection_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_surface_editing_fixtures.py', '--check'),
        ('python3', 'rust/tools/generate_surface_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_real_root_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_root_comparison_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_algebraic_image_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spline_proximity_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spline_plane_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spline_quadric_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_closed_clip_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_affine_parameter_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spline_linear_fixtures.py', '--check'),
    ],
    'brep': [
        ('venv', 'rust/tools/generate_brep_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_identity_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_history_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_split_merge_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_attribute_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_brep_io_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_primitive_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_analytic_intersection_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_procedural_intersection_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_torus_curve_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_ruled_curve_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_curve_surface_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_curve_curve_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_split_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_polyhedral_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_curved_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_tessellation_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_step_fixtures.py', '--check'),
    ],
    'booleans': [
        ('venv', 'rust/tools/generate_procedural_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_turned_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_capped_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_sphere_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spheres_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_spheres_turned_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_cone_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_torus_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_torus_segment_boolean_fixtures.py', '--check'),
        ('venv', 'rust/tools/generate_shallow_lens_fixtures.py', '--check'),
    ],
    'cones': [
        ('venv', 'rust/tools/generate_cones_boolean_fixtures.py', '--check'),
    ],
    'cones-loops': [
        ('venv', 'rust/tools/generate_cones_loops_boolean_fixtures.py', '--check'),
    ],
    'tori': [
        ('venv', 'rust/tools/generate_torus_curved_boolean_fixtures.py', '--check'),
    ],
    'torus-parts': [
        ('venv', 'rust/tools/generate_torus_parts_boolean_fixtures.py', '--check'),
    ],
    'chained': [
        ('venv', 'rust/tools/generate_chained_boolean_fixtures.py', '--check'),
    ],
    'spline-any': [
        ('venv', 'rust/tools/generate_spline_any_boolean_fixtures.py', '--check'),
    ],
    'given': [
        ('venv', 'rust/tools/generate_given_boolean_fixtures.py', '--check'),
    ],
    'given-curved': [
        ('venv', 'rust/tools/generate_given_curved_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'given-met': [
        ('venv', 'rust/tools/generate_given_met_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'spline-parallel': [
        ('venv', 'rust/tools/generate_spline_parallel_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'spline-crossing': [
        ('venv', 'rust/tools/generate_spline_crossing_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'imported': [
        ('venv', 'rust/tools/generate_imported_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'imported-arcs': [
        ('venv', 'rust/tools/generate_imported_arcs_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'spline-sphere': [
        ('venv', 'rust/tools/generate_spline_sphere_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'imported-polyhedra': [
        ('venv', 'rust/tools/generate_imported_polyhedra_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'spline-cone': [
        ('venv', 'rust/tools/generate_spline_cone_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'imported-pieces': [
        ('venv', 'rust/tools/generate_imported_pieces_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'split-pieces': [
        ('venv', 'rust/tools/generate_split_pieces_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'one-sphere': [
        ('venv', 'rust/tools/generate_one_sphere_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'one-sphere-incidence': [
        ('venv', 'rust/tools/generate_one_sphere_incidence_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'piece-forms': [
        ('venv', 'rust/tools/generate_piece_forms_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'piece-trees': [
        ('venv', 'rust/tools/generate_piece_trees_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'imported-joints': [
        ('venv', 'rust/tools/generate_imported_joints_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'primitive-chains': [
        ('venv', 'rust/tools/generate_primitive_chains_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'prism-leaves': [
        ('venv', 'rust/tools/generate_prism_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'plane-parts': [
        ('venv', 'rust/tools/generate_plane_parts_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'polyhedra-curved': [
        ('venv', 'rust/tools/generate_polyhedra_curved_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'deep-trees': [
        ('venv', 'rust/tools/generate_deep_trees_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'two-leaves': [
        ('venv', 'rust/tools/generate_two_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'rim-leaves': [
        ('venv', 'rust/tools/generate_rim_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'flush-leaves': [
        ('venv', 'rust/tools/generate_flush_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'stack-leaves': [
        ('venv', 'rust/tools/generate_stack_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'nested-leaves': [
        ('venv', 'rust/tools/generate_nested_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-profiles': [
        ('venv', 'rust/tools/generate_turned_profiles_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-flats': [
        ('venv', 'rust/tools/generate_turned_flats_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-tilted': [
        ('venv', 'rust/tools/generate_turned_tilted_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-fillets': [
        ('venv', 'rust/tools/generate_turned_fillets_boolean_fixtures.py', '--check', '--workers', '4'),
        ('venv', 'rust/tools/generate_turned_consistency_fixtures.py', '--check'),
    ],
    'turned-cone-fillets': [
        ('venv', 'rust/tools/generate_turned_cone_fillets_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-crossings': [
        ('venv', 'rust/tools/generate_turned_crossings_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'turned-corners': [
        ('venv', 'rust/tools/generate_turned_corners_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'crossed-leaves': [
        ('venv', 'rust/tools/generate_crossed_leaves_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
    'inside-first': [
        ('venv', 'rust/tools/generate_inside_first_boolean_fixtures.py', '--check', '--workers', '4'),
    ],
}


def command_argv(command):
    interpreter, *rest = command
    return [VENV_PYTHON if interpreter == 'venv' else 'python3', *rest]


def relative(path, root=ROOT):
    """The checkout-relative path of an absolute one, or None outside it,
    in an excluded directory or in a bytecode cache (its source instead)."""
    try:
        rel = Path(path).relative_to(root)
    except ValueError:
        return None
    parts = rel.parts
    if not parts or parts[0] in EXCLUDED:
        return None if parts else '.'
    if '__pycache__' in parts:
        i = parts.index('__pycache__')
        if i+1 < len(parts) and parts[i+1].endswith('.pyc'):
            return str(Path(*parts[:i], parts[i+1].split('.')[0]+'.py'))
        return None
    return str(rel)


def listing(path):
    names = sorted(n for n in os.listdir(path) if n != '__pycache__' and not n.endswith('.pyc'))
    return hashlib.sha256('\n'.join(names).encode()).hexdigest()


def describe(rel, root=ROOT, listed=False):
    """What a group could have learnt from a path: its kind, a file's
    contents and, for a listed directory, its names."""
    path = root/rel
    if path.is_symlink() and not path.exists():
        return 'dangling'
    if path.is_dir():
        return 'dir:'+listing(path) if listed else 'dir'
    if path.is_file():
        return 'file:'+hashlib.sha256(path.read_bytes()).hexdigest()
    if not path.exists():
        return 'missing'
    return 'other'


# A root a group's modules name in their path constants is keyed as a whole
# tree, every file's contents, besides the runtime record: when in doubt,
# include. The fixtures stay per file (the record), and target/ holds only
# what the job itself makes.
TREE_ROOTS = {'kernel': 'rust/kernel', 'data': 'data', 'fuzz regressions': 'rust/fuzz/regressions'}
RECORDED_ROOTS = {'fixtures', 'target'}


def tree(rel, root=ROOT):
    digest = hashlib.sha256()
    base = root/rel
    for path in sorted(base.rglob('*')) if base.exists() else []:
        if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc':
            digest.update(f'{path.relative_to(base)}\0'.encode())
            digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def key_trees(group):
    modules = sorted({m for c in GROUPS[group] for m in tool_modules(c[1])})
    roots, _ = static_audit(modules)
    return sorted(TREE_ROOTS[r] for r in roots if r in TREE_ROOTS)


def environment():
    return {'image': f"{os.environ.get('ImageOS','local')}-{os.environ.get('ImageVersion','local')}",
            'python': platform.python_version(), 'system': platform.system()}


def is_python(program):
    return bool(re.fullmatch(r'python[\d.]*', Path(str(program)).name))


def read_audit(directory, root=ROOT):
    """{rel: listed?} of the inputs the audit logs name, and the programs started."""
    inputs, programs = {}, set()
    for log in sorted(Path(directory).glob('*.log')):
        for line in log.read_text(errors='surrogateescape').splitlines():
            kind, _, path = line.partition('\t')
            if kind == 'exec':
                programs.add(path)
                continue
            rel = relative(path, root)
            if rel is None:
                continue
            inputs[rel] = inputs.get(rel, False) or kind == 'list'
    for rel in FIXED_INPUTS:
        inputs.setdefault(rel, False)
    return inputs, sorted(programs)


def build_marker(group, audit_dir, root=ROOT, passed_in=None):
    inputs, programs = read_audit(audit_dir, root)
    foreign = [p for p in programs if not is_python(p)]
    return {'version': MARKER_VERSION, 'group': group, 'commands': [list(c) for c in GROUPS[group]],
            'environment': environment(), 'passed_in': passed_in or {}, 'programs': programs,
            'memoisable': not foreign,
            'trees': {rel: tree(rel, root) for rel in key_trees(group)},
            'inputs': {rel: {'listed': listed, 'state': describe(rel, root, listed)} for rel, listed in sorted(inputs.items())}}


def check_marker(group, marker, root=ROOT):
    """The reasons the group must run, empty if it passed on identical inputs."""
    reasons = []
    if marker.get('version') != MARKER_VERSION or marker.get('group') != group:
        return ['the marker is of another version or group']
    if marker['commands'] != [list(c) for c in GROUPS[group]]:
        reasons.append('the group\'s commands changed')
    if marker['environment'] != environment():
        reasons.append(f"the environment changed: {marker['environment']} then, {environment()} now")
    if not marker.get('memoisable'):
        reasons.append(f"the group started other programs: {marker.get('programs')}")
    if sorted(marker.get('trees', {})) != key_trees(group):
        reasons.append('the group\'s keyed trees changed')
    changed_trees = [rel for rel, digest in marker.get('trees', {}).items() if tree(rel, root) != digest]
    if changed_trees:
        reasons.append(f'{changed_trees} changed')
    changed = [rel for rel, entry in marker['inputs'].items() if describe(rel, root, entry['listed']) != entry['state']]
    if changed:
        reasons.append(f'{len(changed)} input(s) changed, first {changed[:5]}')
    if not any(rel == FIXED_INPUTS[0] for rel in marker['inputs']):
        reasons.append('the marker does not name the workflow file')
    return reasons


def run_group(group, audit_dir, marker_path):
    if audit_dir.exists():
        shutil.rmtree(audit_dir)
    audit_dir.mkdir(parents=True)
    env = dict(os.environ, GENERATOR_AUDIT_DIR=str(audit_dir))
    env['PYTHONPATH'] = os.pathsep.join(filter(None, [str(HOOK), os.environ.get('PYTHONPATH')]))
    for command in GROUPS[group]:
        argv = command_argv(command)
        print(f'::group::{" ".join(argv)}', flush=True)
        code = subprocess.call(argv, cwd=ROOT, env=env)
        print('::endgroup::', flush=True)
        if code != 0:
            print(f'{" ".join(argv)} failed with exit code {code}', flush=True)
            return code
    passed_in = {k: os.environ.get(v) for k, v in [('run_id', 'GITHUB_RUN_ID'), ('run_attempt', 'GITHUB_RUN_ATTEMPT'),
                                                   ('sha', 'GITHUB_SHA'), ('ref', 'GITHUB_REF'), ('event', 'GITHUB_EVENT_NAME')]}
    marker = build_marker(group, audit_dir, passed_in=passed_in)
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    marker_path.write_text(json.dumps(marker, indent=1)+'\n')
    print(f"{group} passed; {len(marker['inputs'])} inputs recorded, programs {marker['programs']}, "
          f"memoisable {marker['memoisable']}", flush=True)
    return 0


def output(**values):
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as f:
            for k, v in values.items():
                f.write(f'{k}={v}\n')


def check_group(group, marker_path):
    if not marker_path.exists():
        print(f'{group}: no earlier pass recorded; running it', flush=True)
        output(skip='false')
        return
    marker = json.loads(marker_path.read_text())
    reasons = check_marker(group, marker)
    if reasons:
        print(f'{group}: running it, since ' + '; '.join(reasons), flush=True)
        output(skip='false')
        return
    run = marker['passed_in']
    note = (f"{group}: skipped; it passed on identical inputs ({len(marker['inputs'])} files and directories, "
            f"the workflow and the environment) in run {run.get('run_id')} attempt {run.get('run_attempt')} "
            f"at {run.get('sha')} ({run.get('event')}, {run.get('ref')})")
    print(note, flush=True)
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as f:
            f.write(note+'\n')
    output(skip='true')


@functools.lru_cache(maxsize=None)
def parsed(name):
    return ast.parse((TOOLS/f'{name}.py').read_text())


@functools.lru_cache(maxsize=None)
def tool_modules(script):
    """The rust/tools modules a script imports, transitively (static)."""
    seen, todo = set(), [Path(script).stem]
    while todo:
        name = todo.pop()
        if name in seen or not (TOOLS/f'{name}.py').exists():
            continue
        seen.add(name)
        for node in ast.walk(parsed(name)):
            if isinstance(node, ast.Import):
                todo += [a.name.split('.')[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                todo.append(node.module.split('.')[0])
    return sorted(seen)


PATH_ROOTS = [('rust/kernel', 'kernel'), ('kernel/', 'kernel'), ('rust/fixtures', 'fixtures'), ('fixtures', 'fixtures'),
              ('data/', 'data'), ('regressions', 'fuzz regressions'), ('fuzz/', 'fuzz regressions'),
              ('rust/fuzz', 'fuzz regressions'), ('target/', 'target')]


@functools.lru_cache(maxsize=None)
def module_audit(name):
    """The roots one module's path-like string constants name, and whether
    it imports subprocess."""
    roots, spawns = set(), False
    for node in ast.walk(parsed(name)):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and re.fullmatch(r'[\w.{}*/-]+', node.value):
            for prefix, root in PATH_ROOTS:
                if node.value.startswith(prefix) or node.value == prefix.rstrip('/'):
                    roots.add(root)
                    break
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [a.name for a in node.names] + ([node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            spawns |= 'subprocess' in names
    return frozenset(roots), spawns


def static_audit(modules):
    """The checkout roots a set of rust/tools modules names in path-like
    string constants, and the modules that can start a program."""
    audits = {m: module_audit(m) for m in modules}
    return sorted(set().union(*(r for r, _ in audits.values()))), sorted(m for m, (_, s) in audits.items() if s)


def audit_table():
    """A static summary per group: its scripts' rust/tools closure, the
    checkout roots it names and the modules that can start a program. The
    marker's runtime record is what a skip relies on; this cross-checks it."""
    rows = []
    for group, commands in GROUPS.items():
        modules = sorted({m for c in commands for m in tool_modules(c[1])})
        roots, spawns = static_audit(modules)
        rows.append((group, len(commands), modules, roots, spawns))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    for action in ('check', 'run'):
        p = sub.add_parser(action)
        p.add_argument('group', choices=sorted(GROUPS))
        p.add_argument('--marker', type=Path, required=True)
    sub.add_parser('audit')
    args = parser.parse_args()
    if args.action == 'check':
        check_group(args.group, args.marker)
    elif args.action == 'run':
        raise SystemExit(run_group(args.group, ROOT/'target/generator-audit'/args.group, args.marker))
    else:
        for group, count, modules, roots, spawns in audit_table():
            print(f'{group}\t{count}\t{len(modules)}\t{",".join(roots)}\t{",".join(spawns)}')


if __name__ == '__main__':
    main()
