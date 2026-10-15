#!/usr/bin/env python3
"""Every Boolean set's comparison with the source-pinned SDK, in shards.

The scheduled workflow "Rust Boolean captures on Linux"
(`.github/workflows/rust-boolean-linux.yml`) runs one shard per job on the
SDK its first job builds: `compare_boolean.py` (S9a.1's prisms, whose Linux
record exists), `compare_boolean.py --splines`, `compare_polyhedral.py` and
every `compare_*_boolean.py`, each with `--allow-unrecorded-platform`, so a
capture without a record of the runner's platform yet reports its
observation pending (its `report.json`'s `unrecorded_platform` and
`platform_reviews_pending`) and a capture with one is compared strictly.
Every comparison runs even after one fails; the shard fails if any did.
Each set's output directory (`target/boolean*-oracle/`) holds its
`native-observed.txt` and `report.json`, from which
`take_platform_records.py` takes the records; `--summary` writes the
shard's commands, exit codes and times.

Shards are balanced by SECONDS, each comparison's local wall time (macOS
arm64, the kernel probe already built). A set missing from SECONDS counts
60 s; `test_platform_records.py` requires every set to have its time and
every shard's estimate on a Linux runner (LINUX_FACTOR times, plus
SHARD_OVERHEAD) to stay within BUDGET, half the job's 60-minute limit.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT/'rust/tools'
SHARDS = 3
# Local wall seconds of each comparison, the kernel probe built (2026-10-14,
# all 46 sets 378 s; the native probe's compilation included; the turned
# profiles', flats', tilted bodies', fillets', cone fillets' and crossings'
# 2026-10-15).
SECONDS = {
    'compare_boolean.py': 3,
    'compare_boolean.py --splines': 9,
    'compare_polyhedral.py': 3,
    'compare_capped_boolean.py': 4,
    'compare_chained_boolean.py': 3,
    'compare_cone_boolean.py': 3,
    'compare_cones_boolean.py': 39,
    'compare_cones_loops_boolean.py': 7,
    'compare_curved_boolean.py': 3,
    'compare_deep_trees_boolean.py': 5,
    'compare_flush_leaves_boolean.py': 11,
    'compare_given_boolean.py': 3,
    'compare_given_curved_boolean.py': 7,
    'compare_given_met_boolean.py': 14,
    'compare_imported_arcs_boolean.py': 5,
    'compare_imported_boolean.py': 5,
    'compare_imported_joints_boolean.py': 7,
    'compare_imported_pieces_boolean.py': 12,
    'compare_imported_polyhedra_boolean.py': 6,
    'compare_nested_leaves_boolean.py': 7,
    'compare_one_sphere_boolean.py': 12,
    'compare_one_sphere_incidence_boolean.py': 8,
    'compare_piece_forms_boolean.py': 10,
    'compare_piece_trees_boolean.py': 12,
    'compare_plane_parts_boolean.py': 27,
    'compare_polyhedra_curved_boolean.py': 10,
    'compare_primitive_chains_boolean.py': 10,
    'compare_prism_leaves_boolean.py': 12,
    'compare_procedural_boolean.py': 5,
    'compare_rim_leaves_boolean.py': 6,
    'compare_sphere_boolean.py': 3,
    'compare_spheres_boolean.py': 5,
    'compare_spheres_turned_boolean.py': 5,
    'compare_spline_any_boolean.py': 10,
    'compare_spline_cone_boolean.py': 8,
    'compare_spline_crossing_boolean.py': 14,
    'compare_spline_parallel_boolean.py': 12,
    'compare_spline_sphere_boolean.py': 12,
    'compare_split_pieces_boolean.py': 8,
    'compare_stack_leaves_boolean.py': 7,
    'compare_torus_boolean.py': 6,
    'compare_torus_curved_boolean.py': 12,
    'compare_torus_parts_boolean.py': 15,
    'compare_torus_segment_boolean.py': 5,
    'compare_turned_boolean.py': 7,
    'compare_turned_cone_fillets_boolean.py': 7,
    'compare_turned_crossings_boolean.py': 8,
    'compare_turned_fillets_boolean.py': 9,
    'compare_turned_flats_boolean.py': 5,
    'compare_turned_profiles_boolean.py': 6,
    'compare_turned_tilted_boolean.py': 5,
    'compare_two_leaves_boolean.py': 6,
}
UNKNOWN_SECONDS = 60
# Linux runners against this host (FUZZING.md's triage rule).
LINUX_FACTOR = 2.6
# A shard's fixed cost on a runner: checkout, the venv, the SDK's download
# and unpacking, and the kernel probe's release build (63 s wall, 240 s CPU
# here; about 5 minutes on a runner).
SHARD_OVERHEAD = 8*60
BUDGET = 30*60


def comparisons():
    """[(name, script arguments)] of every Boolean set, in a fixed order."""
    out = [('compare_boolean.py', ['compare_boolean.py']),
           ('compare_boolean.py --splines', ['compare_boolean.py', '--splines']),
           ('compare_polyhedral.py', ['compare_polyhedral.py'])]
    return out+[(p.name, [p.name]) for p in sorted(TOOLS.glob('compare_*_boolean.py'))]


def cost(name):
    return SECONDS.get(name, UNKNOWN_SECONDS)


def shards(count=SHARDS):
    """The comparisons in `count` shards: longest first onto the least
    loaded shard (ties to the lower index), deterministic."""
    loads = [[0, k, []] for k in range(count)]
    for name, argv in sorted(comparisons(), key=lambda c: (-cost(c[0]), c[0])):
        least = min(loads, key=lambda load: (load[0], load[1]))
        least[0] += cost(name)
        least[2].append((name, argv))
    return [load[2] for load in loads]


def estimate(shard):
    """A shard's estimated seconds on a Linux runner."""
    return SHARD_OVERHEAD+LINUX_FACTOR*sum(cost(name) for name, _ in shard)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--shard', type=int, required=True)
    parser.add_argument('--shards', type=int, default=SHARDS)
    parser.add_argument('--occt-root', type=Path)
    parser.add_argument('--sdk-manifest', type=Path)
    parser.add_argument('--list', action='store_true', help='print the shard and its estimate only')
    parser.add_argument('--summary', type=Path, help="the shard's commands, exit codes and times (JSON)")
    args = parser.parse_args()
    if not 0 <= args.shard < args.shards:
        parser.error('--shard must lie in [0, --shards)')
    plan = shards(args.shards)[args.shard]
    if args.list:
        for name, _ in plan:
            print(f'{name}\t{cost(name)} s')
        print(f'{len(plan)} comparisons, about {estimate(plan)/60:.0f} minutes on a Linux runner')
        return
    if args.occt_root is None or args.sdk_manifest is None:
        parser.error('--occt-root and --sdk-manifest are required to run the comparisons')
    github = os.environ.get('GITHUB_ACTIONS') == 'true'
    results = []
    for name, argv in plan:
        command = [sys.executable, str(TOOLS/argv[0]), *argv[1:], '--allow-unrecorded-platform',
                   '--occt-root', str(args.occt_root), '--sdk-manifest', str(args.sdk_manifest)]
        print(('::group::' if github else '== ')+name, flush=True)
        start = time.monotonic()
        code = subprocess.run(command, cwd=ROOT).returncode
        seconds = round(time.monotonic()-start, 1)
        if github:
            print('::endgroup::')
        print(f'{name}: exit {code}, {seconds} s', flush=True)
        results.append({'comparison': name, 'exit_code': code, 'seconds': seconds})
    failed = [r['comparison'] for r in results if r['exit_code'] != 0]
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps({'shard': args.shard, 'shards': args.shards, 'platform': sys.platform,
                                            'comparisons': results, 'failed': failed}, indent=2)+'\n')
    if failed:
        raise SystemExit('failed: '+', '.join(failed))


if __name__ == '__main__':
    main()
