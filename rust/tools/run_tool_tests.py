#!/usr/bin/env python3
"""The rust/tools unit tests (`unittest discover -s rust/tools -p 'test_*.py'`)
split into shards by test module, for the kernel workflow's bridge-test jobs.

Every shard discovers the whole suite exactly as `unittest discover` does,
then runs the modules the plan gives it: the plan places the modules,
heaviest first by their measured seconds, each on the least loaded shard,
so every discovered module (and every test in it) runs in exactly one
shard. `--list` prints the discovered test ids, or one shard's.
"""
import argparse
from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parent
SHARDS = 4
# Seconds per test module in the upstream job's discover on CI (run
# 39039288617, 1,862 s in all); a module not listed counts as DEFAULT_SECONDS.
WEIGHTS = {
    'test_brep_oracle': 302, 'test_given_curved_boolean_reference': 231, 'test_spline_sphere_boolean_reference': 124,
    'test_torus_parts_boolean_reference': 73, 'test_given_boolean_reference': 65, 'test_torus_curved_boolean_reference': 65,
    'test_one_sphere_boolean_reference': 63, 'test_cones_boolean_reference': 62, 'test_flush_leaves_boolean_reference': 44,
    'test_deep_trees_boolean_reference': 39, 'test_imported_arcs_boolean_reference': 39,
    'test_turned_crossings_boolean_reference': 39, 'test_rim_leaves_boolean_reference': 38,
    'test_shallow_lens_reference': 36, 'test_degree_elevation_reference': 35, 'test_nested_leaves_boolean_reference': 34,
    'test_imported_boolean_reference': 33, 'test_piece_trees_boolean_reference': 33,
    'test_turned_corners_boolean_reference': 30, 'test_given_met_boolean_reference': 28,
    'test_imported_pieces_boolean_reference': 28, 'test_spheres_turned_boolean_reference': 27,
    'test_chained_boolean_reference': 25, 'test_turned_fillets_boolean_reference': 24,
    'test_imported_joints_boolean_reference': 22, 'test_spline_any_boolean_reference': 22,
    'test_spline_proximity_oracle': 22, 'test_draw_bridge': 21, 'test_turned_tilted_boolean_reference': 21,
    'test_torus_boolean_reference': 20, 'test_primitive_oracle': 19, 'test_spline_parallel_boolean_reference': 18,
    'test_turned_cone_fillets_boolean_reference': 18, 'test_crossed_leaves_boolean_reference': 17,
    'test_polyhedra_curved_boolean_reference': 17, 'test_turned_profiles_boolean_reference': 15,
    'test_split_pieces_boolean_reference': 14, 'test_two_leaves_boolean_reference': 12, 'test_fuzz_runner': 10,
    'test_inside_first_boolean_reference': 10, 'test_primitive_chains_boolean_reference': 10,
    'test_piece_forms_boolean_reference': 8, 'test_tessellation_reference': 8, 'test_turned_flats_boolean_reference': 7,
    'test_torus_segment_boolean_reference': 6, 'test_spline_cone_boolean_reference': 5, 'test_spline_oracle': 5,
    'test_capped_boolean_reference': 4, 'test_split_merge_oracle': 4, 'test_cone_boolean_reference': 3,
    'test_history_oracle': 2, 'test_cones_loops_boolean_reference': 1, 'test_imported_polyhedra_boolean_reference': 1,
    'test_sphere_boolean_reference': 1, 'test_spheres_boolean_reference': 1, 'test_spline_crossing_boolean_reference': 1,
    # Under a second each.
    **dict.fromkeys(['test_bezier_editing_oracle', 'test_brep_io_oracle', 'test_degree_elevation_oracle',
                     'test_knot_editing_oracle', 'test_linear_sets_oracle', 'test_one_sphere_incidence_boolean_reference',
                     'test_plane_parts_boolean_reference', 'test_platform_records', 'test_prism_leaves_boolean_reference',
                     'test_procedural_boolean_reference', 'test_proximity_oracle', 'test_revolve_history',
                     'test_spline_jet_reviews', 'test_spline_linear_oracle', 'test_spline_plane_oracle',
                     'test_stack_leaves_boolean_reference', 'test_step_oracle', 'test_surface_editing_oracle',
                     'test_surface_knot_oracle', 'test_turned_boolean_reference', 'test_run_tool_tests'], 0),
    # Measured on the development Mac (about 20 s, 2.6 times that on CI).
    'test_generator_groups': 50,
}
DEFAULT_SECONDS = 30


def discover():
    """The suite `python -m unittest discover -s rust/tools -p 'test_*.py'` runs."""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    return unittest.TestLoader().discover(str(TOOLS), pattern='test_*.py', top_level_dir=str(TOOLS))


def cases(suite):
    for test in suite:
        if isinstance(test, unittest.TestSuite):
            yield from cases(test)
        else:
            yield test


def module_of(test):
    # A module that failed to import is discover's _FailedTest named after it.
    if type(test).__name__ == '_FailedTest':
        return test.id().rsplit('.', 1)[-1]
    return type(test).__module__


def modules_of(suite):
    return sorted({module_of(t) for t in cases(suite)})


def plan(modules, shards=SHARDS):
    """Each module on exactly one shard: heaviest first, each to the least
    loaded shard (the lowest index on a tie)."""
    loads, parts = [0]*shards, [[] for _ in range(shards)]
    for m in sorted(modules, key=lambda m: (-WEIGHTS.get(m, DEFAULT_SECONDS), m)):
        k = min(range(shards), key=lambda k: (loads[k], k))
        loads[k] += WEIGHTS.get(m, DEFAULT_SECONDS)
        parts[k].append(m)
    return [sorted(p) for p in parts]


def shard_suite(suite, shard, shards=SHARDS):
    mine = set(plan(modules_of(suite), shards)[shard])
    return unittest.TestSuite(t for t in cases(suite) if module_of(t) in mine)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shard', type=int, help='run (or list) this shard only')
    parser.add_argument('--shards', type=int, default=SHARDS)
    parser.add_argument('--list', action='store_true', help='print the test ids instead of running them')
    args = parser.parse_args()
    # As under `python -m unittest`, whose sys.path starts with the working directory.
    sys.path.insert(0, '')
    suite = discover()
    if args.shard is not None:
        if not 0 <= args.shard < args.shards:
            parser.error(f'--shard must be in [0,{args.shards})')
        modules = plan(modules_of(suite), args.shards)[args.shard]
        print(f'shard {args.shard} of {args.shards}: {len(modules)} modules, '
              f'{sum(WEIGHTS.get(m, DEFAULT_SECONDS) for m in modules)} s measured: {" ".join(modules)}',
              file=sys.stderr, flush=True)
        suite = shard_suite(suite, args.shard, args.shards)
    if args.list:
        for test in cases(suite):
            print(test.id())
        return
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == '__main__':
    main()
