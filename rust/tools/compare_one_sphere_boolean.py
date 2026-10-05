#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.3c.1 reference: two pieces of one sphere
(bodies OCCT wrote to `.brep` files, two with their rims divided as the
DRAW survey's `so1` and `so4`) given to Booleans
(`generate_one_sphere_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, read from `rust/fixtures`) and the comparison are
S9e.4a's (`compare_imported_boolean.py`, through
`compare_boolean.make_set`) on S9e.4b.3c.1's fixtures, capture
(`fixtures/occt-boolean-one-sphere-preimplementation`) and reviews
(`fixtures/occt-boolean-one-sphere-divergences.json`), the kernel's
enclosures against the reference's totals with S9e.4a's slack of 1e-12
relative. The capture was taken before S9e.4b.3c.1's code existed, while
the kernel refused two inputs on one sphere where either is a piece
(`pieces.rs`'s `OutOfDomain("faces of both inputs on one sphere
(S9e.4b.3c)")`, the refusal the step removes) and a piece whose rim OCCT
divided on import: the capture's `rust_one_sphere_boolean_exists` is false,
and the comparison requires every case `unsupported` until that refusal
is gone; then none of the solid, empty or degenerate cases may stay
`unsupported`.

`--write-bodies` writes the bodies themselves first: the oracle's `write`
blocks of `boolean-one-sphere-bodies.txt` into `rust/fixtures/imported/`.
"""
import json
import sys

import compare_boolean as base
import compare_imported_boolean as imported
import generate_one_sphere_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/pieces.rs'
# The refusal S9e.4b.3c.1 removes.
REFUSAL = 'faces of both inputs on one sphere (S9e.4b.3c)'
FIXTURES = imported.FIXTURES


def rust_one_sphere_boolean_exists():
    """Whether S9e.4b.3c.1's code exists: the pieces' module without the
    refusal of two inputs on one sphere."""
    return KERNEL.exists() and REFUSAL not in KERNEL.read_text()


class OneSphereSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.3c.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/'rust/fixtures/boolean-one-sphere-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-one-sphere-preimplementation'
        self.output = ROOT/'target/boolean-one-sphere-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-one-sphere-divergences.json'
        self.exists_key = 'rust_one_sphere_boolean_exists'
        self.kinds = {c.name: c.kind for c in fixtures.all_cases()}

    def exists(self):
        return rust_one_sphere_boolean_exists()

    def before_code(self):
        return not rust_one_sphere_boolean_exists()

    def must_support(self, case):
        return self.kinds[case.name] != 'unsupported' and rust_one_sphere_boolean_exists()


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


base.make_set = OneSphereSet
base.native_input = native_input
base.case_scale = fixtures.case_size


def write_bodies():
    """`--write-bodies`: OCCT writes every new body's file."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-one-sphere-oracle').resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, _, _ = base.build(prefix, output, base.SOURCE_FILE, 'boolean-oracle',
                                       ('TKBO', 'TKPrim', 'TKShHealing'),
                                       ['TKBO', 'TKPrim', 'TKShHealing']+base.TOOLKITS)
    (FIXTURES/'imported').mkdir(exist_ok=True)
    text = (FIXTURES/fixtures.BODIES_FILE).read_text()
    record = base.run(executable, text, env)
    lines = record['stdout'].split()
    written = [w for w in record['stdout'].splitlines() if w.endswith(' written')]
    wanted = sum(1 for b in fixtures.BODIES.values() if not getattr(b, 'reused', False))
    if record['exit_code'] != 0 or len(written) != wanted or 'failure' in lines:
        raise SystemExit('the bodies were not all written: '+json.dumps(record)[:2000])
    print(f'{len(written)} bodies written')


if __name__ == '__main__':
    if '--write-bodies' in sys.argv:
        write_bodies()
    else:
        base.main()
