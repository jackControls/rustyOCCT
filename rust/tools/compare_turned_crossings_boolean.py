#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.4d.2c.2 reference: imported turned bodies
(fillets and rounds beside planes, cylinders and cones; bodies OCCT wrote to
`.brep` files by BRepPrimAPI_MakeRevol) given to Booleans whose partners
cross their smooth edges, the tangent circles where a fillet or a round
meets its neighbour (`generate_turned_crossings_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, read from `rust/fixtures`) and the comparison are
S9e.4a's (`compare_imported_boolean.py`, through `compare_boolean.make_set`)
on S9e.4b.4d.2c.2's fixtures, capture
(`fixtures/occt-boolean-turned-crossings-preimplementation`) and reviews
(`fixtures/occt-boolean-turned-crossings-divergences.json`), the kernel's
enclosures against the reference's totals with S9e.4a's slack of 1e-12
relative. The capture was taken before S9e.4b.4d.2c.2's code existed, while
`curved/graph.rs`'s arrangement found a plane's spiric section of a torus
once per plane face (its comment saying so: "found once for each plane"),
handing a body's second torus its first torus's sections: the capture's
`rust_turned_crossings_boolean_exists` is false. Before the code a case may
evaluate within the reference or fail (a partner across two tori of one
body, an open Boolean of arcs: the kernel's probe failing, listed
`unsupported`); once the code exists none of the solid, empty or degenerate
cases may stay `unsupported`.

`--write-bodies` writes this step's own body first: the oracle's `write`
block of `boolean-turned-crossings-bodies.txt` into `rust/fixtures/imported/`.
"""
import json
import sys

import compare_boolean as base
import compare_imported_boolean as imported
import generate_turned_crossings_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/graph.rs'
# The meetings found once that S9e.4b.4d.2c.2 keys per torus (its comment
# in the source).
ONCE = 'found once for each plane'
FIXTURES = imported.FIXTURES


def rust_turned_crossings_boolean_exists():
    """Whether S9e.4b.4d.2c.2's code exists: the arrangement's spiric
    sections no longer found once for each plane alone."""
    return KERNEL.exists() and ONCE not in KERNEL.read_text()


class TurnedCrossingsSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.4d.2c.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/f'rust/fixtures/{fixtures.PREFIX}-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-turned-crossings-preimplementation'
        self.output = ROOT/'target/boolean-turned-crossings-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-turned-crossings-divergences.json'
        self.exists_key = 'rust_turned_crossings_boolean_exists'
        self.kinds = {c.name: c.kind for c in fixtures.all_cases()}

    def exists(self):
        return rust_turned_crossings_boolean_exists()

    def before_code(self):
        """Never every case refused: before the code a partner across one
        torus evaluates already."""
        return False

    def must_support(self, case):
        return rust_turned_crossings_boolean_exists() and self.kinds[case.name] != 'unsupported'


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


base.make_set = TurnedCrossingsSet
base.native_input = native_input
base.case_scale = fixtures.case_size


def write_bodies():
    """`--write-bodies`: OCCT writes this step's own body's file."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-turned-crossings-oracle').resolve()
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
    if record['exit_code'] != 0 or len(written) != len(fixtures.OWN) or 'failure' in lines:
        raise SystemExit('the bodies were not all written: '+json.dumps(record)[:2000])
    print(f'{len(written)} body written')


if __name__ == '__main__':
    if '--write-bodies' in sys.argv:
        write_bodies()
    else:
        base.main()
