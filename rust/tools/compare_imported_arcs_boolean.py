#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.1 reference: imported prisms whose arcs'
ends round off their circles in their caps' frames (bodies OCCT wrote to
`.brep` files) given to Booleans
(`generate_imported_arcs_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, the file's one solid read by `BRepTools::Read` from
`rust/fixtures`) and the comparison are S9e.4a's
(`compare_imported_boolean.py`, through `compare_boolean.make_set`) on
S9e.4b.1's fixtures, capture
(`fixtures/occt-boolean-imported-arcs-preimplementation`) and reviews
(`fixtures/occt-boolean-imported-arcs-divergences.json`), the kernel's
enclosures against the reference's totals with S9e.4a's slack of 1e-12
relative. The capture was taken before the kernel's S9e.4b.1 module
(`rust/kernel/src/solid/boolean/curved/snapped.rs`) existed, while S9e.4a
refused every such prism ("an imported prism's arc whose ends round off its
circle in its cap's frame"): the capture's
`rust_imported_arcs_boolean_exists` is false, and the comparison requires
every case `unsupported` until the module exists.

`--write-bodies` writes the bodies themselves first: the oracle's `write`
blocks of `boolean-imported-arcs-bodies.txt` (`MakePrism` of each profile,
`BRepTools::Write` in format version 1) into `rust/fixtures/imported/`.
"""
import json
import sys

import compare_boolean as base
import compare_imported_boolean as imported
import generate_imported_arcs_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/snapped.rs'
FIXTURES = imported.FIXTURES


def rust_imported_arcs_boolean_exists():
    """Whether the kernel's S9e.4b.1 module exists."""
    return KERNEL.exists()


class ArcsSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/'rust/fixtures/boolean-imported-arcs-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-imported-arcs-preimplementation'
        self.output = ROOT/'target/boolean-imported-arcs-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-imported-arcs-divergences.json'
        self.exists_key = 'rust_imported_arcs_boolean_exists'

    def exists(self):
        return rust_imported_arcs_boolean_exists()

    def before_code(self):
        return not rust_imported_arcs_boolean_exists()


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


# S9e.4a's comparison (its case size, slack and oracle environment, set on
# `compare_boolean` by importing it) on this step's set.
base.make_set = ArcsSet
base.native_input = native_input


def write_bodies():
    """`--write-bodies`: OCCT writes every body's file."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-imported-arcs-oracle').resolve()
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
    if record['exit_code'] != 0 or len(written) != len(fixtures.BODIES) or 'failure' in lines:
        raise SystemExit('the bodies were not all written: '+json.dumps(record)[:2000])
    print(f'{len(written)} bodies written')


if __name__ == '__main__':
    if '--write-bodies' in sys.argv:
        write_bodies()
    else:
        base.main()
