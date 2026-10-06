#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.4b.2b.1 reference: imported bodies of several
sphere, cylinder and cone faces whose other plane faces are a primitive's
flat, a prism of one cap or a pocket (the `boolean` fuzz target's cone on a
ball's half, a hexagon under a stadium, a stadium boss on a box, a plate
with a hole less its moved copy's band; bodies OCCT wrote to `.brep` files)
given to Booleans (`generate_plane_parts_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, read from `rust/fixtures`) and the comparison are
S9e.4a's (`compare_imported_boolean.py`, through `compare_boolean.make_set`)
on S9e.4b.4b.2b.1's fixtures, capture
(`fixtures/occt-boolean-plane-parts-preimplementation`) and reviews
(`fixtures/occt-boolean-plane-parts-divergences.json`), the kernel's
enclosures against the reference's totals with S9e.4a's slack of 1e-12
relative. The capture was taken before S9e.4b.4b.2b.1's code existed, while
the import refused such a body as S9e.4b.4b.2a did (`imported.rs`'s
`OutOfDomain("an imported body of several primitives with plane faces other
than their ends or a prism's (S9e.4b.4b.2b)")`, whose text the step
changes): the capture's `rust_plane_parts_boolean_exists` is false, and the
comparison requires every case `unsupported` until that refusal is gone (the
cake's `refused`: S9e.4b.4b.2a took its boss's arcs for primitives tangent to
its line walls, `Degenerate("an imported body of several primitives whose
faces are tangent along an edge")`); then none of the solid, empty or
degenerate cases may stay `unsupported`.

`--write-bodies` writes the bodies themselves first: the oracle's `write`
blocks of `boolean-plane-parts-bodies.txt` into `rust/fixtures/imported/`.
"""
import json
import sys

import compare_boolean as base
import compare_imported_boolean as imported
import generate_plane_parts_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/imported.rs'
# The refusal S9e.4b.4b.2b.1 removes (its text's end in the source).
REFUSAL = 'prism\'s (S9e.4b.4b.2b)"'
FIXTURES = imported.FIXTURES


def rust_plane_parts_boolean_exists():
    """Whether S9e.4b.4b.2b.1's code exists: the import without S9e.4b.4b.2a's
    refusal of plane faces other than the primitives' ends or a prism's."""
    return KERNEL.exists() and REFUSAL not in KERNEL.read_text()


class PlanePartsSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.4b.2b.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/f'rust/fixtures/{fixtures.PREFIX}-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-plane-parts-preimplementation'
        self.output = ROOT/'target/boolean-plane-parts-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-plane-parts-divergences.json'
        self.exists_key = 'rust_plane_parts_boolean_exists'
        self.kinds = {c.name: c.kind for c in fixtures.all_cases()}

    def exists(self):
        return rust_plane_parts_boolean_exists()

    def before_code(self):
        return not rust_plane_parts_boolean_exists()

    def refused_before_code(self, case, rows):
        """Before the code, every case `unsupported` (every body refused on
        import as S9e.4b.4b.2a refused it; the cake's refusal counted so by
        `rust_rows`)."""
        return rows == [['unsupported']]

    def must_support(self, case):
        return rust_plane_parts_boolean_exists() and self.kinds[case.name] != 'unsupported'


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


base.make_set = PlanePartsSet
_rust_rows = base.rust_rows


def rust_rows():
    """The kernel's probe's rows; before the code, a case of the cake
    `refused` (S9e.4b.4b.2a's reading of its boss's arcs as primitives tangent
    to its walls, the body's own tangency) counted `unsupported`."""
    rows, failed = _rust_rows()
    if not rust_plane_parts_boolean_exists():
        for case in fixtures.all_cases():
            cake = any(s == ('imported', 'cake') for s in case.items)
            if cake and rows.get(case.name) == [['refused']]:
                rows[case.name] = [['unsupported']]
    return rows, failed


base.rust_rows = rust_rows
base.native_input = native_input
base.case_scale = fixtures.case_size


def write_bodies():
    """`--write-bodies`: OCCT writes every body's file."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-plane-parts-oracle').resolve()
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
