#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.4c.2a reference: imported bodies of one
cylinder or cone face and plane faces whose pockets hold pockets of their
own (S9e.4b.3c.3b's tooth, a post standing in a pocket, a U island whose
notch is a pocket again; bodies OCCT wrote to `.brep` files) given to
Booleans (`generate_deep_trees_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, read from `rust/fixtures`) and the comparison are
S9e.4a's (`compare_imported_boolean.py`, through `compare_boolean.make_set`)
on S9e.4b.4c.2a's fixtures, capture
(`fixtures/occt-boolean-deep-trees-preimplementation`) and reviews
(`fixtures/occt-boolean-deep-trees-divergences.json`), the kernel's
enclosures against the reference's totals with S9e.4a's slack of 1e-12
relative. The capture was taken before S9e.4b.4c.2a's code existed, while
the import refused such a body (`imported.rs`'s `OutOfDomain("an imported
plane piece other than a Boolean tree of its primitive and its planes'
hulls (S9e.4b.4)")`, the refusal the step removes, its text replaced by one
naming the limits on depth): the capture's `rust_deep_trees_boolean_exists`
is false, and the comparison requires every case `unsupported` until that
refusal is gone, but for the declared `degenerate` cases, whose refusal
(S9's, against the partner) may come first once the body imports; then none
of the solid, empty or degenerate cases may stay `unsupported`.

`--write-bodies` writes this step's own bodies first: the oracle's `write`
blocks of `boolean-deep-trees-bodies.txt` into `rust/fixtures/imported/`
(the tooth is S9e.4b.3c.3b's, the ball S9e.4a's).
"""
import json
import sys

import compare_boolean as base
import compare_imported_boolean as imported
import generate_deep_trees_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/imported.rs'
# The refusal S9e.4b.4c.2a removes (its text, on one line in the source).
REFUSAL = 'an imported plane piece other than a Boolean tree of its primitive and its planes\' hulls \\\n         (S9e.4b.4)'
FIXTURES = imported.FIXTURES


def rust_deep_trees_boolean_exists():
    """Whether S9e.4b.4c.2a's code exists: the import without the refusal of
    a piece other than a Boolean tree of its primitive and its planes'
    hulls, pockets of one level."""
    return KERNEL.exists() and REFUSAL not in KERNEL.read_text()


class DeepTreesSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.4c.2a has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/f'rust/fixtures/{fixtures.PREFIX}-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-deep-trees-preimplementation'
        self.output = ROOT/'target/boolean-deep-trees-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-deep-trees-divergences.json'
        self.exists_key = 'rust_deep_trees_boolean_exists'
        self.kinds = {c.name: c.kind for c in fixtures.all_cases()}

    def exists(self):
        return rust_deep_trees_boolean_exists()

    def before_code(self):
        return not rust_deep_trees_boolean_exists()

    def refused_before_code(self, case, rows):
        """Before the code, every case `unsupported` (each body refused on
        import as S9e.4b.4's)."""
        return rows == [['unsupported']]

    def must_support(self, case):
        return rust_deep_trees_boolean_exists() and self.kinds[case.name] != 'unsupported'


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


base.make_set = DeepTreesSet
base.native_input = native_input
base.case_scale = fixtures.case_size


def write_bodies():
    """`--write-bodies`: OCCT writes this step's bodies' files."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-deep-trees-oracle').resolve()
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
    print(f'{len(written)} bodies written')


if __name__ == '__main__':
    if '--write-bodies' in sys.argv:
        write_bodies()
    else:
        base.main()
