#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9c.1 reference: Booleans of prisms of line, arc
and circle profiles in any relative position
(`generate_curved_boolean_fixtures.py`, `curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: arcs and circles
as `gp_Circ` edges about the profile plane's normal) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9c.1's
fixtures, capture (`fixtures/occt-boolean-curved-preimplementation`) and
reviews (`fixtures/occt-boolean-curved-divergences.json`). The capture was
taken before the kernel's S9c.1 module
(`rust/kernel/src/solid/boolean/curved/`) existed, while it refused arcs in
frames with different axes (the capture's `rust_curved_boolean_exists` is
false; the probe reported `unsupported` on every case). With the module the
kernel's rows are compared as S9a's and S9b's; the cases it leaves to S9c.2
(cylinders meeting in quartics) are `unsupported` rows.
"""
import compare_boolean as base
import generate_curved_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/mod.rs'


def rust_curved_boolean_exists():
    """Whether the kernel's S9c.1 module exists."""
    return KERNEL.exists()


class CurvedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9c.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-curved-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-curved-preimplementation'
        self.output = ROOT/'target/boolean-curved-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-curved-divergences.json'
        self.exists_key = 'rust_curved_boolean_exists'

    def exists(self):
        return rust_curved_boolean_exists()

    def before_code(self):
        return not rust_curved_boolean_exists()


base.make_set = CurvedSet

if __name__ == '__main__':
    base.main()
