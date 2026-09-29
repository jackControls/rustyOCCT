#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.4b.1 reference: Booleans of a torus v-segment or
wedge against a polyhedral prism in any relative position
(`generate_torus_segment_boolean_fixtures.py`,
`torus_segment_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, a part's `torus`
row with its latitudes and turn built by `BRepPrimAPI_MakeTorus(gp_Ax2, R,
r, a1, a2, angle)`, an inside-out v-segment reversed) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on
S9d.4b.1's fixtures, capture (`fixtures/occt-boolean-torus-segment-
preimplementation`) and reviews (`fixtures/occt-boolean-torus-segment-
divergences.json`). The capture was taken before the kernel's S9d.4b.1
module (`rust/kernel/src/solid/boolean/curved/torus_segment.rs`) existed,
while the kernel refused every Boolean of a torus segment or wedge
(`OutOfDomain("a Boolean of a torus segment or wedge (S9d.4b)")`): the
capture's `rust_torus_segment_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the module
exists. With the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_torus_segment_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/torus_segment.rs'


def rust_torus_segment_boolean_exists():
    """Whether the kernel's S9d.4b.1 module exists."""
    return KERNEL.exists()


class TorusSegmentSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.4b.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-torus-segment-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-torus-segment-preimplementation'
        self.output = ROOT/'target/boolean-torus-segment-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-torus-segment-divergences.json'
        self.exists_key = 'rust_torus_segment_boolean_exists'

    def exists(self):
        return rust_torus_segment_boolean_exists()

    def before_code(self):
        return not rust_torus_segment_boolean_exists()


base.make_set = TorusSegmentSet

if __name__ == '__main__':
    base.main()
