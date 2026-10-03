#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9f.2b reference: Booleans of spline prisms against
prisms with arc or circle walls whose axes cross
(`generate_spline_crossing_boolean_fixtures.py`, `curved_boolean_reference.py`
with its crossing walls).

The protocol, the native probe (`occt_boolean_oracle.cpp`: each spline path
segment a `Geom_BSplineCurve` edge through its lifted poles, each arc a
circle's, each prism `BRepPrimAPI_MakePrism` of its profile face in its own
frame) and the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9f.2b's fixtures, capture
(`fixtures/occt-boolean-spline-crossing-preimplementation`) and reviews
(`fixtures/occt-boolean-spline-crossing-divergences.json`). The capture was
taken before the kernel's S9f.2b module
(`rust/kernel/src/solid/boolean/curved/spline_crossing.rs`) existed, while
the kernel refused a spline prism against a prism with arcs on crossing
axes (`OutOfDomain("a spline prism against a prism with arcs on crossing
axes (S9f.2b)")` from the curved engine's `spline_pairs`): the capture's
`rust_spline_crossing_boolean_exists` is false, and until the module
exists the probe must report `unsupported` on every case (a probe failure
or any other row is a failure). With the module the kernel's rows are
compared as S9c.1's; S9f.2b.2's cases (loops) may stay `unsupported`.
"""
import compare_boolean as base
import generate_spline_crossing_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spline_crossing.rs'


def rust_spline_crossing_boolean_exists():
    """Whether the kernel's S9f.2b module exists."""
    return KERNEL.exists()


class SplineCrossingSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9f.2b is one set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spline-crossing-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spline-crossing-preimplementation'
        self.output = ROOT/'target/boolean-spline-crossing-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spline-crossing-divergences.json'
        self.exists_key = 'rust_spline_crossing_boolean_exists'

    def exists(self):
        return rust_spline_crossing_boolean_exists()

    def before_code(self):
        return not rust_spline_crossing_boolean_exists()


base.make_set = SplineCrossingSet

if __name__ == '__main__':
    base.main()
