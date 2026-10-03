#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9f.2a reference: Booleans of spline prisms against
prisms with arc, circle or spline walls whose axes are exactly parallel, in
frames that differ (turned about the axis, rotated about a tilted axis, an
origins' offset that rounds; `generate_spline_parallel_boolean_fixtures.py`,
`curved_boolean_reference.py` with its parallel walls).

The protocol, the native probe (`occt_boolean_oracle.cpp`: each spline path
segment a `Geom_BSplineCurve` edge through its lifted poles, each arc a
circle's, each prism `BRepPrimAPI_MakePrism` of its profile face in its own
frame) and the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9f.2a's fixtures, capture
(`fixtures/occt-boolean-spline-parallel-preimplementation`) and reviews
(`fixtures/occt-boolean-spline-parallel-divergences.json`). The capture was
taken before the kernel's S9f.2a module
(`rust/kernel/src/solid/boolean/curved/spline_parallel.rs`) existed, while
the kernel refused a spline prism against a prism with arcs or splines in
any position (`OutOfDomain("a spline prism against a prism with arcs in any
position (S9f.2)")` or `OutOfDomain("spline walls against spline walls in
any position (S9f.2)")` from the curved engine's `spline_pairs`): the
capture's `rust_spline_parallel_boolean_exists` is false, and until the
module exists the probe must report `unsupported` on every case (a probe
failure or any other row is a failure). With the module the kernel's rows
are compared as S9c.1's.
"""
import compare_boolean as base
import generate_spline_parallel_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spline_parallel.rs'


def rust_spline_parallel_boolean_exists():
    """Whether the kernel's S9f.2a module exists."""
    return KERNEL.exists()


class SplineParallelSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9f.2a is one set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spline-parallel-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spline-parallel-preimplementation'
        self.output = ROOT/'target/boolean-spline-parallel-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spline-parallel-divergences.json'
        self.exists_key = 'rust_spline_parallel_boolean_exists'

    def exists(self):
        return rust_spline_parallel_boolean_exists()

    def before_code(self):
        return not rust_spline_parallel_boolean_exists()


base.make_set = SplineParallelSet

if __name__ == '__main__':
    base.main()
