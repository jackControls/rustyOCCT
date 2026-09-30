#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9f.1 reference: Booleans of spline prisms (lines,
arcs and S8b's nonrational splines in one frame) against polyhedral prisms
in any relative position, same-axis pairs whose offset rounds included
(`generate_spline_any_boolean_fixtures.py`, `curved_boolean_reference.py`
with its spline walls).

The protocol, the native probe (`occt_boolean_oracle.cpp`: each spline path
segment a `Geom_BSplineCurve` edge through its lifted poles, each prism
`BRepPrimAPI_MakePrism` of its profile face in its own frame) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9f.1's fixtures, capture (`fixtures/occt-boolean-spline-any-
preimplementation`) and reviews (`fixtures/occt-boolean-spline-any-
divergences.json`). The capture was taken before the kernel's S9f.1 module
(`rust/kernel/src/solid/boolean/curved/spline_walls.rs`) existed, while the
kernel refused a spline prism in a Boolean in any position
(`OutOfDomain("a Boolean of a solid with curved faces or edges in any
position (S9c)")` from the polyhedral engine, or `OutOfDomain("a spline
profile in a Boolean of prisms in any position (S9c)")` from the curved
engine's model where the profile also holds an arc): the capture's
`rust_spline_any_boolean_exists` is false, and until the module exists the
probe must report `unsupported` on every case (a probe failure or any other
row is a failure). With the module the kernel's rows are compared as
S9c.1's.
"""
import compare_boolean as base
import generate_spline_any_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spline_walls.rs'


def rust_spline_any_boolean_exists():
    """Whether the kernel's S9f.1 module exists."""
    return KERNEL.exists()


class SplineAnySet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9f.1 is one set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spline-any-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spline-any-preimplementation'
        self.output = ROOT/'target/boolean-spline-any-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spline-any-divergences.json'
        self.exists_key = 'rust_spline_any_boolean_exists'

    def exists(self):
        return rust_spline_any_boolean_exists()

    def before_code(self):
        return not rust_spline_any_boolean_exists()


base.make_set = SplineAnySet

if __name__ == '__main__':
    base.main()
