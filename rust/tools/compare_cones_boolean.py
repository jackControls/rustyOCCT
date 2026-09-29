#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.3b reference: Booleans of a cone or frustum
against a prism with arcs and circles (cylindrical walls), a sphere or cap,
and another cone or frustum (`generate_cones_boolean_fixtures.py`,
`cones_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, each cone's
`cone` row built by `BRepPrimAPI_MakeCone` and each sphere's `sphere` row by
`BRepPrimAPI_MakeSphere`, either input or both) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.3b's
fixtures, capture (`fixtures/occt-boolean-cones-preimplementation`) and
reviews (`fixtures/occt-boolean-cones-divergences.json`). The capture was
taken before the kernel's S9d.3b module
(`rust/kernel/src/solid/boolean/curved/cones.rs`) existed, while the kernel
refused a cone against a prism with arcs, a sphere or a cone
(`OutOfDomain("a cone against a prism with arcs, a sphere or a cone
(S9d.3b)")`): the capture's `rust_cones_boolean_exists` is false; the probe
reported `unsupported` on every case, which the comparison requires until
the module exists. With the module the kernel's rows are compared as
S9d.1's.
"""
import compare_boolean as base
import generate_cones_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/cones.rs'


def rust_cones_boolean_exists():
    """Whether the kernel's S9d.3b module exists."""
    return KERNEL.exists()


class ConesSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.3b has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-cones-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-cones-preimplementation'
        self.output = ROOT/'target/boolean-cones-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-cones-divergences.json'
        self.exists_key = 'rust_cones_boolean_exists'

    def exists(self):
        return rust_cones_boolean_exists()

    def before_code(self):
        return not rust_cones_boolean_exists()


base.make_set = ConesSet

if __name__ == '__main__':
    base.main()
