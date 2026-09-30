#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.4c reference: Booleans of a sphere's cap or
zone against a whole torus, and of a torus v-segment or wedge against a
prism with arcs, a sphere (whole or a cap), a cone or a whole torus, in any
relative position (`generate_torus_parts_boolean_fixtures.py`,
`torus_parts_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: a torus's `torus`
row built by `BRepPrimAPI_MakeTorus(gp_Ax2, R, r, low, high, angle)`, a
sphere's or cap's `sphere` row by `BRepPrimAPI_MakeSphere(gp_Ax2, R, low,
high)`, a cone's by `BRepPrimAPI_MakeCone`, a prism's by
`BRepPrimAPI_MakePrism` of its profile face) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.4c's
fixtures, capture (`fixtures/occt-boolean-torus-parts-preimplementation`)
and reviews (`fixtures/occt-boolean-torus-parts-divergences.json`). The
capture was taken before the kernel's S9d.4c module
(`rust/kernel/src/solid/boolean/curved/torus_parts.rs`) existed, while the
kernel refused a sphere's circle of a surd radius against a torus
(`OutOfDomain("a sphere's circle of a surd radius against a torus
(S9d.4b)")`) and a torus segment or wedge against a curved face
(`OutOfDomain("a torus segment or wedge against a curved face (S9d.4b)")`):
the capture's `rust_torus_parts_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the module
exists. With the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_torus_parts_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/torus_parts.rs'


def rust_torus_parts_boolean_exists():
    """Whether the kernel's S9d.4c module exists."""
    return KERNEL.exists()


class TorusPartsSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.4c has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-torus-parts-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-torus-parts-preimplementation'
        self.output = ROOT/'target/boolean-torus-parts-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-torus-parts-divergences.json'
        self.exists_key = 'rust_torus_parts_boolean_exists'

    def exists(self):
        return rust_torus_parts_boolean_exists()

    def before_code(self):
        return not rust_torus_parts_boolean_exists()


base.make_set = TorusPartsSet

if __name__ == '__main__':
    base.main()
