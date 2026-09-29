#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.4b.2 reference: Booleans of a whole torus against
a prism with arcs, a sphere, a cone or another whole torus in any relative
position (`generate_torus_curved_boolean_fixtures.py`,
`torus_curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: a torus's `torus`
row built by `BRepPrimAPI_MakeTorus(gp_Ax2, R, r)`, a sphere's by
`BRepPrimAPI_MakeSphere`, a cone's by `BRepPrimAPI_MakeCone`, a prism's by
`BRepPrimAPI_MakePrism` of its profile face, either input or both) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9d.4b.2's fixtures, capture (`fixtures/occt-boolean-torus-curved-
preimplementation`) and reviews (`fixtures/occt-boolean-torus-curved-
divergences.json`). The capture was taken before the kernel's S9d.4b.2
module (`rust/kernel/src/solid/boolean/curved/torus_curved.rs`) existed,
while the kernel refused every Boolean of a torus against a curved face
(`OutOfDomain("a torus against a curved face (S9d.4b)")`): the capture's
`rust_torus_curved_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the module
exists. With the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_torus_curved_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/torus_curved.rs'


def rust_torus_curved_boolean_exists():
    """Whether the kernel's S9d.4b.2 module exists."""
    return KERNEL.exists()


class TorusCurvedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.4b.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-torus-curved-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-torus-curved-preimplementation'
        self.output = ROOT/'target/boolean-torus-curved-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-torus-curved-divergences.json'
        self.exists_key = 'rust_torus_curved_boolean_exists'

    def exists(self):
        return rust_torus_curved_boolean_exists()

    def before_code(self):
        return not rust_torus_curved_boolean_exists()


base.make_set = TorusCurvedSet

if __name__ == '__main__':
    base.main()
