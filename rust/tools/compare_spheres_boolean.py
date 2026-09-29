#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.2 reference: Booleans of a sphere (whole or a
cap) against a prism with arcs and circles (cylindrical walls), and of two
spheres (`generate_spheres_boolean_fixtures.py`,
`spheres_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, each sphere's
`sphere` row built by `BRepPrimAPI_MakeSphere`, either input or both) and
the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9d.2's fixtures, capture
(`fixtures/occt-boolean-spheres-preimplementation`) and reviews
(`fixtures/occt-boolean-spheres-divergences.json`). The capture was taken
before the kernel's S9d.2 module
(`rust/kernel/src/solid/boolean/curved/spheres.rs`) existed, while the
kernel refused a sphere against a cylinder (`OutOfDomain("a sphere against
a cylinder or a sphere (S9d.2)")`) and two spheres (`OutOfDomain("a
Boolean of a solid with curved faces or edges in any position (S9c)")`):
the capture's `rust_spheres_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the
module exists. With the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_spheres_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spheres.rs'


def rust_spheres_boolean_exists():
    """Whether the kernel's S9d.2 module exists."""
    return KERNEL.exists()


class SpheresSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spheres-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spheres-preimplementation'
        self.output = ROOT/'target/boolean-spheres-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spheres-divergences.json'
        self.exists_key = 'rust_spheres_boolean_exists'

    def exists(self):
        return rust_spheres_boolean_exists()

    def before_code(self):
        return not rust_spheres_boolean_exists()


base.make_set = SpheresSet

if __name__ == '__main__':
    base.main()
