#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.1 reference: Booleans of a sphere (whole, a cap
or a zone) against a polyhedral prism in any relative position
(`generate_sphere_boolean_fixtures.py`, `sphere_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, a sphere's
`sphere` row built by `BRepPrimAPI_MakeSphere`) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.1's
fixtures, capture (`fixtures/occt-boolean-sphere-preimplementation`) and
reviews (`fixtures/occt-boolean-sphere-divergences.json`). The capture was
taken before the kernel's S9d.1 module
(`rust/kernel/src/solid/boolean/curved/sphere.rs`) existed, while the
kernel refused every Boolean with a sphere (`OutOfDomain("a Boolean of a
solid with curved faces or edges in any position (S9c)")`: the capture's
`rust_sphere_boolean_exists` is false; the probe reported `unsupported` on
every case, which the comparison requires until the module exists). With
the module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_sphere_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/sphere.rs'


def rust_sphere_boolean_exists():
    """Whether the kernel's S9d.1 module exists."""
    return KERNEL.exists()


class SphereSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-sphere-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-sphere-preimplementation'
        self.output = ROOT/'target/boolean-sphere-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-sphere-divergences.json'
        self.exists_key = 'rust_sphere_boolean_exists'

    def exists(self):
        return rust_sphere_boolean_exists()

    def before_code(self):
        return not rust_sphere_boolean_exists()


base.make_set = SphereSet

if __name__ == '__main__':
    base.main()
