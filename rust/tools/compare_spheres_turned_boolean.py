#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.2c reference: Booleans of a sphere against a
prism with circles where a frame is turned, a turned cap's circles (of
unequal axes) against a cylinder and a sphere meeting a turned cylinder in a
loop (`generate_spheres_turned_boolean_fixtures.py`,
`spheres_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, each sphere's
`sphere` row built by `BRepPrimAPI_MakeSphere`) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.2c's
fixtures, capture (`fixtures/occt-boolean-spheres-turned-preimplementation`)
and reviews (`fixtures/occt-boolean-spheres-turned-divergences.json`). The
capture was taken before the kernel's S9d.2c module
(`rust/kernel/src/solid/boolean/curved/spheres_turned.rs`) existed, while
the kernel refused a sphere's circle of unequal axes against a cylinder
(`OutOfDomain("a sphere's circle of unequal axes against a cylinder
(S9d.2b)")`) and a sphere meeting a turned cylinder in a loop
(`OutOfDomain("a sphere meeting a turned cylinder in a loop (S9d.2b)")`):
the capture's `rust_spheres_turned_boolean_exists` is false; the probe
reported `unsupported` on every case, which the comparison requires until
the module exists. With the module the kernel's rows are compared as
S9d.1's.
"""
import compare_boolean as base
import generate_spheres_turned_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spheres_turned.rs'


def rust_spheres_turned_boolean_exists():
    """Whether the kernel's S9d.2c module exists."""
    return KERNEL.exists()


class SpheresTurnedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.2c has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spheres-turned-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spheres-turned-preimplementation'
        self.output = ROOT/'target/boolean-spheres-turned-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spheres-turned-divergences.json'
        self.exists_key = 'rust_spheres_turned_boolean_exists'

    def exists(self):
        return rust_spheres_turned_boolean_exists()

    def before_code(self):
        return not rust_spheres_turned_boolean_exists()


base.make_set = SpheresTurnedSet

if __name__ == '__main__':
    base.main()
