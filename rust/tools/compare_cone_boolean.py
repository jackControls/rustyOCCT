#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.3a reference: Booleans of a cone or frustum
against a polyhedral prism in any relative position
(`generate_cone_boolean_fixtures.py`, `cone_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, a cone's `cone`
row built by `BRepPrimAPI_MakeCone`) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.3a's
fixtures, capture (`fixtures/occt-boolean-cone-preimplementation`) and
reviews (`fixtures/occt-boolean-cone-divergences.json`). The capture was
taken before the kernel's S9d.3 module
(`rust/kernel/src/solid/boolean/curved/cone.rs`) existed, while the kernel
refused every Boolean with a cone (`OutOfDomain("a Boolean of a solid with
curved faces or edges in any position (S9c)")`: the capture's
`rust_cone_boolean_exists` is false; the probe reported `unsupported` on
every case, which the comparison requires until the module exists). With
the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_cone_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/cone.rs'


def rust_cone_boolean_exists():
    """Whether the kernel's S9d.3 module exists."""
    return KERNEL.exists()


class ConeSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.3a has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-cone-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-cone-preimplementation'
        self.output = ROOT/'target/boolean-cone-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-cone-divergences.json'
        self.exists_key = 'rust_cone_boolean_exists'

    def exists(self):
        return rust_cone_boolean_exists()

    def before_code(self):
        return not rust_cone_boolean_exists()


base.make_set = ConeSet

if __name__ == '__main__':
    base.main()
