#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.4a reference: Booleans of a whole torus against
a polyhedral prism in any relative position
(`generate_torus_boolean_fixtures.py`, `torus_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, a torus's
`torus` row built by `BRepPrimAPI_MakeTorus`) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.4a's
fixtures, capture (`fixtures/occt-boolean-torus-preimplementation`) and
reviews (`fixtures/occt-boolean-torus-divergences.json`). The capture was
taken before the kernel's S9d.4 module
(`rust/kernel/src/solid/boolean/curved/torus.rs`) existed, while the kernel
refused every Boolean with a torus: the capture's `rust_torus_boolean_exists`
is false; the probe reported `unsupported` on every case, which the
comparison requires until the module exists. With the module the kernel's
rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_torus_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/torus.rs'


def rust_torus_boolean_exists():
    """Whether the kernel's S9d.4 module exists."""
    return KERNEL.exists()


class TorusSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.4a has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-torus-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-torus-preimplementation'
        self.output = ROOT/'target/boolean-torus-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-torus-divergences.json'
        self.exists_key = 'rust_torus_boolean_exists'

    def exists(self):
        return rust_torus_boolean_exists()

    def before_code(self):
        return not rust_torus_boolean_exists()


base.make_set = TorusSet

if __name__ == '__main__':
    base.main()
