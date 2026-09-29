#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9c.2b.1 reference: Booleans of two cylinders in
turned frames (affine models on the stored axes) with crossing axes, their
quartic section crossing no cap's circle
(`generate_turned_boolean_fixtures.py`, `curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9c.2b.1's fixtures, capture (`fixtures/occt-boolean-turned-
preimplementation`) and reviews (`fixtures/occt-boolean-turned-
divergences.json`). The capture was taken before the kernel's S9c.2b.1
module (`rust/kernel/src/solid/boolean/curved/turned.rs`) existed, while
S9c.2a's code refused these pairs (`OutOfDomain(... (S9c.2))`: two
cylinders not circular in a common measure; the capture's
`rust_turned_boolean_exists` is false; the probe reported `unsupported` on
every case, which the comparison requires until the module exists). With
the module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_turned_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/turned.rs'


def rust_turned_boolean_exists():
    """Whether the kernel's S9c.2b.1 module exists."""
    return KERNEL.exists()


class TurnedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9c.2b.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-turned-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-turned-preimplementation'
        self.output = ROOT/'target/boolean-turned-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-turned-divergences.json'
        self.exists_key = 'rust_turned_boolean_exists'

    def exists(self):
        return rust_turned_boolean_exists()

    def before_code(self):
        return not rust_turned_boolean_exists()


base.make_set = TurnedSet

if __name__ == '__main__':
    base.main()
