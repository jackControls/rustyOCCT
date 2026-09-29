#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9c.2b.2 reference: Booleans of two cylinders whose
quartic section crosses a cap's circle within both faces, in exact and in
turned frames (`generate_capped_boolean_fixtures.py`,
`curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9c.2b.2's fixtures, capture (`fixtures/occt-boolean-capped-
preimplementation`) and reviews (`fixtures/occt-boolean-capped-
divergences.json`). The capture was taken before the kernel's S9c.2b.2
module (`rust/kernel/src/solid/boolean/curved/algebraic.rs`) existed, while
S9c.2a's and S9c.2b.1's code refused these pairs (`OutOfDomain(... crossing
a cap's circle ... (S9c.2b))` in exact frames, `(S9c.2b.2)` in turned ones;
the capture's `rust_capped_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the module
exists). With the module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_capped_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/algebraic.rs'


def rust_capped_boolean_exists():
    """Whether the kernel's S9c.2b.2 module exists."""
    return KERNEL.exists()


class CappedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9c.2b.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-capped-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-capped-preimplementation'
        self.output = ROOT/'target/boolean-capped-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-capped-divergences.json'
        self.exists_key = 'rust_capped_boolean_exists'

    def exists(self):
        return rust_capped_boolean_exists()

    def before_code(self):
        return not rust_capped_boolean_exists()


base.make_set = CappedSet

if __name__ == '__main__':
    base.main()
