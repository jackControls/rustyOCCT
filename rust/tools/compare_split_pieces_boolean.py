#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4b.3b reference: the kernel's own plane pieces
(S8's split pieces, `Clipped` and `Half`) given to Booleans with curved faces
(`generate_split_pieces_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: a split piece its
solid's rows and a `split` row, the solid common
`BRepPrimAPI_MakeHalfSpace` of the plane on the kept side) and the
comparison are S9e.4a's (`compare_imported_boolean.py`, through
`compare_boolean.make_set`, its slack of 1e-12 relative: the kernel's plane
is the split frame's stored normal where the reference's half-space box
holds its axes' `x * y`) on S9e.4b.3b's fixtures, capture
(`fixtures/occt-boolean-split-pieces-preimplementation`) and reviews
(`fixtures/occt-boolean-split-pieces-divergences.json`). The capture was
taken before the kernel's S9e.4b.3b module
(`rust/kernel/src/solid/boolean/curved/splits.rs`) existed, while the kernel
refused every split piece against curved faces (`OutOfDomain`, S9e.4's):
the capture's `rust_split_pieces_boolean_exists` is false, and the
comparison requires every case `unsupported` until the module exists.
"""
import compare_boolean as base
import compare_imported_boolean  # noqa: F401 (S9e.4a's slack and oracle environment)
import generate_split_pieces_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/splits.rs'


def rust_split_pieces_boolean_exists():
    """Whether the kernel's S9e.4b.3b module exists."""
    return KERNEL.exists()


class SplitPiecesSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4b.3b has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/'rust/fixtures/boolean-split-pieces-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-split-pieces-preimplementation'
        self.output = ROOT/'target/boolean-split-pieces-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-split-pieces-divergences.json'
        self.exists_key = 'rust_split_pieces_boolean_exists'

    def exists(self):
        return rust_split_pieces_boolean_exists()

    def before_code(self):
        return not rust_split_pieces_boolean_exists()


def native_input():
    """The cases' native rows."""
    return fixtures.native_input()


base.make_set = SplitPiecesSet
base.native_input = native_input
base.case_scale = fixtures.case_size

if __name__ == '__main__':
    base.main()
