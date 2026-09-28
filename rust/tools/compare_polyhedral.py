#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9b reference: Booleans of polyhedral prisms in any
relative position (`generate_polyhedral_fixtures.py`,
`polyhedral_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`) and the
comparison are `compare_boolean.py`'s, run on S9b's fixtures, capture
(`fixtures/occt-boolean-polyhedra-preimplementation`: taken before the
kernel's S9b module `rust/kernel/src/solid/boolean/polyhedra.rs` exists) and
reviews (`fixtures/occt-boolean-polyhedra-divergences.json`), through
`compare_boolean.make_set`.
"""
import compare_boolean as base
import generate_polyhedral_fixtures as fixtures

ROOT = base.ROOT
base.KERNEL_FILE = ROOT/'rust/kernel/src/solid/boolean/polyhedra.rs'


class PolyhedraSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9b has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-polyhedra-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-polyhedra-preimplementation'
        self.output = ROOT/'target/boolean-polyhedra-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-polyhedra-divergences.json'


base.make_set = PolyhedraSet

if __name__ == '__main__':
    base.main()
