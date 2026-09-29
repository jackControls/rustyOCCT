#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9c.2 reference: Booleans of prisms whose
cylindrical walls meet in quartics (perpendicular cylinders of any radii
and offset in exact frames, oblique, skew and parallel ones in turned
frames; `generate_procedural_boolean_fixtures.py`,
`curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9c.2's fixtures, capture (`fixtures/occt-boolean-procedural-
preimplementation`) and reviews (`fixtures/occt-boolean-procedural-
divergences.json`, none so far). The capture was taken before the kernel's
S9c.2 module (`rust/kernel/src/solid/boolean/curved/procedural.rs`) existed,
while S9c.1's code refused these pairs (`OutOfDomain(... (S9c.2))`: the
capture's `rust_procedural_boolean_exists` is false; the probe reported
`unsupported` on every case, which the comparison requires until the module
exists). With the module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_procedural_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/procedural.rs'


def rust_procedural_boolean_exists():
    """Whether the kernel's S9c.2 module exists."""
    return KERNEL.exists()


class ProceduralSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9c.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-procedural-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-procedural-preimplementation'
        self.output = ROOT/'target/boolean-procedural-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-procedural-divergences.json'
        self.exists_key = 'rust_procedural_boolean_exists'

    def exists(self):
        return rust_procedural_boolean_exists()

    def before_code(self):
        return not rust_procedural_boolean_exists()


base.make_set = ProceduralSet

if __name__ == '__main__':
    base.main()
