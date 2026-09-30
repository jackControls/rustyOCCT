#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9d.3c reference: Booleans of a cone against a
cylinder or another cone meeting it in loops (in exact and turned frames,
their curve through infinity where the cones' direction cones cross), of a
turned cone against a sphere in a loop, and of a turned cap against a cone
(`generate_cones_loops_boolean_fixtures.py`, `cones_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`, each cone's
`cone` row built by `BRepPrimAPI_MakeCone`, each sphere's or hemisphere's
`sphere` row by `BRepPrimAPI_MakeSphere`) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9d.3c's
fixtures, capture (`fixtures/occt-boolean-cones-loops-preimplementation`)
and reviews (`fixtures/occt-boolean-cones-loops-divergences.json`). The
capture was taken before the kernel's S9d.3c module
(`rust/kernel/src/solid/boolean/curved/cones_loops.rs`) existed, while the
kernel refused a cone meeting a cylinder or a cone in loops
(`OutOfDomain("a cone meeting a curved face in a loop (S9d.3b.2)")`, a
cylinder along a cone's ruling among them), a turned cone meeting a sphere
in a loop (`OutOfDomain("a sphere meeting a turned cone in a loop
(S9d.3c)")`) and a turned cap's circle against a cone (`OutOfDomain("a
turned cap's circle against a cone (S9d.3c)")`): the capture's
`rust_cones_loops_boolean_exists` is false; the probe reported `unsupported`
on every case, which the comparison requires until the module exists. With
the module the kernel's rows are compared as S9d.1's.
"""
import compare_boolean as base
import generate_cones_loops_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/cones_loops.rs'


def rust_cones_loops_boolean_exists():
    """Whether the kernel's S9d.3c module exists."""
    return KERNEL.exists()


class ConesLoopsSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9d.3c has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-cones-loops-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-cones-loops-preimplementation'
        self.output = ROOT/'target/boolean-cones-loops-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-cones-loops-divergences.json'
        self.exists_key = 'rust_cones_loops_boolean_exists'

    def exists(self):
        return rust_cones_loops_boolean_exists()

    def before_code(self):
        return not rust_cones_loops_boolean_exists()


base.make_set = ConesLoopsSet

if __name__ == '__main__':
    base.main()
