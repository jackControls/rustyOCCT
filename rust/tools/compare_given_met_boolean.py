#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.3b reference: a Boolean's result given to
another Boolean whose solid's faces meet the given result's edges on
meetings of two curved faces (`Rise`, `Meet`, `Toric`) or on a plane's
general section of a cone or of a torus
(`generate_given_met_boolean_fixtures.py`, `chained_curved_boolean_reference.py`
with `given_met_reference.py`'s triple points).

The protocol, the native probe (`occt_boolean_oracle.cpp`: its `then` row,
OCCT running the second Boolean on the first's one solid and the third
solid, each prism `BRepPrimAPI_MakePrism` of its profile face, a sphere,
cone or torus its `BRepPrimAPI` primitive) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9e.3b's
fixtures, capture (`fixtures/occt-boolean-given-met-preimplementation`) and
reviews (`fixtures/occt-boolean-given-met-divergences.json`); the native
rows are the second Boolean's. The capture was taken before the kernel's
S9e.3b module (`rust/kernel/src/solid/boolean/curved/triple.rs`) existed,
while the kernel refused every case (`OutOfDomain`: a given result's meeting
of two curved faces met by another face, naming S9e.3b): the capture's
`rust_given_met_boolean_exists` is false; the probe reported `unsupported`
on every case, which the comparison requires until the module exists. With
the module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import compare_given_curved_boolean as curved
import generate_given_met_boolean_fixtures as fixtures
from identity_reference import native_chained_case

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/triple.rs'


def rust_given_met_boolean_exists():
    """Whether the kernel's S9e.3b module exists."""
    return KERNEL.exists()


class GivenMetSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.3b has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-given-met-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-given-met-preimplementation'
        self.output = ROOT/'target/boolean-given-met-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-given-met-divergences.json'
        self.exists_key = 'rust_given_met_boolean_exists'

    def exists(self):
        return rust_given_met_boolean_exists()

    def before_code(self):
        return not rust_given_met_boolean_exists()


def native_input():
    """The cases' native rows."""
    out = []
    for c in base.SET.cases():
        s = c.solid_cases
        out.append(native_chained_case(s[0], c.op1, s[1], c.op2, s[2], c.swapped, None, []))
    return '\n'.join(out)+'\n'


base.make_set = GivenMetSet
base.native_input = native_input
base.case_scale = curved.case_scale

if __name__ == '__main__':
    base.main()
