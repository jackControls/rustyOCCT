#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.1 reference: a Boolean's result given to
another Boolean, `(A op1 B) op2 C` or, swapped, `C op2 (A op1 B)`, for
prisms of line, arc and circle profiles in any relative position
(`generate_chained_boolean_fixtures.py`, `chained_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: its `then` row,
OCCT running the first Boolean, taking its result's one solid and running
the second with the third prism, each prism `BRepPrimAPI_MakePrism` of its
profile face) and the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9e.1's fixtures, capture
(`fixtures/occt-boolean-chained-preimplementation`) and reviews
(`fixtures/occt-boolean-chained-divergences.json`); the native rows are the
second Boolean's. The capture was taken before the kernel's S9e.1 module
(`rust/kernel/src/solid/boolean/curved/given.rs`) existed, while the kernel
refused a Boolean's result with curved faces as an input
(`OutOfDomain("a Boolean of a solid with curved faces or edges in any
position (S9c)")`): the capture's `rust_chained_boolean_exists` is false;
the probe reported `unsupported` on every case, which the comparison
requires until the module exists. With the module the kernel's rows are
compared as S9c.1's.
"""
import compare_boolean as base
import generate_chained_boolean_fixtures as fixtures
from identity_reference import native_chained_case

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/given.rs'


def rust_chained_boolean_exists():
    """Whether the kernel's S9e.1 module exists."""
    return KERNEL.exists()


class ChainedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.1 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-chained-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-chained-preimplementation'
        self.output = ROOT/'target/boolean-chained-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-chained-divergences.json'
        self.exists_key = 'rust_chained_boolean_exists'

    def exists(self):
        return rust_chained_boolean_exists()

    def before_code(self):
        return not rust_chained_boolean_exists()


def native_input():
    """The chained cases' native rows."""
    return '\n'.join(native_chained_case(c.obj, c.op1, c.tool, c.then, c.third, c.swapped)
                     for c in base.SET.cases())+'\n'


def case_scale(case):
    """The case's size: its three prisms' reach."""
    values = [1.0]
    for c in (case.obj, case.tool, case.third):
        values += [abs(c.start), abs(c.end)]+[abs(x) for x in c.frame[:3]]
        for b in c.boundaries:
            if b.circle is not None:
                values += [abs(b.circle[0])+b.circle[2], abs(b.circle[1])+b.circle[2]]
            else:
                values += [abs(x) for p in b.points for x in p]
    return max(values)


base.make_set = ChainedSet
base.native_input = native_input
base.case_scale = case_scale

if __name__ == '__main__':
    base.main()
