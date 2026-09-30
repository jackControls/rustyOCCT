#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.2 reference: a Boolean's result given to
another Boolean where the given solid is a stack with an arc wall, an S9b.1
result of line prisms given with arcs, or one solid of a first result of
several (`generate_given_boolean_fixtures.py`, `chained_boolean_reference.py`
with its selector).

The protocol, the native probe (`occt_boolean_oracle.cpp`: its `then` row,
OCCT running the first Boolean, taking its result's one solid, or the one
`BRepClass3d_SolidClassifier` finds the pick point inside, and running the
second with the third prism, each prism `BRepPrimAPI_MakePrism` of its
profile face) and the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9e.2's fixtures, capture
(`fixtures/occt-boolean-given-preimplementation`) and reviews
(`fixtures/occt-boolean-given-divergences.json`); the native rows are the
second Boolean's. The capture was taken before the kernel's S9e.2 module
(`rust/kernel/src/solid/boolean/curved/matched.rs`) existed, while the
kernel refused every case (`OutOfDomain`: a stack or an S9b.1 result with
arcs refused by `polyhedra.rs`, a result of several solids by
`curved/given.rs`, each naming S9e.2): the capture's
`rust_given_boolean_exists` is false; the probe reported `unsupported` on
every case, which the comparison requires until the module exists. With the
module the kernel's rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_given_boolean_fixtures as fixtures
from identity_reference import native_chained_case

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/matched.rs'


def rust_given_boolean_exists():
    """Whether the kernel's S9e.2 module exists."""
    return KERNEL.exists()


class GivenSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.2 has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-given-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-given-preimplementation'
        self.output = ROOT/'target/boolean-given-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-given-divergences.json'
        self.exists_key = 'rust_given_boolean_exists'

    def exists(self):
        return rust_given_boolean_exists()

    def before_code(self):
        return not rust_given_boolean_exists()


def native_input():
    """The given cases' native rows."""
    return '\n'.join(native_chained_case(c.obj, c.op1, c.tool, c.then, c.third, c.swapped, c.pick)
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


base.make_set = GivenSet
base.native_input = native_input
base.case_scale = case_scale

if __name__ == '__main__':
    base.main()
