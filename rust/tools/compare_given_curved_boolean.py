#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.3a reference: a Boolean's result given to
further Booleans where the solids are spheres, cones or frusta and whole
tori besides prisms with arcs, deeper chains, and given results against a
sphere, a cone or a torus (`generate_given_curved_boolean_fixtures.py`,
`chained_curved_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: its `then` rows,
OCCT running each Boolean on the previous one's result's one solid and the
next solid, each prism `BRepPrimAPI_MakePrism` of its profile face, a sphere,
cone or torus its `BRepPrimAPI` primitive) and the comparison are
`compare_boolean.py`'s, run through `compare_boolean.make_set` on S9e.3a's
fixtures, capture (`fixtures/occt-boolean-given-curved-preimplementation`)
and reviews (`fixtures/occt-boolean-given-curved-divergences.json`); the
native rows are the last Boolean's. The capture was taken before the
kernel's S9e.3a module (`rust/kernel/src/solid/boolean/curved/chain.rs`)
existed, while the kernel refused every case (`OutOfDomain`: a result of
solids other than prisms, with procedural edges or against a sphere, cone
or torus, each naming S9e.3): the capture's `rust_given_curved_boolean_exists`
is false; the probe reported `unsupported` on every case, which the
comparison requires until the module exists. With the module the kernel's
rows are compared as S9c.1's.
"""
import compare_boolean as base
import generate_given_curved_boolean_fixtures as fixtures
from identity_reference import native_chained_case

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/chain.rs'


def rust_given_curved_boolean_exists():
    """Whether the kernel's S9e.3a module exists."""
    return KERNEL.exists()


class GivenCurvedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.3a has no spline set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-given-curved-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-given-curved-preimplementation'
        self.output = ROOT/'target/boolean-given-curved-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-given-curved-divergences.json'
        self.exists_key = 'rust_given_curved_boolean_exists'

    def exists(self):
        return rust_given_curved_boolean_exists()

    def before_code(self):
        return not rust_given_curved_boolean_exists()


def native_input():
    """The given cases' native rows."""
    out = []
    for c in base.SET.cases():
        s = c.solid_cases
        (op2, sw2) = c.stages[0]
        more = [(op, s[3+k], sw, None) for k, (op, sw) in enumerate(c.stages[1:])]
        out.append(native_chained_case(s[0], c.op1, s[1], op2, s[2], sw2, None, more))
    return '\n'.join(out)+'\n'


def case_scale(case):
    """The case's size: its solids' reach from the origin."""
    values = [1.0]
    for c in case.solid_cases:
        reach = max(abs(x) for x in c.frame[:3])
        if c.sphere is not None:
            values.append(reach+c.sphere[0])
        elif c.cone is not None:
            values.append(reach+max(c.cone[0], c.cone[1])+c.cone[2])
        elif c.torus is not None:
            values.append(reach+c.torus[0]+c.torus[1])
        else:
            values += [reach+abs(c.start), reach+abs(c.end)]
            for b in c.boundaries:
                if b.circle is not None:
                    values += [reach+abs(b.circle[0])+b.circle[2], reach+abs(b.circle[1])+b.circle[2]]
                else:
                    values += [reach+abs(x) for p in b.points for x in p]
    return max(values)


base.make_set = GivenCurvedSet
base.native_input = native_input
base.case_scale = case_scale

if __name__ == '__main__':
    base.main()
