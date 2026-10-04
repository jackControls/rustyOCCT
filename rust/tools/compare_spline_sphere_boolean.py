#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9f.3a reference: Booleans of spline prisms against
spheres and hemispheres in any relative position
(`generate_spline_sphere_boolean_fixtures.py`,
`spline_sphere_boolean_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: each spline path
segment a `Geom_BSplineCurve` edge through its lifted poles, each prism
`BRepPrimAPI_MakePrism` of its profile face in its own frame, each sphere
`BRepPrimAPI_MakeSphere` on its frame between its latitudes) and the
comparison are `compare_boolean.py`'s, run through `compare_boolean.make_set`
on S9f.3a's fixtures, capture
(`fixtures/occt-boolean-spline-sphere-preimplementation`) and reviews
(`fixtures/occt-boolean-spline-sphere-divergences.json`). The capture was
taken before the kernel's S9f.3a module
(`rust/kernel/src/solid/boolean/curved/spline_sphere.rs`) existed, while
the kernel refused a spline prism against a sphere
(`OutOfDomain("a spline prism against a sphere or a cone (S9f.3)")` from
the curved engine's `spline_pairs`): the capture's
`rust_spline_sphere_boolean_exists` is false, and until the module exists
the probe must report `unsupported` on every case. With the module the
kernel's rows are compared as S9c.1's, and none may be `unsupported`
(`rust_unsupported_after_its_code`).
"""
import compare_boolean as base
import generate_spline_sphere_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/boolean/curved/spline_sphere.rs'


def rust_spline_sphere_boolean_exists():
    """Whether the kernel's S9f.3a module exists."""
    return KERNEL.exists()


class SplineSphereSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9f.3a is one set'
        super().__init__(False)
        self.cases = fixtures.cases
        self.expected = ROOT/'rust/fixtures/boolean-spline-sphere-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-spline-sphere-preimplementation'
        self.output = ROOT/'target/boolean-spline-sphere-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-spline-sphere-divergences.json'
        self.exists_key = 'rust_spline_sphere_boolean_exists'

    def exists(self):
        return rust_spline_sphere_boolean_exists()

    def before_code(self):
        return not rust_spline_sphere_boolean_exists()

    def must_support(self, case):
        return rust_spline_sphere_boolean_exists()


base.make_set = SplineSphereSet

if __name__ == '__main__':
    base.main()
