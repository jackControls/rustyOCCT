# Sheets, shells, wires and an acorn (S6): native observations before any kernel code accepts them

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code accepted a body without a solid region:
the validator still required every shell to close around a solid (its
Euler and orientation checks), and nothing built, imported or wrote such a
body. `inputs.txt` holds the explicit OCCT constructions of
`generate_brep_fixtures.sheet_models()`; a `result` row makes the oracle
check a free face, a shell, a wire of `W` edges, a free edge or a vertex
instead of a solid of the shells. `oracle.cpp` is
`occt_brep_check_oracle.cpp` as captured, and `native.txt` its four rows per
case: BRepCheck statuses, distinct subshape counts, tolerances with OCCT's
measured gaps and deviations, and `G`: the area (faces, shells) or length
(wires, edges) from `BRepGProp` with its centre.

## The models

Cut from the neutral prisms (`extract`, `wire_of`): free faces (a
rectangle cap, the bottom cap with its reversed frame, a cap with a round
hole, a stadium, a half disc, a partial cylinder, the spline wall of the
bulge), an open box (five faces in one shell), closed shells without a solid
(the box and the cylinder with its seam), wires (the box's bottom edges, a
closed circle, a line-arc-line run) and free edges, and a free vertex. Three
mutations: a free face's vertex moved by `1e-3`, a pcurve of the open box
shifted by `1e-3`, and a wire of two edges sharing no vertex.

## Observations

* Every constructed body is BRepCheck-valid. Counts: a free face has no
  shell (`4 4 1 1 0 0` for the rectangle), the open box one shell and no
  solid (`8 12 5 5 1 0`), the cylinder shell its seam (`2 3 3 3 1 0`), a
  wire one wire and no face, a free edge no wire (`2 1 0 0 0 0`; the closed
  circle `1 1 0 0 0 0`), the vertex `1 0 0 0 0 0`.
* Measures are the exact areas and lengths to rounding: `6` for the
  rectangle, `12 - 0.5625π = 10.2328541323557` for the holed cap, `π` for
  the half cylinder, `10` for the square wire, `3π` for the circle.
* The moved vertex reports `BRepCheck_InvalidPointOnCurve` (1) on it; the
  shifted pcurve `InvalidCurveOnSurface` and `InvalidSameParameterFlag` (8,
  11) on its edge and `UnorientableShape` (27) on the face, as for the
  prism and spline shifts; the disconnected wire `NotConnected` (29) on the
  wire, its two unused edges absent from the analysed shape.
* Every measured deviation and vertex gap is at most `3.7e-16` (the
  cylinder's seam); tolerances are the requested `1e-7`.

`compare_brep.py --family sheet` requires every later run to reproduce
these rows (statuses and counts exactly, measurements as for M5, measures
within `1e-9` relative) and compares the kernel's verdicts, counts,
enclosures and measures with them.
