# Cones (S3): native observations before any kernel cone code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks (`build_pinned_occt.py`, the Boolean
SDK, which includes TKPrim). `capture.json` records the Rust revision and
that no kernel cone code existed; the kernel tree had no changes.
`inputs.txt` is `primitive-cases.txt` (22 cones from
`generate_primitive_fixtures.py`), `oracle.cpp` is
`occt_primitive_oracle.cpp` as captured, `native.txt` its output.

## OCCT source review

* `BRepPrimAPI_MakeCone(gp_Ax2, R1, R2, H)` builds a `BRepPrim_Cone`, a
  `BRepPrim_OneAxis` revolution of the generatrix from `(R1, 0)` to
  `(R2, H)` about the axis. `myAngle` is `2π` and `myMeridianOffset` 0, so
  the lateral face's seam lies along the frame's x direction (u = 0 and
  u = 2π, `BRepPrim_OneAxis.cxx` lines 426–437), and each end circle runs
  from 0 to `myAngle` from its seam vertex.
* A zero radius end is an apex: one vertex and a degenerated edge whose
  pcurve is the line v = v_apex. A nonzero end is a circle, closed at its
  seam vertex, bounding a planar disc.
* Hence an apex cone has 2 vertices (apex and the base's seam vertex), 3
  edges (base circle, seam, degenerated edge), 2 wires and 2 faces; a
  frustum 2 vertices, 3 edges, 3 wires and 3 faces; one shell, one solid.

## Observations

All 22 solids are BRepCheck-valid. Counts, faces (type, area, centre),
edges (degenerated, closed, length, point at the middle parameter) and
vertices equal `primitive_reference.py`'s, which derives them from the
specification alone. Volume, area, centre and inertia agree to `4.9e-15`
of the case's scale (bound `1e-9`). The degenerated edge reports a length
of about `2.5e-15` from `BRepGProp` (it integrates the pcurve on the
surface), not 0: the bound absorbs it.
