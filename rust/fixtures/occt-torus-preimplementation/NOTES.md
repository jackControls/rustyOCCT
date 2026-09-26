# Tori (S3): native observations before any kernel torus code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0.
`capture.json` records the Rust revision and that no kernel torus code
existed; the kernel tree had no changes. `inputs.txt` is `torus-cases.txt`
(22 whole tori, v-segments and wedges from `generate_primitive_fixtures.py`),
`oracle.cpp` is `occt_primitive_oracle.cpp` as captured, `native.txt` its
output.

## OCCT source review

* `BRepPrimAPI_MakeTorus(gp_Ax2, R, r, angle1, angle2, angle)` builds a
  `BRepPrim_Torus`: the minor circle about `O + R x` in the half-plane of `x`
  and the axis (about `-y`, its parameter the latitude `v`), latitudes
  `angle1..angle2`, revolved by `angle`. Defaults: `0..2π` and a full turn.
* `BRepPrim_OneAxis`: a closed meridian has no top or bottom face; an open
  one has its ends' circles bounding planar discs. A partial turn
  (`2π - angle > Precision::Angular()`) adds planar start and end faces in
  the meridian half-planes. So the whole torus is 1 vertex, 2 edges (the
  meridian circle at `u = 0` and the latitude circle at `v = angle1`, both
  seams), 1 wire, 1 face; a v-segment 2 vertices, 3 edges, 3 wires, 3 faces;
  a wedge 2 vertices, 3 edges (the meridian circles at both ends, closed at
  `v = angle1`, and the latitude arc between them), 3 wires, 3 faces.
* When the meridian's end lies lower than its start (`isHeightInverted`),
  OneAxis reverses only a wedge's start and end faces (lines 546 and 612).

## Observations

All 22 solids are BRepCheck-valid; counts, faces, edges and vertices equal
`primitive_reference.py`'s, and 21 of them its mass properties (worst
`3.2e-15` of the case's scale). `inner_half` (latitudes `π/2..3π/2`, a full
turn) is inside out: `BRepGProp` gives the reference's volume and inertia
negated, its area unchanged. It is the one reviewed difference
(`../occt-primitive-divergences.json`).
