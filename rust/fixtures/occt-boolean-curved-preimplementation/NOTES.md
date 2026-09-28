# Booleans of prisms with arcs in any position (S9c.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_curved_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9c.1's fixtures). `capture.json` records the Rust revision and
that the kernel's S9c code did not exist (`rust_curved_boolean_exists`
false: `rust/kernel/src/solid/boolean/polyhedra.rs` still refuses a
Boolean of prisms with arcs in frames with different axes, `OutOfDomain(...
(S9c))`, and the kernel's probe reports every case `unsupported`).
`inputs.txt` holds the 44 cases of `boolean-curved-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame (arcs
and circles as `gp_Circ` edges about the plane's normal), the two frames'
axes different (`XY`, `SIDE`, `DOWN`, `TURN`, `TILT`, `TILT2`, `TILTX`,
`LEAN`, `R125`).

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 44 (`curved_boolean_reference.py`: the
  regularized Boolean's maximal connected regions), the 4 empty results
  included (a pin through a square hole without touching, a pin filling a
  round hole, and the commons of a tangent plane and of tangent cylinders).
* The 40 nonempty totals agree with the reference within 8.6e-9 relative in
  volume and 7.5e-9 in area (both `across_hole_common`, a tilted
  cylinder's piece bounded by ellipses and the hole's wall), the centre
  within 6.1e-10 of the case's size; 25 of the 40 within 1e-12. The larger
  differences are BRepGProp's integration of faces bounded by ellipses and
  plane sections of cylinders (the exact reference is checked against
  closed forms within 1e-40 in exact frames); all are inside the 2e-8
  allowance, so no review.
* The degenerate fixtures are valid results natively: the Steinmetz cut as
  two solids touching at the cylinders' two tangent points, the fuses of a
  cylinder and a box tangent along a generatrix and of two tangent
  cylinders as two solids touching along it, and the cut of an internally
  tangent cylinder as one solid touching itself along a line (4 faces, 5
  edges, 2 vertices). The decisions refuse them (`Degenerate`).
* Four solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `steinmetz_common` 5/8/5 to
  4/7/5, `parallel_cylinders_fuse` 9/15/8 to 4/6/4 and
  `parallel_cylinders_common` 5/9/6 to 4/6/4 (the two cylinders' coplanar
  caps merged), and the box of `tangent_plane_fuse` 7/15/10 to 6/12/8 (the
  tangency line imprinted on its face). The kernel's results will be
  compared with the unified counts, as S9a's and S9b's.
* Before the capture, the reference's scan for near coincidences moved two
  fixtures off coincidences of its slicing (two distinct events at one
  slice height, parted by the stored axes' rounding: the parallel
  cylinders' tops both at `y = 3`, and a leaning box's corner edge meeting
  the cylinder at the height where the cylinder's rim met the box's face).
