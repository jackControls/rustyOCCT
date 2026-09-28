# Booleans of polyhedral prisms in any position (S9b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_polyhedral.py --capture` (`compare_boolean.py`'s probe and
protocol on S9b's fixtures). `capture.json` records the Rust revision and
that the kernel's S9b module did not exist (`rust_boolean_exists` false:
`rust/kernel/src/solid/boolean/polyhedra.rs` absent; the kernel refuses
frames with different axes). `inputs.txt` holds the 45 cases of
`boolean-polyhedra-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame, the
two frames' axes different (`XY`, `TILT`, `ROT`, `SIDE`, `LEAN`).

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 45 (`polyhedral_reference.py`: the
  regularized Boolean's maximal connected regions, exact convex cells), the
  6 empty results included, and the totals agree with it within 6.1e-16
  relative in volume and 7.9e-16 in area (`across_hole_common`), the centre
  within 3.0e-16 of the cases' size (`caps_same_fuse`). The reference is the
  exact model on the stored axes (a prism's volume its profile's area times
  its height times `det(x, y, n)`); OCCT builds on the same stored values.
* Fuses of prisms touching along an edge (`edge_on_face_fuse`) or at a
  vertex (`vertex_on_face_fuse`) are two valid solids natively; the fixtures
  declare them `degenerate` (the kernel's `Degenerate` until it holds
  non-manifold bodies), as S9a's touching cases.
* Four solids' counts change under `ShapeUpgrade_UnifySameDomain` (coplanar
  faces of both inputs, and collinear edges, merged): the kernel's results
  are compared with the unified counts, as S9a's.
* Before the capture, two fixtures were moved off exact coincidences the
  reference resolves and OCCT's tolerances do not (a bar's side through a
  hole's corner, a face through the other's edge): each gave one solid in
  the exact model and two natively.
* After the capture the kernel met a third near-coincidence the fixtures
  had not declared: in `ell_tilted`'s fuse and cut the tool's top edge
  passes within rounding of the L's corner edge, leaving a wall of two
  triangles joined by a neck thinner than the resolution in the exact
  model (one face) and two faces touching at a point natively. Both are
  now declared `degenerate` (the native inputs, and so this capture, are
  unchanged); the common is unaffected.
