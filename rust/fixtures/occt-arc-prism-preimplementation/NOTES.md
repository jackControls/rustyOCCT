# Arc prisms (S5): native MakePrism observations before any kernel arc code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code built prisms of profiles with arcs: the
kernel's `Boundary` had polygons and circles only. The inputs are the
thirteen cases of `generate_identity_fixtures.arc_cases()` as explicit OCCT
constructions (`identity_reference.native_case`: an `S` wire gives each
stored point and its segment to the next, a line or an arc about the plane's
normal or its reverse); `oracle.cpp` is `occt_history_oracle.cpp` as
captured, and `native.json` holds its complete output per case: the
profile face's validity, every `BRepTools_History`/`MakePrism` query of
every profile vertex, edge and face, every output subshape's coverage, the
distinct subshape counts, the structure-only subshapes and each transform
step.

## Observations

* Every construction is BRepCheck-valid and its history complete: every
  output vertex, edge and face is reached directly by a query.
* An arc segment gives one cylindrical wall, whatever its sweep (over half
  a turn in the pac-man cases, two arcs of one circle in
  `arc_holes_mixed`), and the counts are a polygon's with the same number of
  segments: `2n` vertices, `3n` edges, `n + 2` faces per outer boundary.
* The only structure-only subshapes are the full circle hole's seam edge and
  its two vertices in `arc_holes_mixed`; no boundary with arcs has a seam.

`compare_history.py --family arc` requires every later run to reproduce
these rows and compares the kernel's histories with them.
