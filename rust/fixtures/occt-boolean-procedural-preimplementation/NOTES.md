# Booleans of cylinders meeting in quartics (S9c.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_procedural_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9c.2's fixtures). `capture.json` records the Rust revision and
that the kernel's S9c.2 code did not exist (`rust_procedural_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/procedural.rs` absent,
`curved/meet.rs` refusing two cylinders meeting in curves other than lines
and conics, `OutOfDomain(... (S9c.2))`, and the kernel's probe reports every
case `unsupported`). `inputs.txt` holds the 28 cases of
`boolean-procedural-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame
(circles as `gp_Circ` edges about the plane's normal), the frames S9c.1's
(`XY`, `TURN`, `DOWN` against `SIDE` in S9c.2a; `TILT`, `LEAN`, `R125`
against `XY` in S9c.2b).

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 28 (`curved_boolean_reference.py`: the
  regularized Boolean's maximal connected regions), the pipe's two stubs of
  `rings_crossing_cut` and the two pieces of `hole_wall_common` included.
* The totals of the 23 results with a quartic edge (20 declared solid,
  `inside_tangent`'s 3) differ from the reference by more than the 2e-8
  allowance: up to 7.5e-5 relative in volume and 7.8e-5 in area
  (`equal_offset_common`) and 8.7e-6 of the case's size in centre in exact
  frames, up to 1.1e-4 in volume, 3.0e-4 in
  area (`skew_cut`) and 1.4e-4 in centre (`skew_common`) in turned ones;
  `hole_wall_fuse` and `hole_wall_cut` differ in measure only (centres
  within 6.8e-9). The geometry is right: the same native results measured
  by BRepGProp with an adaptive accuracy (`VolumeProperties` and
  `SurfaceProperties` given `Eps = 1e-10`, a diagnostic build of the probe,
  not captured) agree with the reference within 2.7e-9 in volume, 2.0e-9 in
  area and 4.0e-10 in centre (`bite_common` the largest); the edges'
  tolerances are 1.0e-7 to 1.4e-7. The default measure integrates faces
  bounded by the intersection's approximated B-spline curves at a fixed
  order. The 23 are reviewed in `occt-boolean-procedural-divergences.json`,
  each with both measures.
* The reference itself is checked independently
  (`generate_procedural_boolean_fixtures.py --check`): closed forms within
  7.9e-41 in exact frames (one quadrature in `eta`, Legendre's form for
  crossing axes) and 5.0e-16 in turned ones, inclusion and exclusion
  1.2e-41, the area identity 1.3e-40, Monte Carlo 3.2 standard errors.
* The 5 matches have no quartic edge: `parallel_hole_cut` and
  `parallel_hole_fuse` (parallel cylinders: lines and circles) within
  1.8e-16, `outside_tangent_*` within 1.4e-16.
* The degenerate fixtures are valid results natively: the internally
  tangent pipe (`inside_tangent`) as one solid for each operation, its
  meeting's node not marked (fuse 7 faces, 11 edges, 6 vertices; cut 4, 7,
  4; common 3, 3, 2); the externally tangent one's fuse as two solids
  touching at the tangent point, each with a vertex there (3, 3, 3), its cut
  the thick cylinder with that vertex imprinted (3, 3, 3), its common empty.
  The decisions refuse all six (`Degenerate`).
* Six solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `rings_common` 4/6/4 to 3/5/4,
  `bite_common`, `equal_offset_common` and `skew_common` 3/3/2 to 2/2/2, a
  piece of `hole_wall_common` 4/6/4 to 3/5/4, and `parallel_hole_cut`
  9/21/14 to 8/18/12 (faces of one surface merged with their edges).
  The kernel's results will be compared with the unified counts, as S9c.1's;
  its switch and split vertices on the quartics (the decisions' `Curve3::Meet`
  graphs) will differ from OCCT's intersection splines' splits, reviewed
  counts.
* The reference's scan for near coincidences found none; no fixture was
  moved.
