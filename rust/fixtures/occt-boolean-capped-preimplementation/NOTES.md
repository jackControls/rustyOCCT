# Booleans of cylinders whose sections cross caps' circles (S9c.2b.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_capped_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9c.2b.2's fixtures). `capture.json` records the Rust revision
and that the kernel's S9c.2b.2 code did not exist (`rust_capped_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/algebraic.rs` absent; S9c.2a's
`procedural.rs` refusing the 9 exact-frame cases with
`OutOfDomain("two cylinders' section crossing a cap's circle (S9c.2b)")`
and S9c.2b.1's `turned.rs` the 9 turned-frame ones with
`OutOfDomain("two cylinders' section crossing a cap's circle in turned
frames (S9c.2b.2)")`, so the kernel's probe reports every case
`unsupported`). `inputs.txt` holds the 18 cases of
`boolean-capped-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame
(circles as `gp_Circ` edges about the plane's normal), the frames S9c.1's:
`XY` against `SIDE` (perpendicular, exact) and `XY` against `TILTX`,
`LEAN` and `TILT` pipes (turned), the axes never parallel.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 18 (`curved_boolean_reference.py`: the
  regularized Boolean's maximal connected regions), the two pieces of
  `hole_rim_common` and `hole_rim_tiltx_common` included.
* Every result has a quartic edge through a vertex on a cap's circle. 16
  differ from the reference by more than the 2e-8 allowance: up to 1.4e-5
  relative in exact frames (`rim_bite_common`'s area) and 1.9e-5 in turned
  ones (`hole_rim_tiltx_common`'s volume), 3.1e-6 of the case's size in
  centre (`blind_tilt_common`). The geometry is right: the same native
  results measured by BRepGProp with an adaptive accuracy
  (`VolumeProperties` and `SurfaceProperties` given `Eps = 1e-10`, a
  diagnostic build of the probe, not captured) agree with the reference
  within 1.0e-9 in exact frames and 1.9e-8 in turned ones
  (`blind_tilt_common`: a result of volume 0.19 whose one section edge has
  tolerance 1.1e-7; the same at `Eps = 1e-12`, so the rest is OCCT's
  approximated section within its tolerance, not the integration); the
  edges' tolerances are 1.0e-7 to 1.4e-7. The default measure integrates
  faces bounded by the intersection's approximated B-spline curves at a
  fixed order, as in S9c.2's and S9c.2b.1's captures. The 16 are reviewed
  in `occt-boolean-capped-divergences.json`, each with both measures.
* The 2 matches are `rim_pipe_cut` (9.0e-9 in volume, 1.4e-8 in area) and
  `enter_lean_common` (4.9e-9, 1.3e-8), where the default measure happens
  to fall within the allowance.
* The reference itself is checked independently
  (`generate_capped_boolean_fixtures.py --check`): the exact-frame pairs'
  closed forms (the common of perpendicular cylinders clipped by both
  inputs' caps by one quadrature in `eta`, areas included; the hole as the
  pipe's segment below the box's top face less its common with the hole's
  cylinder) within 6.4e-41, inclusion and exclusion 1.5e-41, the area
  identity 5.9e-41, every face's classes 3.3e-41, Monte Carlo 2.9 standard
  errors; the classes from the intervals as declared (ends at least 0.6
  apart); each pair's caps' circles crossing the other cylinder on its face
  (20 crossings, each at least 2.0 inside the other face along its axis,
  at slopes of at least 0.25; 2 more of `XY`'s bottom circle against
  `enter_lean`'s pipe, 4.6 outside its face).
* Eight solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `rim_pipe_cut` 6/12/8 to 5/11/8,
  `rim_pipe_common` 6/12/8 to 4/9/7, `rim_bite_common` 4/6/4 to 3/4/3,
  `hole_rim_cut` 12/33/22 to 10/29/20, `hole_rim_common` 5/9/6 to 4/7/5 and
  6/12/8 to 4/8/6, `enter_lean_common` 5/9/6 to 4/7/5 and
  `hole_rim_tiltx_common`'s smaller piece 4/6/4 to 3/5/4 (faces of one
  surface merged with their edges). The kernel's results will be compared
  with the unified counts, as S9c.2's; its switch and split vertices on the
  quartics will differ from OCCT's intersection splines' splits, reviewed
  counts.
* The reference's scan for near coincidences found none. Perpendicular
  pairs in stored turned frames (`TILT` against `TILTX`) were tried first
  and dropped: one input's cap plane is then parallel to the other's
  generatrices only within the stored axes' rounding, and the scan found
  breakpoints about 1e-17 apart where the pipe's wall meets that plane in
  two nearly parallel lines. The turned pairs are oblique instead.
