# Booleans of spheres' caps and zones against tori, and of torus v-segments and wedges against curved solids (S9d.4c): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_torus_parts_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9d.4c's fixtures). `capture.json` records the Rust
revision and that the kernel's S9d.4c code did not exist
(`rust_torus_parts_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/torus_parts.rs` absent). The kernel's
probe reports every case `unsupported`: a sphere's cap or zone against a
torus is refused with `OutOfDomain("a sphere's circle of a surd radius
against a torus (S9d.4b)")` (16 cases: `zone_coax`, `cap_coax`,
`coax_cap`, `cap_top`, `zone_side`, `dome_lean`, whose hemisphere's rim is
of a rational radius but on a basis of unequal lengths in `LEAN`,
`zone_tilt` and `rim_touch`), a torus v-segment or wedge against a curved
face with `OutOfDomain("a torus segment or wedge against a curved face
(S9d.4b)")` (37: every `oh_`, `ih_`, `band_`, `qw_`, `hw_` and `tw_` case),
the messages read from `protocol::run` by a throwaway example. No protocol
or oracle change was needed: the native oracle's `torus ... R r LOW HIGH
ANGLE` row (S9d.4b.1) builds each part by `BRepPrimAPI_MakeTorus(gp_Ax2, R,
r, low, high, angle)` and its `sphere R LOW HIGH` row each cap or zone by
`BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)`. `inputs.txt` holds the 53
cases of `boolean-torus-parts-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
cone or frustum is `BRepPrimAPI_MakeCone` in its frame, each prism
`BRepPrimAPI_MakePrism` of its circle in its own frame. Frames: tori and
parts in `XY`, `TILT` and `LEAN`; spheres in `XY`, `SIDE` and `LEAN`;
prisms and cones in `XY`.

## Observations

* Every case is done, 52 of the 53 results valid with no warnings, and the
  solid counts are the reference's in all 53
  (`torus_parts_boolean_reference.py`; `ih_pipe_cut`, `ih_ball_cut` and
  `oh_dome_cut` two solids, the others one).
* 26 cases match within the 2e-8 allowance: every coaxial pair (circles
  within 2.4e-14: `zone_coax`, `cap_coax`, `coax_cap`, `oh_pipe`,
  `ih_pipe`, `ih_cone`, `ih_ball`, `oh_dome`, `band_tilt_pipe` in `TILT`)
  and `oh_ball_cut`, `band_cone_cut` and `band_cone_common`, `qw_pipe` (all
  three) and `hw_ball_cut` (within 1.9e-8).
* 27 cases are reviewed (`occt-boolean-torus-parts-divergences.json`), 25
  of them `measure` (and `centre` for 14): BRepGProp's default
  integration on faces bounded by the Boolean's B-spline sections misses
  by up to 7.4e-6 in volume (`zone_tilt_common`), 4.3e-6 in area
  (`band_torus_common`) and 9.0e-7 of the case's size in centre. The same
  native results measured by BRepGProp with an adaptive accuracy of 1e-10
  (a diagnostic build of the oracle, not committed; unchanged at 1e-12, so
  converged) are within 2.3e-8 in volume (`tw_cone_common`), 2.1e-8 in area
  (`dome_lean_common`) and 4.2e-9 in centre. The geometry agrees with the
  reference within the native tolerances; the default measure does not
  reach the allowance.
* The quarter wedge in `LEAN` against a sphere across its start disc's rim
  (`qw_lean_ball`) is wrong natively: its cut is valid but 3.0e-6 off in
  volume and 4.7e-5 in area, adaptively measured, and its common is
  invalid under `BRepCheck_Analyzer` (the result and its solid) and 1.6e-4
  off in volume; the reference's two families of every face agree within
  1.5e-31 there. Both are reviewed on the reference's evidence.
* The degenerate fixtures are valid results natively: `rim_touch_fuse` (a
  cap's rim within the resolution of tangency to the torus's outer equator)
  and `qw_rim_touch_common` (a quarter wedge's start rim touching a
  sphere). The decisions refuse both (`Degenerate`).
* Fourteen solids' counts change under `ShapeUpgrade_UnifySameDomain`
  (faces, edges, vertices as built, then unified): `cap_top_common` 5/9/6 to
  4/6/4, `zone_side_cut` 5/14/9 to 4/12/8, `zone_side_common` 8/16/10 to
  7/14/9, `dome_lean_cut` 4/11/7 to 3/9/6, `dome_lean_common` 5/9/6 to
  3/6/5, `zone_tilt_common` 5/9/6 to 4/8/6, `oh_ball_common` 4/8/5 to 3/6/4,
  `band_ball_common` 4/6/4 to 3/4/3, `band_cone_common` 5/8/5 to 4/6/4,
  `band_torus_cut` 5/12/8 to 4/9/6, `qw_pipe_cut` 6/12/8 to 5/11/8,
  `qw_pipe_common` 5/9/6 to 4/8/6, `qw_rim_touch_common` 3/3/2 to 2/2/2 and
  `qw_lean_ball_common` 5/13/9 to 4/11/8 (faces of one surface split at
  seams merged). The kernel's results will be compared with the unified
  counts, as S9c's.
* The reference itself is checked independently
  (`generate_torus_parts_boolean_fixtures.py --check`): closed forms of the
  9 coaxial pairs (one quadrature along the torus's axis of both inputs'
  radial intervals, the walls by their own elements, a wedge by its turn and
  its meridian discs) within 6.9e-40 in exact frames and 3.8e-17 in turned
  ones, a half and the other half against the tool as the whole torus
  (S9d.4b.2's reference) and a cap or zone and the rest of its sphere
  against the torus as the whole sphere within 4.3e-42, every operation two
  ways (each face's two families of curves) 1.5e-31, both inputs from their
  faces 8.3e-35, inclusion and exclusion 1.8e-35, the area identity
  8.3e-35, Monte Carlo 2.5 standard errors; its scan for near coincidences
  found none but in the declared degenerate pairs, every other pair's
  surfaces meeting at a sine of at least 0.38, its edges (the rims too)
  crossing the other's surfaces at a sine of at least 0.29, and every plane
  at least 0.0059 of the case's size from tangency to the other's spheres
  and tori. Python 3.9 and 3.12 write the same files.
