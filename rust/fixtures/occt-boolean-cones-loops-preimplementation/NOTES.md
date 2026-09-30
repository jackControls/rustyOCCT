# Booleans of cones meeting cylinders, cones and spheres in loops, and turned caps against cones (S9d.3c): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_cones_loops_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9d.3c's fixtures). `capture.json` records the Rust
revision and that the kernel's S9d.3c code did not exist
(`rust_cones_loops_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/cones_loops.rs` absent). The kernel's
probe reports every case `unsupported`: a cone meeting a cylinder or a cone
in loops, the cylinder exactly along a cone's ruling among them, is refused
with `OutOfDomain("a cone meeting a curved face in a loop (S9d.3b.2)")` (17
cases: `rod_graze`, `rod_bitten`, `tip_graze`, `rod_lean`, `tilt_rod`,
`cones_graze`, `cones_turned`, `cones_asymptotic`, `rod_ruling_exact`), a
turned cone meeting a sphere in a loop with `OutOfDomain("a sphere meeting
a turned cone in a loop (S9d.3c)")` (8: `ball_tilt`, `ball_lean`,
`ball_r125`, and `dome_cone_turned`, whose loop is refused before its
cap's circles), a turned cap's circle against a cone with `OutOfDomain("a
turned cap's circle against a cone (S9d.3c)")` (6: `dome_tilt_cone`,
`dome_lean_frustum`, `rim_tangent_cone`), the messages read from
`protocol::run` by a throwaway example. No protocol or oracle change was
needed. This is the second capture of the set, before the kernel code as
the first: the first had two pairs with a plane tangent to a sphere
(`ball_r125`, `dome_lean_frustum`), which S9d.1's rule makes `Degenerate`
wherever it touches; the generator now checks for them and the pairs moved
clear (the evidence's correction, its own commit). `inputs.txt` holds the 31 cases of `boolean-cones-loops-cases.txt`
as explicit constructions (`identity_reference.native_boolean_case`),
`oracle.cpp` is `occt_boolean_oracle.cpp` as captured and `native.txt` its
results; each cone or frustum is `BRepPrimAPI_MakeCone` in its frame, each
sphere or hemisphere `BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)`, each
prism `BRepPrimAPI_MakePrism` of its circle in its own frame. Frames:
cones in `XY`, `TILT`, `LEAN`, `SIDE` and `R125`; rods in `SIDE`, `LEAN`
and `XY`; hemispheres in `TILT` and `LEAN`.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 31 (`cones_boolean_reference.py`, each
  result one solid).
* 7 cases match within the 2e-8 allowance: the turned hemispheres against
  cones (`dome_tilt_cone`, `dome_lean_frustum`, `dome_cone_turned`, all
  within 1.8e-8) and the declared degenerate `rim_tangent_cone_common`
  (3.2e-9).
* 24 cases are reviewed (`occt-boolean-cones-loops-divergences.json`), all
  `measure` (and `centre` but for `tip_graze_cut`, `tilt_rod_common` and
  `cones_turned_cut`): BRepGProp's default integration on faces bounded by
  the Boolean's approximated loops misses by up to 5.8e-6 in volume
  and 4.0e-6 in area (`rod_graze_common`) and 2.2e-7 of the case's size
  in centre (`rod_lean_common`), and 4.5e-4 in volume on the cylinder along a
  cone's ruling (`rod_ruling_exact_common`, declared degenerate). The same
  native results measured by BRepGProp with an adaptive accuracy of 1e-10
  (a diagnostic build of the oracle, not committed; at 1e-12 the same to
  two digits but one centre, 7.4e-12 for 7.6e-12, so converged) are within 1.7e-8 in volume (`rod_lean_common`), 1.3e-8 in
  area (`tilt_rod_common`) and 4.2e-9 in centre, their edges' tolerances
  at most 2.7e-7 and up to 8 of their edges B-splines. The geometry agrees
  with the reference within the native tolerances; the default measure
  does not reach the allowance.
* The degenerate fixtures are valid results natively:
  `rod_ruling_exact_common` (the cylinder exactly along the cone's ruling,
  its curve one branch running off to infinity) and `rim_tangent_cone_common`
  (the hemisphere's rim touching the cone's wall, the reference's clearance
  9.6e-17 of the radius squared). The decisions refuse both (`Degenerate`).
* Five solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `tilt_rod_common` 6/12/8 to
  5/10/7, `ball_tilt_common` 4/9/6 to 3/8/6, `ball_lean_common` 3/6/4 to
  2/5/4, `dome_cone_turned_cut` 4/9/6 to 3/8/6 and `dome_cone_turned_common`
  5/9/6 to 3/6/5 (faces of one surface split at seams merged). The kernel's
  results will be compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_cones_loops_boolean_fixtures.py --check`): closed forms (a rod
  across a cone's axis by S9d.3b's strip, a whole sphere against a turned
  cone by the lens of two discs along its ideal axis, a hemisphere against
  a coaxial cone by circular segments along the axis, its faces by their
  arcs and the disc by its chords) within 2.2e-40 in exact frames and
  1.3e-16 in turned ones, an input cut in two along its axis against the
  whole within 5.3e-42, a hemisphere and its complement against the whole
  sphere within 1.4e-41, inclusion and exclusion 1.2e-41, the area
  identity and every face's classes 6.9e-41, the cone's wall two ways
  2.7e-41, a second slicing direction 1.7e-41, Monte Carlo 2.5 standard
  errors; its scan for near coincidences found none but in the declared
  degenerate pairs, every other edge and vertex at least 0.0156 (relative)
  from tangency with or incidence on the other's surfaces, every plane at
  least 0.05 of the radius from tangency to the other's spheres and
  parallel cylinders. Python 3.9 and
  3.12 write the same files.
