# Booleans of a whole torus against prisms with arcs, spheres, cones and tori (S9d.4b.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_torus_curved_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9d.4b.2's fixtures). `capture.json` records the Rust
revision (the fixtures and the reference committed, nothing else changed)
and that the kernel's S9d.4b.2 code did not exist
(`rust_torus_curved_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/torus_curved.rs` absent). The
kernel's probe reports every case `unsupported`: the curved Booleans refuse
a torus against anything but a polyhedral prism with `OutOfDomain("a torus
against a curved face (S9d.4b)")` (all 44 cases; every input constructs, so
the kernel's test support reads the rows). `inputs.txt` holds the 44 cases
of `boolean-torus-curved-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`: a torus's block one `torus` row,
a sphere's one `sphere` row, a cone's one `cone` row, a prism's its `plane`,
`wire` and `prism` rows; no row is new), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured (unchanged: every older capture
reproduces) and `native.txt` its results; each torus is
`BRepPrimAPI_MakeTorus(gp_Ax2, R, r)`, each sphere `BRepPrimAPI_MakeSphere`,
each cone `BRepPrimAPI_MakeCone`, each prism `BRepPrimAPI_MakePrism` of its
profile face in its own frame. Frames: `XY` but for the rod and the linked
and crossing tori (`SIDE`), the cone in `LEAN`, and the torus in `TILT` with
a coaxial pipe and with a sphere on its tube.

## Observations

* Every case is done without warnings, every result and solid valid. The
  solid counts are the reference's in all 44
  (`torus_curved_boolean_reference.py`: the torus's normal slices swept,
  each slice's rays from the axis joined where they overlap): the rod's
  common with the torus two solids (its two crossings of the tube) and the
  rod less the torus three (its ends and its middle in the hole), the parallel
  torus's common two, the tori linked apart, the pin and the ball in the
  hole fused two solids each (their commons empty).
* 24 cases match within the 2e-8 allowance: where every section is a
  circle (the coaxial pairs: cylinders, spheres, a cone and a torus about
  the torus's axis, in `XY` and in `TILT`) within 1.6e-13 (BRepGProp's
  integration on the toroidal faces), the pairs apart and the degenerate
  ones within 2.1e-16; three whose results hold B-spline sections within
  1.9e-8 (`bore_cut`, `tori_ring_fuse`, `tori_ring_cut`).
* 20 cases are reviewed (`occt-boolean-torus-curved-divergences.json`),
  all `measure` (and eleven `centre`): BRepGProp's default integration on
  faces bounded by the Boolean's approximated sections misses by up to
  9.8e-6 in volume (`core_ball_cut`), 3.9e-6 in area and 4.5e-6 of the
  case's size in centre (`tori_side_common`). The same native results
  measured by BRepGProp with an adaptive accuracy of 1e-10 (a diagnostic
  build of the oracle, not committed) are within 4.8e-9 in volume, 2.7e-9 in
  area and 3.1e-9 in centre, unchanged at an accuracy of 1e-12: what remains
  is the approximated sections themselves (B-splines within their edges'
  tolerance), well inside the allowance.
* Curves as built: circles where a section is one (a coaxial pair's walls
  and a cap normal to the axis); every other section of the torus and a
  cylinder, sphere, cone, plane or torus is a B-spline approximation, 121 of
  the results' 329 edges. The results' edges' tolerances are 1e-7, 1e-6 on
  the stadium's (its flat walls meeting its arcs tangentially).
* The degenerate fixtures are valid results natively: `pipe_equator_fuse`
  the cylinder (the torus inside it, tangent along the outer equator),
  `ball_touch_fuse` two solids touching at a point, `tori_kiss_fuse` two
  tori touching along a circle. The decisions refuse them (`Degenerate`).
* Thirteen solids' counts (ten cases) change under
  `ShapeUpgrade_UnifySameDomain` (faces, edges, vertices as built, then
  unified): the rod's commons' and ends' pieces (`rod_common` 7/13/8 to
  6/12/8 and 4/6/4 to 3/5/4, `rod_ends_cut` 4/5/3 to 3/4/3 and 4/7/5 to
  3/6/5), `bore_common` and `spike_common` 5/9/6 to 3/7/6, `ball_top_common`
  3/5/3 to 2/4/3, `ball_core_cut` 5/8/5 to 4/7/5, `cone_lean_common` 3/6/4 to
  2/5/4, both of `tori_side_common`'s 4/6/4 to 3/5/4, `pipe_equator_fuse`
  4/5/3 to 3/3/2 and `tori_kiss_fuse`'s outer torus 2/4/2 to 1/0/0 (faces
  split at OCCT's seams merged). The kernel's results will be compared with
  the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_torus_curved_boolean_fixtures.py --check`): closed forms of 23
  of the 27 pairs (one quadrature along the torus's axis of its annulus
  against the other input's sections: lenses of signed discs about a
  parallel axis, a rod's strip by rectangle-in-disc antiderivatives, inputs
  apart by their sums) within 7.5e-39 in exact frames and 4.6e-17 in turned
  ones (the declared degenerate pairs 3.1e-24), every operation two ways
  (each face swept by two families of curves) within 6.6e-35, both inputs'
  measures from their faces 1.4e-39 from their closed forms, inclusion and
  exclusion 1.1e-40, the area identity 1.4e-39, Monte Carlo 2.8 standard
  errors. Its scan for near coincidences found none outside the declared
  degenerate pairs, and every other pair's surfaces meet at a sine of at
  least 0.42 (`tori_side`), keep 0.066 of the case's size apart where they
  do not meet (`pin_hole`), its edges cross the torus at a sine of at least
  0.62 and stay 0.031 from tangency (`sleeve`'s lower rim, coaxial inside the
  tube), and the cone's apex lies 0.13 inside the tube (`cone_lean`).
