# Booleans of a cone against prisms with arcs, spheres and cones (S9d.3b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_cones_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9d.3b's fixtures). `capture.json` records the Rust revision and
that the kernel's S9d.3b code did not exist (`rust_cones_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/cones.rs` absent). The
kernel's probe reports every case `unsupported`: S9d.3a's arrangement
refuses a cone against a prism with an arc or a circle, a sphere or another
cone with `OutOfDomain("a cone against a prism with arcs, a sphere or a cone
(S9d.3b)")` (all 40 cases, the message read from `protocol::run` by a
throwaway example). No protocol or oracle change was needed: either input
of a case was already a `cone` or `sphere` row, and the oracle builds each
side's row on its own. `inputs.txt` holds the 40 cases of
`boolean-cones-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
cone or frustum is `BRepPrimAPI_MakeCone(gp_Ax2, R1, R2, H)`, each sphere or
cap `BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)`, each prism
`BRepPrimAPI_MakePrism` of its profile face (circles, a stadium's lines and
arcs, a square with a circular hole) in its own frame. Frames: `XY`, `SIDE`
(rods and a cone along `x`), `DOWN` (cones pointing down), cones and a
pipe in `TILT`, a cone and a cylinder in `LEAN`.

## Observations

* Every case is done, every result and solid valid; `rod_ruling_common` (a
  cylinder tangent to a cone along a ruling) and `apex_to_apex_fuse` (two
  cones touching at their apexes) report warnings. The solid counts are the
  reference's in 39 of the 40 (`cones_boolean_reference.py`: the
  regularized Boolean's components followed through its slices): the two
  ends of a pipe, a rod and a crossing cone less the frustum
  (`pipe_ends_cut`, `rod_ends_cut`, `cones_ends_cut`), a cone less a rod
  across its tip, its tip and its base (`tip_rod_cut`), two cones apex to
  apex fused (2), and the 2 empty commons (a frustum on another's top disc,
  the cylinder tangent along a ruling). The other is
  `ball_inscribed_cut`, declared degenerate: OCCT's one solid (the cone with
  the ball's cavity touching its wall all round a circle) against the
  reference's two regions meeting along that circle.
* 23 cases match within the 2e-8 allowance. Where every section is a
  circle about the cone's axis (a coaxial pipe, a box with a coaxial hole,
  coaxial spheres and a dome, coaxial cones apex to apex, frusta crossing,
  stacked, equal) within 4.1e-14 (`ball_coax_cut`: the spheres' faces;
  without a sphere within 1.1e-15), the pipe in `TILT` within 1.3e-16;
  `stadium_common` within 2.0e-8 (its flat walls' hyperbolas, its arcs'
  B-splines); `cones_parallel_common` within 6.1e-10: two cones of
  equal half-angle with parallel axes meet in a plane (the difference of
  their equations is linear), and OCCT's section is exact, four
  hyperbolic arcs.
* 17 cases are reviewed (`occt-boolean-cones-divergences.json`). 16 are
  `measure` (and `centre`): BRepGProp's default integration on faces
  bounded by the Boolean's approximated intersection curves misses by up to
  4.1e-6 in volume (`cones_ends_cut`), 2.6e-6 in area and 7.4e-7 of the
  case's size in centre (`bite_common`); every quartic (a cylinder through
  the wall off the axis, rods across a frustum and a cone's tip or apex, a
  sphere off the axis, crossing cones, the stadium's arcs). The same native
  results measured by BRepGProp with an adaptive accuracy of 1e-10 (a
  diagnostic build of the oracle, not committed) are within 6.0e-9 in
  volume (`stadium_cut`), 2.3e-9 in area (`apex_rod_common`) and 1.1e-9 in
  centre; their edges' tolerances at most 1.9e-7 (`apex_rod_common`), 2 to
  8 of their edges B-splines. The geometry agrees with the reference; the
  default measure does not reach the allowance. The seventeenth is
  `ball_inscribed_cut`'s solid count (above), its measures within 1.7e-16.
* Curves as built: lines and circles where the sections are circles; the
  stadium's flat walls meet the frustum in hyperbolas, `cones_lean`'s
  sections include four ellipses and `cones_parallel`'s are hyperbolas;
  every other intersection of two quadrics is a B-spline approximation.
  The cones' apexes are degenerated edges where a result keeps them.
* The degenerate fixtures are valid results natively:
  `rod_ruling_common` empty (the reference empty: the stored frame of
  `LEAN` leaves the cylinder 1e-16 off the ruling), `ball_inscribed_cut`
  above, `apex_to_apex_fuse` two solids touching at a point (2 faces, 3
  edges and 2 vertices each, the apex a degenerated edge), `cones_equal_fuse`
  the cone. The decisions refuse them (`Degenerate`).
* Coplanar discs: `frusta_stack_fuse` (a frustum standing on the frustum's
  top disc, opposite orientations) one solid of 5 faces, the smaller disc
  gone and the larger holed at its rim; the common empty.
  `frusta_cross_common` shares both end planes with the same orientation.
* Eleven solids' counts change under `ShapeUpgrade_UnifySameDomain`
  (faces, edges, vertices as built, then unified): `bite_common` 4/6/4 to
  3/4/3, `rod_common` and `rod_ends_cut`'s end, `cones_cross_common` and
  `cones_ends_cut`'s small end 4/6/4 to 3/5/4, `stadium_cut` and
  `stadium_common` 9/21/14 to 8/19/13, `ball_side_common` 3/6/4 to 2/5/4,
  `ball_tilt_common` 4/9/6 to 3/8/6, `cones_lean_common` 5/9/6 to 4/8/6,
  `cones_parallel_common` 4/6/4 to 3/4/3 (the walls' faces split at OCCT's
  seams merged). The kernel's results will be compared with the unified
  counts, as S9c's.
* The reference itself is checked independently
  (`generate_cones_boolean_fixtures.py --check`): closed forms of 22 of the
  27 pairs (the lens of two discs along the axis for cylinders, spheres
  and cones of parallel axes, coaxial or not; the strip of a rod across the
  axis; S9d.3a's box less the coaxial core) within 1.7e-40 in exact frames
  and 4.4e-17 in turned ones, inclusion and exclusion 1.2e-41, the area
  identity 1.1e-40, every face's classes 5.0e-41, the cone's wall two ways
  3.6e-41, a second slicing direction 4.6e-41, Monte Carlo 2.9 standard
  errors. Its scan for near coincidences found none and every edge and
  vertex is at least 0.0126 (relative) from tangency with or incidence on a
  surface of the other input.
