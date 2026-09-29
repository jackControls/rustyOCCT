# Booleans of a sphere against prisms with arcs and of two spheres (S9d.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_spheres_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9d.2's fixtures). `capture.json` records the Rust revision and
that the kernel's S9d.2 code did not exist (`rust_spheres_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/spheres.rs` absent). The
kernel's probe reports every case `unsupported`: a sphere against a prism
with an arc or a circle is refused by S9d.1's arrangement with
`OutOfDomain("a sphere against a cylinder or a sphere (S9d.2)")` (26
cases), two spheres by S9b.2's stored model with `OutOfDomain("a Boolean of
a solid with curved faces or edges in any position (S9c)")` (7 cases: S9c's
arrangement does not take two spheres, and the sphere's face is not a
plane), the messages read from `protocol::run` by a throwaway example. No
protocol or oracle change was needed: either input of a case was already a
`sphere` row, and the oracle builds each side's row on its own.
`inputs.txt` holds the 33 cases of `boolean-spheres-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
sphere or cap is `BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)`, each prism
`BRepPrimAPI_MakePrism` of its profile face (circles, a stadium's lines and
arcs, a square with a circular hole) in its own frame. Frames: `XY`, `SIDE`
(a pipe along `x`), prisms in `TILT`, whole spheres in `TILTX` and `LEAN`.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 33 (`spheres_boolean_reference.py`: the
  regularized Boolean's components followed through its slices): the two
  ends of a pipe and of a rod less the sphere (`pipe_ends_cut`,
  `rod_ends_cut`), the sphere less a box with a coaxial hole whose walls
  cut it, its core and four bulges (`hole_box_cut`, 5), the sphere less a
  stadium through it (`stadium_cut`, 2), two spheres apart or tangent
  fused (2), and the three empty commons (spheres apart, spheres tangent, a
  cylinder tangent outside). The nested cut is one solid with a cavity.
* 18 cases match within the 2e-8 allowance. Where every section is a
  circle of latitude of the sphere's own frame (a pipe along the sphere's
  axis: `pipe_ring`, `pipe_ends`, `dome_pipe`) within 2.3e-14; two whole
  spheres nested, apart or tangent within 2e-16; crossing within 2.9e-9
  (`spheres_turned_common`: the section circle off both spheres'
  parallels); a pipe along `x` against a sphere of axis `z` (`pipe_side`,
  its circles in planes normal to the sphere's equator) within 2.7e-9 in
  centre;
  `hole_box_cut` within 1.2e-9 and `holed_box_cut` within 1.7e-8 (the
  walls' circles off the sphere's parallels).
* 15 cases are reviewed (`occt-boolean-spheres-divergences.json`), all
  `measure` (and `centre`): BRepGProp's default integration on faces
  bounded by the Boolean's approximated intersection curves misses by up to
  2.4e-6 in volume, 1.5e-6 in area and 3.1e-7 of the case's size in centre
  (`bite_common`); every quartic (a rod through the sphere, a bite, a
  cylinder ending inside it, the stadium's arcs, the tangent rod inside)
  and `hole_box_common` (no B-spline edge: circles off the sphere's
  parallels, 3.1e-8). The same native results measured by BRepGProp with
  an adaptive accuracy of 1e-10 (a diagnostic build of the oracle, not
  committed) are within 3.4e-9 in volume, 1.6e-9 in area and 8.8e-10 in
  centre (`bite_common`, `bite_tilt_cut`); their edges' tolerances at most
  7.5e-7 (`bite_tilt_cut`), 2 to 8 of their edges B-splines. The geometry
  agrees with the reference; the default measure does not reach the
  allowance.
* The degenerate fixtures are valid results natively: `spheres_tangent_fuse`
  two solids touching at a point (each sphere keeps the point as a vertex:
  1 face, 3 edges and 3 vertices, the larger 1, 4 and 3 as built), the
  commons of tangent
  spheres and of a cylinder tangent outside empty, `rod_inside_tangent_cut`
  one solid (the rod's hole touching the sphere at a point; reviewed for its
  measure). The decisions refuse them (`Degenerate`).
* Seven solids' counts change under `ShapeUpgrade_UnifySameDomain`
  (faces, edges, vertices as built, then unified): `spheres_cross_common`
  3/3/2 to 2/4/2, `spheres_turned_common` 3/6/4 to 2/3/2, `pipe_side_common`
  4/6/4 to 3/4/2, `hole_box_cut`'s core 3/3/2 to 2/2/1, `bite_common` 3/3/2
  to 2/2/2, `cap_cross_common` 4/6/4 to 3/4/3 and `spheres_tangent_fuse`'s
  larger sphere 1/4/3 to 1/3/2 (faces of one surface split at seams merged).
  The kernel's results will be compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spheres_boolean_fixtures.py --check`): closed forms (two
  spheres' lens as two caps; a sphere against a coaxial cylinder, the ring
  `pi h^3 / 6`; a box with a coaxial hole as S9d.1's box less the core; an
  off-axis cylinder by the lens of two discs along the axis) within 9.2e-41
  in exact frames and 4.2e-17 in turned ones, inclusion and exclusion
  9.2e-41, the area identity 2.3e-40, every face's classes 1.8e-40, a
  second slicing direction (the prism cut obliquely) 1.8e-40, Monte Carlo
  3.1 standard errors; its scan for near coincidences found none and every
  cap circle is at least 0.48 of the radius squared from tangency.
