# Analytic surface intersections (S7a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel surface/surface intersection code existed
(`rust/kernel/src/intersection/analytic.rs` absent). `inputs.txt` holds the
67 cases of `analytic-intersection-cases.txt` with the kernel's stored frame
axes (from `identity_reference.frame_axes`, as `Frame3::new` stores them);
`oracle.cpp` is `occt_analytic_intersection_oracle.cpp` as captured, and
`native.txt` its rows: `IntAna_QuadQuadGeo`'s result type and number of
solutions per case, then each point, line, circle, ellipse, hyperbola or
parabola, with `Precision::Angular()` and `Precision::Confusion()` as the
tolerances IntPatch passes.

## Observations

* Against the independent reference (`analytic_intersection_reference.py`,
  exact rationals and 80-digit mpmath), `compare_analytic_intersections.py`
  finds 63 cases with the same items and parameters within `1e-9` of each
  item's magnitude, after putting OCCT's results in the reference's
  canonical form (a line's point nearest the origin, unit directions with
  their first nonzero coordinate positive, a hyperbola's branches as one
  item).
* The four other cases are near-degenerate configurations IntAna snaps with
  its tolerances (reviewed in `occt-analytic-intersection-divergences.json`):
  a plane `2^-30` rad off a cylinder's axis direction gives two lines
  instead of a 1.07e9-long ellipse, a plane `2^-30` rad off perpendicular a
  circle instead of an ellipse, a plane parallel to a cone's generatrix up
  to rounding a parabola instead of a hyperbola (an exact parabola cannot
  occur with a binary64 half-angle), and a cylinder axis stored as
  `(9.3e-10, 0, 1.0)` becomes exactly parallel under `RefineDir` (two
  lines instead of two ellipses).
* Every exact degeneracy (parallel and coincident planes, tangent spheres,
  planes and cylinders, a plane containing a cone's axis or through a
  rational apex, coaxial pairs) agrees exactly; general cylinder, cone and
  sphere pairs are `NoGeometricSolution` natively and `not_conic` in the
  reference (S7b).

`compare_analytic_intersections.py` requires every later run to reproduce
these rows (types and counts exactly, numbers within `1e-12` relative) and,
once the kernel intersects surfaces, compares its certified results with
both.
