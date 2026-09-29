# Booleans of cylinders in turned frames (S9c.2b.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_turned_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9c.2b.1's fixtures). `capture.json` records the Rust revision
and that the kernel's S9c.2b.1 code did not exist (`rust_turned_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/turned.rs` absent,
`curved/meet.rs` refusing two cylinders not circular in a common measure,
`OutOfDomain("two cylinders meeting in curves other than lines and conics
(S9c.2)")` on every case, and the kernel's probe reports every case
`unsupported`). `inputs.txt` holds the 15 cases of
`boolean-turned-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame
(circles as `gp_Circ` edges about the plane's normal), the frames S9c.1's,
at least one turned per pair (`XY` against `TILTX` and `TILT2`, `TILT`
against `SIDE`, `LEAN` against `TILTX`, `R125` against `SIDE`), the axes
never parallel.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 15 (`curved_boolean_reference.py`: the
  regularized Boolean's maximal connected regions), the pipe's two stubs of
  `rings_pipe_cut` and the two pieces of `hole_tiltx_common` included.
* The totals of the 13 results with a quartic edge (all declared solid)
  differ from the reference by more than the 2e-8 allowance: up to 3.3e-4
  relative in volume (`skew_lean_common`), 2.4e-4 in area (`skew_lean_cut`)
  and 6.1e-5 of the case's size in centre (`skew_lean_common`), the least
  `hole_tiltx_fuse` and `hole_tiltx_cut` (1.7e-7 to 2.4e-7 in measure, 6.0e-8
  in centre). The geometry is right: the same native results measured by
  BRepGProp with an adaptive accuracy (`VolumeProperties` and
  `SurfaceProperties` given `Eps = 1e-10`, a diagnostic build of the probe,
  not captured) agree with the reference within 8.9e-10 in volume, 9.1e-10
  in area (`rings_pipe_cut`) and 2.0e-10 in centre (`hole_tiltx_common`);
  the edges' tolerances are 1.0e-7 to 1.4e-7. The default measure
  integrates faces bounded by the intersection's approximated B-spline
  curves at a fixed order, as in S9c.2's capture. The 13 are reviewed in
  `occt-boolean-turned-divergences.json`, each with both measures.
* The reference itself is checked independently
  (`generate_turned_boolean_fixtures.py --check`): closed forms in the ideal
  frames within 4.2e-16 (the common of crossing cylinders as the
  perpendicular one over `sin phi`, the perpendicular pair with its cap and
  areas, the hole by Cavalieri, the near nodes by Legendre's form at `k =
  1`), inclusion and exclusion 1.2e-41, the area identity 6.0e-41, every
  face's classes 6.9e-41, Monte Carlo 2.4 standard errors; every cap's
  circle clear of the other cylinder (none crosses its model, the least
  clearance 0.45 of the radius squared).
* The 2 matches are the declared near nodes, within 1.5e-10: `node_lean_common`
  (one solid, 6 faces, 12 edges, 8 vertices) and `node_r125_fuse` (8, 16,
  10). OCCT takes the pair as S9c.1's two ellipses within its tolerance
  (the section's edges are `GeomAbs_Ellipse` arcs, 10 and 8, split at the
  nodes and seams; the diagnostic build's edge types, where the 13 others
  have B-spline section edges), so its faces are bounded by conics and the
  default measure is exact enough; the stored models' intervals' ends are
  2.4e-17 and 8.3e-17 apart (rings about 1e-8 apart at the nodes). The
  decisions refuse both (`Degenerate`).
* Four solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `bite_tiltx_common` and
  `skew_lean_common` 3/3/2 to 2/2/2, `skew_lean_cut` 5/6/4 to 4/5/4, and
  `node_lean_common` 6/12/8 to 4/10/8 (faces of one surface merged with
  their edges). The kernel's results will be compared with the unified
  counts, as S9c.2a's; its switch and split vertices on the quartics will
  differ from OCCT's intersection splines' splits, reviewed counts.
* The reference's scan for near coincidences found only the declared
  nodes' breakpoints (two pairs about 1e-17 apart each); no fixture was
  moved.
