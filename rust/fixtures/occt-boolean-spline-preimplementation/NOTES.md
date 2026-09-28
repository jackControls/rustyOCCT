# Booleans of spline profiles (S9a.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28. `capture.json`
records the Rust revision and that the kernel's Booleans of spline profiles
did not exist (`rust_spline_boolean_exists: false`: `rust/kernel/src/profile/
boolean.rs` still refuses every profile holding a spline segment with
`OutOfDomain("a Boolean of spline profiles (S9a.2)")`). `inputs.txt` holds
the 46 cases of `boolean-spline-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`); `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured (S9a's oracle with S8b's `B` path
segments built as `Geom_BSplineCurve` edges between the path's vertices, as
`occt_split_oracle.cpp` builds them), and `native.txt` its results, in the
shape of S9a's capture.

## Observations

* Every case is done, every result and solid valid, no warnings, and every
  solid count is the reference's (`boolean_reference.SplinePair`; 4 results
  empty). 37 cases agree with the reference within the 2e-8 allowance
  (1.8e-8 in volume at worst, `blobs_common`, 2.0e-8 in area, 2.3e-9 of the
  case's size in the centre); most within 1e-14.
* 9 cases differ beyond it, each reviewed in
  `occt-boolean-spline-divergences.json`: `capsule_parabola_common` (3.8e-8
  in volume), `capsule_parabola_cut_tilted` (1.4e-7), `blob_inside_cut`
  (1.2e-7), `blob_inside_common_tilted` (1.7e-7), `blob_pocket_cut` (2.6e-8),
  `blobs_cut_tilted` (1.8e-7), `waves_fuse_tilted` (5.0e-6), `waves_common`
  (1.3e-4) and `waves_step_fuse_tilted` (8.0e-4). All are BRepGProp's
  integration of faces bounded by B-spline edges, not OCCT's geometry: a
  scratch variant of the oracle (not the captured probe) integrating
  `u dv - v du` over OCCT's own edges of every planar cap face finds the
  caps' areas the reference's slab areas within 4.5e-8 in every case
  (within 2.6e-9 in the eight reviewed cases other than `waves_common`, whose
  cap edges end at wave/wave sections of tolerance 2.6e-5), where
  BRepGProp's own face areas err by up to 1.3e-3; BRepGProp with an adaptive
  accuracy of 1e-13 brings the blob and capsule cases and `waves_fuse_tilted`
  within 8.2e-10, but not `waves_common` and `waves_step_fuse_tilted` (1.9e-4
  and 1.8e-3): their faces are bounded by the wave, the quadratic B-spline of
  three spans whose faces S8b's split capture already reviewed
  (`occt-split-divergences.json`:
  `wave_knot`, `wave_parallel_four`, `wave_tilted`). OCCT is also
  inconsistent with itself there: `waves_common`'s two lenses, mirror images
  about `u = 5`, have volumes 14.534320 and 14.536522; `blob_inside_cut`
  should be 500 - 261.25 = 238.75 exactly and is 238.7500297; OCCT's own cut
  and common of the two blobs sum to 261.2499799, not the blob prism's
  261.25.
* Tangencies (a spline touching a line or a circle): OCCT's results have the
  reference's volume and solid count, and `dome_tangent_disc_cut` (a hole
  tangent to the outer boundary at the dome's apex) is one valid solid, as
  S9a's `tangent_in_cut`; the decisions make it `Degenerate` (the fixture's
  `expect degenerate`). OCCT keeps the tangency point as a vertex and edges
  where the decisions say a tangency cuts nothing: the dome with a line
  touching its apex from outside (`dome_touch_line_cut_tilted`, the dome
  itself) has 4 faces, 10 edges and 9 vertices after unifying, where the
  dome prism has 4, 6 and 4; `dome_tangent_line_fuse_tilted` 12, 37, 32 and
  `dome_tangent_disc_fuse_tilted` (the dome again) 4, 16, 22. The kernel's
  counts will differ there and be reviewed.
* Shared splines: identical profiles and a spline shared in the same or the
  opposite direction give the reference's results within 2.3e-16 (the
  opposite-direction fuse is the rectangle [0, 14] x [0, 6], 6 faces after
  unifying; the lens hole filled is the square, 6 faces).
* Crossings at interior knots of both splines (`waves_at_knots_cut`) and a
  vertex exactly on a spline (`vertex_on_spline_cut`) agree within 4.2e-9
  and 1.6e-14.

`compare_boolean.py --splines` requires every later run to reproduce these
rows (on another platform, its reviewed record; the Linux record is pending
CI). While the kernel refuses spline profiles, the probe
(`rust/kernel/examples/boolean_probe.rs`) must report `unsupported` on all 46
cases, which the comparison lists under `rust_unsupported`; any other row, or
a probe failure, is a failure.
