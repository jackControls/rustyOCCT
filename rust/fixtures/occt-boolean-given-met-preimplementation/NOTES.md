# Given results whose meetings of two curved faces, or cones' and tori's general sections, the partner meets (S9e.3b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-10-03 by
`compare_given_met_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9e.3b's fixtures). `capture.json` records the Rust revision
(`68f51519`, the evidence commit) and that the kernel's S9e.3b code did not
exist (`rust_given_met_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/triple.rs` absent). The kernel's probe
reports every case `unsupported`: its first Boolean evaluates and the
Boolean given the result is refused with `OutOfDomain` by
`curved/meet.rs`'s `edge_surface` ("a given result's meeting of two curved
faces met by another face (S9e.3b)"), on all 50. `inputs.txt` holds the 50
cases of `boolean-given-met-cases.txt` as explicit constructions
(`identity_reference.native_chained_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame
(`XY`, `SIDE`, `TILT` and `TILTX`), each sphere, cone and torus its
`BRepPrimAPI` primitive.

## Observations

* Every case is done, every result valid, and the solid counts are the
  reference's in all 50 (`chained_curved_boolean_reference.py`, its counts
  the declared ones by rays at two resolutions; `torus_cut_pipe_common` two
  solids, the others one). No result reports warnings.
* 12 match within the 2e-8 allowance. 38 are reviewed
  (`occt-boolean-given-met-divergences.json`): BRepGProp's default
  integration on faces bounded by the Boolean's approximated sections (the
  given meetings of two curved faces and the partner's meetings with them)
  misses by up to 4.8e-6 (`peg_pipe_common`'s area); the same native results
  measured with an adaptive accuracy of 1e-10 (a diagnostic build, not the
  capture's probe) are within 3.2e-8 of the reference, unchanged at 1e-12:
  what remains is the approximated sections themselves, on the smallest
  results (`cross_ball_common`, `cone_cut_ball_common`).
* The declared degenerate fixtures are valid results natively: a box whose
  bottom face touches the peg's meeting with the sphere at its lowest point
  (`peg_touch`) and a ball through that point with its normal dependent on
  the two surfaces' (`peg_kiss`), their measures the reference's
  (`peg_touch_common` reviewed with the others).
* Unified counts: 12 results lose faces, edges or vertices when unified
  (OCCT's split faces and seams joined); the others keep their counts.
