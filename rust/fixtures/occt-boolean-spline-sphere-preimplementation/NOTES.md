# Booleans of spline prisms against spheres and hemispheres (S9f.3a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-10-03 by
`compare_spline_sphere_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9f.3a's fixtures) while the kernel's S9f.3a code did not
exist (`rust/kernel/src/solid/boolean/curved/spline_sphere.rs` absent; the
probe reported every case `unsupported`: the curved engine's
`spline_pairs` refused all 33 with `OutOfDomain("a spline prism against a
sphere or a cone (S9f.3)")`). No protocol or oracle change was needed: the
protocol's `B` path segments and the native `wire S` rows carry a spline's
poles lifted by its frame (S8b, S9a.2), and the `sphere` row builds a
sphere or a hemisphere by `BRepPrimAPI_MakeSphere` on its frame between its
latitudes (S9d.1). `capture.json` records the Rust revision (the fixtures,
`576445ac`) and `rust_spline_sphere_boolean_exists` false. `inputs.txt`
holds the 33 cases of `boolean-spline-sphere-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured (unchanged) and `native.txt` its
results. Frames: the spline prisms in `XY` and `TILT`, whole spheres in
`XY` about binary64 centres, hemispheres in `XY` and on their `SIDE`.

## Observations

* Every case is done, every result and solid valid, and every solid count
  is the reference's (`spline_sphere_boolean_reference.py`;
  `dome_touch_common` empty, the sphere touching the dome's arch at a
  point). The declared degenerate `cap_turn` (a loop turning back on the
  bulge's top cap's edge) is valid but reported with warnings on both its
  operations; `dome_touch` (a point contact) and `knot_turn` (a turning
  point at `knot`'s knot of multiplicity two) without. The decisions
  refuse all six (`Degenerate`).
* 4 cases match within the 2e-8 allowance (`dome_ball_fuse`,
  `dome_ball_common` and the declared degenerate `dome_touch`'s two).
* 29 cases are reviewed (`occt-boolean-spline-sphere-divergences.json`),
  all `measure` (and `centre` for 21): BRepGProp's default integration on
  the faces bounded by the Boolean's intersection edges misses by up to
  4.1e-5 in volume, 8.8e-5 in area and 1.2e-5 of the size in the centre
  (`bulge_ball`, `knot_tilt_ball`, `knot_turn`), and by up to 1.0e-3 on
  the wave (`ball_wave`, the quadratic B-spline of three spans whose faces
  every earlier spline capture reviewed). A diagnostic build of the oracle
  (not committed) measures each solid twice more: BRepGProp with an
  adaptive accuracy of 1e-10, within 5.7e-9 in volume but for the wave
  (up to 1.1e-3) and `knot`'s walls (up to 9.2e-4), its surface
  integration erratic on some faces (up to 1.6e-5 in area); and Green's
  theorem over OCCT's own faces and pcurves (the inner integral along each
  surface's `v` by 24-point Gauss-Legendre, the outer along each pcurve by
  2,000 subintervals of 5-point Gauss-Legendre). The better of the two is
  within 6.0e-9 in volume, 7.8e-9 in area and 8.5e-9 of the size in the
  centre on every reviewed case but the wave's cut and common (9.2e-7 in
  volume and 9.3e-7 in the centre, within the result's own edge tolerance
  of 3.7e-5: its intersection edges approximated on the wave's spans): the
  native geometry is the reference's; its default measure is not.
* Two solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `lens_ball_fuse` 6/10/7 to
  5/9/7, `ball_wave_common` 3/4/3 to 2/3/3. The kernel's results will be
  compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spline_sphere_boolean_fixtures.py --check`): each region's
  volume and first moments by two slicings (along the caps' normal and
  along `(2, -3, 5)`) within 1.8e-41 of the case's size, both inputs
  against their closed forms within 5.7e-42, every face's classes within
  2.8e-41 of its area (the sphere's by Archimedes' area element), the area
  identity within 8.0e-41, Monte Carlo within 2.8 standard errors; margins
  outside the declared pairs at least 0.023 (a vertex of the prism from the
  sphere), the declared pairs' tangency gap, turning point's distance from
  the knot and loop's depth below the cap zero. Python 3.9 and 3.12 write
  the same files.

`compare_spline_sphere_boolean.py` requires every later run to reproduce
these rows (on another platform, its reviewed record; the Linux record is
pending CI). While `spline_sphere.rs` is absent every case must be
`unsupported`; once it exists none may be
(`rust_unsupported_after_its_code`).
