# Spline models of the certified quadrature (F8): native observations before the kernel integrates them by it

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code integrated mass properties by the
certified quadrature of F8 (REVIEW_NOTES.md): the kernel then enclosed these
models by the first-order strips and Green integrals of S4d. The inputs
(`inputs.txt`) are `generate_brep_fixtures.quadrature_models()`, three
spline models whose routes the S4 models do not exercise; `oracle.cpp` is
`occt_brep_check_oracle.cpp` as captured, run with `BREP_ORACLE_PROPERTIES`
set, and `native.txt` its four rows per case: the `BRepCheck_Analyzer`
statuses, the subshape counts, the tolerances and measured deviations, and
`BRepGProp::VolumeProperties` and `SurfaceProperties` at a requested
relative precision of `1e-12` (volume, its relative error estimate, area,
its relative error estimate, centre of mass, matrix of inertia about it,
row-major).

* `spline_stadium_parallel`: the stadium of half width `0.25` (tolerance
  `1e-6`) whose first cylinder wall's top pcurve, along the parallel
  `v = 1`, is a degree-1 spline: the Green integral along a spline pcurve
  with `du ≠ 0`.
* `spline_bulge_split_wall`: the bulge whose ruled spline wall has a knot
  at half its height, so the wall's upper patch has one below it in its
  column and the vertical pcurves cross the knot line.
* `spline_rounded_corner_split_wall`: the same for the rational corner.

## Observations

* All three are valid natively (`R 1`); the counts are those of the models
  they vary.
* The volumes and areas are the unsplit models' (`1.5 + π/16` and
  `9 + 5π/8` for the stadium, whose arcs are binary64 half turns). OCCT's
  area of the split bulge is again `23.6289205238958`, `1.7e-9` relative
  above the closed form `40/3 + 8 + √2 + asinh 1` while estimating `2e-16`,
  and its volume of the split rounded corner again `5.7e-12` relative
  below the reference, as in `occt-spline-properties`.
* The stadium's `I_yy` is `0.17568042331343026` where the reference's
  quadrature gives `0.17568042331321812` (`1.2e-12` relative), again
  beyond OCCT's own estimate (`5e-16`).

The bridge (`compare_brep.py --family spline`) requires these rows to
reproduce as the S4d rows do, and every kernel enclosure to contain OCCT's
value within OCCT's error estimate plus `1e-8` relative.
