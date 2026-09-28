# Tessellation of splines (T-b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks and TKMesh (the SDK of T-a's
capture, the same manifest digest). `capture.json` records the Rust
revision (the reference and fixtures committed, `394f1491`) and that no
kernel code tessellated spline geometry (`rust/kernel/src/tessellation/spline.rs`
absent; `tessellate` still returned `OutOfDomain` for every spline cell).
`inputs.txt` holds the twelve bodies of `tessellation-spline-cases.txt` as
the kernel's existing `.brep` writer produces them (`brep_io_probe parts`,
records 7 and 9 for the splines), then the 24 mesh rows, a coarse (a
hundredth of the case's scale, 0.5 rad) and a fine setting (a thousandth,
0.3 rad) per body; `oracle.cpp` is T-a's `occt_tessellation_oracle.cpp`,
unchanged (it already projects onto B-spline surfaces and curves with
`GeomAPI_ProjectPointOnSurf` and `GeomAPI_ProjectPointOnCurve`), and
`native.txt` its `S` rows.

## Observations

Checked against the independent reference (`tessellation_reference.py`) on
the welded native meshes (the probe's `dump` mode, not part of the
capture): its largest sampled distances equal the probe's own projections
within `1.6e-16` of each case's scale on every row but `spline_bulge_far`'s
(`6.5e-14`, its coordinates two thousand times its scale), which checks the
reference's projections as much as OCCT.

* Every solid is watertight once nodes are joined through OCCT's edge
  polygons (no free, non-manifold or misoriented mesh edge), with the
  boundary's Euler characteristic (`spline_wave_holes` -2, the others 2);
  each face body's free edges are exactly its boundary loops. No triangle
  is degenerate.
* All eight spline prisms, the spline ring face and the rational sheet stay
  within the requested deflection at both settings, and the dome and the
  trimmed sheet at the coarse one.
* The requested deflection is exceeded on two rows: `spline_dome` fine
  (`0.00385` against `0.003`) and `sheet_spline_hole` fine (`0.00630`
  against `0.005375`). On the dome OCCT's own
  `Poly_Triangulation::Deflection` understates its measured distance at both
  settings (`0.02829` against `0.02986` coarse, `0.00377` against
  `0.00385` fine).
* OCCT's recorded deflections exceed the request on five rows whose
  measured distances do not (`spline_sharp_rational` records `0.266` for
  `0.0225` measured, the sheets up to `0.087`): its per-face value is not
  the distance of the final mesh.

`compare_tessellation.py --family spline` requires every later run to
reproduce these rows and the `.brep` texts, reviews OCCT's differences
with fingerprints (`occt-tessellation-divergences.json`), and requires the
kernel's meshes to pass every check of the reference.
