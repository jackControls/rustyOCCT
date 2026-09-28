# STEP files read natively (STEP-b): native observations before the importer

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks and the data exchange toolkits
(`build_pinned_occt.py --toolkit TKDESTEP`). `capture.json` records the Rust
revision and that the kernel's translation of STEP-b's entities did not
exist (`rust/kernel/src/step/spline.rs` absent): the kernel's importer of
that revision rejected all seven files by their first construct outside
STEP-a (`ELLIPSE` twice, `B_SPLINE_CURVE_WITH_KNOTS` twice,
`B_SPLINE_SURFACE_WITH_KNOTS` twice, `RATIONAL_B_SPLINE_SURFACE`).
`inputs.txt` holds STEP-b's seven fixture files of `fixtures/step`
(`generate_step_fixtures.py`, `STEP_B`), framed as for STEP-a; `oracle.cpp`
is `occt_step_oracle.cpp` as captured, unchanged since STEP-a, and
`native.txt` its rows: class, OCCT's counts, BRepCheck's verdict, volume,
area, centre and largest tolerance. STEP-a's capture
(`occt-step-preimplementation`) is unchanged and still reproduced.

## Observations

* Every file reads and every body is valid, with the independent
  reference's class and counts (`step-expected.tsv`); every tolerance stays
  `1e-7` mm.
* Measures, against the reference's closed forms, exact integrals and
  quadratures: the half-ellipse sheet and the B-spline plate agree within
  `2.2e-16` relative; the B-spline prism's area within `2.2e-12`, the
  trimmed patch's area and centre within `2.0e-10` (both inside the bridge's
  `1e-9`). Three bodies differ beyond it, all through `BRepGProp`'s default
  integration (fixed Gauss orders), not the reader: the oblique cylinder
  (volume `3.5e-6`, area `3.3e-6` relative, centre `1.5e-7` of the size),
  the two-span patch (area `5.7e-7`, centre `3.2e-7`) and the rational
  cylinder (volume `2.0e-3`, area `8.2e-4`). A scratch probe of the same
  shapes with `BRepGProp`'s adaptive integration (`Eps = 1e-12`, not part of
  the capture) gives the oblique cylinder's volume within `4e-16` of
  `250 pi`, the rational cylinder's volume and area within `6.5e-12` of
  `200 pi` and `130 pi`, and the patch's area within `3.8e-13` of the
  quadrature.

`compare_step.py` requires every later run to reproduce these rows and
compares the kernel's bodies with both the reference and these
observations; the three differences above need fingerprinted reviews.
