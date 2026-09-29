# Prisms of spline profiles split by a plane (S8b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28. `capture.json`
records the Rust revision and that no kernel code for the spline split existed
(`rust/kernel/src/solid/split/spline.rs` absent; the kernel's profiles hold
no spline segments yet). `inputs.txt` holds the 17 cases of
`split-spline-cases.txt` as explicit constructions
(`identity_reference.native_case`: the profile face at the start offset, a
spline segment a `Geom_BSplineCurve` edge through the lifted poles, with the
knots and multiplicities as given, between the path's vertices; its prism
vector) with the splitting plane; `oracle.cpp` is `occt_split_oracle.cpp` as
captured, and `native.txt` its solids: side, volume, area, centre, face, edge
and vertex counts, validity.

The profiles: a rectangle with a quadratic bulge, a blob of four cubics
(given clockwise), a rectangle under a quadratic wave of three spans
(interior knots 1 and 2), a stadium whose left end is a cubic of two spans
(a double interior knot) beside the right end's arc, and a square with a
lens-shaped hole of two cubics (given clockwise). The planes: normal to the
axis, parallel to it (crossing a spline once, twice and four times, through
the wave's interior knot's point, through two joins of the blob, tangent to
the bulge's extreme, beyond the bulge inside its control polygon), oblique
(across spline walls and caps, across every wall and no cap, across a spline
beside an arc, touching a join only) and in a cap, in both frames.

## Observations

* Every piece is valid, and the sides agree with the reference
  (`split_reference.py`); the tangent plane, the plane inside the control
  polygon only, the plane touching a join and the plane in a cap leave one
  solid (the tangent plane with an extra edge along the touching ruling: 6
  faces, 14 edges, 10 vertices); the wave cut above its dip gives one piece
  below and two above.
* 10 cases agree with the reference within 2e-8, most to 1e-15 (planes
  parallel to the axis cut the spline walls along rulings).
* Where a plane cuts a single-span spline wall across its rulings
  (`bulge_height`, `bulge_oblique`, `blob_walls_only`, `capsule_steep`),
  BRepGProp's default integration of the trimmed extrusion faces errs by
  3.6e-8 to 2.6e-7 relative; with an adaptive accuracy of 1e-13 the same
  pieces agree with the reference within 1e-10.
* BRepGProp integrates faces bounded by the wave (three spans) badly: the
  uncut prism's volume is 329.093 by default and 329.623 adaptively, not
  328.75, though Green's theorem over OCCT's own edges of the profile gives
  65.7500000397 (the reference's 65.75) and the same curve given as three
  Bezier spans gives the same numbers. Pieces holding more than one span
  (`wave_knot`, `wave_parallel_four`, `wave_tilted`) err by 5e-6 to 1.8e-3;
  those within one span agree to 1e-15.
* All seven are reviewed native differences in
  `occt-split-divergences.json`; the adaptive and Green's-theorem figures
  came from a scratch variant of the oracle, not from the capture.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record: `platform-linux/`, from CI run 36556520655). Until
the kernel reads spline segments the probe cannot parse the `B` rows and the
comparison lists every case under `rust_unsupported` and `rust_probe_failed`.
