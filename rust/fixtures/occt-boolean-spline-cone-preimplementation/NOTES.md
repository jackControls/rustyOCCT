# Booleans of spline prisms against cones and frustums (S9f.3b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-10-04 by
`compare_spline_cone_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9f.3b's fixtures) while the kernel's S9f.3b code did not
exist (`rust/kernel/src/solid/boolean/curved/spline_cone.rs` absent; the
probe reported every case `unsupported`: the curved engine's
`spline_pairs` refused all 27 with `OutOfDomain("a spline prism against a
cone (S9f.3b)")`). No protocol or oracle change was needed: the protocol's
`B` path segments and the native `wire S` rows carry a spline's poles
lifted by its frame (S8b, S9a.2), and the `cone` row builds a cone or a
frustum by `BRepPrimAPI_MakeCone` on its frame between its radii (S9d.3a).
`capture.json` records the Rust revision (the fixtures, `a5818a39`) and
`rust_spline_cone_boolean_exists` false. `inputs.txt` holds the 27 cases of
`boolean-spline-cone-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured (unchanged) and `native.txt` its
results. Frames: the spline prisms in `XY` and `TILT`, the cones in exact
frames (`XY`, `DOWN`, `SIDE` and one along `y`) about binary64 points.

## Observations

* Every case is done, every result and solid valid, and every solid count
  is the reference's (`spline_cone_boolean_reference.py`). The declared
  degenerate `apex_wall` (the cone's apex on the dome's arch) is valid but
  reported with warnings on both its operations; `ruling_tilt` (`A = 0`)
  and `dome_touch` (a tangency) without. The decisions refuse all six
  (`Degenerate`).
* 4 cases match within the 2e-8 allowance (`bulge_frustum`'s three, a
  frustum on the bulge's axis, and the declared degenerate
  `ruling_tilt_cut`).
* 23 cases are reviewed (`occt-boolean-spline-cone-divergences.json`), all
  `measure` (and `centre` for 15): BRepGProp's default integration on the
  faces bounded by the Boolean's intersection edges misses by up to 5.5e-6
  in volume, 5.0e-6 in area and 9.5e-7 of the size in the centre on the
  quadratic and cubic walls (`apex_wall_common`), by up to 7.4e-5 on
  `knot` in the tilted frame and by up to 1.0e-3 on the wave (`cone_wave_fuse`, the quadratic
  B-spline of three spans every earlier spline capture reviewed). A
  diagnostic build of the oracle (not committed) measures each solid twice
  more: BRepGProp with an adaptive accuracy of 1e-10 (within 2.0e-8 in
  volume on the quadratic and cubic walls, but up to 8.0e-4 on `knot` and
  1.1e-3 on the wave, its surface integration up to 1.7e-5 off on faces
  bounded by the meeting); and Green's theorem over OCCT's own faces and
  pcurves (the inner integral along each surface's `v` by 6-point
  Gauss-Legendre, exact for its planes, cones and linear extrusions, the
  outer along each pcurve by 2,000 subintervals of 5-point Gauss-Legendre).
  The better of the two is within 3.2e-8 in volume, 7.0e-9 in area and
  9.0e-9 of the size in the centre on every reviewed case (the largest
  volume and centre `knot_tilt_cone_common`'s, the area
  `ruling_tilt_common`'s, their intersection edges approximated on the
  splines' spans): the native geometry is the
  reference's; its default measure is not.
* Three solids' counts change under `ShapeUpgrade_UnifySameDomain`
  (faces, edges, vertices as built, then unified): `cone_wave_common`
  5/9/6 to 4/6/4, `knot_tilt_cone_cut` 9/21/14 to 8/20/14 and
  `knot_tilt_cone_common` 5/9/6 to 4/8/6. The kernel's results will be
  compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spline_cone_boolean_fixtures.py --check`): each region's
  volume and first moments by two slicings in directions cutting the cone
  in ellipses within 1.1e-39 of the case's size, both inputs against their
  closed forms within 1.0e-39, every face's classes within 1.4e-40 of its
  area (the cone's wall by its area element against `pi (b + t)` times its
  slant), the area identity within 9.2e-41, Monte Carlo within 2.1 standard
  errors; margins outside the declared pairs at least 0.0156 (`|A|` of the
  frustum on the bulge's axis), the declared pairs' apex distance, `|A|`
  and tangency gap zero. Python 3.9 and 3.12 write the same files.

`compare_spline_cone_boolean.py` requires every later run to reproduce
these rows (on another platform, its reviewed record; the Linux record is
pending CI). While `spline_cone.rs` is absent every case must be
`unsupported`; once it exists none may be
(`rust_unsupported_after_its_code`).
