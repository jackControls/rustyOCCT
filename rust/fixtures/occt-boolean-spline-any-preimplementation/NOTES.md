# Booleans of spline prisms against polyhedral prisms in any position (S9f.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_spline_any_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9f.1's fixtures). `capture.json` records the Rust
revision and that the kernel's S9f.1 code did not exist
(`rust_spline_any_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/spline_walls.rs` absent). The kernel's
probe reports every case `unsupported`: 32 cases are refused by the
polyhedral engine with `OutOfDomain("a Boolean of a solid with curved faces
or edges in any position (S9c)")` and the 6 `capsule_turned` and
`capsule_stand` cases (the capsule's profile also holds an arc, so they
reach the curved engine) by its model with `OutOfDomain("a spline profile
in a Boolean of prisms in any position (S9c)")`, the messages read from
`protocol::run` by a throwaway example. No protocol or oracle change was
needed: the protocol's `B` path segments and the native `wire S` rows
already carry a spline's poles lifted by its frame (S8b, S9a.2), and the
oracle builds each as a `Geom_BSplineCurve` edge and each prism by
`BRepPrimAPI_MakePrism` of its profile face in its own frame. `inputs.txt`
holds the 38 cases of `boolean-spline-any-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results. Frames:
spline prisms in `XY` and `TILT`; the polyhedral prisms in `XY`, `SIDE`,
`TURN`, `R125`, `LEAN`, `TILT`, `TILT2` and `TILTX`.

## Observations

* Every case is done, every result and solid valid with no warnings, and
  every solid count is the reference's
  (`curved_boolean_reference.py` with S9f.1's spline walls; `lens_pin_common`
  and `kink_lean_cut` two solids, `dome_tangent_fuse` two touching along
  the dome's apex generatrix, `capsule_stand_common` and
  `dome_tangent_common` empty).
* 23 cases match within the 2e-8 allowance (at worst 1.7e-8 in volume,
  1.8e-8 in area, 6.9e-9 of the size in the centre), the declared
  degenerate ones and the coplanar caps among them (`capsule_stand`, a
  face shared with the opposite orientation, within 1.6e-15).
* 15 cases are reviewed (`occt-boolean-spline-any-divergences.json`), all
  `measure` (and `centre` for 10): BRepGProp's default integration on the
  faces bounded by the Boolean's B-spline edges misses by up to 6.7e-7 in
  volume and area and 2.4e-7 of the size in the centre (`dome_side`,
  `lens_pin`, `dome_turned`, `kink_lean`, `blob_tool`), and by 1.0e-3,
  4.8e-4 and 3.1e-4 on the wave (`wave_tiltx`, the quadratic B-spline of
  three spans whose faces S8b's and S9a.2's captures reviewed). With an
  adaptive accuracy of 1e-10 (a diagnostic build of the oracle, not
  committed; unchanged at 1e-12) the volumes are within 6.0e-9 and the
  centres 2.6e-9 but for the wave, where BRepGProp does not converge (2.5e-3
  at 1e-10, 2.6e-3 at 1e-12); the adaptive surface integration is erratic
  on these faces (up to 2.2e-5 at 1e-10). The same diagnostic build
  integrating each native solid by Green's theorem over OCCT's own faces
  and pcurves (the inner integral along each surface's `v` by 6-point
  Gauss-Legendre, exact for its planes and linear extrusions, the outer
  along each pcurve by 2,000 subintervals of 5-point Gauss-Legendre) gives
  the reference within 1.3e-8 in volume, 6.8e-9 in area and 5.8e-9 of the
  size in the centre on all 38 cases (the wave within 3.2e-9): the native
  geometry is the reference's; its default measure is not.
* The degenerate fixtures are valid results natively: `dome_tangent` (a
  plane tangent to the dome's wall along its apex generatrix; its fuse two
  solids, its common empty), `kink_knot_cut` (a plane tangent to the kink's
  wall along the generatrix of its knot of multiplicity two) and
  `blob_rounding` (a `TILTX` cap plane within rounding, 8.9e-17, of the
  `TILT` blob's axis). The decisions refuse all five (`Degenerate`).
* Two solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `capsule_turned_fuse` 17/39/24 to
  13/33/22 and `capsule_stand_cut` 8/19/13 to 6/15/11 (faces of one plane
  or one cylinder split at the coplanar caps' imprints merged). The
  kernel's results will be compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spline_any_boolean_fixtures.py --check`): every operation two
  ways (the slicing, and the divergence theorem over the classified face
  pieces) within 4.2e-41 of the case's size, inclusion and exclusion
  4.5e-41, the area identity 1.2e-40, every face's classes 3.5e-41, a
  second slicing direction for parallel axes 2.1e-41; each spline
  profile's prism against a half-space box, oblique (creases) and parallel
  to its axis (generatrices), as S8b's split reference within 1.4e-41 in
  exact frames and 1.9e-18 in turned ones; the rounding-offset pair
  `bulge_offset` and S9a.2's own 17 one-spline fixtures as S9a.2's
  `SplinePair` within 2.2e-17 (the 10 in `XY` to every printed digit); its
  margins outside the declared pairs at least 0.019 (vertices from the
  other's faces), 0.12 (creases' extreme heights from the caps), 0.16
  (parallel planes from tangency to a span), 0.37 (sines of edges crossing
  spline walls) and 0.48 (sines of oblique planes to the axis). Python 3.9
  and 3.12 write the same files.

`compare_spline_any_boolean.py` requires every later run to reproduce these
rows (on another platform, its reviewed record; the Linux record is
pending CI). While `spline_walls.rs` is absent the probe
(`rust/kernel/examples/boolean_probe.rs`) must report `unsupported` on all
38 cases, which the comparison lists under `rust_unsupported`; any other
row, or a probe failure, is a failure.
