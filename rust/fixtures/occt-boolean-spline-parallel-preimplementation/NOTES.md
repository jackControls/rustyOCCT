# Booleans of spline prisms against prisms with arc, circle or spline walls on exactly parallel axes (S9f.2a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-10-03 by
`compare_spline_parallel_boolean.py --capture` (`compare_boolean.py`'s
probe and protocol on S9f.2a's fixtures). `capture.json` records the Rust
revision and that the kernel's S9f.2a code did not exist
(`rust_spline_parallel_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/spline_parallel.rs` absent). The
kernel's probe reports every case `unsupported`: the curved engine's
`spline_pairs` refuses the 20 cases of a spline prism against a prism with
arcs or circles with `OutOfDomain("a spline prism against a prism with arcs
in any position (S9f.2)")` and the 23 of two spline prisms with
`OutOfDomain("spline walls against spline walls in any position (S9f.2)")`,
the messages read from `protocol::run` by a throwaway example. No protocol
or oracle change was needed: the protocol's `B` path segments and the
native `wire S` rows already carry a spline's poles lifted by its frame and
an arc's circle (S5, S8b, S9a.2), and the oracle builds each spline as a
`Geom_BSplineCurve` edge and each prism by `BRepPrimAPI_MakePrism` of its
profile face in its own frame. `inputs.txt` holds the 43 cases of
`boolean-spline-parallel-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results.
Frames: the spline prisms in `XY` and `TILT`, their partners in `R125`,
`TURN`, `FLIP` (a half turn about the `XY` normal), `TILT2` (about
`TILT`'s axis) and `TILT` (rounding offsets).

## Observations

* Every case is done, every result and solid valid with no warnings, and
  every solid count is the reference's (`curved_boolean_reference.py` with
  S9f.2a's parallel walls; `wave_wave_flip_common` two solids,
  `dome_disc_touch_fuse` and `dome_flip_touch_fuse` two touching along the
  dome's apex generatrix, `dome_disc_touch_common` empty).
* 31 cases match within the 2e-8 allowance (at worst 1.8e-8 in volume,
  6.5e-9 in area, 1.2e-9 of the size in the centre), the declared
  degenerate ones among them.
* 12 cases are reviewed (`occt-boolean-spline-parallel-divergences.json`),
  all `measure` (and `centre` for the waves): BRepGProp's default
  integration on the faces bounded by the Boolean's B-spline edges misses
  by up to 1.5e-7 in volume, 1.9e-7 in area and 8.7e-9 of the size in the
  centre (`bulge_dome_turn`, `hump_lens_turn`, `capsule_disc_turn_common`,
  `blob_lens_r125_common`, `bulge_disc_offset_common`), and by 8.6e-4,
  2.4e-4 and 2.0e-4 on the waves (`wave_wave_flip`, the quadratic B-spline
  of three spans whose faces S8b's, S9a.2's and S9f.1's captures
  reviewed). With an adaptive accuracy of 1e-10 (a diagnostic build of the
  oracle, not committed) the volumes are within 2.3e-9 and the centres
  3.1e-10 but for the waves, where BRepGProp does not converge (6.3e-4);
  the adaptive surface integration is erratic on these faces (up to
  6.9e-6). The same diagnostic build integrating each native solid by
  Green's theorem over OCCT's own faces and pcurves (the inner integral
  along each surface's `v` by 6-point Gauss-Legendre, exact for its planes,
  cylinders and linear extrusions, the outer along each pcurve by 2,000
  subintervals of 5-point Gauss-Legendre) gives the reference within
  9.5e-9 in volume, 1.0e-8 in area and 1.1e-8 of the size in the centre on
  all 39 cases with a result (the waves within 1.1e-8): the native
  geometry is the reference's; its default measure is not.
* The degenerate fixtures are valid results natively: `dome_disc_touch` (a
  disc's cylinder tangent to the dome's wall along its apex generatrix;
  its fuse two solids, its common empty) and `dome_flip_touch` (two domes,
  one turned half a turn, tangent along their apexes' generatrix; the
  fuse two solids). The decisions refuse all four (`Degenerate`).
* Four solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `blob_ring_turn_cut` 13/33/22
  to 12/30/20, `blob_ring_turn_common` 9/21/14 to 8/18/12,
  `capsule_disc_turn_fuse` and `bulge_disc_offset_fuse` 10/21/14 to
  9/18/12 (a cylinder's or a plane's faces split by the Boolean merged).
  The kernel's results will be compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_spline_parallel_boolean_fixtures.py --check`): every
  operation two ways (the slicing, and the divergence theorem over the
  classified face pieces) within 1.3e-40 of the case's size, inclusion and
  exclusion 1.3e-40, the area identity 9.8e-41, every face's classes
  3.0e-41, a second slicing direction 1.3e-40; S9a.2's `SplinePair` on the
  33 cases whose map between the frames is an exact turn (`TURN`, `FLIP`,
  `TILT2` against `TILT`) or an offset in one frame (the tool's profile
  turned exactly into the object's frame, the offset taken exactly)
  within 2.4e-39 in `XY` and 2.2e-17 in `TILT`; its margins outside the
  declared pairs at least 0.06 (crossings' sines 0.47, near misses 0.37,
  vertices from the other's faces 0.13, creases' heights from the caps
  0.06), the declared pairs' crossing sines below 5.7e-21. Python 3.9 and
  3.12 write the same files.

`compare_spline_parallel_boolean.py` requires every later run to reproduce
these rows (on another platform, its reviewed record; the Linux record is
pending CI). While `spline_parallel.rs` is absent the probe
(`rust/kernel/examples/boolean_probe.rs`) must report `unsupported` on all
43 cases, which the comparison lists under `rust_unsupported`; any other
row, or a probe failure, is a failure.
