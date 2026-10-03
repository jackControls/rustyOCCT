# Booleans of spline prisms against prisms with arc or circle walls on crossing axes (S9f.2b): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, first captured 2026-10-03 by
`compare_spline_crossing_boolean.py --capture` (`compare_boolean.py`'s
probe and protocol on S9f.2b's fixtures) while the kernel's S9f.2b code did
not exist (`rust/kernel/src/solid/boolean/curved/spline_crossing.rs`
absent; the probe reported every case `unsupported`: the curved engine's
`spline_pairs` refused all 34 with `OutOfDomain("a spline prism against a
prism with arcs on crossing axes (S9f.2b)")`, the message read from
`protocol::run` by a throwaway example). No protocol or oracle change was
needed: the protocol's `B` path segments and the native `wire S` rows
already carry a spline's poles lifted by its frame and an arc's circle (S5,
S8b, S9a.2), and the oracle builds each spline as a `Geom_BSplineCurve`
edge and each prism by `BRepPrimAPI_MakePrism` of its profile face in its
own frame.

Taken again, whole, on 2026-10-03 for S9f.2b.2's evidence (17 cases added,
51 in all) before S9f.2b.2's code: `capture.json` records the Rust revision
(the extended fixtures, `2db283c4`) and that S9f.2b.2's code did not exist
(`rust_spline_loops_boolean_exists` false: `spline_crossing.rs` still
refusing a turning point inside both faces and a cap circle in a plane
holding the wall's axis with its two `OutOfDomain(... (S9f.2b.2))`, every
S9f.2b.2 case `unsupported`, the 28 S9f.2b.1 cases the kernel's as before).
The 34 earlier rows are the first capture's to the bit. `inputs.txt` holds
the 51 cases of `boolean-spline-crossing-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured (unchanged) and `native.txt` its
results. Frames: the spline prisms in `XY`, their partners in `TILT`,
`LEAN`, `SIDE` and `STEEP` (the normal (0, 5, 12) / 13), the steep cylinder
the object against the blob.

## Observations

* Every case is done, every result and solid valid with no warnings, and
  every solid count is the reference's (`curved_boolean_reference.py` with
  S9f.2b's crossing walls; `dome_touch_fuse` two solids touching at the
  dome's apex, `dome_touch_common` empty).
* 13 cases match within the 2e-8 allowance, the declared degenerate
  `dome_touch` among them.
* 21 cases are reviewed (`occt-boolean-spline-crossing-divergences.json`),
  all `measure` (and `centre` for 16): BRepGProp's default integration on
  the faces bounded by the Boolean's B-spline intersection edges misses by
  up to 1.5e-6 in volume, 7.9e-7 in area and 1.7e-7 of the size in the
  centre (`dome_side`, `capsule_tilt`, `steep_blob`, `bulge_side_loop`,
  `lens_tilt_loop`, `bulge_tilt_common`, `knot_turn_common`), by 1.3e-5 on
  `knot_turn_cut` and by up to 2.9e-4 on the waves (`wave_lean`, the
  quadratic B-spline of three spans whose faces S8b's, S9a.2's, S9f.1's
  and S9f.2a's captures reviewed). A diagnostic build of the oracle (not
  committed) measures each solid twice more: BRepGProp with an adaptive
  accuracy of 1e-10, within 2.6e-9 in volume but for the waves (2.3e-3,
  no better) and `knot_turn_cut` (1.5e-4), its surface integration erratic
  on some faces (up to 1.3e-5 in area); and Green's theorem over OCCT's own
  faces and pcurves (the inner integral along each surface's `v` by
  6-point Gauss-Legendre, exact for its planes, cylinders and linear
  extrusions, the outer along each pcurve by 2,000 subintervals of 5-point
  Gauss-Legendre), within 4.0e-9 in volume, 1.4e-9 in area and 8.6e-9 of
  the size in the centre on every reviewed case but `bulge_tilt_common`
  (4.6e-8 in volume and 2.0e-8 in the centre, a sliver of volume 0.5 where
  the approximated intersection edges' own deviation shows) and
  `lens_tilt_loop_common` (5.1e-9 in area, 2.2e-8 in the centre). The better of the two is within 4.0e-9 in
  volume, 1.4e-9 in area and 8.6e-9 in the centre on all 21: the native
  geometry is the reference's; its default measure is not.
* The degenerate fixtures are valid results natively: `dome_touch` (a rod
  tangent to the dome's wall at its apex, a point contact; its fuse two
  solids, its common empty) and `knot_turn` (a rod whose meeting with
  `knot`'s wall turns back at its knot of multiplicity two). The
  decisions refuse all four (`Degenerate`).
* Five solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `capsule_tilt_cut` 9/21/14 to
  8/20/14, `capsule_tilt_common` 7/15/10 to 6/14/10, `ring_dome_fuse` and
  `ring_dome_cut` 12/30/20 to 11/27/18, `lens_tilt_loop_fuse` 9/15/10 to
  8/15/10 (a plane's, a cylinder's or a spline wall's faces split by the
  Boolean merged). The kernel's results will be compared with the unified
  counts, as S9c's.

## S9f.2b.2's cases (the recapture)

* All 17 done, every result and solid valid with no warnings, every solid
  count the reference's: towers alone (`dome_cap_tower`, a rod ending
  inside the dome, its cap circle on the arch's generatrix at `x = 5/2`;
  `lens_cap_tower`, a wide rod ending inside the lens, its cap circle on
  both cubics' generatrices), loops cut by caps (`bulge_cap_loop`, a tower
  cap circle near the bulge's apex; `bulge_side_cap_loop`, the bulge's top
  cap; `lens_tilt_cap_loop`, a tilted rod's cap circle), and the declared
  degenerate `cap_turn` (a turning point on the rod's cap rim, valid
  natively; the decisions refuse it).
* 4 match within the 2e-8 allowance (`bulge_cap_loop` and `cap_turn_cut`),
  13 are reviewed, all `measure` (and `centre` for 12): BRepGProp's default
  integration misses by up to 2.0e-6 in volume (`lens_tilt_cap_loop_common`,
  a sliver of volume 0.35), 9.7e-7 in area and 3.5e-7 of the size in the
  centre. The diagnostic build's adaptive BRepGProp (up to 1.7e-5 off in
  area on some faces, as before) and Green's theorem over OCCT's own faces
  and pcurves: the better of the two within 2.1e-9 in volume, 3.7e-9 in
  area and 1.8e-10 in the centre on all 13: the native geometry is the
  reference's; its default measure is not.
* Five more solids' counts change when unified: `bulge_cap_loop_fuse`
  13/24/16 to 12/22/15, `bulge_side_cap_loop_cut` 9/21/14 to 8/19/13,
  `bulge_side_cap_loop_common` 5/9/6 to 4/7/5, `lens_tilt_cap_loop_cut`
  7/12/8 to 6/10/7, `lens_tilt_cap_loop_common` 4/6/4 to 3/4/3.
* The reference itself is checked independently
  (`generate_spline_crossing_boolean_fixtures.py --check`): every operation
  two ways (the slicing, and the divergence theorem over the classified
  face pieces) within 2.0e-41 of the case's size, inclusion and exclusion
  2.0e-41, the area identity 9.4e-41, every face's classes 4.5e-41; on the
  four perpendicular pairs the common as the product of the profile's
  chord and the disc's height chord within 7.1e-44; margins outside the
  declared pairs at least 0.011 (vertical edges piercing curved walls at a
  sine's complement of 0.011, vertices 0.037 from the other's faces,
  creases 0.073 from the caps, meetings crossing caps at 0.11, turning
  points 0.56 outside a face, the loops' 1.5 inside both), the declared
  pairs' tangency gap and turning point's distance from the knot zero.
  Python 3.9 and 3.12 write the same files.

`compare_spline_crossing_boolean.py` requires every later run to reproduce
these rows (on another platform, its reviewed record; the Linux record is
pending CI). While S9f.2b.2's two refusals stand in `spline_crossing.rs`
its 23 cases may stay `unsupported` (listed under `rust_unsupported`);
once they are gone an `unsupported` S9f.2b.2 case is a failure
(`rust_unsupported_after_its_code`).
