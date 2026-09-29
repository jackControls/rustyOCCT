# Booleans of torus v-segments and wedges against polyhedral prisms (S9d.4b.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_torus_segment_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9d.4b.1's fixtures). `capture.json` records the Rust
revision (the fixtures and the reference committed, nothing else changed)
and that the kernel's S9d.4b.1 code did not exist
(`rust_torus_segment_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/torus_segment.rs` absent). The
kernel's probe reports every case `unsupported`: S9d.4a's torus model
refuses anything but a whole torus with `OutOfDomain("a Boolean of a torus
segment or wedge (S9d.4b)")` (all 29 cases, the message read from
`protocol::run` by a throwaway example; every input constructs, so the
kernel's test support reads the rows). `inputs.txt` holds the 29 cases of
`boolean-torus-segment-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`: a part's block one `torus` row,
its frame's origin, normal and x axis, major and minor radii, latitudes and
turn), `oracle.cpp` is `occt_boolean_oracle.cpp` as captured (its `torus`
row taking the latitudes and turn new with this capture; a whole torus's
row and every other row unchanged, every older capture reproducing) and
`native.txt` its results; each part is `BRepPrimAPI_MakeTorus(gp_Ax2, R, r,
a1, a2, angle)`, reversed when its BRepGProp volume is negative (OCCT builds
a v-segment whose meridian ends lower than it starts inside out: S3's
capture, here the inner half), each prism `BRepPrimAPI_MakePrism` of its
profile face in its own frame. Frames: `XY` throughout, but the outer half
in `TILT` and the half turn in `TILT2` against prisms in `XY`.

## Observations

* Every case is done without warnings, every result and solid valid. The
  solid counts are the reference's in all 29
  (`torus_segment_boolean_reference.py`: the meridian half-planes swept,
  each section's components by bands of height): the spool cut by a slab
  through its waist two solids, the column less the spool six (its ends
  above and below and four corners through the waist), a bar along the half
  turn's end discs two chunks (`hw_bar_common`), a box beyond the quarter
  turn fused two solids and a box against its start disc one
  (`qw_back_fuse`, the common empty), the tilted half turn less a slab two
  (`hw_tilt_cut`); the degenerate `ub_top_fuse` and `qw_touch_fuse` two
  solids touching along the top circle and at a point, `oh_lid_fuse` one
  (the box on the outer half's end disc, its plane tangent to the wall along
  the disc's ring).
* 16 cases match within the 2e-8 allowance: where every section is a
  circle or a line (the slabs normal to the axis through the spool and the
  half turn, the lid and the box on the planes of the end discs, the parts
  touching or apart) within 2.6e-16; the others, whose results hold
  B-spline sections, within 1.1e-8 (`oh_bar_cut`).
* 13 cases are reviewed (`occt-boolean-torus-segment-divergences.json`),
  all `measure` (and eight `centre`): BRepGProp's default integration on
  faces bounded by the Boolean's approximated sections misses by up to
  1.7e-5 in volume, 2.2e-5 in area and 3.2e-6 of the case's size in centre
  (`hw_tilt_common`). The same native results measured by BRepGProp with an
  adaptive accuracy of 1e-10 (a diagnostic build of the oracle, not
  committed) are within 5.7e-9 in volume, 3.8e-9 in area and 1.6e-9 in
  centre, but for `oh_bar_common` (3.7e-8 in volume, 1.5e-8 in area) and
  `ub_lid_common` (3.2e-8 in volume), unchanged at an accuracy of 1e-12:
  there the approximated sections themselves (B-splines within their edges'
  tolerance, 1e-6) bound a region that far from the exact one, while the
  reference agrees with the closed form of `oh_bar` within 5.5e-40.
* Curves as built: lines and circles where a section is one (a plane normal
  to the axis cuts the wall in circles, a plane through the axis in the
  tube's circles or arcs, the end discs' planes in lines); every other
  section of a plane and the wall is a B-spline approximation, 80 of the
  results' 414 edges. The results' edges' tolerances are 1e-6 where they
  hold a B-spline and at most 1.5e-7 otherwise.
* The degenerate fixtures are valid results natively: `oh_lid_fuse` one
  solid (the box's face on the outer half's upper end disc, its plane
  tangent to the wall along the ring: the disc and the face merged),
  `ub_top_fuse` two solids touching along the band's top circle,
  `qw_touch_fuse` two solids touching at the quarter turn's outermost point
  of its start disc's circle. The decisions refuse them (`Degenerate`).
* Nine solids' counts (seven cases) change under
  `ShapeUpgrade_UnifySameDomain` (faces, edges, vertices as built, then
  unified): `oh_side_common` and `ub_side_common` 5/9/6 to 4/6/4,
  `oh_tilt_common` and `hw_tilt_cut`'s two solids 4/6/4 to 3/4/3,
  `qw_box_common` 6/12/8 to 5/10/7, `slab_hw_cut` 11/27/18 to 10/24/16,
  `ub_top_fuse`'s band 4/5/3 to 3/3/2 and box 7/13/9 to 6/12/8 (the
  toroidal face split at OCCT's seams merged). The kernel's results will be
  compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_torus_segment_boolean_fixtures.py --check`): closed forms of
  17 of the 19 pairs (the part's sections by hand along the axis: an
  aligned box by S9d.1's rectangle-in-disc antiderivatives, a wedge's box
  clipped to its sector, a half-space by circular segments, inputs meeting
  on a face or at a point by their sums) within 8.2e-41 in exact frames,
  1.3e-16 in turned ones and 4.3e-16 where a wedge's rounded end cuts the
  box (the declared degenerate pairs within 4.2e-41), every operation two
  ways (normal slices, meridian half-planes) within 6.4e-41, inclusion and
  exclusion 5.8e-41, the wall two ways 7.0e-37, the end discs two ways
  2.9e-41, the area identity and every face's classes 2.8e-40, Monte Carlo
  3.3 standard errors. Its scan for near coincidences found none outside
  the declared degenerate pairs, and every other pair keeps its vertices at
  least 0.072 from the part's faces and rims, its faces' planes at least
  0.063 from the wall's tangent planes (points of contact on the face and
  the wall), its edges' crossings of the wall or an end plane at a sine of
  at least 0.28, and the rims' crossings of its faces at a sine of at least
  0.17 (distances relative to the case's size).
