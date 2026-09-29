# Booleans of a whole torus against polyhedral prisms (S9d.4a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_torus_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9d.4a's fixtures). `capture.json` records the Rust revision and
that the kernel's S9d.4 code did not exist (`rust_torus_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/torus.rs` absent). The
kernel's probe reports every case `unsupported`: S9b.2's stored model
refuses every Boolean with a torus with `OutOfDomain("a Boolean of a solid
with curved faces or edges in any position (S9c)")` (all 35 cases, the
message read from `protocol::run` by a throwaway example), its wall not a
plane. `inputs.txt` holds the 35 cases of `boolean-torus-cases.txt` as
explicit constructions (`identity_reference.native_boolean_case`: a torus's
block one `torus` row, its frame's origin, normal and x axis, major and
minor radii), `oracle.cpp` is `occt_boolean_oracle.cpp` as captured (its
`torus` row new with this capture; the prism, sphere and cone rows
unchanged) and `native.txt` its results; each torus is
`BRepPrimAPI_MakeTorus(gp_Ax2, R, r)`, each prism `BRepPrimAPI_MakePrism` of
its profile face in its own frame. Frames: `XY` and `SIDE` (a prism along
`x`), tori in `TILT` and `TILTX` against prisms in `XY`, a bar in `LEAN`
against a torus in `XY`.

## Observations

* Every case is done without warnings, every result and solid valid. The
  solid counts are the reference's in all 35
  (`torus_boolean_reference.py`: components followed about the axis through
  the meridian half-planes): the two chunks of the bar's common and the
  torus less the bar two C-shaped halves, the strip less the torus three
  (the hole's piece and two ends), the torus less a band two and the slab
  less the torus two (the hole's disc and the outside), the four caps of
  the torus less a box through its hole (`frame_cut`), the four caps cut
  from the tube's inside by a box's vertical edges (`corner_caps_common`),
  the two chunks of `tiltx_band` (a slab between the saddle levels of a
  torus in `TILTX`) and of `lean_bar`, a box in the hole fused two solids
  (its common empty), the cavities of `in_tube_cut` and `in_box_cut` one
  solid each, and the degenerate `top_tangent_fuse` two solids touching
  along the top circle.
* 16 cases match within the 2e-8 allowance. Where every section is a circle
  (a band or half-space normal to the axis, a half-space through the axis,
  a box inside the tube, the torus inside a box, a box in the hole, the
  face tangent along the top circle) within 4.7e-16; `corner_caps_cut`
  within 1.1e-9 (its eight B-splines).
* 19 cases are reviewed (`occt-boolean-torus-divergences.json`), all
  `measure` (and most `centre`): BRepGProp's default integration on faces
  bounded by the Boolean's approximated spiric sections misses by up to
  6.3e-6 in volume, 4.3e-6 in area and 2.4e-6 of the case's size in centre
  (`tiltx_band_common`). The same native results measured by BRepGProp with
  an adaptive accuracy of 1e-10 (a diagnostic build of the oracle, not
  committed) are within 1.8e-9 in volume (`corner_caps_common`), 1.7e-9 in
  area (`strip_cut`) and 5.8e-10 in centre (`tiltx_band_common`), but for
  `wedge_cap_common` (1.2e-8 in volume, 5.4e-9 in area: a small cap, 0.12 of
  volume, on two B-splines) and the degenerate `villarceau_common` (2.0e-8
  in volume, 2.1e-8 in area, 8.1e-9 in centre: the two Villarceau circles
  crossing at the points of tangency, approximated by seven B-splines). The
  geometry agrees with the reference; the default measure does not reach the
  allowance.
* Curves as built: lines and circles where a section is a circle (a plane
  normal to the axis cuts the torus in two circles about it, a plane through
  the axis in two of the tube's circles); every other section of a plane
  and the torus, the spiric sections (loops about the tube, caps' loops,
  the Villarceau circles, the figure eight of a plane tangent to the inner
  equator), is a B-spline approximation, 115 edges in all (147 lines, 142
  circles, no conic other than circles). The results' edges' tolerances are
  1e-6 where they hold a B-spline and 1e-7 otherwise.
* The degenerate fixtures are valid results natively: `top_tangent_fuse`
  two solids touching along the top circle, `inner_tangent_common` one
  solid (the wall tangent to the inner equator, its section a figure eight
  as four B-splines), `villarceau_common` one solid (half the torus, by
  symmetry), `vertex_on_cut` one solid (the box's vertex 2e-17 of the case's
  size from the torus). The decisions refuse them (`Degenerate`).
* Nine solids' counts (five cases) change under
  `ShapeUpgrade_UnifySameDomain` (faces, edges, vertices as built, then
  unified): `strip_cut`'s two ends 7/15/10 to 6/14/10 and 6/12/8,
  `frame_cut`'s three caps 3/3/2 to 2/2/2 (the fourth, crossing OCCT's seam,
  5/8/5 both ways), `wedge_cap_common` 4/4/2 to 3/3/2,
  `top_tangent_fuse`'s torus 2/4/2 to 1/0/0 and box 7/13/9 to 6/12/8,
  `villarceau_common` 5/9/5 to 4/9/5 (the torus's faces split at OCCT's
  seams merged). The kernel's results will be compared with the unified
  counts, as S9c's.
* The reference itself is checked independently
  (`generate_torus_boolean_fixtures.py --check`): closed forms of 19 of the
  20 pairs (an aligned box by rectangles inside both circles of each normal
  slice; the torus less half-spaces by circular segments of both circles)
  within 3.3e-40 in exact frames and 1.5e-16 in turned ones (the declared
  degenerate pairs within 2.3e-29), inclusion and exclusion 1.1e-40, the
  common two ways 4.6e-41, the area identity and every face's classes
  3.0e-35 (1.3e-39 but for the vertex on the torus), a second slicing
  direction (the meridian half-planes) 1.6e-40, the wall two ways 2.1e-37,
  Monte Carlo 2.8 standard errors. Its scan for near coincidences found none
  outside the declared degenerate pairs, and every other pair keeps its
  vertices at least 0.18, its faces' planes at least 0.031 (from tangent
  planes whose points of contact lie on the face) and its edges' crossings
  at a sine of at least 0.54 from tangency with the torus (relative to the
  case's size).
