# Booleans of a cone or frustum against polyhedral prisms (S9d.3a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-29 by
`compare_cone_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9d.3a's fixtures). `capture.json` records the Rust revision and
that the kernel's S9d.3 code did not exist (`rust_cone_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/cone.rs` absent). The kernel's
probe reports every case `unsupported`: S9b.2's stored model refuses every
Boolean with a cone with `OutOfDomain("a Boolean of a solid with curved
faces or edges in any position (S9c)")` (all 30 cases, the message read
from `protocol::run` by a throwaway example), its wall not a plane.
`inputs.txt` holds the 30 cases of `boolean-cone-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`: a cone's block one
`cone` row, its frame's origin, normal and x axis, bottom and top radii and
height), `oracle.cpp` is `occt_boolean_oracle.cpp` as captured (its `cone`
row new with this capture; the prism and sphere rows unchanged) and
`native.txt` its results; each cone or frustum is
`BRepPrimAPI_MakeCone(gp_Ax2, R1, R2, H)`, each prism
`BRepPrimAPI_MakePrism` of its profile face in its own frame. Frames:
`XY` and `SIDE` (prisms along `x`), cones in `TILT`, `LEAN` and `TILTX`
against prisms in `XY`.

## Observations

* Every case is done, every result and solid valid; `apex_plane_common` (a
  box's wall through the cone's apex, the common a half cone) reports
  warnings. The solid counts are the reference's in all 30
  (`cone_boolean_reference.py`: the regularized Boolean's maximal connected
  regions, by convexity): the two pieces of `slab_cut`, `bar_cut` and
  `tilt_slab_cut`, the cavities of `cone_in_box_cut`, `box_in_frustum_cut`
  and `inside_tilt_cut` one solid each, and the 2 empty commons (a box on
  the frustum's top disc, a face tangent along a ruling).
* All 30 match within the 2e-8 allowance, so no review. The 28 nonempty
  totals agree with the reference within 1.4e-8 relative in volume and
  8.8e-9 in area, the centre within 1.4e-9 of the case's size
  (`ellipse_cut`: the frustum above an oblique face, its elliptic
  section); in turned frames within 2.9e-9 in volume and 5.5e-9 in area
  (`lean_apex_common`, an oblique cone on an ellipse). Where every section
  is a circle of the cone's own frame (`normal_cut`, a face normal to the
  axis; the cavities; `on_top_fuse`) within 4.1e-16, and a wall through the
  apex (`apex_plane_common`) within 3.3e-16. Every section
  is an exact conic: the results' edges are lines, circles, ellipses,
  parabolas and hyperbolas (no B-spline), their tolerances at most 9.9e-7
  (`parabola`) and 1.5e-7 otherwise. The same native results measured by
  BRepGProp with an adaptive accuracy of 1e-10 (a diagnostic build of the
  oracle, not committed) are within 5.0e-9 in volume, 2.3e-9 in area and
  4.0e-10 in centre (`apex_corner_common`, two hyperbolic walls near the
  apex), the ellipses within 2.9e-11: the default integration dominates
  the larger differences.
* The degenerate fixtures are valid results natively: `apex_plane_common` a
  half cone (4 faces, 8 edges, 4 vertices as built: the wall's section two
  rulings through the apex), `ruling_tangent_fuse` two solids touching
  along the ruling (the frustum keeps it as an edge: 4 faces, 6 edges, 4
  vertices), the common of the tangent face empty, `vertex_on_cut` one
  solid, the box's vertex on the wall. The decisions refuse them
  (`Degenerate`).
* Coplanar discs: `on_top_fuse` (a box on the frustum's top disc,
  opposite orientations) one solid of 8 faces, the disc gone and the box's
  bottom holed at its rim; `base_wall_common` (a box standing on the base
  plane, the same orientation, a hyperbolic wall) one solid of 5 faces.
* Thirteen solids' counts (12 cases) change under
  `ShapeUpgrade_UnifySameDomain` (faces, edges, vertices as built, then
  unified): `parabola_cut` 5/9/6 to 4/6/4, `hyperbola_cut` 9/18/12 to
  8/16/11, `hyperbola_common` 4/6/4 to 3/4/3, `apex_corner_common` 5/9/6 to
  4/7/5, `slab_cut`'s piece 5/9/6 to 4/6/4, `bar_cut`'s 7/15/10 to 6/12/8,
  `base_wall_common` 5/9/6 to 4/6/4, `tilt_ellipse_common` 4/6/4 to 3/4/3,
  `tilt_slab_cut`'s pieces 4/6/4 to 3/4/3 and 5/9/6 to 4/7/5,
  `apex_plane_common` 4/8/4 to 3/5/3 and `ruling_tangent_fuse`'s frustum
  4/6/4 to 3/3/2 (faces of the wall split at OCCT's seam, along the
  frame's x, merged). The kernel's results will be compared with the
  unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_cone_boolean_fixtures.py --check`): closed forms of every pair
  (an aligned box by rectangles inside the sections' discs, a half-space by
  circular segments, a slab as two half-spaces, along the axis) within
  1.2e-40 in exact frames and 1.2e-16 in turned ones, inclusion and
  exclusion 9.2e-41, the common two ways 2.2e-42, the area identity
  1.8e-40, every face's classes 9.2e-41, a second slicing direction (every
  section an ellipse) 1.8e-40, Monte Carlo 2.7 standard errors. Its scan
  for near coincidences found none.
