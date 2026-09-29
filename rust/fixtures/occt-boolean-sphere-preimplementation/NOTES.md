# Booleans of a sphere against polyhedral prisms (S9d.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28 by
`compare_sphere_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9d.1's fixtures). `capture.json` records the Rust revision and
that the kernel's S9d.1 code did not exist (`rust_sphere_boolean_exists`
false: `rust/kernel/src/solid/boolean/curved/sphere.rs` absent; S9b.2's
stored model refusing every Boolean with a sphere with
`OutOfDomain("a Boolean of a solid with curved faces or edges in any
position (S9c)")`, its spherical face not a plane, so the kernel's probe
reports every case `unsupported`). `inputs.txt` holds the 30 cases of
`boolean-sphere-cases.txt` as explicit constructions
(`identity_reference.native_boolean_case`: a sphere's block one `sphere`
row, its frame's origin, normal and x axis, radius and latitudes),
`oracle.cpp` is `occt_boolean_oracle.cpp` as captured and `native.txt` its
results; each sphere, cap or zone is `BRepPrimAPI_MakeSphere(gp_Ax2, R,
low, high)`, each prism `BRepPrimAPI_MakePrism` of its profile face in its
own frame. Frames: `XY` and `SIDE` (caps and zones, centred at the origin),
boxes in `TILT` and `LEAN`, a whole sphere in `TILTX`.

## Observations

* Every case is done, every result and solid valid; `octant_common` (the
  box's corner at the sphere's centre, three faces through it) reports
  warnings. The solid counts are the reference's in all 30
  (`sphere_boolean_reference.py`: the regularized Boolean's maximal
  connected regions, by convexity), the two pieces of `bar_through_cut`,
  `slab_cut` and `lshape_cut` and the 3 empty results included (the common
  of a dome standing on a box's face, and of a face and of an edge tangent
  to the sphere).
* The 27 nonempty totals agree with the reference within 3.3e-9 relative in
  volume and 2.5e-9 in area (`corner_box_common`: a box's corner, its walls
  cutting the sphere in circles off its axis), the centre within 4.5e-10 of
  the case's size (`edge_wedge_common`); in turned frames within 3.1e-9
  (`bar_tilt_common`). Where no wall cuts the sphere, every section a
  circle of latitude of its own frame (`face_cap`, `zone_post`, the slab,
  the cavities, the coplanar discs, `inscribed_cut`), within 4e-15;
  `octant_common`, its walls cutting along meridians, within 1e-14. The
  larger differences, 1e-12 to 3.3e-9, are BRepGProp's
  default integration of spherical faces bounded by circles off their
  parallels (the reference's closed forms agree within 9.3e-40 in exact
  frames); all are inside the 2e-8 allowance, so no review.
* The degenerate fixtures are valid results natively: `face_tangent_fuse`
  two solids touching at a point (the box keeps the point as a vertex on
  its face: 6 faces, 12 edges, 9 vertices), `inscribed_cut` one solid, the
  box's cavity touching the sphere at its eight vertices, and the commons
  of a tangent face and a tangent edge empty. The decisions refuse them
  (`Degenerate`).
* Coplanar discs: `dome_on_box_fuse` (a hemisphere's disc on a box's top
  face, opposite orientations) one solid of 7 faces, the disc gone and the
  top face holed at the dome's rim; `hemisphere_in_box` (the disc on the
  box's bottom face, the same orientation) cut to a cap and common to a
  zone of 2 and 3 faces.
* Three solids' counts change under `ShapeUpgrade_UnifySameDomain` (faces,
  edges, vertices as built, then unified): `zone_side_common` 7/15/10 to
  6/12/8, `hemisphere_box_common` 6/12/8 to 5/9/6 and `lshape_cut`'s
  larger piece 7/15/10 to 6/12/8 (faces of one surface split by the
  zone's or the profile's seams merged). The kernel's results will be
  compared with the unified counts, as S9c's.
* The reference itself is checked independently
  (`generate_sphere_boolean_fixtures.py --check`): closed forms of every
  pair (a rectangle's part of the sphere's sections integrated along the
  box's axis, both surfaces' parts inside the other) within 9.3e-40 in
  exact frames and 2.2e-16 in turned ones, inclusion and exclusion 9.2e-41,
  the common two ways 1.4e-42, the area identity 3.3e-40, every face's
  classes 1.8e-40, a second slicing direction 1.6e-40, Monte Carlo 2.7
  standard errors. Its scan for near coincidences found none.
