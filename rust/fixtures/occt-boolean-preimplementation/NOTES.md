# Booleans of two prisms (S9a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28. `capture.json`
records the Rust revision and that no kernel code for the Booleans existed
(`rust/kernel/src/solid/boolean.rs` absent; `Solid` has no `fuse`, `cut` or
`common`). `inputs.txt` holds the 45 cases of `boolean-cases.txt` as explicit
constructions (`identity_reference.native_boolean_case`: the object's
`native_case` rows, a `boolean fuse|cut|common` row, the tool's rows);
`oracle.cpp` is `occt_boolean_oracle.cpp` as captured, and `native.txt` its
results.

Each prism is `BRepPrimAPI_MakePrism` of its profile face (the frame's plane
at the start offset, the boundaries' wires, holes reversed), as the split
oracle builds them; the tool's frame is its own (its origin offset from the
object's, the axes bitwise equal). The object and the tool go to
`BRepAlgoAPI_Fuse`, `BRepAlgoAPI_Cut` or `BRepAlgoAPI_Common` with default
options. Per case: the number of solids, the result's validity
(`BRepCheck_Analyzer`) and whether the operation reported warnings; per
solid, ordered by volume and centre: volume, area and centre of mass
(`BRepGProp`), its face, edge and vertex counts as built and after
`ShapeUpgrade_UnifySameDomain` (coplanar faces and collinear edges merged),
and its validity.

## Observations

* Every case is done, every result and solid valid, no warnings. The solid
  counts are the reference's in all 45 (`boolean_reference.py`: the
  regularized Boolean's maximal connected regions), and the totals agree
  with it within 6.1e-15 relative in volume (`lens_common`), 4.9e-15 in
  area and 4.5e-16 of the case's size in the centre: all 45 match within the
  2e-8 allowance, no review needed.
* The empty results (`identical_cut`, `disjoint_common`, `edge_touch_common`,
  `tangent_out_common`, `cap_touch_common`: a common or cut leaving at most a
  face, an edge or a point) are empty compounds: 0 solids, as the reference
  reports them.
* Touching: prisms sharing a wall fuse into one solid (`edge_touch_fuse`, a
  20 x 10 box after unifying; `edge_partial_fuse`); a cut by a tool sharing
  a wall keeps the object (`edge_touch_cut_tilted`, 6 faces); stacked prisms
  with caps at equal heights fuse into one solid (`cap_touch_fuse`,
  `stacked_fuse_tilted`: one box after unifying).
* Where the recorded decisions differ from OCCT: a fuse of prisms touching
  along a vertical edge (`vertex_touch_fuse`: squares meeting at a corner;
  `tangent_out_fuse`: circles touching outside) is two valid solids, and the
  cut leaving a hole tangent to the outer circle (`tangent_in_cut`) one valid
  solid (4 faces, 5 edges, 2 vertices) whose boundary touches itself along
  the tangent line; the decisions make all three `Degenerate` until the
  kernel holds non-manifold bodies (the fixtures' `expect degenerate`).
* Counts as built are not the kernel's: OCCT keeps each input's cap and
  wall images apart (`rects_fuse`: 14 faces, 32 edges, 20 vertices; 10, 24
  and 16 after unifying) where the decisions merge coplanar caps and
  coincident walls (`Merged`). After unifying, OCCT also merges collinear
  edges of the two inputs into one (`collinear_fuse`: 8 faces, 18 edges, 12
  vertices, the bottom wall one face), which the decisions do not state for
  the traced profile. OCCT's cylinders keep a seam edge and seam vertices
  (`disc_inside_common`: 3 faces, 3 edges, 2 vertices, where the kernel's
  seamless circle prism has 3 faces, 2 edges, 0 vertices), and a wall
  crossed by a seam is two faces as built (`lens_common`: 5 faces, 9 edges,
  6 vertices; 4, 6 and 4 after unifying). The kernel's counts will be
  compared with the unified ones and reviewed.

`compare_boolean.py` requires every later run to reproduce these rows (on
another platform, its reviewed record; the Linux record is pending CI).
Until `rust/kernel/examples/boolean_probe.rs` exists the comparison lists
every case under `rust_unsupported`.
