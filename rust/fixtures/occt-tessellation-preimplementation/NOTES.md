# Tessellation (T-a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks and TKMesh
(`build_pinned_occt.py --toolkit TKTopAlgo --toolkit TKMesh`).
`capture.json` records the Rust revision and that no kernel tessellation
code existed (`rust/kernel/src/tessellation.rs` absent). `inputs.txt` holds
the 28 bodies of `tessellation-cases.txt` as the kernel's existing `.brep`
writer produces them (`brep_io_probe prisms` and `bodies`), then the 56
mesh rows, a coarse (a hundredth of the case's scale, 0.5 rad) and a fine
setting (a thousandth, 0.3 rad) per body; `oracle.cpp` is
`occt_tessellation_oracle.cpp` as captured, and `native.txt` its `S` rows:
`BRepMesh_IncrementalMesh` with an absolute deflection, sequential, on a
clean copy of each body (the row's fields are listed in the probe's header).

## Observations

Checked against the independent reference (`tessellation_reference.py`) on
the welded native meshes (the probe's `dump` mode, not part of the capture):
its largest sampled distances equal the probe's own ElSLib measurements
within `2e-16` of each case's size on every row, which checks the reference
as much as OCCT.

* Every solid is watertight once nodes are joined through OCCT's edge
  polygons (no free, non-manifold or misoriented mesh edge), with the
  boundary's Euler characteristic; each face body's free edges are exactly
  its boundary loops. OCCT duplicates seam nodes (the `nodes` and `welded`
  columns differ on every periodic face).
* All twelve prisms, the box, the frustum, the inverted and tilted cones,
  the hemisphere, the whole and tilted tori, the torus wedge and both face
  bodies stay within the requested deflection at both settings.
* The requested deflection is exceeded on seven rows: `cone_apex` coarse
  (`0.045` against `0.03`), `sphere` fine (`0.0103` against `0.005`),
  `sphere_far` fine (`0.0041` against `0.002`), `sphere_zone` coarse
  (`0.094` against `0.05`) and fine (`0.0096` against `0.005`),
  `torus_segment` fine (`0.035` against `0.008`, with inward normals) and
  `torus_inner_half` fine (`0.035` against `0.008`). OCCT's own recorded
  `Poly_Triangulation::Deflection` states the excess on each of them.
* Apexes and poles give degenerate triangles (two nodes joined into one):
  one per apex or pole on `cone_apex`, `cone_inverted`, `hemisphere`, two
  on `sphere` and `sphere_far`.

`compare_tessellation.py` requires every later run to reproduce these rows
and the `.brep` texts, reviews OCCT's differences with fingerprints
(`occt-tessellation-divergences.json`), and requires the kernel's meshes to
pass every check of the reference.
