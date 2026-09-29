# Sheets and wires split by a plane (S8e): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-28. `capture.json`
records the Rust revision and that no kernel code for the split existed
(`rust/kernel/src/body/split.rs` absent; `Body` has no `split_by_plane`).
`inputs.txt` holds the 25 cases of `split-sheet-cases.txt` as explicit
constructions (`identity_reference.native_case`: the frame's plane, the
boundaries' wires, as for the prisms of S8a and S8b, and a `make face` or
`make wire` row in place of the prism vector) with the splitting plane;
`oracle.cpp` is `occt_split_oracle.cpp` as captured, and `native.txt` its
pieces.

A sheet is `BRepBuilderAPI_MakeFace` of the wires on the frame's plane, a
wire the one `BRepBuilderAPI_MakeWire`; `BRepAlgoAPI_Splitter` splits either
by a planar face 2,000 wide. A sheet's pieces are the result's faces: side by
the centre of mass, area (`SurfaceProperties`), perimeter
(`LinearProperties` over the face's edges), centre, face, edge and vertex
counts, validity. The splitter keeps a split wire one wire with its edges
split (its images); the probe groups them: the input's edges in stored order,
each replaced by its images in order along it, each image on the side of its
midpoint (an image on the plane, within 1e-9, on the side of the image
before it), each maximal cyclic run on one side a piece: its length
(`LinearProperties`), centre, 0 faces, its edges and vertices, and whether a
wire made of them in order is valid. The images cover the result's edges in
every case.

The cases: sheets of a square (planes at an angle to the sheet, through two
corners, touching a corner, parallel off it and containing it), a U cut
along its notch's bottom edge (two pieces above, one below), a stadium
touched by a plane tangent to its arc, a disc crossed and missed, a square
with a circular hole touched by the plane, a square with a slot hole in the
tilted frame, and S8b's bulge (crossed twice), wave (through its knot's
point) and square with a lens hole (tilted frame); wires of a square (tilted
frame, through two corners), an L cut along an edge between runs on both
sides, a circle and a stadium crossed, a square parallel to the plane and
one lying in it (tilted frame, within binary64 of parallel: their traces lie
at `v = 1.1e16` and along the bottom edge), and S8b's blob (through two
joins), wave (four crossings: two runs on each side), capsule (tilted frame)
and bulge (tangent).

## Observations

* Every piece is valid and every case done; the sides and, for wires, the
  number of runs on each side agree with the reference (`split_reference.
  planar_rows`, `Planar.wire_runs`).
* 23 cases agree with the reference within 2e-8, most to 1e-15; the wires
  with spline edges to 1.1e-8 (`LinearProperties`' default integration).
* `sheet_wave_knot`: BRepGProp's default integration of the face bounded by
  two spans of the wave errs by 1.3e-3 in area (41.8202 against 41.875);
  with an adaptive accuracy of 1e-13 it gives 41.874999999999993 and the
  reference's centre, and Green's theorem in exact Fractions over the
  piece's boundary gives 335/8 and centre (15887/5025, 32561/10050), the
  reference's.
* `sheet_lens_tilted`: `LinearProperties` measures the lens's whole upper
  cubic 6.3908444 instead of 6.3909761 (`GCPnts_AbscissaPoint::Length` at
  1e-13, and the reference's quadrature), so the perimeters are off by
  3.5e-6 (above) and 5.7e-8 (below); the edges' abscissa lengths sum to the
  reference's perimeters within 1e-15. Areas and centres agree to 1e-15.
* Both are reviewed native differences in `occt-split-divergences.json`;
  the adaptive, abscissa and Green figures came from a scratch variant of
  the oracle and a scratch computation, not from the capture.
* Where the kernel's recorded decisions will differ in counts: OCCT splits
  an edge where the plane touches it (the stadium's arc, the bulge's spline:
  one edge and vertex more) and keeps a circle's seam vertex (the disc and
  circle pieces holding angle 0); the sheet with the tangent hole is split,
  its lower face holding the hole touching its chord at (5, 7) (7 edges, 6
  vertices, valid), where S8e's decisions make the pinch `Degenerate`. A
  sheet lying in the plane is reported on the side of its centre (above);
  a whole wire along the plane below, by the probe's rule.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record: `platform-linux/`, from CI run 36556520655). Until
`Body::split_by_plane` exists the probe cannot build a body case: the
comparison lists every case under `rust_unsupported` and `rust_probe_failed`.
