# Cone history (S3): native MakeRevol observations

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
(`build_pinned_occt.py`). `inputs.txt` is `primitive-cases.txt`, the 22
cones of the S3 capture; `oracle.cpp` is `occt_revolve_oracle.cpp` as
captured; `native.txt` its output.

## When it was taken (REVIEW_NOTES.md R11)

The S3 cone capture (`../occt-primitive-preimplementation`) observed
`BRepPrimAPI_MakeCone` before any kernel cone code, but no history:
`MakeCone` reports none and DRAW `pcone` saves none. This capture of
`BRepPrimAPI_MakeRevol` was taken afterwards, at revision `1a76d29e` with
the kernel's cone builder and its derivations written but uncommitted
(`capture.json` lists the worktree). It was taken before the Python
enumeration of cone entities (`identity_reference.cone_entities`) and before
the role correspondences of `compare_revolve_history.py` were written, and
the builder was not changed after it. R11 records the deviation from R7 for
the user's decision.

## OCCT source review

* `BRepPrimAPI_MakeRevol` delegates to `BRepSweep_Revol` and
  `BRepSweep_Rotation`. `Generated(S)` returns the swept shape of a meridian
  subshape when `IsUsed(S)`; `IsDeleted(S)` is `!IsUsed(S)`
  (`BRepPrimAPI_MakeRevol.cxx`). `BRepSweep_NumLinearRegularSweep::IsUsed`
  reports a generated shape used when it became a subshape of another; a
  full, closed revolution leaves the solid and the discs swept by the radial
  edges unused, so their generators are reported deleted though both are in
  the result.
* A point on the axis is invariant: it generates nothing, or a degenerated
  edge when it closes a face (the apex). An edge on the axis generates
  nothing (`HasShape`).
* The swept face's surface comes from
  `GeomAdaptor_SurfaceOfRevolution::GetType`: a revolved line is a cone only
  when the cosine of its semi-angle is at most `1 - Precision::Confusion()`;
  otherwise the face is a `Geom_SurfaceOfRevolution` of the same line.
* `FirstShape(S)`/`LastShape(S)` return the meridian copy at the start and
  end of the sweep: for a closed revolution the seam and its vertices, the
  apex vertex, or shapes the result leaves out (the radial edges and the
  meridian face).

## Observations

Every revolved solid is BRepCheck-valid with the counts `MakeCone` gives
(2/3/2/2/1/1 with an apex, 2/3/3/3/1/1 for a frustum). The rim points
generate the circles (or, at an apex, a degenerated edge), the slant edge
generates the lateral face; the axis, the radial edges and the meridian face
are reported deleted and the discs and solid are reached by no query.
`nearly_cylinder` (1 - cos a = 2.98e-8) gets a surface of revolution: the
one reviewed difference (`../occt-revolve-history-divergences.json`).
