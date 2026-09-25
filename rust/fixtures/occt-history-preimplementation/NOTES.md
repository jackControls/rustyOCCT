# Extrusion and transform history: preflight, before any Rust identity or history

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`. Rust revision at
capture: `427cecb3`. The uncommitted work was only the new Python references,
generator and this probe; `rust/kernel` had no identity, history or attribute
code.

## OCCT source review

* `BRepPrimAPI_MakePrism` delegates to `BRepSweep_Prism`. `Generated(S)`
  returns `myPrism.Shape(S)` for a profile subshape: the swept shape one
  dimension up (vertex to edge, edge to face, face to solid).
  `FirstShape(S)` and `LastShape(S)` return the copies of `S` at the start
  and end of the sweep; `FirstShape()` and `LastShape()` return the caps.
  `IsDeleted(S)` is false for every profile subshape. There is no `Modified`.
* `BRepTools_History(arguments, algorithm)` asks the algorithm for
  `Generated`, `Modified` and `IsDeleted` of every supported subshape
  (vertex, edge, face, solid) of the arguments, and stores lists keyed by
  shape. It does not include `FirstShape`/`LastShape`, and it never checks
  that every output is reached. Whether history is complete depends on each
  algorithm.
* `BRepBuilderAPI_Transform` with copy reports `Modified(S)` through
  `BRepBuilderAPI_ModifyShape`: one image per subshape of the transformed
  shape. DRAW's `ttranslate`/`trotate` change the location instead and
  record no history.
* DRAW `prism` (`BRepTest_SweepCommands.cxx`) calls
  `BRepTest_Objects::SetHistory({base}, prism)`. `savehistory`, `generated`,
  `modified` and `isdeleted` (`BRepTest_HistoryCommands.cxx`) query that
  saved `BRepTools_History`.

## Probe

`occt_history_oracle.cpp` reads rows from `identity_reference.native_case`:
the profile face at the start offset, the prism vector and each rigid
transform as a 3x4 matrix. It prints each query result as a geometric
signature (vertex point; edge curve type and points at first, middle and last
parameter; face surface type, area and centroid; solid volume and centroid).
It also classifies every output vertex, edge and face as reached directly by
a query, reached through a reached shape's subshapes, or not reached.

## Observations (98 cases: 34 explicit, 64 corpus)

* Every profile face built from the stored counter-clockwise outer wire and
  reversed hole wires is valid, and every prism and transformed copy is
  valid.
* Every profile vertex generates exactly one line edge; every edge generates
  exactly one face (plane or cylinder); the face generates the solid.
  `FirstShape`/`LastShape` give the start and end vertex, edge and cap.
* Every output vertex, edge and face is reached directly. For prisms, OCCT's
  history plus First/Last shapes is complete.
* Every subshape of every transform step has exactly one `Modified` image.
* All cases took at most 0.023 seconds.

Expected contract differences, to review against the Rust history:

* OCCT reports the solid as `Generated` from the face; Rust histories relate
  vertices, edges and faces only (bodies are referenced, not related).
* OCCT's start/end copies are `FirstShape`/`LastShape`, outside
  `Generated`; Rust reports them as `Generated` from the same profile element
  with roles `BottomVertex`/`TopVertex`/`BottomEdge`/`TopEdge`/caps.
* OCCT's caps come from the face; Rust caps are generated from every boundary
  label (the profile face has no label of its own).
