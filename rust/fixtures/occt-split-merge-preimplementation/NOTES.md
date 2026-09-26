# Height split and stacked fuse: preflight, before any Rust M3 code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`. Rust revision at
capture: `8807667c` (T1 accepted). The uncommitted work was only the new
Python reference (`split_merge_reference.py`), its generator, new operation
kind and role codes in `identity_reference.py`, new encoding vectors and this
probe; `rust/kernel` had no split or fuse code.

## OCCT source review

* `BRepAlgoAPI_Splitter` runs `BOPAlgo_Splitter`: the arguments are split by
  the tools and only argument pieces are kept. `BOPAlgo_Builder::PrepareHistory`
  records, for every source subshape, the splits kept in the result as
  `Modified`, section edges and vertices created by face/face and edge/face
  interferences as `Generated` (`LocGenerated`), and `Remove` for a shape that
  has no split in the result and is not itself in it. An untouched shape has
  no record and stays in the result as the same `TShape`.
* A solid split by a face is two solids that share the cut face (one
  `TShape` in both, opposite orientations) and its section edges and vertices.
  The cut face is the tool face's split; the tool is not an argument, so
  nothing in the arguments' history reaches it.
* `BRepAlgoAPI_Fuse` of two solids that touch along a common face removes the
  common faces; the coplanar walls stay separate faces and the shared cap
  edges stay as edges between them. `ShapeUpgrade_UnifySameDomain(shape,
  true, true, true)` then merges same-domain faces and edges and records its
  own `History()`: every rebuilt shape is `Modified` into its new `TShape`,
  merged faces and edges into the one merged shape, and removed ones `Remove`.
  `BRepTools_History::Merge` chains the fuse and unification histories.
* DRAW: `bapisplit` and `bapibop` (`BOPTest_APICommands.cxx`) and
  `unifysamedom` (`SWDRAW_ShapeUpgrade.cxx`) store their history with
  `BRepTest_Objects::SetHistory` when history is enabled. The original
  history cases that use them (`bugs/modalg_7/bug29333_*`, `bug21264`,
  `bug28113_*`) operate on faces, edges and lines from `plane`, `mkface`,
  `line` and `mkedge`, not on prisms.

## Probe

`occt_split_merge_oracle.cpp` reads rows from
`generate_split_merge_fixtures.native_rows`: named prisms (profile face at the
start offset, prism vector, rigid transforms, as for the history probe) and one
operation: `split` by a square planar face through the transformed axis point
at the split offset, `fuse` then unify, or `splitfuse`, the split followed by
the fuse and unification of its two solids. For every stage it prints, per
argument vertex, edge, face and solid (`TopExp::MapShapes` order), the counts
and signatures of `Modified` and `Generated` images and `IsRemoved`; every
result solid's signature, validity and distinct subshape counts; and every
result vertex, edge and face as `kept` (an argument subshape), `image` or
`none`. The `composed` stage queries the original prism in the merged
split-fuse-unify history.

## Observations (160 scenarios: 26 split-fuse-split, 6 fuse, 128 corpus)

* Every input prism, split piece and fused solid is valid; every split gives
  two solids, every fuse one. All cases took at most 0.043 seconds.
* Split: every wall, vertical (and seam) edge and the solid is `Modified`
  into exactly two pieces; each wall `Generates` one section edge and each
  vertical edge one section vertex. Caps, cap edges and cap vertices have no
  record and are kept. Each piece has the counts of a prism (8, 12, 6, 6, 1,
  1 for a square). The shared cut face is the only result shape no query
  reaches (`none`, 154 times, one per split).
* Fuse and unify: every argument subshape on the outer sides is `Modified`
  into one rebuilt shape with equal geometry; walls and vertical edges of both
  bodies into the merged wall or edge; the shared caps, their edges and
  vertices and both solids are `Remove`d.
* Composed: the original walls, vertical edges and seams are `Modified` into
  the merged shape and at the same time reported `IsRemoved` (2,030 edges,
  2,031 faces): `Merge` keeps the removal of the intermediate pieces. The
  original solid is `Remove`d.
* `split_labels_inserted_vertex` (a hexagon with a collinear inserted vertex):
  unification also merges the two coplanar walls on either side of that
  vertex and removes its vertical edge and both its vertices, which the
  profile keeps.

Expected contract differences, to review against the Rust histories:

* OCCT reports a split parent as `Modified` into its pieces; Rust reports
  `Split` (by design, `IDENTITY_AND_HISTORY.md` contract 3).
* OCCT shares one cut face, and its edges and vertices, between the two
  pieces; Rust bodies are separate, each with its own `Generated` cut face,
  edges and vertices. OCCT's cut face is the tool's image, reached by no
  argument query; Rust generates it from the walls.
* OCCT relates solids; Rust splits and merges the solid region and relates
  bodies by id.
* Unification's rebuilt `TShape`s are `Modified` with equal geometry; Rust
  keeps those entities `Unchanged`.
* Unification merges coplanar walls of the profile itself; Rust keeps every
  profile segment.
