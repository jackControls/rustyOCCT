# Generic B-rep validator: preflight, before any Rust implementation

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871` (identical in this
checkout). Rust revision at preflight start: `12098ee3`.

## OCCT source review

`BRepCheck_Analyzer` runs `BRepCheck_{Vertex,Edge,Wire,Face,Shell,Solid}`:
`Minimum` checks each shape alone, and `InContext` checks each subshape inside
each ancestor. Statuses come from `BRepCheck_Status` (37 values).

* Edge/face agreement (`BRepCheck_Edge::InContext`): the edge must be flagged
  SameParameter and SameRange, otherwise the geometric check is skipped and the
  flags are reported. `BRepLib_ValidateEdge` then compares the 3D curve with
  `surface(pcurve)` at 23 uniform parameters by default (`myControlPointsNumber
  = 22`). The optional "exact" method uses `GeomLib_CheckCurveOnSurface`
  numerical maximization. Neither is a certified bound over the whole range.
* Vertices (`BRepCheck_Vertex`): the vertex point against each curve's
  endpoint parameter and each pcurve's endpoint on the surface, within the
  vertex tolerance.
* Wires (`BRepCheck_Wire`): Closed (3D connectivity), Closed2d, Orientation
  (each vertex joins exactly one incoming and one outgoing use), SelfIntersect
  (2D pcurve intersections, via `Geom2dInt`).
* Faces (`BRepCheck_Face`): IntersectWires, ClassifyWires (imbrication by 2D
  classification), OrientationOfWires, NoSurface, RedundantWire.
* Shells (`BRepCheck_Shell`): Closed (every edge used exactly twice by the
  shell, seams counted in their face), Orientation (the two uses of each edge
  must have opposite senses), NotConnected, RedundantFace, and
  InvalidMultiConnexity.
* Solids (`BRepCheck_Solid::Minimum`): a face in two shells gives
  InvalidImbricationOfShells. More than one non-hole closed shell gives
  EnclosedRegion. A shell outside another gives SubshapeNotInShape.
  `BRepCheck_ToolSolid::IsHole` classifies an infinite point with
  `BRepClass3d`, and `IsOut` classifies a point of one shell against another.

## Rust contract for milestone 1

A `Topology` is built only through a public validating constructor. Its report
is a complete list of typed issues with the offending entity, not a first-error
string. Exact combinatorial decisions:

* references in range; every vertex, edge and face used; each face in exactly
  one shell; no empty loop, face or shell
* each loop closes in 3D vertex order
* within each shell, every edge has exactly two uses with opposite senses. One
  use is free; more than two is non-manifold. A seam is two uses in one face.
  An edge used by two shells is an error.
* each shell is face-connected. Each vertex link (the fan of face corners
  around a vertex) is a single cycle, so pinch vertices are rejected.
  Euler-Poincaré per shell, `V - E + 2F - L = 2 - 2g`, gives an even
  characteristic `<= 2`.

Certified geometric decisions, three-valued. A check passes only on a
certified upper bound `<= tol`, fails only on a certified lower bound
`> tol`, and otherwise reports an explicit `Uncertified` issue. Never sample.

* curve and surface definitions: finite values, positive radius, nonzero
  line length, `0 < |sweep| <= 2π`, and valid frames
* edge curve endpoints within tolerance of their vertices
* `sup_t |C(t_edge) - S(P(t))|` over the WHOLE use, with the edge parameter
  reversed for a reversed use
* loop winding in UV: the outer loop is positive and inner loops negative with
  respect to the face orientation, and inner loops lie inside the outer loop
* shell orientation: the outer shell's signed volume is positive, each cavity's
  negative, and each cavity lies inside the outer shell and outside the other
  cavities

Out of scope for milestone 1: 2D self-intersection of loops, face/face
intersection (OCCT's BRepCheck does not do the latter either), tolerance
healing, and operation history.

## Native bridge encoding

Rust curves use a normalized fraction `t in [0,1]` for both 3D curves and
pcurves. OCCT needs SameParameter, so the generator computes OCCT curves whose
native parameter matches at corresponding points:

* a 3D line uses `Geom_Line` over arc length `[0,L]`
* a 3D arc uses `Geom_Circle` over the angle range, with its axis flipped for
  negative sweep
* a plane pcurve line uses a `Geom2d_Line` over the same arc length (the plane
  frame is orthonormal); a plane pcurve arc uses a `Geom2d_Circle` whose axis
  angle and sense make its parameter equal the 3D angle
* cylinder pcurves are `Geom2d_Line`s whose origin is shifted so that their
  parameter is the 3D angle (circles) or the axial length (seams)

A reversed use maps the pcurve by the reversed edge parameter. Seams use
`BRep_Builder::UpdateEdge(E, C_forward, C_reversed, F, tol)`. Pairs whose
speeds cannot match (e.g. a mutated pcurve radius) are still sent; OCCT's
range/flag response is itself an observation. The native probe reads explicit
OCCT constructions produced by the Python generator, so it contains no
geometry mapping of its own.

## Pre-implementation observations (capture.json, native.json)

Captured at Rust revision `12098ee3`. The uncommitted work was only the new
Python tools, fixtures and this probe; no Rust validator existed. The 54
generated cases are independent of Rust and OCCT: 19 valid solids and 35
mutations. The two reference-error cases cannot be built as OCCT shapes. On the
other 52, OCCT 8.1.0 (pinned headless SDK) completes every case in at most
0.02 s:

* All 19 valid solids are valid natively. They include holes, arcs (convex
  and concave), full circles with seams, one and two cavities, rotated and
  far-translated copies, and a millimetre-scale box at tolerance 1e-9.
* Every mutation the oracle rejects is also rejected natively. The status
  classes usually correspond: free edge / NotClosed, same sense /
  BadOrientationOfSubshape, vertex off curve / InvalidPointOnCurve, pcurve /
  InvalidCurveOnSurface, UV seam gap / wire NotClosed, non-manifold /
  InvalidMultiConnexity, loop outside / InvalidImbricationOfWires, cavity not
  inverted / EnclosedRegion, cavity outside or nested / SubshapeNotInShape.
* Contract differences to review later:
  * Unused vertices and edges are outside the analyzed shape, so OCCT
    reports the solid valid.
  * An inverted outer shell gives SubshapeNotInShape.
  * A 1e-3 UV gap left by a shifted pcurve is not flagged by OCCT's Closed2d.
  * Faces listed twice give NotConnected plus InvalidImbricationOfShells.
  * A zero-length edge or zero-radius cylinder gives several geometric
    statuses.
* Three mutated pcurves cannot be encoded with matching parameter speed. They
  are marked `inexact_pcurve_encodings`, so an encoding approximation is never
  read as an OCCT verdict.
