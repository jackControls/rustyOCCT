# Tessellation

A parallel track of `REVIEW_NOTES.md`: deflection-controlled, watertight
triangle meshes of bodies, for display and for downstream consumers that
need triangles. T-a, described here, covers faces on planes, cylinders,
cones, spheres and tori bounded by lines, circles and arcs: profiles and
prisms (S5), the primitives (S3), sheets and wires (S6) and imported bodies
of those. T-b adds spline edges, pcurves and faces; T-c procedural
intersection edges (D13 of `TOPOLOGY_MODEL.md`) once faces carry them. The
decisions were recorded before any code (`REVIEW_NOTES.md`, parallel
tracks).

## Contract

`tessellation::tessellate(&Topology, Parameters)` (also `Solid::tessellate`
and `Body::tessellate`) returns a `Mesh`: `nodes`, `triangles` grouped by
face (`faces`, one range per face in face order), each triangle's certified
`Bound`, one polyline per edge (`edges`, with each segment's `Bound`), the
node of each vertex, and the largest bounds.

* **Parameters.** An absolute linear deflection `δ > 0` and an angle
  `0 < θ <= π/2`, as `IMeshTools_Parameters::Deflection` and `Angle` of
  OCCT's BRepMesh; anything else is `NonFinite` or `OutOfDomain`.
* **Watertight by construction.** Every edge is discretized once, uniformly
  in its fraction; every face using it takes exactly that polyline's nodes
  and segments. A vertex is one node at its stored position. There is no
  welding by distance, no T-junction and no seam: a periodic face is meshed
  in a chart of its whole domain, so nothing is duplicated along a seam.
* **Orientation.** Every triangle's normal (right-hand rule) leaves the
  solid region behind its face; a face between void regions follows its
  oriented normal. A closed shell gives a closed, consistently oriented
  2-manifold mesh (every mesh edge in exactly two triangles, once each way)
  whose Euler characteristic is the shell's; an open sheet a manifold with
  its boundary polylines; a wire body polylines only; an acorn a node. The
  kernel checks these on every result before returning it
  (`ComputationLimit("tessellation contract")` otherwise, never observed).
* **Deflection.** Each triangle comes with an explicit map onto its face's
  surface and each segment with the map onto its edge's arc between the same
  fractions; the reported bound certifies that no point moves farther under
  the map. So every mesh point is within `δ` of its face's surface (as
  BRepMesh measures deflection) and every polyline point within `δ` of its
  edge. The parameter triangles cover the face's domain up to the slivers
  between its pcurves and their chords (none along a line pcurve; on a plane
  no wider than the edge's own deflection), where the map lands just outside
  the face, so every mesh point is within `2δ` of the face itself.
* **Angle.** The surface normal turns by at most `θ` over each triangle's
  parameter triangle, and an edge's tangent by at most `θ` over each segment
  (a `1e-12` relative allowance covers the stored frames' departure from
  orthonormality, a few units in the last place).
* **Determinism.** The same topology and parameters give the same mesh bit
  for bit: fixed orders, no hashing, no threads, a fixed-seed walk.
* **Limits.** A node farther than `δ/2` from its face at its pcurve point
  (only an imported or invalid body can have one) is `InvalidTopology`;
  boundary chords that cross between line edges, or still cross after
  twelve doublings of a face's curved edges, `InvalidTopology` and
  `ComputationLimit`; more than 4,000,000 nodes `ComputationLimit`; a loop
  wound twice or both ways, a vertex loop that is not a pole, and spline
  geometry `OutOfDomain` (T-b). Tessellation does not validate the body.

## Certified bounds

For a parameter triangle `p_i` with extents `U`, `V` and bounds `a >=
|S_uu|`, `b >= |S_uv|`, `c >= |S_vv|` over its box, the linear map
`Σ λ_i X_i ↦ S(Σ λ_i p_i)` moves no point more than `(a U² + 2 b U V + c V²)
/ 8` plus the nodes' gaps `|X_i - S(p_i)|` (`MATHEMATICS.md`). Two sharper
maps are certified beside it and the smallest bound is kept:

* **Tangential correction.** On every analytic surface `S_uv` is parallel to
  `S_u`, so it has no normal component. The normal part of the linear map's
  deviation is at most `(a U² + 2 t b U V + c V²) / 8`, `t` the triangle's
  normal turn; the tangential part is taken up by moving the parameter
  point, leaving a second-order remainder. This removes the cone's and the
  torus's mixed term, which dominates elongated triangles.
* **Fan map at a pole.** A triangle with one vertex at a sphere's pole or a
  cone's apex maps each ray from the pole to the meridian through the
  interpolated `u` of the opposite edge; it deviates by at most the opposite
  edge's interpolation plus a meridian chord, `c V² / 8`. On a cone the
  meridian is a ruling and the second term vanishes: an apex fan meets the
  bound of its base segment.

The coefficients are closed forms times a bound on the stored frame's
spectral norm (plane 0; cylinder `r`, 0, 0; cone `|R + v sin α|`, `|sin α|`,
0; sphere `R |cos v|`, `R |sin v|`, `R`; torus `R + r cos v`, `r |sin v|`,
`r`), and everything is evaluated in the outward-rounded binary64 tier of
`certified.rs`, its trigonometry included, as are the node positions
(`bounds.rs`). An arc segment deviates by at most `r φ² / 8` for its sweep
`φ` plus its end nodes' gaps, a line segment by its end gaps.

## Algorithm

* **Edges first.** An edge's segment count is the least meeting its curve's
  `δ` and `θ` (with a 0.9 share) and, for each curved face using it, the
  thin-triangle condition at `0.45 δ` and `0.45 θ`: a triangle on a boundary
  segment cannot have smaller extents than the segment, so the segment must
  leave the face room. Ring edges and closed edges get at least three.
* **Charts.** Each face in a planar chart of its domain on the cover: a
  plane in its frame coordinates; a periodic face whose loops do not wind in
  the sinusoidal chart `((u - u_c) |S_u|(v), s(v))`, `s` the arc length
  along `v`, which collapses a pole line to a point (a sphere lune); a face
  wound in `u` in an annulus chart `P(v) (cos u, sin u)`, `P` exponential in
  `s(v)` with its logarithmic range capped at 4, or the arc length from the
  pole when a vertex loop closes the band (a hemisphere, an apex cone); a
  torus face wound in `v` the same with `u` and `v` exchanged. Every chart
  is a homeomorphism of the face's domain onto a planar region; its
  orientation sign relates chart and parameter orientation.
* **Triangulation.** A constrained Delaunay triangulation of the chart
  polygons (`cdt.rs`): exact orientation predicates for every combinatorial
  decision, Shewchuk's filtered in-circle test that flips only on a certain
  violation, Sloan's flips to recover constraints, the domain by the parity
  of constraints crossed. A boundary chord crossing another, a node on
  another's segment, a repeated constraint or a loop on the wrong side of
  its chords (nesting changed by coarse chords, checked against the loops'
  traversal) doubles the face's curved edges and starts again.
* **Refinement.** On curved faces a triangle is split while its certified
  deflection or turn exceeds the request, its lifted parameter triangle does
  not have the chart's orientation, or its flat normal disagrees with the
  surface's: at the midpoint of the edge contributing most to its bounds
  (an edge at a pole has no `u` extent), unless that edge is a boundary
  segment, else at its centroid; both lie inside the domain. Planar faces
  need no interior node. Steiner points are mapped back through the chart's
  inverse.
* **Whole sphere and torus.** A face without loops is a structured grid in
  `(u, v)` (a sphere's poles single nodes with fans), shrunk until every
  triangle's certified bounds hold.

## Evidence

* **Independent reference.** `tessellation_reference.py` derives each
  fixture body's boundary from the identity case alone: a prism's caps and
  walls (each profile piece times the height range), a solid of
  revolution's meridian pieces in the point's own half-plane, a torus
  wedge's tube and end discs, with exact point distances to each patch and
  to its whole surface; exact area and volume (mpmath, primitive_reference
  for the revolved solids), the Euler characteristic and a face body's
  loops. It checks a mesh: closedness and opposite orientations (a face
  body's boundary loops), the Euler characteristic, every node on the
  boundary, 12 barycentric samples of every triangle within the request
  and within the triangle's reported bound of its face's whole surface and
  within `2δ` of the boundary, edge polylines within `δ`, normals leaving
  the solid (a centroid offset by twice the larger of `δ` and the
  triangle's own deviation lands outside), and the enclosed volume positive
  and within `δ (A + A_mesh)` of the exact one; `test_tessellation_reference.py`
  requires it to name every deliberate failure of hand-made box and
  octahedron meshes (a flipped or missing triangle, a moved node, a
  deflection or bound understated, inward normals). `generate_tessellation_fixtures.py`
  writes 28 bodies (twelve prisms with polygons, circles, arcs, notches,
  lenses, scallops, holes, a tilted and a far frame; a box; four cones; four
  spheres; five tori, among them an inner half, a wedge and a tilted one;
  two face bodies) with a coarse (`scale / 100`, 0.5 rad) and a fine
  (`scale / 1000`, 0.3 rad) setting (`tessellation-cases.txt`,
  `-expected.tsv`).
* **Native capture before code.** `occt_tessellation_oracle.cpp` reads each
  body as the kernel's existing `.brep` writer produces it and runs
  `BRepMesh_IncrementalMesh` at each setting: counts, OCCT's own
  deflections, the deflection measured against each face's surface,
  watertightness after joining nodes through OCCT's edge polygons,
  orientation, area and volume. The capture
  (`fixtures/occt-tessellation-preimplementation`) was committed before any
  kernel tessellation code; the probe's own measurements agree with the
  reference's within `2e-16` of each case's size.
* **Kernel tests.** `tessellation.rs` checks all 56 meshes against the
  reference's expectations with its own closed-form surface distances
  (samples within each triangle's bound), closedness, orientation, the
  Euler characteristic, the volume bound, determinism and convergence; the
  parameters' errors; wire bodies' polylines and a typed budget error; the
  spline cases of `brep-cases.txt` out of domain; and every one of the 54
  certified `data/occ` solids without spline geometry (planes, cylinders,
  cones, spheres and tori with general loops, merged seams and holes),
  whose volume must lie within the deflection times the areas of its
  certified mass enclosure.
* **Native bridge.** `compare_tessellation.py` reproduces the capture and
  the `.brep` texts, requires the kernel's meshes (`tessellation_probe`) to
  pass every check of the reference on all 56 (they do) and compares OCCT's
  welded meshes with the reference: 42 matches and 14 reviewed differences
  (`occt-tessellation-divergences.json`). OCCT exceeds the requested
  deflection on seven rows (an apex cone coarse, the whole spheres fine, the
  sphere zone at both settings, the torus v-segment and inner half fine; up
  to 4.4 times the request, which its own `Poly_Triangulation::Deflection`
  records) and leaves one or two degenerate triangles at every apex and
  pole, where its degenerated edge's polygon joins into one node. Every
  OCCT solid is watertight after the join.
* **Fuzzing.** The `tessellation` target builds polygon and arc prisms,
  arc and polygon face bodies, cones, spheres and zones, tori, v-segments
  and wedges, rigidly moved, at a deflection from an eighth to a 2048th of
  the body's size and an angle from 0.2 rad to `π/2`, and requires a mesh
  meeting the whole contract (`FUZZING.md`). A clean local 600-second
  campaign at `c9dbf4f6` (AddressSanitizer, standard 20-second/2 GiB
  limits) ran 4,855 mutation executions after replaying 350 inputs in
  209 s (14,000 coverage edges, 961 MB peak) without an artifact; a
  120-second smoke campaign before it found nothing either.

Triangle counts against OCCT on the fixtures (kernel/OCCT): prisms 1.0–1.2
(the kernel's arcs meet the thin-triangle condition), the apex cones
0.2–0.3 (fans certified by their fan map), frusta and tilted cones 1.4–2.6,
spheres and zones 1.1–2.5, tori 0.9–2.5, face bodies 0.6–0.8 (a planar
face is its boundary's triangulation only). Where the ratio exceeds one the
kernel pays for a bound that holds everywhere, not only at sampled points;
OCCT's excess deflections above are the other side of that trade.

## Acceptance

T-a's acceptance needs kernel CI green at the accepted revision (the
tessellation bridge on macOS and Linux), the schedule run's full fuzz
replay green there, and a clean local 600-second `tessellation` campaign.
Recorded here once met.
