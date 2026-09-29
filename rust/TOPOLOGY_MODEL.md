# Topology model: decisions, migration and evidence

This is the decision record for the kernel's topology model and the delivery
guide for migrating to it. The decisions were taken on 2026-09-25 after
comparing OCCT, Parasolid and CGM; T1 and T2 are implemented (see
[Implementation notes](#implementation-notes-t1) and
[T2](#implementation-notes-t2)). It sits
beside `IDENTITY_AND_HISTORY.md`, whose contracts are model-agnostic and whose
milestones M0–M2 are accepted on the current prism model at `34efd36c`. The
migration here, **T1**, must land before M3, because M3's first splits and
merges would otherwise be built on seams that the decided model removes.

**Status:** decided; T1 accepted at `e4adb869` (see
[Acceptance](#acceptance)); T2 accepted at `0913b5e2`; the surfaces of
revolution (S3) and spline cells (S4) accepted by `ac92ec89`.

## Why this model

OCCT's shape is a reference into a tree of shared sub-shapes: adjacency is
derived, loops are unordered, space is not partitioned, periodic faces need
seam edges, identity is a pointer and tolerances only grow. Parasolid's body
is a cellular partition of space with stored adjacency, no seams, identifiers
and tracking. CGM's body is a cell complex of matter with a mandatory journal,
versioned operators and one fixed resolution with measured gaps. The kernel
takes Parasolid's skeleton, CGM's contracts, and keeps three things of its own:
stored certified pcurves, one resolution per body with computed enclosures,
and exact arithmetic. The long-term target is a Parasolid rival; the model
must therefore be at least as expressive as Parasolid's, and every structural
claim must be checkable by our own validator, because no other kernel's tests
exercise our structure.

## Decisions

Each decision states what it forbids, so a reviewer can reject code that
violates it without re-deriving the reasoning.

* **D1. Cellular skeleton.** A body is a set of regions, shells, faces, loops,
  fins, edges and vertices. A shell is one connected boundary component of one
  region and lists oriented face sides, wire edges and acorn vertices. A face
  has two sides, and each side belongs to exactly one shell. A loop is an
  ordered cycle of fins, or a single vertex. A fin is one loop's oriented use
  of an edge. *Forbids:* faces with only one side, adjacency reconstructed
  from a tree, unordered loops.
* **D2. Regions from day one, including the infinite void.** Space is fully
  partitioned. Region 0 of every body is the infinite void; it is structure
  and has no id. Bounded regions are solid or void and are entities with
  identity. A cavity is a bounded void region, not "a later shell". A prism
  has one solid region. *Forbids:* an "outer shell by convention", cavity
  classification by ray casting at validation time.
* **D3. Stored adjacency.** Each edge stores its fins in radial order. Each
  face stores its two shells. Each shell stores its region. The validator
  certifies stored adjacency against geometry; algorithms may rely on it.
  *Forbids:* scanning all faces to find an edge's neighbours.
* **D4. No seams.** Periodic faces have no seam edge. A loop on a periodic
  surface closes modulo the period and carries an integer winding number per
  periodic direction; pcurves live in the universal cover and may cross the
  period. Closed curves without vertices are ring edges. A closed surface may
  be a face with zero loops. A pole that an edge ends on is a vertex loop.
  *Forbids:* seam edges, degenerate edges, seam vertices.
* **D5. Stored, certified pcurves per fin.** Every fin carries its pcurve and
  the validator certifies agreement with the edge curve, as today. Planar
  faces may declare a derived pcurve, computed exactly from the edge curve and
  the plane, instead of storing one. *Forbids:* on-demand projection onto
  non-planar surfaces, a pcurve without a certified agreement check.
* **D6. Sense on entities, matter from regions.** A face has a sense against
  its surface, an edge against its curve, a fin against its edge. Which side
  of a face is material follows from the kinds of the regions its two shells
  bound; there is no separate matter flag to keep consistent. *Forbids:*
  orientation on references composed down a tree, internal or external
  orientation states.
* **D7. One resolution per body, computed enclosures.** Contract 5 of the
  identity guide. Imported tolerant entities become enclosures with
  provenance `Imported`. *Forbids:* a global session precision, growing a
  tolerance to make a step succeed.
* **D8. No instancing inside bodies.** A placed copy is a new body with its
  own ids related by `Modified`. Assemblies are a later layer. *Forbids:*
  location chains on entities.
* **D9. Body type is computed.** Validation classifies a body as solid,
  sheet, wire, acorn or general from its regions and domains. Operations
  state their required class as a precondition and check it. *Forbids:* a
  stored type that algorithms branch on.
* **D10. Immutable bodies, sharing deferred.** Operations never modify an
  input body. Sharing untouched entities between input and output bodies is an
  optimisation for later; the arena must not make it impossible, which value
  ids and per-body slots already guarantee. *Forbids:* in-place mutation.
* **D11. C1 within a cell.** Within an edge or a face interior, geometry is
  C1 in its parameters; a tangent discontinuity is a vertex or an edge. The
  exact spline modules already detect interior-knot discontinuities; R4 of
  `REVIEW_NOTES.md` implements the check (`edge_not_c1`, `pcurve_not_c1`,
  `face_not_c1`). *Forbids:* a single edge across a corner.
* **D12. Blend surfaces are an open decision with a constraint.** Procedural
  blend surfaces evaluated through their supports are preferred, as in
  Parasolid and CGM, but undecided. The `Surface` representation must remain
  extensible with a certified evaluation contract. *Forbids:* approximating a
  fillet to a B-spline as the only representation.
* **D13. Intersection curves are procedural** (U7 of `REVIEW_NOTES.md`,
  answered 2026-09-27). An intersection edge whose curve is neither a conic
  nor a spline stores its two surfaces and a certified parameterisation, and
  evaluates by a certified iteration with an enclosure, as CGM's edge
  curves and Parasolid's SP-curves do. It lies on both faces by definition,
  so its pcurves are exact projections. It is approximated by a spline only
  for tessellation and interchange, with the approximation's bound recorded
  beside it. Conics and splines stay explicit where the intersection is one.
  *Forbids:* a spline approximation as the kernel's own representation of
  an intersection curve.

Two contracts from CGM are adopted into the identity guide rather than here:
every operation carries an algorithm level and replays at the recorded level
(requirement H8), and every operation's history is complete and checked.

## The four models

| Aspect | OCCT | Parasolid | CGM | This kernel |
| --- | --- | --- | --- | --- |
| Space model | none, outer shell by classification | full partition, solid and void regions | matter only, volumes in lumps | full partition, region 0 the infinite void |
| Non-manifold, mixed dimension | internal and external orientations, compsolid, compound | general bodies | cell complex, any configuration | general bodies, type computed |
| Adjacency | derived | stored, fins in radial order | stored bounding relations | stored, fins in radial order, certified |
| Loops | unordered wires | ordered, vertex loops, zero-loop faces | ordered, outer, inner, full | ordered, vertex loops, zero-loop faces, winding numbers |
| Periodic faces | seam edges | no seams | seams forbidden | no seams, pcurves cross the period |
| Edge geometry | curve, pcurves, polygons, flags | curve, optional SP-curves | edge curve with mapped representations and max gap | exact curve, certified pcurve per fin, enclosure per fin |
| Tolerance | per entity, only grows | fixed precision, tolerant entities | fixed resolution, measured gaps | one resolution per body, computed enclosures |
| Identity | pointer plus location | identifiers | tags plus generic naming | structural value ids |
| History | optional, unverified | tracking records | mandatory journal, checked | mandatory history, all entities, checked |
| Replay | none | not public | software configuration levels | algorithm level per operation |
| Immutability | mutable, shared | mutable, rollback marks | frozen, shared cells | immutable, sharing deferred |
| Instancing | in-model locations | assembly instances | assembly level | none in a body |
| Native file | `.brep`, no ids | XT | proprietary | own format, XT as interchange |
| Exact decisions | none | not public | high precision, not exact | exact predicates, certified enclosures |

## Data model

The sketch refines the current arena in `topology.rs`; names that already
exist keep their meaning. Sketches are not normative; the decisions are.

```rust
index_type!(RegionId); index_type!(ShellId); index_type!(LoopId); index_type!(FinId);

pub enum RegionKind { Solid, Void }

/// regions[0] is the infinite void: kind Void, no id, no derivation.
pub struct Region { pub kind: RegionKind, pub shells: Vec<ShellId> }

pub struct Shell {
    pub region: RegionId,
    pub sides: Vec<(FaceId, Side)>,   // Side::Front | Side::Back
    pub wire_edges: Vec<EdgeId>,      // edges with no fins
    pub acorn_vertices: Vec<VertexId>,
}

pub struct Face {
    pub surface: Surface,
    pub sense: Orientation,           // against the surface normal
    pub loops: Vec<LoopId>,           // may be empty for a closed surface
    pub front: ShellId,
    pub back: ShellId,
}

pub enum Loop {
    Edges { fins: Vec<FinId>, winding: [i32; 2] }, // per periodic direction
    Vertex(VertexId),                               // a pole or an immersed vertex
}

pub struct Fin {
    pub edge: EdgeId,
    pub sense: Orientation,           // against the edge curve
    pub pcurve: Pcurve,               // Stored(Curve2) | DerivedOnPlane
    pub enclosure: Enclosure,         // M5
}

pub struct Edge {
    pub curve: Curve3,
    pub start: Option<VertexId>,      // None, None: a ring edge
    pub end: Option<VertexId>,
    pub fins: Vec<FinId>,             // radial order about the curve tangent
}
```

`Coedge` becomes `Fin`, `Face.loops: Vec<Vec<Coedge>>` becomes loop ids,
`Topology.shells: Vec<Vec<FaceId>>` becomes `Shell`s owned by regions, and
`Identity` gains a `Slot::Region`. The prism builder produces one solid region
whose single shell lists every face's front side, and the infinite void whose
single shell lists every back side.

### Entities with identity

Bodies, bounded regions, faces, edges and vertices have ids and appear in
histories. Shells, loops, fins and the infinite void are structure: they are
recomputed by validation, never named, never journaled. `EntityKind` gains
`Region` (dimension 3) with an appended encoding code, so every existing
identity vector stays valid. Extrusion adds `Role::Region` for its solid
region, generated from every boundary label.

## Validator

`Topology::check` keeps its exact combinatorial and certified geometric tiers
and its complete-issue-list contract (`BREP_VALIDATION.md`). The cell-complex
invariants below replace the shell-set invariants; kinds that no longer apply
(`same_sense_uses` for seams, `cavity_outside` and `nested_cavity` by ray
parity) are retired or re-expressed through regions.

Combinatorial, exact:

* every face side is listed by exactly one shell, and that shell's region is
  the one the face records for that side: `side_without_shell`,
  `side_in_two_shells`, `side_region_mismatch`
* every shell belongs to one region and is face-connected; region 0 exists and
  is the only infinite region: `region_without_shell`, `no_infinite_region`
* around each edge, the fins in radial order alternate the regions of the
  faces' sides consistently, and a manifold edge has exactly two fins:
  `radial_order_inconsistent`; free and non-manifold edges are classified,
  not rejected, unless the operation requires a solid
* each loop of edges closes in vertex order, or modulo the period on a
  periodic surface with the recorded winding numbers: `open_loop`,
  `winding_mismatch`
* a ring edge has no vertices and a closed curve; a vertex loop lies on the
  surface: `ring_edge_with_vertex`, `vertex_loop_off_surface`
* no cell is bounded twice by the same lower cell: `double_bounding`
* Euler–Poincaré per shell and per body class; the body class itself is a
  result, not an input

Geometric, certified with the two arithmetic tiers:

* pcurve against edge curve per fin, closure modulo the period, loop winding
  and inner-loop containment, all as today but on the universal cover
* C1 at interior knots of every edge and face (`continuity`)
* enclosures within the resolution (M5)

## Failure modes this model introduces

These cannot occur in OCCT's model, so OCCT's tests cannot find them. Each
has a checker kind above, a fuzz mutation and, where possible, a metamorphic
law that needs no oracle.

| Failure mode | Mutation that must be detected | Law |
| --- | --- | --- |
| Period wrap error in a pcurve or loop | shift one fin's pcurve by a period; flip a winding number | loop closure modulo period; face area equals the unwrapped integral |
| Radial order or region assignment error | swap two fins around an edge; point a side at the wrong shell | fins alternate regions; volume of solid regions equals the signed-shell integral |
| Vertex loop or zero-loop face mishandled | add a vertex loop off the surface; remove the only loop of an open face | area and centroid unchanged by a valid vertex loop |
| Id derivation drift | rebuild after a no-op refactor | ids of every fixture byte-identical |
| Incomplete history | drop, duplicate or retype a relation | checker clean; `then` associative |
| Unsound enclosure | store a bound below the truth | recomputation never exceeds the stored bound |

## Effects on accepted work

T1 changes structure, not geometry or ids, and the following are
requirements, not hopes:

* **Polygon prisms keep every id.** Every row of `identity-expected.tsv` and
  every relation of `history-expected.tsv` for polygon cases must be
  byte-identical after T1.
* **Circle prisms lose only seam entities.** Roles `Seam` and `SeamVertex`
  retire. The circular edges keep roles `BottomEdge` and `TopEdge` with
  ordinal 0 and become ring edges with no vertices. The wall face keeps its id.
  The fixture generator drops the seam rows by rule; nothing is hand-edited.
* **Native captures are untouched.** OCCT's encoding keeps its seams. The
  history bridge classifies native seam edges and seam vertices as
  structure-only by rule (an edge closed on its face, and vertices used only
  by such edges), requires every other native entity to map one to one, and
  keeps its 98 matches with zero failures. The rule is verified on every case
  by count synthesis: the synthesized OCCT counts of the Rust body must equal
  the native counts.
* **Mass properties, bounds and classification are bitwise unchanged** for
  every prism fixture.
* **The three original and two derived DRAW cases keep passing.** A third
  derived case exercises `pcylinder` counts through the synthesizer.

## OCCT structure interop (T2)

OCCT is one encoding of the neutral fixture model, and the upstream suite is
evidence for geometry, properties and validity verdicts, never for our
structure. T2 makes that usable in both directions:

* **Converter** from `.brep` (the format in `dox/specification/brep_format.md`)
  into the cell model: merge the two uses of a seam into one periodic loop
  with winding numbers, drop degenerate edges into poles, turn per-entity
  tolerances into `Imported` enclosures, compounds into bodies with regions,
  and classify every unsupported construct as `Unsupported` by name. Never
  approximate.
* **Writer** back to `.brep`, inserting seams and degenerate edges by rule,
  so a round trip through native `checkshape`, `nbshapes` and `props` must
  reproduce the original verdict, counts and properties.
* **Count synthesizer** for `checknbshapes`: from a body, the counts OCCT
  would report: one seam per periodic direction of every face whose loops wind
  that direction, one degenerate edge per pole loop, their vertices, loops as
  wires. Implemented once in the kernel and used by the DRAW adapter's
  `nbshapes`, so original assertions run unchanged.
* **Native selector** for index-based picks such as `s_5` after `explode`:
  run the construction in DRAWEXE, take that sub-shape's geometry, select the
  Rust entity by exact match. The kernel never reproduces OCCT's exploration
  order.
* **Coverage ledger** in the bridge report: every upstream assertion carries
  one of three statuses, `model-independent`, `mapped-and-verified` or
  `lost`, and the aggregate is a gate in `PRODUCTION_READINESS.md`. A
  mapping counts as verified only on a case where native DRAW confirms it.

Baseline from the 2026-09-25 survey of `tests/`, to be recomputed by the
ledger tool: 17,879 cases; 11,599 load external data; 7,658 touch the viewer;
property assertions in 8,415 files, sub-shape counts in 4,689, validity in
2,655, tolerance maxima in 860.

## Milestones

Each milestone ships the full pattern of `IDENTITY_AND_HISTORY.md`: kernel
code, an independent Python oracle generator run with `--check`, in-crate
fixture tests, a fuzz target registered in `run_fuzz.py` and the fuzz
workflow matrix, a source-pinned native bridge where OCCT can observe the
behaviour, and doc updates. Native observations are captured before
implementation.

### T1 — migrate the prism model to the cell complex

* Replace the arena per the data model: regions, shells with sides, loops,
  fins, ring edges, vertex loops, winding numbers. Rebuild the prism builder
  seamlessly: circle boundaries yield ring edges; polygon boundaries are
  unchanged. `Topology::from_parts` accepts the new parts; `Solid::extrude`
  and `transformed` generate the solid region.
* Extend `Topology::check` with the invariants above and retire the seam
  kinds. Extend `brep_reference.py` with the same invariants on a seamless
  model, keeping its OCCT-row output seamed.
* **Fixtures:** regenerate `brep-expected.tsv`, `identity-expected.tsv` and
  `history-expected.tsv`; the polygon subsets must be byte-identical. Add
  seamless circle, torus-like and vertex-loop cases to the neutral model as
  the builders allow.
* **Fuzz:** extend `brep_validation` with the six mutations of the failure
  table and the `identity` and `history` targets with circle profiles in both
  directions. Clean local 600-second campaigns for all three.
* **Native bridges:** `compare_brep.py` and `compare_history.py` gain the
  structure-only classification and the count-synthesis verification; match
  counts must not fall and failures must stay zero.
* **Upstream tests:** `nbshapes` in the DRAW adapter reports synthesized
  counts; add a derived `pcylinder` case checked on both backends.
* **Docs:** `BREP_VALIDATION.md` invariants, `VALIDATION.md`, `FUZZING.md`,
  `SOURCE_MAP.md`, the acceptance record here and the milestone rows in
  `PORTING.md` and `PRODUCTION_READINESS.md`.

Gate: kernel CI, fuzz CI, native bridges on macOS and Linux, and a clean
600-second campaign at one revision, plus the five requirements under
[Effects on accepted work](#effects-on-accepted-work).

### T2 — OCCT structure interop and the coverage ledger

* Converter, writer, count synthesizer and native selector as specified
  above, in a kernel module for the format and bridge code for the selector
  and ledger.
* **Fixtures:** round trips of every prism fixture and of the 37 `.brep`
  files under `data/occ`: import, `Topology::check`, export, native
  `checkshape`, `nbshapes` and `props` equal to the original; unsupported
  constructs reported by name, with the count of each recorded.
* **Fuzz target `brep_io`:** structure-aware mutations of `.brep` text;
  malformed input is rejected with a typed error, never a panic; valid input
  round-trips to identical synthesized counts and properties.
* **Native bridge:** the round trip itself, on both platforms.
* **Upstream tests:** the ledger statuses in `run_upstream_tests.py` and
  `upstream-draw.json`; the first data-dependent original cases once the
  public test data directory is configured.
* **Docs:** `UPSTREAM_TESTS.md` ledger section, `PRODUCTION_READINESS.md`
  gate, acceptance record here.

T2 touches tools, the format module and the bridge; it can run in parallel
with M4 by a second agent after M3 has landed.

## Program order

M0, M1 and M2 are accepted at `34efd36c`. The remaining order is fixed by
dependencies:

| Step | Depends on | Unlocks |
| --- | --- | --- |
| T1 cell model | M2 | M3 on the decided model, seamless circle histories |
| M3 height split and stacked fuse | T1 | first `Split`, `Merged`, `Deleted` on real bodies, region relations |
| T2 OCCT interop and ledger | T1, M3 for split fixtures | data-dependent upstream cases, measurable retained evidence |
| M4 attributes on operations | M3 | application colours, names, PMI links through Booleans |
| M5 enclosures as data | T1 | Contract 5 complete, `enclosure_*` kinds |

## Implementation notes (T1)

Where the implementation refines the sketches above, the refinement and its
reason:

* **Arenas.** `TopologyParts` holds vertices, edges, fins, loops, faces,
  shells and regions; faces list loop ids, loops list fin ids, edges list
  their fins. Pcurves are stored (`DerivedOnPlane` and fin enclosures wait for
  M5; M5 added enclosures on vertices, fins and faces, while `DerivedOnPlane`
  is still not implemented). A loop's winding is `[u, v]`; only `u` on a cylinder may be nonzero.
* **Prism layout.** Shell 0 lists every face's front side and belongs to the
  solid region 1; shell 1 lists every back side and belongs to the infinite
  void. A cavity adds its material shell to region 1 and a bounded void
  region whose shell lists the opposite sides. A shell listing exactly the
  opposite sides of an earlier shell is its twin; shell-level checks (face
  connectivity, Euler, orientation, containment) run once per surface.
* **Region orientation** is a flux sign per shell: `-∮ v f(u) du` per face in
  closed form, negated for back sides, positive for a bounded region's first
  shell and negative for the others. Cavity containment keeps the ray parity
  of the earlier validator, re-expressed per region, and gains a `+v`
  cover-crossing test on cylinders that counts every period alias of a wound
  loop.
* **Retired kinds.** The seam cases of `same_sense_uses` became `seam_edge`,
  since this model forbids seams; `cavity_outside` and `nested_cavity` stay,
  re-expressed through regions.
* **Structure-only vertices.** The rule "vertices used only by such edges"
  would keep OCCT's seam vertex, which its circles also use. The bridges
  classify a vertex as structure-only when every edge using it is a seam or
  an edge closed on that vertex, and at least one is a seam. Count synthesis
  verifies the rule on every case.
* **Count synthesis in T1.** `Topology::occt_counts` (one seam per wound face
  and periodic direction, one seam vertex per ring edge used by a wound loop,
  a wire per unwound loop plus one per wound face, shells and solids from
  solid regions) landed with T1, since the bridges verify with it and the
  derived `pcylinder` case needs it. Degenerate edges for poles join it with
  the surfaces that have poles (T2).
* **Fixture layout.** Region rows of `identity-expected.tsv` and region
  relations of `history-expected.tsv` are listed beside the corpus digests
  rather than inside them, so the vertex, edge and face digests of polygon
  corpora are byte-identical across T1 and prove the first requirement.
* **Algorithm levels (H8).** `AlgorithmLevel::FIRST` is the only level. The
  `_with` constructors run at `AlgorithmLevel::CURRENT`; `Solid::extrude_at`
  and `Solid::transform_at` replay a recorded level and reject any other with
  `Error::UnknownAlgorithmLevel`. `History::then` first gave a composite the
  level of its last step; since `REVIEW_NOTES.md` S1 a composite records
  every step's operation, kind and level (H8 in `IDENTITY_AND_HISTORY.md`).
* **Bitwise properties per host.** `prism-properties-baseline.tsv` was
  recorded on macOS/aarch64 at `26fc457f`, before any T1 code. Frames and
  rotations use the platform's trigonometry, so Linux and Windows differ from
  it in the last bits for rotated copies. There, CI checks out `26fc457f`,
  writes that host's rows with the pre-migration code and requires the
  migrated kernel to reproduce them bitwise; the requirement is unchanged
  properties on the same host, not equal libm across hosts.
* **Radial order law.** "Swap two fins around an edge" cannot be detected on a
  manifold edge: a cyclic order of two is unchanged, and the fuzz target
  requires exactly that reversal to stay valid. The detectable form, fins
  exchanged between two edges, is the mutation.

## Implementation notes (T2)

* **Order of evidence.** Unlike every earlier milestone, the native
  observations of the upstream corpus were captured after the reader,
  converter and writer existed, and the independent Python reader was
  written by the same author afterwards. The break is recorded in
  `fixtures/occt-brep-io-capture/NOTES.md`, with what it does and does not
  affect; it is not repaired.
* **Converter scope.** `occt_brep::import` converts a solid when every face
  is a plane or a cylinder (indirect cylinders flip the face sense), every
  edge a line or circle with both vertices, every orientation forward or
  reversed and every location rigid. A trimmed line or circle is its basis
  over a range, and a plane without a pcurve gets one from the edge
  (`CurveOnPlane`). Each solid becomes its own body: a solid region, the
  shell with the largest bounding box as its outer shell, and the others as
  cavities. Compounds are flattened. Free shapes, compsolids, degenerate
  edges, mesh-only faces and every other curve or surface are reported by
  name and counted. Degenerate edges belong with surfaces that have poles,
  which the kernel does not represent yet, so none is dropped into a pole.
* **Tolerances.** An imported body's resolution is its largest vertex,
  edge or face tolerance (`ImportedSolid::tolerance`). Since M5, each
  vertex's and fin's OCCT tolerance is its `Imported` enclosure (see
  `IDENTITY_AND_HISTORY.md`), which the checker verifies.
* **Writer.** Version 1 text. Each wound face gets one seam at a `u` where
  both of its wound loops have a vertex, or where a ring loop starts. A ring
  edge gets a seam vertex there and becomes a closed edge. A face whose wound
  loops share no such position is `Unwritable`: splitting an edge would
  change the counts the synthesizer promises. Text is not a fixed point of a
  round trip, because import numbers entities in discovery order and
  re-derives plane frames. Round trips compare counts, cells, bit-identical
  vertices and native properties instead.
* **Near-frequency bound.** OCCT prints `2π` as `6.28318530717959`. A ring
  pcurve over that span against a circle swept by `2π` gives two harmonic
  terms a few ulps apart in frequency, which the validator bounded apart
  (`2r`) and could not certify. Both validators now pair terms adjacent by
  frequency (`BREP_VALIDATION.md`), and the importer snaps such spans to
  `2π`.
* **Native selector.** Picks are matched by measure and centre of gravity,
  never by order, and only faces and edges so far. A seam pick is lost.
* **Ledger.** The statuses and their rules are in `UPSTREAM_TESTS.md`. No
  original assertion is mapped-and-verified yet; the adapter's `restore`
  waits for an original data-dependent case that a configured data
  directory makes runnable on both backends, because `data/occ` alone makes
  none runnable.
* **Fuzzing.** The first `brep_io` smoke run found reader references beyond
  their tables (edge representation and face locations). The reader now
  range-checks every reference; the input is a retained regression.

## Surfaces of revolution (S3)

S3 of `REVIEW_NOTES.md` adds the surfaces of revolution in the order cone,
sphere, torus, one per accepted revision (R6).

### Cone

* **Surface.** `Surface::Cone { frame, radius, half_angle }` is OCCT's
  `Geom_ConicalSurface`: `S(u, v) = O + (R + v sin a)(cos u x + sin u y) +
  v cos a n`, `v` arc length along a ruling, the apex at `v = -R / sin a`.
  A valid cone has a finite `R ≥ 0` and `0 < |a| < π/2`. Pcurves live on the
  universal cover in `(u, v)`, as on a cylinder.
* **Pole.** The apex is a vertex loop of the wall. On a cone face whose edge
  loops wind once in total, the first vertex loop is the pole: it closes the
  band, so the windings need not balance, and its vertex must lie within
  tolerance of the apex (`pole_off_apex`) as well as on the surface. The
  periodic area and the orientation flux take the pole into account
  (`MATHEMATICS.md`). The count synthesizer adds one degenerate edge per
  pole, as OCCT has.
* **Builder.** `Solid::cone_with`, `cone_at` and `cone_in` build
  `BRepPrimAPI_MakeCone(gp_Ax2, bottom, top, height)`: discs, then the wall
  whose loops are the rings (bottom `+u` at `v = 0`, top `-u` at the slant)
  and, for a zero radius, the pole. Every entity is generated from its
  meridian element (`IDENTITY_AND_HISTORY.md`); the history's kind is
  `Revolve`. Mass properties are the general certified ones (`U2`), and the
  cone's enclosures are measured like every builder's. Split and fuse
  rebuild prisms: on a cone they return `OutOfDomain`, and `Solid::profile`
  is `None`.
* **Interop.** The reader takes surface record 3 (with an indirect axis
  becoming the cone about `-N` with `v` and the angle negated) and accepts
  degenerated edges. In the seam merge `du` is scaled by `|R + v sin a|`,
  and a run that is exactly one degenerated edge on a cone face becomes the
  pole's vertex loop. A degenerated edge anywhere else is `DegeneratedEdge`,
  unsupported. The writer inserts OCCT's structure: the seam from the base
  ring's seam vertex to the apex, and the degenerated edge (no 3D curve, the
  pcurve `v = v_apex` over one turn from the seam's `u`, the apex its vertex
  at both ends).
* **DRAW.** `pcone name R1 R2 H`, and `vprops`/`sprops` of restored bodies and
  cone faces through the certified mass properties. The derived case
  `pcone_counts` checks an apex cone, a frustum and a base apex on both
  backends.
* **Evidence.** The native `MakeCone` capture came before any kernel cone
  code (`fixtures/occt-primitive-preimplementation`). The `MakeRevol` history
  capture came after the builder was written, before it was committed
  (`R11`, `fixtures/occt-revolve-history-capture/NOTES.md`). The independent
  references are `primitive_reference.py` (counts, mass, faces, edges,
  vertices), the cone rules of `brep_reference.py` and `cell_reference.py`,
  and `identity_reference.cone_entities`.

### Sphere

* **Surface.** `Surface::Sphere { frame, radius }` is OCCT's
  `Geom_SphericalSurface`: `S(u, v) = O + R (cos v (cos u x + sin u y) +
  sin v n)`, `v` the latitude, the poles at `v = ±π/2`. `S_u × S_v = R² cos v`
  times the outward radial direction.
* **The first zero-loop face.** A whole sphere is one face without edge
  loops, a closed surface; `empty_face` is reported only on other surfaces.
  Its count synthesis is OCCT's one wire: a seam from pole to pole, a
  degenerated edge at each pole and their two vertices.
* **Poles.** On a sphere face whose edge loops wind once in total, the first
  vertex loop is the pole on the band's material side: the material lies
  left of a loop's traversal, so a band winding `+u` on a forward face
  closes at the north pole (`v = π/2`), otherwise at the south pole. The
  pole's vertex must be within tolerance of that pole (`pole_off_apex`, as
  for a cone's apex). A loop may also pass through a pole: consecutive fins
  meet at the pole vertex with different `u`, a UV gap of no length there
  (`|R cos v| = 0`) closed by a chord along the pole's line; OCCT has a
  degenerated edge there, and the count synthesizer counts it.
* **Builder.** `Solid::sphere_with`, `sphere_at` and `sphere_in` build
  `BRepPrimAPI_MakeSphere(gp_Ax2, R, low, high)` for latitudes
  `-π/2 <= low < high <= π/2`, an end at `±π/2` (the binary64 value) being a
  pole: discs at the other ends, then the wall with its rings and, when
  exactly one end is a pole, that pole. Ids follow the meridian as for the
  cone, the arc being segment 1 and a pole's vertex having role `pole`.
  Split and fuse return `OutOfDomain`.
* **Mass properties.** On a sphere the integrands are trigonometric in both
  `u` and `v`; `-∮ F du` is exact along parallels (`F(u, v0)` collapses to a
  trigonometric polynomial in `u`), zero along meridians, and enclosed along
  any other line (a chord) by the interval of `F` over its bounding box. The
  whole sphere is bounded on the cover by the north pole's line
  (`MATHEMATICS.md`).
* **Interop.** The reader takes surface record 4; a sphere bounded only by
  its poles imports as the whole sphere, and a degenerated edge among other
  uses is dropped (the loop passes through the pole). The writer writes the
  whole sphere's seam as a meridian circle, its poles' vertices and
  degenerated edges, and a degenerated edge at every pole pass.
* **DRAW.** `psphere name R [angle1 angle2]` (degrees) and the derived case
  `psphere_counts`. The original case `bugs/modalg_6/bug27264_2` (a restored
  whole sphere) now passes on both backends.
* **Evidence.** Both native captures, `MakeSphere` and `MakeRevol`, came
  before any kernel sphere code (`fixtures/occt-sphere-preimplementation`,
  `fixtures/occt-sphere-revolve-capture`).

### Torus

* **Surface.** `Surface::Torus { frame, major, minor }` is OCCT's
  `Geom_ToroidalSurface`: `S(u, v) = O + (R + r cos v)(cos u x + sin u y) +
  r sin v n`, a ring torus (`R - r > tol`); a horn or spindle torus is
  `degenerate_surface`, and the reader reports it as
  `NonRingToroidalSurface`, unsupported. `S_u × S_v = r (R + r cos v)` times
  the tube's outward normal.
* **Windings in `v`.** A torus is the one surface periodic in both
  parameters: a loop winds `[wu, wv]` and closes shifted by `(2πwu, 2πwv)`.
  Windings in `v` are `winding_mismatch` on every other surface, and on a
  torus they must balance. A face wound in `v` (a wedge's wall, bounded by
  two meridian circles) has no outer/inner order: the periodic area `∮ u dv`
  over all loops has the face's sign, each unwound loop the opposite one. A
  face wound in both directions (a torus knot) is
  `uncertified_loop_winding`. UV gaps scale `du` by `R + r cos v` and `dv`
  by `r`.
* **Whole torus.** One face without loops, as the whole sphere. Its count
  synthesis is OCCT's one wire of two seams (the meridian and the latitude
  circle through one vertex).
* **Builder.** `Solid::torus_with`, `torus_at` and `torus_in` build
  `BRepPrimAPI_MakeTorus(gp_Ax2, R, r, low, high, angle)` for the whole
  torus (`high - low` and `angle` the binary64 `2π`), a v-segment (a full
  turn of latitudes `low..high`: rings at both ends bounding discs, the wall
  wound in `u`) or a wedge (the whole tube from `v = 0` revolved by
  `0 < angle < 2π`: the tube's circles at both ends bounding discs, the wall
  wound in `v`). A segment of a partial turn is `OutOfDomain`, as is a
  segment whose meridian region (between the arc and the axis) is not
  simple. A v-segment is the revolved region between its arc and the axis,
  so its wall is reversed when the arc bulges toward the axis: OCCT's
  `inner_half` has the wall forward and is inside out (a reviewed native
  difference). Ids follow the meridian as for the cone; a wedge's discs come
  from the boundary and its circles from the arc (segment 1's start and end
  copies, as a prism's cap edges). Split and fuse return `OutOfDomain`.
* **Mass properties.** The integrands are trigonometric in `u` and `v`.
  Loops wound in `u`, or none, integrate `-∮ F du` with `F` from `v = 0`;
  loops wound in `v` integrate `∮ G dv` with `G` from `u = 0`, the same
  routine with `u` and `v` exchanged. The whole torus is bounded on the
  cover by both periods (`MATHEMATICS.md`).
* **Interop.** The reader takes surface record 5 (an indirect axis negates
  `v`); the seam merge joins seams in `u` and in `v`, and a torus bounded
  only by its two seams imports as the whole torus. The writer writes a
  wound face's seam in each wound direction (a latitude circle for `u`, a
  meridian circle for `v`) and the whole torus's two seams through one
  vertex.
* **DRAW.** `ptorus name R r [angle1 angle2] [angle]` (degrees) and the
  derived case `ptorus_counts`.
* **Evidence.** Both native captures, `MakeTorus` and `MakeRevol`, came
  before any kernel torus code (`fixtures/occt-torus-preimplementation`,
  `fixtures/occt-torus-revolve-capture`).

## Arcs in profiles (S5)

S5 of `REVIEW_NOTES.md` lets a profile boundary mix lines and circular arcs,
so the Extrude job takes sketches with fillets and slots.

* **Input.** `Boundary::path(points, segments, tolerance)`: segment `i` from
  point `i` to `i + 1` is `Segment::Line` or `Segment::Arc { center, radius,
  ccw }`. Both ends of an arc lie within tolerance of its circle; its sweep
  is the turn between them in its direction. A path of lines only is the
  polygon, so polygon and circle ids are unchanged. The path is stored
  counter-clockwise by its area with arcs' bulges; a clockwise input is
  reversed as a polygon is, arcs' directions flipped, labels mapped.
* **Validity** (`decide/arcs.rs`, exact where the distance is a sum of input
  terms, conservative in rational intervals otherwise): chords and radii
  above the tolerance; non-adjacent pieces farther apart than it (endpoint
  distances, the interior pair through the foot of a centre or along the
  line of centres, crossings); adjacent pieces neither doubling back nor
  meeting again (the reflection of the shared point) nor bringing a far
  end within tolerance. Holes against paths use the same pieces.
* **Prism.** An arc segment gives circular-arc bottom and top edges about
  the profile normal, the vertical lines at its ends and a partial cylinder
  wall: one unwound loop of line pcurves on the cover, the rectangle
  `[a, a + sweep] x [0, h]`, sensed so its normal leaves the material
  (forward on a counter-clockwise outer arc). Caps take circular-arc
  pcurves. Roles, ids, layout and history are a line segment's. Profile
  moments add each arc's Green integrals in closed form; the certified mass
  enclosure contains them. Points classify by exact ray parity over the
  arcs' parts monotone in `y`. The height split and stacked fuse work
  unchanged.
* **Interop and DRAW.** The writer and reader already carried arc edges and
  partial cylinder walls; every arc prism round-trips and OCCT reads it
  valid with equal counts and properties. The adapter runs upstream's
  `profile` sketch command (lines, `C` arcs, `F`, `O`, `P`, `X`, `Y`, `L`,
  `T`, `R`, `D`, `I`, `W`) and `prism` of its faces; the derived case
  `profile_arcs` passes on both backends.
* **Evidence.** The native `MakePrism` capture of thirteen arc prisms came
  before any kernel arc code (`fixtures/occt-arc-prism-preimplementation`);
  `compare_history.py --family arc` matches all thirteen. The identity
  reference enumerates their ids (`identity-arc-cases.txt`,
  `identity-arc-expected.tsv`); the builder reproduces the neutral
  generator's stadium, notch and half disc.

## Sheets, wires and acorns (S6)

S6 of `REVIEW_NOTES.md` admits bodies without a solid region: faces and
shells, wires and single vertices. The cell model already had their slots
(shells list wire edges and acorn vertices); S6 validates, builds, measures,
imports and writes them.

* **Classes (D9).** `Topology::class()` computes `Solid` (faces only, each
  between a solid and a void region), `Sheet` (faces only, no solid region),
  `Wire` (wire edges only), `Acorn` (one vertex) or `General`; nothing stores
  it. The independent reference computes the same class
  (`cell_reference.body_class`) for every valid fixture
  (`brep-classes.tsv`).
* **Validation.** An open sheet's faces have both sides in one shell of the
  infinite void: such two-sided faces carry no orientation flux and their
  shell no Euler condition, but a sound shell still needs one one-sided face
  more than it has two-sided ones, so a solid's shell never degenerates to
  a sheet. A closed shell without a solid bounds a bounded void region and
  is checked as a solid's shell. A wire edge may have no fin
  (`wire_edge_with_fins` otherwise) and an acorn vertex no edge
  (`acorn_vertex_used`); their curves and vertex gaps are checked as edges'.
  Every face-level check is unchanged.
* **Counts.** OCCT closes every closed edge at a vertex, so every ring edge
  counts one vertex, whatever it bounds (a cylinder's seam, a disc, a circle
  wire). Without a solid, a sheet of one face is a free face (no shell), of
  several one shell per face-bearing shell, twins counted once; a wire of
  one edge is a free edge (no wire), of several one wire. An imported
  one-face shell or one-edge wire would count differently; none occurs in
  the corpus.
* **Measures.** `Topology::measure_enclosure()` is a sheet's area, a wire's
  length (zero for an acorn) and the centre, certified: faces through the
  mass module's face integrals, lines and arcs in closed form, a full circle
  with the certified `pi`. A wire with a spline edge is not measured.
* **Builders.** `Body::face_from_profile_with` makes one planar face on the
  profile's frame bounded by every boundary; `Body::wire_from_boundary_with`
  one closed wire of a boundary's edges in stored order. Their operation
  kinds are `MakeFace` and `MakeWire`, the roles `Face` (from every
  boundary), `Edge` (from its segment) and `Vertex` (from its point), so a
  face body's entities derive from the same profile elements as a prism's
  start cap, bottom edges and bottom vertices. Rigid motions rebuild in the
  moved frame and keep every id.
* **Interop.** The reader imports every shell, face, wire, edge and vertex
  reached from the root outside a solid (`Import::free`): a shell is closed
  when every edge with a curve is used exactly twice, else a sheet. On a
  spline face, as on a plane, the loop enclosing the largest UV area comes
  first (OCCT stores wires in any order; a hammer face had its hole first).
  The writer writes a body without a solid as its closed shells, its face or
  open shell, its edge or wire (each edge oriented to continue the last) or
  its vertex; a ring edge no seam meets is closed at the start of its
  record. A spline with more poles than `MAX_POLES` is reported as
  `BSplineControlDataLimit`, and the independent reader models the limit.
* **Evidence.** Native `BRepCheck`, counts and `BRepGProp` of eighteen
  sheet, shell, wire and acorn models (`fixtures/occt-sheet-preimplementation`)
  and of the corpus's 6,223 free shapes (`fixtures/occt-free-shape-capture`)
  were captured before any kernel code accepted them. `compare_brep.py
  --family sheet` gives 17 matches and 1 reviewed difference (the open
  box's shifted pcurve, the prism shifts' Closed2d rule) with 15 counts and
  15 measures verified; `compare_brep_io.py` imports 6,188 free shapes with
  OCCT's counts and certified measures containing OCCT's, reads all of them
  back natively valid with equal properties, and reads the eleven face and
  wire bodies of `identity-sheet-cases.txt` back. The hammer's 31 faces the
  validator rejects though BRepCheck accepts them (C0 spline surfaces,
  pcurves certified farther from their edges than OCCT's stored tolerance,
  uses the exact tiers leave undecided) are pinned in `occt_brep.rs`.

## Spline cells (S4)

S4 of `REVIEW_NOTES.md` adds B-spline edges, pcurves and faces: R4's
continuity check first, then the geometric checks, mass properties and
interop.

* **Cells.** `Curve3::BSpline` and `Curve2::BSpline` hold a `SplineSpan`:
  a `BSplineCurve3` or `BSplineCurve2`, a closed range inside its domain
  (the whole domain when the kernel builds it, a sub-range when a file
  trims it) and a `reversed` flag, the fraction mapped affinely onto the
  range, backwards when flagged. `Surface::BSpline` is a
  `BSplineSurface3` in its own `(u, v)`; its loops do not wind.
* **Validation.** C1 by exact knot removal (R4); deviations by exact
  rational composition or second-order Taylor enclosures; vertex gaps from
  exact ends; UV gaps on a spline surface in 3D; signed and periodic areas,
  containment parity and orientation fluxes with spline geometry
  (`BREP_VALIDATION.md`, `MATHEMATICS.md`). Periodic spline surfaces are
  not certified.
* **Mass properties.** Enclosed with every spline geometry but a periodic
  spline surface (S4d): exact volume and moments on nonrational spline
  surfaces; since F8 certified Gauss–Legendre quadrature within `1e-12` of
  each property's scale on the fixtures (first-order strips and Green
  integrals only as a fallback) for areas, for
  rational surfaces and along spline pcurves on curved surfaces.
* **Interop.** The reader takes curve record 7, 2D curve record 7 and
  surface record 9; a trimmed curve (record 8) whose basis is a line,
  circle or B-spline reads as its basis, the edge's range bounding the part
  used. A range within `1e-12` of a domain end is snapped to it, and a
  periodic curve's range exceeding one period by at most `1e-12` of it is
  shortened to one period. A periodic spline surface
  (`PeriodicBSplineSurface`) and a seam on a spline surface
  (`SeamOnBSplineSurface`) are unsupported: their windings would need the
  domain's period. The writer writes records 7 and 9 with the span's
  range; a pcurve against its spline's parameter is written mirrored
  (knots `a + b - k`, only when exact), and an edge against its curve
  (S8b: a spline profile segment given clockwise, a hole's spline in an
  outer boundary after a split) is written along it, used the other way in
  every wire (its pcurves, stored in their uses' directions, are unchanged).
* **Evidence.** The native observations of the spline models came before
  any kernel code certified spline geometry
  (`fixtures/occt-spline-preimplementation`), and `BRepGProp`'s properties
  before any kernel code integrated them (`fixtures/occt-spline-properties`);
  `compare_brep.py --family spline` requires both to reproduce and compares
  verdicts, status classes, counts, enclosures and mass properties (7
  matches, 3 reviewed differences: BRepCheck has no continuity status, and
  a shifted pcurve's gap is masked by its edge's status). The independent references are `spline_cell_reference.py`
  (spans, jets, exact and interval evaluation), the spline rules and
  `mass_properties` of `cell_reference.py`, and `brep_io_reference.py` for
  the records.

## Plane sections (S8a.2)

S8a.2 of `REVIEW_NOTES.md` splits prisms by planes oblique to their axis
(`SPLIT.md`), which adds the cells of a plane's section of a cylinder.

* **Cells.** `Curve3::EllipseArc` is `origin + major cos a x + minor sin a
  y`, `a = start + sweep f` (OCCT's `Geom_Ellipse`); a ring edge's sweep is
  a full turn. `Curve2::EllipseArc` is the same in a plane's `(u, v)` with
  its axes along `u` and `v` (the split's cut face shares its frame's axes
  with every section ellipse). `Curve2::Sinusoid` is `u = start + sweep f`,
  `v = a0 + a1 cos u + a2 sin u`: the graph of a plane's height over a
  cylinder's angle, on the cylinder's universal cover.
* **Validation.** The harmonic bound certifies an ellipse edge against
  either pcurve: both are an affine term plus a rotating one at the edge's
  own frequency, so the sweeps a split stores equal cancel exactly and
  their phases' roundings remain. Degeneracy (both semi-axes above the
  resolution), signed and periodic areas, containment parity (an ellipse
  is a circle scaled about its centre), boundary clearance, orientation
  fluxes and mass properties are closed forms in `cos` and `sin`
  (`MATHEMATICS.md`); a wire's ellipse edge has no certified length (an
  elliptic integral).
* **Tessellation.** An ellipse's segment deflects at most its frame's norm
  times its major radius times the angle step squared over eight, and turns
  at most `major / minor` times the step.
* **Interop.** The writer prints `Geom_Ellipse` (3D record 3, its major
  radius first) and `Geom2d_Ellipse` (2D record 3) on a plane whose axes
  are the ellipse's, the edge's parameter its angle; a sinusoid pcurve has
  no OCCT record and is `Unwritable`. The reader keeps skipping ellipses.

## Conic sections and projections (S8d.2)

S8d.2 of `REVIEW_NOTES.md` splits cones, frusta and sphere zones by planes
neither normal to their axis nor containing it (`SPLIT.md`): D13's
procedural pcurves with explicit conic edges.

* **Cells.** `Curve3::HyperbolaArc` is `origin + major cosh t x + minor sinh
  t y`, `t = start + sweep f` (OCCT's `Geom_Hyperbola`); `Curve3::ParabolaArc`
  is `origin + t^2 / (4 focal) x + t y` (`Geom_Parabola`). `Curve2::Projection`
  is the exact inverse of its face's surface map applied to its fin's edge
  (at fraction `f` the edge's point at `f`, or `1 - f` for a reversed use),
  continuous on the surface's cover: `u` (and a torus's `v`) is the
  representative nearest the linear interpolation of lifts recorded every
  thirty-second of the edge, consecutive lifts within a quarter period.
* **Validation.** A projection of the fin's own edge onto the face's own
  stored surface deviates by zero by definition; its ends, gaps and
  degeneracy come from certified evaluation, and its twice-areas, periodic
  areas, orientation fluxes and mass moments from certified Taylor
  quadrature (`MATHEMATICS.md`). `+u` ray crossings of a projection are not
  decided: a split's projection loops are outer loops, or holes whose own
  first point is tested against the other loops.
* **Tessellation.** A hyperbola's segment deflects at most `σ max(major,
  minor) cosh(T)` times the step squared over eight (`T` the largest
  `|t|`) and turns at most `major / minor` times the step; a parabola's at
  most `σ / (2 focal)` times the step squared over eight, turning `1 / (2
  focal)` times the step. A projection's `(u, v)` range is the union of its
  certified enclosures on 64 pieces.
* **Torus sections (S8d.3).** `Curve3::Section` is a plane's section of a
  torus as a graph over one of its angles (`v = psi +- acos(-C / W)` over
  `u`, `u = atan2(b, a) +- acos(q / |(a, b)|)` over `v`), each edge's range
  clear of the turning points; a whole turn is a ring. Its jets come from
  the formula; on its own torus its projection is its own angles, lifted.
  Its segments' bounds come from interval jets of order 2 on 64 pieces
  (the largest `|C''|` and `|C''| / |C'|`). A face on a torus whose loops
  all run as holes is the torus less them.
* **Cylinders' meetings (S9c.2).** `Curve3::Meet` is two cylinders'
  meeting as a graph over the first's angle (S7b's ruled parameterisation):
  the ruling `frame.point(radius (cos u, sin u), v)` at `u = start + sweep
  f` meets the other cylinder (`|(w . x2, w . y2)| = other_radius`) where
  `a v^2 + 2 b v + c = 0`, and `v = (-b + sign sqrt(b^2 - a c)) / a`; an
  edge's range keeps the discriminant positive (no turning point), a whole
  turn is a ring. Its jets come from the formula; on its own cylinder its
  projection is its own angle and height, lifted. Its segments' bounds are
  a section's (interval jets on 64 pieces); a rigid motion moves both
  frames. The writer refuses it (`Unwritable`) until D13's interchange
  approximation. A sphere may be its other quadric (`other_sphere`,
  S9d.2).
* **A cylinder's and a sphere's meeting over the height (S9d.2b).**
  `Curve3::Rise`: at `w = start + sweep f` the cylinder's circle meets the
  sphere where `alpha cos u + beta sin u = g(w)` (`g` quadratic), `u = phi +
  sign acos(g(w) / rho)`, `rho` and `phi` the binary64 `hypot` and `atan2` of
  `(alpha, beta)`; an edge's range keeps `|g / rho| < 1` strictly. Jets,
  quadrature, tessellation rates, rigid motion and the writer's refusal as
  `Curve3::Meet`'s.
* **A cone's sections in Boolean results (S9d.3a).** `HyperbolaArc` and
  `ParabolaArc` edges bound Boolean results too; a rigid motion moves their
  frames, and a body's bounds hold their coordinates' extremes inside their
  arcs and every vertex (an apex's vertex loop among them).
* **Interop.** The writer prints `Geom_Hyperbola` and `Geom_Parabola` (3D
  records 5 and 6), and refuses a torus section; a projection pcurve has no record until D13's
  interchange approximation (a B-spline whose certified bound becomes the
  edge's tolerance), so every S8d.2 piece is `Unwritable` for now. The
  reader keeps skipping the conics.

## Acceptance

Record each milestone's clean revision, kernel CI job count, fuzz CI target
count, native bridge match, review and failure counts on macOS and Linux, and
the clean local 600-second campaign, in the style of `BREP_VALIDATION.md`.

* **T1 — accepted at `e4adb869`** (the migration is `77401b00`; `e4adb869`
  changes only the property-baseline test, the workflow and docs).
  * Rust kernel workflow: all eleven jobs passed, including the B-rep and
    history comparisons against OCCT 8.1.0 built from the pinned source, the
    original-test bridge, Rust 1.85 and WebAssembly. The first run at
    `77401b00` failed only the bitwise property test on Linux, Windows and
    Rust 1.85: the baseline was a macOS observation (see the notes above).
  * Fuzzing workflow: all twenty-one targets passed. On Linux,
    `brep_validation` replayed 1,417 inputs in 229 seconds, then ran 60.09
    seconds of mutation (813 executions, 7,406 coverage edges, 796 MB RSS
    peak); `identity` replayed 800 in 121 seconds, then 60.08 seconds (699
    executions, 7,210 edges, 675 MB); `history` replayed 880 in 148 seconds,
    then 60.09 seconds (740 executions, 7,976 edges, 712 MB). None produced
    an artifact.
  * Native B-rep bridge: 44 matches, 8 reviewed differences, 0 failures on
    macOS and Linux, unchanged from M1; one native seam status set aside as
    structure-only; synthesized counts equal native on all 21 cases valid on
    both sides.
  * Native history bridge: 98 matches, 0 reviewed differences, 0 failures on
    macOS and Linux; 38 circle cases classify OCCT's seam and its two
    vertices as structure-only, and every case's synthesized counts equal
    native's. The pre-implementation captures are unchanged.
  * DRAW: the three original and three derived cases (with `pcylinder_counts`)
    pass on the Rust adapter and native DRAW (OCCT 7.9.3 locally, the Ubuntu
    24.04 package in CI).
  * Effects on accepted work: every row, relation and corpus digest of the
    284 polygon-only cases is byte-identical to `26fc457f` apart from the
    added region rows; the explicit circle cases differ only by their removed
    seam rows and relations and the added region. Mass properties, bounds and 5×5×5 classifications are
    bitwise equal to the pre-migration revision on macOS, Linux and Windows.
  * Fixtures: 64 B-rep reports (23 valid), 46,494 entities and 114,270
    relations over 546 extrusions, all regenerated byte-identically on
    Python 3.9 and 3.12.
  * Clean local 600-second campaigns at `77401b00` (AddressSanitizer,
    standard 20-second/2 GiB limits): `brep_validation` completed 600.09
    seconds of mutation after 105.6 seconds of replay (18,173 executions,
    7,351 edges, 1,222 MB peak); `identity` 600.03 seconds after 148.1
    (6,004 executions, 7,387 edges, 864 MB); `history` 600.03 seconds after
    103.2 (9,334 executions, 7,965 edges, 853 MB). There were no crash,
    timeout, OOM, slow-unit or disagreement artifacts.

  What remains: faces without loops, windings in `v`, unwound loops on wound
  faces (only `uncertified_containment` today), free and non-manifold
  results, and degenerate-edge counts for poles wait for the surfaces and
  operations that need them; the `.brep` interop, native selector and
  coverage ledger are T2. The `continuity` check of D11 (C1 at interior
  knots of every edge and face) is not implemented: lines, arcs, planes and
  cylinders are C1 by construction, so nothing is unchecked today, but it
  must land, with its fixtures and a mutation, before the first spline edge
  or face enters a topology. Composed histories record only their last
  step's level; see H8 in `IDENTITY_AND_HISTORY.md`.
* **T2 — accepted at `0913b5e2`** (first pushed at `cb57697b`, whose kernel
  workflow also passed; `0913b5e2` adds M5, which gives imported vertices
  and fins their OCCT tolerances as enclosures).
  * Rust kernel workflow: all twelve jobs passed. On Linux the interop
    bridge certified the independent reader on all 77 corpus solids. OCCT
    read back all 546 prisms and all 29 imported corpus solids as valid,
    with equal counts and properties (worst difference `9.5e-15`, bound
    `1e-11`); 575 matches, no reviewed differences, no failures. On macOS the
    worst difference was `2.5e-14`. The original-test bridge passed
    `explode_selector` on the Rust adapter and on Ubuntu's DRAW. The
    coverage ledger equals the recorded one.
  * Fuzzing workflow: all twenty-four targets passed. On Linux `brep_io`
    replayed 423 inputs in 110.0 seconds, then ran 60.08 seconds of
    mutation (539 executions, 9,540 coverage edges, 664 MB RSS peak) without
    an artifact.
  * Fixtures: `brep-io-expected.tsv` (37 files, 77 solids, 29
    representable) regenerates byte-identically on Python 3.9 and 3.12, and
    `occt_brep.rs` agrees with it on every file and solid.
  * Order of evidence: the native corpus observations were captured after
    implementation, as recorded in `fixtures/occt-brep-io-capture/NOTES.md`.
  * Clean local 600-second campaign at `0913b5e2` (AddressSanitizer,
    standard 20-second/2 GiB limits): `brep_io` completed 600.03 seconds of
    mutation after 40.6 seconds of replay, with 15,855 executions, 10,163
    coverage edges and a 997 MB RSS peak, and no crash, timeout, OOM,
    slow-unit or disagreement artifacts. The earlier smoke campaign's crash
    (a reference beyond its table) is fixed and retained as a regression.

  What remains: no original upstream assertion is mapped-and-verified yet
  (0 of 35,766). Data-dependent original cases wait for OCCT's external
  test data, and the adapter's `restore` for the first of them. The selector
  handles faces and edges only. Surfaces with poles, and so degenerate
  edges, are not representable, and neither are splines, cones, spheres and
  tori.

* **S3 cone — accepted at `166fc905`** (the cone is `0c94aa53`; its native
  `MakeCone` capture `dde086c1` and validation rules `1a76d29e` came before
  it). The first push at `0c94aa53` failed only the revolve bridge's Linux
  fingerprint (reviewed in `166fc905`) and the `identity` and
  `surface_knots` replay budgets (R12); `166fc905` also shares the mass
  engine's trigonometric moments.
  * Both workflows passed at `166fc905`: twelve kernel jobs (the primitive
    bridge 22 of 22 cones against the independent reference, worst
    `4.9e-15`; the revolve history bridge 21 matches and one reviewed
    difference, `nearly_cylinder`; the B-rep bridge 44/8/0) and twenty-four
    fuzz targets.
  * Clean local 600-second campaigns at `166fc905` (AddressSanitizer,
    standard 20-second/2 GiB limits): `brep_validation` 9,644 mutation
    executions after 140 s of replay (10,112 edges, 684 MB peak),
    `identity` 2,013 after 24 s (10,074, 734 MB), `history` 7,493 after
    8 s (11,173, 815 MB), `brep_io` 10,091 after 6 s (11,996, 767 MB); no
    artifact.
  * Order of evidence: the `MakeRevol` history capture came after the
    builder was written (R11 of `REVIEW_NOTES.md`, accepted by the user).
* **S3 sphere — accepted at `059616f9`** (the sphere is `d7677e4b`, after
  its native captures `3641a2f7`; `059616f9` adds only the torus's native
  captures and tools, no kernel or fuzz source).
  * `d7677e4b`'s kernel workflow passed and its fuzzing run was cancelled by
    the next push; both passed at `059616f9` (twelve and twenty-four), with
    the sphere primitive bridge 22 of 22 (worst `6.7e-15`) and the sphere
    revolve bridge 22 of 22.
  * Clean local 600-second campaigns at `d7677e4b` (the same kernel and fuzz
    sources): `brep_validation` 12,537 executions after 62 s of replay
    (11,064 edges, 767 MB), `identity` 3,645 after 16 s (10,948, 829 MB),
    `history` 7,700 after 9 s (12,086, 880 MB), `brep_io` 11,797 after 5 s
    (13,151, 805 MB); no artifact.
  * Upstream: `bugs/modalg_6/bug27264_2` (a restored whole sphere) passes on
    both backends and is registered; the ledger has 2 mapped-and-verified
    assertions.
* **S3 torus — accepted at `ac92ec89`, with S4** (the torus is `e2e7c200`,
  after its native captures `059616f9`). `e2e7c200`'s kernel workflow
  failed only on the Linux observation of OCCT's inside-out inner half,
  reviewed in `896606bc`, where both workflows passed (twelve and
  twenty-four; the torus primitive and revolve bridges 22 of 22 each).
  `896606bc` already carried R4, so the torus's clean campaigns at
  `e2e7c200` (`brep_validation` 17,019 executions, `identity` 2,964,
  `history` 6,290, `brep_io` 9,345; no artifact) predate that code, and the
  gate closes at `ac92ec89` with S4's evidence below.
* **R4 (S4) — accepted at `ac92ec89`.** The continuity check is `395bc3e8`;
  both workflows passed at its descendant `896606bc`, which adds only the
  S4a spline capture and a torus review. Its fixtures (nineteen spline cells,
  mutation 31) and campaigns are part of S4's record.
* **S4 — accepted at `ac92ec89`** (S4a capture `85b324e6` and the S4d
  `BRepGProp` capture `fc389e8d`, each before the kernel code they observe;
  S4b-d `38725497`, `e7c2a682`; S4b-e `79422a58`; spline mass `ac92ec89`).
  * Rust kernel workflow: all twelve jobs passed (run 36302492680). On Linux
    the spline bridge had 7 matches, 3 reviewed differences and no failures,
    reproduced the ten S4a and ten S4d rows, and found OCCT's properties in
    the kernel's enclosures for all five models valid on both sides; the
    prism bridge kept 44/8/0 and the interop bridge 622/0/0 (58 written
    corpus solids, worst property difference `4.4e-14`).
  * Fuzzing workflow: all twenty-four targets passed (run 36302492715). On
    Linux `brep_validation` replayed 2,371 inputs in 920 s, then ran 60.08 s
    of mutation (221 executions, 17,559 edges, 1,001 MB); one replayed input
    was a 14-second slow-unit diagnostic, triaged and retained after this
    revision (`fuzz/regressions/README.md`: about 4.1 s locally under
    AddressSanitizer once the mass check reuses the report's topology).
    `brep_io` replayed 1,659 in 470 s, then 60.09 s (300 executions, 18,573
    edges, 856 MB).
  * Fixtures: 152 B-rep reports (63 valid), the spline masses of the sixteen
    valid spline cases and `brep-io-expected.tsv` (77 solids, 60
    representable, 58 certified) regenerate byte-identically on Python 3.9
    and 3.12.
  * Clean local 600-second campaigns at `ac92ec89` (AddressSanitizer,
    standard 20-second/2 GiB limits, `brep_validation` with the allocator
    purge): `brep_validation` 4,266 mutation executions after 249 s of
    replay (17,856 edges, 1,285 MB), `identity` 2,730 after 218 s (12,254,
    990 MB), `history` 6,224 after 110 s (12,985, 938 MB), `brep_io` 8,318
    after 109 s (18,805, 1,137 MB); no artifact. Earlier smoke campaigns
    found the far, small rational corner (fixed in `ac92ec89`, fixture
    `spline_rounded_corner_far`) and the RSS growth that moved
    `brep_validation` to the allocator targets.
  * Decisions taken during S4 (in `REVIEW_NOTES.md`): periodic spline
    surfaces and seams on spline surfaces stay unsupported; spline mass
    enclosures are first order where not exact.

## Delivery rules

The rules in `IDENTITY_AND_HISTORY.md` apply unchanged: native observations
before implementation, subject-only commits without attribution footers or
absolute paths, ask before pushing, AddressSanitizer seeds under 20 s ÷ 2.6,
mpmath rounding in fixture generators, no stripped upstream tests counted as
passes, no widened tolerances, no reviewed timeouts, and no behaviour change
to an existing operation without a new algorithm level.
