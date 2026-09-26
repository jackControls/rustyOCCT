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
[Acceptance](#acceptance)); T2 accepted at `0913b5e2`.

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
  exact spline modules already detect interior-knot discontinuities.
  *Forbids:* a single edge across a corner.
* **D12. Blend surfaces are an open decision with a constraint.** Procedural
  blend surfaces evaluated through their supports are preferred, as in
  Parasolid and CGM, but undecided. The `Surface` representation must remain
  extensible with a certified evaluation contract. *Forbids:* approximating a
  fillet to a B-spline as the only representation.

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
  `Error::UnknownAlgorithmLevel`. `History::then` gives a composite the level
  of its last step, which is exact only while every step ran at one level
  (open item under H8).
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

## Delivery rules

The rules in `IDENTITY_AND_HISTORY.md` apply unchanged: native observations
before implementation, subject-only commits without attribution footers or
absolute paths, ask before pushing, AddressSanitizer seeds under 20 s ÷ 2.6,
mpmath rounding in fixture generators, no stripped upstream tests counted as
passes, no widened tolerances, no reviewed timeouts, and no behaviour change
to an existing operation without a new algorithm level.
