# Identity, history, attributes and error bounds

This is a design contract and delivery guide, not an acceptance record. Nothing
in it is implemented yet. It defines five kernel contracts that must exist
before any topology-changing operation (Booleans, fillets, shells, splits) is
written, and it tells the implementing agent what to build, in what order, and
what evidence promotes each milestone. Read `PORTING.md`, `PRODUCTION_READINESS.md`
and `BREP_VALIDATION.md` first; this guide follows the same shipping pattern.

The long-term goal of the kernel is to rival Parasolid, not to reproduce
OCCT. Parasolid's advantage is not feature count but the contracts around its
features: every operation reports how entities were created, split, merged,
deleted or transferred; attributes ride on entities with declared behaviour;
identifiers survive transmit and receive. OCCT lacks all three at kernel level:
a `TopoDS_Shape` is the same shape only while it holds the same `TShape` pointer
and location (`TopoDS_Shape::IsSame`), history is optional per algorithm and
never verified for completeness, and naming lives in the application framework
as a heuristic re-selection layer. These are retrofit failures. The kernel has
one builder and no Booleans, so the decisions below are still cheap.

**Status:** M0, M1 and M2 accepted at `34efd36c` (value ids, complete
checked histories and the attribute checker, on extrusions and rigid
transforms). The topology model is decided in `TOPOLOGY_MODEL.md`; its
migration, T1, precedes M3. M3–M5 pending. Acceptance evidence per milestone
is under [Acceptance](#acceptance).

## What the kernel promises, and what it does not

The kernel provides the raw material for persistent naming. It does not solve
persistent naming, and neither does Parasolid; SolidWorks, NX and Onshape each
build naming policy on top of kernel tracking. The boundary:

| Kernel promises | Application decides |
| --- | --- |
| Ids are values, deterministic, never reused within a body's lineage | Which of several split children a user selection now means |
| Every operation returns a complete, independently checkable history | Whether a broken reference is an error or a re-selection prompt |
| Split and merge are reported as relations, never as a chosen survivor | Feature-level naming, display names, undo semantics |
| Attributes propagate by a declared policy, and conflicts are reported | What attributes mean (colour, material, PMI, sketch links) |
| Every constructed coordinate carries a computed error bound | Presentation tolerances, units, rounding for display |

## Definitions

* **Body** — one `Solid` (later: any B-rep body). **Lineage** — the sequence
  of bodies produced from an origin construction through operations. Ids are
  unique within a lineage, not merely within one body.
* **Entity** — a vertex, edge, face, or later a loop, shell or region. Every
  entity has an **id** (a value, see below) and a **slot** (its arena index in
  the owning `Topology`). Slots are an implementation detail and remain
  body-local and dense; ids are the public identity.
* **Operation** — a kernel call that produces one or more output bodies from
  zero or more input bodies plus parameters. Constructions (extrude) have no
  input bodies; their inputs are labelled profile entities.
* **Derivation** — the structural record from which an id is computed:
  operation, parent ids, role and canonical ordinal.
* **History** — the complete list of relations between the input entities and
  output entities of one operation.
* **Entities with identity** — bodies, bounded regions, faces, edges and
  vertices. Shells, loops, fins and the infinite void are structure: never
  named, never journaled (`TOPOLOGY_MODEL.md`).
* **Algorithm level** — the version of an operation's behaviour. A history
  records the level it ran at, and the kernel replays an operation at a
  recorded level. Levels never enter a derivation, so upgrading an algorithm
  never changes ids.
* **Resolution** — a body's modelling tolerance, today `Tolerance::linear` and
  `Tolerance::angular` (`math.rs`). It is a contract about what counts as
  distinct. **Enclosure** — a computed bound on where the true geometry lies
  relative to its representable value. Enclosures are outputs of arithmetic;
  resolutions are inputs to decisions.

## Contract 1: value identity

### Requirements

* **I1. Ids are values.** `EntityId` is `Copy`, `Eq`, `Ord`, `Hash`,
  serializable, and meaningful without access to the process that created it.
  No pointers, no reference counting, no interior mutability.
* **I2. Ids are derived structurally, not positionally.** An id is a fixed,
  documented 128-bit digest of its `Derivation`. Two runs, on any platform,
  of the same operation on inputs with the same ids and parameters produce the
  same ids. No binary64 value enters a derivation except through an exact
  canonical ordering (I5).
* **I3. Caller labels are the roots.** Construction inputs accept optional
  caller labels (`InputLabel(u64)`) on profile boundaries, segments and
  vertices. When present they are the parents in derivations; when absent,
  indices are. A sketch edit that keeps labels but moves points keeps every id.
  Inserting a segment changes only ids derived from it.
* **I4. Ids are never reused.** Within a lineage, an id names one entity for
  ever. A deleted id is never reassigned. A collision between two derivations
  is an operation error (`Error::InvalidTopology("id collision")`), never
  resolved silently by ordinal.
* **I5. Canonical ordinals are exact.** When one parent yields several
  children of the same role (a face split into three), the ordinal comes from
  an exact comparison of rational quantities defined per operation (for the
  height split in M3: the axial order of the pieces). The documentation of
  each operation states its ordinal rule. Ordinals are deterministic for
  identical inputs; they may change under perturbation, and that is reported
  through the history, not hidden.
* **I6. Rigid transforms and copies keep ids.** `Solid::transformed` reports
  every entity as `Modified` with the same id.
* **I7. Slots stay.** `Topology` keeps its arena and adds `id_of(slot)` and
  `slot_of(id)` maps. Existing index-based APIs continue to work.

### Data model sketch

```rust
/// 128-bit digest of a Derivation, computed by a fixed in-crate hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId(pub [u8; 16]);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InputLabel(pub u64);

/// Why an entity has its id. Stored beside the id, never recomputed from geometry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Derivation {
    pub operation: OperationId,      // caller-supplied, e.g. a feature id
    pub kind: OperationKind,         // Extrude, Transform, HeightSplit, StackedFuse, ...
    pub parents: Vec<Parent>,        // Label(InputLabel) | Entity(EntityId), canonical order
    pub role: Role,                  // StartCap, EndCap, Wall, BottomEdge, TopEdge, Vertical, Seam, CutFace, ...
    pub ordinal: u32,                // exact canonical ordinal among same-role siblings
}
```

The hash is FNV-1a with a 128-bit state, implemented in the crate (the
standard library's `DefaultHasher` is not stable across Rust versions). Its
weak collision resistance is acceptable because inputs are short structured
records and collisions are detected (I4). Document the byte encoding of
`Derivation` beside the hash and pin both with fixed test vectors. The
`Derivation` is retained in the body so that "why does this face have this id"
is always answerable, and so that collisions are detectable (I4).

## Contract 2: complete, validated history on every operation

### Requirements

* **H1. Every operation returns a `History`.** No operation may construct an
  output body without one. The signature pattern is
  `fn op(...) -> Result<(Outputs, History)>`.
* **H2. Every input entity appears exactly once as a source.** Its relation is
  one of `Unchanged`, `Modified`, `Split`, `Merged` or `Deleted`.
* **H3. Every output entity appears at least once as a target.** Targets of
  `Generated` may have several parents (an intersection edge has two faces).
  An output with no relation is a contract failure, never "unknown".
* **H4. Relations respect dimension and support.** `Unchanged` requires
  identical geometry and identical loop structure. `Modified` keeps the id and
  the surface or curve family. `Split` children and `Merged` parents share the
  parent's support geometry up to the operation's rigid motion. `Generated`
  parents are of lower or equal dimension. `Deleted` ids appear in no output.
* **H5. Histories compose.** `History::then(&next) -> History` is associative
  and yields a valid history from the first inputs to the last outputs. A
  split followed by a merge of the same pieces composes to `Modified`, not to
  `Unchanged`, because geometry may have changed.
* **H6. An independent checker validates every history.**
  `history::check(inputs, outputs, &history) -> Vec<HistoryIssue>` follows the
  `Topology::check` style: sorted, duplicate-free, typed kinds, each with its
  offending id. It knows nothing about the algorithm that produced the history.
  Every operation runs it in debug builds and in every test.
* **H7. The application resolves through the history only.** There is no API
  that answers "which output face is the old face" other than
  `History::resolve(id) -> Resolution` (Contract 3).
* **H8. Every operation is versioned.** A `History` records the
  `AlgorithmLevel` of the operation that produced it. Changing an existing
  operation's behaviour requires a new level; the previous level stays
  callable so a stored history replays identically. Levels are not part of
  `Derivation` (adopted from CGM's software configurations; pending, added in
  T1).

### Types

```rust
pub enum Relation {
    Unchanged { id: EntityId },
    Modified  { id: EntityId },
    Generated { from: Vec<Parent>, to: EntityId, role: Role },
    Split     { from: EntityId, into: Vec<EntityId> },   // ids ordered by canonical ordinal
    Merged    { from: Vec<EntityId>, into: EntityId },
    Deleted   { id: EntityId },
}

pub struct History {
    pub operation: OperationId,
    pub kind: OperationKind,
    pub input_bodies: Vec<BodyRef>,
    pub output_bodies: Vec<BodyRef>,
    pub relations: Vec<Relation>,   // sorted by (source ids, kind, targets)
    pub attributes: Vec<AttributeOutcome>,  // Contract 4
}

pub enum HistoryIssueKind {
    MissingSource, DuplicateSource, MissingTarget, DanglingId, DimensionMismatch,
    UnchangedGeometryDiffers, ModifiedFamilyDiffers, SplitSupportDiffers,
    MergedSupportDiffers, DeletedStillPresent, OrdinalNotCanonical, CompositionInvalid,
}
```

### History of the existing operations

`Solid::extrude` (construction, `OperationKind::Extrude`), with the profile's
labelled entities as parents:

| Output | Relation |
| --- | --- |
| Start cap, end cap | `Generated { from: [every boundary label], role: StartCap / EndCap }` |
| Wall face of segment *s* of boundary *b* | `Generated { from: [Label(b, s)], role: Wall }` |
| Bottom and top edge of segment *s* | `Generated { from: [Label(b, s)], role: BottomEdge / TopEdge }` |
| Vertical edge at vertex *v* | `Generated { from: [Label(b, v)], role: Vertical }` |
| Bottom and top vertex at *v* | `Generated { from: [Label(b, v)], role: BottomVertex / TopVertex }` |
| Circle boundary: circular edges, seam edge, cylinder wall, seam vertices | Same pattern with roles `BottomEdge`, `TopEdge`, `Seam`, `Wall`, `SeamVertex` |

`Solid::transformed` (`OperationKind::Transform`): every entity `Modified`
with the same id, in slot order.

After T1 (`TOPOLOGY_MODEL.md`) the circle row changes: the circular edges are
ring edges with roles `BottomEdge` and `TopEdge`, there is no seam and no seam
vertex, and the solid region is `Generated` from every boundary label with
role `Region`. Every polygon id and relation is unchanged by T1.

This replaces `FaceOrigin` as the source of provenance. Keep `FaceOrigin`
until M1 lands, then derive it from the `Derivation` and remove the field.

## Contract 3: explicit ambiguity as a type

* **A1. No survivor selection in the kernel.** A `Split` never names a
  primary child and a `Merged` never names a primary parent. There is no
  "modified" record for a split parent.
* **A2. Resolution is a closed enum.**

```rust
pub enum Resolution {
    Same(EntityId),                 // Unchanged or Modified
    Split(Vec<EntityId>),           // in canonical ordinal order
    Merged { into: EntityId, with: Vec<EntityId> },
    Deleted,
    Unknown,                        // the id was not an input of this history: caller error
}
```

  `Unknown` is for ids that were never inputs; it is a caller error, not a
  kernel answer. A composed history resolves through every step.
* **A3. Attribute conflicts are typed, not dropped silently.** See Contract 4.
* **A4. Geometric coincidence is never identity.** Two faces with equal
  surfaces, equal areas and equal centroids are two ids. Any API that would
  match entities by geometric equality lives in the application, with its own
  documented tolerance, not in the kernel.

## Contract 4: attribute propagation policy

* **P1. Attributes are opaque to the kernel.** An `Attribute` is an
  `AttributeKey(u64)` and a byte value. The kernel never interprets values.
* **P2. Every key has a declared policy.** The operation context carries a
  `PolicyTable`. A key without a policy is an operation error.

```rust
pub struct Policy {
    pub on_modify:    Keep | Drop,
    pub on_transform: Keep | Drop | Recompute,   // Recompute: application callback, kernel records it
    pub on_split:     CopyToAll | Drop,
    pub on_merge:     KeepIfAllEqual | Drop,
    // on_delete is always Drop; on_generate never creates attributes
}
```

  Parasolid's attribute definitions declare comparable behaviour per class;
  the names above are ours.
* **P3. Outcomes are recorded.** Every attribute on an input entity yields one
  `AttributeOutcome { key, from, to: Vec<EntityId>, result: Kept | Copied | Dropped | Conflict }`
  in the history. `Conflict` is reported when `KeepIfAllEqual` sees unequal
  values, and the attribute is then absent from the output.
* **P4. Deterministic.** Outcomes depend only on the history and the table.
* **P5. Checked independently.** `attributes::check(inputs, outputs, &history, &table)`
  verifies every output attribute is traceable to an input attribute and an
  outcome, no attribute survives on a deleted id, and no outcome contradicts
  its policy.

## Contract 5: tolerance as computed error bounds

* **T1. One resolution per body.** `Tolerance` is fixed at construction and
  stored on the body. Operations on bodies with different resolutions are an
  error until an explicit, documented combination rule exists.
* **T2. Enclosures are data.** Each vertex carries a radius bound; each edge a
  bound on the distance between its curve and every use's pcurve image; each
  face a bound on its uses' closure gaps. They are computed with the certified
  tiers in `certified.rs` from the construction's exact inputs and the chosen
  representable outputs, exactly as the intersection modules already return
  minimal binary64 enclosures.
* **T3. Valid means enclosures fit the resolution.** `Topology::check` gains
  `enclosure_exceeds_resolution` and `enclosure_missing`. The certified checks
  it already performs (`vertex_off_curve`, `pcurve_off_edge`, `uv_gap`) then
  verify the stored bounds rather than recomputing from scratch, and a stored
  bound that the recomputation shows to be too small is `enclosure_unsound`.
* **T4. Never widen to succeed.** No operation may increase a resolution or
  an enclosure to make a step pass. If a construction cannot be represented
  within the resolution, it fails with `Error::PrecisionLoss` or
  `Error::Unrepresentable`, which already exist. Review every comparison of a
  binary64 value against `tolerance.linear()` and route it through a certified
  tier or document why it is exact.
* **T5. Enclosures propagate.** An output enclosure is a function of input
  enclosures plus the construction's own error, never a reset to zero. Rigid
  transforms add the transform's rounding error.
* **T6. OCCT tolerances are observations, not targets.** OCCT's per-entity
  tolerances (`BRep_Tool::Tolerance`) are requested upper bounds that only
  grow. The native bridge records them for the same constructions as
  observations; a Rust enclosure larger than OCCT's tolerance for an exactly
  representable construction is a bug to investigate, never a reason to
  loosen the bound.

## Topology model: decided

The model is decided in `TOPOLOGY_MODEL.md` (D1–D12): a cellular partition of
space with regions including an infinite void, stored adjacency with fins in
radial order, ordered loops with vertex loops and zero-loop faces, no seams,
stored certified pcurves per fin, sense on entities with matter derived from
regions, one resolution per body with enclosures, no instancing inside
bodies, computed body types, immutable bodies, C1 within a cell, and blend
surfaces left open with an extensibility constraint. The contracts here are
unchanged by it: ids attach to bodies, bounded regions, faces, edges and
vertices; relations are generic over those kinds; shells, loops and fins have
no identity. Its migration, T1, precedes M3.

## Milestones

Each milestone ships the full pattern: kernel code, an independent Python
oracle generator (`rust/tools/generate_*_fixtures.py --check`, run with
`target/math-oracle-venv/bin/python`), in-crate fixture tests, a fuzz target
registered in `run_fuzz.py` `TARGETS` with seeds and added to the
`.github/workflows/rust-fuzz.yml` matrix, a source-pinned native bridge
(`compare_*.py` plus a reviewed divergence JSON) where OCCT can observe the
behaviour, and doc updates (this file, `PORTING.md`, `PRODUCTION_READINESS.md`,
`VALIDATION.md`, `FUZZING.md`, `SOURCE_MAP.md`, and `MATHEMATICS.md` when
arithmetic changes). Native observations are captured **before** implementation
and checked in under `rust/fixtures/occt-*-preimplementation/`.

### M0 — types and checkers, no behaviour change

* Add `identity.rs` (`EntityId`, `InputLabel`, `Derivation`, hash, encoding),
  `history.rs` (`Relation`, `History`, `Resolution`, `check`, `then`) and
  `attributes.rs` (`Attribute`, `Policy`, `PolicyTable`, `AttributeOutcome`,
  `check`).
* Unit tests for hash stability (fixed vectors checked into
  `rust/fixtures/identity-vectors.tsv`), encoding round trips, checker issue
  kinds on hand-written histories, composition laws on hand-written chains.
* No operation is changed. Gate: kernel CI green, clippy clean at MSRV 1.85,
  `unsafe_code = "forbid"` untouched.

### M1 — value ids on extrude and transform

* `Profile` accepts optional labels; `Topology` stores ids, derivations and
  the slot maps (I7); `extrude` and `transformed` derive ids per the tables
  above; `FaceOrigin` becomes derived.
* **Fixtures:** `generate_identity_fixtures.py` builds the same 256-prism
  corpus as `invariants` plus labelled variants and writes expected ids from
  its own implementation of the encoding and hash. Cases must include:
  polygons of 3 to 12 sides, circles, 0 to 3 holes, both extrusion
  directions, all eight sign combinations of `box_at`, rotated and translated
  copies, labels absent, labels present, labels permuted.
* **Fuzz target `identity`:** structure-aware (`arbitrary` input: profile,
  labels, frame, offsets, transform). Invariants: building twice gives
  identical ids; a transform keeps ids; permuting label values permutes ids
  bijectively; no two entities share an id; every id's derivation re-hashes to
  the id; slot maps are inverse.
* **Native bridge:** not applicable. OCCT has no value ids; record this in
  `VALIDATION.md` rather than inventing a comparison.

### M2 — complete history on extrude and transform

* `extrude` and `transformed` return `(Solid, History)`; the old signatures
  become thin wrappers that discard the history, marked deprecated.
* **Fixtures:** `generate_history_fixtures.py` writes the complete expected
  relation list per case from an independent enumeration of the prism
  structure, not from Rust. Include every case of M1.
* **Fuzz target `history`:** the M1 input plus a mutation of the produced
  history that must be detected: drop a relation, duplicate one, change a kind,
  point at a dangling id, swap split ordinals, compose out of order. A clean
  history must report no issues; each mutation must report at least its
  predicted kind on the mutated relation.
* **Native bridge:** `occt_history_oracle.cpp` builds the profile as a face
  with `BRep_Builder`, runs `BRepPrimAPI_MakePrism`, wraps it in
  `BRepTools_History`, and prints `Generated`, `Modified`, `IsRemoved`,
  `FirstShape` and `LastShape` per input subshape in a stable order.
  `BRepBuilderAPI_Transform::ModifiedShape` covers transforms. `compare_history.py`
  maps native records to Rust relations (native `Generated(vertex)` is our
  `Vertical` edge, `Generated(edge)` our `Wall` face, first and last shapes
  our caps). Review differences in `occt-history-divergences.json` with the
  usual fingerprint rules; timeouts and crashes are never reviewable; deadline
  120 s.
* **Upstream tests:** DRAW's `prism` records history (`BRepTest_SweepCommands.cxx`),
  so `savehistory`, `generated`, `modified` and `isdeleted` adapters can
  delegate to the kernel history. Add the four commands to the Rust DRAW
  adapter. Where an original self-contained case exercises them, register it
  in `upstream-draw.json`; otherwise add derived cases labelled as such, never
  stripped originals.

### T1 — topology model migration

Specified in `TOPOLOGY_MODEL.md`. It replaces seams with periodic loops,
adds regions, shells with sides, fins and ring edges, and must keep every
polygon id and relation byte-identical. M3 is built on the migrated model.

### M3 — first topology-changing operations: height split and stacked fuse

These are deliberately narrow. They exist to exercise `Split`, `Merged`,
`Generated` and `Deleted` on real bodies with an OCCT oracle, before any
general Boolean.

* `Solid::split_at_height(h)` splits a prism by the plane at offset `h`
  strictly between its offsets: the solid region and each wall face and
  vertical edge `Split` into two (ordinal by axial order), caps and ring
  edges `Unchanged` in their body, two cut faces `Generated` from the plane
  and the wall set, cut edges `Generated` from walls (a ring edge for a
  cylinder wall, with no vertices), cut vertices `Generated` from vertical
  edges, the input body deleted and two bodies output.
* `Solid::fuse_stacked(other)` merges two prisms that share a cap exactly
  (same profile, same frame, adjacent offsets): the two solid regions and,
  pairwise, walls and vertical edges `Merged` (parents in body order,
  ordinal 0), the shared caps `Deleted`, the shared cap edges and vertices
  `Deleted`, outer caps and outer ring edges `Unchanged`.
* **Fixtures:** independent enumeration of expected relations and of the
  resulting `Topology` structure; every result must pass `Topology::check`.
* **Fuzz target `split_merge`:** random prism, random split height, then
  fuse back; invariants: histories check clean; `split.then(fuse)` resolves
  every wall to `Same`; volumes are conserved within the certified enclosure;
  `Resolution` for every input id is exhaustive.
* **Native bridge:** `BRepAlgoAPI_Splitter` with a plane face tool, and
  `BRepAlgoAPI_Fuse` followed by `ShapeUpgrade_UnifySameDomain` with
  `History()`. Compare relation kinds and counts; OCCT's `Modified` for a split
  parent is our `Split` (reviewed divergence, by design).
* **Upstream tests:** the self-contained history cases in `tests/bugs/modalg_6`
  and `modalg_7` that use `bsplit`, `bfuse`, `unifysamedom` with history
  commands are the candidates; register those the adapter can run.

### M4 — attributes on every operation

* Attribute storage on `Topology`, `PolicyTable` in the operation context,
  outcomes in `History`, for extrude, transform, split and fuse.
* **Fixtures:** the full policy matrix (each `on_*` variant × each relation
  kind) on the M3 corpus, expected outcomes from `generate_attribute_fixtures.py`.
* **Fuzz target `attributes`:** random attributes with random policies through
  random operation chains; invariants: `attributes::check` clean; conflicts
  appear exactly when unequal values merge under `KeepIfAllEqual`; no
  attribute on a deleted id; outcomes are deterministic under re-run.
* **Native bridge:** not applicable at kernel level. Record in `VALIDATION.md`.

### M5 — enclosures as entity data

* Add enclosure fields, compute them in every builder and operation, extend
  `Topology::check` per T3, and audit every tolerance comparison per T4.
* **Fixtures:** the B-rep validation corpus with expected enclosure bounds from
  `brep_reference.py`, and perturbation cases just inside and just outside the
  resolution.
* **Fuzz target:** extend `brep_validation` with enclosure mutations (a bound
  set too small must be `enclosure_unsound`; a construction placed so its
  enclosure exceeds the resolution must fail construction, never validate).
* **Native bridge:** extend `occt_brep_check_oracle.cpp` to print
  `BRep_Tool::Tolerance` for each subshape; record as observations per T6.

## Implementation decisions

The contracts above are normative; the data model sketches are not. Where the
implementation refines a sketch, the refinement and its reason are recorded
here.

* **Encoding.** `Derivation` also carries the entity kind, so a vertex and an
  edge derived from the same parent with the same role and ordinal can never
  share an id. Version 1 of the encoding is documented in `identity.rs` and
  pinned by `fixtures/identity-vectors.tsv` (11 vectors, plus the published
  FNV-1a-128 vectors for `""`, `"a"` and `"foobar"`).
* **Unlabelled parents.** Without labels, a parent is
  `Parent::Profile { boundary, element }` with the boundary's *stored* index:
  polygons are stored counter-clockwise with the first input point kept, so
  labels given in the caller's input order are mapped through that reversal.
* **Start and end.** Roles `BottomEdge`, `TopEdge`, `BottomVertex` and
  `TopVertex` mean the extrusion's start and end sides, so swapping the offsets
  keeps every id. The two seam vertices of a circle share role `SeamVertex`,
  ordinal 0 on the start side and 1 on the end side.
* **`Modified { from, to }`.** H5 requires a split followed by a merge of the
  same pieces to compose to `Modified`, but the merged entity has a new id (I4
  forbids reusing the split parent's). So `Modified` names both ids. A single
  operation must keep the id (`modified_id_changed` otherwise); only a
  composed history (`OperationKind::Composite`) may change it.
* **Generated dimension rule.** H4 says `Generated` parents have lower or equal
  dimension, but M3 generates cut edges from walls and cut vertices from
  vertical edges. The checker enforces the rule both cases satisfy: a
  generated entity's dimension is at most one more than each parent's (a
  sweep adds one dimension; an intersection lowers it).
* **Checker input.** `history::check` takes `EntitySet`s (per body: its id,
  resolution, and each entity's kind, ordinal, geometry and loop structure),
  so hand-written histories can be checked before any operation produces
  them. Split and merge supports are compared exactly: equal surfaces, equal
  circles, or line pieces whose endpoints lie within the resolution of the
  whole's infinite line (exact rational arithmetic). The split ordinal check
  applies to single operations only; composed children keep their own
  ordinals.
* **Extra issue kinds.** `modified_id_changed`, `invalid_arity` (a split into
  fewer than two, a merge of fewer than two, a generation from nothing),
  `body_mismatch`, `duplicate_id` and `relation_order` (relations not in
  canonical order) join the contract's list.
* **Composition limits.** `History::then` maps each input through both
  histories. A many-to-many relation (for example, a split piece merged with
  another input) has no single-relation form and is `composition_invalid`.
  Attribute outcomes are not composed.
* **Attributes.** `OutcomeResult::Recomputed` records `OnTransform::Recompute`.
  A merge keeps an attribute only when every parent carries the key with the
  same value; a missing value counts as unequal. The checker adds
  `missing_attribute` (an outcome that keeps or copies has no value on its
  target) and `attribute_on_dropped` (a value survives where every outcome
  dropped it).

## Acceptance

Promote a milestone only when all applicable gates pass at one clean
revision, and record that revision, the kernel CI job count, the fuzz CI
target count, the native bridge match/review/failure counts on macOS and
Linux, and a clean local 600-second campaign here, in the style of
`BREP_VALIDATION.md`. A finite corpus does not make the contract complete;
say what remains.

* **M0 — accepted at `d388158a`.** `identity.rs`, `history.rs` and
  `attributes.rs` with their checkers; no operation changed. The Rust kernel
  workflow passed all ten jobs (Linux, macOS and Windows debug/release, Rust
  1.85, WebAssembly, the original-test bridge and four source-pinned native
  comparisons) and the fuzzing workflow all nineteen targets. Evidence: the
  11 fixed vectors in `identity-vectors.tsv` (produced by the independent
  `identity_reference.py`) plus the published FNV-1a-128 vectors, decoding and
  malformed-input tests, and `history_contracts.rs`, which triggers every
  history issue kind on hand-written histories, checks composition laws
  (identity steps, associativity, split-then-merge is `Modified`, many-to-many
  is rejected) and the full attribute policy matrix. No native bridge applies.
* **M1 — accepted at `34efd36c`**, together with M2 (the ids are unchanged by
  M2; that revision carries both). `identity.rs` in the test suite matches
  the independent ids, roles, ordinals and parents of 46,734 entities in 546
  extrusions, each entity located from geometry alone, before and after every
  rigid motion. The `identity` fuzz target recomputes every id with its own
  encoder. No native bridge applies: OCCT has no value ids (recorded in
  `VALIDATION.md`).
* **M2 — accepted at `34efd36c`.**
  * Rust kernel workflow: all eleven jobs passed, including the new
    source-pinned history comparison.
  * Fuzzing workflow: all twenty-one targets passed. On Linux, `identity`
    replayed 97 inputs in 70 seconds, then ran 60.08 seconds of mutation
    (881 executions, 6,178 coverage edges, 641 MB RSS peak); `history`
    replayed 113 inputs in 71 seconds, then ran 60.07 seconds (2,730
    executions, 6,933 edges, 630 MB). Neither produced an artifact.
  * Native history bridge (`compare_history.py`, OCCT 8.1.0 built from the
    pinned source): 98 matches, 0 reviewed differences and 0 failures on
    both macOS and Linux; the slowest Linux case took 0.044 seconds. The
    observations were captured at `427cecb3`, before any Rust identity or
    history code existed. OCCT's own prism history plus First/Last shapes
    reached every output entity.
  * Independent fixtures: 114,890 relations over the same 546 cases
    (constructions, every transform step and their compositions), each
    history also passing the independent checker.
  * Derived DRAW cases: `prism_history_rectangle` and
    `prism_history_reversed_triangle` pass on the Rust adapter and on native
    DRAW (OCCT 7.9.3 locally, the Ubuntu 24.04 package in CI). They are
    derived, not original upstream tests.
  * Clean local 600-second campaigns at `34efd36c` (AddressSanitizer,
    standard 20-second/2 GiB limits): `history` completed 600.05 seconds of
    mutation after 26.17 seconds of replay, with 22,614 executions, 6,968
    coverage edges and a 794 MB RSS peak; `identity` completed 600.09 seconds
    after 106.83 seconds of replay, with 6,726 executions, 6,592 edges and an
    848 MB peak. There were no crash, timeout, OOM, slow-unit or disagreement
    artifacts. (An earlier `identity` campaign at the M1 commit `bc8f7015`
    also passed, with 14,350 executions.)

  What remains: histories cover only constructions and rigid motions, so
  `Split`, `Merged` and `Deleted` have been exercised only on hand-written and
  fuzz-synthesized histories, not produced by an operation (M3). Attribute
  storage and outcomes on real operations are M4, and enclosures M5. Ids are
  not yet persisted outside the process (no native format).
* T1 — pending; recorded in `TOPOLOGY_MODEL.md`
* M3 — pending
* M4 — pending
* M5 — pending

## Native format and Parasolid XT

**Decision: the kernel's native format is its own, and XT is an interchange
format, import first.**

The native format must carry everything the contracts above create and
nothing XT can hold: exact rational geometry where the kernel has it, value
ids with their derivations, histories, attributes with outcomes, and
enclosures. XT stores binary64 geometry and Parasolid's tolerant-entity
semantics, so writing our bodies to XT is lossy by construction and reading
XT gives upper bounds rather than enclosures. A format that loses the
kernel's own guarantees cannot be its native format.

XT is still the most valuable interchange target for a Parasolid rival: NX,
Solid Edge, SolidWorks and Onshape exchange it losslessly among themselves,
and its schema is the cheapest public description of Parasolid's data model
(regions, shells, faces, loops, fins, edges, vertices, tolerant entities,
attributes with definitions, per-node identifiers). Siemens publishes the
Parasolid XT Format Reference, and the XT B-rep segment of the ISO 14306 JT
format refers to it. Confirm the document's terms of use before implementing
against it, and do not copy its text into this repository.

Feasibility with OCCT as the oracle:

* **OCCT cannot read or write XT.** There is no differential oracle for an XT
  writer, and an unverified writer is worse than none: a file that NX or
  SolidWorks rejects fails in ways we cannot observe. Do not ship XT export
  until a Parasolid-based validator is part of the evidence. Onshape imports
  Parasolid v10.0–v38.1, is built on Parasolid, and exports both x_t and STEP;
  it is the practical validator for export and the source of twin corpora.
* **An XT reader is checkable.** Read the x_t twin, read the STEP twin
  through OCCT (and later through our own STEP reader), compare face and edge
  counts, surface and curve families, mass properties within enclosures, and
  sampled points; run every imported body through `Topology::check`. Accept a
  subset of node types (points, lines, circles, ellipses, B-curves, planes,
  cylinders, cones, spheres, tori, B-surfaces, the topology nodes, attributes)
  and classify every other node as `Unsupported` with its name, never as an
  approximation. Intersection curves stored as B-curve charts import with the
  chart's tolerance as an `Imported` enclosure.
* **Ordering.** The XT reader belongs after T1 and after the OCCT `.brep`
  converter of T2 (`TOPOLOGY_MODEL.md`), because both readers need the same
  generic topology builder and the `.brep` reader has a native oracle for
  every file.

## Delivery rules for the implementing agent

* Capture native observations before writing kernel code, and check them in
  under `rust/fixtures/occt-<capability>-preimplementation/`.
* Commit subject-only, no AI attribution footers, no absolute paths; store
  paths relative to `$REPO` or `$SDK`. Ask before pushing.
* The Linux CI runner is about 2.6× slower than the development Mac; keep
  every AddressSanitizer fuzz seed under 20 s ÷ 2.6 locally.
* Python fixture generators must round trigonometry and norms once through
  mpmath (`cos_rn`, `sin_rn`, `atan2_rn`, `hypot_rn` in `brep_reference.py`);
  platform libm differs in the last bit and fixture bytes must be identical
  on every host.
* Never count a stripped upstream test as an unchanged pass, never widen a
  tolerance to pass a case, never accept a timeout as a review.
* Update this file's [Acceptance](#acceptance) section and the roadmap rows in
  `PORTING.md` and `PRODUCTION_READINESS.md` in the same commit that claims a
  milestone.
* Never change an existing operation's behaviour without a new algorithm
  level (H8); the old level stays callable and tested.
