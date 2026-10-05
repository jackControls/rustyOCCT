# Handoff: the Rust kernel's S9 Booleans

State of the work for whoever continues it. `REVIEW_NOTES.md` remains the
plan of record: its decisions, evidence, implemented, survey and campaign
bullets are authoritative, and this page only summarizes them and says
where to pick up. Written 2026-10-01, brought up to date 2026-10-03 for
S9f.2a, S9f.2b.1 and S9e.3b (branch `s9e3b`, `s9c2-kernel` merged in),
their DRAW survey (branch `s9-draw-3`), then for S9e.4a (branch `s9e4`,
over `s9c2-kernel` at `e80e2fd9`) and S9f.2b.2 (branch `s9f2b2`, over
`s9c2-kernel` at `f3362b50`), then for S9e.4b.1 (branch `s9e4b`, over
`s9c2-kernel` at `7199e06a`), then for the DRAW survey of S9f.2b.2 and
S9e.4b.1 (branch `s9-draw-4`, over `s9c2-kernel` at `93e6fcd0`),
then for S9f.3a (branch `s9f3`, over `s9c2-kernel` at `7199e06a`) and a
fix of two cylinders within rounding of parallel (branch `fix-conic-eval`),
then for the near-parallel audit (branch `near-parallel-audit`), S9e.4b.2
(branch `s9e4b2`) and S9f.3b (branch `s9f3b`), each over `s9c2-kernel` at
`b903f3fa`, then for the DRAW survey of S9e.4b.2, S9f.3a and S9f.3b (branch
`s9-draw-5`) and the loops' certified integrals (branch `loop-integrals`),
each over `s9c2-kernel` at `12b6c176`, and a Linux fix of prism walls'
frames (branch `linux-imported-arcs`), then for S9e.4b.3a (branch
`s9e4b3`, over `s9c2-kernel` at `12b6c176`), the validator's curve points
guarded (branch `validate-curve-at`) and the arithmetic of fields of
degree eight (branch `given-met-speed`), both over `s9c2-kernel` at
`abebe31c`,
then for S9e.4b.3b (branch `s9e4b3b`, over `s9c2-kernel` at `b99799ec`),
its pole lift combined with the meridian loop's fix (branch `fix-uv-gap`),
then for the DRAW survey of S9e.4b.3a, S9e.4b.3b and the fuzz fixes
(branch `s9-draw-6`, over `s9c2-kernel` at `d665df29`),
then for S9e.4b.3c.1 (branch `s9e4b3c`, over `s9c2-kernel` at `d665df29`),
then for the DRAW survey of S9e.4b.3c.1 and the switches' speed-up (branch
`s9-draw-7`, over `s9c2-kernel` at `f4b584f7`),
then for S9e.4b.3c.2 (branch `s9e4b3c2`, over `s9c2-kernel` at `f4b584f7`),
then for S9e.4b.3c.3a (branch `s9e4b3c3`, over `s9c2-kernel` at `8a55a3e6`).

## Where things stand

- **S1–S8 done.** S1–S4 accepted, S5–S8 implemented with their evidence and
  campaigns (`REVIEW_NOTES.md`, sections S5–S8).
- **S9a–S9d done.** Booleans of prisms in one frame (S9a), polyhedra in any
  position (S9b), arcs and cylinders in any position (S9c), spheres, cones
  and tori (S9d.1–S9d.4c). Every sub-step has its refined decisions, an
  independent reference, fixtures, a native capture taken before its kernel
  code, the kernel, a DRAW survey and a clean 600 s campaign.
- **S9e.1, S9e.2 and S9e.3a done.** A Boolean's result given to another
  Boolean: results of prisms with lines, arcs and circles (S9e.1,
  `curved/given.rs`), stacks with arc walls, polyhedral results given with
  arcs, and one solid of a result of several (S9e.2, `curved/matched.rs`),
  results of spheres, cones and tori, deeper chains and given results
  against spheres, cones and tori (S9e.3a, `curved/chain.rs`), each with its
  DRAW survey (the last full survey, that of S9e.4b.3c.1 and the
  switches' speed-up at `f4b584f7`, branch `s9-draw-7`: the Boolean
  group's 987 cases and 41 restore cases registered, none failing or
  timing out, native DRAW's audited values those of `d665df29`'s bit for
  bit, the Rust adapter's too but 4 centre coordinates of two restore
  cases moved within rounding by the prism caps' frames) and a clean campaign (S9e.3a's and S9f.1's at `b2765f20`, 989
  runs).
- **S9f.1 done.** Spline prisms against polyhedral prisms in any position
  (`curved/spline_walls.rs`), with its survey and campaign as above.
- **S9f.2a implemented** (on `s9c2-kernel`): spline walls against arc,
  circle and spline walls of a prism whose axis is exactly parallel (turned
  about the axis, rotated about a tilted axis, offset by an amount that
  rounds), generatrices over the profiles' exact crossings
  (`curved/spline_parallel.rs`): decisions ("S9f.2a refined"), 43 cases
  referenced and captured before the kernel code, the kernel within the
  reference on all 39 results and refusing the 4 degenerate tangencies; the
  validator's mass integrals take a steep line by its box. Campaign clean
  at `089c5fbe` (860 runs, the slowest input 44 s under AddressSanitizer on
  a loaded host); its DRAW survey done (that of S9e.3b, S9f.2a and
  S9f.2b.1: no case of the group reaches it).
- **S9f.2b.1 implemented** (on `s9c2-kernel`): spline walls against
  cylinder walls of a prism whose axis crosses theirs (leaning, tilted or
  side rods, rings and stadiums), the meeting's branches as graphs over the
  spline's parameter between turning points (`curved/spline_crossing.rs`)
  and a new procedural edge, `Curve3::WallMeet`, certified per knot span
  (`topology/validate/wall_meet.rs`): decisions ("S9f.2b refined"), 34
  cases referenced and captured before the kernel code, the kernel within
  the reference on all 24 S9f.2b.1 results, refusing the 4 degenerate
  cases and S9f.2b.2's 6 (loops round a turning point inside the faces,
  tower fields). Campaign clean at `6c77655a` (902 runs, the slowest input
  50 s under AddressSanitizer at load 6 to 8); its DRAW survey done (as
  S9f.2a's: no case of the group reaches it).
- **S9f.2b.2 implemented** (branch `s9f2b2`, merged into `s9c2-kernel`): a spline wall's meeting with
  a crossing cylinder turning back inside both faces, a graph over the
  spline prism's height about each such turning point (the run parameter
  the one root of `F(., w)` in a rational window of one arc, verified
  exactly; `Curve3::WallMeet` with a `window`), switched at rational
  parameters; a cap circle on a wall in a plane holding its axis, its
  points in one field of the circle's half-angle tangent instead of a tower
  (`curved/spline_crossing.rs`'s `meeting`, `height_piece`,
  `tower_points`). Decisions ("S9f.2b.2 refined"), 17 more cases referenced
  and the whole set captured again before the kernel code, the kernel
  within the reference on all 45 results of S9f.2b and refusing the 6
  degenerate cases. Campaign clean at `7199e06a` (with S9e.4a's); its DRAW
  survey done (that of S9f.2b.2 and S9e.4b.1: no case reaches its loops or
  towers).
- **S9f.3a implemented** (branch `s9f3`, merged into `s9c2-kernel`):
  spline prisms against spheres, caps and zones in any position, either
  the object (`curved/spline_sphere.rs`): S9f.2b's meeting along the wall's
  rulings over the sphere's three rows (`A = n . n` for every ruling), its
  branches over the run, loops' graphs over the height and switches
  unchanged; a sphere's circle (a rim, the split's great circle) on a wall
  along its crease or, in a plane holding the wall's axis, by a primitive
  element from the arc's implicit equation reduced by the circle's;
  `Curve3::WallMeet`'s `other_sphere` and `wall_meet.rs`'s spans over the
  partner's rows. Decisions ("S9f.3 refined": spheres S9f.3a, cones
  S9f.3b), 33 cases referenced by an independent reference of its own
  (`spline_sphere_boolean_reference.py`) and captured before the kernel
  code, the kernel within the reference on all 27 results and refusing the
  6 degenerate cases. Its fuzz switch `SPLINE_SPHERE` is on since the
  loops' certified integrals were sped up (the corpus's `d7599dbe` 472 s
  under AddressSanitizer before, 24 to 38 s on a host at load 8 to 22).
  Campaign clean at `788f8861` with the switch off and at `cfeab65d`
  with it on; its DRAW survey done (that of S9e.4b.2, S9f.3a and S9f.3b:
  no case reaches it, the group's spline solids boxes by `nurbsconvert`).
- **S9f.3b implemented** (branch `s9f3b`, merged into `s9c2-kernel`): spline prisms
  against cones and frustums in any position, either the object
  (`curved/spline_cone.rs`): S9f.2b's meeting with the cone's radius row of
  negative sign (`spline_crossing`'s rows signed, one code for cylinders,
  spheres and cones); `A = q_u^2 + q_v^2 - k^2 q_w^2` one constant per
  pair, `A > 0` S9f.3a's loops and branches, `A < 0` both nappes over the
  whole run (the other nappe's branch outside the cone's face), `A = 0`
  and the apex (real or virtual) on a spline wall `Degenerate`; a rim on a
  wall along its crease in its own elliptic cylinder or in S9f.2b.2's
  half-angle chart; `Curve3::WallMeet`'s `other_half_angle` with
  `wall_meet.rs`'s polynomials kept per power of the stored half angle's
  tangent. Decisions ("S9f.3b refined"), 27 cases referenced by an
  extension of S9f.3a's reference (`spline_cone_boolean_reference.py`) and
  captured before the kernel code, the kernel within the reference on all
  21 results and refusing the 6 degenerate cases. Its fuzz switch
  `SPLINE_CONE` is on since the loops' certified integrals were sped up
  (the slowest variant 97 s under AddressSanitizer before, 23 to 33 s on a
  host at load 8 to 22). Campaign clean at `788f8861` with the switch off
  and at `cfeab65d` with it on; its DRAW survey done (as S9f.3a's: no
  case reaches it).
- **S9e.3b implemented** (branch `s9e3b`, over `507b8054`, with
  `s9c2-kernel` at `c8e37abe` merged in; not yet pushed): a given result's
  meeting of two curved faces (`Meet`, `Rise`, `Toric`) or a cone's or
  torus's general section met by the partner (`curved/triple.rs`), its
  decisions, evidence and capture committed before the kernel; the DRAW
  survey of the chained cases done (none reaches it), and the full DRAW
  survey (that of S9e.3b, S9f.2a and S9f.2b.1: no status or refusal
  changes, the volume audit's values bit for bit). Campaign clean at
  `b0b9adc6`.
  A given `WallMeet` edge never reaches it (a given result with spline
  walls is refused before, S9f). Its fuzz switch `GIVEN_MET` is on since
  the arithmetic of fields of degree eight was sped up (branch
  `given-met-speed`; REVIEW_NOTES.md's "The degree-eight arrangement's
  arithmetic"), and a given meeting's points on two rulings at equal
  heights are separated by a sheared projection there.
- **S9e.4a implemented** (branch `s9e4`, pushed at `b032c1be`): imported solids
  (a `.brep` or STEP body without a construction, `Solid::imported_with`)
  whose stored topology is a prism of lines, arcs and circles, a sphere,
  cap or zone, a cone or frustum or a whole torus, decided on that
  construction read off their stored surfaces and matched to the stored
  topology (`solid/imported.rs`), standing in a Boolean for the body with
  the history over its stored ids: decisions ("S9e.4 refined"), 69 cases
  on 13 bodies OCCT wrote (`rust/fixtures/imported/`) referenced and
  captured before the kernel, the kernel within the reference on all 57
  solid cases, the 9 degenerate refused, the 3 turned profiles S9e.4b's;
  the DRAW adapter's restored solids reach it: of the 1,814 cases restoring
  a shape for a Boolean (never surveyed before), the import reaches 171
  and 16 evaluate on both backends with native DRAW's volumes, registered
  (1,089 cases; the ledger's `checknbshapes` 2 to 4 mapped and verified).
  Campaign clean at `7199e06a` (845 runs, the slowest 47 s under ASan).
- **S9e.4b.1 implemented** (branch `s9e4b`, merged into `s9c2-kernel`): imported prisms
  whose arcs' ends round off their circles in their caps' frames (every
  CTO-like part in a turned frame, `dee_turn`), each end taken onto its
  circle in the exact model (`curved/snapped.rs`: the circle's rational
  point at the end's half-angle tangent rounded once; a joint of a line
  and an arc or of two arcs of one circle shared): decisions ("S9e.4b
  refined", with the sub-steps S9e.4b.1 to S9e.4b.4), 36 cases on 5 bodies
  OCCT wrote in turned frames referenced and captured before the kernel,
  the kernel within the reference on all 27 solid cases, the 6 degenerate
  refused, the lens (two circles at a joint) S9e.4b.4's; S9e.4a's
  `dee_turn` cases now degenerate under S9's rules (a corner on the box's).
  Its DRAW survey (that of S9f.2b.2 and S9e.4b.1, branch `s9-draw-4`): of
  the 31 such restore cases 7 evaluate on both backends with every check
  and native DRAW's volumes (`bcut_complex/H3`, `K8`, `bfuse_complex/C9`,
  `E9`, `I6`, `N1`, `N9`), registered (1,096 cases), 18 refused by S9's
  rules, 6 S9e.4b.4's; nothing else moves. Campaign clean at `788f8861`.
- **S9e.4b.2 implemented** (branch `s9e4b2`, merged into `s9c2-kernel`): imported
  polyhedra other than prisms (pyramids, frustums of them, slanted wedges,
  results of boxes in different frames) decided on their stored vertices,
  not their planes' common points (a vertex of four planes has none once
  rounded; two files sharing a face, the survey's `buc60803a` and `b`, meet
  at their planes' points 5.3e-15 to 7.2e-15 apart where their stored
  points agree): S9b.2's stored model, each face's polygon of stored
  vertices in exactly planar triangles, an exactly coplanar face joined
  with a partner's on its plane, membership by parity, the history over
  the stored ids (`solid/boolean/polyhedra/imported.rs`). Decisions
  ("S9e.4b.2 refined"), 48 cases on 15 bodies OCCT wrote (wedges, sewn
  polyhedra, results of boxes, the survey's shapes from exact points)
  referenced by an exact reference of convex cells
  (`imported_polyhedra_boolean_reference.py`) and captured before the
  kernel (47 matching, the flush fuse reviewed), the kernel within the
  reference on all 44 solid and empty cases, the flush fuse refused, the
  cavity S9e.4b.4's; against curved faces S9e.4b.4's. The fuzz target's
  `IMPORTED` stage imports the chained stage's planar first result too.
  A trial of the survey's 7 restore cases: 4 evaluate on both backends
  with every check (`buc60803`, `bug102_1`, `bug102_2`, `bopfuse_complex/K5`,
  the last since an imported prism of lines against an imported
  polyhedron is decided on its stored vertices too), `bug578_1` and `_2`
  refused (the frustums' bases 6.6e-7 to 2.0e-6 apart: their fuse two
  solids) and `bfuse_complex/D9` (a shared corner stored 1e-13 apart).
  Campaign clean at `788f8861`. Its DRAW survey (that of S9e.4b.2, S9f.3a
  and S9f.3b, branch `s9-draw-5`) confirms the trial: the 4 evaluate with
  native DRAW's volumes and are registered (1,100 cases), nothing else
  moves.
- **S9e.4b.3a implemented** (branch `s9e4b3`, merged into `s9c2-kernel`): imported plane
  pieces, a body of one sphere, cylinder or cone face and plane faces that
  is none of S9e.4a's constructions, decided as its primitive common its
  planes' half-spaces: S9e.1's given model of that Boolean against a hull
  leaf model of the stored planes, matched to the stored topology on
  import, the history over the stored ids (`curved/pieces.rs`). Decisions
  ("S9e.4b.3 refined", with the split S9e.4b.3a to S9e.4b.3c), 45 cases on
  8 sphere pieces OCCT wrote in turned rational frames
  (`generate_imported_pieces_boolean_fixtures.py`, S9e.3a's chained
  reference with each piece's closed form) captured before the kernel (35
  matching, 10 reviewed), the kernel within the reference on all 33 solid
  cases, the 6 degenerate refused, the 6 one-sphere and bitten cases
  S9e.4b.3c's. Only spheres' pieces come from `.brep` files (a cylinder's
  or cone's oblique section is an ellipse record the reader does not
  read; a sphere's section off its frame's meridians and parallels an
  uncertified pcurve); the kernel's own split pieces are tested as pieces.
  The fuzz target's `IMPORTED` stage imports such first results too. A
  trial of the survey's 39 S9e.4b.3 restore cases: none evaluates (18 a
  piece other than its primitive common its planes, 20 two pieces of one
  sphere, both S9e.4b.3c's; `buc60926` a plane through a cone's apex).
  Its DRAW survey (that of S9e.4b.3a, S9e.4b.3b and the fuzz fixes,
  branch `s9-draw-6`) confirms the trial, native DRAW evaluating all 39;
  two CTO parts reach the recognition too, `bcut_complex/G4` refused as a
  piece other than its primitive common its planes and `bcut_complex/I6`
  as "a tangency between the inputs", raised in its notched tool's own
  first arrangement (the reason names the wrong pair); nothing registered.
- **S9e.4b.3b implemented** (branch `s9e4b3b`): the kernel's own split
  pieces (S8's `Clipped` and `Half`) against curved faces, on S9e.4b.3a's
  model read off the split: its primitive (the prism, the cone, the whole
  torus; a zone's or cap's whole sphere and its ends' parallels' planes)
  common the half-space of the plane it was built on, taken exactly into
  the world through the frame's axes, matched to the piece's stored
  topology, the history over its ids (`curved/splits.rs`). Decisions
  ("S9e.4b.3b refined"), 48 cases on 8 split pieces
  (`generate_split_pieces_boolean_fixtures.py`, a `split` row in both
  protocols, the chained reference with each piece its solid common a
  half-space box) captured before the kernel (33 matching, 15 reviewed),
  the kernel within the reference on all 39 solid cases, the 6 degenerate
  refused (a frustum's half through its axis at its apex, a box on a cut
  plane), the 3 one-sphere cases S9e.4b.3c's. A cone's half through its
  axis stays `Degenerate`, a torus's spiric piece and a spline prism's
  piece refused. The `split` fuzz target's new stage (`PIECE_BOOLEANS`)
  found three assembly rules, each an imported piece's too, now fixed: a
  loop through a stored pole lifted from the pole on (its holes with it),
  a stored circle and a cone's stored section matched as rings run the
  result's way. Its campaigns found two more, both fixed: a torus cap's
  round end lifted the wrong way (S8d.3's, latent) and a plane within
  rounding of a frustum's virtual apex reaching the arrangement (now
  `Degenerate` before its edges enter it). At `cc7ea7ff`, with those
  fixes and `CONE_PAIRS` on, both campaigns are clean (`boolean` 996
  runs, the slowest input 13 s; `split` 1,731 runs, the slowest 19 s).
  Its DRAW survey (with S9e.4b.3a's and the fuzz fixes', branch
  `s9-draw-6` at `d665df29`) moves no case.
- **S9e.4b.3c.1 implemented** (branch `s9e4b3c`): two inputs on one sphere
  in general position and sphere pieces whose rims OCCT divided (the
  survey's `so1` and `so4`): faces on one sphere as S9c.1's faces on one
  surface, the circles of both crossing there met exactly (a whole
  sphere's split great circle too), a piece's sphere split at the second
  arrangement's seam, a stored vertex on a piece's ring splitting it in its
  own arrangement (`curved/graph.rs`, `curved/pieces.rs`); a body whose one
  curved face's material lies outside its quadric (a groove) refused on
  import before its model. Decisions ("S9e.4b.3c refined", with the split
  S9e.4b.3c.1 to S9e.4b.3c.3), 22 cases on 5 pieces OCCT wrote
  (`generate_one_sphere_boolean_fixtures.py`, a `write` block's `divide`
  row) captured before the kernel (8 matching, 14 reviewed), the kernel
  within the reference on all 17 solid and empty cases, the 2 degenerate
  refused, the 3 exact incidences S9e.4b.3c.2's; S9e.4b.3a's and S9e.4b.3b's
  one-sphere cases solid. `PIECE_BOOLEANS` gives a zone's or cap's piece to
  a ball of its sphere too. A trial of the survey's 38 `so` cases: 14
  evaluate on both backends (`so1` and `so4`, `so4` and `so2`, `so2` and
  `so6`'s common), 13 S9e.4b.3c.2's, 11 `Degenerate` (`so6` and `so7`'s
  corners within rounding of the partner's); `bcut_complex/I6` and `G4`
  S9e.4b.3c's. Its DRAW survey (with the switches' speed-up, branch
  `s9-draw-7` at `f4b584f7`) confirms the trial: the 14 evaluate with
  native DRAW's volumes and are registered (1,114 cases), the other 24
  and `I6` are refused as it found; no other status or reason moves but
  two restores ending within the survey's 120 s by the host's load.
- **S9e.4b.3c.2 implemented** (branch `s9e4b3c2`): exact incidences of two
  inputs on one sphere (the survey's `so1` and `so2`, `so2` and `so3`,
  `so5` and `so2`): model vertices of both at one exact point one vertex
  (`VKey::Both`), a vertex of one inside a line or circle edge of the other
  splitting it, a part of B's edge alike a part of A's one arrangement edge
  in both inputs' faces (`Arr::shared`), sections and edges along it taken
  by it, plane faces on one plane as faces on one surface, the splits'
  great circles alike (`curved/graph.rs`). Decisions ("S9e.4b.3c.2
  refined"), 45 cases on 11 pieces OCCT wrote on the world's axes (in
  turned frames an exact incidence along a line holds only between planes
  stored bit for bit alike) with an exact reference of their own (the ball
  cut into cells by the inputs' heights and half-planes about the common
  axis, `generate_one_sphere_incidence_boolean_fixtures.py`) captured
  before the kernel (42 matching, 3 reviewed), the kernel within the
  reference on all 33 solid and empty cases, the 12 degenerate refused by
  S9's rules; S9e.4b.3c.1's `hemi_octant` solid. A panic on importing a
  piece whose stored pole lies exactly above its rim's centre (S9e.4b.3c.1's
  `Arr::split_at`) fixed before the capture. `PIECE_BOOLEANS` gives a
  zone's or cap's piece to the solid it was split from too. A trial of the
  survey's restore cases: the 13 evaluate on both backends with every
  check (27 of the 38 `so` cases in all; `so6` and `so7`'s 11 `Degenerate`
  as before).
- **S9e.4b.3c.3a implemented** (branch `s9e4b3c3`): an imported body of one
  sphere, cylinder or cone face and plane faces that is another Boolean of
  its primitive and the convex hull of its other planes than a common
  (`solid/imported.rs`'s `Form`, `curved/pieces.rs`): the hull less the
  primitive (a hole, groove, slot or dimple: its curved face's material
  outside its quadric), the primitive less the hull of its planes turned
  over (a bite) and their fuse (a boss: the survey's `bcut_complex/G4`
  part), tried in turn with S9e.2's match the arbiter; a cylinder's or
  cone's primitive over its curved face's range ending at its caps, two
  faces on one plane facing one way one plane of the hull, OCCT's vertex
  loops at a sphere's poles left unmatched, a tangency in the piece's own
  arrangement named for the body (`I6`'s tool). Decisions ("S9e.4b.3c.3
  refined", with the split S9e.4b.3c.3a and S9e.4b.3c.3b), 39 cases on 9
  bodies OCCT wrote (`generate_piece_forms_boolean_fixtures.py`, the
  chained reference with each body's closed form) captured before the
  kernel (26 matching, 13 reviewed), the kernel within the reference on
  all 33 solid cases, the 6 degenerate refused (the notch's tangency, the
  three-quarter frustum's apex); S9e.4b.3a's `bitten_box` solid. A trial
  of the survey's restore cases: `G4` evaluates on both backends with every
  check, `I6` refused for its tool's own tangency (its part tangent to the
  tool besides), `buc60926` at its frustum's apex as before, the 38 `so`
  cases as before.
- **CI.** Both workflows ("Rust kernel", "Rust geometry fuzzing") were green
  at `6c221525`. They had been red from S7 until 2026-09-29, unnoticed;
  check them after every push (see "Working rules").

## What is open, in order

1. **`GIVEN_MET` is on** (branch `given-met-speed`): the arithmetic of
   fields of degree eight was sped up (REVIEW_NOTES.md's "The degree-eight
   arrangement's arithmetic"), the corpus's slowest inputs too ("The
   boolean target's slowest inputs": `e36969f1` from 67 to 20 s under
   AddressSanitizer), and a meridian loop through one pole of a sphere
   fixed (S9e.4b.3a's amendment, found by the replay with the switch on).
   The campaign with `GIVEN_MET`, `SPLINE_SPHERE` and `SPLINE_CONE` on is
   clean at `a32d256f` (895 runs, the slowest input 31 s under
   AddressSanitizer). Open levers: libFuzzer's leak-check rerun (an input
   whose run keeps the kernel's caches runs twice under the sanitizer),
   the torus meetings' root sampling (the certified integrals along a cone
   carrier's meetings are narrowed since the off switches' track:
   `487e8cac` 39 to 20 G instructions, 16.5 to 9.7 s under the sanitizer).
2. **The Linux records and reviews** of the recaptured
   `occt-boolean-spline-crossing-preimplementation` and of the
   spline-sphere, spline-cone, imported-polyhedra, imported-arcs,
   imported-pieces, split-pieces, one-sphere, one-sphere-incidence and
   piece-forms captures from CI's run, as every capture's. `SPLINE_SPHERE` and
   `SPLINE_CONE` are on since the loops' certified integrals were sped up
   (branch `loop-integrals`; REVIEW_NOTES.md's "Certified integrals beside
   spline walls' loops"), and the campaign with both on is clean at
   `cfeab65d` (803 runs, the slowest input 52 s under AddressSanitizer, an
   existing corpus input reaching neither). The DRAW survey of S9e.4b.2,
   S9f.3a and S9f.3b is done at `12b6c176` (branch `s9-draw-5`: S9e.4b.2's
   4 restore cases registered, no case reaching the spline walls), that of
   S9e.4b.3a, S9e.4b.3b and the fuzz fixes at `d665df29` (branch
   `s9-draw-6`: no status moving, 41 restore cases' reasons, no case
   registered), that of S9e.4b.3c.1 and the switches' speed-up at
   `f4b584f7` (branch `s9-draw-7`: the 14 `so` cases registered).
3. **S9e.4b**, split in "S9e.4b refined": S9e.4b.1 (arcs rounded off
   their circles) and S9e.4b.2 (polyhedra other than prisms, on their
   stored vertices) are implemented with their DRAW surveys and
   campaigns; **S9e.4b.3**, split in "S9e.4b.3 refined": S9e.4b.3a (plane
   pieces as their primitive common their planes' half-spaces) and
   S9e.4b.3b (S9e.2's deferred `Clipped` and `Half` against curved faces
   on that model) are implemented with clean campaigns and their DRAW
   survey (`s9-draw-6` at `d665df29`: the 39 refused as the trial found,
   nothing registered; pending their captures' Linux records);
   **S9e.4b.3c**, split in "S9e.4b.3c refined": S9e.4b.3c.1 (two inputs on
   one sphere in general position, rims split by stored vertices: 14 of the
   survey's 38 `so` cases evaluate; a groove such as `bcut_complex/I6`'s
   tool or `G4`'s boss refused on import as S9e.4b.3c's) is implemented
   with its DRAW survey (`s9-draw-7` at `f4b584f7`: the 14 registered,
   1,114 cases) and a clean campaign at `59d0c57b` (pending its capture's
   Linux record); **S9e.4b.3c.2**, exact incidences on one sphere (a
   vertex of both, a circle of both, a line of both, plane faces on one
   plane with overlapping edges: `so1` and `so2`, `so2` and `so3`, `so5`
   and `so2`, 13 cases, all evaluating in its trial) is implemented
   (campaigns clean at `8a55a3e6`; pending its capture's Linux record and
   the DRAW survey); **S9e.4b.3c.3**, split in "S9e.4b.3c.3 refined":
   **S9e.4b.3c.3a** (a body of one curved face and planes as one Boolean of
   its primitive and the hull of its other planes: grooves, bites, bosses,
   two faces on one plane; `bcut_complex/G4` evaluating in its trial) is
   implemented (campaigns clean at `e8940c22`; pending its capture's Linux
   record and the DRAW survey); next **S9e.4b.3c.3b**, bodies of one curved face needing
   more than one hull (a groove in a body not convex in its planes, a
   primitive bitten twice, a sphere's ends with another form: the fuzz
   target's 4 remaining refusals); then
   **S9e.4b.4**, the S9e text's plan in full
   (joints of two circles, prisms with walls of two directions, bodies of
   several curved surfaces: 21 cases; an imported polyhedron against
   curved faces, or with a cavity). The reader's header check refuses
   OCCT 8.1's version-3 `.brep` (`(c) Open Cascade`; 27 dataset files,
   none among the surveyed restore cases): widening it needs a survey of
   the restores it opens.
4. **S9's acceptance** (U6): kernel and fuzz CI green at the accepted
   revision, the schedule run's full replays green (boolean and
   `degree_elevation` are sharded across four jobs plus a completeness check,
   `REPLAY_SHARDS`), and a clean local 600 s campaign; record it in
   `BOOLEAN.md` and mark S9 done in `REVIEW_NOTES.md`.

Refused by design and staying refused (each documented): a tangency between
the inputs, a cavity beside several solids, spline segments along one curve
of different forms, an arc ending off its circle, a section through a
sphere's pole off its meridians, a torus's tube circle on the other surface,
a result touching itself, splines against tori, spline walls against spline
walls with crossing axes, coincident spline walls in different frames, a
torus among three curved surfaces or two tori in one plane at a given
meeting (S9e.3b, by cost).

Parallel tracks (`REVIEW_NOTES.md`, "Parallel tracks"): the degree-eight
arrangement's arithmetic is done (`GIVEN_MET` on); its open levers (the
leak-check rerun, the torus meetings' root sampling); the boolean target's
off switches are measured (branch `fuzz-switches`: `CONE_PAIRS` on and
its campaign clean at `cc7ea7ff`, its
three kinds of failure fixed and two quadrics' meetings' integrals
narrowed, among them a cone carrier's); the last two off switches are on
(branch `torus-turned-speed`: `TORUS_PAIRS` and `TURNED_PARTS`, their
slowest variants 22 s a run under the sanitizer, 30 and 37 to 53 s before
side by side; the campaigns with every switch on are clean at `59d0c57b`:
`boolean` 1,080 runs, the slowest input 13 s; `split` 1,818 runs); the
fuzz switches below.

## Open user decisions

- **U9**: tessellation certifies every triangle's deflection; keep that as
  the only mode, or add a display mode with sampled control and no bound.
  Until answered, only the certified mode exists.
- **U10**: R4's homogeneous C1 test refuses the usual rational NURBS circle
  (STEP-b's `rational_cylinder`): (a) keep refusing, (b) test C1 of the
  rational curve exactly at such knots (recommended), (c) split at such
  knots on import. Until answered, (a) holds.

## How a step is done

Each sub-step follows the same order, and the history shows it in its
commits:

1. **Refined decisions**, a `REVIEW_NOTES.md` bullet "… refined, before its
   code": why the case is refused today, the representation, the degrees and
   fields, what is `Degenerate`, what stays refused.
2. **Evidence**: an independent reference in `rust/tools/` (mpmath, rounding
   once; never `math.cos`/`sin`/`hypot`, whose results differ between Python
   versions and platforms), with two-way checks, closed forms and margins;
   fixtures from a generator with `--check` (identical under Python 3.9 and
   3.12); a reference unit test; the generator in a CI fixture-generator
   group of `.github/workflows/rust-kernel.yml`.
3. **Native capture** before the kernel file exists: the compare script is
   keyed on that file, the capture's metadata records it absent; reviews in a
   divergences JSON only with independent evidence (an adaptive BRepGProp
   build, a Green's-theorem integration, a DRAW check), fingerprinted by the
   native row.
4. **Kernel**, with a test file: every fixture within `1e-9` enclosures of the
   reference, degenerate cases refused, histories complete, results
   deterministic and moved rigidly.
5. **Verification** (below), then merge and push.
6. **DRAW survey** of the Boolean group, registering cases that now evaluate
   within the 30 s contract with a volume audit; then a **campaign**.

Agents worked each step in its own git worktree, merged into `s9c2-kernel`,
which was verified and pushed to `rust-kernel` (a fast-forward).

## Verification

Run from the repository root, with the pinned SDK built at
`target/spline-linear-preflight/pinned-sdk` and the math venv at
`target/math-oracle-venv`:

```bash
cargo +stable fmt -p rusty-occt --check
cargo +stable clippy -q -p rusty-occt --all-targets --release -- -D warnings
cargo +1.85 check -q -p rusty-occt --all-targets
cargo +stable test -q -p rusty-occt --release --no-fail-fast
```

Each comparison takes `--occt-root target/spline-linear-preflight/pinned-sdk/install --sdk-manifest target/spline-linear-preflight/pinned-sdk/build-manifest.json`:

| Comparison | Matches / reviewed |
|---|---|
| `compare_boolean.py` | 45 / 0 |
| `compare_boolean.py --splines` | 33 / 13 |
| `compare_polyhedral.py` | 43 / 2 |
| `compare_curved_boolean.py` | 42 / 2 |
| `compare_procedural_boolean.py` | 4 / 24 |
| `compare_turned_boolean.py` | 2 / 13 |
| `compare_capped_boolean.py` | 0 / 18 |
| `compare_sphere_boolean.py` | 30 / 0 |
| `compare_spheres_boolean.py` | 12 / 21 |
| `compare_cone_boolean.py` | 25 / 5 |
| `compare_cones_boolean.py` | 21 / 19 |
| `compare_torus_boolean.py` | 11 / 24 |
| `compare_torus_segment_boolean.py` | 15 / 14 |
| `compare_torus_curved_boolean.py` | 15 / 29 |
| `compare_spheres_turned_boolean.py` | 0 / 18 |
| `compare_cones_loops_boolean.py` | 5 / 26 |
| `compare_torus_parts_boolean.py` | 16 / 37 |
| `compare_chained_boolean.py` | 24 / 6 |
| `compare_given_boolean.py` | 36 / 0 |
| `compare_given_curved_boolean.py` | 25 / 23 |
| `compare_spline_any_boolean.py` | 22 / 16 (the kernel within the reference on all 33 solid and empty cases, the 5 degenerate refused) |
| `compare_spline_parallel_boolean.py` | 31 / 12 (the kernel within the reference on all 39 results, the 4 degenerate refused) |
| `compare_spline_crossing_boolean.py` | 14 / 37 (the kernel within the reference on all 45 results, the 6 degenerate refused, none `unsupported`) |
| `compare_given_met_boolean.py` | 8 / 42 (the kernel within the reference on all 42 solid cases, the 8 degenerate refused) |
| `compare_imported_boolean.py` | 54 / 15 (the kernel within the reference on all 57 solid cases, the 12 degenerate refused: since S9e.4b.1 the 3 turned profiles' corner on the box's) |
| `compare_imported_arcs_boolean.py` | 27 / 9 (the kernel within the reference on all 27 solid cases, the 6 degenerate refused, the 3 lens cases `unsupported`, S9e.4b.4's) |
| `compare_spline_sphere_boolean.py` | 4 / 29 (the kernel within the reference on all 27 results, the 6 degenerate refused, none `unsupported`) |
| `compare_imported_polyhedra_boolean.py` | 47 / 1 (the kernel within the reference on all 44 solid and empty cases, the flush fuse refused, the cavity `unsupported`, S9e.4b.4's) |
| `compare_spline_cone_boolean.py` | 4 / 23 (the kernel within the reference on all 21 results, the 6 degenerate refused, none `unsupported`) |
| `compare_imported_pieces_boolean.py` | 30 / 15 (the kernel within the reference on all 39 solid cases, the 6 degenerate refused; the one-sphere cases solid since S9e.4b.3c.1, the bitten ball's since S9e.4b.3c.3a, two of them reviewed for OCCT's seam and pole edges) |
| `compare_split_pieces_boolean.py` | 26 / 22 (the kernel within the reference on all 42 solid and empty cases, the 6 degenerate refused; the one-sphere cases solid since S9e.4b.3c.1) |
| `compare_one_sphere_boolean.py` | 8 / 14 (the kernel within the reference on all 20 solid and empty cases, the 2 degenerate refused; the 3 exact incidences solid since S9e.4b.3c.2) |
| `compare_one_sphere_incidence_boolean.py` | 42 / 3 (the kernel within the reference on all 33 solid and empty cases, the 12 degenerate refused) |
| `compare_piece_forms_boolean.py` | 20 / 19 (the kernel within the reference on all 33 solid cases, the 6 degenerate refused) |

Every one must report 0 failures; since S9f.1 also `compare_split.py`
(72 / 56) and `compare_brep.py --family spline` (10 / 3), which share its
spline prisms, and `compare_brep_io.py` (6,835 / 7); since S9e.4
`compare_step.py` (23 / 6, with STEP-b's SDK at `target/step-b-sdk`: the
pinned SDK lacks the STEP toolkits' headers), whose bodies imported solids
may be. Where `curved_boolean_reference.py`,
`chained_curved_boolean_reference.py` or `torus_curved_boolean_reference.py`
changes, every generator importing it runs with `--check --workers 4`. Then:

- the fuzz crate: `cd rust/fuzz && cargo +nightly-2026-09-22 fmt --check && cargo +nightly-2026-09-22 check`;
- a replay with debug assertions of every boolean corpus input and every
  `rust/fuzz/regressions/boolean*` file, one process per file: build with
  `CARGO_INCREMENTAL=0 RUSTFLAGS="-C debug-assertions" cargo +nightly-2026-09-22 build --release --example replay_boolean` in `rust/fuzz`;
  for a switch that is off (none is now; `TORUS_PAIRS`, `CONE_PAIRS`,
  ... when one is), replay once more with it switched on in the source,
  which reaches the paths it gates;
- `target/math-oracle-venv/bin/python -m unittest discover -s rust/tools -p 'test_*.py'` (the system Python lacks mpmath);
- `python3 rust/tools/run_upstream_tests.py --ledger`;
- a campaign: `python3 rust/tools/run_fuzz.py --target boolean --seconds 600 --toolchain nightly --replay sample` (the full local corpus exceeds the startup hour under AddressSanitizer; the sampled replay plus the debug-assertion replay above stand in for it, and CI's schedule runs the full replay in shards).

## Working rules and lessons

- **Watch CI after every push.** `gh run list -R jackControls/rustyOCCT -L 6`,
  and `--log-failed` on any failure. Recurring causes so far:
  - a new native capture without its Linux record (`platform-linux/`, taken
    from the run's artifact and reviewed against the macOS capture);
  - reviews fingerprinted with macOS native rows only (the Linux rows differ
    in their last bits and need their own reviews);
  - job time limits as generators and fixtures grow;
  - Linux's correctly rounded `hypot` against macOS's: a frame normalized
    again differs by an ulp between hosts, so exact parallel or
    perpendicular tests must tolerate it (reproduce on macOS by turning a
    normal by an ulp in a test).
- **Ubuntu's packaged OCCT is 7.6.3**, older than the pinned 8.1: a case may
  state that release's own outcome (`expected_occt_by_version`).
- **Fuzz timing.** AddressSanitizer costs about twelve times the plain time
  on the exact arithmetic, and Linux runners about 2.6 times this host; the
  boolean target's limit is 60 s an input. Configurations too slow for it are
  switched off in `rust/fuzz/src/boolean.rs`, below the 35 s a run under the
  sanitizer here the switches on keep to; none is off now (`TORUS_PAIRS` and
  `TURNED_PARTS` on since REVIEW_NOTES.md's "The boolean target's last two
  off switches", their slowest variants 22 s; `CONE_PAIRS`, `GIVEN_CURVED`,
  `GIVEN_ROUND`, `GIVEN_BALL`, `GIVEN_MET`, `SPLINE_WALLS`,
  `SPLINE_PARALLEL`, `SPLINE_CROSSING`, `SPLINE_SPHERE`, `SPLINE_CONE` and
  `IMPORTED` too). A heavily loaded host makes campaigns time out
  spuriously; run them on a quiet machine.
- **Keep debug-assertion tests optimized.** CI runs them with
  `CARGO_PROFILE_DEV_OPT_LEVEL=2`; time new test files that way.
- **Commits** are subject-only and long and descriptive, with no AI
  attribution; a staged diff must contain no machine paths. Push only after
  every check passes, and chain the push after the checks.
