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
each over `s9c2-kernel` at `12b6c176`, a Linux fix of prism walls'
frames (branch `linux-imported-arcs`), and the arithmetic of fields of
degree eight (branch `given-met-speed`, over `s9c2-kernel` at
`abebe31c`).

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
  DRAW survey (the last full survey, that of S9e.4b.2, S9f.3a and S9f.3b
  at `12b6c176`, branch `s9-draw-5`: the Boolean group's 987 cases and 27
  restore cases registered, none failing or timing out, the volume audit's
  values those of `93e6fcd0`'s bit for bit) and a clean campaign (S9e.3a's and S9f.1's at `b2765f20`, 989
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
- **CI.** Both workflows ("Rust kernel", "Rust geometry fuzzing") were green
  at `6c221525`. They had been red from S7 until 2026-09-29, unnoticed;
  check them after every push (see "Working rules").

## What is open, in order

1. **A campaign with `GIVEN_MET` on** (branch `given-met-speed`, not yet
   merged or pushed). The switch is on: the arithmetic of fields of degree
   eight was sped up (REVIEW_NOTES.md's "The degree-eight arrangement's
   arithmetic": signs from a certified dyadic point before Sturm-Tarski,
   numbers as integers over one denominator, fewer inverses, rational
   reductions by the crate's Lehmer gcd, coprimality modulo a prime), the
   corpus's slowest chained operations reaching a given meeting 24 to 42 G
   instructions with debug assertions (92 to 123 G before), under
   AddressSanitizer 21 to 35 s at load 6 to 7 (60 to 71 s before), each two
   runs (libFuzzer's leak check); the corpus and its regressions replay
   with debug assertions and the switch on without a failure. The "points
   not separated by a projection" limits were two cylinders with parallel
   axes and a sphere crossing their two common rulings at equal heights; a
   sheared projection separates them now (4 evaluate, 4 degenerate). Next:
   merge, push, CI, and the boolean campaign with the switch on. Open
   levers, in REVIEW_NOTES.md: libFuzzer's leak-check rerun (an input whose
   run keeps the kernel's cached arrangements runs twice under the
   sanitizer), the certified integrals along a cone carrier's meetings, the
   torus meetings' root sampling (`e36969f1`, the corpus's slowest input,
   43.6 G).
2. **The Linux records and reviews** of the recaptured
   `occt-boolean-spline-crossing-preimplementation` and of the
   spline-sphere, spline-cone, imported-polyhedra and imported-arcs
   captures from CI's run, as every capture's. `SPLINE_SPHERE` and
   `SPLINE_CONE` are on since the loops' certified integrals were sped up
   (branch `loop-integrals`; REVIEW_NOTES.md's "Certified integrals beside
   spline walls' loops"), and the campaign with both on is clean at
   `cfeab65d` (803 runs, the slowest input 52 s under AddressSanitizer, an
   existing corpus input reaching neither). The DRAW survey of S9e.4b.2,
   S9f.3a and S9f.3b is done at `12b6c176` (branch `s9-draw-5`: S9e.4b.2's
   4 restore cases registered, no case reaching the spline walls).
3. **S9e.4b**, split in "S9e.4b refined": S9e.4b.1 (arcs rounded off
   their circles) is implemented (above, its DRAW survey and campaign done); **S9e.4b.2**, polyhedra other than prisms, is
   implemented on their stored vertices (above; its DRAW survey and
   campaign done, 4 restore cases registered); next **S9e.4b.3**, a
   plane's pieces of a sphere, cylinder or cone and S9e.2's deferred
   `Clipped` and `Half` against curved faces (39 cases), then
   **S9e.4b.4**, the S9e text's plan in full (joints of two circles, prisms
   with walls of two directions, bodies of several curved surfaces: 21
   cases; an imported polyhedron against curved faces, or with a cavity).
   The reader's header check
   refuses OCCT 8.1's version-3 `.brep` (`(c) Open Cascade`; 27 dataset
   files, none among the surveyed restore cases): widening it needs a
   survey of the restores it opens.
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
leak-check rerun, the integrals along a cone carrier's meetings, the torus
meetings' root sampling); the fuzz switches below.

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
  for a switch that is off (`TORUS_PAIRS`, `CONE_PAIRS`, ...), replay once
  more with it switched on in the source, which reaches the paths it gates;
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
  switched off in `rust/fuzz/src/boolean.rs` (`TORUS_PAIRS`, `CONE_PAIRS`,
  `TURNED_PARTS`; `GIVEN_CURVED`, `GIVEN_ROUND`, `GIVEN_BALL`, `GIVEN_MET`,
  `SPLINE_WALLS`, `SPLINE_PARALLEL`, `SPLINE_CROSSING`, `SPLINE_SPHERE`,
  `SPLINE_CONE` and `IMPORTED` are on; the three off switches and
  `GIVEN_MET` on at once make corpus input `6fab9d41` fail
  `vertex_off_curve`, each alone clean, open). A heavily loaded
  host makes campaigns time out spuriously; run them on a quiet machine.
- **Keep debug-assertion tests optimized.** CI runs them with
  `CARGO_PROFILE_DEV_OPT_LEVEL=2`; time new test files that way.
- **Commits** are subject-only and long and descriptive, with no AI
  attribution; a staged diff must contain no machine paths. Push only after
  every check passes, and chain the push after the checks.
