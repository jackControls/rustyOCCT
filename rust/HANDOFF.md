# Handoff: the Rust kernel's S9 Booleans

State of the work for whoever continues it. `REVIEW_NOTES.md` remains the
plan of record: its decisions, evidence, implemented, survey and campaign
bullets are authoritative, and this page only summarizes them and says
where to pick up. Written 2026-10-01; the code it describes is that of
`6c221525` on `rust-kernel`, with the documentation commits after it.

## Where things stand

- **S1–S8 done.** S1–S4 accepted, S5–S8 implemented with their evidence and
  campaigns (`REVIEW_NOTES.md`, sections S5–S8).
- **S9a–S9d done.** Booleans of prisms in one frame (S9a), polyhedra in any
  position (S9b), arcs and cylinders in any position (S9c), spheres, cones
  and tori (S9d.1–S9d.4c). Every sub-step has its refined decisions, an
  independent reference, fixtures, a native capture taken before its kernel
  code, the kernel, a DRAW survey and a clean 600 s campaign.
- **S9e.3a implemented** (branch `s9e3`, not yet pushed): given results of
  spheres, cones and tori, deeper chains and given results against spheres,
  cones and tori (`curved/chain.rs`), with its evidence captured before the
  kernel; its DRAW chained cases surveyed (`G9` and `H3` refused as a
  tangency); pending the full DRAW survey and the campaign.
- **S9e.1 and S9e.2 done.** A Boolean's result given to another Boolean:
  results of prisms with lines, arcs and circles (S9e.1, `curved/given.rs`),
  stacks with arc walls, polyhedral results given with arcs, and one solid of
  a result of several (S9e.2, `curved/matched.rs`), each with its DRAW
  survey (the last full survey: 987 cases registered, none failing or
  timing out) and a clean campaign (S9e.2's at `6c221525`, 1,230 runs).
- **S9f decided, S9f.1 implemented** (branch `s9f1`, not yet pushed).
  Splines in any position: decisions recorded ("S9f refined"), S9f.1's 38
  cases (spline prisms against polyhedral prisms in any position)
  referenced and captured before any kernel code, and its kernel
  (`curved/spline_walls.rs`): spline walls in the curved engine, R4's knot
  of multiplicity `p` removed exactly before every lift (the extrusion in a
  turned frame failed before), the validator's exact Green path and
  quadrature taking a crease's pieces across a knot line by a sliver and
  past a cap by their control polygon; `boolean-spline-any-r4-*` adds the
  rounded knot and exact parallels after the capture. Pending its DRAW
  survey and campaign.
- **S9f.2a implemented** (branch `s9f2a`, not yet pushed): spline walls
  against arc, circle and spline walls of a prism whose axis is exactly
  parallel (turned about the axis, rotated about a tilted axis, offset by
  an amount that rounds), generatrices over the profiles' exact crossings
  (`curved/spline_parallel.rs`): decisions ("S9f.2a refined"), 43 cases
  referenced and captured before the kernel code, the kernel within the
  reference on all 39 results and refusing the 4 degenerate tangencies; the
  validator's mass integrals take a steep line by its box. Pending its DRAW
  survey and campaign.
- **S9f.2b.1 implemented** (branch `s9f2b`, not yet pushed): spline walls
  against cylinder walls of a prism whose axis crosses theirs (leaning,
  tilted or side rods, rings and stadiums), the meeting's branches as
  graphs over the spline's parameter between turning points
  (`curved/spline_crossing.rs`) and a new procedural edge,
  `Curve3::WallMeet`, certified per knot span
  (`topology/validate/wall_meet.rs`): decisions ("S9f.2b refined"), 34
  cases referenced and captured before the kernel code, the kernel within
  the reference on all 24 S9f.2b.1 results, refusing the 4 degenerate
  cases and S9f.2b.2's 6 (loops round a turning point inside the faces,
  tower fields). Pending its DRAW survey and campaign.
- **CI.** Both workflows ("Rust kernel", "Rust geometry fuzzing") are green
  at `6c221525`. They had been red from S7 until 2026-09-29, unnoticed; check
  them after every push (see "Working rules").

## What is open, in order

1. **S9e.3a's DRAW survey and campaign** (its decisions, evidence, capture
   and kernel done on branch `s9e3`: given results of spheres, cones and
   tori, with procedural edges the partner does not reach, deeper chains to
   three Booleans, given results against spheres, cones and tori; DRAW's
   `G9` and `H3` refused as the tangency they are). Then **S9e.3b**,
   evidence first: a given edge that is a meeting of two curved faces
   (`Meet`, `Rise`, `Toric`, a cone pair's), a cone's section against a
   curved face or a torus's general plane section against a curved face,
   met by the partner (three surfaces, two curved: resultants and certified
   isolation along the procedural curve); 156 of the fuzz corpus's 1,238
   chained operations are refused as S9e.3b's. The old local branch
   `s9e3-wip` is superseded: do not merge it.
2. **S9f.1's DRAW survey and campaign** (its kernel done on branch `s9f1`, over
   `s9c2-kernel` with S9e.3a: spline prisms against prisms of lines in any
   position, `REVIEW_NOTES.md`'s "S9f.1 implemented"). The survey should
   find upstream cases of spline prisms in turned frames newly evaluating;
   the campaign runs the `boolean` target with `SPLINE_WALLS` on and R4's
   knot decoded (the byte at index 14). The old local branch `s9f1-kernel`
   is superseded: do not merge it.
3. **S9e.4**: imported bodies decided on their stored surfaces.
4. **S9f.2a's DRAW survey and campaign** (its kernel done on branch
   `s9f2a`, over `s9c2-kernel`: `REVIEW_NOTES.md`'s "S9f.2a implemented";
   the campaign runs the `boolean` target with `SPLINE_PARALLEL` on and
   the tilted offset decoded, the byte at index 15). Then **S9f.2b.1's
   DRAW survey and campaign** (its kernel done on branch `s9f2b`, over
   `s9c2-kernel`: `REVIEW_NOTES.md`'s "S9f.2b.1 implemented"; the campaign
   runs the `boolean` target with `SPLINE_CROSSING` on), **S9f.2b.2**
   (loops round a turning point inside the faces, each piece about it a
   graph over the cylinder's angle in a verified window, and the tower
   fields of a cap circle on a wall along its axis), **S9f.3** (spheres
   and cones; "S9f refined" gives the degrees).
5. **S9's acceptance** (U6): kernel and fuzz CI green at the accepted
   revision, the schedule run's full replays green (boolean and
   `degree_elevation` are sharded across four jobs plus a completeness check,
   `REPLAY_SHARDS`), and a clean local 600 s campaign; record it in
   `BOOLEAN.md` and mark S9 done in `REVIEW_NOTES.md`.

Refused by design and staying refused (each documented): a tangency between
the inputs, a cavity beside several solids, spline segments along one curve
of different forms, an arc ending off its circle, a section through a
sphere's pole off its meridians, a torus's tube circle on the other surface,
a result touching itself, splines against tori, spline walls against spline
walls with crossing axes, coincident spline walls in different frames.

Parallel tracks (`REVIEW_NOTES.md`, "Parallel tracks"): the degree-eight
arrangement arithmetic (the next lever for speed); the fuzz switches below.

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
| `compare_spline_crossing_boolean.py` | 13 / 21 (the kernel within the reference on all 24 S9f.2b.1 results, the 4 degenerate refused, S9f.2b.2's 6 `unsupported`) |

Every one must report 0 failures; since S9f.1 also `compare_split.py`
(72 / 56) and `compare_brep.py --family spline` (10 / 3), which share its
spline prisms. Then:

- the fuzz crate: `cd rust/fuzz && cargo +nightly-2026-09-22 fmt --check && cargo +nightly-2026-09-22 check`;
- a replay with debug assertions of every boolean corpus input and every
  `rust/fuzz/regressions/boolean*` file, one process per file: build with
  `CARGO_INCREMENTAL=0 RUSTFLAGS="-C debug-assertions" cargo +nightly-2026-09-22 build --release --example replay_boolean` in `rust/fuzz`;
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
  `TURNED_PARTS`; `GIVEN_CURVED`, `GIVEN_ROUND`, `GIVEN_BALL`, `SPLINE_WALLS`,
  `SPLINE_PARALLEL` and `SPLINE_CROSSING` are on). A heavily loaded
  host makes campaigns time out spuriously; run them on a quiet machine.
- **Keep debug-assertion tests optimized.** CI runs them with
  `CARGO_PROFILE_DEV_OPT_LEVEL=2`; time new test files that way.
- **Commits** are subject-only and long and descriptive, with no AI
  attribution; a staged diff must contain no machine paths. Push only after
  every check passes, and chain the push after the checks.
