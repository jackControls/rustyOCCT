# Handoff: the Rust kernel's S9 Booleans

State of the work for whoever continues it. `REVIEW_NOTES.md` remains the
plan of record: its decisions, evidence, implemented, survey and campaign
bullets are authoritative, and this page only summarizes them and says
where to pick up. Written 2026-10-01, brought up to date 2026-10-03 for
S9f.2a, S9f.2b.1 and S9e.3b (branch `s9e3b`, `s9c2-kernel` merged in)
and their DRAW survey (branch `s9-draw-3`).

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
  DRAW survey (the last full survey, that of S9e.3b, S9f.2a and S9f.2b.1
  at `b0b9adc6`, branch `s9-draw-3`: 987 cases registered, none failing or
  timing out, the volume audit's values those of `507b8054`'s bit for
  bit) and a clean campaign (S9e.3a's and S9f.1's at `b2765f20`, 989
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
- **S9e.3b implemented** (branch `s9e3b`, over `507b8054`, with
  `s9c2-kernel` at `c8e37abe` merged in; not yet pushed): a given result's
  meeting of two curved faces (`Meet`, `Rise`, `Toric`) or a cone's or
  torus's general section met by the partner (`curved/triple.rs`), its
  decisions, evidence and capture committed before the kernel; the DRAW
  survey of the chained cases done (none reaches it), and the full DRAW
  survey (that of S9e.3b, S9f.2a and S9f.2b.1: no status or refusal
  changes, the volume audit's values bit for bit). Pending the campaign.
  A given `WallMeet` edge never reaches it (a given result with spline
  walls is refused before, S9f).
- **CI.** Both workflows ("Rust kernel", "Rust geometry fuzzing") were green
  at `6c221525`. They had been red from S7 until 2026-09-29, unnoticed;
  check them after every push (see "Working rules").

## What is open, in order

1. **S9e.3b's campaign** (`REVIEW_NOTES.md`'s "S9e.3b implemented"; the
   DRAW survey of S9e.3b, S9f.2a and S9f.2b.1 is done, branch
   `s9-draw-3`: no status or refusal changes, the volume audit's 987 bit
   for bit). The fuzz target's `GIVEN_MET` is off: the corpus's slowest
   chained operations reaching a given meeting take 60 to 71 s an input
   under AddressSanitizer. The campaign runs with it off; the next lever is the degree-eight
   arrangement arithmetic (the second arrangement's predicates on vertices
   in fields of degree eight), after which `GIVEN_MET` can be switched on
   and a campaign run with it. Of the chained operations the debug replay
   with it on reaches, the "points not separated by a projection" limits
   (`triple.rs`'s retries exhausted: a fibre of two points under every
   shear tried) are the one open refusal worth a closer look.
2. **S9f.2b.2** (loops round a turning point inside the faces, each piece
   about it a graph over the cylinder's angle in a verified window, and the
   tower fields of a cap circle on a wall along its axis), then **S9f.3**
   (spheres and cones; "S9f refined" gives the degrees).
3. **S9e.4**: imported bodies decided on their stored surfaces.
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
| `compare_given_met_boolean.py` | 8 / 42 (the kernel within the reference on all 42 solid cases, the 8 degenerate refused) |

Every one must report 0 failures; since S9f.1 also `compare_split.py`
(72 / 56) and `compare_brep.py --family spline` (10 / 3), which share its
spline prisms, and `compare_brep_io.py` (6,835 / 7). Where `curved_boolean_reference.py` or
`chained_curved_boolean_reference.py` changes, every generator importing it
runs with `--check --workers 4`. Then:

- the fuzz crate: `cd rust/fuzz && cargo +nightly-2026-09-22 fmt --check && cargo +nightly-2026-09-22 check`;
- a replay with debug assertions of every boolean corpus input and every
  `rust/fuzz/regressions/boolean*` file, one process per file: build with
  `CARGO_INCREMENTAL=0 RUSTFLAGS="-C debug-assertions" cargo +nightly-2026-09-22 build --release --example replay_boolean` in `rust/fuzz`;
  for a switch that is off (`GIVEN_MET`, `TORUS_PAIRS`, ...), replay once
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
  `TURNED_PARTS`, `GIVEN_MET`; `GIVEN_CURVED`, `GIVEN_ROUND`, `GIVEN_BALL`,
  `SPLINE_WALLS`, `SPLINE_PARALLEL` and `SPLINE_CROSSING` are on). A heavily loaded
  host makes campaigns time out spuriously; run them on a quiet machine.
- **Keep debug-assertion tests optimized.** CI runs them with
  `CARGO_PROFILE_DEV_OPT_LEVEL=2`; time new test files that way.
- **Commits** are subject-only and long and descriptive, with no AI
  attribution; a staged diff must contain no machine paths. Push only after
  every check passes, and chain the push after the checks.
