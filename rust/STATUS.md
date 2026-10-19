# Status of the Rust kernel

The plan of record (R16 of `REVIEW_2026-10-18.md`): start here. The review
file holds the findings, decisions and next steps; the per-capability
documents linked below are the source of truth for contracts and evidence.
`REVIEW_NOTES.md` is a frozen archive. How a step is done, verified and
committed: [CONTRIBUTING-KERNEL.md](CONTRIBUTING-KERNEL.md). The CI state
of the head: [HANDOFF.md](HANDOFF.md)'s header (R15).

## What exists, accepted at which revision

| Work | Accepted at | Record |
|---|---|---|
| M0–M2 ids and histories on extrude and transform | `34efd36c` | [IDENTITY_AND_HISTORY](IDENTITY_AND_HISTORY.md) |
| T1 cell complex; M3 split and fuse; M4 attributes | `e4adb869`; `00f0034c`; `979cf939` | [TOPOLOGY_MODEL](TOPOLOGY_MODEL.md), IDENTITY_AND_HISTORY |
| M5 enclosures; T2 OCCT interop and the ledger | `0913b5e2` | the same, [UPSTREAM_TESTS](UPSTREAM_TESTS.md) |
| Exact root ordering; degree elevation | `de1d1d43`; `27647fad` | [VALIDATION](VALIDATION.md), [DEGREE_ELEVATION](DEGREE_ELEVATION.md) |
| Spline proximity; spline/linear preimages; B-rep validation | `d1206b15`; `30c5a247`; `dff912e5` | [SPLINE_PROXIMITY](SPLINE_PROXIMITY.md), [SPLINE_LINEAR_INTERSECTIONS](SPLINE_LINEAR_INTERSECTIONS.md), [BREP_VALIDATION](BREP_VALIDATION.md) |
| S1 closures; S2 `restore` and public data | `7e463cb2`; `1a76d29e` | IDENTITY_AND_HISTORY; UPSTREAM_TESTS |
| S3 cone, sphere, torus; S4 C1 check and spline cells | `166fc905`, `059616f9`, `ac92ec89`; `ac92ec89` | TOPOLOGY_MODEL, [MATHEMATICS](MATHEMATICS.md) |
| S5 arcs, S6 sheets and wires, S7 analytic intersections, S8 plane split, S9 Booleans | `54030783` (2026-10-14) | [BOOLEAN](BOOLEAN.md)'s Acceptance, [ANALYTIC_INTERSECTIONS](ANALYTIC_INTERSECTIONS.md), [SPLIT](SPLIT.md) |
| T-a, T-b tessellation; STEP-a, STEP-b import | implemented, own acceptance not recorded | [TESSELLATION](TESSELLATION.md), VALIDATION |

Since `54030783` (merge revisions on `rust-kernel`; detail in BOOLEAN.md,
UPSTREAM_TESTS.md and the review file):

- Boolean captures' Linux records on CI: the workflow `14a326de`, the 46
  records of run 38670188123 `512cdb5e`, each later set's after it; U10 (b)
  `b18183a8`.
- S9e.4b.4d.2 turned bodies (the user's word of 2026-10-14): 2a
  `7b43d0b4`, 2b and moved bodies `61162084`, 2c.1a `712115c6`, the DRAW
  survey `ba6186d5`, 2c.1b `5aaf69ce`, the consistency fix `a6f394cb`,
  2c.2 `679fdf39`, 2c.1c `12a6545e`; complete with campaigns at `5a7acb4c`.
- Near misses of other curves `5f50b30c`; the import tail (S9e.4b.4b.2b.2b.3b,
  S9e.4b.4c.2b.2) `37949a87`; the shallow lens `6e3423f9`; a Boolean's
  resolution `36cdd03b`; faces' near misses `75b925ae`; the torus meeting's
  verification budget `d86e8b74`; the tilted torus's validation `f2976ac8`;
  the DRAW survey at `1bc0c327` `216ca888` (1,132 registered cases hold).
- The review of 18 October: X1 `6d4c521f`; X2 (R13, correctly rounded
  constructions) `dacad00d` and `cb991b48`, X3 (R14) `3d0c89af`, both done
  at `36f83779` with every workflow green; X6 (R17) implemented `5e861738`;
  X4 and X5 (R15, R16: this document, the handoff's cut) pending the gate.

## Open items, in order

1. **CI at the head**: green at `50afeef4` (X6's CI changes; HANDOFF.md's
   header has the runs). Red is always the next task.
2. **X6's gate**: no scheduled job hitting its budget for a week, the
   Sunday minimisation of 25 October first ([FUZZING](FUZZING.md)'s
   sharded minimisation; the review's Status).
3. **S10** cancellation and work budgets; **S11** Revolve; **S12** Hole;
   **S13** STEP export (U11's order); **S14** mirror and
   orientation-reversing transforms; **S15** chamfers, then constant-radius
   fillets as procedural rolling-ball blends (U12). Each step's scope and
   evidence: the review's "Next development steps"; the jobs: [PORTING](PORTING.md).
4. **Acceptance records** of T-a, T-b, STEP-a and STEP-b: their gates
   (TESSELLATION.md's Acceptance, VALIDATION.md's STEP section).
5. **`rational_f64`'s ties to even**: the history's `Unchanged` hangs on
   tie coincidences (the review's F11; BOOLEAN.md's "Open items and refusals").
6. **Boolean items beyond acceptance** (imports not built, refusals
   diagnosed, the shallow lens's remainders, fuzz levers): BOOLEAN.md's
   "Open items and refusals".

Deliberate limits, each staying until a step takes it: the Booleans'
refusals by design (BOOLEAN.md); the DRAW adapter's command subset
(UPSTREAM_TESTS.md's "Deliberate limits") and `bug476_1` to `_8`
unregistered (its "Evaluating but not registered"); periodic spline
surfaces and seams on them unsupported, spline mass enclosures first order
where not exact (REVIEW_NOTES.md's S4 answer); BinTools' binary `.brep`
refused (UPSTREAM_TESTS.md, the version-3 header); procedural edges in
tessellation (T-c) pending (TESSELLATION.md); the public dataset local
only (U1, UPSTREAM_TESTS.md); no cancellation, WASM runtime or
performance gates yet ([PRODUCTION_READINESS](PRODUCTION_READINESS.md)).

## User decisions

None is open. Answered:

- **U9** (2026-10-14): tessellation's certified deflection stays the only mode.
- **U10** (2026-10-14): (b), C1 of the rational curve tested exactly.
- **U11** (2026-10-18): the default order S11, S12, S13.
- **U12** (2026-10-18): procedural rolling-ball blends for S15b.
- **U13** (2026-10-18): keep the Linux host's toolchain and sparse clone.

## How to run the gates

**CI** ([FUZZING](FUZZING.md) for the fuzz rules): "Rust kernel"
(`rust-kernel.yml`; every push and nightly: Ubuntu, macOS, Rust 1.85 and
Windows, a smoke subset per push; the fixture generators; the bridges), "Rust geometry fuzzing"
(`rust-fuzz.yml`; every push a sampled replay and campaign, daily the full
replays in shards, Sundays the sharded minimisation), "Rust Boolean
captures on Linux" (`rust-boolean-linux.yml`; daily and on dispatch, all 55
comparisons strictly against their Linux records). Check with `gh run list
-R jackControls/rustyOCCT -L 6`.

**Locally** (the commands in CONTRIBUTING-KERNEL.md): fmt, clippy, the
Rust 1.85 check, the release tests; the 55 Boolean comparisons by
`rust/tools/run_boolean_comparisons.py --shard 0 --shards 1` or each
script (counts in BOOLEAN.md's "Comparisons"), with `compare_split`,
`compare_brep`, `compare_brep_io`, `compare_step` and
`compare_tessellation`; the corpus replays with debug assertions
(`rust/fuzz/examples/replay_*.rs`); the `rust/tools` unit tests; campaigns
by `rust/tools/run_fuzz.py --target T --seconds 600` (FUZZING.md); the DRAW
contract run, `rust/tools/run_upstream_tests.py --backend both` with the
dataset and `--ledger` ([UPSTREAM_TESTS](UPSTREAM_TESTS.md)).

**The green-gate rule (R15).** A red workflow is the next task, always.
No "done", "verified" or "accepted" record is written while a workflow is
red at that revision; the handoff's header states the head's CI state; a
push that turns a workflow red is followed only by commits that address it
or revert it. A fixture-generator group skipped on identical inputs counts
as green (FUZZING.md).

**CI wall times** (X6's measurements of 11 to 18 October, before its
shards; FUZZING.md's "CI wall times (R17)" has each target's): a push ran
the fuzz workflow's 31 jobs in 250 to 274 job-minutes (39 to 79 minutes of
wall time) and the kernel workflow's 56 to 58 jobs in 910 to 1,045
job-minutes, 561 to 713 of them fixture generators, which X6 now skips on
identical inputs (some 550 job-minutes on a kernel-only push), its tool
tests in four shards of about 486 s; a daily fuzz run took 65 jobs, 1,284
to 1,320 job-minutes (91 to 119 minutes); the Sunday minimisations 792 and
886 job-minutes (115 and 125 minutes, `surface_knots`' 106-minute job the
longest). The times with X6's shards are measured from `50afeef4`'s runs on.
