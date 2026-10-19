# Working on the Rust kernel

How a step is done, how it is verified and the working rules learned so
far. The plan of record is `STATUS.md`; contracts and evidence are in the
per-capability documents it links.

## How a step is done

Each sub-step follows the same order, and the history shows it in its
commits:

1. **Refined decisions** before its code: why the case is refused today,
   the representation, the degrees and fields, what is `Degenerate`, what
   stays refused. Since R16 they go in the capability document (or the
   current review file), not in the frozen `REVIEW_NOTES.md`.
2. **Evidence**: an independent reference in `rust/tools/` (mpmath,
   rounding once; never `math.cos`, `sin` or `hypot`, whose results differ
   between Python versions and platforms), with two-way checks, closed
   forms and margins; fixtures from a generator with `--check` (identical
   under Python 3.9 and 3.12); a reference unit test; the generator in a
   fixture-generator group of `.github/workflows/rust-kernel.yml`
   (`rust/tools/generator_groups.py`).
3. **Native capture** before the kernel file exists: the compare script
   is keyed on that file, the capture's metadata records it absent;
   reviews in a divergences JSON only with independent evidence (an
   adaptive BRepGProp build, a Green's-theorem integration, a DRAW check),
   fingerprinted by the native row. Its Linux record follows from CI
   (`BOOLEAN.md`, "Linux records (CI)").
4. **Kernel**, with a test file: every fixture within `1e-9` enclosures
   of the reference, degenerate cases refused, histories complete, results
   deterministic and moved rigidly.
5. **Verification** (below), then merge and push, then CI green at the
   pushed revision before anything is called done (R15).
6. **DRAW survey** of the affected group (`UPSTREAM_TESTS.md`),
   registering cases that now evaluate within the 30 s contract with a
   volume audit; then a **campaign** (`FUZZING.md`).

Each step is worked in its own git worktree and branch, merged into the
integration branch `s9c2-kernel`, verified there and pushed to
`rust-kernel` as a fast-forward.

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

Then:

- the Boolean comparisons, each with 0 failures and the counts of
  `BOOLEAN.md`'s "Comparisons" table (all 55 by
  `rust/tools/run_boolean_comparisons.py --shard 0 --shards 1` with the
  pinned SDK's two arguments); since S9f.1 also `compare_split.py`
  (72 / 56) and `compare_brep.py --family spline` (10 / 3), which share
  its spline prisms, and `compare_brep_io.py` (6,835 / 7; with
  `--v3-dataset target/occt-test-data/opencascade-dataset-7.9.0` also 11
  version-3 matches, locally only); since S9e.4 `compare_step.py` (23 / 6,
  with STEP-b's SDK at `target/step-b-sdk`: the pinned SDK lacks the STEP
  toolkits' headers); `compare_tessellation.py` (42 / 14, and 21 / 3 with
  `--family spline`). Where `curved_boolean_reference.py`,
  `chained_curved_boolean_reference.py` or
  `torus_curved_boolean_reference.py` changes, every generator importing
  it runs with `--check --workers 4`;
- the fuzz crate: `cd rust/fuzz && cargo +nightly-2026-09-22 fmt --check
  && cargo +nightly-2026-09-22 check`;
- a replay with debug assertions of every corpus input and regression of
  the targets the change reaches (`boolean`, `split`, `identity`,
  `brep_validation` at least for a Boolean), one process per file: build
  with `CARGO_INCREMENTAL=0 RUSTFLAGS="-C debug-assertions" cargo
  +nightly-2026-09-22 build --release --example replay_boolean` in
  `rust/fuzz`; for a switch that is off (none is now), replay once more
  with it switched on in the source;
- `target/math-oracle-venv/bin/python -m unittest discover -s rust/tools
  -p 'test_*.py'` (the system Python lacks mpmath; CI runs the same suite
  in four shards by `run_tool_tests.py`, `FUZZING.md`);
- `python3 rust/tools/run_upstream_tests.py --ledger`, and for a change a
  DRAW case can reach the full contract run (`UPSTREAM_TESTS.md`'s "Run":
  `--backend both`, the dataset present, 30 s a case);
- a campaign: `python3 rust/tools/run_fuzz.py --target boolean --seconds
  600 --toolchain nightly --replay sample` (the full local corpus exceeds
  the startup hour under AddressSanitizer; the sampled replay plus the
  debug-assertion replay above stand in for it, and CI's schedule runs the
  full replay in shards).

## Working rules and lessons

- **Watch CI after every push** (R15): `gh run list -R
  jackControls/rustyOCCT -L 6`, and `gh run view RUN --log-failed` on any
  failure. A red workflow is the next task; record the head's state in
  `HANDOFF.md`'s header. Both workflows were red from S7 until 2026-09-29
  unnoticed, and six pushes went out red before the review of 18 October.
  Recurring causes so far:
  - a new native capture without its Linux record (`platform-linux/`,
    taken from the run's artifact and reviewed against the macOS capture);
  - reviews fingerprinted with macOS native rows only (the Linux rows
    differ in their last bits and need their own reviews);
  - job time limits as generators, fixtures and corpora grow (R17: any
    replay or minimisation over 45 minutes is sharded);
  - platform libm differences. Since R13 the kernel's constructions call
    no platform libm: every elementary function goes through `cr::Rn`
    (`MATHEMATICS.md`, "Correctly rounded elementary functions") and
    clippy denies `f64`'s libm methods in the library. Keep it so; a
    host-dependent bit in a fixture is a defect, not an allowance
    (`VALIDATION.md`'s cross-platform table holds native observations
    only).
- **Ubuntu's packaged OCCT is 7.6.3**, older than the pinned 8.1: a case
  may state that release's own outcome (`expected_occt_by_version`).
- **Fuzz timing.** AddressSanitizer costs about twelve times the plain time
  on the exact arithmetic, and Linux runners about 2.6 times the
  development Mac; the boolean target's limit is 60 s an input.
  Configurations too slow for it are switched off in
  `rust/fuzz/src/boolean.rs`, below the 35 s a run under the sanitizer the
  switches on keep to; none is off now (`TORUS_PAIRS`, `TURNED_PARTS`,
  `CONE_PAIRS`, `GIVEN_CURVED`, `GIVEN_ROUND`, `GIVEN_BALL`, `GIVEN_MET`,
  `SPLINE_WALLS`, `SPLINE_PARALLEL`, `SPLINE_CROSSING`, `SPLINE_SPHERE`,
  `SPLINE_CONE` and `IMPORTED` all on). A heavily loaded host makes
  campaigns time out spuriously; run them on a quiet machine.
- **Keep debug-assertion tests optimized.** CI runs them with
  `CARGO_PROFILE_DEV_OPT_LEVEL=2`; time new test files that way.
- **Commits** are subject-only, long and descriptive, with no AI
  attribution; a staged diff must contain no machine paths, host names or
  personal names. Push only after every check passes, and chain the push
  after the checks.
