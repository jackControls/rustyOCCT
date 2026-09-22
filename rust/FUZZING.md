# Sustained geometry fuzzing

`rust/fuzz` is a separate cargo-fuzz workspace. Its fuzz-specific dependencies and C++
libFuzzer runtime are test-only. The production kernel remains native Rust.
Both workspaces have checked-in dependency locks. The campaigns run optimized
code with AddressSanitizer, debug assertions, overflow checks, and libFuzzer
coverage feedback. They are distinct from the deterministic invariant tests.

| Target | Generated inputs | Assertions |
| --- | --- | --- |
| `predicates` | Full binary64 bit patterns, scaled integers, coplanarity, sphere-boundary points, non-finite coordinates | Exact 2D/3D orientation and insphere against independent rational matrix elimination; sphere ordering; typed rejection |
| `intersections` | Lines/segments, three-point planes/triangles, full binary64 exponents, coplanarity and degeneracies | All three intersection APIs against a rational barycentric oracle; exact classification; minimal finite coordinate/parameter enclosures; explicit unrepresentable results |
| `modeling` | Valid radial polygons/circles, optional holes, 49 scales, up to eight operations, plus raw invalid input | Repeated rigid transforms, reversed winding/offsets and planar splits; mass/first-moment conservation, topology, classification, bounds and finite positive properties |
| `curved` | Full binary64 coefficients, centers/radii/axes and line endpoints; scaled integers; exact/neighboring tangencies; generator and point segments | Quadratic root count, multiplicity, exact comparisons, minimal enclosures; circle/sphere/cylinder hits against independent polynomial-sign and axial-projection oracles; endpoint clipping and typed failures |
| `splines` | Raw binary64 poles/weights/knots, scaled geometry, degree 1..25, repeated knots, endpoint/interior parameters, explicit sides and derivative requests | Exact basis-function derivatives and closed quotient formulas independently check homogeneous pole interpolation; minimal position/derivative bounds; discontinuity, domain, nonfinite and overflow errors |

The rational oracle uses `num-rational` with Gaussian elimination and
barycentric coordinates, plus polynomial-sign/vertex comparisons and cylinder
axial projection. Production uses a fixed binary64 integer lattice, determinant
expansion, edge half-planes, exact radical comparisons and cylinder cross
products. Spline production uses differentiated de Boor pole interpolation;
its fuzz oracle instead evaluates basis functions and their derivative identity.
The Rust oracles and production share `num-bigint`/`num-rational`;
the checked-in Python `Fraction` fixtures provide a separate integer runtime.
Coverage counts include the oracle and dependencies: they are not kernel-only
coverage percentages or evidence of exhaustive input coverage.

## Continuing campaigns

[Rust geometry fuzzing](../.github/workflows/rust-fuzz.yml) runs all five targets
for 60 seconds of mutation each on relevant pushes/PRs, and 600 seconds each every day at
06:23 UTC on the default branch. Manual runs accept 1–3,600 seconds per target.
GitHub can delay scheduled jobs. The schedule must remain enabled on the fork.

Every run restores the previous evolving corpus, adds seeds derived from the
exact fixtures, and retains newly discovered inputs for the next run. PRs read
the corpus but do not publish corpus caches. All runs upload logs, JSON metadata,
corpora, and crash/timeout/OOM artifacts for 30 days, including failed runs.
Cache eviction does not remove the checked-in fixtures or regressions.

Each input has a 20-second limit, a 2 GiB process RSS limit, and at most 256
bytes. The modeling harness bounds geometry to 24 vertices and eight operations;
spline inputs have at most 51 poles (the kernel allows 4096). An
outer deadline kills the build/fuzzer process group. Corpus replay without any
subsequent mutation is an incomplete run. An incomplete run, crash,
timeout, OOM, changed dependency lock or mathematical disagreement fails CI.
These are harness limits, not kernel production performance guarantees.

The mutation timer starts when the pinned libFuzzer reports `INITED`, after
corpus replay. Its `max_total_time` flag includes initialization and previously
allowed a growing corpus to consume the entire short campaign; this was caught
as an incomplete CI run. The runner now creates a nonempty `stop_file` after
the full requested mutation budget. LibFuzzer stops normally and emits final
statistics. The outer deadline still bounds startup/build plus mutation to
`requested_seconds + 600`; a missing initialization marker or ignored stop
request cannot hang indefinitely. Early exits never count as completed budgets.

## Local reproduction

On a cargo-fuzz-supported Unix host with a C++ compiler:

```sh
rustup toolchain install nightly-2026-09-22 --profile minimal --component rustfmt
cargo install cargo-fuzz --version 0.13.1 --locked
python3 rust/tools/run_fuzz.py --seconds 300
# Reproduce a particular random campaign start (corpus state also matters):
python3 rust/tools/run_fuzz.py --target intersections --seconds 600 --seed 314159
```

The nightly toolchain is only for fuzzing; normal builds retain Rust 1.85 as
their minimum. `--toolchain` permits deliberate local toolchain experiments.
`target/fuzz-reports/summary.json` records the compiler, fuzzer, lock hash,
revision/dirty state, command, startup/mutation timing, budget completion,
replay/mutation execution counts, coverage edges and artifact
names. The raw log records the random seed. Cargo-fuzz 0.13.1 has no `--locked`
run option: the runner fetches with `--locked`, builds offline, and verifies
that the lock remains unchanged.

To replay and minimize a discovered failure (replace the example path):

```sh
cargo +nightly-2026-09-22 fuzz run intersections rust/fuzz/artifacts/intersections/crash-HASH --fuzz-dir rust/fuzz
cargo +nightly-2026-09-22 fuzz tmin intersections rust/fuzz/artifacts/intersections/crash-HASH --fuzz-dir rust/fuzz -- -max_total_time=120
```

Keep the original artifact. Investigate whether the defect is in the kernel,
the oracle, or its input contract. Add the minimized bytes under
`rust/fuzz/regressions/<target>/*.bin`, with a short explanation and a readable
ordinary Cargo regression. The runner automatically seeds those bytes on every
campaign. Fix the underlying defect before changing expected values or bounds.
Minimization is an explicit triage step, not a claim that CI automatically
understands or repairs failures. No lifetime reliability guarantee follows
from any finite campaign.

If corpus replay starts exhausting the outer startup allowance, compact it with
`cargo +nightly-2026-09-22 fuzz cmin <target> --fuzz-dir rust/fuzz` and retain
the uncompressed artifact until the compacted corpus is validated. A stalled
campaign must be repaired, not counted as a successful mutation run.
