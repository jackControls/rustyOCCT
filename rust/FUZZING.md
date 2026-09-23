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
| `linear_sets` | All 25 linear primitive pairings, raw binary64, coplanar/shared-vertex modes, full exponent range, invalid definitions and collapsed segments | Complete exact sets from independent boundary crossings and a gift-wrapping hull; canonical operand/winding symmetry; minimal construction bounds; zero-distance consistency on scaled-integer modes |
| `proximity` | All 25 point/line/segment/plane/triangle pairings; arbitrary binary64 coordinates, scaled/coplanar/shared-vertex cases, collapsed segments and invalid inputs | Independent exact membership and supporting-plane certificates prove each returned pair globally minimal; minimal output enclosures, threshold comparisons, operand symmetry, coordinate permutation/vertex reversal and deterministic witnesses |
| `modeling` | Valid radial polygons/circles, optional holes, 49 scales, up to eight operations, plus raw invalid input | Repeated rigid transforms, reversed winding/offsets and planar splits; mass/first-moment conservation, topology, classification, bounds and finite positive properties |
| `curved` | Full binary64 coefficients, centers/radii/axes and line endpoints; scaled integers; exact/neighboring tangencies; generator and point segments | Quadratic root count, multiplicity, exact comparisons, minimal enclosures; circle/sphere/cylinder hits against independent polynomial-sign and axial-projection oracles; endpoint clipping and typed failures |
| `splines` | Raw binary64 poles/weights/knots, scaled geometry, degree 1..25, repeated knots, periodic seams, full-range wrapped parameters, explicit sides and derivative requests | Exact basis-function derivatives and closed quotient formulas independently check homogeneous pole interpolation; minimal position/derivative bounds; discontinuity, domain, nonfinite and overflow errors |
| `bezier_editing` | Raw binary64 and scaled rational controls, degree 1..25, clamped/unclamped/periodic splines, subnormal spans, multi-period extraction and composed edits with rational cuts | Complete homogeneous polynomial identities from independent Cox basis coefficients and affine substitution; exact jets/minimal bounds, positive weights, shared endpoints, reversal, elevation/split commutation and preflight limits |
| `knot_editing` | Degrees 1..25, raw binary64 and scaled controls, rational cuts, unclamped inactive controls, repeated knots, periodic seams/origin changes, refinement/removal sequences and invalid data | Entire raw-support homogeneous identities; independent coefficient-equation removal feasibility, exact complete controls, rational jets/extraction, minimal enclosures, batch order/duplicate behavior, round trips and unchanged rejection limits |
| `exact_spline_intersections` | Rational controls/weights and knots, degree 1..25, rational trims, close roots through 2^-2048 spacing, huge/subnormal domains, periodic seams and edit sequences | Complete known-factor contacts, rational secants and exact circle overlaps; independent full polynomial edit identity; exact parameter/coordinate comparisons, contact orders through 50, minimal finite bounds or typed conversion failure, traversal limits and malformed-rational rejection |
| `surface_editing` | Tensor degrees 1..25 in both directions, raw binary64 and scaled controls, independent periodicity, unclamped knots, subnormal domains, full low-degree multi-period queries and high-degree selected spans, composed edits and rational cuts | Complete tensor polynomial identities; all exact partials through order two and minimal bounds; exact isocurves/shared boundaries, reversals, transposition, elevation/split commutation and Cartesian output limits |
| `surfaces` | Rational tensor grids, independently periodic U/V, repeated knots, high degree in either direction, full-range wrapped parameters and malformed data | Independent tensor basis plus closed bivariate quotient formulas; exact mixed-partial and quadrant continuity decisions; minimal enclosures and typed failures |
| `roots` | Products of rational/irrational/complex factors through degree 25, repeated roots, closed-domain clipping, power-of-two coefficient scaling and arbitrary binary64 query polynomials | Complete expected root list from known factors; algebraic signs reduce independently in Q(sqrt(d)); multiplicities, exact comparisons, minimal enclosures and nonfinite rejection |
| `spline_intersections` | Rational Bézier curves with known factored plane numerators and squared sphere/cylinder contact equations; weighted rational lines against quadrics; nonperiodic/periodic rational polylines; varying weights, parameter/space scales, oblique planes, tangencies, knots and zero spans; explicit trim bounds, neighboring floats and large periodic offsets | Complete parameter/position bounds from rational Bernstein evaluation, affine span equations, or independent geometric quadratic roots with rational weight-parameter conversion; one-sided orders, crossing/tangent/boundary classification, maximal clipped overlaps, closed seam events and repeated turns |

The rational oracle uses `num-rational` with Gaussian elimination and
barycentric coordinates, plus polynomial-sign/vertex comparisons and cylinder
axial projection. Production uses a fixed binary64 integer lattice, determinant
expansion, edge half-planes, exact radical comparisons and, for analytic lines,
cylinder cross products. Spline production uses differentiated de Boor pole interpolation;
its fuzz oracle instead evaluates basis functions and their derivative identity.
Spline/quadric line fuzzing independently solves the physical line's quadratic
and converts its parameter through the rational weight map; production instead
isolates the spline's homogeneous implicit polynomial. Separate Python quadric
fixtures use cylinder cross products to check production's dot-product form.
The proximity checker uses convex supporting-plane inequalities, independently
of the production face enumeration and exact normal-system solve. Python
fixtures additionally use closed analytic projection/cross-product formulas.
Complete linear-set fuzzing uses cross-product line/plane formulas and boundary
edge candidates with a gift-wrapping hull; production solves affine equalities
and enumerates feasible halfspace vertices before a monotone-chain hull.
Editing production uses spline blossoms (curves), exact boundary knot insertion
(tensor extraction) and homogeneous de Casteljau;
its checker uses basis polynomials, binomial affine substitution and Bernstein
coefficient expansion. Clearing common denominators before those linear
transforms reduces repeated GCD work without changing the assertions.
Knot editing uses exact Boehm insertion and inverse insertion in production;
its independent checker solves Cox power coefficient equations with integer
elimination, including inactive unclamped support and failed removals.
The Rust oracles and production share `num-bigint`/`num-rational`;
the checked-in Python `Fraction` fixtures provide a separate integer runtime.
Coverage counts include the oracle and dependencies: they are not kernel-only
coverage percentages or evidence of exhaustive input coverage.

## Continuing campaigns

[Rust geometry fuzzing](../.github/workflows/rust-fuzz.yml) runs all fourteen targets
for 60 seconds of mutation each on relevant pushes/PRs, and 600 seconds each every day at
06:23 UTC on the default branch. Manual runs accept 1–3,600 seconds per target.
GitHub can delay scheduled jobs. The schedule must remain enabled on the fork.

Every run restores the previous evolving corpus, adds seeds derived from the
exact fixtures, and retains newly discovered inputs for the next run. PRs read
the corpus but do not publish corpus caches. All runs upload logs, JSON metadata,
corpora, and crash/timeout/OOM artifacts for 30 days, including failed runs.
Cache eviction does not remove the checked-in fixtures or regressions.

Each input has a 20-second limit and a 2 GiB process RSS limit. Surface editing
permits 4096 bytes to populate full tensor grids; knot editing permits 512 bytes;
other targets permit 256
bytes. The modeling harness bounds geometry to 24 vertices and eight operations;
curve spline/Bézier-editing inputs have at most 51 poles; knot editing starts with
up to 101 poles and can insert up to 50 more; surfaces have up to 108 poles with
degree 25 in either direction. The surface-editing target allows up to 784
poles, including degree 25 in both axes (the kernel allows 4096 total poles).
It queries one selected knot rectangle for high-degree grids and also full
domains/multiple turns for degrees whose product is at most 25. Full high-degree
decompositions and all edits remain in the complete native coefficient bridge;
ordinary fixtures retain extraction and composed high-degree edits for all four
periodicity combinations. A Bézier
known-factor intersection input has degree at most eight, and a rational polyline
has up to eight spans. Power-curve quadric inputs reach degree 25 (degree-50
equations), checked by independent monotonic rational powers. Independent exact
fixtures additionally cover dense degree-25 intersections.
The exact edited-curve target starts with at most 26 poles and refines to at
most 52. It checks degree-25 known-factor plane numerators and squared quadric
contacts through order 50, rational secants, periodic line crossings and
exact periodic circle overlaps. Parameter domains include rational offsets and
widths through 2^2048 and 2^-2048; positive rational weights, common homogeneous
scaling and coordinate permutations exercise conditioning independently of
root identity. These cases retain exact results when floating conversion fails.
General root products reach degree 25. An
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

Slow-unit diagnostics also receive triage even when a campaign passes. The
checked-in [regression notes](fuzz/regressions/README.md) link saved inputs to
their originating run and ordinary exact fixtures. A release replay example
can time complete spline-intersection oracle checks separately from sanitizer
instrumentation; these timings are diagnostic, not a production latency gate.

If corpus replay starts exhausting the outer startup allowance, compact it with
`cargo +nightly-2026-09-22 fuzz cmin <target> --fuzz-dir rust/fuzz` and retain
the uncompressed artifact until the compacted corpus is validated. A stalled
campaign must be repaired, not counted as a successful mutation run.
