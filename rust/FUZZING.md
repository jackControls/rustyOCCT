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
| `spline_proximity` | Positive-weight rational space curves through degree 25 with known zero and nonzero minima, nonplanar cylinder curves, rational/irrational factor roots, composed exact knot/degree edits, trimmed/singleton ranges, near-tie parabolas, rational circles and periodic polylines | Complete minimum parameter sets and whole intervals from independent factor and radical signs, a monotone cubic, circle identities and exhaustive rational segment projections; tight parameter/coordinate/distance bounds; exact ties and periodic aliases |
| `spline_linear` | Known-factor rational Bézier curves through degree 25 with sub-float root pairs, rational polylines with periodic turns, constant spans and retracing, a rational circle against rational slopes, a quadratic retrace, lines/collapsed/reversed segments, invertible integer affine maps, parameter domains through 2^±2048 and exact knot/degree edits | Complete point/interval sets from constructed factors, per-span rational linear inequalities and closed-form quadratic sign functions; exact parameter, coordinate and line-parameter identity within 2^-96 independent brackets; tight or typed unrepresentable views; reversal symmetry, malformed-rational and work-limit rejection |
| `brep_validation` | Valid star-outline prisms with no hole, round or square holes, or an inverted box cavity merged by fuzz-crate code, at scales through 2^±10 with random frames; one of sixteen mutations: extra vertex/edge, empty shell/loop, bad references, dropped or repeated faces, flipped uses and faces, moved vertices, shifted pcurves, inverted shells and swapped loops | Clean report for every base; exact issue lists for local mutations and a required issue on the mutated entity otherwise; a `10·tol` pcurve shift clean at `100·tol`; deterministic, duplicate-free reports identical to `from_parts` |
| `identity` | Structure-aware (`arbitrary::Unstructured`) profiles: 3–12-sided or circular outlines, clockwise or counter-clockwise input, up to three square or round holes, labels absent or present, random operation ids, frames, directions, scales through 2^±8 and up to two rigid motions | An independent version-1 encoder and FNV-1a-128 recompute every id from its derivation; ids unique and slot maps inverse; rebuilding, rigid motion, reversed direction and moved labelled points keep ids; permuting label values permutes parents bijectively; unlabelled builds share no id; entity counts and roles follow the profile |
| `modeling` | Valid radial polygons/circles, optional holes, 49 scales, up to eight operations, plus raw invalid input | Repeated rigid transforms, reversed winding/offsets and planar splits; mass/first-moment conservation, topology, classification, bounds and finite positive properties |
| `curved` | Full binary64 coefficients, centers/radii/axes and line endpoints; scaled integers; exact/neighboring tangencies; generator and point segments | Quadratic root count, multiplicity, exact comparisons, minimal enclosures; circle/sphere/cylinder hits against independent polynomial-sign and axial-projection oracles; endpoint clipping and typed failures |
| `splines` | Raw binary64 poles/weights/knots, scaled geometry, degree 1..25, repeated knots, periodic seams, full-range wrapped parameters, explicit sides and derivative requests | Exact basis-function derivatives and closed quotient formulas independently check homogeneous pole interpolation; minimal position/derivative bounds; discontinuity, domain, nonfinite and overflow errors |
| `bezier_editing` | Raw binary64 and scaled rational controls, degree 1..25, clamped/unclamped/periodic splines, subnormal spans, multi-period extraction and composed edits with rational cuts | Complete homogeneous polynomial identities from independent Cox basis coefficients and affine substitution; exact jets/minimal bounds, positive weights, shared endpoints, reversal, elevation/split commutation and preflight limits |
| `knot_editing` | Degrees 1..25, raw binary64 and scaled controls, rational cuts, unclamped inactive controls, repeated knots, periodic seams/origin changes, refinement/removal sequences and invalid data | Entire raw-support homogeneous identities; independent coefficient-equation removal feasibility, exact complete controls, rational jets/extraction, minimal enclosures, batch order/duplicate behavior, round trips and unchanged rejection limits |
| `surface_knots` | Tensor degrees 1..25 in both axes, independent periodicity and raw unclamped support, arbitrary binary64 and scaled controls, rational insertions, batch/transpose identities, inverse edits and periodic origin deletion | Every transverse homogeneous Cox coefficient; shared independent fraction-free removal equations; full extracted patch coefficients and quotient jets; retained isocurves, parameter wrapping, invalid rational rejection and atomicity |
| `degree_elevation` | Curve/tensor degrees 1..25; clamped, unclamped, inactive and periodic axes; raw binary64 and extreme rational atoms; invalid targets and full 4096-control boundaries | Complete independent Cox-equation control reconstruction; original domains, positive weights, staged edits, refinement composition, tensor order/transposition, exact jets and typed degree/count rejection |
| `exact_spline_intersections` | Rational controls/weights and knots, degree 1..25, rational trims, close roots through 2^-2048 spacing, huge/subnormal domains, periodic seams and edit sequences | Complete known-factor contacts, rational secants and exact circle overlaps; independent full polynomial edit identity; exact parameter/coordinate comparisons, contact orders through 50, minimal finite bounds or typed conversion failure, traversal limits and malformed-rational rejection |
| `surface_editing` | Tensor degrees 1..25 in both directions, raw binary64 and scaled controls, independent periodicity, unclamped knots, subnormal domains, full low-degree multi-period queries and high-degree selected spans, composed edits and rational cuts | Complete tensor polynomial identities; all exact partials through order two and minimal bounds; exact isocurves/shared boundaries, reversals, transposition, elevation/split commutation and Cartesian output limits |
| `surfaces` | Rational tensor grids, independently periodic U/V, repeated knots, high degree in either direction, full-range wrapped parameters and malformed data | Independent tensor basis plus closed bivariate quotient formulas; exact mixed-partial and quadrant continuity decisions; minimal enclosures and typed failures |
| `roots` | Products of rational/irrational/complex factors through degree 25, repeated roots, closed-domain clipping, power-of-two coefficient scaling and arbitrary binary64 query polynomials | Complete expected root list from known factors; algebraic signs reduce independently in Q(sqrt(d)); multiplicities, cross-equation root ordering and equality against independently constructed factors, minimal enclosures and nonfinite rejection |
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

[Rust geometry fuzzing](../.github/workflows/rust-fuzz.yml) runs all twenty targets
for 60 seconds of mutation each on relevant pushes/PRs, and 600 seconds each every day at
06:23 UTC on the default branch. Manual runs accept 1–3,600 seconds per target.
GitHub can delay scheduled jobs. The schedule must remain enabled on the fork.

Every run restores the previous evolving corpus, adds seeds derived from the
exact fixtures, and retains newly discovered inputs for the next run. PRs read
the corpus but do not publish corpus caches. All runs upload logs, JSON metadata,
corpora, and crash/timeout/OOM artifacts for 30 days, including failed runs.
Cache eviction does not remove the checked-in fixtures or regressions.

Each input has a 20-second limit, except complete surface-knot and degree-elevation verification,
which has 60 seconds. All targets retain the 2 GiB process RSS limit. The new
tensor target checks both axes at degree 25 and every raw-support homogeneous
equation; the densest retained input takes approximately 23 seconds with
instrumentation on the development machine. Its larger verification budget is
explicit in each campaign report. Surface editing
and surface knot/degree editing permit 4096 bytes to populate full tensor grids; curve knot editing permits 512 bytes;
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
Spline minimum-distance inputs additionally offset known-factor curves so the
positive squared minimum is exactly `h²`. Their stationary equations can reach
degree 73, while an attained affine-hull bound can certify minima directly. Rational
parameter ties reach degree 25, while positive-distance irrational ties and
composed edits start at degree at most eight, with optional elevation to nine.
The independent factor identities still certify every minimizing parameter;
all original zero-distance degree-25 cases remain. Nonplanar rational cylinder
curves additionally have `D=h²+scale²*s²*(s²-a)²`, with `0<a<1`. Their full 3D
affine hull has an unattainable zero lower bound, so non-singleton ranges retain
the general stationary solver. Independently known roots and coordinates
certify all minima through parameter changes, translations, axis permutations,
degree-five-to-six elevation, knot insertion and closed trims. The input and
RSS budgets are unchanged; the original degree-24 timeout seed remains included.
Spline/linear inputs use at most 128 bytes and the standard 20-second/2 GiB
limits. Their known-factor curves reach degree 25; polylines can be elevated to
degree 25 after optional knot insertion. Each input also runs the reversed
query and typed rejection checks. Early development replay of 3,000 random
inputs found known-factor sub-float root pairs taking up to 65 seconds. Profiling
traced this to schoolbook integer gcds during polynomial content removal and to
late rational-root recognition. The binary gcd and rational-root-theorem
candidate (see [the capability](SPLINE_LINEAR_INTERSECTIONS.md)) reduced that
input to 0.79 seconds in release, without changing any limit or assertion.
General root products reach degree 25. An
outer deadline kills the build/fuzzer process group. Corpus replay without any
subsequent mutation is an incomplete run. An incomplete run, crash,
timeout, OOM, changed dependency lock or mathematical disagreement fails CI.
These are harness limits, not kernel production performance guarantees.

B-rep validation inputs use at most 256 bytes and the standard 20-second/2 GiB
limits. Every prism is built by `Solid::extrude` and validated before mutation,
so the base's validity is also re-checked on every input.

The surface-knot runner enables the test-only `asan-allocator` feature together
with AddressSanitizer. At most once per second, after a complete input and all
its mathematical checks have returned, the harness invokes
`__sanitizer_purge_allocator`. The pinned libFuzzer already invokes this API
during mutation, but omits it during seed replay. The additional call covers
replay as well. It drains freed-allocation quarantine and releases unused pages;
it does not free live geometry or change the RSS limit. This does shorten the
quarantine across independent inputs, as the runtime's existing purge does;
the normal sanitizer checks remain active throughout each input. Ordinary
replay binaries and the production kernel have no sanitizer FFI dependency.
See the pinned runtime's `FuzzerLoop.cpp::{ReadAndExecuteSeedCorpora,PurgeAllocator}`
and [LLVM's allocator implementation](https://github.com/llvm/llvm-project/blob/main/compiler-rt/lib/asan/asan_allocator.cpp).

This target also sets `ASAN_OPTIONS=quarantine_size_mb=64`, retaining a 64 MiB
freed-block quarantine while keeping the same 2 GiB process gate. Other targets
use their existing sanitizer settings. Default quarantine plus allocator
cleanup still exhausted macOS RSS during retained-corpus replay despite roughly
38 MB of live allocations at the failure. Ordinary Rust replay of the same 197
inputs twice peaked at approximately 62 MiB. The smaller quarantine is a
documented instrumentation tradeoff: it can miss a stale-pointer access after
its freed allocation leaves quarantine sooner. It does not suppress a reported
error, remove a mathematical assertion, or raise the process memory cap.
Campaign manifests record the effective sanitizer options. These measurements
diagnose this corpus and runtime; they do not establish a general kernel memory
bound.

`spline_linear` uses the same two settings, for the same measured reason. Its
second clean local 600-second campaign at `809f8b4b` stopped with an OOM at
2,055 MB. RSS was already about 1.9 GB when the 535-input corpus replay ended.
The saved OOM input alone peaks at 65 MB under AddressSanitizer and 3.6 MB
without it. The campaign-wide growth is therefore freed-block quarantine and
allocator retention from exact BigInt churn, not live data. The input limit,
2 GiB gate and all assertions are unchanged.

The mutation timer starts when the pinned libFuzzer reports `INITED`, after
corpus replay. Its `max_total_time` flag includes initialization and previously
allowed a growing corpus to consume the entire short campaign; this was caught
as an incomplete CI run. The runner now creates a nonempty `stop_file` after
the full requested mutation budget. LibFuzzer stops normally and emits final
statistics. Build and corpus replay have a separate deadline of
`min(3600, 600 + input_limit_seconds * (initial_corpus_files + 1))` seconds.
This reserves build time and each input's configured allowance, including the
initial empty input, subject to a one-hour cap. An earlier average-cost estimate
exhausted startup on a 385-input surface corpus. Every saved input remains in
replay, with its own timeout, and mutation still receives its full separate
budget. A late `INITED` marker cannot borrow mutation or
shutdown time. The pinned libFuzzer checks `stop_file` between mutation batches,
each containing up to five callbacks. The runner explicitly keeps
`mutate_depth=5` and allows five times the per-input timeout plus five seconds
for that final batch and its statistics: 105 seconds normally, 305 for surface
knots. The process group is bounded by the startup allowance plus the requested
mutation duration plus that final-batch grace. Individual input limits stay at
20/60 seconds, and ignored stop requests are killed.
An early exit, startup overrun, missing final statistics or absent mutations
still fails the campaign. Reported peak RSS or slowest input exceeding its limit
also fails, even if the runtime exits successfully: a Linux seed took 64 seconds
against a 60-second alarm without a nonzero exit. Separating the phases fixes a reproduced Linux run
that completed its mutation budget but was killed before its last input and
final statistics finished; the failed evidence remains retained.

Surface-knot coefficient checks use exact continuity induction across the full
common raw-support partition. Equality on the preceding span and the known
knot multiplicities prove the lower coefficients; explicit integer comparisons
prove the rest. Periodic curves require all first-span coefficients because
they have no exterior zero polynomial. This remains a complete identity proof,
with adversarial comparisons against the exhaustive coefficient checker; it
does not substitute sampled points or approximate comparisons.

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

For the surface-knot target, set `ASAN_OPTIONS=quarantine_size_mb=64` and include
`--sanitizer address --features asan-allocator` before the final `--` when
reproducing the campaign's allocator behavior.

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


Surface knot editing starts with up to 784 controls, preserving every complete
transverse coefficient identity after each edit. Its independent checker shares
basis construction and equation factorization across transverse fields rather
than repeating them per row. Every field and residual equation is retained.
Pairs of identical complete before/after control columns share one projection
only after comparing every integer in both columns. Any changed field creates
its own pair; tests corrupt every control component to verify rejection.
The fifteen-target daily workflow includes this campaign with the same startup,
mutation and memory limits; its larger per-input limit is documented above.
