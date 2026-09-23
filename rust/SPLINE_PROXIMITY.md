# Complete point-to-spline minimum distances

The `proximity::closest_points_on_spline*` and
`closest_points_on_exact_spline*` APIs minimize squared Euclidean distance on a
closed rational B-spline parameter range. The supported curve contract remains
degree 1..25, positive weights and continuous positions, including unclamped
active domains and periodic ranges spanning multiple turns. Rendering and C++
runtime dependencies are absent.

## Result and failure contract

`SplineClosestSet` contains **every** isolated minimizing parameter and every
maximal whole minimizing interval. Shared knot endpoints appear once; an
interval includes its endpoints, which are omitted from the isolated list.
Different periodic parameters remain different even when their points coincide.
A singleton query gives one point. Exact roots and curve equations remain the
result's identity; binary64 bounds are optional views and never decide equality.

The API provides exact parameter, coordinate, squared-distance and inter-result
distance comparisons. Finite bounds are tight (one exact float or adjacent
floats); an unrepresentable view returns `Unrepresentable` without destroying
the exact result. Invalid/nonfinite queries, malformed rationals and reversed or
out-of-domain ranges return errors. Negative distance thresholds compare below
every result. No geometric snapping or approximate tie threshold is applied.

`SplineProximityOptions` limits span traversal, candidate count, shared root
isolation subdivisions, image coefficient updates and intermediate image integer size.
Defaults are 4,096 spans, 65,536 candidates, 65,536 subdivisions, 2,000,000 image
updates and 65,536 image bits. Exhaustion returns `ComputationLimit` and no partial
minimum set. These structural limits are **not hard CPU or RSS deadlines**;
BigInt cost, coefficient construction and view conversion also matter.

## Exact construction

On each nonempty normalized knot span, write the homogeneous curve as
`(X,Y,Z,W)` with `W>0`, the query as `Q`, and
`N = Σ(X_i - Q_i W)^2`. Then squared distance is `N/W²`, and its stationary
equation is `F = N'W - 2NW'`, of degree at most `3p-2` (73 at degree 25).

First, exact coefficient vectors determine the curve span's affine hull. Let
`c` be the orthogonal projection of the query onto that hull, expressed relative
to the query, and let `R_i = X_i-Q_i W-c_i W`. Exact orthogonality gives
`N/W² = ||c||² + ΣR_i²/W²`. Thus `||c||²` is a certified lower bound. Every real
common root of the residual polynomials attains it, and these are all local
minima when at least one lies in the closed query range. Their common polynomial
has degree at most `p`; the solver can omit the larger stationary equation and
retain the certified rational distance even at irrational parameters. The hull
may be a point, line, plane or all of 3D; no coordinate alignment is assumed.

An unattained hull bound never becomes a result. In that case, every real
stationary root, nonsmooth knot and range endpoint is considered. `F=0`
identically gives a constant-distance interval. The equation never divides by
curve speed, so singular stationary points are retained. A positive local bound
does not prune later spans, which may contain smaller minima. After an attained
global zero is found, only further zero minima can survive. Shared rational
knot hits merge exactly, while periodic aliases and whole intervals survive.
Both isolation paths use the same subdivision budget.

Exact interval bounds filter candidate distances. Zero distance uses the
sum-of-squares identity. When refined distance intervals still overlap, a
continued-fraction candidate proposes a rational squared distance. The candidate
is accepted only if `N-cW²` vanishes exactly at a selected parameter; an exact
threshold comparison then orders the other distance. This also resolves rational
distance ties at irrational parameters without constructing their images.
A failed proposal retains the general exact fallback, including arbitrarily
close unequal values. Unresolved comparisons use a rational-function image
in the square-free parameter polynomial's quotient ring. Pole factors shared
with the denominator are removed; they cannot contain the selected parameter.
Extended Euclid gives the denominator's inverse. Successive powers of the image
and primitive integer row elimination produce its first linear dependence,
an annihilating polynomial. Common-denominator scales and complete dependence
witnesses are retained through every operation. Exact threshold signs identify
the correct image root, and exact root ordering resolves equality. Images are
reused within a query for candidates sharing a span's defining equation.

Output conversion reuses a bounded, refined isolator across exact threshold
queries. A certified zero returns its zero enclosure immediately. Neither
optimization weakens the exact fallback or replaces it with a tolerance.

The sustained checker includes positive known minima of
`h² + scale² P(t)²(1+t^8)/W(t)²`. Positive weights and the independently known
factors of `P` prove the complete minimum parameter set. Rational parameter
ties reach degree 25; positive-distance irrational ties and composed edits reach
degree eight before optional elevation to nine. These planar offset cases
exercise the hull shortcut, including the unchanged degree-24 timeout seed.
Another family uses the nonplanar rational curve
`C-Q = (scale*s(s²-a), h*(1-s²)/(1+s²), 2h*s/(1+s²))` for `0<a<1`.
Its distance is `h² + scale²*s²*(s²-a)²`, so the complete minimum set is
`{-sqrt(a),0,sqrt(a)}`. Its affine hull is all of 3D and its zero lower bound is
unattainable, retaining general stationary-solver coverage. Positive degree-five
weights, parameter changes, coordinate permutations, translations, knot
insertion, optional elevation to six and closed trims are varied. Expected
roots and coordinates use independent rational radical signs. Separate regressions retain
irrational minimum distances and image-budget failures, and distinguish
parabola minima under a `2^-180` query displacement.

## Reference and independent evidence

OCCT reference revision: `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Both legacy `GeomAPI_ProjectPointOnCurve`/`Extrema_ExtPC` and the newer
`ExtremaPC_Curve` family were read and captured before the Rust implementation.
The committed [original captures](fixtures/occt-spline-proximity-preimplementation)
retain both 29-case native families plus the separate degree-73 legacy case,
source/input/observation hashes and loaded library provenance. Ten input shapes
come from original extrema GTests; four Béziers are represented as equivalent
clamped B-splines. **Their original assertions are not executed.** This adds no
passing DRAW scripts to the existing unchanged-test bridge.

The independent Python oracle derives complete homogeneous Cox coefficients,
factors stationary equations over QQ, isolates roots with VAS and compares
resultant images. Certified interval separation avoids a costly image when a
unique minimum is already proved. Thirty cases contain 66 isolated minima and
three minimizing intervals. Tests check exact parameter identities against
independent irreducible equations and isolators, every curve coefficient, and
every reported parameter/coordinate/distance enclosure. A separate 94-case
SymPy resultant oracle checks every coefficient of the production image
construction, including relative scales, reducible equations and complex poles.

The native bridge uses a source-pinned, headless `TKGeomAlgo` SDK with exception
checks enabled. It validates source/library hashes and actual loaded paths.
Legacy observations include the reported trimmed endpoints; the newer family
uses `PerformWithEndpoints` and retains both variants' raw output. Each exact
minimum needs a distinct compatible native witness under fixed
`1e-6 + 1e-10*abs(expected)` numerical budgets. Extra stationary candidates are
permitted because the native APIs return extrema, not complete global-minimum
sets. Infinite-solution flags cannot represent arbitrary partial intervals.

The local macOS capture has 50 matching family/case pairs and ten separately
[reviewed differences](fixtures/occt-spline-proximity-divergences.json): missing
whole-interval representations, omitted periodic aliases and degree-25 witness
accuracy. The first 29 observations reproduce the pre-implementation stdout
byte-for-byte for both native families. Native exceptions, crashes, timeouts,
malformed output, unreviewed fingerprints and any Rust mathematical disagreement
fail the gate. Native differences never relax Rust's expected result.

## Current verification and remaining scope

The initial local development campaign completed 600.04 seconds of mutation
after 94.87 seconds of corpus replay, with 562 mutations, a 1,405 MiB reported
RSS peak and no artifacts. This was a dirty development tree, not clean-revision
release acceptance. The retained corpus grew from 73 to 311 files. The new
`spline_proximity` target runs with the unchanged 20-second input/2 GiB limits
and joins the daily retained-corpus workflow.

At `e70225a0`, all seventeen CI fuzz targets passed their full 60-second mutation
budgets after replay. The clean local spline campaign completed 600.03 seconds
of mutation after 455.47 seconds of replay, with 584 mutations and a 1,453 MiB
RSS peak. There were no crash, timeout, OOM or mathematical-disagreement
artifacts. Two local slow inputs and one CI slow input were retained; the CI
input reached the 20-second threshold. The subsequent exact zero-set path
reduced that input's local sanitizer callback from about 9.3 to 2.8 seconds in
an isolated experiment. These timings are observations under concurrent work,
not a platform-independent speed guarantee.

The optimized implementation at `fc022e82` completed another clean local run:
600.04 seconds of mutation after 227.91 seconds of replay, 1,908 mutations,
and an 898 MiB RSS peak. No new artifacts were produced; the two previously
retained slow inputs remain archived and both pass the optimized checker.
All seventeen CI fuzz targets passed on this revision. Its Linux spline campaign
completed its full 60-second mutation budget with no artifacts and a 625 MiB RSS
peak. Linux, macOS and Windows each passed 132 debug and 132 release tests;
Rust 1.85 passed all 132 tests, and the WebAssembly library check passed.

The first probe retained 120-second failures for two high-degree bound views.
The zero and shared-refinement fixes now complete all 30 cases locally; the
hard rational degree-73 case takes about 24 seconds including all tight bounds.
These are workload observations, not general latency guarantees.

The first Linux native run independently certified all thirty Rust results,
then exposed an incorrect bridge dependency check: `ExtremaPC_Curve` uses
`TKGeomBase`, while the legacy wrapper also uses `TKGeomAlgo`. Linux correctly
omits that unused library from the newer probe. The corrected bridge requires
the actual family dependencies and still validates every loaded OCCT path and
hash. Separate Linux degree-25 output fingerprints are reviewed with the same
eight-of-25 legacy and ten-of-25 newer-family witness coverage as macOS. Twelve
stored fingerprints cover ten family/case differences across the two platforms;
neither comparison budgets nor Rust expectations were relaxed.

At `14fb4ea3`, the fresh Linux native gate passed with 50 matches, ten reviewed
differences and zero failures, and all seventeen CI fuzz targets passed again.
Only the review data and this document changed from `fc022e82`; its kernel,
tests, tools and dependencies are identical to the fully tested platform source.
The earlier workflow remains recorded as failed because its new Linux native
fingerprint had not yet been reviewed. The final rerun also passed all eight
jobs, including 132 debug and release tests on each of Linux, macOS and Windows,
132 tests on Rust 1.85, and the WebAssembly check.

The subsequent rational-distance filter passed 134 release tests, nine focused
debug tests, both lint checks, ten bridge self-tests, sixteen fuzz-runner
self-tests, and another fresh native comparison with the same 50/10/0 result.
The expanded checker completed 600.05 seconds of mutation after 80.79 seconds
of replay in a scratch development tree: 1,062 mutations, a 1,502 MiB RSS peak,
and no crash, timeout, OOM or mathematical disagreement. One 10-second slow seed
was retained. Clean-revision acceptance is separate from that experiment and
from the accepted zero-distance optimization above.

At `31f378b0`, all eight kernel CI jobs passed: 134 debug and release tests on
each native platform, 134 tests on Rust 1.85, WebAssembly, the original-test
bridge and the source-pinned native comparisons. The clean local spline run
completed 600.04 seconds of mutation after 810.18 seconds of replay, with 1,238
mutations and a 1,674 MiB RSS peak. Sixteen CI fuzz targets passed, but the spline
target timed out during replay of its deterministic degree-24 positive-distance
seed. That failure remains recorded; a local pass does not excuse it.

The affine-hull optimization preserves that input and all comparison, timeout
and RSS limits. In paired local sanitizer replays, the Linux timeout seed went
from 8.51 to 0.67 seconds; three other retained slow inputs also improved.
These timings are local observations, not portable latency guarantees. The
expanded development checker passed 39 separate sanitizer replays, including
35 nonplanar seeds or retained mutations and the four retained slow/timeout
inputs. Development checks passed 136 release tests, eleven focused debug
tests, both lint checks, and a fresh 50/10/0 native comparison. Clean-revision
fuzzing and fresh Linux CI are still required to close the published failure.

This capability
does not implement curve/curve or curve/surface minimum distance, general
intersections, Boolean topology, persistent history, STEP, cancellation or hard
per-operation resource ceilings. Full production kernel parity remains open.

```sh
target/math-oracle-venv/bin/python rust/tools/generate_algebraic_image_fixtures.py --check
target/math-oracle-venv/bin/python rust/tools/generate_spline_proximity_fixtures.py --check
cargo +stable test --locked --release --lib proximity::spline::tests
cargo +stable test --locked --release --test spline_proximity
target/math-oracle-venv/bin/python rust/tools/test_spline_proximity_oracle.py
python3 rust/tools/build_pinned_occt.py --output target/spline-proximity-sdk --toolkit TKGeomAlgo --jobs 2
target/math-oracle-venv/bin/python rust/tools/compare_spline_proximity.py --occt-root target/spline-proximity-sdk/install --sdk-manifest target/spline-proximity-sdk/build-manifest.json
python3 rust/tools/run_fuzz.py --target spline_proximity --seconds 600
```
