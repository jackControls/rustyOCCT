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
Every real stationary root, every nonsmooth knot and both range endpoints are
considered. `F=0` identically gives a constant-distance interval. The equation
never divides by curve speed, so singular stationary points are retained.

First, the common polynomial of `X_i-Q_i W` is checked for every exact
zero-distance parameter on the span. It has degree at most `p`. If it has real
roots in the queried range, they are all local minima with distance zero; other
stationary points cannot improve them. Three identically zero deltas give a
whole zero-distance interval. After any zero is found, subsequent spans need
only this complete zero-set check. Earlier positive-distance winners are
discarded. When there are no real zeros, the stationary equation remains the
fallback. Shared rational knot hits merge exactly, while periodic aliases and
whole intervals survive. Both isolation paths use the same subdivision budget.

Exact interval bounds filter candidate distances. Zero distance uses the
sum-of-squares identity. Unresolved comparisons use a rational-function image
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

The first probe retained 120-second failures for two high-degree bound views.
The zero and shared-refinement fixes now complete all 30 cases locally; the
hard rational degree-73 case takes about 24 seconds including all tight bounds.
These are workload observations, not general latency guarantees.

The first Linux native run independently certified all thirty Rust results,
then exposed an incorrect bridge dependency check: `ExtremaPC_Curve` uses
`TKGeomBase`, while the legacy wrapper also uses `TKGeomAlgo`. Linux correctly
omits that unused library from the newer probe. The corrected bridge requires
the actual family dependencies and still validates every loaded OCCT path and
hash. A separate Linux degree-25 output fingerprint is reviewed with the same
eight-of-25 witness coverage; neither comparison budgets nor Rust expectations
were relaxed.

Clean-revision native/platform/fuzz acceptance of these fixes is still pending. This capability
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
