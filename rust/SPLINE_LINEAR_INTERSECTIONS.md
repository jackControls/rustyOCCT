# Complete spline/line and spline/segment intersections

The `intersection::spline_line*`, `spline_segment*`, `exact_spline_line*` and
`exact_spline_segment*` APIs return the complete preimage of an infinite line
or closed segment on a closed rational B-spline parameter range. Binary64 and
`ExactBSplineCurve3` inputs share one exact engine; binary64 inputs are
converted exactly before use. The supported curve contract is unchanged:
degree 1..25, positive weights, continuous positions, unclamped active domains
and periodic ranges spanning multiple turns. No rendering or C++ runtime
dependency is involved.

## Result and failure contract

`SplineLinearIntersection` contains **every** isolated curve parameter and every
maximal closed parameter interval whose points lie on the query. Interval
endpoints are omitted from the isolated list. Adjacent spans merge only at an
exactly equal shared knot. Different parameters remain different when their
points coincide: periodic aliases, backtracking curves and constant spans all
keep their separate parameter sets. A constant curve can therefore produce a
whole interval whose geometric image is one point. A singleton range gives at
most one point.

The line or segment is `A + t(B-A)`, in the caller's endpoint order. Lines
require `A != B` and return `Degenerate` otherwise. A collapsed segment is the
point `A`, with deterministic linear parameter zero. Reversing a segment keeps
every curve parameter and maps `t` to `1-t`.

`SplineLinearPoint` retains an algebraic curve parameter. It supports exact
`compare_parameter`, `compare_coordinate`, `compare_linear_parameter` and
cross-result `parameter_cmp`, including results on different spans, curves,
queries or periodic turns. None of these uses a binary64 enclosure to decide
equality. `parameter_bounds`, `coordinate_bounds`, `point_bounds` and
`linear_parameter_bounds` are optional views. They are tight when finite (one
exact float or two adjacent floats); otherwise they return `Unrepresentable`
and leave the exact result usable. Overlap endpoints are ordinary points.

Nonfinite binary64 inputs, malformed rationals, reversed ranges and
nonperiodic ranges outside the active domain are errors.
`SplineLinearOptions` bounds span traversal (default 4,096), root-isolation
subdivisions (65,536) and retained boundary candidates (65,536). Exhaustion
returns `ComputationLimit` and no partial set. These are structural limits,
**not hard CPU or RSS deadlines**: BigInt size and later view conversions also
contribute to cost.

## Exact construction

On each nonempty normalized span, write the homogeneous curve as `H=(X,Y,Z,W)`
with `W>0`, and the query direction as `D=B-A`. Let `Δ = (X,Y,Z) - A W`.

* **Line or nondegenerate segment.** The curve point lies on the line exactly
  when `Δ × D = 0`. The primitive integer gcd of the three cross-product
  components has degree at most the curve degree. Its real roots in the span
  are the isolated candidates. If all three components vanish identically,
  the whole span lies on the line.
* **Collapsed segment.** The three point-preimage equations `Δ = 0` are used
  directly instead of dividing by `D·D = 0`.
* **Segment clipping.** The linear parameter is `L/((D·D)W)` with
  `L = Δ·D`. Every cleared denominator is positive, so membership is
  `L >= 0` and `(D·D)W - L >= 0`. An isolated candidate is kept when both
  signs, evaluated exactly at its algebraic root, are nonnegative.

A contained span needs a complete one-variable sign decomposition. The engine
isolates every root of both boundary polynomials and orders them exactly. The
sign immediately to the right of a boundary is the sign of its first nonzero
derivative: the Taylor factorial and positive displacement powers do not
change it. Since every boundary root is isolated, that sign holds on the whole
next open cell. No approximate midpoint, repeated separation or binary64
decision is needed. Closed inequalities imply that an accepted cell contains
both endpoints. Adjacent accepted cells merge; an accepted boundary with no
accepted neighboring cell is an isolated contact. An identically zero boundary
polynomial is satisfied everywhere rather than being treated as a finite root
list.

Results from different spans are ordered in original parameter units by an
exact comparison of positive affine images `a + b·α` and `c + d·β`. First, `α`
is compared with the rational images of `β`'s isolator endpoints; a rational
`β` needs one such comparison. Otherwise `α` maps strictly inside `β`'s
unique-root isolator. The sign at `α` of `d^n p((a-c+bt)/d)`, where `p` is `β`'s
square-free defining polynomial, then decides the order or proves equality.
This reuses the existing algebraic sign engine and avoids constructing an
image resultant.

## Shared arithmetic changes

Two changes to the shared exact root arithmetic were driven by this
capability's fuzz inputs. Neither changes a mathematical result.

* Integer content during primitive polynomial remainders now uses Stein's
  binary gcd from `num-integer`, which `num-bigint` and `num-rational` already
  require. The previous schoolbook Euclid performed a full multiprecision
  division at every step.
* Rational-root recognition now also tries the unique candidate given by the
  rational root theorem: a root `a/b` in lowest terms has `b | lead(p)`, so
  once the isolator is narrower than `1/|lead(p)|`, `ceil(lead·lower)/lead` is
  its only rational candidate. As before, a candidate is accepted only after
  an exact zero check inside the unique-root isolator. The continued-fraction
  candidate needs width about `1/b²`; this one needs `1/|lead|`.
* Before an exact sign query falls back to a gcd or Sturm-Tarski chain, the
  isolator is bisected up to 64 further steps if that reaches width
  `1/|lead|`. Bisection evaluates only the small defining polynomial, while the
  fallback gcd can involve a much larger query polynomial. An uncapped version
  slowed spline/plane tests by 46%; with the cap, every suite stays within
  noise of the baseline.

On a degree-25 known-factor input with a sub-float root pair near `1/3`,
queries at the large-denominator root previously built a multi-thousand-bit
integer gcd each time. Its complete release fuzz-oracle replay fell from 65.5
seconds to 7.6 seconds with the binary gcd, then to 0.79 seconds with both
changes.
The first clean campaign then found a degree-25 query exceeding the 20-second
sanitizer limit. The capped lead-bound bisection cut its release check from
4.48 to 1.25 seconds (11.7 seconds under AddressSanitizer). The complete release
kernel suite shows no suite-level regression, and kernel unit tests fell from
23.5 to 12.1 seconds on the development machine.

## Independent evidence

* `generate_spline_linear_fixtures.py --check` reproduces 60 complete
  preimages: 92 isolated parameters and 20 maximal intervals. Its answers come
  from independent Cox basis power coefficients, a QQ gcd, irreducible
  factorization and VAS isolation. Contained spans are clipped at exact
  rational **sample points** between separated boundary roots, not by
  production's derivative signs. Each case records every parameter's
  irreducible equation, unique-root isolator and a rigorous linear-parameter
  bracket. Cargo proves exact membership in that equation and isolator, and
  checks the linear parameter. The first 31 cases are the native preflight
  inputs, and each of their counts matches that separate SymPy-inequality
  preimage. The remaining 29 add rational trims, singletons, collapsed
  segments on and off curves, periodic aliases, constant curves, oblique 3D
  secants, irrational contacts, parameters beyond binary64 range, `2^-1100`
  coordinates, `2^-180` gaps, contact order six, unclamped knots, periodic
  quadratic turns and weighted retracing. The 43 binary64-representable cases
  also run through the binary64 entry points. Every case is re-queried after
  exact knot insertion and degree elevation.
* `generate_closed_clip_fixtures.py --check` checks the clipping primitive on
  79 SymPy inequality sets (79 isolated points, 69 maximal intervals). These
  have equations through degree 25, gaps through `2^-600` and bounds through
  `2^2048`.
* `generate_affine_parameter_fixtures.py --check` checks 851 cross-equation
  affine comparisons in both directions (79 equal). They come from 69 equation
  pairs with degree-25 clusters, `2^-600` separations and offsets and scales
  through `2^2048`.
* Nine public tests cover secants, tangencies, backtracking with irrational
  interval endpoints, constant curves, periodic aliases, exact edits,
  sub-float gaps, unrepresentable views and atomic failure.

The eighteenth fuzz target, `spline_linear`, constructs complete answers
algebraically. It covers four families:

1. Known-factor rational Bézier curves through degree 25, including sub-float
   root pairs, with rational segment clipping.
2. Rational polylines with rational interval endpoints, periodic multi-turn
   ranges, constant spans, retracing and collapsed segments.
3. A rational quarter circle against rational-slope lines, whose contacts are
   irrational and are described by exact sign functions.
4. A quadratic retrace whose overlap endpoints are `(1 ± √g)/2`.

Every input also applies an invertible integer affine map to the curve and
query, a common homogeneous scale, parameter domains through `2^±2048`,
optional exact knot insertion and elevation to degree 25, and query reversal.
Algebraic values are bracketed to `2^-96` in construction units by
independent bisection. The kernel must place each result strictly inside its
bracket and produce tight, or typed unrepresentable, views.

A local development campaign ran on a dirty tree based on `d1206b15`, with
AddressSanitizer and the standard 20-second/2 GiB limits. It completed 300.01
seconds of mutation after 36.94 seconds of seed replay, with 1,011 mutation
executions, 6,185 coverage edges and a 1,603 MB RSS peak. It produced no crash,
timeout, OOM or disagreement artifact, and the corpus grew from 100 to 388
inputs. The two formerly slow random inputs are retained as regressions. This
is not clean-revision acceptance. Its RSS headroom should be watched in longer
daily campaigns.

Clean local 600-second campaigns followed. At `da7be8be`, one input exceeded
the 20-second sanitizer limit; it is now a retained regression. At
`809f8b4b`, accumulated allocator retention caused an OOM at 2,055 MB. The
`d6b2105c` run used the allocator settings described in `FUZZING.md` and
passed. It completed 600.07 seconds of mutation after 291.76 seconds of
replay, with 1,488 mutation executions, 6,288 coverage edges and a 1,035 MB
RSS peak. There were no crash, timeout, OOM or disagreement artifacts; one
slow-unit input was reported. Linux CI and platform acceptance are still
pending.

## Native OCCT observations

Before any Rust implementation existed, the source review read
`IntTools_EdgeEdge::{Perform,Prepare,FindSolutions,MergeSolutions,AddSolution,
FindBestSolution,ComputeLineLine,IsIntersection,IsCoincident}` and
`IntTools_CommonPrt` at `3d097a0328e71b826377d4814ab05ec3c3d23871`. A
source-pinned, headless TKBO SDK captured the first 31 fixture inputs with
explicit curve ranges, zero extra fuzzy value and quick coincidence checks
disabled. Thirty completed. The degree-25 Chebyshev case exceeded the
20-second process deadline and remains a recorded failure, not a match. Its
32-window dyadic partition finished in 17.7 seconds and reported eight distinct
verified witnesses of 25 exact roots.

These pre-implementation observations show contract differences, not
necessarily native defects:

* `IntTools_EdgeEdge` is a finite-edge algorithm, so both unbounded lines
  returned empty results. A finite pole-hull cover recovered the expected
  result types.
* `BRepLib_MakeEdge` adjusts periodic ranges to the fundamental period. The
  multi-turn query was therefore outside the constructed edge. Querying each
  span inside that range and retaining integer period offsets reproduced the
  exact multi-turn preimage for that fixture.
* `MergeSolutions` classifies a result as an edge when either range is fully
  covered, clears earlier common parts and stops. It is not designed to list
  every disjoint parameter preimage of a retracing curve.

## Native comparison bridge

`compare_spline_linear.py` first verifies the pinned SDK manifest, every loaded
toolkit, the unchanged pre-implementation capture in
`fixtures/occt-spline-linear-preimplementation/`, and byte-identical
regeneration of the native inputs. That capture's library paths are stored
relative to the SDK prefix, with their original hashes. The bridge then
certifies all 31 Rust probe rows against the independent generator. Each
tight parameter bound must contain its irreducible root, and each line
parameter must meet the rigorous bracket. Native observations are compared
only after this, under the existing 1e-6 absolute plus 1e-10 relative budget.

`occt_spline_linear_oracle.cpp` is the pre-implementation probe with two
structural adapters, and neither changes a tolerance or a result:

* An unbounded line uses the finite range `[-2e, 2e]` of the same `gp_Lin`,
  where `e` bounds every pole's L1 distance from `A`. Positive weights keep the
  curve in its pole hull, and every hull point's unit-direction projection has
  magnitude below `e`, so no intersection is excluded.
* A trailing mode can restrict the curve with `Geom_BSplineCurve::Segment`,
  which keeps original parameters.

The driver adds one fundamental-period window per requested periodic turn,
with its exact integer offset. Only after a recorded whole-range timeout does
it run 32 uniform `Segment` windows.

On the pinned macOS SDK (OCCT 8.1.0), 23 cases match and eight are reviewed
contract or tolerance differences, with no failures. The eight are:

* two tolerance contacts on exact misses by `2^-40` and `2^-54`
* two exact crossings `2^-19` apart merged into one vertex
* the second disjoint retrace interval omitted by `MergeSolutions`
* three collinear pieces reported as single vertices
* the degree-25 whole-range timeout, whose segmented observation witnesses
  all 25 exact roots

The timeout is reviewable only together with that complete segmented
observation, and it remains listed in every report.
`occt-spline-linear-divergences.json` fingerprints each observation and gives
its independent evidence. `test_spline_linear_oracle.py` checks deliberately
wrong rows and native classifications.

Linux CI builds the same pinned TKBO SDK. Its first native fingerprints need
their own review, as spline proximity's did. Platform acceptance and a
clean-revision fuzz campaign remain **pending**. General curve/curve and
curve/surface intersections beyond planes, spheres, cylinders, lines and
segments remain separate work.
