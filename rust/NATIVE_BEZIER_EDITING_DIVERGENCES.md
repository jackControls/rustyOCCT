# Native Bézier extraction and editing observations

Source reference: `3d097a0328e71b826377d4814ab05ec3c3d23871`. All 546
native inputs were captured before the Rust implementation. They exercise
full/clipped clamped, unclamped and periodic B-spline extraction, Bézier
trimming, splitting, reversal, elevation and combined edits. Source inputs from
`Geom_BezierCurve_Test` include the polynomial cubic and rational increase/reverse
examples. The additional quarter-conic case uses `sqrt(0.5)`, one binary64 ULP
above the source RationalSegment test's `1/sqrt(2)` weight; it is not a verbatim
input replay of that test. This is an API bridge, not execution of the unchanged
GoogleTest file or an additional unchanged DRAW pass.

Rust keeps original parameter ranges; OCCT's individual Bézier arcs use local
`[0,1]`. The bridge applies that explicit affine map. It compares every pole and
weight, degree, arc count and domain. Native poles are Euclidean and weights
are normalized by the first weight to allow a common homogeneous scale. The
budget is `1e-10 + 2e-12 * abs(exact_component)`. Every Rust homogeneous control
must first equal the exact independent basis-polynomial result without any
normalization or numerical budget.

## OCCT 7.9.3 on macOS

- 546 independently verified Rust cases, producing 1,672 exact edited arcs.
- 506 native cases match the unchanged numerical budget.
- 40 cases differ: 13 trims, 14 splits, 13 combined trim/elevate/reverse/splits.
  All involve degree 25. Every native count, degree and domain agrees.
- All extraction-only, reversal-only and elevation-only cases match, including
  degree 25. All lower-degree edit cases match.

For `d25_r0_s1_w0_op2`, the unweighted degree-25 Bézier curve's last input pole
is `(22,4,0)`. Splitting at 3/8 must leave the right arc's endpoint exactly there.
OCCT reports z = `4.0531125140574886e-05`; Rust retains z = 0. The same observation
occurs on the corresponding first span of `d25_r0_s3_w0_op2`. This endpoint
check is independent of both implementations' interior evaluations. It also
demonstrates that these discrepancies affect geometry, not just a harmless
common rescaling of homogeneous controls.

The read source route is `Geom_BezierCurve::Segment` → `BSplCLib::BuildCache`
→ `PLib::Trimming` → `PLib::CoefficientsPoles`. It converts to power coefficients,
substitutes the parameter interval and converts back using floating arithmetic.
The observed errors are consistent with cancellation in that degree-25 route;
the individual floating instruction responsible has not been isolated. Rust
instead uses exact positive homogeneous de Casteljau interpolation. Independent
Python and Rust basis-polynomial checks preserve the complete parameterized
curve, including endpoints and exact derivatives.

## OCCT 7.6.3 on Linux

Linux's OCCT 7.6.3 reports 508 matching cases and 38 differences: 12 trims,
14 splits and 12 combined edits. These are a subset of the same macOS degree-25
cases; `periodic_d25_m1_w0_op1` and `periodic_d25_m1_w0_op5` fit the budget on
Linux. Native values differ, so the 38 observations have separate full bit/field
pins. The unweighted split example's Linux endpoint is
`(22,4.000002201702663,5.0594034521456166e-05)`. All 546 complete exact Rust
outputs are byte-identical across the platforms, and all native count/degree/
domain comparisons agree. The registry contains 78 version-specific records.

## Review guards

`fixtures/occt-bezier-editing-divergences.json` pins runtime version, input hash,
every native control/domain bit, exact-output hash and differing field indices.
A changed observation, input or unreviewed runtime fails; `--strict-native`
rejects reviewed differences too. Reviews never bypass the exact Rust check.
Corruption tests cover changed interior controls with intact endpoints, missing
arcs, reversed order, incorrect ranges, malformed values and every review pin.

Additional exact-only cases cover full-exponent data, subnormal and adjacent
knots, multiple periodic turns and distinct rational cuts with indistinguishable
floating values. They do not claim native agreement outside the native corpus.
