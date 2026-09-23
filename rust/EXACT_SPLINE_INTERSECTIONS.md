# Exact edited-curve intersection contract

The plane, sphere and infinite-cylinder intersection engine accepts
`ExactBSplineCurve3`, including arbitrary rational knots and positive homogeneous
controls, without converting the curve or query interval to binary64.

The exact result retains every isolated algebraic contact and maximal closed
overlap interval in the original parameter units. Rational comparisons remain
available even when a parameter or position cannot be represented by a finite
binary64 interval. Enclosures are an explicit fallible conversion. Existing
binary64 curve entry points retain their current return types and failure rules.

Nonperiodic query intervals must be positive and inside the curve domain.
Periodic queries can cross seams and span multiple turns; there is no snapping,
automatic endpoint adjustment or reversal. Contacts at a shared knot are merged
by exact parameter identity. Distinct roots remain distinct when their floating
enclosures coincide. Left/right contact orders, crossing/tangent/boundary
classification and exclusion of overlap endpoints retain the existing contract.

The existing span and root-subdivision limits apply to the complete operation;
exhaustion is an error, never partial success. These limits do not impose a bound
on arbitrary-precision operand size or a hard wall-clock deadline. Malformed
rationals, including zero denominators, are rejected before use.

Reference OCCT `GeomAPI_IntCS::{Perform,Parameters,Segment}`,
`IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`, and
`IntSurf_Quadric::Distance` at `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Their curve parameters, isolated-point/interval result model and implicit
surface equations inform the contract. The Rust engine instead uses exact
homogeneous span polynomials and complete algebraic root isolation.

Before implementation, fresh native captures retained 214 plane and 92
sphere/cylinder cases on the original represented curves. Exact refinement and
removal preserve their parameter-to-point functions; comparisons reuse
those native inputs and report their existing reviewed differences separately.
Arbitrary rational data requires independent exact oracles, because rounding
it into OCCT would change the input. No new unchanged upstream GTest or DRAW
passes are implied. General curve/surface intersections remain separate work.

## Using retained results

`exact_spline_plane`, `exact_spline_sphere` and `exact_spline_cylinder` query the
fundamental domain. Their `_in` variants accept rational trim bounds, and
`_with_options` / `_in_with_options` variants expose the existing work limits.
An exact refinement or removal result can be passed directly to these functions.

The returned `ExactSplineSurfaceIntersection` exposes ordered `points()` and
maximal `overlaps()`. A point supports `compare_parameter(&BigRational)` and
`compare_coordinate(component, &BigRational)` without requesting floating bounds.
Its `parameter_bounds()` and `coordinate_bound(component)` request individual
minimal binary64 enclosures; `coordinate_bounds()` requests all three. Calling
`enclosed()` on a point or the complete result requests the corresponding legacy
bounded representation. An `Unrepresentable` error leaves the original exact
result available. Overlap endpoints are available directly as rational values,
with analogous comparisons and fallible enclosure conversion.

## Independent evidence

The 402 existing spline/surface intersection fixtures check the original curve,
its exact representation, a rationally refined representation and an exact
refinement/removal round trip against the same independently derived complete
answers. The native bridge repeats these four representations for its 306
geometries; this is 306 distinct native cases, not 1,224.

An additional 274 Python `Fraction` fixtures construct complete known-factor
equations, contained arcs and periodic crossings. They check
exact contacts, coordinates, one-sided orders, maximal overlaps and minimal
finite enclosures or typed conversion failure. They include degree 25, quadric
contact order 50, separations of 2^-2048 and domains outside binary64 range.
The dedicated fuzz target mutates these families and exact knot edits, checking
complete polynomial preservation with a separate Cox coefficient oracle.
