# Exact B-spline knot editing contract

This milestone changes curve representation without changing its homogeneous
piecewise polynomial. It does not implement general degree elevation, surface
knot editing, approximation, topology, or tolerance-based simplification.
The `exact_spline_plane`, `exact_spline_sphere` and `exact_spline_cylinder`
intersection families accept edited curves directly, including rational query
intervals. [Their exact results](EXACT_SPLINE_INTERSECTIONS.md) retain contacts
outside binary64 range and provide optional finite enclosures.

## Representation and operations

`ExactKnotVector` accepts rational, strictly increasing distinct knots, degree
1 through 25, and at most 4096 poles, with more poles than the degree. Its
multiplicity and periodic pole-order conventions match `KnotVector`. Rational
inputs are normalized; a zero denominator is an error. `ExactBSplineCurve3`
stores positive homogeneous controls `(w*x, w*y, w*z, w)`. Converting an existing
binary64 curve preserves its input values exactly. Nothing is rounded back to
binary64 between operations.

Refinement raises the **total** multiplicity at each requested knot. Requests
may be unordered; duplicates and the two periodic seam aliases use the maximum
target. Zero and targets below the existing multiplicity are no-ops. Targets
above the allowed multiplicity, knots outside the closed fundamental domain,
invalid rationals, more than 4096 requests, and output exceeding 4096 poles are
errors, even when another request would be a no-op. No epsilon-based snapping
is performed. Existing unclamped exterior knots and inactive end controls are
preserved. Nonperiodic physical end knots can have multiplicity degree+1;
other knots and periodic seam endpoints can have at most degree.

Removal requests a lower total multiplicity at an existing knot. A target at
or above the current multiplicity returns an unchanged curve. Nonperiodic
removal is restricted to knots strictly inside the parameter domain. Periodic
endpoints identify one seam. Removing its last occurrence shifts the origin to
the next distinct knot and retains the same period, following OCCT. Requests
that leave the supported representation family are errors. A valid request
returns no curve if exact inverse insertion cannot preserve all four
homogeneous components, or if the result would contain a nonpositive weight.
This stronger contract can reject removals which OCCT accepts within a tolerance;
it does not claim to find every equivalent Euclidean rational representation.
All operations are immutable and atomic.

Evaluation returns exact position and derivatives through order two in original
parameter units, with the existing explicit side/continuity rules. Optional
binary64 enclosures can fail when the exact answer is outside binary64 range.
Extraction accepts rational intervals and counts periodic spans before allocation.
The existing default limit of 4096 Bézier arcs applies.

## Evidence and reference algorithms

Source baseline: `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Read `Geom_BSplineCurve::{InsertKnots,RemoveKnot}` and
`BSplCLib::{PrepareInsertKnots,InsertKnots,RemoveKnot,BoorScheme,AntiBoorScheme}`.
Refinement uses their local homogeneous insertion relation with exact rational
arithmetic. Removal inverts that relation and checks the opposite anchor by
exact equality. Periodic curves are temporarily unrolled into a bounded five
period window; the new extended knot sequence selects the canonical pole block.

The initial 665 native captures precede the Rust implementation; two additional
nonconstant seam cases were added during validation. They preserve every operation's
success flag, complete knot/multiplicity arrays, control points and weights,
degree, periodicity, and domain. The installed OCCT SDK is an observation oracle,
not necessarily a runtime build of the source baseline. Tests include adaptations
of `Geom_BSplineCurve_Test.cxx` insertion/removal cases, all supported degrees,
unclamped inactive end controls, repeated knots, periodic seam aliases, origin
changes, round trips and failed removal. These are not unchanged upstream GTest
or DRAW passes.

Independent exact Cox basis power coefficients check complete homogeneous
polynomial identity, rather than only sampled points. Exact fixtures also cover
rational cuts, sub-ULP separations, extreme scales, invalid inputs, operation
sequences and resource limits. Fuzzing must retain this identity check and
removal/reinsertion assertions; native output alone is not proof of correctness.

Count limits bound allocations and iteration counts, not the cost of arbitrary
precision arithmetic. Arbitrarily large rational inputs have no wall-clock or
bit-length guarantee. The sustained fuzz campaign uses its existing process
timeout and memory limits.
