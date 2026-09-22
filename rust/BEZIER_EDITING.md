# Exact Bézier extraction and editing contract

This capability supplies geometry-preserving curve operations for later edge
splitting and modeling. It does not build B-rep edges or topology/history maps,
edit surfaces, remove knots, join curves, or implement arbitrary B-spline
knot-vector editing. Source reference: `3d097a0328e71b826377d4814ab05ec3c3d23871`.

`BSplineCurve3` can be decomposed into exact rational Bézier arcs over its
fundamental domain or an explicit positive-length parameter interval. All
existing degree 1..25, positive-weight, clamped/unclamped and periodic inputs
are supported. Periodic queries may cross seams or span multiple periods,
subject to an explicit preflight arc-count limit. Nonperiodic extrapolation is
an error. Knots are never snapped together.

Every arc retains positive homogeneous rational control points `(wx,wy,wz,w)`
and an exact increasing parameter interval `[a,b]`. The polynomial is expressed
in local `t=(u-a)/(b-a)`, but public evaluation and derivatives use the original
parameter `u`. Floating bounds are requested separately. Exact evaluations and
control data survive unrepresentable derivatives or other output conversions.
Each arc supplies its own endpoint limits; adjacent arcs can have different
one-sided derivatives at a knot. No claim of differentiability is inferred by
splitting the curve.

Immutable operations on an arc:

- Subdivide at a strictly interior rational parameter, returning two closed
  arcs with exactly the same shared endpoint and the original parameter units.
- Trim to a positive-length closed subinterval, including the full interval.
  No extrapolation or zero-length curve is silently substituted.
- Reverse while retaining the same increasing domain: the new value at `u`
  equals the old value at `a+b-u`. First derivatives change sign.
- Elevate to a requested degree between the current degree and 25, preserving
  the complete parameterized rational function. Equal degree is an identity;
  degree reduction is an error.

Operations do not round their control data back into `BSplineCurve3`'s f64 pole
storage. Rational parameters permit repeated editing without a float round trip.
Input rational parameters are normalized; a zero denominator is rejected.
Resource limits bound arc enumeration and degree, not a hard time deadline for
arbitrarily large user-supplied rational parameters or edit sequences.

Read entry points and called algorithms together:
`GeomConvert_BSplineCurveToBezierCurve` copies, segments and raises interior
multiplicities; `Geom_BSplineCurve::Segment` clips knots/control data and handles
periodic origin changes; `BSplCLib::InsertKnots` performs control refinement.
`Geom_BezierCurve::{Segment,Reverse,Increase}` reference `PLib::Trimming`
(polynomial parameter substitution) and `BSplCLib::IncreaseDegree` (degree
increment averaging). Rust uses exact blossom evaluation, de Casteljau splits
and Bernstein elevation. It preserves the supported geometry rather than the
native floating tolerance policy. In particular, Rust's original parameter
ranges differ from OCCT Bézier's reset local `[0,1]` range; the bridge accounts
for that map, and multi-period extraction exceeds native Segment's one-period
limit without pretending it is a matching native query.

The 546 native extraction/editing inputs were captured before Rust implementation,
including upstream `Geom_BezierCurve_Test` polynomial and rational edit inputs.
The bridge retains every native output, independently verifies complete exact
homogeneous polynomial identities, and distinguishes numerical mismatches from
agreement. The 636 exact fixtures extend this corpus to extreme inputs; an
independent coefficient oracle also checks sustained fuzz campaigns. These API
observations do not count as executing unchanged GoogleTest or DRAW files.
