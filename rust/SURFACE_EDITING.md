# Exact tensor-product patch editing

The capability contract is defined before implementation and comparison with
the native oracle. Inputs are every currently supported `BSplineSurface3`:
degrees 1–25 in each direction, at most 4096 source poles, finite coordinates,
positive finite weights, clamped or unclamped knots and independent U/V
periodicity. This is rectangular parameter-space geometry, not arbitrary face
trimming, sewing, topology, or a general surface intersector.

Extraction preserves the homogeneous tensor polynomial exactly on each
nonempty knot rectangle. Output controls are reduced rationals `(wx,wy,wz,w)`
in U-major order and retain the original increasing U/V parameter domains.
The fundamental domain and explicit finite positive rectangles are supported;
periodic axes can cross seams and cover multiple turns. Patch order is
increasing U, then increasing V. Both axis traversals and their Cartesian
product are bounded before constructing control grids. The default limits are
4096 patches and 1,048,576 output controls. Caller-supplied larger limits permit
more work; these limits are not deadlines or rational bit-size limits.

Immutable operations are strict-interior U/V subdivision, positive rectangular
restriction, reversal in either axis, exchange of U and V (including domains),
degree elevation up to 25 independently in each axis, and constant-U or
constant-V curves including boundary curves. Reversal retains increasing
domains and maps `u` to `a+b-u` (similarly for V). A constant-U curve varies in
the original V units. Every edit preserves exact homogeneous data, without a
binary64 control-point round trip or projective rescaling. Invalid rational
parameters (zero denominator), domains and degrees fail explicitly. Exact
rational cuts are supported; nonrational algebraic cuts need a future
representation extension.

Exact position and all requested partials through total order two retain
rational values, including the mixed derivative. Optional binary64 enclosures
are minimal and fail atomically if any requested value is unrepresentable.
Each patch owns its one-sided boundary limits; continuity across adjoining
patches is not inferred from this representation.

OCCT references at `3d097a0328e71b826377d4814ab05ec3c3d23871`:
`GeomConvert_BSplineSurfaceToBezierSurface` (segmentation, multiplicity raising,
patch controls and boundaries), `Geom_BSplineSurface::Segment`,
`Geom_BezierSurface::{Segment,Increase,UReverse,VReverse,ExchangeUV,UIso,VIso}`,
`BSplSLib::{Iso,IncreaseDegree,BuildCache}` and
`BSplCLib::{InsertKnots,BoorScheme}` and
`PLib::{UTrimming,VTrimming,CoefficientsPoles}`. Native segmentation has
tolerance snapping and a one-period restriction; Rust's exact domain contract
is deliberately independent of those restrictions. Original GTest input
families come from `Geom_BezierSurface_Test.cxx`; the bridge is an adaptation,
not an additional unchanged upstream test pass.

Validation must compare complete homogeneous tensor polynomials, exact jets,
shared boundaries, independent axis operations and explicit resource errors.
Native controls are compared after Euclidean conversion and a single common
weight normalization. Native exceptions and numerical differences remain
visible; they cannot excuse a Rust mismatch with the independent exact oracle.

There are 744 complete-control fixtures, including exact extremes not accepted
by the native SDK. A shared denominator per output item stores all rational
controls without repeating denominator strings; no values are sampled or
omitted by this lossless representation. Every one of the 691 live native
cases also recomputes every Rust control independently. See
[validation](VALIDATION.md), [native differences](NATIVE_SURFACE_EDITING_DIVERGENCES.md)
and [sustained fuzzing](FUZZING.md).

The second Rust oracle independently forms Cox tensor power coefficients,
uses binomial substitutions, and evaluates direct monomial derivative sums.
It carries an integer grid over one positive denominator across transforms;
complete polynomial equality uses exact cross multiplication, with no sampled
or omitted coefficients. Linux release CI checks this oracle on all 744
fixtures, beyond the representative checks in ordinary debug tests. Reproduce:

```sh
RUSTY_VERIFY_ALL_SURFACE_ORACLES=1 cargo test --locked --release \
  --test surface_editing complete_tensor_controls_and_independent_quotient_jets
```
