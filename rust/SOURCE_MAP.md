# Source-driven implementation

Use the original implementation, its callers and its tests together when
porting a capability. The source baseline is OCCT commit
`3d097a0328e71b826377d4814ab05ec3c3d23871`, retained in this repository.
The runtime oracle version is recorded separately: local OCCT 7.9.3 and the
Ubuntu 24.04 distribution build in CI are not builds of that source commit.

## Current traceability

The table distinguishes behavior ported from source from an independent
analytic implementation checked against OCCT. It does not claim a line-for-line
translation of the current Rust kernel.

| Rust behavior | Original source and contract | Implementation / evidence |
| --- | --- | --- |
| `Tolerance::default` | [Precision.hxx](../src/FoundationClasses/TKernel/Precision/Precision.hxx), `Confusion()` and `Angular()` | Same `1e-7` and `1e-12` defaults. Rust additionally validates supplied tolerances and rejects unresolvable coordinates. |
| `predicates::orient2d`, polygon crossings/classification | [CSLib_Class2d.cxx](../src/FoundationClasses/TKMath/CSLib/CSLib_Class2d.cxx) for ray/boundary behavior; [Shewchuk's predicate filter](https://www.cs.cmu.edu/~quake/robust.html) for numerical error bounds | Exact sign of the represented finite-f64 inputs. Conservative filtered path and independently implemented bounded integer fallback; not an OCCT algorithm port. The OCCT normalization/grid is not reproduced. Independent rational/integer oracles and integration invariants; see `MATHEMATICS.md`. |
| `Solid::box_at` | [BRepPrimAPI_MakeBox.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeBox.cxx), `pmin` and point/signed-size constructor | Negative dimensions move the minimum corner; absolute dimensions define positive volume. All eight sign combinations checked in Rust and through native DRAW. Zero/sub-tolerance sizes remain errors. Existing `cuboid` still takes positive dimensions. |
| `orient3d`, `in_sphere` | [Shewchuk's predicates](https://www.cs.cmu.edu/~quake/robust.html), `orient3d` / `o3derrboundA`; exact determinant mathematics | Independently implemented finite-binary64 integer determinants. Above-plane orientation is opposite Shewchuk's sign. Sphere classification normalizes tetrahedron orientation and rejects coplanar defining points. Python rational matrix and continuous rational fuzz oracles; no OCCT predicate-port claim. |
| `line_plane`, `segment_plane` | [IntAna_IntConicQuad.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_IntConicQuad.cxx), `Perform(gp_Lin,gp_Pln)`, distance/direction substitution and parallel/contained branches | Retain affine intersection semantics with endpoint parameterization. Exact three-point plane and rational construction replace rounded normalized coefficients and angular-tolerance decisions. 72 live native OCCT cases cover the well-conditioned common domain; rational fixtures cover extremes and deliberate exact/tolerance divergence. |
| `segment_triangle` and construction enclosures | Exact projected half-plane membership/clipping and binary64 ordering; independently implemented | Includes coplanar overlap, tangency and degenerate segments. Every coordinate and parameter has a minimal finite enclosure; unrepresentable infinite-line constructions fail explicitly. Python and fuzz rational barycentric oracles; no claim to have ported OCCT's general intersectors or Boolean algorithms. |
| `polynomial::quadratic_roots` | [math_DirectPolynomialRoots.cxx](../src/FoundationClasses/TKMath/math/math_DirectPolynomialRoots.cxx), `Solve(A,B,C)`, degree reduction, discriminant, double-root multiplicity and cancellation behavior | Algebraic roots use exact integer discriminants and radical comparisons instead of the source's coefficient threshold, discriminant uncertainty band and floating evaluation/refinement. Independent polynomial-sign/vertex oracles define the exact contract. Native observations preserve multiplicity; no cubic/quartic implementation. |
| Circle/sphere/cylinder line and segment intersections | [IntAna_IntConicQuad.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_IntConicQuad.cxx), line/quadric substitution; [IntAna_Quadric.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_Quadric.cxx); [gp_Sphere.cxx](../src/FoundationClasses/TKMath/gp/gp_Sphere.cxx) and [gp_Cylinder.cxx](../src/FoundationClasses/TKMath/gp/gp_Cylinder.cxx), `Coefficients`; [IntAna2d_AnaIntersection_3.cxx](../src/ModelingData/TKGeomBase/IntAna2d/IntAna2d_AnaIntersection_3.cxx), line/circle empty/tangent/two-point cases | Equivalent implicit equations with exact local differences and unnormalized cylinder axes avoid rounded global coefficient formation. Closed segments clip exact roots before construction. Circle-plane candidates are checked exactly. 174 native polynomial/curved cases, 1,540 independent exact fixtures and a fourth sustained fuzz target. Near-degenerate exact decisions intentionally differ from OCCT tolerance policies; general quadric/curve/surface intersections remain pending. |
| `Bounds3::intersects` | [Bnd_Box.cxx](../src/FoundationClasses/TKMath/Bnd/Bnd_Box.cxx), finite/non-open branch of `IsOut` | Ported separation comparisons, sum of both gaps, inclusive touching. No infinite, void or open-bound box representation. This is broad-phase overlap, not exact collision. |
| `RigidTransform::rotation` | [gp_Trsf.cxx](../src/FoundationClasses/TKMath/gp/gp_Trsf.cxx), axis rotation: translation `origin - R*origin` | Same pivot semantics, independently expressed with Rust vectors. No mirror or scale. Transform tests and native DRAW checks. |
| `Solid::extrude` | [BRepPrimAPI_MakePrism.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakePrism.cxx) and [BRepSweep_Prism.cxx](../src/ModelingAlgorithms/TKPrim/BRepSweep/BRepSweep_Prism.cxx) | Specialized finite normal extrusion of a material profile. Reject thickness at/below linear tolerance, matching the finite prism's vector-length guard. General swept subshapes, oblique/infinite prisms and OCCT history APIs are not implemented. |
| Cylinder seam topology | [BRepPrim_OneAxis.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_OneAxis.cxx), lateral-face pcurves | Rust retains one shared seam edge with two opposed uses and distinct surface pcurves. This is a specialized builder, not a port of the complete one-axis builder. Topology, normals and 66-solid oracle tests. |
| Mass, bounds, classification | [BRepGProp.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepGProp/BRepGProp.cxx), [BRepClass3d_SolidClassifier.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepClass3d/BRepClass3d_SolidClassifier.cxx), `BRepBndLib::AddOptimal` | Independently derived exact prism formulas and profile classification, checked against these native OCCT APIs. The general algorithms have not been ported. |
| DRAW adapter | [BRepTest_PrimitiveCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_PrimitiveCommands.cxx) and [BRepTest_BasicCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_BasicCommands.cxx) | Small command/signature subset delegates to Rust geometry; Tcl remains the real interpreter. `isbbinterf` uses the finite AABB logic above for box solids. |
| Original test judgments | [CheckCommands.tcl](../resources/DrawResources/CheckCommands.tcl), [TestCommands.tcl](../resources/DrawResources/TestCommands.tcl), group/grid `begin`, `end`, `parse.rules` and selected original tests | Execute original files without rewriting expected values. SHA-256 manifest pins every participating upstream file. Known failures, unsupported behavior and absent data are not passes. |

The original copyright/license notices remain in the inherited source. Rust
additions use LGPL-2.1-only WITH OCCT-exception-1.0; algorithmic ports should
retain source attribution here and beside the relevant implementation.

## Spline traceability

Spline evaluation references are
[`Geom_BSplineCurve.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BSplineCurve.cxx)
(`CheckCurveData`, constructors),
[`Geom_BezierCurve.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BezierCurve.cxx),
[`BSplCLib_CurveComputation.pxx`](../src/FoundationClasses/TKMath/BSplCLib/BSplCLib_CurveComputation.pxx)
(`PrepareEval_T`, `BSplCLib_D0/D1/D2`),
[`BSplCLib.cxx`](../src/FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx)
(`NbPoles`, `KnotSequence`, `PoleIndex`, `Eval`, `Bohm`), and
[`PLib.cxx`](../src/FoundationClasses/TKMath/PLib/PLib.cxx) (`RationalDerivative`).
Rust's `curve` module retains degree/multiplicity and homogeneous quotient
semantics, with an independent exact, differentiated de Boor implementation.
It preserves all positive represented weights and strictly distinct knots;
OCCT's weight-resolution tests, knot snapping and extrapolation are deliberately
not applied. Derivatives at interior knots need exact two-sided agreement or an
explicit side. The original
[`Geom_BSplineCurve_Test.cxx`](../src/ModelingData/TKG3d/GTests/Geom_BSplineCurve_Test.cxx)
`SetUp` cubic supplies a native-comparison input; it is not counted as a passing
unchanged upstream test. Evidence: 456 native curve observations, 1,059 independent exact curve fixtures,
invariance tests and the `splines` fuzz target. `spline::KnotVector` follows the
OCCT periodic end-multiplicity, knot-extension and cyclic pole conventions;
parameter normalization and extended knots use exact arithmetic. Supported
periodic directions require more poles than their degree.

Surface references are
[`Geom_BSplineSurface.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BSplineSurface.cxx)
(`CheckSurfaceData`, constructors),
[`Geom_BSplineSurface_1.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BSplineSurface_1.cxx)
(`LocalD0/D1/D2`), and
[`BSplSLib.cxx`](../src/FoundationClasses/TKMath/BSplSLib/BSplSLib.cxx)
(`PrepareEval`, `D0/D1/D2`, `RationalDerivative`). The Rust `surface` module
retains the tensor-product and multivariate quotient semantics with independent
U/V periodicity. Exact pole interpolation and quadrant continuity checks replace
floating evaluation; no tolerance-based weight simplification is applied.
Evidence: 210 native surface jets, 791 exact tensor-basis fixtures, parameter
transpose/weight invariance tests and the sixth fuzz target, `surfaces`.
Curve/surface editing and spline B-rep integration remain pending. High-degree
native discrepancies are [reviewed explicitly](NATIVE_SPLINE_DIVERGENCES.md).
The OCCT 7.6.3 degree-25 second-derivative discrepancy is documented in
[`VALIDATION.md`](VALIDATION.md) and pinned in the reviewed-divergence registry.
It is counted separately from matching native cases, with the exact independent
answer recomputed before accepting the review; the comparison budget is unchanged.

## Spline/plane intersection traceability

The source reference remains `3d097a0328e71b826377d4814ab05ec3c3d23871`:

- [`GeomAPI_IntCS.cxx`](../src/ModelingAlgorithms/TKGeomAlgo/GeomAPI/GeomAPI_IntCS.cxx),
  `Perform`, `Point` and `Segment`: adapted curves/surfaces and point/interval results.
- [`IntCurveSurface_QuadricCurveExactInterUtils.pxx`](../src/ModelingAlgorithms/TKGeomAlgo/IntCurveSurface/IntCurveSurface_QuadricCurveExactInterUtils.pxx)
  and `IntCurveSurface_TheQuadCurvExactHInter.cxx`: substitute the curve into the
  quadric equation and search smooth parameter intervals for roots/zero intervals.
- `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter.cxx`: implicit value
  and tangent/gradient derivative of that scalar equation.
- [`math_FunctionAllRoots.cxx`](../src/FoundationClasses/TKMath/math/math_FunctionAllRoots.cxx):
  sampled zero-interval classification and root refinement. Its numerical
  tolerance policy is deliberately replaced with exact polynomial root counts
  and identically-zero span decisions.

Rust uses exact homogeneous span polynomials, degree-25 Sturm isolation and
Sturm–Tarski signs at algebraic roots. These are independently implemented
mathematical algorithms, not translations of OCCT's sampled root finder. See
`MATHEMATICS.md` for proof references and the subdivision-budget contract.
The native 132-case oracle was captured before Rust implementation. All native
events remain visible; [reviewed divergences](NATIVE_SPLINE_PLANE_DIVERGENCES.md)
require independently verified complete results. The 94 root and 185 spline/plane
fixtures use exact factorization/continued fractions and basis polynomials;
separate sustained campaigns mutate roots and spline intersections.

Explicit interval queries additionally reference
[`Geom_TrimmedCurve.cxx`](../src/ModelingData/TKG3d/Geom/Geom_TrimmedCurve.cxx),
`SetTrim`, preserving positive-length explicit bounds with periodic adjustment
disabled, and
[`IntCurveSurface_InterUtils.pxx`](../src/ModelingAlgorithms/TKGeomAlgo/IntCurveSurface/IntCurveSurface_InterUtils.pxx),
`ComputeAppendPoint`, for native periodic-hit normalization. Rust's query keeps
all parameter events across the supplied interval; it does not implement curve
sense reversal or OCCT tolerance-based extrapolation. Exact shifted knots and
preflight traversal limits cover huge periodic parameters without float wrapping
or uncontrolled turn enumeration. Evidence adds 82 native observations and 98
independent exact fixtures; all prior fixture rows and review pins are retained.

## Spline/quadric intersection traceability

`spline_sphere` and `spline_cylinder` use the same source baseline and
`GeomAPI_IntCS` entry point. Read alongside the sphere/cylinder branches of
`IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection` and
[`IntSurf_Quadric.cxx`](../src/ModelingAlgorithms/TKGeomAlgo/IntSurf/IntSurf_Quadric.cxx),
`Distance`, `Gradient`, and `ValAndGrad`: OCCT evaluates signed distance from
the curve to each primitive and refines roots over C1 intervals. Rust uses
equivalent exact homogeneous implicit polynomials (degree at most 50) and the
shared certified spline intersection engine, including exact contained spans.
The cylinder axis is never rounded to a unit vector. The mathematical root
isolator remains an independent implementation, not an OCCT numerical port.

All 92 native inputs were captured before the Rust quadric implementation.
The 118 independent fixtures add dense/high-degree, full-exponent and
high-multiplicity cases; the existing plane fixtures are unchanged. The shared
fuzz target adds factored squared contact equations and a weighted-line oracle
using independent rational parameter conversion. Native omissions of overlap
intervals and later periodic events are [reviewed explicitly](NATIVE_SPLINE_QUADRIC_DIVERGENCES.md).
The source regression `tests/bugs/modalg_4/bug23076` also exercises curve/surface
intersection, but its external DRAW geometry is unavailable; it is not counted
as an executed or passing upstream regression. General curve/surface
intersection and spline B-rep integration remain pending.

## Linear proximity traceability

The source baseline is unchanged. Read entry points and their implementations:

- [`Extrema_ExtPElC.cxx`](../src/ModelingData/TKGeomBase/Extrema/Extrema_ExtPElC.cxx),
  line projection and accepted parameter ranges;
  [`Extrema_ExtPElS.cxx`](../src/ModelingData/TKGeomBase/Extrema/Extrema_ExtPElS.cxx),
  plane projection and surface parameters.
- [`Extrema_ExtElC.cxx`](../src/ModelingData/TKGeomBase/Extrema/Extrema_ExtElC.cxx),
  line/line normal equations, parallel branch and angular/resolution tests;
  `Extrema_ExtElCS.cxx` and `Extrema_ExtElSS.cxx`, parallel line/plane and
  plane/plane extrema.
- [`BRepExtrema_DistShapeShape.cxx`](../src/ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_DistShapeShape.cxx),
  `Perform` and vertex/edge/face traversal;
  [`BRepExtrema_DistanceSS.cxx`](../src/ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_DistanceSS.cxx),
  elementary pairs, boundary classification, tolerance filtering and infinite-edge handling.
- [`gp_Lin.cxx`](../src/FoundationClasses/TKMath/gp/gp_Lin.cxx), `Distance`;
  [`gp_Pln.hxx`](../src/FoundationClasses/TKMath/gp/gp_Pln.hxx), distance and
  parallel/normal decisions, alongside installed OCCT 7.9.3's earlier inline
  `Distance` implementation; `gp_Dir::Angle`, `IsParallel` and `gp::Resolution`.

`proximity::closest_points` generalizes the line/line stationarity equations to
exact rational face pairs. It retains geometric projection/minimum-distance
semantics but replaces tolerance-based rank/containment and rounded unit normals.
The complete convex face solver and its supporting-plane checker are independent
implementations, not a port of the general B-rep extrema engine. There is one
deterministic Rust minimum witness, not OCCT's list of reported extrema.

The original
[`BRepExtrema_DistShapeShape_Test.cxx`](../src/ModelingAlgorithms/TKTopAlgo/GTests/BRepExtrema_DistShapeShape_Test.cxx)
`BUC60870_EdgeToVertexMinimumDistance` supplies the edge/point input. The native
harness uses its standard default deflection; the original GoogleTest's loose
deflection/EXPECT_NEAR assertion is not replayed. No additional unchanged DRAW
pass is claimed. All 370 native inputs were captured before Rust implementation;
554 independent exact fixtures and the ninth fuzz target validate all ordered
linear-primitive pairings. Unbounded B-rep nonresults and native affine
parallelism differences are [reviewed explicitly](NATIVE_PROXIMITY_DIVERGENCES.md).
General B-rep/curved distances, extrema enumeration and topological attachment
of closest points remain pending.

## Complete linear-set intersections

`intersection::linear_intersection` references the same source baseline and:

- `IntAna_QuadQuadGeo.cxx`, `Perform(gp_Pln,gp_Pln)`: cross-normal intersection
  direction, coincident/parallel cases and the near-parallel origin refinement.
- `IntTools_EdgeEdge.cxx`, line/line branch: incidence, finite-range clipping and
  overlap. `IntTools_EdgeFace.cxx` calls `IntCurveSurface_HInter` and clips its
  surface intersections to the original edge range.
- `BRepAlgoAPI_Common.cxx` / `BOPAlgo_BOP.cxx`, `BuildRC`: explicit filtering
  below the minimum input dimension. `BRepAlgoAPI_Section.cxx` /
  `BOPAlgo_Section.cxx`, `PerformInternal1` and `BuildSection`: contact results.

Rust preserves geometric incidence and bounded overlap semantics with an
independent exact affine/halfspace solver. It returns the entire closed set,
including contacts omitted by COMMON and coplanar polygons with three to six
vertices. Canonical rational constructions are distinct from rounded positions.
This is not a port of the general Boolean engine or its tolerance/history policy.

The 433 native input geometries were captured before implementing the Rust
solver. The initial plane pair from `tests/lowalgos/intss/buc60815` is reused;
its subsequent swept surfaces and original DRAW assertions are not executed.
No new unchanged upstream pass is claimed. Evidence includes 684 independent
fixtures, separate Python/Rust boundary oracles, defining-point permutations,
minimal enclosures and the tenth sustained fuzz target. Native dimension
semantics and [reviewed differences](NATIVE_LINEAR_INTERSECTION_DIVERGENCES.md)
remain explicit. General curved intersections and B-rep splitting are pending.

## Exact Bézier extraction and editing

Source baseline: `3d097a0328e71b826377d4814ab05ec3c3d23871`.

- `GeomConvert_BSplineCurveToBezierCurve.cxx`, both constructors and `Arc`:
  segment a copy, raise interior knot multiplicities, extract consecutive poles.
- `Geom_BSplineCurve.cxx`, `Segment`: locate/insert endpoints, reset periodic
  origin, remove periodicity and enforce the native one-period limit.
- `BSplCLib.cxx`, `InsertKnots`, `Bohm`, `IncreaseDegree`: local refinement,
  derivative interpolation and repeated degree-increment averaging.
- `Geom_BezierCurve.cxx`, `Segment`, `Reverse`, `Increase`: parameter
  substitution, pole/weight reversal and elevation constraints; `PLib.cxx`,
  `Trimming` and `CoefficientsPoles`, plus `BSplCLib::BuildCache` for the
  native floating power-coefficient editing route. `BSplCLib_3.cxx` forwards
  `BuildCache` to `BSplCLib_CurveComputation.pxx`, where `Bohm` derivatives are
  scaled by factorials and span length to produce the power coefficients.

Rust preserves exact shape and parameterization using rational blossom
extraction, de Casteljau subdivision and Bernstein elevation. It keeps original
parameter ranges rather than resetting every result to local `[0,1]`, and
supports explicitly budgeted multi-period extraction. It does not port native
tolerance snapping, out-of-domain Bézier extrapolation, general B-spline knot
editing or topology/history. Read [the full contract](BEZIER_EDITING.md).

The original `Geom_BezierCurve_Test.cxx` cubic and rational Segment/Increase/
Reverse inputs are included in the 546 native captures made before Rust
implementation. The bridge maps the different parameter conventions and
preserves native numerical differences. Source test inputs are reused, not
executed as unchanged GoogleTests. There are no new unchanged DRAW passes.
The 636 independent exact fixtures, coefficient identities, endpoint/jet and
edit-commutation checks, and eleventh sustained fuzz target validate this scope.

## Rule for the next capability

1. Define the standalone kernel input/output, numerical, topology/history and
   failure contracts. Use noBS-CAD only to inform capability scope; application
   integration is deferred. Exclude rendering and unrelated OCCT infrastructure.
2. Read the OCCT entry point **and the algorithm it calls**, including degenerate
   branches, orientation, periodicity, tolerance propagation and history maps.
   Record files, symbols and source revision in this table.
3. Select upstream regressions and native oracle observations before writing the
   Rust implementation. Capture a native result independently of Rust output.
4. Implement idiomatic Rust ownership and errors while preserving the supported
   behavior. Record deliberate restrictions or divergences; never silently
   substitute a mesh or bounding box for an exact operation.
5. Add independent mathematical, invariant and standalone operation-sequence checks. A test becoming
   supported requires a manifest review, not automatic baseline regeneration.
6. Treat discovered OCCT defects as reviewed divergences with an analytic or
   independent reference, not defects we must reproduce for the sake of parity.
