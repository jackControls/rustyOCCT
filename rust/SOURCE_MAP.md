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
| `Solid::extrude` | [BRepPrimAPI_MakePrism.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakePrism.cxx) and [BRepSweep_Prism.cxx](../src/ModelingAlgorithms/TKPrim/BRepSweep/BRepSweep_Prism.cxx) | Specialized finite normal extrusion of a material profile (lines, arcs and, since S8b, nonrational spline segments screened exactly for separation, with exact degree-(p, 1) spline walls: `decide/splines.rs`). Reject thickness at/below linear tolerance, matching the finite prism's vector-length guard. General swept subshapes, oblique/infinite prisms and OCCT history APIs are not implemented. |
| Cylinder topology | [BRepPrim_OneAxis.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_OneAxis.cxx), lateral-face pcurves and seam | Deliberately different: the wall has no seam. Two ring edges bound it through ring loops with winding numbers ±1, pcurves on the universal cover (`TOPOLOGY_MODEL.md`). OCCT's seam, seam vertices and their counts are synthesized by `Topology::occt_counts` for comparison only. Topology, normals, 66-solid oracle tests and the native count checks. |
| Mass, bounds, classification | [BRepGProp.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepGProp/BRepGProp.cxx), [BRepClass3d_SolidClassifier.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepClass3d/BRepClass3d_SolidClassifier.cxx), `BRepBndLib::AddOptimal` | Independently derived exact prism formulas and profile classification, checked against these native OCCT APIs. The general algorithms have not been ported. |
| DRAW adapter | [BRepTest_PrimitiveCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_PrimitiveCommands.cxx) and [BRepTest_BasicCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_BasicCommands.cxx) | Small command/signature subset delegates to Rust geometry; Tcl remains the real interpreter. `isbbinterf` uses the finite AABB logic above for box solids. `nbshapes` reports the kernel's synthesized OCCT counts, as `TopExp::MapShapes` would count the seamed body; `lprops` sums per edge use, as `BRepGProp::LinearProperties` explores edges. |
| `Body::face_from_profile_with`, `Body::wire_from_boundary_with` (S6) | [BRepBuilderAPI_MakeFace.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_MakeFace.cxx), [BRepLib_MakeFace.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepLib/BRepLib_MakeFace.cxx), [BRepBuilderAPI_MakeWire.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_MakeWire.cxx), [BRepBuilderAPI_MakeEdge.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_MakeEdge.cxx) | Independently built cell topology: one planar face with both sides in the void, or a closed wire of the boundary's edges; a circle is a ring edge where OCCT closes the edge at a vertex. Counts follow OCCT's (a free face has no shell, a free edge no wire). The DRAW adapter's `plane`, `cylinder`, `line`, `circle` (with gp_Ax2's automatic X direction, [gp_Ax2.cxx](../src/FoundationClasses/TKMath/gp/gp_Ax2.cxx)), `mkface` and `mkedge` follow [GeomliteTest_SurfaceCommands.cxx](../src/Draw/TKTopTest/GeomliteTest/GeomliteTest_SurfaceCommands.cxx), [GeomliteTest_CurveCommands.cxx](../src/Draw/TKTopTest/GeomliteTest/GeomliteTest_CurveCommands.cxx), [BRepTest_SurfaceCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_SurfaceCommands.cxx) and [BRepTest_CurveCommands.cxx](../src/Draw/TKTopTest/BRepTest/BRepTest_CurveCommands.cxx) for their 3D, bounded forms. |
| `intersection::surface_surface` (S7a) | [IntAna_QuadQuadGeo.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_QuadQuadGeo.cxx) | Independently derived closed forms on exact surfaces with exact degeneracy predicates and certified enclosures; IntAna's tolerance snapping (its angular and linear tolerances, `RefineDir`) is deliberately not reproduced, and its parabola never occurs for a binary64 half-angle. Compared natively on 67 cases. |
| `intersection::ProceduralCurve` (S7b.1) | [GeomInt_IntSS.cxx](../src/ModelingAlgorithms/TKGeomAlgo/GeomInt/GeomInt_IntSS.cxx), IntPatch walking lines | Independently derived: a ruled parameterisation with the discriminant in closed form and exact classes, certified points instead of approximated walking lines; cones on their rulings through the apex or through the cylinder's, with certified root isolation where no closed form exists; a torus with a plane or a sphere on its meridians (S7b.3a), with exact classes and closed-form loop ends, and its special and coaxial pairs as exact circles. Compared natively on 93 cases. |
| `intersection::TracedCurve` (S7b.3b) | [GeomInt_IntSS.cxx](../src/ModelingAlgorithms/TKGeomAlgo/GeomInt/GeomInt_IntSS.cxx), IntPatch walking lines | Independently derived: a torus's curve with a cylinder, a cone or another torus as a certified graph on its meridians (folds by subdivision and the Krawczyk operator, branches by certified windows, tangencies exactly from the pipes' spines and axes), where IntPatch walks and approximates lines and misses isolated tangencies and tiny loops; two cones and a cone's apex on another surface on the cone's projective rulings, unbounded components counted through infinity (S7b.4). Compared natively on 50 cases. |
| `intersection::curve_surface` (S7c.1) | [GeomAPI_IntCS.cxx](../src/ModelingAlgorithms/TKGeomAlgo/GeomAPI/GeomAPI_IntCS.cxx), [IntAna_IntConicQuad.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_IntConicQuad.cxx) | Independently derived: a line's exact polynomial along a plane, cylinder, sphere or torus and a circle's resultant against them, roots isolated exactly with their multiplicities (tangencies) and containment decided exactly, a cone's in certified intervals; GeomAPI_IntCS never reports a contained curve as a segment. S7c.2: `conic_surface` (`Geom_Ellipse`, `Geom_Hyperbola` data) by the same resultant, `spline_torus` exactly and `spline_cone` with exact overlaps and apexes (the cone's `tan` transcendental). Compared natively on 71 cases. |
| `intersection::curve_curve` (S7d.1) | [IntTools_EdgeEdge.cxx](../src/ModelingAlgorithms/TKBO/IntTools/IntTools_EdgeEdge.cxx), [IntAna_IntConicQuad.cxx](../src/ModelingData/TKGeomBase/IntAna/IntAna_IntConicQuad.cxx) | Independently derived: lines by exact linear algebra, a line and a conic by its plane and equation, two conics by S7c's resultant in one plane or by the gcd of their equations along two planes' common line, all exact with tangencies and coincidence; IntTools_EdgeEdge samples and refines edge ranges with tolerances. S7d.2: `spline_curve`, a spline's span polynomials in a conic's plane and equation, overlaps and tangencies exact. Compared natively on 54 cases. |
| `Solid::split_by_plane` (S8) | [BRepAlgoAPI_Splitter.cxx](../src/ModelingAlgorithms/TKBO/BRepAlgoAPI/BRepAlgoAPI_Splitter.cxx), `BOPAlgo_Splitter` | Independently derived: exact sides of the plane, the profile's section by the plane's line with exact crossings and tangencies, pieces rebuilt as prisms and renamed by provenance (S8a.1); oblique planes by two sections (the plane's traces on the caps' planes) lifted to general bodies with ellipse edges and sinusoid pcurves (S8a.2, `solid/split/oblique.rs`); cones, frusta and zones by planes normal to (S8c.1) or containing their axis (S8c.2) and by any other plane in explicit conics with projection pcurves (S8d.2, `solid/split/conic.rs`, D13's certified jets in `jet.rs` and `topology/validate/projection.rs`), whole tori by any plane (S8d.1 normal to or containing the axis, S8d.3 in spiric sections: `solid/split/spiric.rs`), prisms of spline profiles by exact spline/line roots and knot insertion (S8b.3, `solid/split/spline.rs`); OCCT's general builder splits with tolerances, keeps seams, splits walls along tangent rulings and projects a cylinder's section onto it as a B-spline pcurve. Compared natively on 103 cases (26 prisms, 23 cones and spheres, 12 conic sections, 12 tori, 13 spiric tori, 17 spline prisms). |
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
Exact curve/surface and curve knot editing are described below; surface knot editing and spline
B-rep integration remain pending. High-degree
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
tolerance snapping, out-of-domain Bézier extrapolation or topology/history.
General curve knot editing is described separately below. Read [the Bézier contract](BEZIER_EDITING.md).

The original `Geom_BezierCurve_Test.cxx` cubic and RationalIncrease/
RationalReverse inputs are included in the 546 native captures made before Rust
implementation. The bridge maps the different parameter conventions and
preserves native numerical differences. Source test inputs are reused, not
executed as unchanged GoogleTests. There are no new unchanged DRAW passes.
The 636 independent exact fixtures, coefficient identities, endpoint/jet and
edit-commutation checks, and eleventh sustained fuzz target validate this scope.

## Exact tensor patch extraction and editing

The [surface editing contract](SURFACE_EDITING.md) references
[`GeomConvert_BSplineSurfaceToBezierSurface.cxx`](../src/ModelingData/TKGeomBase/GeomConvert/GeomConvert_BSplineSurfaceToBezierSurface.cxx),
constructors, `Patch`, `UKnots` and `VKnots`;
[`Geom_BSplineSurface.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BSplineSurface.cxx),
`segment` and `Segment`;
[`Geom_BezierSurface.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BezierSurface.cxx),
`Segment`, `Increase`, `UReverse`, `VReverse`, `ExchangeUV`, `UIso`, `VIso`;
[`BSplSLib.cxx`](../src/FoundationClasses/TKMath/BSplSLib/BSplSLib.cxx),
`Iso`, `IncreaseDegree`, `BuildCache`; and
[`PLib.cxx`](../src/FoundationClasses/TKMath/PLib/PLib.cxx),
`UTrimming`, `VTrimming`, `CoefficientsPoles`.

Rust raises the clipped boundary multiplicities in each axis using exact local
knot insertion on unit controls, referencing `BSplCLib::{InsertKnots,BoorScheme}`.
The resulting Bernstein block supplies a reusable extraction map; integer dot
products avoid recomputing it for every grid row. This is the same exact map
as spline blossoming, with fewer interpolation stages. Immutable curve
Bernstein operations act on homogeneous tensor rows;
isocurves retain the other original parameter domain. Exact quotient jets
include mixed partials. Positive weights are preserved without OCCT's rational
flag simplification, parameter snapping or one-period segmentation restriction.

The native 691-case corpus was captured before Rust implementation and includes
the `Geom_BezierSurface_Test.cxx` SetUp, RationalSegment, RationalIncrease,
RationalSurface_UIso and VIso_Rational input families. Reused GTest inputs do
not count as unchanged upstream passes. All 4,347 outputs undergo independent
exact tensor coefficient checks. The 745 ordinary complete-control fixtures,
boundary/commutation/overflow/budget tests and twelfth fuzz target add independent
mathematical evidence. Native degree-25 segmentation differences are
[reviewed separately](NATIVE_SURFACE_EDITING_DIVERGENCES.md).

## Exact curve knot refinement and removal

Read at source baseline `3d097a0328e71b826377d4814ab05ec3c3d23871`:

- [`Geom_BSplineCurve.cxx`](../src/ModelingData/TKG3d/Geom/Geom_BSplineCurve.cxx),
  `InsertKnot`, `InsertKnots`, `RemoveKnot`, and their declarations/contracts.
- [`BSplCLib.cxx`](../src/FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx),
  `PrepareInsertKnots`, `InsertKnots`, `BoorScheme`, `RemoveKnot`, `AntiBoorScheme`:
  multiplicity and count preflight, local homogeneous interpolation/inversion,
  cyclic indexing and periodic origin changes.
- [`Geom_BSplineCurve_Test.cxx`](../src/ModelingData/TKG3d/GTests/Geom_BSplineCurve_Test.cxx),
  SetUp, InsertKnot, RemoveKnot, InsertKnots_Multiple and periodic inputs.

The first 665 native operation sequences were captured before the Rust
implementation; [the capture record](fixtures/occt-knot-editing-capture.json)
pins the input/output hashes and records the implementation's absence. Two
supplemental nonconstant seam cases were added during validation, for 667 total.
These reuse source inputs and API operations; they are not unchanged GTest or
DRAW executions. Source and installed runtime versions remain distinct.

`ExactKnotVector` and `ExactBSplineCurve3` preserve rational knots and homogeneous
controls throughout immutable edits. Exact inverse insertion proves removal or
returns no result, with final positive-weight validation. A five-period unroll
handles cyclic indexing and canonical seam origin changes. The complete
supported family, limits and stricter failure rules are in [KNOT_EDITING.md](KNOT_EDITING.md).
This does not add general degree elevation, surface knot editing or spline B-rep.

The 723 complete-control fixtures solve independent Cox power coefficient
equations. A second Rust equation solver checks all full-support identities and
removal feasibility in fuzzing; selected Python cases also use independent
Greville collocation. Native differences are
[reviewed separately](NATIVE_KNOT_EDITING_DIVERGENCES.md). The thirteenth fuzz
target retains a degree-25 timeout regression without relaxing its checks.

## Exact edited-curve intersection traceability

`exact_spline_plane`, `exact_spline_sphere` and `exact_spline_cylinder` retain
the implicit equations and point/interval model of the source entry points
listed above: `GeomAPI_IntCS::{Perform,Parameters,Segment}`,
`IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection`, and
`IntSurf_Quadric::Distance`. Their adapters, tolerance-based roots and interval
branches were read before this extension. `GeomAPI_IntCS_Test.cxx`'s OCC26979
general-surface regressions remain outside the analytic surface family and are
not newly claimed as passes.

[Fresh captures](fixtures/occt-exact-spline-intersection-capture.json) preserve
the existing 214 plane and 92 sphere/cylinder observations before the new Rust
API was implemented. Exact conversion and rational edits retain those inputs'
parameter-to-point functions. All four representations must agree with the
independent complete certificate; the native corpus count remains 306.

The engine now accepts rational knot data directly, preserves exact algebraic
contacts and overlap endpoints, and computes finite output enclosures only on
request. Exact comparisons and positive-weight control hulls do not require
finite binary64 poles or parameters. The existing binary64 entry points use the
same engine and retain their fallible enclosure contract. A separate 274-case
Python Fraction oracle and the fourteenth sanitizer target cover arbitrary
rational inputs. [Limits and output semantics](EXACT_SPLINE_INTERSECTIONS.md)
remain explicit; no general surface intersector or tolerance-based equivalence
is implied.

## Exact B-spline surface knot editing

Source revision remains `3d097a0328e71b826377d4814ab05ec3c3d23871`.
`Geom_BSplineSurface_1.cxx::{InsertUKnots,InsertVKnots,RemoveUKnot,RemoveVKnot}`
and `BSplSLib.cxx::{SetPoles,GetPoles,InsertKnots,RemoveKnot}` were read before
implementation, including transverse homogeneous packing and delegated
`BSplCLib` insertion/removal. The Rust surface uses the existing exact curve
insertion rule as a sparse shared map and its inverse across all transverse rows;
combined refinement preflights both axes
and the complete Cartesian control count before control arithmetic.

[Native capture](fixtures/occt-surface-knot-capture.json) pins 1,748 observations
made before the Rust file existed, including adapted SetUp and U/V insertion/
removal inputs from `Geom_BSplineSurface_Test.cxx`. These are native API
comparisons, not claims that the original GTests run unchanged. The Python
oracle solves every Cox power-coefficient equation over full raw support and
verifies 1,756 complete grids. A separate Rust fraction-free equation solver
checks tensor sequences and removal feasibility during fuzzing. The source pin
and installed OCCT runtime version remain separately reported.

Exact partials reuse differentiated de Boor interpolation weights and integer
dot products across the tensor grid. The 791 independent surface fixtures also
compare against the original de Boor path. Isocurves preserve complete rational
B-splines and compose with the existing certified analytic intersections.
[Contracts and limits](SURFACE_KNOT_EDITING.md) exclude general surface/surface
intersection, arbitrary face trims, general degree elevation and B-rep topology.
[Native discrepancies](NATIVE_SURFACE_KNOT_DIVERGENCES.md) are recorded separately
from actual matches.

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

## General exact degree elevation

Source revision `3d097a0328e71b826377d4814ab05ec3c3d23871`:
`Geom_BSplineCurve::IncreaseDegree`, `Geom_BSplineSurface::IncreaseDegree`,
`BSplCLib::{IncreaseDegreeCountKnots,IncreaseDegree}` and delegated `BSplSLib`
homogeneous packing were read before implementation. The Prautzsch rank-average
construction is implemented as a shared exact sparse map. Unclamped working
padding/cropping and canonical periodic origins are explicit in
[the degree-elevation contract](DEGREE_ELEVATION.md).

The 142 curve and 125 tensor inputs were captured before Rust implementation,
including three adapted original GTests. A pristine source-pinned headless SDK
reproduces those observations. Complete independent Python coefficient systems
and a separate Rust integer solver check all 267 control grids. The native gate
reports 187 matches and 80 fingerprinted discrepancies separately; their range
errors and invalid control counts are not emulated. Fuzz and platform acceptance
passed at Rust revision `27647fad`; the linked contract retains exact campaign
and workflow evidence. This capability adds no rendering or C++ runtime dependency.

## Point-to-rational-spline global minima

At reference revision `3d097a0328e71b826377d4814ab05ec3c3d23871`, the source
review covered `GeomAPI_ProjectPointOnCurve`, `Extrema_ExtPC`,
`Extrema_GGExtPC`/`Extrema_GGenExtPC`/`Extrema_GFuncExtPC`, `math_FunctionRoots`,
and the newer `ExtremaPC_Curve`, `ExtremaPC_BSplineCurve` and `ExtremaPC::Result`
endpoint/status semantics. Native observations preceded implementation.
Inputs were adapted from `ExtremaPC_ExtendedGeometry_Test.cxx` and
`ExtremaPC_Comparison_Test.cxx`; their original assertions were not executed.

Rust preserves the stationary equation, knot partition and closed-domain
endpoints, using exact polynomial isolation and algebraic distance comparisons
instead of sampled searches and numerical root merging. Complete constant
intervals and periodic aliases are explicit. Native contract/numerical
differences have exact observation fingerprints and separate independent proofs.
See [spline proximity](SPLINE_PROXIMITY.md) for the supported mathematical
contract and its clean-revision acceptance at `d1206b15`.

## Spline/line and spline/segment intersections

At reference revision `3d097a0328e71b826377d4814ab05ec3c3d23871`, the source
review covered `IntTools_EdgeEdge::{Perform,Prepare,FindSolutions,
MergeSolutions,AddSolution,FindBestSolution,ComputeLineLine,IsIntersection,
IsCoincident}` and `IntTools_CommonPrt`. Native observations of 31 inputs
preceded implementation. Rust keeps OCCT's result model of isolated vertices
and common parameter ranges. It replaces tolerance-based range searches and
merging with exact common roots, closed polynomial inequalities and exact
parameter identity. Native finite-edge, periodic-edge-range and
`MergeSolutions` coverage differences are recorded in
[spline/linear intersections](SPLINE_LINEAR_INTERSECTIONS.md). A reviewed
native comparison bridge passes on macOS and Linux with reviewed
differences; acceptance passed at `30c5a247`.

## Generic B-rep validation

The source review covered `BRepCheck_Analyzer` and
`BRepCheck_{Vertex,Edge,Wire,Face,Shell,Solid}` (with `BRepCheck_ToolSolid`
and `BRepLib_ValidateEdge`) at `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Native observations of 52 cases preceded implementation. Rust keeps OCCT's
per-subshape status model as typed issues on explicit entities. It replaces
sampled curve/pcurve agreement (23 points by default) with a certified
whole-range bound and replaces classifier tolerances with certified signs and
margins. It also validates every supplied entity, not only subshapes of the
solid. Reviewed differences, including `Closed2d`'s first/last-junction check,
are recorded in [B-rep validation](BREP_VALIDATION.md). The bridge passes on
macOS and Linux with the same reviewed differences; acceptance passed at
`dff912e5`.

## Cell-complex topology

The model (`TOPOLOGY_MODEL.md`) is decided against OCCT's `TopoDS`
(orientations, compounds, seams via `BRep_Tool::IsClosed`), Parasolid's
regions and fins and CGM's cells, and ports none of them. Regions partition
space, faces have front and back sides listed by shells, fins sit in radial
order on their edges, and periodic surfaces carry no seams. The native
bridges keep OCCT's seamed encoding: a seam is an edge closed on its face
(`BRep_Tool::IsClosed(edge, face)`), and its vertices are those used only by
seams and by edges closed on them. Those native entities are classified
structure-only by that rule, and every comparison verifies the rule with the
synthesized counts (`TopExp::MapShapes` over vertices, edges, wires, faces,
shells and solids).

## Value identity and operation history

The source review for history covered `BRepPrimAPI_MakePrism`
(`Generated`, `FirstShape`, `LastShape`, `IsDeleted`), `BRepSweep_Prism`,
`BRepTools_History`, `BRepBuilderAPI_Transform::Modified` and DRAW's
`BRepTest_SweepCommands.cxx` (`prism`) and `BRepTest_HistoryCommands.cxx`
(`savehistory`, `generated`, `modified`, `isdeleted`) at
`3d097a0328e71b826377d4814ab05ec3c3d23871`. Native observations of 98 cases
preceded implementation. OCCT has no value ids: a shape's identity is its
`TShape` pointer and location. Rust ids are digests of derivations, and every
operation returns a complete, independently checked history. Rust reports
OCCT's first/last shapes as `Generated` with start/end roles. The solid
OCCT generates from the face corresponds to the solid region, generated from
every boundary label; bodies are related by id. See
[identity and history](IDENTITY_AND_HISTORY.md).

## Height split and stacked fuse

The source review covered `BRepAlgoAPI_Splitter` (`BOPAlgo_Splitter`),
`BOPAlgo_Builder::PrepareHistory`, `LocGenerated` and `LocModified`
(`BOPAlgo_Builder_4.cxx`), `BRepAlgoAPI_Fuse`,
`ShapeUpgrade_UnifySameDomain::History`, `BRepTools_History::Merge` and DRAW's
`bapisplit`, `bapibop` (`BOPTest_APICommands.cxx`) and `unifysamedom`
(`SWDRAW_ShapeUpgrade.cxx`) at `3d097a0328e71b826377d4814ab05ec3c3d23871`.
Native observations of 160 scenarios preceded implementation. Rust does not
port the general Boolean builder: `Solid::split_at_height` and
`Solid::fuse_stacked` rebuild prisms of one profile and frame and name every
entity from the inputs by rule. OCCT reports a split parent as `Modified`
into its pieces, shares one cut face between them and relates solids; Rust
reports `Split`, gives each piece its own generated cut face, edges and
vertices, and splits and merges the solid region. Unification also merges
coplanar walls of the profile itself, which the kernel keeps (one reviewed
difference). See [identity and history](IDENTITY_AND_HISTORY.md).

## Enclosures

The source review for M5 covered `BRep_Tool::Tolerance`, `BRep_Builder`'s
`MakeVertex`, `MakeEdge`, `MakeFace` and `UpdateEdge` with
`BRep_TEdge::UpdateTolerance` (tolerances are stored as requested and only
grow), `BRepCheck_Vertex` (curve and curve-on-surface ends against
`max(vertex, edge)` tolerance), `BRepCheck_Edge` and `BRepLib_ValidateEdge`
(23 control points by default, or `BRepLib_CheckCurveOnSurface`
maximization) at `3d097a0328e71b826377d4814ab05ec3c3d23871`. Native
observations of 52 cases preceded implementation
(`fixtures/occt-enclosure-preimplementation`). OCCT tolerances are requested
bounds that algorithms grow; kernel enclosures are certified bounds computed
from the stored geometry, never grown to make a step pass. Imported OCCT
tolerances become `Imported` enclosures with `BRepCheck`'s acceptance rules,
and the checker verifies them. See
[enclosures](BREP_VALIDATION.md#enclosures-m5).

## OCCT `.brep` interop

The source review covered `dox/specification/brep_format.md`,
`BRepTools_ShapeSet` (geometry and shape records, edge representations,
`CN` continuity glued to a closed surface's second pcurve),
`TopTools_ShapeSet` (backward record numbering, flags, orientations),
`TopTools_LocationSet` (matrix and composite location records),
`GeomTools_CurveSet`, `GeomTools_Curve2dSet` and `GeomTools_SurfaceSet`
(every curve and surface record, trimmed curves over a basis, B-spline
rational and periodic flags), `BRep_Tool::CurveOnSurface` and
`CurveOnPlane` (the second pcurve serves a reversed seam use; planes may omit
pcurves) and DRAW's `restore` and `explode` (`DBRep.cxx`: the type from its
first letter, duplicates skipped in `TopExp_Explorer` order) at
`3d097a0328e71b826377d4814ab05ec3c3d23871`. Rust's `occt_brep` reads every
record of versions 1 to 3 and converts only what the cell model represents:
plane and cylinder faces (direct or indirect), line and circle edges, rigid
locations. Seams merge into periodic loops by the rule of
`cell_reference.to_cell`; everything else is reported by name and counted,
never approximated. The writer emits version 1 and inserts a seam per wound
face at a `u` where both wound loops have a vertex, with a seam vertex on a
ring edge. Native observations of the unmodified `data/occ` corpus were
captured after implementation, unlike every earlier milestone
(`fixtures/occt-brep-io-capture/NOTES.md`). OCCT's exploration order is never
reproduced: the DRAW adapter selects picks by native geometry. See
[validation](VALIDATION.md#occt-brep-interop) and
[the coverage ledger](UPSTREAM_TESTS.md#structure-mapping-and-the-coverage-ledger).

## Cones (S3)

The source review covered `BRepPrimAPI_MakeCone`, `BRepPrim_Cone` and
`BRepPrim_OneAxis` (the generatrix from `(R1, 0)` to `(R2, H)`, the lateral
face's seam at angle 0 along the frame's x, an apex as one vertex with a
degenerated edge), `Geom_ConicalSurface` and `ElSLib::ConeValue` (the
parameterization `O + (R + v sin a)(cos u X + sin u Y) + v cos a N`),
`GeomTools_SurfaceSet` (record 3: location, axis, X, Y, radius,
semi-angle), `BRepPrimAPI_MakeRevol`, `BRepSweep_Revol`, `BRepSweep_Rotation`
and `BRepSweep_NumLinearRegularSweep::IsUsed` (the revolve history:
`IsDeleted` is `!IsUsed`, invariant axis points and edges generate nothing
or a degenerated edge) and `GeomAdaptor_SurfaceOfRevolution::GetType` (a
revolved line is a cone only when `|cos a| <= 1 - Precision::Confusion()`)
at `3d097a0328e71b826377d4814ab05ec3c3d23871`.

| Rust | OCCT | Notes |
| --- | --- | --- |
| `Surface::Cone` | [Geom_ConicalSurface.cxx](../src/ModelingData/TKG3d/Geom/Geom_ConicalSurface.cxx) | Same parameterization, `v` along the ruling; the apex is a pole (a vertex loop), not a degenerated edge. |
| `Solid::cone_with` | [BRepPrimAPI_MakeCone.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeCone.cxx), [BRepPrim_Cone.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_Cone.cxx) | Full turn only; equal radii (a cylinder) and two apices are refused. History as revolving the meridian ([BRepPrimAPI_MakeRevol.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeRevol.cxx)), with the caps and the region reported, which `MakeRevol` reports deleted. |
| `occt_brep` cones | [GeomTools_SurfaceSet.cxx](../src/ModelingData/TKGeomBase/GeomTools/GeomTools_SurfaceSet.cxx) | The degenerated apex edge becomes the pole and is written back by rule. |
| Mass properties | [BRepGProp.cxx](../src/ModelingAlgorithms/TKTopAlgo/BRepGProp/BRepGProp.cxx) | Exact integration in certified intervals instead of Gauss quadrature. |

## Spheres (S3)

The source review covered `BRepPrimAPI_MakeSphere` and `BRepPrim_Sphere`
(the meridian circle about `-y` through `x`, offset `2π`; latitudes
`angle1..angle2`; an end at `±π/2` a pole with a degenerated edge, others a
circle bounding a disc), `Geom_SphericalSurface` and `ElSLib::SphereValue`
(`O + R (cos v (cos u X + sin u Y) + sin v N)`), `GeomTools_SurfaceSet`
record 4, and DRAW's `psphere` (angles in degrees, converted by
`* M_PI / 180`), at `3d097a0328e71b826377d4814ab05ec3c3d23871`.

| Rust | OCCT | Notes |
| --- | --- | --- |
| `Surface::Sphere` | [Geom_SphericalSurface.cxx](../src/ModelingData/TKG3d/Geom/Geom_SphericalSurface.cxx) | Same parameterization; a whole sphere is a face without loops, a pole closing a band a vertex loop. |
| `Solid::sphere_with` | [BRepPrimAPI_MakeSphere.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeSphere.cxx), [BRepPrim_Sphere.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_Sphere.cxx) | Full longitude only. History as revolving the meridian, the arc generating the wall. |
| `occt_brep` spheres | [GeomTools_SurfaceSet.cxx](../src/ModelingData/TKGeomBase/GeomTools/GeomTools_SurfaceSet.cxx) | Poles and passes through them are OCCT's degenerated edges, rebuilt by rule when writing. |

## Tori (S3)

The source review covered `BRepPrimAPI_MakeTorus` and `BRepPrim_Torus`
(the minor circle about `O + R x` in the half-plane of `x` and the axis,
latitudes `angle1..angle2` revolved by `angle`), `BRepPrim_OneAxis` (a
closed meridian has no top or bottom face; a partial turn adds the end
faces; `isHeightInverted` reverses only a wedge's end faces),
`Geom_ToroidalSurface` and `ElSLib::TorusValue`, `GeomTools_SurfaceSet`
record 5, and DRAW's `ptorus`, at
`3d097a0328e71b826377d4814ab05ec3c3d23871`.

| Rust | OCCT | Notes |
| --- | --- | --- |
| `Surface::Torus` | [Geom_ToroidalSurface.cxx](../src/ModelingData/TKG3d/Geom/Geom_ToroidalSurface.cxx) | Same parameterization; ring tori only; a whole torus is a face without loops, loops may wind in `v`. |
| `Solid::torus_with` | [BRepPrimAPI_MakeTorus.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeTorus.cxx), [BRepPrim_Torus.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_Torus.cxx), [BRepPrim_OneAxis.cxx](../src/ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_OneAxis.cxx) | The whole torus, v-segments (full turn) and wedges (whole tube from `v = 0`); a segment of a partial turn is refused. The wall faces out: OCCT's inner half is inside out. History as revolving the meridian. |
| `occt_brep` tori | [GeomTools_SurfaceSet.cxx](../src/ModelingData/TKGeomBase/GeomTools/GeomTools_SurfaceSet.cxx) | Seams in `u` and `v` merged on import and rebuilt by rule when writing; horn and spindle tori unsupported. |

## Spline continuity (R4)

| Rust | OCCT | Notes |
| --- | --- | --- |
| `topology::validate::continuity` | [BSplCLib.cxx](../src/FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx) `RemoveKnot`, [Geom_BSplineCurve.cxx](../src/ModelingData/TKG3d/Geom/Geom_BSplineCurve.cxx) `Continuity` | OCCT measures continuity by knot multiplicity and removes knots within a tolerance; the kernel decides C1 by one exact homogeneous removal with zero residual at each knot of multiplicity equal to the degree. |

## Tessellation (T-a)

The source review covered `BRepMesh_IncrementalMesh` and
`IMeshTools_Parameters` (`Deflection`, `Angle`, `Relative`, `InParallel`,
`ControlSurfaceDeflection`, `MinSize`), `BRepMesh_EdgeDiscret` and
`BRepMesh_CurveTessellator` (an edge discretized once, its polygon on each
face's triangulation, a seam's two polygons), `BRepMesh_Deflection`, the
range splitters of cylinders, cones, spheres and tori (the parameter grid
of a periodic face, meshed between its seam's two uses, so the seam's
nodes appear twice in the face's triangulation),
`BRepMesh_DelaunayBaseMeshAlgo`, `BRepMesh_Delaun` and
`BRepMesh_DelaunayDeflectionControlMeshAlgo` (a Delaunay triangulation in
scaled parameters, then up to 11 passes testing each triangle's centre and
its links' middles against the surface and inserting nodes; frontier links
are not split; the largest sampled deviation is stored as the face's
`Poly_Triangulation::Deflection`), at
`3d097a0328e71b826377d4814ab05ec3c3d23871`.

| Rust | OCCT | Notes |
| --- | --- | --- |
| `tessellation::Parameters` | [IMeshTools_Parameters.hxx](../src/ModelingAlgorithms/TKMesh/IMeshTools/IMeshTools_Parameters.hxx) | Absolute deflection and angle only; no relative mode, minimum size or parallel faces yet. The angle is limited to `(0, π/2]`. |
| edge polylines | [BRepMesh_EdgeDiscret.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_EdgeDiscret.cxx), [BRepMesh_CurveTessellator.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_CurveTessellator.cxx) | Uniform in the edge's fraction, the count from a certified chord bound and each curved face's thin-triangle condition; no seam polygons, no degenerated edges (a pole is one node). |
| face charts | [BRepMesh_CylinderRangeSplitter.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_CylinderRangeSplitter.cxx), [BRepMesh_ConeRangeSplitter.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_ConeRangeSplitter.cxx), [BRepMesh_SphereRangeSplitter.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_SphereRangeSplitter.cxx), [BRepMesh_TorusRangeSplitter.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_TorusRangeSplitter.cxx) | OCCT meshes a periodic face between its seam's two uses; the kernel maps the whole domain homeomorphically into the plane (sinusoidal and annulus charts), so there is no seam and nothing to join. |
| `tessellation/cdt.rs` | [BRepMesh_Delaun.cxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_Delaun.cxx) | Exact orientation predicates, a filtered in-circle that never flips on doubt, Sloan's constraint recovery; crossing chords refine the edges instead of being repaired. |
| refinement and bounds | [BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx](../src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx) | OCCT samples each triangle's centre and links' middles for at most 11 passes and records what remains; the kernel certifies every triangle's deviation over the whole triangle (`MATHEMATICS.md`) and refines until it holds. |
| `Mesh` | [Poly_Triangulation.hxx](../src/FoundationClasses/TKMath/Poly/Poly_Triangulation.hxx), [Poly_PolygonOnTriangulation.hxx](../src/FoundationClasses/TKMath/Poly/Poly_PolygonOnTriangulation.hxx) | One node array for the whole body, triangles grouped by face and oriented out of the material; OCCT stores one triangulation per face in its surface's orientation and joins them through the edges' polygons. |

## STEP import (STEP-a)

The source review covered OCCT's Part 21 lexer and grammar (`step.lex`,
`step.yacc`, `StepData_StepReaderData`: simple and complex instances,
typed parameters, `\S\` and the other string directives, the header),
`STEPControl_Reader` and `STEPControl_ActorRead` (roots, units through
`StepData_Factors`, the uncertainty as the healing precision while entities
are built at `Precision::Confusion`), `StepToTopoDS_Builder` (a manifold
solid, voids reversed when their oriented shell is `.F.`),
`StepToTopoDS_TranslateShell`, `StepToTopoDS_TranslateFace` (a wire reversed
when its bound's orientation differs from `same_sense`, the face reversed
when `same_sense` is false, a vertex loop alone on a sphere the natural
bounds), `StepToTopoDS_TranslateEdgeLoop`, `StepToTopoDS_TranslateEdge`
(vertices swapped when an edge curve's `same_sense` is false),
`StepToTopoDS_TranslateVertexLoop` (a degenerated edge) and `StepToGeom`
(placements, lines, circles, elementary surfaces), with `TopoDSToStep`'s
writer for the fixtures' conventions (degenerated edges dropped, voids
written reversed), at `3d097a0328e71b826377d4814ab05ec3c3d23871`.

| Rust | OCCT | Notes |
| --- | --- | --- |
| `step::part21` | [step.lex](../src/DataExchange/TKDESTEP/StepFile/step.lex), [step.yacc](../src/DataExchange/TKDESTEP/StepFile/step.yacc), [StepData_StepReaderData.cxx](../src/DataExchange/TKDESTEP/StepData/StepData_StepReaderData.cxx) | Schema-free; strings kept raw; a dangling reference, a duplicate entity number or a nesting deeper than 64 is an error, where OCCT reports it and goes on. |
| `step::import` | [STEPControl_Reader.cxx](../src/DataExchange/TKDESTEP/STEPControl/STEPControl_Reader.cxx), [STEPControl_ActorRead.cxx](../src/DataExchange/TKDESTEP/STEPControl/STEPControl_ActorRead.cxx) | Every solid and surface-model shell in entity order, in its representation's units; no product structure, no placements (`AssemblyPlacement`), no shape healing. The uncertainty is each entity's imported tolerance, which the validator verifies; OCCT uses it as its healing precision instead. |
| body, shell and face translation | [StepToTopoDS_Builder.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_Builder.cxx), [StepToTopoDS_TranslateShell.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateShell.cxx), [StepToTopoDS_TranslateFace.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateFace.cxx) | The same orientation rules, into an `occt_brep::Document` that the `.brep` converter turns into cells. |
| edges, loops and vertices | [StepToTopoDS_TranslateEdge.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdge.cxx), [StepToTopoDS_TranslateEdgeLoop.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdgeLoop.cxx), [StepToTopoDS_TranslateVertex.cxx](../src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateVertex.cxx) | An edge runs along its curve between its vertices' projections. OCCT projects a missing pcurve and repairs with ShapeFix; the kernel derives rulings, parallels and meridians exactly as straight segments on the universal cover and adds the degenerated edge at a pole on the face's side. |
| geometry | [StepToGeom.cxx](../src/DataExchange/TKDESTEP/StepToGeom/StepToGeom.cxx) | Placements by ISO 10303-42's `build_axes` (OCCT's `gp_Ax2` picks another `x` when the reference direction is absent, which moves no point); a trimmed curve as its basis with its sense; a torus with a negative major radius is `NegativeMajorRadius` (OCCT flips the face). |
