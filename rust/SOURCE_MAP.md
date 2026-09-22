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
