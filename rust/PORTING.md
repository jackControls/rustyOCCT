# Kernel scope and implementation order

This is a capability-driven Rust implementation, not a translation of every
OCCT package or a C++ wrapper. Prioritize reliability, then performance, then
ease of integration. Mathematical correctness and standalone kernel robustness
come first. noBS-CAD's product directions inform scope, not acceptance criteria;
application adapters and saved-project replay are deferred.

## Historical capability-scope reference

The initial scope was informed by noBS-CAD commit
`da92416309debfe52f63530804683c837e3b7c35`. This remains useful context for which
geometry families to include, without making its project architecture a kernel
dependency or making application integration the next milestone:

- [`crates/solid/src/dto.rs`](https://github.com/jackControls/noBS-CAD/blob/da92416309debfe52f63530804683c837e3b7c35/crates/solid/src/dto.rs):
  all 14 `KernelJobDto` variants, profiles/curves, scene and topology contracts.
- [`crates/occt/src/native.rs`](https://github.com/jackControls/noBS-CAD/blob/da92416309debfe52f63530804683c837e3b7c35/crates/occt/src/native.rs)
  and [`shim.cpp`](https://github.com/jackControls/noBS-CAD/blob/da92416309debfe52f63530804683c837e3b7c35/crates/occt/src/shim.cpp):
  actual modeling, measurement, tessellation, drawing and STEP behavior.
- [`crates/occt/build.rs`](https://github.com/jackControls/noBS-CAD/blob/da92416309debfe52f63530804683c837e3b7c35/crates/occt/build.rs):
  19 linked modeling/data-exchange toolkits; no OCCT visualization toolkits.
- [`docs/goals.md`](https://github.com/jackControls/noBS-CAD/blob/da92416309debfe52f63530804683c837e3b7c35/docs/goals.md),
  `docs/2D_DRAWINGS.md`, `docs/ASSEMBLIES.md`, `docs/CAM_EDGE_CHAINS.md` and
  `docs/cam/README.md`: accepted current and future product needs.

The upstream source reference is pinned to
`3d097a0328e71b826377d4814ab05ec3c3d23871`. The first *runtime oracle* is installed
OCCT 7.9.3, matching the currently used SDK; it is not a build of that fork point.
Record the source and runtime versions separately on every future comparison.
Read [the source map and porting rule](SOURCE_MAP.md) before implementing each
capability; original OCCT code is the behavioral reference, alongside its tests.

## Current job coverage

"Partial foundation" below means callable Rust geometry exists. It does not
mean that the complete noBS-CAD job or its DTO adapter has been implemented.

| noBS-CAD job | Required geometry | Rust status / remaining work |
| --- | --- | --- |
| Extrude | Profile and exact planar-face extrusion, offsets, taper, holes, join/cut/intersect | **Partial foundation:** polygon/circle New Body geometry with holes and signed offsets. Need arcs, source-face adapter, taper and Booleans. |
| Revolve | Trimmed revolutions about sketch axes, partial/full turns, Booleans | Planned. |
| Sweep | Analytic/spline paths, frames, orientation modes, transitions and guide rail | Planned. |
| Loft | Ruled/smooth sections, continuity, centerline and guide rail | Planned. |
| Rib | Thin profile construction and Boolean attachment | Extrusion foundation only; job planned. |
| Fillet | Selected edges, constant-radius blend, tangent chains, topology history | Planned. |
| Chamfer | Selected edges, edge/face adjacency, tangent chains | Planned; explicit adjacency exists. |
| Hole | Blind/through, stepped/countersunk bottoms, modeled internal threads | Through voids exist **only as input profile holes**. Cutting an existing solid, drilling geometry and thread operations are planned. |
| ExternalThread | Selected cylinder, exact helical profiles, handedness and fit geometry | Planned. |
| Shell | Offset faces, remove selected faces, sew/validate result | Planned. |
| Transform | Translation, rotation, rigid poses, mirror, patterns | **Partial foundation:** translate/rotate supported for current solids. Mirrors and body/job mapping planned; pattern meaning stays in noBS-CAD. |
| Combine | Fuse/cut/common, many tools, keep-tools semantics | Planned. |
| SplitBody | Exact planar splitting into retained solids | Planned. |
| ImportStep | STEP topology, units, analytic and spline geometry | Planned; no mesh substitution. |

## Required query and interchange coverage

| Contract / use | Required kernel work | Status |
| --- | --- | --- |
| Scene recompute | Owned bodies, feature-local errors, atomic operation results and repeatable replay | Immutable standalone solids and typed errors exist; scene/job adapter planned. |
| Topology references | Face/edge identities, membership, geometry signatures, generated/modified/deleted mappings | Shared topology and extrusion face origins exist. Body-local IDs are not persistent naming. History mappings required before feature migration. |
| Face/edge metadata | Plane, cylinder, cone, circle, curvature, edge lengths and face signatures | Planes/cylinders/circles retained; DTO signatures and additional queries planned. |
| Measurement | Bounds, mass/area/centroid/inertia, point classification, extrema/closest points | Bounds, mass properties and point classification supported for current prisms; certified linear-set minimum distances exist. Complete point-to-rational-spline minimum sets are implemented and accepted at `d1206b15`. General B-rep and other curved-pair distance/extrema remain planned. |
| Exact interference | Occurrence transforms, minimum clearance, closest points, overlap volume | Transform, classifier and certified linear-set distance foundations exist. Solid clearance, common-solid and multi-body queries remain planned. |
| Tessellation | Deflection-controlled watertight triangles, normals, face ranges, shared edge samples | Planned geometry output, independent of any renderer. Preserve f64 exact geometry; f32 output conversion belongs at the adapter boundary. |
| STL / 3MF | Feed trustworthy tessellation to `nbcad-export` | Keep existing Rust writers; do not port OCCT mesh file writers or build a second export stack. |
| STEP export | Exact AP242 geometry, units, names, long thread metadata and assembly placements | Planned; round-trip with imports and independent readers is a release gate. |
| Drawings | Exact orthographic visible/hidden curves, section/slab cuts, circles and anchors | Planned geometry algorithms (including HLR), not a drawing UI or raster renderer. |
| CAM | Curve/face intersections, offsets, exact analytic boundaries, surface normals, closest/contact queries | Planned as geometry services. Tool libraries, postprocessors, feeds, motion planning and stock simulation remain in noBS-CAD. |
| Future analysis | Reliable closed geometry, geometric repair diagnostics, suitable meshing boundary and stable selections | Preserve these interfaces; a tetrahedral mesher, FEA solver, materials and dynamics are separate validated systems. |

## Implementation order and acceptance gates

1. **Mathematical foundation — current priority.** The prism implementation
   provides a working test subject. Exact finite-f64 2D orientation, independent
   rational/integer oracles and generated geometric invariants are implemented,
   along with 3D orientation/insphere, exact quadratic roots, certified
   line/segment intersections with planes/triangles/circles/spheres/cylinders and
   sustained structured fuzzing. Exact point/line/segment/plane/triangle minimum
   distances retain rational witnesses and certified comparisons. Certified rational Bézier/B-spline curves and surfaces, including periodic
   directions and derivatives through order two, also exist. Exact Bézier
   extraction and subdivision/trimming/reversal/degree elevation preserve
   homogeneous rational controls and original parameter units. Tensor patch
   extraction through degree 25, independent U/V edits and isocurves also exist,
   including exact mixed partials and Cartesian output limits. Curve knot
   refinement and removal through degree 25 preserve homogeneous functions,
   including periodic seams and unclamped controls, with rational knots
   throughout. Surface knot refinement and exact removal apply the same
   contract to the complete tensor grid, with independent U/V periodicity and
   atomic Cartesian limits. General exact curve/tensor degree elevation now
   has complete independent coefficient reconstruction, source-pinned native
   checks, sustained fuzzing and published platform acceptance at revision
   `27647fad`; see `DEGREE_ELEVATION.md` for the evidence. Exact degree-25
   root isolation and certified spline/plane/sphere/cylinder intersections cover
   isolated contacts
   and whole overlap intervals, including explicit trimmed ranges and multiple
   periodic turns with exact clipping and traversal limits; quadratic surfaces
   use internal equations through degree 50. The same complete queries accept
   exact edited curves and rational trim bounds, retaining algebraic contacts
   and overlaps independently of optional binary64 enclosures. Complete
   spline/line and spline/segment preimages add exact closed clipping and
   cross-span parameter ordering; their native comparison gate is pending. Extend certified comparisons and constructions
   to curved geometry and propagate uncertainty through topology changes.
   Continue minimizing fuzz failures. See `MATHEMATICS.md` and `FUZZING.md`.
2. **Topology invariants and operation history.** Strengthen generic B-rep
   validation beyond sampled prism checks: connectivity, orientation, shells,
   cavities, seams and curve/pcurve/surface consistency. Define generated,
   modified and deleted mappings with explicit split/merge ambiguity before
   topology-changing operations. Body-local indices are not persistent names.
3. **Reliable geometry and intersections.** Circular arcs, trimmed curves,
   arbitrary face trimming; extend the existing certified
   curve/surface jets with
   projection, general intersections, trimming and sewing. Every
   numerical algorithm needs a declared domain, degeneracy behavior and error
   evidence. Introduce acceleration after correctness and workload benchmarks.
4. **Booleans and mechanical features.** Split/classify/assemble solids and
   preserve source history, then fuse/cut/common, many-tool operations, plane
   splits, drilling, revolutions and rib attachment. Cover tangency, coincident
   faces, tiny edges, periodic seams and disconnected results. Never replace a
   failed exact operation with a triangle or bounding-box approximation.
5. **Interchange, meshing and operational guarantees.** Deflection-controlled
   watertight tessellation, STEP import/export and diagnosis. Begin each as soon
   as its geometry exists; compare units, topology and geometric error against
   independent readers. Bound memory/work, support cancellation for long
   operations, test atomic failure and deterministic standalone operation
   sequences. Add real WASM execution and kernel performance budgets.
6. **Advanced modeling and geometric queries.** Sweeps, lofts, offsets/shells,
   chamfers, fillets, modeled helices/threads, exact section/HLR and proximity
   queries. Each capability passes the same mathematical, topology, interchange
   and operational gates. Application migration is separate future work.

No dates or full-kernel equivalence are implied by this ordering. A robust
CAD kernel is a substantial continuing effort, and each stage needs its own
failure corpus and measurable acceptance evidence.

## Code boundaries

Start with one coherent `rusty-occt` crate instead of empty crates for every
future subsystem. Split crates when there is a working boundary to isolate.

- `math`: f64 points/vectors, orthonormal frames, proper rigid transforms,
  dimensional tolerance and bounds. No global mutable precision state.
- `predicates`: exact signs for represented coordinates, separate from
  approximate constructions and tolerance-band decisions.
- `intersection`: validated three-point planes/triangles, exact intersection
  classifications, analytic circles/spheres/infinite cylinders, and rational
  or algebraic constructions with binary64 enclosures. Complete linear-set
  intersections additionally retain canonical rational points, segments, filled
  convex polygons, infinite lines and planes without eager output rounding.
  Rational spline curves against lines and closed segments return complete
  algebraic parameter preimages, including whole overlap intervals.
- `polynomial`: specialized exact quadratics and general real-root isolation
  through degree 25, retained algebraic identity/multiplicity, polynomial signs
  at exact roots and certified comparisons/enclosures.
- `proximity`: complete minimum-distance queries among points, lines, closed
  segments, planes and closed triangles; one exact rational closest pair,
  original affine parameters and independently requested output enclosures.
  Point-to-rational-spline queries additionally retain all algebraic minimum
  parameters and whole minimizing intervals, including periodic aliases.
- `curve`: validated periodic/nonperiodic rational Bézier/B-spline data, exact homogeneous
  evaluation, derivative continuity decisions and minimal output enclosures;
  exact Bézier extraction and immutable editing without rounding control data;
  exact rational B-spline knot refinement/removal and evaluation.
- `spline`: validated binary64 and rational knot vectors, exact periodic extension and parameter wrapping.
- `surface`: immutable rational tensor-product patches and certified partial derivatives;
  independent U/V periodicity and side selection; exact patch extraction,
  rectangular editing and isocurves; exact B-spline surface refinement/removal
  across complete U/V grids. These are not yet B-rep faces.
- `profile`: validated material boundaries, containment and planar moments.
  Polygons are normalized CCW; holes are assigned orientation in the B-rep.
- `topology`: immutable owned vertices/edges/faces, opposite oriented uses,
  exact surfaces and per-face pcurves. A periodic seam is one edge with two
  uses on the same face. Extrusion provenance is explicit.
- `solid`: validated normal prisms and geometry-derived queries. The current
  retained profile is an exact construction model for these solids, not a
  general-purpose post-Boolean shape representation.
- `tools/occt_oracle.cpp`: optional test executable only. Never FFI-linked or
  silently invoked by the Rust library.

All current constructors return `Result` before publishing a solid. No
renderer, OS calls or native SDK dependency is required by the library, and the
kernel crate forbids `unsafe` code.
`num-bigint` supplies exact integer arithmetic for predicates/constructions;
`num-rational` supplies exact rational spline arithmetic. The libFuzzer runtime
is isolated to the fuzz workspace. Dependencies must have an actual geometry need; dependency count
alone is not a correctness or speed metric.

## OCCT families to study, not mechanically translate

- `TKernel` / `TKMath`: precision, linear algebra, bounds, numerical geometry.
- `TKG2d`, `TKG3d`, `TKGeomBase`, `TKGeomAlgo`: analytic and spline geometry,
  conversion, parameterization, extrema and intersections.
- `TKBRep`, `TKTopAlgo`, `TKPrim`: shared topology, building/checking,
  classification, properties and primitives.
- `TKBO`, `TKBool`, `TKShHealing`: splitting, Boolean assembly and repair.
- `TKFillet`, `TKOffset`: blends, offsets, sweeps, lofts and shells.
- `TKMesh`: surface/edge tessellation for export and consumer display data.
- `TKHLR`: technical drawing geometry.
- `TKXSBase`, `TKDE`, `TKDESTEP`: required STEP exchange semantics.

Do not port `Visualization` (AIS/V3d/OpenGL/Vulkan/viewers), `Draw`/Tcl test UI,
windowing, shaders, OCAF/TDF document/history frameworks, or unused formats such
as IGES/VRML. Existing upstream copies remain reference material only. Do not
duplicate noBS-CAD's sketch solver, feature planner, history/undo, assemblies,
MCP, scripting, UI, rendering, file schema or manufacturing application logic.
