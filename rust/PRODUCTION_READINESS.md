# Kernel production readiness

The kernel must be dependable independently of any application. noBS-CAD's
requirements help choose capability scope; its documents, DTOs and feature
history are not the definition of geometric correctness or a prerequisite for
kernel work. Application integration is deferred.

The current implementation is an early prism kernel, not a production-ready
replacement for OCCT. The priority is **mathematical correctness and numerical
robustness, then topology/history, then interchange and operational reliability**.
Test and fuzz every layer as it is built. Passing upstream tests is evidence,
not proof of correctness for all inputs or a reason to reproduce an OCCT defect.

## Release gates, by capability

Promote kernel capabilities only after all applicable gates pass. Use
`PORTING.md` for capability scope and `MATHEMATICS.md` for numerical contracts.

| Gate | Evidence required | Current state |
| --- | --- | --- |
| Mathematical decisions | Correct predicate signs over declared domains, independent exact references, documented arithmetic assumptions | Exact finite-f64 2D/3D orientation, insphere and linear-set distance decisions exist; incircle and curved/B-rep distances pending |
| Numerical constructions | Error/residual bounds, conditioning, explicit degeneracy, tolerance propagation through projections/intersections | Exact degree-25 real roots, complete exact linear sets, certified analytic and spline/plane/sphere/cylinder intersections, rational spline curve/surface jets, exact curve/tensor-patch editing, curve/surface knot refinement/removal and exact edited-curve intersection results exist; general degree elevation passed independent complete-grid/native, sustained fuzz and platform gates at 27647fad; arbitrary face trimming, intersections and topology tolerance propagation pending |
| Geometric invariants | Analytic cases, conservation, covariance, winding independence, positive mass/inertia, classification consistency | 256 generated prisms, 2,417 rational 2D fixtures, 1,648 spatial predicate fixtures, 963 linear intersections, 448 quadratic roots, 1,092 curved intersections, 1,059 curve and 791 surface spline evaluations, 95 general roots, 284 spline/plane, 118 spline/quadric, 554 exact proximity cases, 684 complete linear sets, 636 exact curve edits, 745 exact surface edits, 723 exact curve knot edits, 1,756 exact surface knot edits, 274 arbitrary rational intersection cases and 10,000 integer cases run in Cargo tests |
| Structural validity | Closed oriented shells, curve/pcurve/surface consistency, self-intersection, disconnected regions, cavities and invalid-import diagnosis | Specialized prism checks exist; generic B-rep validation/healing pending |
| Stable topology history | Generated/modified/deleted face and edge mappings across kernel operations; explicit split/merge ambiguity | Only body-local IDs and extrusion provenance exist; general operation history pending |
| Fuzzing and differential geometry | Structured operation sequences, invalid input, minimized failures, independent oracles and original upstream tests | Sixteen sanitizer/coverage-guided targets (all passed at 27647fad) with daily campaigns and retained corpora, 66 OCCT prisms, 72 native line/plane cases, 174 polynomial/curved cases, 456 curve, 210 surface, 214 spline/plane, 92 spline/quadric, 370 proximity, 433 complete linear-set, 546 Bézier curve, 691 tensor-patch, 667 curve knot-editing and 1,748 surface knot-editing cases (reviewed numerical/contract differences and nonresults reported separately) and three original DRAW cases; a separate source-pinned degree-elevation runtime oracle passes with 187 matches and 80 reviewed discrepancies; full algorithm coverage remains pending |
| Source and behavior traceability | Pinned source symbols, supported domains, tolerance/error contracts, reviewed divergences | `SOURCE_MAP.md` records OCCT references and independently implemented mathematical algorithms |
| Failure containment | Atomic operations, bounded work/memory, cancellation, useful error context, no silent approximation | Immutable solid construction, typed errors and profile size limits exist; long-operation cancellation pending |
| Interchange | STEP unit/topology round trips, independent reader checks, watertight deflection-controlled export meshes | Not implemented |
| Platform/performance | Native and WASM execution, concurrency determinism, kernel workloads with time/memory budgets | Native debug/release Cargo CI and WASM compile check exist; WASM runtime and benchmark gates pending |

## Work in priority order

**Establish mathematical contracts first.** Separate exact combinatorial
decisions from approximate geometric construction and from modeling tolerance.
Extend the existing spatial predicates and certified analytic intersections with
distance/comparison predicates as needed. Define parameter domains, units,
residual/error bounds and degeneracy behavior for curve evaluation, projection,
root finding and intersections before building Booleans on top. Never increase
a tolerance merely to hide an unstable calculation. Exact predicates do not
make constructed coordinates exact.

**Test the mathematics adversarially.** Generate valid profiles and operation
sequences, then perturb gaps, angles and sizes around degeneracies. Mix in
malformed inputs to test rejection. Check rigid-motion invariance, inverse
transforms, winding independence, conservation, positive mass/inertia, outward
normals and bounds containing geometry. Later, check Boolean identities and
volume conservation with explicit dimensional error budgets.

Use [structure-aware fuzzing](https://rust-fuzz.github.io/book/cargo-fuzz/structure-aware-fuzzing.html),
retain seeds and minimize every failure into a deterministic regression. Fuzz
for crashes, non-finite results, inconsistent topology and excessive work, as
well as disagreement with OCCT. [The fuzz workflow](FUZZING.md) now runs bounded
smoke campaigns and longer daily mutation against retained corpora; it is
separate from deterministic seeded tests. Crashes and oracle disagreements fail
CI and preserve reproducers. Minimization/diagnosis remains an explicit triage
step; broader algorithm coverage grows with the kernel.

**Make topology and operation history kernel contracts.** Define operation
results with generated, modified and deleted mappings before Booleans, fillets
and shelling. A split or merge need not have a unique surviving face identity;
return that relationship or ambiguity explicitly. Array positions, approximate
geometric equality and coincident-looking faces are not persistent identities.
Check manifoldness and curve/pcurve/surface agreement independently of the
algorithm that constructed the shape. Sample checks are useful regressions but
do not prove continuous surface agreement.

**Require independent references.** OCCT is a valuable reference, not a proof.
Keep exact arithmetic oracles, analytic cases and mathematical invariants beside
differential tests. Investigate disagreements against the geometric contract
and preserve reviewed differences with their rationale. As STEP and meshing
arrive, test units/topology round trips with another reader, approximation
error, shared-edge consistency and watertightness.

**Bound operational failure and cost.** Record inputs, kernel/source versions,
platform, tolerance and operation context. Benchmark standalone kernel sequences
and worst-case geometric arrangements with CPU/memory budgets. Add cancellation
and recovery checks when operations become long-running. Add WASM execution
tests, not just compilation, before claiming that runtime is supported. Define
per-operation budgets from measured workloads; do not promise generic speedups
over C++ or trade reliability for a faster average.

Continue with stronger B-rep invariants, topology history, certified curved
geometry and general intersections. Interchange and performance gates become
applicable as those capabilities arrive. A finite test suite and certified
analytic primitives do not make the entire kernel infallible.
