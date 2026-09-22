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
| Mathematical decisions | Correct predicate signs over declared domains, independent exact references, documented arithmetic assumptions | Exact finite-f64 2D orientation exists; 3D orientation, incircle/insphere and certified distance decisions pending |
| Numerical constructions | Error/residual bounds, conditioning, explicit degeneracy, tolerance propagation through projections/intersections | Basic tolerance guards exist; certified construction and general intersection work pending |
| Geometric invariants | Analytic cases, conservation, covariance, winding independence, positive mass/inertia, classification consistency | 256 generated prism cases, 2,417 rational predicate fixtures and 10,000 integer predicate cases run in Cargo tests |
| Structural validity | Closed oriented shells, curve/pcurve/surface consistency, self-intersection, disconnected regions, cavities and invalid-import diagnosis | Specialized prism checks exist; generic B-rep validation/healing pending |
| Stable topology history | Generated/modified/deleted face and edge mappings across kernel operations; explicit split/merge ambiguity | Only body-local IDs and extrusion provenance exist; general operation history pending |
| Fuzzing and differential geometry | Structured operation sequences, invalid input, minimized failures, independent oracles and original upstream tests | Deterministic generated tests, 66 OCCT prism cases and three original DRAW cases exist; sustained coverage-guided fuzzing and a source-pinned runtime oracle pending |
| Source and behavior traceability | Pinned source symbols, supported domains, tolerance/error contracts, reviewed divergences | `SOURCE_MAP.md` records OCCT references and independently implemented mathematical algorithms |
| Failure containment | Atomic operations, bounded work/memory, cancellation, useful error context, no silent approximation | Immutable solid construction, typed errors and profile size limits exist; long-operation cancellation pending |
| Interchange | STEP unit/topology round trips, independent reader checks, watertight deflection-controlled export meshes | Not implemented |
| Platform/performance | Native and WASM execution, concurrency determinism, kernel workloads with time/memory budgets | Native debug/release Cargo CI and WASM compile check exist; WASM runtime and benchmark gates pending |

## Work in priority order

**Establish mathematical contracts first.** Separate exact combinatorial
decisions from approximate geometric construction and from modeling tolerance.
Extend reliable orientation to the 3D decisions required by intersections, then
add distance/comparison predicates as needed. Define parameter domains, units,
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
well as disagreement with OCCT. The bounded seeded tests now run in CI; they are
not a coverage-guided fuzzer. Sustained fuzz campaigns and automatic failure
minimization remain to be implemented.

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

The next milestone is a documented, adversarially tested mathematical layer and
stronger B-rep invariants, followed by topology history and reliable
intersections. Interchange and performance gates become applicable as those
capabilities arrive. A finite test suite and one exact predicate do not make
the entire kernel infallible.
