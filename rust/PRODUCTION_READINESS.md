# Production readiness for noBS-CAD

The current kernel is an early, validated prism implementation. It is not yet
ready to replace OCCT in production. Upstream regression compatibility is one
layer of evidence; it cannot establish correctness for all inputs, preserve
application selections by itself, or justify copying a known upstream defect.

## Release gates, by capability

Promote individual noBS-CAD operations only after all applicable gates pass.
Use `PORTING.md` for the actual required jobs, queries and future capabilities.

| Gate | Evidence required | Current state |
| --- | --- | --- |
| Source and behavior traceability | Pinned source symbols, supported domains, tolerance/error contracts, reviewed divergences | `SOURCE_MAP.md` and initial implementation tests exist; expand for every operation |
| Differential geometry | Original upstream tests plus tighter independent volume, area, centroid, inertia, bounds, classification, geometry and topology comparisons | 66 prism cases and first three original DRAW tests; exact application-version oracle build still needed |
| Application replay | Saved noBS-CAD documents, feature edits/suppression/reordering, transformed instances, export and deterministic replay | Not implemented; required before switching any application feature |
| Stable selections | Generated/modified/deleted face and edge mappings; edits preserve downstream references or return an explicit broken-reference error | Only body-local IDs and extrusion provenance exist |
| Numerical robustness | Reviewed tolerance propagation, periodic seams, tangency, coincident faces, slivers, tiny/huge scales; robust predicates where comparisons require them | Basic tolerance guards and boundary cases exist; general intersection/Boolean work pending |
| Structural validity | Closed oriented shells, curve/pcurve/surface consistency, self-intersection, disconnected regions, cavities and invalid-import diagnosis | Specialized prism checks exist; generic B-rep validation/healing pending |
| Failure containment | Atomic operations, bounded work/memory, cancellation, useful error context, no silent approximation | Immutable solid construction, typed errors and profile size limits exist; long-operation cancellation and adapter transaction semantics pending |
| Interchange | STEP unit/topology round trips, trusted independent reader checks, watertight deflection-controlled export meshes | Not implemented |
| Platform/performance | Native and WASM runtime replay, concurrency determinism, realistic operation time/memory regressions | Native Cargo CI and WASM compile check exist; WASM runtime and benchmark gates pending |

## Highest-value next work

**First, make application replay the acceptance test.** Build a headless adapter
for one supported New Body extrusion job. Feed identical real noBS-CAD job DTOs
to OCCT and Rust and compare resulting measurements, face/edge geometry,
classification, selections and feature errors. Save every disagreement as a
small permanent fixture. Run this comparison in tests/development before opting
an application operation into Rust. Keep switching explicit; a hidden fallback
to OCCT would conceal missing Rust behavior.

**Add structured fuzzing and mathematical invariants.** Generate valid profiles
and operation sequences so tests reach geometric algorithms, then perturb gaps,
angles and sizes around degeneracies. Mix in malformed inputs to test rejection.
Check rigid-motion invariance, inverse transforms, winding independence, positive
mass, outward normals and bounds containing sampled geometry. Later, check
Boolean identities and volume conservation with explicit tolerance budgets.
Use [structure-aware fuzzing](https://rust-fuzz.github.io/book/cargo-fuzz/structure-aware-fuzzing.html),
retain seeds and minimize every failure into a deterministic regression. Fuzz
for crashes, non-finite results, inconsistent topology and excessive work, as
well as disagreement with OCCT. This campaign is proposed, not yet running.

**Treat topology history and tolerance as architecture.** Stable feature
references and explicit tolerance propagation are harder to retrofit than
another primitive constructor. Define operation results with generated,
modified and deleted mappings before Booleans, fillets and shelling. Avoid
equality based on array order or transient face numbers. Track accumulated
tolerance through intersections; do not inflate tolerances until a test passes.

**Require independent checks.** OCCT is a valuable reference, not a proof.
Keep analytic cases and mathematical invariants alongside differential tests.
Where OCCT and Rust disagree, investigate which behavior satisfies the geometric
contract. Preserve reviewed differences with their rationale. For interchange,
use another STEP reader and inspect exact geometry as well as tessellated output.

**Make performance and failures reproducible.** Record inputs, kernel/source
versions, platform, tolerance and operation context. Benchmark realistic part
histories and worst-case geometric arrangements with bounded CPU/memory budgets.
Include cancellation and recovery tests once operations become long-running.
Add WASM execution tests, not just successful compilation, before a browser path
is considered supported. Define per-operation release budgets from measured
workloads rather than promising a generic speedup over C++.

The next milestone should be one production-quality extrusion path with replay
and stable references, followed by the geometry/intersection foundation needed
for exact Booleans. Broad feature coverage without these gates would not make a
dependable CAD kernel.
