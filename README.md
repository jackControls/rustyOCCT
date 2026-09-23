# rustyOCCT

A native Rust geometry kernel forked from
[Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT).

Prioritize mathematical correctness, numerical robustness, topology/history,
interchange and operational reliability. [noBS-CAD](https://github.com/jackControls/noBS-CAD)
provides a capability-scope reference for mechanical CAD, CAM, additive
manufacturing and analysis; kernel correctness is application-independent.
No renderer, windowing, GPU integration, viewer, or duplicate application framework.

**Status: first working kernel milestone, not a replacement for OCCT yet.**
The Rust implementation uses `num-bigint` and `num-rational` for exact arithmetic, has no C++
bindings, and forbids unsafe code in the kernel. The inherited C++ source remains available for reference and comparison;
Cargo does not build it. Rust work lives on the `rust-kernel` branch, while
`master` retains the upstream fork point.

## Implemented

- Validated simple concave/convex polygon and circle profiles with multiple
  disjoint polygonal or circular holes.
- Exact normal extrusions on arbitrary planes, with signed start/end offsets;
  convenience constructors for boxes and cylinders.
- Analytic lines, circles, planes and cylinders; shared vertices/edges, oriented
  face loops, and per-face parameter curves, including cylinder seams.
- Closed-shell connectivity validation, face provenance, and hole-aware Euler
  characteristic. Body-local indices are not persistent names across edits.
- Volume, surface area, centroid, central inertia tensor, exact bounds, point
  classification, rigid translation and rotation.
- Explicit tolerances, finite-input and coordinate-resolution checks, and errors
  for self-intersections, touching/nested holes and degenerate geometry.
- Exact 2D/3D orientation and sphere-membership predicates for finite `f64`
  inputs, including subnormals and values whose floating-point products overflow.
- Certified line–plane, segment–plane and segment–triangle intersections,
  including coplanar overlap. Decisions are exact; new coordinates and parameters
  have minimal binary64 enclosures. Rounded representatives retain those bounds.
- Complete closed-set intersections of all point/line/segment/plane/triangle
  pairings, including coplanar polygons and coincident infinite lines/planes.
  Canonical rational constructions survive finite floating-point overflow.
- Exact minimum distances and closest points across all point/line/segment/plane/
  triangle pairings. Rational witnesses and threshold comparisons remain available
  independently of floating-point output limits; non-unique minima are supported.
- Exact quadratic root classification and algebraic root comparisons; certified
  line/segment intersections with circles, spheres and infinite cylinders,
  including exact tangencies, segment clipping and preserved root multiplicity.
- Rational Bézier and B-spline curves and surfaces through degree 25, with
  certified positions and derivatives through order two, including mixed surface
  partials. Clamped, unclamped and periodic directions have explicit knot/seam sides.
- Exact rational Bézier extraction from those splines, including clipped and
  multi-period ranges; subdivision, trimming, reversal and degree elevation
  retain exact controls and original parameter units. See the
  [editing contract](rust/BEZIER_EDITING.md).
- Exact tensor Bézier patch extraction and editing, including independent U/V
  periodicity, rectangular trimming, subdivision, parameter exchange, degree
  elevation and exact boundary/isoparametric curves. Mixed partials and original
  parameter domains are retained. See [surface editing](rust/SURFACE_EDITING.md).
- Exact B-spline knot refinement and removal, including unclamped end controls
  and periodic seams. Rational knots and positive homogeneous controls survive
  edit sequences without rounding. See [knot editing](rust/KNOT_EDITING.md).
- Exact B-spline surface knot refinement and removal in both directions, with
  independent periodicity, atomic grid edits, exact partials, rational patch
  extraction and retained isocurves. See [surface knot editing](rust/SURFACE_KNOT_EDITING.md).
- Exact real-root isolation through degree 25 and polynomial signs at algebraic
  roots. Certified spline intersections with planes, spheres and infinite cylinders;
  crossings, tangencies and maximal overlap intervals
  retain distinct parameters even when their floating enclosures coincide. Explicit
  parameter ranges include periodic seam crossings and multiple turns, with
  exact clipping and traversal budgets. Quadric substitution supports equations
  through degree 50.
- Coverage-guided fuzzing of predicates, intersections, splines and modeling sequences,
  with independent mathematical oracles and retained corpora/failure inputs.

```sh
cargo test --workspace --locked
cargo run --locked --example plate
cargo check --workspace --lib --target wasm32-unknown-unknown
```

The plate example constructs an 80 x 40 x 6 mm plate with four exact circular
through holes and prints its geometry properties. It is entirely headless.
Rust 1.85 or newer is required; the optional WebAssembly check requires that
target to be installed.

## Scope and next work

[The kernel scope and porting plan](rust/PORTING.md) maps needed capabilities to
OCCT families. [The mathematical foundation](rust/MATHEMATICS.md) is the current
priority, followed by topology/history and general curve/surface intersections. General
Booleans, arc profiles, spline topology, revolutions, sweeps/lofts, fillets/chamfers, shelling,
modeled threads, STEP, drawing HLR, and tessellation remain to be implemented.

Tessellation and hidden-line geometry belong in the kernel because noBS-CAD
needs export meshes and technical drawings. The application continues to own
rendering, sketches, feature history, assemblies, CAM planning, and document I/O.

Application adapters and migration are deferred. Kernel acceptance uses its own
contracts, independent mathematical references, geometric invariants, adversarial
inputs and OCCT differential tests.

## Validation

The mathematical tests include **2,417 exact rational 2D orientation fixtures**
(all six permutations), **1,648 spatial predicate fixtures**, **963 certified
linear-intersection fixtures**, **448 quadratic-root fixtures**, **1,092 curved
intersection fixtures**, **1,059 curve and 791 surface spline fixtures**, **95 general
root, 284 spline–plane and 118 spline–quadric fixtures**, **554 exact proximity fixtures**, **684 complete linear-set fixtures**, **636 exact Bézier curve, 745 tensor-patch, 723 curve knot-editing and 1,756 surface knot-editing fixtures**,
**10,000 generated integer predicate cases**, and
**256 generated prism invariant cases**. Debug and optimized native builds run
the same checks. See [the numerical contracts](rust/MATHEMATICS.md) for limits.

The checked-in OCCT 7.9.3 reference corpus covers **66 solids and 2,292 point
classifications**, comparing volume, area, centroid, bounds, inertia and topology
counts. Ordinary Cargo tests run this corpus without an OCCT SDK.

Native OCCT also checks 72 line–plane cases, 174 polynomial/curved cases,
456 spline-curve cases, 210 spline-surface cases, 214 spline–plane and
92 spline–quadric cases. [Quadric differences](rust/NATIVE_SPLINE_QUADRIC_DIVERGENCES.md) and
Spline–plane [contract and numerical differences](rust/NATIVE_SPLINE_PLANE_DIVERGENCES.md)
are independently checked and version-pinned. High-degree native numerical
differences are [reviewed separately](rust/NATIVE_SPLINE_DIVERGENCES.md); they
are not counted as parity matches.
[Sustained fuzzing](rust/FUZZING.md)
runs sixteen instrumented targets on pushes/PRs and daily, restoring the evolving
corpus and retaining crashes, timeouts and mathematical disagreements.

[Exact edited-curve intersections](rust/EXACT_SPLINE_INTERSECTIONS.md) preserve
rational knots, homogeneous controls and trim endpoints through plane, sphere
and cylinder queries. An additional **274 exact rational fixtures** check
complete contacts and overlaps, including values outside binary64 range. All
402 existing spline intersection fixtures also run after exact conversion,
refinement and removal; the 306 native cases retain their original input count
and separately verify all four representations.

Another 370 inputs compare linear-set minimum distances with two native OCCT
APIs and independently certify Rust's returned witnesses. Unbounded native
nonresults and affine parallelism differences are [recorded separately](rust/NATIVE_PROXIMITY_DIVERGENCES.md).

A further 433 inputs compare complete linear intersections with native COMMON
and SECTION; [intentional dimension filtering, unbounded nonresults and three
extra-edge observations](rust/NATIVE_LINEAR_INTERSECTION_DIVERGENCES.md) are counted separately.

Another 546 inputs compare spline extraction and Bézier editing, including
complete control data. [Degree-25 native editing differences](rust/NATIVE_BEZIER_EDITING_DIVERGENCES.md)
are retained separately; Rust must match independent exact polynomials first.

A further 691 inputs compare tensor patch extraction, edits and isocurves.
All complete homogeneous tensor coefficients are checked independently;
[native high-degree differences](rust/NATIVE_SURFACE_EDITING_DIVERGENCES.md)
remain separate from matches.

Another 667 inputs compare complete knot refinement/removal representations and
operation success flags. Independent coefficient equations check all Rust
results; [native round-trip and seam correspondence differences](rust/NATIVE_KNOT_EDITING_DIVERGENCES.md)
remain separate from matches.

Another 1,748 inputs compare complete surface knot grids, domains and operation
flags. All Rust results must match independent coefficient equations; OCCT
round-trip differences are [reviewed separately](rust/NATIVE_SURFACE_KNOT_DIVERGENCES.md).

```sh
# Optional live differential comparison against an installed OCCT SDK:
python3 rust/tools/compare_occt.py --occt-root /path/to/occt
```

See [validation details and limits](rust/VALIDATION.md). CI checks Rust on Linux,
macOS and Windows, plus the minimum Rust version and a WebAssembly library build.

The [upstream test bridge](rust/UPSTREAM_TESTS.md) runs unchanged OCCT DRAW tests
and assertions against Rust and native OCCT. Three original AABB regressions
currently pass on both; unsupported commands and missing data are explicit
non-passes. CI runs both backends and preserves JSON/JUnit results and logs.

```sh
python3 rust/tools/run_upstream_tests.py --backend both --draw-exe /path/to/DRAWEXE
```

The [source map](rust/SOURCE_MAP.md) records implementation references and
deliberate limitations. [Production release gates](rust/PRODUCTION_READINESS.md)
cover mathematical correctness, topology/history, numerical robustness, fuzzing,
interchange, failure containment and performance.

## Upstream and license

Fork point: `Open-Cascade-SAS/OCCT` commit
`3d097a0328e71b826377d4814ab05ec3c3d23871`. The
[original upstream README](rust/UPSTREAM_README.md) is preserved.

This fork retains the [GNU LGPL 2.1](LICENSE_LGPL_21.txt) with the
[OCCT exception](OCCT_LGPL_EXCEPTION.txt). The Rust additions use the same
license expression: `LGPL-2.1-only WITH OCCT-exception-1.0`. Existing upstream
copyright and license notices remain in place.
