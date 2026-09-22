# rustyOCCT

A native Rust geometry kernel for [noBS-CAD](https://github.com/jackControls/noBS-CAD),
forked from [Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT).

Port the modeling and geometry capabilities noBS-CAD needs today and for its
accepted mechanical CAD, CAM, additive manufacturing, and analysis directions.
No renderer, windowing, GPU integration, viewer, or duplicate application framework.

**Status: first working kernel milestone, not a replacement for OCCT yet.**
The Rust implementation has no dependencies, no C++ bindings, and forbids unsafe
code. The inherited C++ source remains available for reference and comparison;
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

[The noBS-CAD porting plan](rust/PORTING.md) maps every current kernel job and
query to the required OCCT families and a staged Rust implementation. General
Booleans, arcs/B-splines, revolutions, sweeps/lofts, fillets/chamfers, shelling,
modeled threads, STEP, drawing HLR, and tessellation remain to be implemented.

Tessellation and hidden-line geometry belong in the kernel because noBS-CAD
needs export meshes and technical drawings. The application continues to own
rendering, sketches, feature history, assemblies, CAM planning, and document I/O.

The noBS-CAD application has not been switched to this kernel. Migration requires
operation-by-operation comparison and saved-project replay, including selected
face/edge identity; successful primitive construction alone is insufficient.

## Validation

The checked-in OCCT 7.9.3 reference corpus covers **66 solids and 2,292 point
classifications**, comparing volume, area, centroid, bounds, inertia and topology
counts. Ordinary Cargo tests run this corpus without an OCCT SDK.

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
cover application replay, stable selections, numerical robustness, fuzzing,
interchange, failure containment and performance.

## Upstream and license

Fork point: `Open-Cascade-SAS/OCCT` commit
`3d097a0328e71b826377d4814ab05ec3c3d23871`. The
[original upstream README](rust/UPSTREAM_README.md) is preserved.

This fork retains the [GNU LGPL 2.1](LICENSE_LGPL_21.txt) with the
[OCCT exception](OCCT_LGPL_EXCEPTION.txt). The Rust additions use the same
license expression: `LGPL-2.1-only WITH OCCT-exception-1.0`. Existing upstream
copyright and license notices remain in place.
