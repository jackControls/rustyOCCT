# Kernel validation

For original test reuse, see [the DRAW bridge](UPSTREAM_TESTS.md): unchanged Tcl
tests, assertion helpers and result rules execute against both kernels. Its
three passing upstream AABB cases and explicit capability/data gaps are separate
from the 66-solid corpus below. Implementation provenance is in [SOURCE_MAP.md](SOURCE_MAP.md);
remaining release requirements are in [PRODUCTION_READINESS.md](PRODUCTION_READINESS.md).

## Current evidence

The [mathematical foundation](MATHEMATICS.md) includes exact 2D/3D orientation
and insphere, 2,417 independent 2D rational fixtures in all six permutations,
1,648 spatial predicate fixtures, 963 certified linear-intersection fixtures,
448 quadratic-root fixtures, 1,092 curved-intersection fixtures,
10,000 integer-oracle predicate cases and 256 generated prism invariant cases.
These tests do not depend on OCCT or an application and run in native debug and
release CI. These deterministic generated tests are separate from the four
[coverage-guided fuzz targets and daily retained-corpus campaigns](FUZZING.md).

`compare_intersections.py` executes native OCCT `IntAna_IntConicQuad` and Rust's
line/plane primitive on 72 shared well-conditioned cases, comparing intersection
type, affine parameter and point coordinates. The test budget is
`1e-10 + 2e-12*abs(expected)`. Local OCCT 7.9.3 used at most 0.000128 of that
budget. CI runs the same comparison against its recorded distribution runtime.
Exact parallelism intentionally differs from OCCT's angular-tolerance policy;
near-degenerate/extreme cases use the independent exact oracles instead.

`compare_curved.py` first captures independent native observations, then checks
Rust against 174 cases of polynomial roots, line/sphere, line/cylinder and
coplanar XY line/circle intersections. Repeated polynomial roots retain their
multiplicity; exactly duplicate native hit parameters are merged only when
comparing geometric hit counts. Cases include exact tangencies and degree
reduction. Local OCCT 7.9.3 used at most 0.000078 of the same comparison budget.
Its source/runtime versions remain separate. Full-exponent, near-tangent,
arbitrary-plane circle and closed-segment behavior is tested by independent
polynomial-sign/axial-projection oracles, not inferred from this native corpus.

The first implementation has deterministic analytic tests plus a recorded,
independently evaluated OCCT 7.9.3 corpus. The live comparison was run on Apple
Silicon macOS with Rust 1.96.0 and the installed OCCT 7.9.3 SDK.

`fixtures/prisms.txt` contains 66 solids in 11 families: boxes, triangles,
concave profiles, polygonal frames, cylinders, tubes, four-hole plates, mixed
polygon/circle holes, round profiles with square holes, collinear boundary
segments, and concave profiles with circular holes. Each family has six
variations covering scale, reversed profile winding, reversed offsets, arbitrary
plane orientation, spatial translation and large local coordinates.

The independent OCCT executable constructs real B-reps with MakeWire/MakeFace/
MakePrism, verifies them with `BRepCheck_Analyzer`, and evaluates:

- volume, surface area and centroid;
- exact geometric axis-aligned bounds without tessellation;
- central inertia tensor in world axes;
- unique vertex, edge and face counts;
- 2,292 inside/boundary/outside point classifications, including caps, seams,
  hole rims, random interior/exterior samples and out-of-height samples.

Both executables consume the same input file. The expected values in
`fixtures/occt-baseline.tsv` come from OCCT, never from the Rust result.
`fixtures/occt-baseline.json` records the oracle version and comparison summary.

The numerical acceptance budget is `absolute + 2e-9 * abs(expected)`, with
absolute `1e-7` for volume/area/coordinates and `1e-6` for inertia. Topology
counts and classifications must match exactly. This is a test budget for this
corpus, not a universal geometric error guarantee. The greatest observed error
was approximately 0.000255 of that budget.

OCCT's default global-origin inertia accumulation lost precision on translated
small parts. The oracle therefore centers a geometry copy before calculating
the *central* tensor; it uses the original solid for all other queries. This
does not change central inertia or relax the comparison tolerance. Separate
analytic box/tube tests and translation-invariance tests check the Rust tensor.

## Reproduce

The usual tests require only Rust:

```sh
cargo fmt --all --check
cargo test --workspace --locked
cargo test --workspace --locked --release
python3 rust/tools/generate_predicate_fixtures.py --check
python3 rust/tools/generate_spatial_fixtures.py --check
python3 rust/tools/generate_curved_fixtures.py --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo check --workspace --lib --locked --target wasm32-unknown-unknown
cargo run --locked --example plate
```

The live comparison requires a C++17 compiler and OCCT's modeling SDK:

```sh
python3 rust/tools/compare_occt.py --occt-root /path/to/occt
python3 rust/tools/compare_intersections.py --occt-root /path/to/occt
python3 rust/tools/compare_curved.py --occt-root /path/to/occt
```

On macOS the default SDK is `/opt/homebrew/opt/opencascade`; on Linux it is
`/usr`. `OCCT_ROOT` and `CXX` can override discovery. This helper currently
supports Unix-style SDK layouts; Windows still runs all ordinary Rust and
recorded-corpus tests. It does not install software or change noBS-CAD.

Results and optional native executables go in ignored `target/occt-oracle/`
and `target/intersection-oracle/`, `target/curved-oracle/`.
The live oracle does not link any rendering or data-exchange toolkit.

To intentionally refresh fixture inputs or reference data:

```sh
python3 rust/tools/generate_fixtures.py
python3 rust/tools/compare_occt.py --occt-root /path/to/occt --write-baseline
```

Review all changes to fixture inputs, expected results, version metadata and
tolerances together. Never update a baseline merely to make a failing test pass.
The generator uses a fixed seed and writes inputs only.

## Representation and numeric limits

- Only finite normal extrusions of simple polygon/circle profiles are supported.
  No partial arcs, freeform surfaces, taper, oblique extrusion, general trimming,
  Boolean edits, closed cavities or multiple disconnected regions per solid.
- Holes must be strictly inside, non-touching, disjoint and non-nested. Material
  islands can become separate solids; this crate does not discover sketch loops.
- Linear tolerance defaults to `1e-7 mm`, angular tolerance to `1e-12 rad`.
  Callers may supply validated tolerances. Inputs are revalidated if boundaries
  built under one tolerance are assembled with another.
- Coordinates that cannot resolve the requested tolerance are rejected using
  a 16-ULP-scale budget; geometry is not rescaled to bypass this check.
  This is a conservative solid-construction guard. The separate exact
  predicates accept all finite `f64` inputs, and implemented linear/curved
  intersections use certified coordinate/parameter enclosures. Solid distance, area and moment
  calculations remain ordinary floating-point arithmetic.
- Boundary validation is quadratic; input is bounded to 4,096 total profile
  edges and 128 holes. This is an initial implementation limit.
- Topology is immutable and body-local. Indices and `FaceOrigin` do not claim
  persistent naming through arbitrary feature edits or Booleans.
- The topology validator checks connectivity, opposite uses and analytic curve/
  surface agreement at fixed parameter samples. Validity also depends on the
  profile builders. It is not a general B-rep validator, healer or importer.
- Mass properties are geometric at unit density; density and material metadata
  belong to the application. No strength or manufacturing result is implied.

These comparisons establish the implemented prism subset. They do not establish
parity with the fork-point OCCT development revision or all of OCCT. Kernel
readiness requires the mathematical, topology/history, interchange and
operational gates in `PRODUCTION_READINESS.md`. Application migration is
separate, deferred work.
