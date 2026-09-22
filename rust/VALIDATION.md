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
448 quadratic-root fixtures, 1,092 curved-intersection fixtures, 1,059 curve and 791 surface spline fixtures,
95 general polynomial-root, 284 spline/plane, 118 spline/quadric, 554 proximity and 684 complete linear-set fixtures,
10,000 integer-oracle predicate cases and 256 generated prism invariant cases.
These tests do not depend on OCCT or an application and run in native debug and
release CI. These deterministic generated tests are separate from the ten
[coverage-guided fuzz targets and daily retained-corpus campaigns](FUZZING.md).

`compare_intersections.py` executes native OCCT `IntAna_IntConicQuad` and Rust's
line/plane primitive on 72 shared well-conditioned cases, comparing intersection
type, affine parameter and point coordinates. The test budget is
`1e-10 + 2e-12*abs(expected)`. Local OCCT 7.9.3 used at most 0.000128 of that
budget. CI runs the same comparison against its recorded distribution runtime.
Exact parallelism intentionally differs from OCCT's angular-tolerance policy;
near-degenerate/extreme cases use the independent exact oracles instead.

`compare_proximity.py` captures 370 point/line/segment/plane/triangle inputs
through native B-rep and lower-level affine distance APIs. Exact rational Rust
witnesses pass independent analytic formulas and supporting-plane certificates
before native comparison. Local OCCT 7.9.3 gives 316 B-rep matches, 54 B-rep
nonresults, and 120 affine matches with six numerical differences; those two
observation sets overlap. All 412 native witness pairs are also checked for
shape ownership and distance. The [reviewed differences](NATIVE_PROXIMITY_DIVERGENCES.md)
are version/input/value pinned, and never weaken the independent Rust checks.
Linux OCCT 7.6.3 has the same 316 B-rep matches and 54 separately reviewed
nonresults; all 126 affine distances match. Rust's complete rational results
are byte-identical between the two platforms.
The 554 Fraction fixtures add full-exponent, near-degenerate and unrepresentable
cases outside this well-scaled native corpus. No general B-rep distance or
unchanged upstream extrema-test coverage is implied.

`compare_linear_sets.py` checks all 25 ordered linear primitive pairings over
433 native inputs. The initial native geometries were captured before Rust
implementation. All exact results equal an independent boundary-crossing oracle.
Local OCCT 7.9.3 yields 384 complete combined COMMON/SECTION matches, 46 empty
results for unbounded intersections, and three extra-edge observations. The
[reviews](NATIVE_LINEAR_INTERSECTION_DIVERGENCES.md) pin each native result;
COMMON's intentional omission of lower-dimensional contacts is counted separately.
The 684 exact fixtures include five/six-vertex polygons, full-exponent and
subnormal geometry, invalid inputs and unrepresentable constructions. Rust
and Python boundary oracles differ from production affine/halfspace enumeration.
Every defining-point permutation, operand order and minimal construction bound
is checked, and `linear_sets` continuously mutates the corpus. General face
splitting, B-rep Boolean operations and topology/history remain pending.

`compare_curved.py` first captures independent native observations, then checks
Rust against 174 cases of polynomial roots, line/sphere, line/cylinder and
coplanar XY line/circle intersections. Repeated polynomial roots retain their
multiplicity; exactly duplicate native hit parameters are merged only when
comparing geometric hit counts. Cases include exact tangencies and degree
reduction. Local OCCT 7.9.3 used at most 0.000078 of the same comparison budget.
Its source/runtime versions remain separate. Full-exponent, near-tangent,
arbitrary-plane circle and closed-segment behavior is tested by independent
polynomial-sign/axial-projection oracles, not inferred from this native corpus.

`compare_splines.py` captures 456 native OCCT observations before evaluating
Rust positions and first/second derivatives. It covers degree 1..25,
rational/polynomial Bézier curves, clamped/unclamped/periodic B-splines and
explicit left/right derivatives at repeated knots and seams. The original
`Geom_BSplineCurve_Test.cxx` cubic fixture is included. This is a native API
comparison, not execution of those original C++ GoogleTest assertions.
`compare_surfaces.py` adds 210 native rational surface jets, including mixed
partials, all U/V periodicity combinations and degree 25 in both directions.

The 1,850 independent Python `Fraction` fixtures use basis-function recursion
and closed quotient formulas; production differentiates homogeneous pole
interpolation and uses recursive quotient rules. They cover subnormal spans,
huge/tiny weights, extreme periodic wrapping, position versus derivative
overflow and exact side/quadrant continuity decisions. All preexisting 385
curve fixture rows were preserved byte-for-byte when expanding the corpus.

The enlarged native corpus exposed [degree-25 numerical discrepancies](NATIVE_SPLINE_DIVERGENCES.md).
On local OCCT 7.9.3, curves have 427 matches and 29 reviewed divergent cases;
surfaces have 169 matches and 41 reviewed divergent cases. Each review pins
version, input and the complete native jet, and recomputes the independent exact
jet before accepting the Rust result. The budget remains `1e-10 + 2e-12*abs(expected)`.
These observations do not establish that OCCT promises that error bound.

On Linux OCCT 7.6.3, the complete curve corpus has 426 matches and 30 reviewed
differences (the 29 periodic cases plus the original case below); surfaces have
169 matches and 41 reviewed differences. Rust observations are identical on
macOS and Linux. Native numerical differences remain explicit and version-pinned.

Ubuntu's OCCT 7.6.3 has one **reviewed numerical divergence**, not a matching
case: `d25_r1_s3_5`, second derivative X, is `-0.47077685134240305` instead of
the exact value enclosed by `[-0.47077685157054794, -0.4707768515705479]`.
The difference is 2.2602 times the unchanged native comparison budget; OCCT
7.9.3 stays within that budget. The independent basis/quotient calculation
agrees with Rust. This records an observed numerical difference, not a claim
that OCCT promises our error bound.

[`occt-spline-divergences.json`](fixtures/occt-spline-divergences.json) pins the
input hash, native version/value bits, exact rational answer and original CI
evidence. The comparison recomputes the independent rational answer and checks
Rust against its minimal enclosure before accepting that specific review.
Changed inputs, native results, versions or incorrect Rust outputs still fail.
Reports count matched, reviewed and unexpected cases separately. In the original
214-case subset, OCCT 7.6.3 has **213 matching, 1 reviewed divergence, 0 unexpected**. Use
`compare_splines.py --strict-native` to fail on reviewed differences too. The
ordinary exact tests and fuzz comparisons have no such exception.

`compare_spline_plane.py` checks 214 native `GeomAPI_IntCS` observations against
complete independent exact certificates. Local OCCT 7.9.3 has 155 matching cases
and 59 reviewed differences. Extra near-tangent native points, missing overlap
intervals, periodic endpoint conventions and numerical discrepancies are
[documented separately](NATIVE_SPLINE_PLANE_DIVERGENCES.md). Every review pins
version, input, complete native values and independently recomputed certificate;
it never exempts a Rust result from exact checks. `--strict-native` fails on
reviewed differences too. SymPy is hash-pinned and test-only; Cargo tests use
checked-in exact certificates without requiring Python or native OCCT.
Linux OCCT 7.6.3 has 158 matches and 56 separately reviewed differences in the
combined 214-case corpus. All Rust certificates agree across macOS and Linux.
The additional 82 cases exercise explicit parameter ranges, including shifted
periods, with native observations reviewed independently for each version.

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
python3 rust/tools/generate_proximity_fixtures.py --check
python3 rust/tools/generate_linear_sets_fixtures.py --check
python3 rust/tools/generate_curved_fixtures.py --check
python3 rust/tools/generate_spline_fixtures.py --check
python3 rust/tools/generate_surface_fixtures.py --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo check --workspace --lib --locked --target wasm32-unknown-unknown
cargo run --locked --example plate
```

The live comparison requires a C++17 compiler and OCCT's modeling SDK:

```sh
python3 rust/tools/compare_occt.py --occt-root /path/to/occt
python3 rust/tools/compare_intersections.py --occt-root /path/to/occt
python3 rust/tools/compare_proximity.py --occt-root /path/to/occt
python3 rust/tools/compare_linear_sets.py --occt-root /path/to/occt
python3 rust/tools/compare_curved.py --occt-root /path/to/occt
python3 rust/tools/compare_splines.py --occt-root /path/to/occt
python3 rust/tools/compare_surfaces.py --occt-root /path/to/occt
```

On macOS the default SDK is `/opt/homebrew/opt/opencascade`; on Linux it is
`/usr`. `OCCT_ROOT` and `CXX` can override discovery. This helper currently
supports Unix-style SDK layouts; Windows still runs all ordinary Rust and
recorded-corpus tests. It does not install software or change noBS-CAD.

Results and optional native executables go in ignored `target/occt-oracle/`
and `target/intersection-oracle/`, `target/curved-oracle/`, `target/spline-oracle/`, `target/surface-oracle/`.
The live oracle does not link any rendering or data-exchange toolkit.

To intentionally refresh fixture inputs or reference data:

```sh
python3 rust/tools/generate_fixtures.py
python3 rust/tools/compare_occt.py --occt-root /path/to/occt --write-baseline
```

Review all changes to fixture inputs, expected results, version metadata and
tolerances together. Never update a baseline merely to make a failing test pass.
The generator uses a fixed seed and writes inputs only.

## Spline/quadric comparisons

`compare_spline_quadric.py` captures 92 sphere/cylinder observations with the
same headless `GeomAPI_IntCS` helper. OCCT 7.9.3 and Linux OCCT 7.6.3 each have
66 direct matches and 26
[reviewed differences](NATIVE_SPLINE_QUADRIC_DIVERGENCES.md), primarily missing
overlap intervals and later periodic events. Every Rust result is independently
recomputed, including complete point counts, minimal enclosures, contact orders
and overlap endpoints. Review hashes and native budgets use the same fail-closed
checks as the plane bridge. The ordinary 118 exact fixtures also exercise dense
degree-25 equations, order-50 contacts and extreme radii/axis magnitudes.

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
