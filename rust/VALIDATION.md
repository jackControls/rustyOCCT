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
95 general polynomial-root, 284 spline/plane, 118 spline/quadric, 554 proximity, 684 complete linear-set, 636 exact Bézier curve, 745 tensor-patch and 723 knot-editing fixtures,
10,000 integer-oracle predicate cases and 256 generated prism invariant cases.
These tests do not depend on OCCT or an application and run in native debug and
release CI. These deterministic generated tests are separate from the seventeen
[coverage-guided fuzz targets and daily retained-corpus campaigns](FUZZING.md).

Exact root-to-root ordering was accepted at revision `de1d1d43`: 122 debug and
122 release tests on Linux/macOS/Windows, Rust 1.85, WASM compilation and all
sixteen then-existing fuzz targets passed. Its local roots campaign completed
600.06 seconds of mutation with 10,483 mutations. Full kernel parity remains
incomplete. The new [point-to-spline minimum-set capability](SPLINE_PROXIMITY.md)
adds 30 independent complete-minimum cases, 94 whole image-polynomial equations,
a source-pinned native bridge and a seventeenth fuzz target; its clean-revision
platform/native/fuzz acceptance is pending.

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
Local OCCT 7.9.3 and Linux OCCT 7.6.3 each yield 384 complete combined
COMMON/SECTION matches, 46 empty results for unbounded intersections, and three
extra-edge/off-line-point observations. All 433 Rust result rows are identical
between platforms. The
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
python3 rust/tools/generate_bezier_editing_fixtures.py --check
python3 rust/tools/generate_surface_editing_fixtures.py --check
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
python3 rust/tools/compare_bezier_editing.py --occt-root /path/to/occt
python3 rust/tools/compare_surface_editing.py --occt-root /path/to/occt
python3 rust/tools/compare_surfaces.py --occt-root /path/to/occt
```

On macOS the default SDK is `/opt/homebrew/opt/opencascade`; on Linux it is
`/usr`. `OCCT_ROOT` and `CXX` can override discovery. This helper currently
supports Unix-style SDK layouts; Windows still runs all ordinary Rust and
recorded-corpus tests. It does not install software or change noBS-CAD.

Results and optional native executables go in ignored `target/occt-oracle/`
and `target/intersection-oracle/`, `target/curved-oracle/`, `target/spline-oracle/`, `target/surface-oracle/`, `target/bezier-editing-oracle/`, `target/surface-editing-oracle/`.
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

## Exact Bézier extraction and editing

The bridge captures 546 native inputs before running Rust and compares all
1,672 returned arcs' control data, domains and degrees. Independent exact
Cox basis polynomials plus affine substitutions must reproduce every Rust
homogeneous control exactly. The native comparison permits common weight
normalization and the documented floating budget. OCCT 7.9.3 has 506 matching
cases and 40 [reviewed degree-25 edit differences](NATIVE_BEZIER_EDITING_DIVERGENCES.md);
Linux OCCT 7.6.3 has 508 matches and 38 separately pinned differences. All 546
complete exact Rust outputs are byte-identical between the platforms.
Source GoogleTest inputs are reused, not executed unchanged.

The 636 exact fixtures add subnormal/adjacent knots, full-exponent geometry and
weights, multi-period intervals and shifted knots below a floating ULP. Separate
tests retain cuts 2^-2048 apart, prove edit commutation, check distinct one-sided
derivatives, invalid rational parameters and traversal exhaustion. Fuzzing
compares complete coefficients as well as exact jets and minimal bounds.
No fixture or review exempts Rust from the exact mathematical contract.

## Exact edited-curve intersections

All 402 existing spline/plane/sphere/cylinder fixtures run on the original
binary64 curve, exact conversion, rational refinement, and exact removal back
to the original. Independent Cox coefficients verify the refined homogeneous
function, and every contact count, enclosure, order and overlap must equal the
independent original certificate. The native bridges likewise check all four
representations against independently recomputed certificates for the same
214 plane and 92 quadric inputs. Existing reviewed OCCT differences remain
separate from native matches; these are not 1,224 distinct native geometries.

The 274 new Python `Fraction` fixtures manufacture known-factor polynomials,
verify the complete implicit equation in Bernstein/power form, and provide
exact rational contacts and minimal binary64 bounds (or explicit conversion
failure). They cover degrees through 25, contact orders through 50, rational
trims, roots separated by 2^-2048, mixed representable/unrepresentable results,
huge/subnormal parameter domains, contained arcs and periodic crossings. Ordinary
tests also check exact seam-origin removal, malformed rational input, atomic
work-limit failure and huge control values cancelling to a finite contact.

See [the exact result contract](EXACT_SPLINE_INTERSECTIONS.md). Enclosure
failure does not discard the exact contacts. The fourteenth fuzz target checks
complete known-factor, secant and periodic answers after exact edits, including
independent full polynomial identity and output-bound checks.

## Exact curve knot refinement and removal

The native bridge retains all operation flags and complete curve representations
for 667 inputs, including 665 captures made before Rust implementation and two
supplemental nonconstant seam cases. The 723 ordinary fixtures independently
solve every homogeneous Cox coefficient equation; selected cases also use a
separate Greville collocation reconstruction. Rust's coefficient solver checks
the entire unclamped raw support and full periodic function, so inactive end
controls and seam correspondence cannot disappear behind point samples.

Debug/release tests cover rational cuts separated by `2^-2048`, exact derivative
overflow, discontinuities, positive-weight rejection, output preflight, periodic
origin changes and a retained degree-25 solver timeout. Linux release CI runs
the second exact equation oracle over every fixture, including failed removal.
The thirteenth sanitizer target mutates edit sequences and independently solves
removal feasibility. [Native differences](NATIVE_KNOT_EDITING_DIVERGENCES.md)
remain separate from parity matches. No new unchanged upstream DRAW passes are
claimed; surface knot editing, general degree elevation and topology remain out
of this milestone's scope.

```sh
python3 rust/tools/generate_knot_editing_fixtures.py --check
python3 rust/tools/compare_knot_editing.py --occt-root /path/to/occt
```

## Exact tensor patches and isocurves

The 691-case native bridge captures independently before Rust execution and
compares all 4,023 patches and 324 isocurves. The complete rational homogeneous
outputs must agree with independent Cox tensor polynomials and binomial affine
substitutions. The macOS OCCT 7.9.3 capture has 637 matches and 54
[reviewed degree-25 Segment differences](NATIVE_SURFACE_EDITING_DIVERGENCES.md).
Linux OCCT 7.6.3 has 635 matches and 56 separately pinned differences. All
complete exact Rust outputs are byte-identical on the two platforms.
All counts, degrees and original U/V domains agree. Native differences cannot
exempt Rust from the exact oracle.

The 745 ordinary cases retain every complete homogeneous control using a
lossless shared-denominator text format. They include all native input families,
subnormal/adjacent knots, extreme homogeneous values, multiple turns, parameter
rectangles whose interior knots round to the same float, and seeded rational
surfaces. Full degree-25 x 25 native operations are recomputed in the bridge;
Cargo fixtures retain extraction and composed edits for all four periodicity
combinations. Exact boundary curves, mixed derivatives, axis commutation,
rational cuts 2^-512 apart, invalid raw rational parameters and Cartesian
budget exhaustion have additional tests. The twelfth fuzz target checks complete
tensor coefficient identities, exact jets, enclosures and operation sequences.

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

The Bézier second oracle also has an exhaustive gate:

```sh
RUSTY_VERIFY_ALL_BEZIER_ORACLES=1 cargo test --locked --release \
  --test bezier_editing exact_controls_match_independent_polynomials_and_original_parameter_jets
```

Linux release CI verifies all 636 Bézier and 745 tensor fixtures against both
the independent Python expectations and the separate Rust coefficient/jet
checker. Shared integer denominators reduce checker replay cost without
changing its equations or assertions. Retained fuzz inputs and failed startup
campaigns are documented in [the regression record](fuzz/regressions/README.md).

## Exact surface knot representations

`generate_surface_knot_fixtures.py --check` independently regenerates 1,756
complete U/V knot edits, success flags and homogeneous grids using Cox power
coefficient equations over full raw support. Ordinary Rust tests consume all
fixtures. The original 791 surface jet fixtures additionally verify the retained
exact representation and reusable differentiated de Boor evaluation.
The Linux job also sets `RUSTY_VERIFY_ALL_SURFACE_KNOT_ORACLES=1` to run every
surface-knot fixture through the separate Rust coefficient/feasibility solver.
`compare_surface_knots.py` compares 1,748 complete native observations;
[native round-trip differences](NATIVE_SURFACE_KNOT_DIVERGENCES.md) remain
separate from matches. Six deliberate-corruption bridge tests protect control,
weight, flag, input-set and versioned-review comparisons.

The `surface_knots` sanitizer target checks every transverse coefficient,
independent exact removal feasibility, whole patch polynomials, quotient jets,
isocurves, periodic wrapping, malformed data and operation sequences. Its
retained degree-25 timeout reproducer is also an ordinary Cargo regression.
The fifteen-target daily workflow retains evolving corpora and artifacts.

## General degree-elevation acceptance

All 142 curve and 125 tensor fixture grids are reconstructed with Python
Fraction/Cox equations and independently with the Rust integer equation solver.
`RUSTY_VERIFY_ALL_DEGREE_ORACLES=1 cargo test --release --test degree_elevation`
checks every second-oracle result, including unclamped inactive controls.

The source-pinned OCCT SDK is built from commit `3d097a0328e71b826377d4814ab05ec3c3d23871`
with exception checks enabled. The native comparison runner reports 187 matches,
46 reviewed native exceptions and 34 reviewed invalid grids. Fingerprints cover
every input, exact output, native output/diagnostic and source revision; crashes
remain failures. Deliberate-corruption tests reject altered controls, flags,
domains, case sets, stale captures and wrong runtime libraries. See
[degree elevation](DEGREE_ELEVATION.md) for the native comparison commands.

Revision `27647fadbb82ca8e516b66b606687d65144e0cf7` passed all sixteen Linux
sanitizer targets, each with a full 60-second mutation budget after replay,
and a separate 600-second local degree-elevation mutation campaign. Its
platform gate passed 120 tests in debug and release on Linux, macOS and Windows,
Rust 1.85, WebAssembly compilation and all five Linux second-oracle checks.
The retained periodic tensor timeout is covered by ordinary tests as well as
fuzzing. [The acceptance evidence](DEGREE_ELEVATION.md#accepted-revision-and-campaign-evidence)
pins the code revision, workflow runs, observed limits and remaining scope.

## Exact ordering across algebraic equations

`generate_root_comparison_fixtures.py --check` verifies 130 equation pairs using
irreducible-factor identity and continued-fraction isolation. Cargo checks all
pairwise root comparisons in both directions, including exact ties, degree-25
clusters, unrepresentable roots and shifts below binary64 resolution. The
existing 95 real-root fixtures also check every within-equation ordering.
The `roots` sanitizer harness constructs independent one- and two-root factor
equations to check cross-equation equality and ordering throughout campaigns.
This foundation does not yet implement rational-function image comparison or
complete curved-distance minima.
