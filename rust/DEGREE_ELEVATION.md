# Exact B-spline degree elevation

The general curve and tensor-surface implementation passed its mathematical,
source-pinned native, sustained fuzz and platform gates at code revision
`27647fadbb82ca8e516b66b606687d65144e0cf7`. The evidence below covers this
capability; full kernel parity and production readiness remain unfinished.

`ExactBSplineCurve3::elevated(degree)` and
`ExactBSplineSurface3::elevated(u_degree, v_degree)` preserve the four
homogeneous functions on the original active parameter domain. Equal degree
returns the same data. Degree reduction, degree above 25, and output counts
above 4096 are errors. Surface requests preflight both axes and the complete
Cartesian control count before evaluating controls. All results are immutable.

Knots and controls remain rational. The supported domain includes nonuniform
positive weights, periodic directions, unclamped axes, and rational atoms
beyond binary64's exponent range. Evaluation and derivative units retain the
original parameters. Positive weights are preserved by positive combinations
of the original homogeneous controls.

## Unclamped axes and periodic origins

Let the original degree be p, target degree q, and delta = q-p. Let f and l
be the zero-based indices of the first and last active distinct knots, and K
the original distinct-knot count. The output pole count is
`old_count + delta * (l-f)`.

For a nonperiodic axis, add delta to every multiplicity, remove exactly
`delta*f` flat knots at the beginning, and remove `delta*(K-1-l)` at the end.
Removal can reduce an exterior multiplicity without discarding its distinct
knot. This preserves the active interval: the count of knots strictly before
its start is unchanged while its endpoint multiplicity rises by delta.

For example, degree 2 with simple knots 0 through 7 has five controls and
active interval [2,5]. Elevation to degree 3 produces eight controls, knots
1 through 6 with multiplicity 2, and the same active interval [2,5].

Inactive controls follow a fully clamped working extension with temporary
zero homogeneous padding. The elevated working controls are cropped by the
original padding plus the flat-knot removals above. Removed basis functions
vanish on the original active interval. The complete resulting control array
is independently checked; the contract does not assert equality outside the
original active interval, where the available exterior support changes.

For a periodic axis, keep every distinct knot and increase each multiplicity
by delta. The pole count rises by `delta*(K-1)`. The cyclic origin and period
are unchanged. Its canonical extension length is also unchanged because
`q+1-(old_seam_multiplicity+delta) = p+1-old_seam_multiplicity`.

## Implementation and independent mathematics

The production map references `BSplCLib::IncreaseDegree` and its Prautzsch
averaging: duplicate controls by congruence class, increase the selected knot
multiplicities, refine the remaining knots, and average p+1 rank curves.
Sparse maps are reused for every transverse homogeneous field of a surface.
Identical source and destination axes also share the map between U and V.
Periodic working data spans three periods: the fundamental interval and one
neighbor on either side. With n cyclic poles and degree p, validation requires
n > p. Thus the unclamped working active interval contains the whole central
period; zero padding at the two outer ends does not change its function.
After elevation, n' > q and the canonical extension length e remains unchanged.
In particular, n'-e is at least the new seam multiplicity, so the first and
last canonical extension knots lie strictly inside the outer period endpoints.
The complete canonical extended knot sequence selects the cyclic controls,
preserving their original origin. Two more outer periods would add computation
without contributing support to this result.

The Python oracle in `tools/degree_elevation_reference.py` instead solves
every Cox power-coefficient equation on every support interval using exact
fractions, with full-rank and residual checks. Its nonperiodic working
extension and cropping are explicit. Its periodic system covers a complete
fundamental interval. It does not use production insertion or averaging.
The tensor oracle also verifies exact equality of both axis orders.
The separately implemented Rust oracle uses primitive integer elimination
on the full coefficient system, including temporary zero homogeneous fields.
It reconstructs every working control before the independently defined crop.
Zero pivot columns are skipped; nonzero eliminations remove common integer
content with checked exact divisions. Full rank and every residual equation
remain required. This avoids inflating sparse degree-changing systems with
unrelated minors across hundreds of transverse fields.
For wide tensor grids the solver reconstructs every source unit column once,
proves its complete coefficient identity, and applies that independent map to
every homogeneous field. Coordinate magnitudes therefore do not inflate the
matrix factorization. No production transform is used by either oracle.

The initial corpus has 142 curve and 125 surface cases: no-ops, single- and
both-axis elevation, degrees through 25, clamped and unclamped support,
independent periodicity, and unit/nonuniform weights. It includes the original
curve `IncreaseDegree` and `RationalCurveIncreaseDegree` GTests and the original
surface `IncreaseDegree` GTest. Rust tests compare every resulting knot,
multiplicity, domain, and homogeneous control. Additional tests cover extreme
rational values, staged elevation, transposition, and atomic count rejection.
The tightest periodic support (only p+1 cyclic poles) is checked separately
for every raising source degree, low/middle/full seam multiplicity, a single
increment and a jump to degree 25. These tests pack every source unit column
into tensor coordinate fields and reconstruct the complete map with the
independent global periodic coefficient solver. The saved degree-24 periodic
tensor timeout also reconstructs all 54-by-54 resulting controls and checks
transposition in the ordinary test suite.

Regenerate or verify the independent fixtures with:

```sh
python3 rust/tools/generate_degree_elevation_fixtures.py --check
cargo +stable test --locked --release --test degree_elevation
RUSTY_VERIFY_ALL_DEGREE_ORACLES=1 cargo +stable test --locked --release --test degree_elevation
```

## Native evidence before implementation

Reference source: OCCT commit
`3d097a0328e71b826377d4814ab05ec3c3d23871`. The native geometry libraries were
built from a pristine archive of that commit, with release exception checks
enabled and all default modules disabled, requesting only TKG3d and its
dependencies. The runtime reports OCCT 8.1.0. Loaded-library paths and library
hashes were checked locally. Git suffix generation was disabled because an
archive inside the Rust checkout would otherwise inherit the enclosing
repository's revision in its version string.

`tools/build_pinned_occt.py` reproduces the source archive and headless SDK in
a fresh output directory. It records the exact build commands, source archive,
CMake configuration, compiler-command database, generated version header, and
installed library hashes. Existing output is preserved. CMake and a C++17
compiler must already be available; the build does not install system packages.

```sh
python3 rust/tools/build_pinned_occt.py --output target/pinned-degree-sdk --jobs 2
```

Use `--cmake /path/to/cmake` when it is outside the executable search path.
The initial builder was checked on macOS with CMake 3.31.10 and Apple Clang;
its 267 native observations reproduced the preimplementation capture exactly.

The capture occurred at Rust commit
`54ce761f39b09b1d8697e00f4a39cfc5fb427364`, before a Rust degree-elevation
implementation existed. `fixtures/occt-degree-elevation-capture.json` retains
the input/output hashes and every per-case outcome.

| Corpus | Complete matches | Native exceptions | Invalid native axes/grids |
| --- | ---: | ---: | ---: |
| Curves, 142 cases | 106 | 24 | 12 |
| Surfaces, 125 cases | 81 | 22 | 22 |

These exceptions and invalid outputs also occur in the installed OCCT 7.9.3
probe. They remain explicit discrepancies, not passing matches. In the unclamped example above, native OCCT reports
eight controls but returns knots and multiplicities requiring ten controls
and a changed domain. The Rust target layout is derived from the coefficient
contract and the original interval, rather than reproducing that malformed
result.

`tools/compare_degree_elevation.py` reconstructs all expected controls from
fresh Python equations, checks Rust exactly, then compares the complete native
result. It validates axis multiplicities, pole counts and original domains;
native Cartesian poles and weights are compared as homogeneous coefficients.
The comparison budget is `1e-10 + 2e-12*abs(exact)` per homogeneous field.

`fixtures/occt-degree-elevation-divergences.json` reviews the 46 native range
exceptions and 34 invalid axes/grids. Every record pins the input, exact result,
native stdout/stderr, native source revision, classification and differences.
The review records include the independently preserved domain and, for
malformed outputs, the declared versus required control counts. A new result
needs a new review; a native crash cannot receive a behavior exemption.
Linux x86_64 with GCC 13.3.0 reproduces the same match/exception/invalid-grid
counts. Twenty-two malformed outputs have different serialized control fields
from the macOS capture, with identical knots, domains and declared/required
counts. Their additional fingerprints and component-change summaries are
reviewed separately. They were captured after implementation; they do not add
new test cases or matches to the preimplementation corpus.

```sh
python3 rust/tools/compare_degree_elevation.py \
  --occt-root target/pinned-degree-sdk/install \
  --sdk-manifest target/pinned-degree-sdk/build-manifest.json
```

`--strict-native` disables all reviews and currently fails on those 80 explicit
differences. `--reuse-capture` accepts only unchanged corpus/probe/source/SDK
and observation hashes. The runner also verifies every resolved OCCT runtime
library against the selected SDK. CI builds the pinned headless SDK separately
and retains observations, reports, library hashes and build logs.

## Accepted revision and campaign evidence

The [platform workflow](https://github.com/jackControls/rustyOCCT/actions/runs/35904050143)
passed at the revision above: 120 tests in each debug and release suite on
Linux, macOS and Windows, 120 tests on Rust 1.85, Clippy, formatting and a
WebAssembly library compile. Linux also ran all five complete second-oracle
checks, including the seven degree-elevation tests. The source-pinned native
job reproduced 187 matches and 80 reviewed discrepancies with no new failures.
The original DRAW bridge passed three unchanged tests on both kernels;
three Rust cases remain unsupported and one lacks its external data fixture.

All sixteen [Linux sanitizer targets](https://github.com/jackControls/rustyOCCT/actions/runs/35904050060)
completed at least 60 seconds of mutation after corpus replay. The degree
target performed 43 mutations; its slowest input took 55 seconds and its peak
resident memory was 400 MiB. A separate clean-revision local degree campaign
completed 600.03 seconds of mutation and 359 mutations after 397.79 seconds of
startup/replay, with a 21-second slowest input and 646 MiB peak resident memory.
It retained 278 corpus files. These timings describe observed campaigns, not
worst-case bounds for arbitrary rational inputs.

The two earlier Linux degree campaigns failed on the same periodic tensor
timeout during replay. Their saved input remains a regression. Sharing equal
U/V transforms and reducing the periodic working support from five to three
periods resolved that input without relaxing assertions or resource limits.
Later documentation commits do not change the revision to which this evidence
applies. General intersections, curved distances, B-rep operations and full
production guarantees remain separate work.

The `degree_elevation` fuzz target reconstructs every successful control grid
independently. It covers degrees 1..25, all five axis families in the fixture
corpus, raw binary64 atoms, rational exponents beyond binary64, staged elevation,
refinement composition, tensor axis order/transposition, exact jets, invalid
degrees and atomic 4096-control limits. Its instrumented resource gates remain
60 seconds per complete input and 2 GiB RSS; corpus replay does not count as
mutation time. No complete campaign is claimed until its manifest passes.
