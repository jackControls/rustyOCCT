# Native complete linear intersection observations

The source reference is `3d097a0328e71b826377d4814ab05ec3c3d23871`.
The initial 433 native input geometries were captured before writing the Rust
implementation. These are native API comparisons, not additional passing
unchanged DRAW tests. The initial plane pair in `tests/lowalgos/intss/buc60815`
is included; that test's later swept-surface operations remain unsupported.

The native harness preserves both `BRepAlgoAPI_Common` and
`BRepAlgoAPI_Section` outputs. It reports finite vertices, whole finite edges,
finite planar faces, and any returned infinite lines/planes. Faces are compared
as filled regions, not just samples or boundary points. Every Rust result must
first equal the complete independent Python `Fraction` boundary-intersection
answer, including all extreme points and canonical unbounded constructions.

## OCCT 7.9.3 on macOS and 7.6.3 on Linux

| Observation | Count |
| --- | ---: |
| Inputs / independently verified Rust sets | 433 |
| Native COMMON sets matching their dimension-filtered contract | 394 |
| Of those, intentional lower-dimensional contact omissions | 154 |
| Native SECTION sets contained in the exact intersection | 430 |
| Combined COMMON/SECTION sets matching the complete intersection | 384 |
| Separately reviewed cases | 49 |
| Native components retained across both operations | 597 |

Both runtime observations have these counts. Their 49 divergent cases have
identical input hashes, native coordinate bits and exact expected results;
the review registry has separate version pins (98 entries). All 433 complete
Rust results are byte-identical across the two platforms.
These counts overlap and must not be added as distinct case counts.
The positional budget is `1e-10 + 2e-12 * case_scale`, where scale comes from
input/exact geometry, never an erroneous native observation. Direction/normal
comparison uses `2e-12`. Classification and required region dimension are exact.

COMMON deliberately filters out contacts below the minimum argument dimension
(`BOPAlgo_BOP::BuildRC`). For example, a triangle/triangle point contact normally
has empty COMMON and a point SECTION. The bridge records these omissions and
requires the combined result to cover the complete contact. It does not call
an empty COMMON a match to the nonempty Rust closed-set intersection.

### Empty results for unbounded operands: 46 cases

Both native operations return empty for `FL`, `coincident_lines`,
`coincident_planes`, `contained_line`, and `FF`, each in nine scale/frame
combinations, plus the initial BUC60815 plane pair. Rust's exact intersection
is a line or plane. Thirty-six also differ from COMMON's dimension-filtered
expectation; ten plane/plane crossings correctly have empty COMMON but no
SECTION line. These observations establish a limitation of this native B-rep
comparison route, not that OCCT's lower-level analytical APIs cannot compute
these intersections. They are never counted as complete-set matches.

### Extra edge / off-line point in the oblique triangle-edge case: 3 cases

`triangle_edge_line` with frame 2 returns two COMMON edges and one SECTION
point, at scales 1/8, 1 and 32. At scale 1 the triangle vertices are `(2,-3,5)`,
`(10,1,5)`, `(2,5,9)`; the line runs through `(-2,-5,5)` and `(14,3,5)`.
The true intersection is the first triangle edge, from `(2,-3,5)` to `(10,1,5)`.
The additional COMMON edge ends at `(2,5,9)`, which cannot lie on a line whose
z coordinate is identically 5. SECTION returns approximately
`(2,1.3055050463303892,7.1527525231651943)`, also off that line. Both errors
are far outside the unchanged budget.

Exact affine constraints, the independent boundary algorithm, and defining-point
permutation checks agree. The immediate native algorithmic cause has not been
isolated; the observation is pinned without attributing an unverified cause.

## Review guards

[`fixtures/occt-linear-sets-divergences.json`](fixtures/occt-linear-sets-divergences.json)
pins runtime version, input SHA-256, every native component's coordinate bits,
complete exact result, and discrepancy categories. A changed observation or
unreviewed runtime fails. `--strict-native` rejects even reviewed differences.
No review can bypass independent Rust correctness checks. Deliberate-corruption
tests cover missing area, missing segment interiors, extra/wrong vertices,
malformed observations and all review pins.

The test harness reads vertex positions from `BRep_TVertex` and applies the
vertex location exactly as `BRep_Tool::Pnt` does. This avoids an unrelated
triangulation include whose transitive `NCollection_AliasedArray.hxx` header
is missing from the packaged Ubuntu 7.6 modeling development SDK. The kernel
and native oracle need no visualization SDK. Native build diagnostics are
retained in `native-build.log`, including failed compilations.
