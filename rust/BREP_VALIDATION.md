# Generic B-rep validation

`Topology::from_parts` builds a cell complex (`TOPOLOGY_MODEL.md`) from
caller-supplied vertices, edges, fins, loops, faces, shells and regions only
when the complete validation contract holds. Otherwise it returns **every**
issue found, each with its offending entity. `TopologyParts::check` and
`Topology::check` return the same list without constructing anything, and
`Topology::validate` reports the first issue as an `InvalidTopology` error.
The prism builders use the same validator.

Supported geometry is what the kernel's topology can represent: line segments,
circular arcs and full circles (`Curve3::Circle`, the curve of a ring edge) in
3D, line and arc pcurves, and plane and cylinder surfaces. Space is
partitioned into regions: region 0 is the infinite void, and every bounded
region, solid or void, lists its shells, the first outer and the others
cavities. Each face has a front and a back side; its oriented normal points
from the front side's region to the back side's. Each side is listed by one
shell. Cylinders carry no seams: a wall that closes around the axis is bounded
by ring loops with winding numbers.

## Result contract

An `Issue` is an `IssueKind` and an `Entity`. It displays as `kind:entity`,
e.g. `pcurve_off_edge:use 2.0.1` for face 2, loop 0, fin 1. Loops and fins
are indexed within their face; a fin or loop that no face reaches is named by
its arena index (`fin 7`, `loop slot 3`). Shells and regions are named by
index (`shell 1`, `region 2`). The list is sorted and free of duplicates.
Issue classes:

| Class | Kinds |
| --- | --- |
| References and usage | `reference`, `unused_vertex`, `unused_edge`, `face_without_shell`, `face_reused`, `empty_face`, `empty_loop`, `empty_shell`, `fin_without_loop`, `fin_reused`, `loop_without_face`, `loop_reused`, `edge_fins_mismatch` |
| Sides and regions | `side_without_shell`, `side_in_two_shells`, `side_region_mismatch`, `region_without_shell`, `no_infinite_region`, `region_shell_mismatch`, `double_bounding` |
| Loops and shells | `open_loop`, `winding_mismatch`, `ring_edge_with_vertex`, `ring_edge_open`, `seam_edge`, `free_edge`, `non_manifold_edge`, `same_sense_uses`, `radial_order_inconsistent`, `edge_across_shells`, `disconnected_shell`, `non_manifold_vertex`, `euler` |
| Definitions | `degenerate_vertex`, `degenerate_curve`, `degenerate_surface`, `degenerate_pcurve` |
| Certified geometry | `vertex_off_curve`, `vertex_loop_off_surface`, `pcurve_off_edge`, `uv_gap`, `loop_winding`, `inner_loop_outside`, `shell_orientation`, `cavity_outside`, `nested_cavity` |
| Undecided geometry | `uncertified_vertex_off_curve`, `uncertified_vertex_loop`, `uncertified_pcurve_off_edge`, `uncertified_uv_gap`, `uncertified_loop_winding`, `uncertified_containment`, `uncertified_shell_orientation` |

An out-of-range reference stops validation after the reference pass, since
nothing else can be indexed safely. Other failures gate only the checks that
depend on them. A degenerate curve skips its vertex and pcurve checks. A face
with a structural, side or geometric issue skips winding. A shell with an
issue skips orientation and containment. A shell listing exactly the opposite
sides of an earlier shell (the void's view of the same surface) is that
shell's twin, and its shell-level checks are not repeated. This gating is part of the contract, so the
independent reference validator reproduces complete issue lists, not just
verdicts.

Conventions: every curve and pcurve uses a normalized fraction `t` in `[0,1]`.
A fin's pcurve follows the oriented face, so a reversed fin pairs pcurve
fraction `t` with edge fraction `1-t`. Outer loops wind counter-clockwise
about the oriented face normal and inner loops clockwise. Cylinder UV is
(angle, axial height), and pcurves live on its universal cover: a loop with
winding number `w` closes with its end `2πw` in `u` after its start.

## Exact combinatorial checks

* Every vertex, edge and face is used, every fin is in exactly one loop and
  every loop in exactly one face, and each edge lists exactly the fins that
  name it. No loop or shell is empty; a face with no loops is `empty_face`
  (a closed surface without boundary is not yet representable).
* Each face side is listed by exactly one shell, the one the face records for
  that side. Each shell is listed by its region, no region lists a shell twice,
  every bounded region has a shell, and region 0 is a void.
* Each loop of edges closes in 3D vertex order; a ring edge alone closes a
  loop. A ring edge has no vertices and a closed curve. Winding numbers are
  zero on planes and in `v`; the windings of a cylinder face sum to zero.
* Around each edge, the regions ahead of and behind its fins, in the stored
  radial order, must alternate. A single fin is `free_edge`, two same-sense
  fins are `same_sense_uses`, two opposite fins out of order are
  `radial_order_inconsistent`, more are `non_manifold_edge`, and fins between
  different shell pairs are `edge_across_shells`. Two fins of one edge in one
  face are a seam, which this model forbids (`seam_edge`).
* Each shell is face-connected. At each vertex, the graph of edges joined by
  consecutive fins around that vertex must be connected, so pinch vertices are
  rejected. Each shell's Euler characteristic `V - E + 2F - L`, with ring
  edges excluded and every loop counted (a vertex loop adds its vertex and
  its loop), must be even and at most 2.

## Certified geometric checks

A geometric check passes only on a certified upper bound `<= tol`. It fails
only on a certified lower bound `> tol`. Otherwise it reports the matching
`uncertified_*` issue. No decision rests on unverified floating point.

* **Definitions.** Finite values, nonzero line length (`|b-a| > tol`), radius
  `> tol`, and `0 < |sweep| <= 2π`.
* **Vertices.** Each edge end is within tolerance of its vertex, and each
  vertex loop's vertex within tolerance of its face's surface.
* **Curve on surface.** Over a whole use, `D(t) = C(t_edge) - S(P(t))` is a
  harmonic sum: `A0 + A1 t + Σ_ω (C_ω cos ωt + S_ω sin ωt)`, with exact
  rational frequencies. This covers line and arc edges against plane
  line/arc pcurves and cylinder line pcurves. Equal frequencies are combined
  first, so a correctly matched arc cancels exactly. Then
  `sup |D| <= sqrt(|A0|² + |A0+A1|²) + Σ sqrt(λmax(Gram(C_ω, S_ω)))`, since the
  affine part is convex and each harmonic term is an ellipse. Terms at
  nearby but different frequencies `ω1 < ω2` (a circle swept by `2π` against
  a pcurve spanning OCCT's printed `6.28318530717959`) would each count at
  full size although they almost cancel. With `z = C - iS`, their sum is
  `Re(z1 e^{iω1 t} + z2 e^{iω2 t})`, and since `|e^{i(ω2-ω1)t} - 1| <=
  (ω2-ω1) t`, on `[0, 1]` it is at most the ellipse bound of `z1 + z2` plus
  `|z2| (ω2 - ω1)`. Terms adjacent by frequency are bounded that way whenever
  it is smaller. A failure is certified when one of 33 sample points lies
  beyond tolerance. An arc pcurve
  on a cylinder is not harmonic and can only be certified as a failure.
* **UV closure.** Consecutive fins meet in UV within tolerance, with cylinder
  angle differences scaled by the radius; the last fin meets the first shifted
  by `2πw` in `u`.
* **Loop winding.** Loops are closed exactly by straight chords between
  consecutive fins (gaps certified to be within tolerance). Twice the signed
  area is the closed form of `∮ u dv - v du`, including the chords. Its sign
  must be positive for an outer loop and negative for an inner loop, relative
  to the face orientation. On a wound cylinder face the loops have no
  outer/inner order: the total periodic area `-∮ v du` over all loops must
  have the face's sign, and each unwound loop the opposite one.
* **Inner loops.** The first point of each inner loop must lie inside the outer
  loop. A `+u` ray uses a half-open crossing rule. Arcs are split at their
  `v` extrema `π/2 + kπ`, using a certified `π`, so each piece is monotone.
  An unwound loop on a wound face is `uncertified_containment` for now.
* **Region orientation.** The flux of `x/3` through a face is
  `-∮ v f(u) du` over its loops, closed chords included, with
  `f = S·(S_u × S_v)` integrated in `v` from 0: for a plane `f = o·(x×y)`, for
  a cylinder `f = -r(o·(x×n)) sin u + r(o·(y×n)) cos u + r²·det(x,y,n)`. A
  shell's flux sums its faces' fluxes, negated for back sides. A bounded
  region's first shell must have positive flux and every other shell
  negative; the infinite void's shells negative.
* **Cavities.** A cavity's point (its first fin's start, or the surface point
  there for a ring fin) must lie inside its region's outer shell and outside
  every other cavity. Rays from that exact point in up to eight fixed integer
  directions are intersected with each face by exact Cramer solves (planes)
  or an exact quadratic (cylinders). A direction counts only when every hit is
  certifiably more than twice the tolerance from its face's boundary, and a
  start point on a face plane is certified outside that face. Adjacent faces
  only meet within tolerance, so without this margin, a ray through a shared
  edge could fall between the two faces' closed UV regions. With it, the
  parity is the same for any watertight surface within tolerance of the
  faces. On a cylinder, a hit is inside the face when the `+v` ray from it
  on the universal cover crosses the face's loops an odd number of times,
  counting every period alias of a wound loop.

## Two arithmetic tiers

All geometry is generic over `certified::Real`. Each decision first runs in
`Fast`, a binary64 interval, and falls back to `Interval`, a rational interval
on a `2^-192` outward grid, only when the first tier cannot decide. Both are
enclosures; neither decides from an unverified float.

* `Fast` rounds outward only when an operation was inexact. The error sign of
  each sum comes from TwoSum and of each product from a fused multiply-add
  (outside the underflow range). Exact results stay points, so exact vertex
  coordinates, exactly shared pcurve endpoints and the angle zero keep exact
  signs. Rationals convert to the tightest binary64 bracket: a leading-bits
  quotient, corrected by exact comparisons.
* `Fast` cosine and sine reduce by the binary64 bracket
  `[FRAC_PI_2, next_up(FRAC_PI_2)]` of `π/2` (checked against the Machin
  enclosure) and sum 13 Taylor terms with the alternating-series remainder.
* `Interval` uses a Machin `π`, argument-reduced Taylor series for cosine,
  sine and `atan`, and an integer square root at `2^-200`, all rounded
  outward to the grid. Rational trigonometry is
  memoized per exact angle.
* A line-length decision below the rational grid compares exact rationals.

On the development machine, the 256-prism invariant suite took 101.9 seconds
with rational intervals only. With the fast tier it takes 0.19 seconds, and no
decision falls back to rationals. Before this validator, with sampled checks,
it took about 0.05 seconds.

## Independent evidence

* `brep_reference.py` is a separate Python validator in mpmath. It
  reimplements the harmonic deviation bound from its own curve and surface
  formulas, certifies failures on 257 samples instead of 33, integrates loop
  areas by Gauss–Legendre quadrature, decides inner loops by a winding-angle
  integral and uses its own ray parity for cavities. Its sign decisions keep
  explicit margins far above mpmath's working precision. `cell_reference.py`
  extends it to the cell model: `to_cell` converts a seamed model by rule
  (a seam pair on a cylinder merges only when its two fins run opposite ways,
  lie exactly one period apart and continue their neighbours in UV; any other
  seam is kept and rejected as `seam_edge`), and `validate` implements the
  side, region, radial, winding, vertex-loop and region-flux invariants
  independently, with its own `+v` cover-crossing parity on cylinders.
  `generate_brep_fixtures.py --check` rebuilds 66 cases. 54 come from an
  independent seamed prism builder (19 valid solids and 35 mutations; the
  valid solids include holes, convex and concave arcs, full circles, one and
  two cavities, rotated and far-translated copies and a millimetre-scale
  box), converted to cells. Twelve are cell-model cases with no seamed form: a
  cylinder parametrized from `π`, a box with a valid vertex loop, a ring
  pcurve spanning OCCT's printed period `6.28318530717959` (valid only
  through the near-frequency bound; without it the reference cannot decide),
  the same pcurve off by `10^-6` of a period, and the model's own failure modes
  (a fin shifted by a period, a flipped winding, fins swapped between edges, a
  side at the wrong shell, a vertex loop off its surface, a face without
  loops, a ring edge with one vertex and a shell its region does not list).
  24 cases are valid. `brep_validation.rs` requires
  Rust's complete sorted issue list to equal the reference's for every case.
* The existing prism suites (`invariants`, `occt_regression`, `modeling`)
  build every solid through the new validator.

The nineteenth fuzz target, `brep_validation`, builds valid star-outline
prisms with no hole, a round (seamless) hole, a square hole or an inverted box
cavity. Scales range over `2^±10` with random frames and offsets. The cavity
is merged by fuzz-crate code, not by a kernel builder: its material shell
joins the body's solid region and its twin bounds a new void region. The
base must be valid. Then one of 24 mutations must produce its predicted
issues:

* an exact report for local changes: an extra vertex, edge, empty shell or
  empty loop, a bad edge reference, an inverted outer shell or an uninverted
  cavity
* a required issue on the mutated entity for the others: a face whose sides
  leave their shells or repeated in one, a flipped fin, a moved vertex, a
  shifted pcurve, a flipped face or swapped outer and inner loops
* the cell model's failure modes (`TOPOLOGY_MODEL.md`): a wall fin shifted by
  one to three periods (on the stadium fixtures, the only prisms with
  multi-fin loops on a cylinder: `uv_gap` at both of its ends), a flipped ring
  loop winding (`uv_gap` and `winding_mismatch`), two edges' fin lists
  exchanged (`edge_fins_mismatch` on both), a front side recorded at the back
  shell (`side_region_mismatch`), a vertex loop off its surface, and an open
  face's only loop removed (`empty_face` and `loop_without_face`)
* laws that must stay valid: a vertex loop on a cap, and a two-fin edge's
  radial order reversed (a cyclic order of two is unchanged)
* for a pcurve shifted by `10·tol`, a clean report at `100·tol`

Every report must be deterministic, duplicate-free and identical to
`from_parts`. A local 300-second development campaign (dirty tree based on
`12098ee3`, AddressSanitizer, standard 20-second/2 GiB limits) ran 13,807
mutation executions with 4,408 coverage edges and a 762 MB RSS peak. It
produced no artifacts. The clean-revision campaign is recorded under
[Acceptance](#acceptance).

## Native OCCT observations

Before any Rust implementation existed, the source review read
`BRepCheck_Analyzer` and `BRepCheck_{Vertex,Edge,Wire,Face,Shell,Solid}` at
`3d097a0328e71b826377d4814ab05ec3c3d23871` (see `SOURCE_MAP.md`). A
source-pinned headless SDK (OCCT 8.1.0) then analyzed the 52 cases that OCCT
can represent; the two out-of-range references cannot be built. That capture
is checked in under `fixtures/occt-brep-preimplementation/`, with its notes,
inputs and observations. Every case completed within 0.02 seconds.

`occt_brep_check_oracle.cpp` builds each case with `BRep_Builder` from explicit
rows written by `brep_reference.native`. It computes no geometry, tolerance or
orientation itself. Rust fractions become OCCT parameters with matching speed:

* lines use arc length
* circles use the angle, with the axis flipped for negative sweeps
* plane pcurves use the same arc length or angle
* cylinder pcurves are lines whose parameter is the angle or axial length

Seams use `UpdateEdge(E, C_forward, C_reversed, F)`. Where a mutated pcurve
cannot match its edge's speed, the use is marked inexact, and that case is
never read as a native verdict. The native rows stay seamed: they are built
from the seamed models, and the Rust side validates their converted cells. The
oracle's second row per case counts the solid's distinct vertices, edges,
wires, faces, shells and solids, as DRAW's `nbshapes` does.

## Native comparison bridge

`compare_brep.py` first verifies the SDK manifest, loaded toolkits and the
unchanged pre-implementation capture. The generator computes trigonometry and
norms with mpmath and rounds each result once to binary64, so every host and
Python version produces the same bytes. Platform libm results differ in the
last bit, and CPython 3.10 changed `math.hypot`: the first Linux CI runs
regenerated different inputs from the macOS capture. The capture itself used
macOS libm and Python 3.9, so the regenerated native inputs must match it
token for token, with numbers within `2^-50`. Only components of the triangle,
hexagon and two-hole plate moved, by at most `2.3e-16`, and all 54 expected
reports are unchanged. The bridge also checks byte-identical regeneration of the
fixtures, and that the Rust probe's issue lists equal
the independent reference's. A case matches when the verdicts agree and every
Rust issue class with a BRepCheck counterpart has one of its corresponding
native statuses. For example, `pcurve_off_edge` corresponds to
`InvalidCurveOnSurface` and `shell_orientation` to `EnclosedRegion` or
`SubshapeNotInShape` (see `CORRESPONDING` in `brep_reference.py`). Each native
process has a 120-second deadline. Timeouts, crashes and malformed output can
never be reviewed.

Two rules bridge the seamed and seamless encodings (`TOPOLOGY_MODEL.md`).
Native seam edges (used twice by one loop) and the vertices only seams and
edges closed on them use are **structure-only**
(`brep_reference.structure_only`): their statuses cannot stand for a Rust
issue class. Every case valid on both sides is then checked by **count
synthesis**: the kernel's `Topology::occt_counts` of the Rust cell (one seam
per wound face, one seam vertex per ring edge, one wire per wound face plus
one per unwound loop) must equal the native count row exactly. Reviews
fingerprint the status row, which the count row follows, so the eight
existing reviews still bind the same native observations.

On macOS and Linux, 44 cases match and 8 are reviewed differences, with no
failures. After the cell-model migration (T1) the same counts hold. The only
structure-only status is the seam's `InvalidCurveOnSurface` on the zero-radius
cylinder, whose Rust issue (`degenerate_surface`) compares by verdict; all 21
cases valid on both sides agree on the synthesized counts.
`occt-brep-divergences.json` fingerprints each observation with its reason
and independent evidence:

* **Unused vertex and edge.** `BRepCheck_Analyzer` only visits subshapes of
  the solid, so these are absent and OCCT reports the solid valid. The Rust
  contract validates every supplied entity.
* **1e-3 and 5e-4 UV gaps.** `BRepCheck_Wire::Closed2d` compares only the wire's
  first/last junction. `IsDistanceIn2DTolerance` also accepts any gap below 1%
  of the face's UV extent before consulting tolerances. OCCT still rejects
  both cases through `InvalidCurveOnSurface`.
* **A face listed twice.** OCCT reports `NotConnected` and
  `InvalidImbricationOfShells`. The Rust contract also reports the edges of
  the excluded duplicate face as free.
* **Three inexact encodings.** A zero-length edge, a changed arc sweep and a
  cavity face using an outer edge. Both verdicts are invalid.

`test_brep_oracle.py` checks malformed native rows, the difference
classification, the structure-only rule (a seam status cannot match an issue;
a cylinder's seam and seam vertices are structure-only, a box has none) and
the review fingerprint rules. Linux CI builds the same
pinned SDK (OCCT 8.1.0). Its native output was byte-identical to macOS for
every case, so the same eight fingerprinted reviews apply and no Linux-only
record was needed. The slowest Linux native case took 0.007 seconds.

## Acceptance

Revision `dff912e5` passed every gate:

* **Rust kernel workflow:** all ten jobs passed. They cover:
  * Linux, macOS and Windows debug/release tests
  * Rust 1.85 and WebAssembly
  * the original-test bridge, including `generate_brep_fixtures.py --check`
  * four source-pinned native comparisons

  The B-rep comparison certified all 54 Rust reports: 44 matches, eight
  reviewed differences and no failures.
* **Geometry fuzzing workflow:** all nineteen targets passed. The Linux
  `brep_validation` campaign replayed 775 corpus inputs in 106 seconds, then
  completed 60.08 seconds of mutation (1,020 executions, 4,499 coverage
  edges) with a 640 MB RSS peak and no artifacts.
* **Clean local 600-second campaign** at `24244b2e`, whose kernel and fuzz code
  are unchanged at `dff912e5` (only the Python generator's rounding,
  regenerated fixture text and docs changed). It completed 600.09 seconds of
  mutation after 55.78 seconds of replay, with 23,390 mutation executions,
  4,533 coverage edges and a 789 MB RSS peak. There were no crash, timeout,
  OOM, slow-unit or disagreement artifacts.

The two commits after `24244b2e` fix host-dependent fixture generation. The
first Linux runs regenerated different bytes because of platform libm
trigonometry and CPython's changed `math.hypot`; see
[the bridge](#native-comparison-bridge). They are observations of this
revision and these runners, not portable latency guarantees.

The cell-complex migration (T1) was accepted at `e4adb869` with the same
bridge counts, the synthesized-count checks above, 66 fixture reports and a
clean 600-second `brep_validation` campaign; the record is in
`TOPOLOGY_MODEL.md`.

## Limits

2D loop self-intersection, face/face intersection (which OCCT's BRepCheck does
not check either), tolerance healing, per-entity tolerances and operation
history are out of scope. Validation cost is not bounded by explicit work
limits. The curve checks are linear in the number of fins, but the
combinatorial passes use ordered maps, and containment is linear in faces per
ray. This validator certifies the supplied boundary; it does not make an
invalid import valid. Free and non-manifold edges, wire edges and acorn
vertices are representable but rejected: every current operation requires a
solid. Faces without loops (closed surfaces), unwound loops on wound faces
(certified only as `uncertified_containment`), windings in `v` and
per-entity enclosures wait for the surfaces and operations that need them.
