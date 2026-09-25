# Generic B-rep validation

`Topology::from_parts` builds a topology from caller-supplied vertices, edges,
faces and shells only when the complete validation contract holds. Otherwise
it returns **every** issue found, each with its offending entity.
`TopologyParts::check` and `Topology::check` return the same list without
constructing anything, and `Topology::validate` reports the first issue as an
`InvalidTopology` error. The prism builders now use the same validator in
place of the earlier sampled curve/pcurve checks.

Supported geometry is what the kernel's topology can represent: line segments
and circular arcs (`Curve3::CircularArc`, with the full circle as a `TAU`
sweep) in 3D, line and arc pcurves, and plane and cylinder surfaces. Any
number of shells is accepted. The first bounds the solid and later ones are
cavities.

## Result contract

An `Issue` is an `IssueKind` and an `Entity`. It displays as `kind:entity`,
e.g. `pcurve_off_edge:use 2.0.1` for face 2, loop 0, use 1. Loops and uses are
indexed within their face. The list is sorted and free of duplicates. Issue
classes:

| Class | Kinds |
| --- | --- |
| References and usage | `reference`, `unused_vertex`, `unused_edge`, `face_without_shell`, `face_reused`, `empty_face`, `empty_loop`, `empty_shell` |
| Loops and shells | `open_loop`, `free_edge`, `non_manifold_edge`, `same_sense_uses`, `edge_across_shells`, `disconnected_shell`, `non_manifold_vertex`, `euler` |
| Definitions | `degenerate_vertex`, `degenerate_curve`, `degenerate_surface`, `degenerate_pcurve` |
| Certified geometry | `vertex_off_curve`, `pcurve_off_edge`, `uv_gap`, `loop_winding`, `inner_loop_outside`, `shell_orientation`, `cavity_outside`, `nested_cavity` |
| Undecided geometry | `uncertified_vertex_off_curve`, `uncertified_pcurve_off_edge`, `uncertified_uv_gap`, `uncertified_loop_winding`, `uncertified_containment`, `uncertified_shell_orientation` |

An out-of-range reference stops validation after the reference pass, since
nothing else can be indexed safely. Other failures gate only the checks that
depend on them. A degenerate curve skips its vertex and pcurve checks. A face
with a structural or geometric issue skips winding. A shell with an issue
skips orientation and containment. This gating is part of the contract, so the
independent reference validator reproduces complete issue lists, not just
verdicts.

Conventions: every curve and pcurve uses a normalized fraction `t` in `[0,1]`.
A use's pcurve follows the oriented face, so a reversed use pairs pcurve
fraction `t` with edge fraction `1-t`. Outer loops wind counter-clockwise
about the oriented face normal and inner loops clockwise. Cylinder UV is
(angle, axial height).

## Exact combinatorial checks

* Every vertex, edge and face is used; each face belongs to exactly one shell;
  no face, loop or shell is empty.
* Each loop closes in 3D vertex order.
* Within each shell, every edge has exactly two uses with opposite senses. One
  use is `free_edge`, and more than two is `non_manifold_edge`. A seam is two
  uses in one face. An edge used by two shells is `edge_across_shells`.
* Each shell is face-connected. At each vertex, the graph of edges joined by
  consecutive uses around that vertex must be connected, so pinch vertices are
  rejected. Each shell's Euler characteristic `V - E + 2F - L` must be even
  and at most 2.

## Certified geometric checks

A geometric check passes only on a certified upper bound `<= tol`. It fails
only on a certified lower bound `> tol`. Otherwise it reports the matching
`uncertified_*` issue. No decision rests on unverified floating point.

* **Definitions.** Finite values, nonzero line length (`|b-a| > tol`), radius
  `> tol`, and `0 < |sweep| <= 2π`.
* **Vertices.** Each edge end is within tolerance of its vertex.
* **Curve on surface.** Over a whole use, `D(t) = C(t_edge) - S(P(t))` is a
  harmonic sum: `A0 + A1 t + Σ_ω (C_ω cos ωt + S_ω sin ωt)`, with exact
  rational frequencies. This covers line and arc edges against plane
  line/arc pcurves and cylinder line pcurves. Equal frequencies are combined
  first, so a correctly matched arc cancels exactly. Then
  `sup |D| <= sqrt(|A0|² + |A0+A1|²) + Σ sqrt(λmax(Gram(C_ω, S_ω)))`, since the
  affine part is convex and each harmonic term is an ellipse. A failure is
  certified when one of 33 sample points lies beyond tolerance. An arc pcurve
  on a cylinder is not harmonic and can only be certified as a failure.
* **UV closure.** Consecutive uses meet in UV within tolerance, with cylinder
  angle differences scaled by the radius.
* **Loop winding.** Loops are closed exactly by straight chords between
  consecutive uses (gaps certified to be within tolerance). Twice the signed
  area is the closed form of `∮ u dv - v du`, including the chords. Its sign
  must be positive for an outer loop and negative for an inner loop, relative
  to the face orientation.
* **Inner loops.** The first point of each inner loop must lie inside the outer
  loop. A `+u` ray uses a half-open crossing rule. Arcs are split at their
  `v` extrema `π/2 + kπ`, using a certified `π`, so each piece is monotone.
* **Shell orientation.** Green's theorem gives the signed volume. Each use
  contributes `∫ G dv` with `∂G/∂u = S·(S_u × S_v)`. For a plane,
  `G = u (o·(x×y))`. For a cylinder,
  `G = r[a_y sin u + b_x cos u + r·det(x,y,n)·u]`, with `a_y = o·(y×n)` and
  `b_x = o·(x×n)`. The outer shell's volume must be positive and each
  cavity's negative.
* **Cavities.** A cavity vertex must lie inside the outer shell and outside
  every other cavity. Rays from that exact point in up to eight fixed integer
  directions are intersected with each face by exact Cramer solves (planes)
  or an exact quadratic (cylinders). A direction counts only when every hit is
  certifiably more than twice the tolerance from its face's boundary, and a
  start point on a face plane is certified outside that face. Adjacent faces
  only meet within tolerance, so without this margin, a ray through a shared
  edge could fall between the two faces' closed UV regions. With it, the
  parity is the same for any watertight surface within tolerance of the
  faces.

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
  explicit margins far above mpmath's working precision. `generate_brep_fixtures.py
  --check` rebuilds 54 cases from an independent prism builder: 19 valid
  solids and 35 mutations. The valid solids include holes, convex and concave
  arcs, full circles with seams, one and two cavities, rotated and
  far-translated copies and a millimetre-scale box. Two mutations remain valid
  (a shift below tolerance, and a looser tolerance absorbing a larger shift),
  so 21 cases are valid. `brep_validation.rs` requires Rust's complete sorted
  issue list to equal the reference's for every case.
* The existing prism suites (`invariants`, `occt_regression`, `modeling`)
  build every solid through the new validator.

The nineteenth fuzz target, `brep_validation`, builds valid star-outline
prisms with no hole, a round hole, a square hole or an inverted box cavity.
Scales range over `2^±10` with random frames and offsets. The cavity is merged
and inverted by fuzz-crate code, not by a kernel builder. The base must be
valid. Then one of sixteen mutations must produce its predicted issues:

* an exact report for local changes: an extra vertex, edge, empty shell or
  empty loop, a bad edge reference, an inverted outer shell or an uninverted
  cavity
* a required issue on the mutated entity for the others: a face dropped from
  or repeated in its shell, a flipped use, a moved vertex, a shifted pcurve, a
  flipped face or swapped outer and inner loops
* for a pcurve shifted by `10·tol`, a clean report at `100·tol`

Every report must be deterministic, duplicate-free and identical to
`from_parts`. A local 300-second development campaign (dirty tree based on
`12098ee3`, AddressSanitizer, standard 20-second/2 GiB limits) ran 13,807
mutation executions with 4,408 coverage edges and a 762 MB RSS peak. It
produced no artifacts.

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
never read as a native verdict.

## Native comparison bridge

`compare_brep.py` first verifies the SDK manifest, loaded toolkits and the
unchanged pre-implementation capture. The generator uses correctly rounded
trigonometry (mpmath, rounded once to binary64) so that every host produces
the same bytes. Platform libm results differ in the last bit: the first Linux
CI run regenerated different inputs from the macOS capture. The capture itself
used macOS libm, so the regenerated native inputs must match it token for
token, with numbers within `2^-50`. Only components of the triangle and
hexagon moved, by at most `1.2e-16`, and all 54 expected reports are
unchanged. The bridge also checks byte-identical regeneration of the
fixtures, and that the Rust probe's issue lists equal
the independent reference's. A case matches when the verdicts agree and every
Rust issue class with a BRepCheck counterpart has one of its corresponding
native statuses. For example, `pcurve_off_edge` corresponds to
`InvalidCurveOnSurface` and `shell_orientation` to `EnclosedRegion` or
`SubshapeNotInShape` (see `CORRESPONDING` in `brep_reference.py`). Each native
process has a 120-second deadline. Timeouts, crashes and malformed output can
never be reviewed.

On macOS, 44 cases match and 8 are reviewed differences, with no failures.
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
classification and the review fingerprint rules. Linux observations need
their own review records if their output differs.

## Limits

2D loop self-intersection, face/face intersection (which OCCT's BRepCheck does
not check either), tolerance healing, per-entity tolerances and operation
history are out of scope. Validation cost is not bounded by explicit work
limits. The curve checks are linear in the number of uses, but the
combinatorial passes use ordered maps, and containment is linear in faces per
ray. This validator certifies the supplied boundary; it does not make an
invalid import valid.
