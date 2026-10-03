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
| Loops and shells | `open_loop`, `winding_mismatch`, `ring_edge_with_vertex`, `ring_edge_open`, `seam_edge`, `free_edge`, `non_manifold_edge`, `same_sense_uses`, `radial_order_inconsistent`, `edge_across_shells`, `disconnected_shell`, `non_manifold_vertex`, `euler`, `wire_edge_with_fins`, `acorn_vertex_used` |
| Definitions | `degenerate_vertex`, `degenerate_curve`, `degenerate_surface`, `degenerate_pcurve` |
| Certified geometry | `vertex_off_curve`, `vertex_loop_off_surface`, `pole_off_apex`, `pcurve_off_edge`, `uv_gap`, `loop_winding`, `inner_loop_outside`, `shell_orientation`, `cavity_outside`, `nested_cavity` |
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
winding number `w` closes with its end `2πw` in `u` after its start. Cone UV
is (angle, arc length along the ruling), as OCCT's `Geom_ConicalSurface`,
on the same cover; the apex is at `v = -R / sin a` (S3 of `REVIEW_NOTES.md`).
Sphere UV is (longitude, latitude), the poles at `v = ±π/2`; a sphere face
without edge loops is the whole sphere. Torus UV is (longitude, the tube's
angle from the outer equator), periodic in both: a loop winding `[wu, wv]`
closes shifted by `(2πwu, 2πwv)`; a torus face without edge loops is the
whole torus. A spline edge (`Curve3::BSpline`) or pcurve (`Curve2::BSpline`)
spans its whole domain, the fraction mapped affinely onto it; a spline face
(`Surface::BSpline`) uses its own `(u, v)` and its loops do not wind yet
(R4 of `REVIEW_NOTES.md`). A spline edge or pcurve is a `SplineSpan`: the
curve over a range, its whole domain when the kernel builds it, a sub-range
when a file trims it.

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
  zero on planes and in `v`; the windings of a cylinder or cone face sum to
  zero, except on a cone or sphere face whose edge loops wind once in total
  and that has a vertex loop: its first vertex loop is the pole, closing the
  band at the apex, or at the sphere's pole on the band's material side (a
  band winding `+u` on a forward face closes at the north pole). Windings in
  `v` are allowed only on a torus, where they must sum to zero. A face
  without loops is `empty_face` unless it is a sphere or a torus.
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
* **Sheets, wires and acorns (S6).** A face with both sides in one shell is
  two-sided: an open sheet in the infinite void. Its shell has no Euler
  condition and the face no orientation flux, and a sound shell lists more
  one-sided face sides than twice its two-sided faces, so a solid's shell
  never passes as a sheet. A wire edge a shell lists has no fin
  (`wire_edge_with_fins`), an acorn vertex no edge (`acorn_vertex_used`);
  a shell listing wire edges must connect them, and a shell with neither
  faces, edges nor a vertex is `disconnected_shell`. Wire edges' curves and
  vertex gaps are checked as every edge's.

## Certified geometric checks

A geometric check passes only on a certified upper bound `<= tol`. It fails
only on a certified lower bound `> tol`. Otherwise it reports the matching
`uncertified_*` issue. No decision rests on unverified floating point.

* **Definitions.** Finite values, nonzero line length (`|b-a| > tol`), radius
  `> tol`, and `0 < |sweep| <= 2π`.
* **Vertices.** Each edge end is within tolerance of its vertex, and each
  vertex loop's vertex within tolerance of its face's surface. A pole's
  vertex must also be within tolerance of the apex or pole (`pole_off_apex`).
  The distance to a cone is the smaller of the distances to the two rulings
  of the meridian half-plane (both nappes); to a sphere, `||p - O| - R|`;
  to a torus, the distance to its central circle less `r`.
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
  on a cylinder is not harmonic and can only be certified as a failure. On a
  cone, a line pcurve is harmonic along a ruling (`du = 0`) or a parallel
  (`dv = 0`); any other pcurve certifies only as a failure. The same holds
  on a sphere, where a meridian (`du = 0`) is a great-circle arc in `v` and
  a parallel a circle in `u`, and on a torus, where a meridian is a circle
  of the tube in `v`. A ring torus (`R - r > tol`) is required
  (`degenerate_surface`). A torus's meeting with a quadric
  (`Curve3::Toric`, S9d.4b.2) needs a ring torus, a positive other radius
  (a cone's zero at its origin allowed), a window of the other angle under
  a turn and a nonzero sweep of at most a turn (`degenerate_curve`); with
  another torus (S9d.4b.2b, `other_minor` positive) a ring torus there too,
  neither a sphere nor a cone; its
  points and jets are its root inside the window (interval Newton, then
  the implicit function theorem), undecided where Newton's step does not
  close inside the window.
* **UV closure.** Consecutive fins meet in UV within tolerance, with cylinder
  angle differences scaled by the radius, cone ones by `|R + v sin a|` and
  sphere ones by `|R cos v|` and torus ones by `R + r cos v` at the larger
  end (zero at an apex or pole, so a loop may pass through a pole), torus
  `v` differences by `r`; the last fin meets the first shifted by `2πw` in
  `u` (and `2πw_v` in `v` on a torus).
* **Loop winding.** Loops are closed exactly by straight chords between
  consecutive fins (gaps certified to be within tolerance). Twice the signed
  area is the closed form of `∮ u dv - v du`, including the chords. Its sign
  must be positive for an outer loop and negative for an inner loop, relative
  to the face orientation. On a wound cylinder or cone face the loops have
  no outer/inner order: the total periodic area `-∮ v du` over all loops
  must have the face's sign, and each unwound loop the opposite one. A
  pole adds `2π · W · v_pole`, `W` the winding it closes. On a torus face
  wound in `v` the same holds for `∮ u dv`; a face wound in both directions
  is `uncertified_loop_winding`.
* **Inner loops.** The first point of each inner loop must lie inside the outer
  loop. A `+u` ray uses a half-open crossing rule. Arcs are split at their
  `v` extrema `π/2 + kπ`, using a certified `π`, so each piece is monotone.
  An unwound loop on a wound face (S8d.2: a hole in a band) lies inside
  when the signed crossings of the `+v` ray from its first point with the
  other loops' lines, sinusoids, projections and chords over every `u`
  alias, plus one for a north pole, equal the face's sign (a sinusoid's
  height at the alias; a projection's from certified pieces, each counted
  where `v` lies above the point all along it, `u` is monotone and its ends
  straddle the alias, others bisected at most 40 times: S9c.2); other
  pcurves on those loops are `uncertified_containment`. On a torus or a
  sphere a face whose first loop runs as a hole (and none winds in `u`) is
  the surface less its loops (S8d.3, S9d.1): every loop then has the inner
  sign, no outer loop holds them, and its integrals are the whole
  surface's less the loops'. An inner loop's containment is tried
  at its start and at two points along its first fin (two rings may start
  on one ray).
* **Region orientation.** The flux of `x/3` through a face is
  `-∮ v f(u) du` over its loops, closed chords included, with
  `f = S·(S_u × S_v)` integrated in `v` from 0: for a plane `f = o·(x×y)`, for
  a cylinder `f = -r(o·(x×n)) sin u + r(o·(y×n)) cos u + r²·det(x,y,n)`, and
  for a cone `f = ρ(v) h(u)`, integrated in `v` from the apex on a face with
  a pole and from 0 otherwise, and on a sphere `f = R² cos v (O·q + R)`,
  from the pole or the south pole, a whole sphere adding its north pole's
  line; on a torus `f = r (R + r cos v)(A(u) cos v + (O·n) sin v + R cos v
  + r)`, `A(u) = O·x cos u + O·y sin u`, integrated in `v` from 0, or for
  loops wound in `v` integrated in `u` from 0 and taken `∮ · dv`; a whole
  torus's flux is `6π²Rr²` with the face's sense (`MATHEMATICS.md`). A
  shell's flux sums its faces' fluxes, negated for back sides. A bounded
  region's first shell must have positive flux and every other shell
  negative; the infinite void's shells negative.
* **Cavities.** A cavity's point (its first fin's start, or the surface point
  there for a ring fin) must lie inside its region's outer shell and outside
  every other cavity. Rays from that exact point in up to eight fixed integer
  directions are intersected with each face by exact Cramer solves (planes)
  or an exact quadratic (cylinders, whole spheres; cones certified). A direction counts only when every hit is
  certifiably more than twice the tolerance from its face's boundary, and a
  start point on a face plane is certified outside that face. Adjacent faces
  only meet within tolerance, so without this margin, a ray through a shared
  edge could fall between the two faces' closed UV regions. With it, the
  parity is the same for any watertight surface within tolerance of the
  faces. On a cylinder, a hit is inside the face when the `+v` ray from it
  on the universal cover crosses the face's loops an odd number of times,
  counting every period alias of a wound loop. On a cone (S9d.3a) the ray
  meets `X^2 + Y^2 = (R + Z tan a)^2` in the frame at a quadratic's roots
  with certified coefficients, each hit on the surface's nappe (`R + Z tan
  a > 0`) tested as a cylinder's, its `u` distances scaled by its radius,
  the apex counting as one crossing more when the face closes at it (its
  loops wind once) and it lies up the ray (`sin a < 0`). A closing chord of
  no length within rounding (a ring's one fin) is clear when its start is
  farther than the margin by more than its length. A ray against a sphere
  with loops or a torus face with loops is not yet solved, so a body with a
  cavity and such a face reports `uncertified_containment` unless the ray
  misses the whole sphere or torus (no positive root: no crossing, S9d.4b.1;
  a Boolean whose result reports only this is `ComputationLimit`); a whole
  torus (no loops, S9d.4a) is crossed at its quartic's positive roots,
  counted exactly (a root at the start or a repeated one undecided). A cavity whose first
  face is a whole sphere takes its point at `(0, 0)`.
* **Continuity (R4).** Every spline edge, pcurve and face must be C1 in its
  own parameter: `edge_not_c1`, `pcurve_not_c1` and `face_not_c1`, decided
  exactly and never uncertified. A binary64 B-spline's interior knots have
  multiplicity at most its degree `p`, so only a knot of multiplicity `p`
  (every knot of degree 1, a periodic seam included) is tested, by one
  exact homogeneous removal with zero residual (`MATHEMATICS.md`); a face
  is tested along every knot line in `u` and `v`. Only knots strictly inside
  a nonperiodic domain count. A spline ring edge must be periodic
  (`ring_edge_open`). A spline edge is degenerate when every pole lies
  within tolerance of the first, a spline pcurve when every pole is the
  first.
* **Spline geometry (S4).** The vertex gaps at a spline edge's ends are
  certified from its exact end points, and consecutive spline pcurves' UV
  gaps from theirs. A use whose factors are all rational (a line or spline
  edge, a line or spline pcurve, a plane or spline surface) is decided by
  exact composition (S4b, `MATHEMATICS.md`), provided each common Bézier
  piece of a spline pcurve on a spline surface lies in one patch by its
  control points, or but for slivers at most `2^-20` of the patch across
  its knot lines, whose departure is bounded exactly and added (S8b.3), its
  control polygon past the surface's domain only where the curve certainly
  keeps to it (its homogeneous coordinate against the bound exactly
  nonnegative: a crease nearly touching a cap, S9f.1), and the composed
  degree is at most 96. Every other use
  with spline geometry (an arc, an analytic curved surface, a pcurve across
  a knot line) is decided by second-order Taylor enclosures on halved
  pieces; only a periodic spline surface leaves it
  `uncertified_pcurve_off_edge`. A spline pcurve's share of a loop's signed
  area and periodic areas are enclosed, and its crossings of a containment
  ray are counted by parity on exact halvings (S4c). The orientation flux of
  a face with spline geometry is enclosed (S4d): on a plane exactly; along
  a spline pcurve on a cylinder, cone, sphere or torus as `-∫ F du` over
  halved pieces, the surface's exact antiderivative `F` evaluated over each
  piece's box; on a nonrational spline surface whose pcurve pieces each lie
  in one patch exactly, by Green's theorem with column antiderivatives
  (S9f.1: a piece past the surface's domain by its control polygon alone,
  as above, boxed by the bound; a piece across a `u` knot line by a sliver
  at most `2^-20` of the patch, a crease's rounded identity in `u` a step
  past the line, integrated on the patch's polynomial and the slivers'
  error added, at most `2 (deg + 1) δ` of `u` travel there times `|G|` of
  both patches over the sliver's box); and
  on any other nonperiodic spline surface by strips, the antiderivative in
  `v` enclosed over 32 strips of the domain from the surface's jets.
  Otherwise the shell's orientation is `uncertified_shell_orientation`. A
  UV gap on a spline surface is measured in 3D between the exact surface
  points at the two ends (a pcurve end within `1e-9` of the domain's width
  outside it is on the boundary patch's polynomial extension, as OCCT
  evaluates it); an inexact end leaves it `uncertified_uv_gap`. A vertex
  loop on a spline surface is `uncertified_vertex_loop`. Mass properties
  are enclosed with every spline geometry except a periodic spline surface
  (S4d, F8): exactly along nonrational spline pcurves on planes, and the
  volume and moments of a nonrational spline surface exactly (their
  integrands are tensor polynomials); everything else (rational pcurves on
  planes, spline pcurves on cylinders, cones, spheres and tori, the area
  terms of a nonrational spline surface and every term of a rational one)
  by a certified eight-node Gauss–Legendre quadrature whose remainders are
  bounded by interval Taylor series (F8, `MATHEMATICS.md`; S9f.1: its
  pieces past the domain or across a knot line by a sliver taken as the
  exact path takes them, the slivers' error bounded by `|f̄|` over their
  columns), with S4d's first-order Green integrals and strips as the
  fallback where it cannot run. A line pcurve or closing chord on a
  cylinder or cone whose `du` is at most `2^-20` of its `dv` (a gap
  between two fins' rounded ends, their angles an ulp apart) is enclosed
  as `-du F` over its box rather than by the expansion in its slope's
  powers, which lost its smallness (S9f.2a: a cylinder face of a `TILT2`
  stadium in a common, its moments 10^27 wide). On the fixtures every enclosure lies within `1e-12` of its
  property's scale (at worst `7.6e-15` of the volume, `5.6e-15` of the
  area, `2.3e-14` of `V^(1/3)` for the centroid and `1.6e-13` of
  `V^(5/3)` for the inertia), where S4d's first-order enclosures were
  0.2–0.6% on the nonrational walls' areas and 5% on the rational corner's
  volume. Validation's orientation fluxes need signs only and keep the S4d
  routes.

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

## Enclosures (M5)

Every vertex, fin and face stores an enclosure: a certified upper bound on
its gaps (Contract 5 of `IDENTITY_AND_HISTORY.md`). A vertex's covers the
ends of its edges' curves and the surface of a vertex loop. A fin's covers
the deviation between the edge curve and the pcurve's image. A face's
covers the UV gaps between its consecutive fins. Each certified relation
above is decided against the stored bound first and against the resolution
only if that fails:

| Against the bound | Against the resolution | Report |
| --- | --- | --- |
| within | (not needed) | nothing |
| beyond or undecided | beyond | the geometric kind (`vertex_off_curve`, `pcurve_off_edge`, `uv_gap`, `vertex_loop_off_surface`) |
| beyond or undecided | undecided | its `uncertified_*` kind |
| beyond | within | `enclosure_unsound` |
| undecided | within | `uncertified_enclosure` |

A missing bound is `enclosure_missing`; one outside `[0, resolution]`
(including NaN) is `enclosure_exceeds_resolution`. Either way, the relation
is decided against the resolution alone. Measured bounds
(`TopologyParts::with_measured_enclosures`, and every builder) use the same
expressions in the same tiers. Each is stored one binary64 step above its
interval, and at least `2^-80`, so the checker's comparison is strict and
decides without falling back.

Evidence:

* `cell_reference.py` declares every fixture's bounds independently: twice
  its own highest certain gap plus `2^-20` of the tolerance, rounded up and
  capped at the tolerance, unless a case sets one. It implements the table
  above. The ten enclosure cases move a box vertex by 0.5 and 1.5
  resolutions, shift a cylinder cap pcurve or a box side pcurve by half a
  resolution (valid with generous bounds, unsound with a quarter-resolution
  fin bound or a tenth-resolution face bound), and remove, double or negate
  a bound. Rust's issue lists equal the reference's on all 76 cases.
* `brep-enclosure-lows.tsv` holds the reference's certain lower value of
  every gap of the 27 valid cases. The kernel's measured bounds on the same
  parts, with every declared bound removed, must lie between that value and
  the declared bound, within `2^-46` of the case's size. The reference builds
  frame axes exactly while the kernel stores them rounded, and binary64
  intervals widen with the coordinates; the allowance covers both. It is
  vacuous for rounding-level gaps and strict for the gaps the cases place
  deliberately (half a resolution, or a pcurve shifted below it).
* The `brep_validation` fuzz target removes a bound, sets one outside
  `[0, resolution]` (including NaN and infinity), or moves a vertex by half
  the resolution: exactly `enclosure_missing`, `enclosure_exceeds_resolution`
  or `enclosure_unsound` must follow. A circle prism placed up to `2^48`
  from the origin, at a tolerance of a few ulps there, must either fail
  construction or enclose every entity within the resolution. The
  `split_merge` target and `enclosures.rs` require no continued entity's
  bound to fall through a transform, split or fuse.
* Native observations (`fixtures/occt-enclosure-preimplementation`) were
  captured before any enclosure code. On every case valid on both sides,
  the kernel's measured vertex and fin bounds are not below OCCT's own
  measurements of the same gaps (`BRepLib_ValidateEdge` exact method, curve
  ends against vertices; seams excluded). They also never exceed the
  tolerance OCCT stores (T6). Where OCCT's representation (a line as point,
  direction and range) has rounding-level gaps that the kernel's segments
  do not, the same `2^-46` allowance applies. The two cases with real gaps
  compare strictly: OCCT measures `1.414e-9` and `1.414e-3`, and the kernel
  bounds are `2.0e-9` and `2.0e-3` (the harmonic bound's affine part is
  `sqrt(|A0|² + |A0+A1|²)`).

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
  `generate_brep_fixtures.py --check` rebuilds 152 cases. 54 come from an
  independent seamed prism builder (19 valid solids and 35 mutations; the
  valid solids include holes, convex and concave arcs, full circles, one and
  two cavities, rotated and far-translated copies and a millimetre-scale
  box), converted to cells. Twenty-two are cell-model cases with no seamed form: a
  cylinder parametrized from `π`, a box with a valid vertex loop, a ring
  pcurve spanning OCCT's printed period `6.28318530717959` (valid only
  through the near-frequency bound; without it the reference cannot decide),
  the same pcurve off by `10^-6` of a period, and the model's own failure modes
  (a fin shifted by a period, a flipped winding, fins swapped between edges, a
  side at the wrong shell, a vertex loop off its surface, a face without
  loops, a ring edge with one vertex and a shell its region does not list).
  Ten more place gaps just inside and just outside the resolution and
  declare enclosures that are missing, out of range or unsound (see
  [Enclosures](#enclosures-m5)). Twelve are cones (S3), built as cells by an
  independent cone builder: an apex at the top or the base, frusta
  narrowing and widening, a rotated frame and a far one, and six mutations
  (the pole moved along a ruling, off the surface or removed, a ring loop
  winding twice, a right angle and a ring pcurve shifted in `v`). Fourteen
  are spheres from an independent sphere builder: the whole sphere, both
  hemispheres, a zone, rotated and far copies, a whole sphere with an
  immersed vertex (valid), and seven mutations (the pole at the other pole
  or off the surface or removed, a ring winding twice, a zero radius, a
  whole sphere turned inside out and a ring pcurve shifted), and one
  sphere loop winding in `v`. Thirteen are tori from an independent torus
  builder: the whole torus, the outer and inner halves, a segment, two
  wedges and a far copy, and six mutations (a meridian loop winding twice or
  unbalanced in `v`, a meridian pcurve shifted, a whole torus inside out, a
  spindle torus and the inner half with its wall forward, as OCCT builds
  it). Twenty-one are spline cells (R4) on the box and cylinder: C1 and
  broken knots on a quadratic, a linear, a rational and an unclamped edge
  (whose knots beyond its domain are not tested), periodic ring edges with
  a C1 seam or a broken one, a small periodic basis whose removal needs a
  refinement first, a nonperiodic ring edge, spline pcurves, spline faces
  (C1, broken, with a vertex loop) and degenerate spline edges and
  pcurves, three spline ranges (an edge and a pcurve trimmed from a
  longer line spline, one with a corner knot outside its range), a pcurve
  running against its spline's parameter, two
  holes against a spline side (inside and outside), and the rational
  rounded corner shrunk by `2^-10` far from the origin (a fuzz finding,
  S4d). Ten are the spline models of S4 (`spline_models()`, with OCCT
  rows for `compare_brep.py --family spline`): prisms with a spline side
  and a ruled spline wall, and a stadium with spline geometry on its
  cylinder, with their mutations. Three are F8's (`quadrature_models()`,
  with OCCT rows): a narrow stadium with a spline pcurve along a cylinder's
  parallel, and the bulge and the rounded corner with a knot in their
  walls. 63 cases are valid. `brep_validation.rs`
  requires Rust's complete sorted issue list to equal the reference's for
  every case. For the nineteen valid spline cases, `cell_reference.mass_properties`
  integrates volume, area, centroid and inertia by Green's theorem in UV
  with nested Gauss–Legendre quadrature of the exact surface jets
  (`brep-spline-mass.tsv`, twenty digits). F8: generation integrates every
  case again with each quadrature interval halved on the faces whose
  integrands 24 nodes do not integrate exactly, and requires agreement
  within `1e-20` of each property's scale (the bulge's area differs by
  `1.9e-22` relative, the rest by `1e-30` or less), and the bulges' volume
  and area to equal their closed forms `20/3` and `40/3 + 8 + √2 + asinh 1`
  to the same bound. `spline_mass_encloses_the_reference` requires every
  kernel enclosure to contain the reference within that bound and the
  rounding of the printed digits, and to be within `1e-12` of its
  property's scale; `spline_parallels_integrate_as_their_lines` writes the
  parallels of a cone, a sphere zone and two tori as degree-1 spline
  pcurves and requires the quadrature's enclosures to overlap the closed
  forms' along the lines, as narrowly.
* The existing prism suites (`invariants`, `occt_regression`, `modeling`)
  build every solid through the new validator.

The nineteenth fuzz target, `brep_validation`, builds valid star-outline
prisms with no hole, a round (seamless) hole, a square hole or an inverted box
cavity. Scales range over `2^±10` with random frames and offsets. The cavity
is merged by fuzz-crate code, not by a kernel builder: its material shell
joins the body's solid region and its twin bounds a new void region. The
base must be valid. Then one of 33 mutations must produce its predicted
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
* a cone or frustum from `Solid::cone_with` (mutation 28, S3), valid with
  enclosures within the resolution, certified volume within `1e-9` of
  `πh(R² + Rr + r²)/3` and synthesized counts `2/3/2/2/1/1` or
  `2/3/3/3/1/1`, then its pole moved `1000·tol` along a ruling (exactly
  `pole_off_apex`) or along the axis (`vertex_loop_off_surface`), the pole
  removed (`loop_without_face` and `winding_mismatch`), a right semi-angle
  (`degenerate_surface`) or a ring pcurve shifted in `v`
  (`pcurve_off_edge`); a moved pole is declared at the resolution, as the
  fixtures declare it
* a whole sphere, hemisphere or zone (mutation 29) with the same checks and
  its pole moved over the sphere or along the axis, removed, a zero radius,
  a ring pcurve shifted, or a whole sphere turned inside out
  (`shell_orientation`)
* a whole torus, v-segment or wedge (mutation 30) with the same checks, a
  volume within `1e-9` of Pappus's for a full tube, and a ring or meridian
  loop winding twice (`winding_mismatch`), a ring pcurve shifted off the
  tube (`pcurve_off_edge`), the tube reaching the axis
  (`degenerate_surface`) or a whole torus turned inside out
* a star prism whose line edge, line pcurve or plane becomes a spline of
  degree 2 or 3 with its knot `1/2` repeated to the degree (mutation 31,
  R4): exactly C1 there (poles on a dyadic grid, the knot's pole their
  midpoint), it reports no continuity issue; with that pole moved half a
  grid step, exactly one more issue, `edge_not_c1`, `pcurve_not_c1` or
  `face_not_c1` (the uses' deviations set aside: the grid poles leave the
  prism's geometry)
* a valid spline fixture (the spline prisms, by exact composition, and
  the stadiums with spline geometry on their cylinder, by Taylor
  enclosures) moved by `p -> s p + t`, `s` a power of two and `t` dyadic
  (mutation 32, S4b-d), which must stay valid with a volume and centroid
  enclosing the reference's, moved (the rational corner included), then its
  first spline use's
  pcurve shifted by `1000·tol` (`pcurve_off_edge`) or by `0.001·tol` (at
  most `enclosure_unsound`), or a ruled wall reversed (`loop_winding`)

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

`compare_brep.py --family spline` runs the ten spline models of S4 against
their pre-implementation captures: the S4a rows (statuses, counts,
tolerances and deviations) and the S4d `BRepGProp` rows
(`occt-spline-properties`), both reproduced on every run; and the three
models of F8 against theirs (`occt-spline-quadrature-preimplementation`,
all four rows, captured before the certified quadrature). It has 10 matches
and 3 reviewed differences: `spline_c0_bulge` (BRepCheck has no continuity
status; the kernel reports the C0 knot), and the two shifted pcurves (the
edge's `InvalidCurveOnSurface` makes `BRepCheck_Analyzer` skip the wire
checks and mark the face `UnorientableShape`, and `Closed2d` accepts the
gaps under its 1% rule, as for the prism shifts). For the eight models
valid on both sides, the kernel's enclosures are at least OCCT's measured
deviations and at most its tolerances, and its certified volume, area,
centroid and inertia contain `BRepGProp`'s values within OCCT's error
estimate plus `1e-8` relative: OCCT's area of `spline_bulge` is `1.7e-9`
relative from the closed form while it estimates `2e-16`. With F8 the
kernel's absolute widths there are `7e-15`–`1.2e-13` for volumes and areas
and at most `4.1e-12` for the inertia's entries, so OCCT's own integration error, not the
kernel's enclosure, is what the allowance absorbs. The sheet bridge allows
the bulge's lone spline wall (`sheet_spline_wall`) `2e-8` of its row's
scale for the same reason: OCCT's area is `1.8e-8` relative above the
closed form `√2 + asinh 1`, and its centre `4.6e-9` off.

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
bridge counts, the synthesized-count checks above, 64 fixture reports and a
clean 600-second `brep_validation` campaign; the record is in
`TOPOLOGY_MODEL.md`.

## Limits

2D loop self-intersection, face/face intersection (which OCCT's BRepCheck does
not check either), tolerance healing and operation
history are out of scope. Validation cost is not bounded by explicit work
limits. The curve checks are linear in the number of fins, but the
combinatorial passes use ordered maps, and containment is linear in faces per
ray. This validator certifies the supplied boundary; it does not make an
invalid import valid. Free and non-manifold edges, wire edges and acorn
vertices are representable but rejected: every current operation requires a
solid. Faces without loops (closed surfaces), unwound loops on wound faces
whose other loops are not lines (`uncertified_containment`) and windings in
`v` wait for the surfaces and operations that need them. A `Projection`
pcurve (S8d.2) of its fin's own edge on its face's own surface deviates by
zero by definition; its integrals come from certified Taylor quadrature and
its `+u` ray crossings from certified pieces (none where `u` or `v` keeps
clear of the ray, one where `v` is monotone with its ends on either side,
others bisected; `TOPOLOGY_MODEL.md`).
