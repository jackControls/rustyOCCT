# Splitting a solid by a plane

S8 of `REVIEW_NOTES.md` splits a solid by an arbitrary plane: the first
general topology-changing algorithm, and the rehearsal for S9's Booleans.
This document describes what is implemented; the decisions are in
`REVIEW_NOTES.md` (S8).

## Contract

`Solid::split_by_plane(operation, plane)` takes a `Frame3` (the plane through
its origin normal to its normal) and returns the pieces with their sides
(`Side::Below`, against the normal, first) and the operation's history. A
plane missing the solid, touching it (a vertex, an edge, a ruling of an arc
wall) or lying in one of its faces returns the solid itself, every entity
`Unchanged`. A side may hold several pieces (a U's two prongs). Each piece is
a validated solid; an error is one of:

* `OutOfDomain`: a solid of a later sub-step (a split piece split again),
  or a spline profile segment lying along the plane's trace.
* `Degenerate`: the split would leave an edge or a piece thinner than the
  resolution (a crossing within the resolution of a vertex, a plane within
  binary64 of a cap, a vertex within the resolution of the plane, a plane
  whose traces on the two caps' planes lie within the resolution of each
  other: parallel to the axis over the prism's height to binary64), or
  pinch a piece (a plane tangent to a hole, or touching the profile at a
  vertex, inside the solid; an oblique plane touching a cap's arc or spline
  edge between its ends; a plane tangent to a spline where it crosses it,
  or at a knot): a profile cannot hold a hole touching its boundary.
* `ComputationLimit`: a certified comparison it could not decide (an arc's
  crossing at its end).
* `PrecisionLoss`: a new point the coordinates cannot resolve (the centre
  of an ellipse where a steep plane meets a wall's axis far from the
  solid).

## Prisms (S8a.1)

In the prism's frame coordinates the plane is `a u + b v + c w + d = 0`, the
coefficients exact rationals of the stored data (the frame's axes as stored).

* **Whole.** The solid lies on one side exactly when the function's extremes
  over its caps' outer boundaries (at each vertex exactly, at each arc's
  support point `f(centre) +- r |(a, b)|` by the exact sign of a quadratic
  surd, when that point lies inside the arc) do not change sign.
* **Normal to the axis** (`a = b = 0`): M3's height split at `w = -d / c`
  rounded to binary64, its history `HeightSplit`.
* **Parallel to the axis** (`c = 0`): the pieces are the prisms of the
  profile's pieces cut by the line `a u + b v + d = 0`. Every vertex's side is
  exact; a line segment crosses where its ends' signs differ (the rational
  crossing rounded); an arc or a circle meets the line twice, tangentially or
  not by the exact comparison of `f(centre)^2` with `r^2 |(a, b)|^2`, its
  crossings' places on the arc decided by certified angles. The pieces of the
  boundary take their sides from their ends or their midpoints; the chords
  of the line inside the profile (classified at an off-centre point) join
  them; each side's cycles are traced (the next piece at a vertex the first
  clockwise from the incoming one), counter-clockwise ones becoming outer
  boundaries and clockwise ones holes. Each piece's prism is renamed from the
  input by provenance: a whole wall, edge or vertex keeps its id
  (`Unchanged`, or `Modified` when its stored frame or bounding ids changed:
  a hole's wall now on an outer boundary, a wall whose vertical edge was
  split in two), parts are `Split` children with canonical ordinals, a
  vertex on the plane is split into one copy per piece, cut walls are
  `Generated` from the caps and the walls their chord meets, cut edges from
  the cap they cut, new vertical edges from the wall they cut and new cap
  vertices from the cap edge they cut (history `PlaneSplit`).
* **Oblique** (S8a.2, `solid/split/oblique.rs`): over a profile point the
  material under the plane along the axis spans `[low, min(high, g)]`,
  `g = -(a u + b v + d) / c`, and the material over it `[max(low, g),
  high]`. A lower piece's footprint is a piece of the profile's exact
  section by the plane's trace on the bottom cap's plane (where its height
  vanishes: those chords are edges where its bottom cap meets its cut
  face), and its top is that footprint's section by the trace on the top
  cap's plane: top-cap faces where `g >= high`, cut faces on the plane where
  `g <= high`; an upper piece is the same with the ends exchanged. Walls run
  over the footprint's boundary from its flat end to its creased one: a
  planar wall's top is a line on the plane, a cylindrical wall's an ellipse
  arc (`Curve3::EllipseArc`, OCCT's `Geom_Ellipse`: its centre where the
  plane meets the wall's axis, its major axis the plane's steepest ascent,
  `major = r |m| / |c|`, `minor = r`, its angle the circle's less the angle
  of `(a, b)`) whose pcurve is the graph `v = a0 + a1 cos u + a2 sin u`
  (`Curve2::Sinusoid`), and on the cut face an axis-aligned
  `Curve2::EllipseArc`. A plane touching a cap's circle at a point leaves a
  vertex there, where the wall's height vanishes (its loop runs once round
  the cylinder from it and back, not winding). Every side, crossing and
  tangency is decided exactly by the two sections; new vertices, heights,
  ellipses and pcurves are rounded from exact values, and each piece is a
  general body (`Topology::from_parts_named`: measured, then the whole
  contract checked). Names follow provenance: an input entity whole in one
  piece keeps its id (`Unchanged`, or `Modified` when its stored geometry or
  bounding ids changed: a cap circle given its touch vertex), one in several
  pieces or in parts is `Split` (vertices and edges lying in the plane into
  copies), cut vertices, edges and faces are `Generated` from the input
  entities they cut. A piece keeps its construction (the prism's profile,
  the plane in its frame, its index), so a rigid motion rebuilds it exactly
  and it classifies points by its footprint, the height range and the
  plane's side. It tessellates within the request (the ellipse's curvature
  bounds its segments); `.brep` writes its planar pieces and its ellipses on
  planes, not a cylinder's section, which OCCT has no analytic pcurve for
  (`Unwritable`).

## Spline prisms (S8b.3)

`solid/split/spline.rs` decides a spline profile segment against a plane's
trace `a u + b v + d = 0`; the prism builders above then take its pieces.

* **Meetings.** On each exact Bézier arc of the (nonrational) spline the
  trace's function is a polynomial with rational coefficients in the arc's
  parameter; its roots in the arc are isolated exactly with their
  multiplicities (`polynomial::real::isolate`). A simple root is a
  crossing, at the curve parameter rounded to binary64 from its refined
  isolator (a root at an interior knot at the knot itself); an even one a
  touch (a tangency: the segment keeps its side; a touch inside the solid
  pinches a piece as an arc's does); an odd one of multiplicity three or
  more, or a tangency at a knot, is `Degenerate`. Each piece's side is the
  exact sign at a rational parameter between two distinct roots, never at a
  rounded crossing. Whether the prism lies on one side takes, per spline,
  whether the function at the caps' heights takes each strict sign: at the
  arcs' ends and between their roots, exactly.
* **Pieces.** A piece of a spline segment is its curve restricted to the
  rounded parameters by Boehm's knot insertion in rationals, its poles
  rounded: it keeps the curve's parameter, and its ends are the rounded
  crossing points exactly (the section's new vertices). A traced piece's
  boundary is then a path of lines, arcs and spline pieces, validated as
  every profile is (S8b.1's screen), and its prism has exact spline walls
  (S8b.2).
* **Oblique.** A wall over a footprint's spline piece is the degree-(p, 1)
  wall of that piece. Its flat end and a crease part on a cap's side are
  the piece's lifts; a crease part on the plane is the piece's affine image
  (each pole at the plane's height over it; the part's ends at their
  vertices' heights), and its pcurve on the wall is exact: `u` the curve's
  own parameter (the spline whose poles are the Greville abscissae, rounded)
  and `v` the height over the wall's low end at each pole. The rounded
  Greville abscissae leave the pcurve's hull across the wall's knot lines by
  a rounding step: the exact deviation bound takes the patch holding the
  piece but for such a sliver, plus a certified bound of the neighbouring
  patch's departure from it there (both re-expressed exactly over the
  sliver's box); Green's exact integrals take a patch whose box the piece's
  curve provably keeps to (each coordinate against each bound, exactly
  nonnegative), where a crease nearly touching a cap takes its control
  polygon out of it.
* **Names and records.** Roles and ids are the arcs': a restricted spline
  edge or wall is a `Split` child whose support is its parent's (the
  independent check compares a spline edge's exact Bézier arcs with its
  parent's over its range, a wall's rows with its parent's extrusion, a
  planar piece's spline boundary by its poles). An edge whose span runs
  against its curve (a profile given clockwise, a hole's spline now on an
  outer boundary) is written to `.brep` along its curve and used the other
  way.

## Sheets and wires (S8e)

`Body::split_by_plane(operation, plane)` (`body/split.rs`) splits S6's
bodies, planar sheets and closed wires, by the plane's trace
`a u + b v + d = 0` on the body's plane (exact rationals of the stored
frame and plane).

* **Parallel.** A plane parallel to the body (`a = b = 0`) returns it, on
  the side it lies on, `Below` in the plane. A plane within the resolution
  of the body's plane over the whole body (`|d| + |(a, b)| R` within it,
  `R` the profile's reach from the frame's origin: a tilted frame's normal
  normalized again) that would split it is `Degenerate`, its pieces
  thinner than the resolution along the plane's normal.
* **Sheets.** The profile's section by the trace (the prisms' exact
  `Section`, lines, arcs, circles, holes and splines) gives the pieces;
  each is a planar sheet of its profile on the body's frame. Whole edges
  and vertices keep their ids (`Unchanged`, or `Modified` when their
  geometry or bounding ids changed), parts are `Split` children with
  canonical ordinals, a stored vertex on the trace reached from both sides
  is split into a copy per piece, each chord is a cut edge `Generated` from
  the face and each crossing's vertex from the edge it cuts; the face is
  `Split` into one child per piece. A trace missing or touching the face
  returns it; one pinching it (tangent to a hole inside it) is `Degenerate`.
* **Wires.** The boundary's exact arrangement with the trace (its pieces
  and sides, without chords or pinch checks) is cut into maximal runs on
  one side, an edge along the trace joining the run it continues in stored
  order; each run is an open wire (a path of segments with two free ends,
  `Body::path`, built by `Topology::open_wire`). Whole edges keep their
  ids, parts are `Split`, a stored vertex where two runs meet is split into
  one copy per run, and each crossing gives a cut vertex in each run
  `Generated` from the edge it cuts. A wire touching the trace without
  crossing it returns it. An open wire split again is `OutOfDomain`.
* **Measures.** A spline wire's length and centre are certified by
  quadrature over each exact Bézier arc (the speed as a Taylor jet's square
  root, D13's order, relative widths and depth), so a wire's pieces' lengths
  add up.

## Cones and spheres (S8c.1)
## Cones and spheres (S8c.1)

`solid/split/revolved.rs`. In the solid's frame the plane is `a u + b v + c w
+ d`; which side the solid lies on is exact (the function's extremes over a
cone are over its end circles, `c w + d +- r |(a, b)|`, an apex at `c w +
d`; over a zone at the sphere's extreme points `d +- R |m|` when their
heights lie between its ends, otherwise over its end circles: each a
quadratic surd's sign).

* **Normal to the axis** (`a = b = 0`, or within a quarter of the
  resolution over the solid's widest circle: a tilted frame's stored axes):
  the cut at `w = -d / c` rounded; the pieces are the same primitive, a cone
  or frustum between exact heights with the radius there rounded, a zone
  between latitudes with `asin(w / R)` rounded, so they keep the
  primitives' exact queries.
* **A whole sphere, any plane**: two caps on a frame whose axis is the
  plane's normal, cut at the plane's signed distance from the centre.
* **History**: the wall and the region `Split` (below first), the ends a
  piece keeps the input's (`Unchanged`, or `Modified` by a rounding), the
  cut disc and ring `Generated` from the wall, a whole sphere's caps' poles
  `Generated` from it (role `Pole`). The history checker's support test
  knows cones (axes parallel, half-angles equal, the piece on the whole's
  surface within the tolerance) and spheres (centres and radii within it).
* **Containing the axis** (S8c.2, `solid/split/meridian.rs`; `c = d = 0`,
  or within a quarter of the resolution over the height): each piece is a
  half over the angles on its side, from where the plane's trace leaves the
  axis to its opposite. Its wall lies on the input's own surface, bounded
  by the end circles' halves (arcs on the input's rings) and two meridians
  (a cone's rulings, a sphere's great-circle arcs, `u`-constant pcurves);
  the end discs' halves close with chords along the trace; the cut face in
  the plane is bounded by the meridians and the chords, and an apex or pole
  ends both meridians (the wall's loop passes it). Each half is a general
  body (`Construction::Half`: the primitive, the plane in its frame, its
  index, so a rigid motion rebuilds it and it classifies points). Names:
  the wall, region, discs and rings `Split` into one child per half, an
  apex or pole into a copy per half; chord ends `Generated` from their
  ring, chords from their disc, meridians from the wall, the cut face from
  the wall and the discs.
* Other planes cut a cone in a conic and a zone in a circle (S8d.2, below).

## Conic sections (S8d.2)

`solid/split/conic.rs`: a cone, frustum, zone or cap by a plane neither
normal to its axis nor containing it. Exact on the stored data: whether the
plane misses, touches or crosses each end circle (`(c w + d)^2` against `r^2
(a^2 + b^2)` on the stored radius), whether it passes through an apex, a
pole or a frustum's virtual apex, which conic it cuts from a cone, and
whether a zone's closed circle winds round the axis (`d^2 < R^2 c^2`).

* **Sections.** An ellipse, hyperbola or parabola on a cone (`EllipseArc`,
  `HyperbolaArc`, `ParabolaArc`), the two rulings through a frustum's virtual
  apex (or within the resolution of it), a circle on a sphere
  (`CircularArc`), all explicit (D13); cut at the rims' crossings into the
  arcs inside the solid (`MATHEMATICS.md`). Their pcurves are `Projection`s
  on the wall and exact on the cut plane (an ellipse's `EllipseArc` on a
  frame sharing its axes, a circle's arc, a hyperbola's or parabola's
  projection).
* **Pieces.** Each is a general body (`Construction::Half`, rebuilt the same
  way under a rigid motion): its part of the wall, bounded by its rim arcs
  or rings and the section's arcs (loops chained through their vertices and
  walked on the cover, a wall wound once about the axis when it holds a
  ring, round an apex or pole on its side as a vertex loop, with a hole where
  a zone's circle does not wind round the axis), its end discs or their
  parts closed by chords, and the cut face.
* **Configurations.** A closed section (no rim crossed); one crossing one
  rim twice (a tongue); two arcs joining the rims; a closed section touching
  one or both rims keeps a vertex at each touch (the wall a bigon, as S8a.2's
  touch). `Degenerate`: a plane touching one rim while crossing the other or
  touching a rim with a zone's side circle (a wall pinched mid-loop), within
  the resolution of a rim's tangent, of an apex or pole, or of a frustum's
  virtual apex where the conic's axes collapse, a cone section within
  binary64 of a parabola (its axes over a million times the solid), or a
  piece thinner than the resolution. `OutOfDomain`: a plane through an apex
  or pole off the axis (rulings through the apex, a circle through the
  pole).
* **History.** The wall, the region, and each crossed rim and disc `Split`
  into one child per piece (below first); a rim, disc, apex or pole whole in
  one piece keeps its id (`Modified` when a touch gives it a vertex,
  otherwise `Unchanged`); crossing and touch vertices `Generated` from their
  rim, chords from their disc, section edges from the wall, the cut face
  from the wall and the discs it crosses.
* Pieces are not written to `.brep` yet (a projection pcurve is
  `Unwritable` until D13's interchange approximation).

## Tori (S8d.1)

`solid/split/torus.rs`, whole tori. The plane misses or touches the tube
when its distance from the core circle is at least the tube's radius
everywhere, `|d| - R |(a, b)| >= r |m|`, decided by squares.

* **Normal to the axis** (within a quarter of the resolution over the outer
  equator): the cut at `h = -d / c` meets the tube in two parallels at the
  tube's angles `asin(h / r)` and `pi - asin(h / r)` (rounded). Each piece is
  a band of the tube between them, wound once about the axis, and the planar
  annulus between the parallels: a general body on the input's own torus
  surface (`Construction::Half` with a torus), its pcurves `v`-constant
  lines. A torus v-segment of the S3 construction is the revolved region
  between the tube's arc and the axis, so it is not such a piece.
* **Containing the axis**: two half-turn wedges of the S3 construction, on
  frames whose x axes point along the plane's trace and against it.
* **History**: the wall and the region `Split` (below first), the cut
  faces and their circles `Generated` from the wall. The history checker's
  support test knows tori (axes parallel, centres and radii within the
  tolerance).
* **Any other plane** (S8d.3, `solid/split/spiric.rs`): on the torus `F = C
  + W cos(v - psi)` (`alpha = a cos u + b sin u`, `C = R alpha + d`, `W =
  r |(alpha, c)|`), so the section exists over `u` where a quadratic in
  `alpha` is non-negative; its roots against `+-|(a, b)|` decide exactly
  (surds on the stored data) between two loops about the axis (the pieces
  tube bands, their cut face a planar annulus), one contractible loop (a cap
  and the torus less that disc, each with a planar disc) and two loops about
  the tube (C-shaped pieces, two discs each); two caps on one side are
  `OutOfDomain`, a tangency or a plane within the resolution of one
  `Degenerate`. The loops about the axis are `Curve3::Section` graphs over
  `u` for a whole turn, those about the tube over `v`; a cap's loop is two
  graphs over `u` and two over `v` round its turning points, joined where
  the slope is one, so every edge is analytic. Pieces are general bodies
  (`Construction::Half` with the torus) with `Projection` pcurves; the wall
  and region `Split`, everything else `Generated` from the wall.
* Tori other than whole ones: `OutOfDomain`.

## Evidence

* **Independent reference.** `split_reference.py` (mpmath) gives each side's
  volume and centre by slicing the profile at each `u`, integrating the
  clipped height exactly between breaks (Simpson on the polynomial pieces)
  and over `u` with `mp.quad` between every break, and its area as its caps',
  walls' and cut faces' by the same slicing; `generate_split_fixtures.py
  --check` writes 26 cases (planes across caps and walls at angles, normal
  and parallel to the axis, through vertices and edges, tangent to an arc,
  in a cap, missing; holes, a U, a circle's whole wall) and each prism's
  stored axes (`stored_axes`, checked bit for bit).
* **Native.** `occt_split_oracle.cpp` splits each prism with
  `BRepAlgoAPI_Splitter` and a planar face; its capture was taken before any
  kernel code (`fixtures/occt-split-preimplementation`): all 26 within 2e-8
  of the reference (BRepGProp's accuracy on elliptic faces). Against the
  kernel three reviewed count differences: OCCT splits an arc wall along a
  tangent ruling, and keeps a seam in a cut cylinder (parallel to the axis,
  and through a plane touching both caps' circles).
* **Revolved solids (S8c).** `split_reference.revolved_rows` slices each
  cone, frustum, sphere or zone into discs cut by the plane's line and
  integrates their areas and moments along the axis (areas: the ends'
  parts, the lateral surface's angle below, the cut face's chords);
  `split-primitive-cases.txt` holds 23 cases and
  `occt-split-primitive-preimplementation` their `BRepPrimAPI_MakeCone` and
  `MakeSphere` splits captured before any kernel code, all within 2e-8 but
  two reviewed BRepGProp errors on spheres cut at an angle. The kernel is
  inside the reference on the 16 cases S8c.1 splits, its counts OCCT's.
* **Conic sections (S8d.2).** `split-conic-cases.txt` holds 12 more cases (a
  closed ellipse, a parabola, hyperbolas across one rim and both, a tongue at
  the top, the rulings through a virtual apex, an ellipse touching both
  rims, tilted cones and zones, a zone's side cap, a cap with its pole,
  circle arcs across both rims) beside S8c's three captured before any code;
  `occt-split-conic-postimplementation` holds their native splits, captured
  after `conic.rs` existed and recorded as such. The kernel is inside the
  reference on all 15; three match OCCT and twelve are reviewed count
  differences (OCCT's seam at `u = 0` splits a strip, tongue, bigon or side
  cap into two faces, or a rim or section arc it crosses), one of them also
  BRepGProp's error on a cap cut at an angle.
* **Spiric sections (S8d.3).** `split-spiric-cases.txt` holds 13 tori (bands,
  caps and C-shaped pieces, both frames) with the reference
  (`split_reference.torus_rows`) and `occt-split-spiric-preimplementation`,
  captured before any kernel code for them; with S8d's four spiric cases the
  kernel is inside the reference on all 17. OCCT approximates the sections by
  B-splines (BRepGProp's volumes and areas off by 1.6e-7 to 7.6e-4) and keeps
  the torus's seams: 17 reviewed differences.
* **Spline prisms (S8b), before the code.** A path's segment may be a
  nonrational B-spline (`B` in an `S` row). `split_reference.py` slices its
  exact Bezier pieces (crossings on `x`-monotone runs; breaks at the pieces'
  ends and extremes and the clip lines' roots), checked against polygons
  (straight splines) and Green's theorem in exact Fractions, within 1e-40;
  `split-spline-cases.txt` holds 17 prisms (a bulge, a blob of four cubics,
  a three-span wave, a stadium with a cubic end, a lens-shaped hole; planes
  normal, parallel through a knot's point, two joins, tangent, oblique, in a
  cap, both frames) and `occt-split-spline-preimplementation` their splits
  with `Geom_BSplineCurve` edges: 10 within 2e-8, seven reviewed BRepGProp
  errors (trimmed spline walls to 2.6e-7; faces bounded by the three-span
  wave to 1.8e-3, though OCCT's edges enclose the reference's area). The
  kernel (S8b.3) is inside the reference on all 17; one more reviewed
  difference: OCCT splits the bulge's cap edges where a plane touches them
  (`bulge_tangent`), the kernel returns the prism.
* **Sheets and wires (S8e), before the code.** A case may be a face or wire
  body (`make face` or `make wire` in place of `offsets`; natively a `make`
  row in place of the prism vector). `split_reference.planar_rows` gives a
  sheet's sides as a prism of height one cut parallel to its axis (area and
  first moments by the same slicing) and its perimeter as its boundary on
  the side plus the trace bounding it, and a wire's sides from its
  boundary's pieces cut exactly where the trace crosses or touches them (a
  piece along the trace on the side of the piece before it in stored
  order): rows `side S area perimeter cx cy cz` and `side S length 0 cx cy
  cz`. Checked within 1e-38 against closed forms (a square's halves, circle
  arcs and segments, a U's and a holed square's perimeters, a quadratic's
  length), exact half-plane clipping of polygons, straight splines, Green's
  theorem and the sides' sums. `split-sheet-cases.txt` holds 25 cases
  (crossing, through vertices and a knot's point, along an edge between
  prongs and between runs, tangent to an arc, a hole and a spline, missing,
  parallel, containing; lines, arcs, circles, holes and S8b's splines; both
  frames) and `occt-split-sheet-preimplementation` their
  `BRepAlgoAPI_Splitter` splits, the probe grouping a split wire's edges
  into runs (OCCT keeps it one wire): every piece valid, sides and runs the
  reference's, 23 within 2e-8 and two reviewed BRepGProp errors (a face
  bounded by two spans of the wave, 1.3e-3 in area; `LinearProperties` on a
  cubic, 3.5e-6 in perimeter). The kernel (S8e) is inside the reference on
  all 25; five more reviews: OCCT splits edges where the plane touches a
  sheet's arc and a wire's spline, keeps circles' seam vertices, and splits
  the sheet whose hole the plane touches, which the kernel refuses.
* **Kernel.** `tests/split.rs`: every side's sums of the kernel's enclosures
  contain the reference's volume, area and moments; histories pass the
  independent check, cover every input entity and repeat exactly; oblique
  pieces move rigidly with their ids, classify, tessellate and write where
  OCCT has records; spline pieces (parallel, knot, oblique and holed) move,
  tessellate and round-trip through `.brep`; sheets and wires match the
  reference and their histories check; `compare_split.py`: all 128 cases
  inside the reference, 70 matching the native counts and 58 reviewed.
* **Fuzzing.** The `split` target cuts rectangles, regular polygons,
  stadiums, U shapes and holed squares in two frames with planes chosen
  degenerate on purpose (oblique ones through a cap's vertex or touching a
  cap's circle or arc); volumes add up, pieces lie on their sides and move
  rigidly with their ids, and the split's own history check runs in its
  debug build. Cones, spheres, zones, caps and tori (a first byte of 224 on)
  are cut by nine plane modes, S8d.2's among them (touching a rim, parallel
  to a ruling, through a frustum's virtual apex). Spline prisms (a first
  byte of 192 to 223: a bulge, a three-span wave, a lens hole given either
  way round) take the same modes, the tangent ones at a spline's apex.
