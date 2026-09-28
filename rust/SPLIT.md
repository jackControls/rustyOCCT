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

* `OutOfDomain`: a solid of a later sub-step (S8b: spline prisms; S8c:
  primitives; a split piece split again).
* `Degenerate`: the split would leave an edge or a piece thinner than the
  resolution (a crossing within the resolution of a vertex, a plane within
  binary64 of a cap, a vertex within the resolution of the plane, a plane
  whose traces on the two caps' planes lie within the resolution of each
  other: parallel to the axis over the prism's height to binary64), or
  pinch a piece (a plane tangent to a hole, or touching the profile at a
  vertex, inside the solid; an oblique plane touching a cap's arc edge
  between its ends): a profile cannot hold a hole touching its boundary.
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
* Other planes cut a cone in a conic and a zone in a circle whose pcurves
  are transcendental graphs over the angle, which wait for D13's procedural
  edges (S8d): `OutOfDomain`.

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
* **Kernel.** `tests/split.rs`: every side's sums of the kernel's enclosures
  contain the reference's volume, area and moments; histories pass the
  independent check, cover every input entity and repeat exactly; oblique
  pieces move rigidly with their ids, classify, tessellate and write where
  OCCT has records; `compare_split.py`: all 26 cases inside the reference,
  23 matching the native counts and 3 reviewed.
* **Fuzzing.** The `split` target cuts rectangles, regular polygons,
  stadiums, U shapes and holed squares in two frames with planes chosen
  degenerate on purpose (oblique ones through a cap's vertex or touching a
  cap's circle or arc); volumes add up, pieces lie on their sides and move
  rigidly with their ids, and the split's own history check runs in its
  debug build.
