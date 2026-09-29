# Booleans of solids

S9 of `REVIEW_NOTES.md` fuses, cuts and intersects solids (the Combine
job): the faces intersected (S7), split (S8), classified by regions and
assembled into shells and regions, with complete histories. This document
describes what is implemented; the decisions are in `REVIEW_NOTES.md` (S9).
S9a is implemented (`profile/boolean.rs` with its spline meetings in
`profile/boolean/splines.rs`, `solid/boolean.rs`, S9a.2's stacks in
`solid/boolean/stack.rs`), and S9b's polyhedral Booleans in any relative
position (`solid/boolean/polyhedra.rs`: prisms, and any planar solid with
straight edges as an input, a Boolean's result or a plane's piece among
them), S9c.1's prisms with arcs in any relative position
(`solid/boolean/curved/`), S9c.2a's perpendicular cylinders in exact
frames meeting in quartics (`curved/procedural.rs`), S9c.2b.1's
cylinders with crossing axes in turned frames (`curved/turned.rs`) and
S9c.2b.2's algebraic vertices where sections cross caps' circles, with
parallel cylinders in turned frames (`curved/algebraic.rs`); S9d is not.

## Contract

`Solid::fuse(operation, other)`, `Solid::cut(operation, tool)` and
`Solid::common(operation, other)` return the result's solids (each a
maximal connected solid region, ordered by the lowest input face id it
keeps, none for an empty result) and the operation's history (`Fuse`,
`Cut`, `Common`). A result that is itself a prism of one of the inputs'
frames (one profile over one height range) is built as one, keeping the
prism's exact queries; any other is a general body built through
`TopologyParts` and validated before it is returned. An error is one of:

* `OutOfDomain`: a pair of a later sub-step (a prism with arcs against a
  solid other than a prism, a spline profile or an arc whose ends lie off
  its circle in any position, S9c), two spline segments along one curve in
  different forms or a spline span along a line, or a cavity in a result
  of several solids.
* `Degenerate`: a crossing within the resolution of a vertex, two crossings
  within it of each other, a piece thinner than the resolution, or a result
  touching itself at a point or along an edge (two solids sharing an edge,
  a hole tangent to the outer boundary) until the kernel holds non-manifold
  bodies.
* `ComputationLimit`: two cuts of a segment, or a crossing at a segment's
  end, whose order its enclosures leave undecided.
* `InvalidLabel`: inputs sharing an entity id (built by one operation, one
  solid twice, or a result and an input whose entities it keeps), which
  the history could not tell apart.

### Prisms in one frame (S9a)

Two prisms whose frames have bitwise-equal axes and whose origins differ by
a vector with binary64 coordinates in them, of lines, arcs and circles, with
holes. The Boolean is a stack of height slabs between the four caps'
heights, each slab's profile the 2D Boolean of the profiles present in it
(an absent profile is empty), consecutive slabs with equal profiles merged
(M3's fuse). S9a.1: results that are one prism (every common; a cut whose
tool spans the object's heights or misses them; a fuse of equal height
ranges, or whose slabs all have one profile). S9a.2: the other stacks,
general bodies whose caps between slabs are the regions where consecutive
profiles differ.

The 2D Boolean (`profile/boolean.rs`) arranges the two boundaries exactly:
two lines cross at their rational crossing, or overlap along one line and
cut each other at their ends; a line and a circle meet where the exact sign
of their discriminant says, at a quadratic surd's enclosure; two circles on
their radical line, or overlap as one circle. A tangency cuts nothing (the
pieces on either side of it lie on one side of the other boundary); a
meeting within the resolution of a segment's stored end is that vertex,
which then cuts the other segment; a stored vertex of one exactly on the
other's line, or within the resolution of its circle, cuts it; B's vertices
equal to A's are A's. Pieces along one line or one circle between the same
ends are shared (in the same or the opposite direction); every other piece
is classified against the other profile at an off-centre point (fraction
`0.4453125` along it, the certified `Profile::classify`; one within the
resolution of the other boundary without being shared is `Degenerate`).
The operation keeps: for a fuse each profile's pieces outside the other and
the shared pieces of the same direction; for a cut the object's pieces
outside the tool, the tool's inside the object reversed and the shared
pieces of opposite directions; for a common each profile's pieces inside
the other and the shared pieces of the same direction. The kept pieces,
directed with the region on their left, are traced into cycles: each end
starts exactly one kept piece (two would make the result touch itself
there: `Degenerate`). A cycle keeps no vertex where it does not turn
(consecutive collinear lines and consecutive arcs of one circle in one
sense are joined, each result segment listing the input pieces it holds; a
cycle left as one arc is its whole circle), counter-clockwise cycles are
outer boundaries and clockwise ones holes, each result a validated profile;
separate results must not touch.

In 3D (`solid/boolean.rs`) the tool's profile is translated exactly into
the object's frame; a common is the 2D common over the ranges'
intersection (empty when they only touch); a cut whose tool spans the
object's heights is the 2D cut over the object's, one that misses them the
object; a fuse of equal ranges the 2D fuse over them, of meeting ranges and
identical profiles one prism over their union, of an input inside the other
(in 2D and in height) that input, of profiles apart both inputs.

S9a.2 (`solid/boolean/stack.rs`) builds every other result, a stack, as
general bodies from one arrangement of both profiles: each piece knows
whether the region on its left and on its right lies in the object and in
the tool, and in each slab the result's region is a set function of the
two (`A ∪ B`, `A`, `B`, `A - B` or nothing). A piece is a wall in a slab
where the result holds one side and not the other; at each slab height the
upward faces (the result below and not above) and the downward ones (above
and not below) are traced from the arrangement as a result profile is,
without joining. Walls on one line facing one way, or on one circle facing
one way, join across a slab height and across a piece's end where only
the two meet there; a horizontal face joins nothing. An edge is a chain of
fine edges (a piece at a height, or a point's vertical over a slab) running
straight on (collinear lines, arcs of one circle, one vertical line)
between the same two faces, its ends vertices; a closed chain of arcs is a
circle. Four faces at an edge, or a face meeting itself at an edge or a
vertex, is `Degenerate` (the result touching itself). Each connected set of
faces is a solid; one whose lowest face faces up bounds a cavity, a second
shell of the solid around it and a bounded void inside it (`OutOfDomain`
beside several solids). A plane wall's frame is the prism's (x its
region-left direction, its normal leaving the material), a cylinder wall's
its circle's centre at the face's lowest height, its pcurves lines in
(angle, height) continuing round each loop and every loop on the sheet of
the first (the outer, first on any face, by its signed area). Each solid is
validated as it is built; its mass properties are the general certified
enclosure's, its classification the slabs' set function at the point's
height over both profiles' classifications (a point on a profile's
boundary taking both memberships, on a slab height both slabs), and a
rigid motion rebuilds it in the moved frame (ids kept, mass moved).

History: each result entity continues the input entities it is a part of
(a cap the caps at its height of the inputs whose material it holds, a wall
or cap edge the walls or cap edges its segment joins, a vertical edge or cap
vertex its vertex's, the region the inputs' regions); a cut's tool faces
the other way where it bounds the result, so what it gives is generated,
not continued. An input continued by one result entity alone keeps its id
(`Unchanged`, or `Modified` when its geometry or bounding ids changed), by
several each continuing it alone is `Split`; several inputs continued by one
result alone are `Merged` into it; where several inputs reach several
results (coplanar caps over several results, a wall shared in part) no
single relation says it: each such result is `Generated` from its parents
and those inputs are `Deleted`. An entity continuing nothing (a vertical
edge where two walls cross, a cap edge inside a wall where the other
input's cap meets it) is `Generated` from what it lies on; an input
continued by nothing is `Deleted`. Inputs returned whole keep every id. The
independent history check runs in debug builds and every test.

In a stack a wall continues the input walls it lies on facing their way
(any other it lies on, and a cut's tool's, it touches); a horizontal face
the input caps at its height facing its way whose region it overlaps (an
exact 2D common), touching a cut's tool's cap facing the other way; an edge
or vertex at an input's end height continues its cap edges or vertices,
and one inside an input's height range lies on that input's wall or
vertical edge and on the caps its faces lie on; a vertical edge continues
the vertical edges of its point's input vertices over their ranges and
lies on the walls through the point; a solid continues the regions of the
inputs whose faces it continues (a cut's: the object's), and a cavity's
void is generated from the tool's region.

### Polyhedral prisms in any position (S9b.1)

Two prisms of line profiles (polygons with polygon holes) whose frames'
axes differ, or whose offset or tool profile would round in S9a's frame,
are decided on their constructions' exact models (`solid/boolean/
polyhedra.rs`): a point `o + u x + v y + w n` in rationals from the stored
binary64 origin, axes, profile and heights, a cap's plane normal to `x * y`
(the stored axes are not exactly orthogonal), a wall's the plane of its
segment's direction and `n`, so every model vertex lies exactly on its
faces' planes. Each boundary face, as convex pieces (the cap's trapezoids,
the wall's rectangle), is split by every plane of the other prism's faces;
each fragment's centroid, pushed an infinitesimal step along its normal
and against it, is classified exactly against the other prism's convex
cells (the profile's trapezoids swept, closed half-spaces), and the
fragment kept, oriented to leave the result's material, where the set
function differs across it (the object's fragment once where both
boundaries lie). Kept fragments are made conforming (each edge split at
every kept vertex on it, found exactly), joined across shared edges into
maximal faces of one oriented plane, their edges joined where they run
straight on between the same two faces. A face with a vertex within the
resolution of another vertex or of an edge not ending there (on the exact
model), an edge with four faces, or two solids sharing a vertex is
`Degenerate`. Each connected set of faces is a shell: an outer one
(positive exact volume) or a cavity of the outer shell holding it (exact
ray parity). Vertices and planes are rounded once, each solid validated as
it is built. A result's rigid motion moves its stored geometry (vertices,
lines and plane frames; the pcurves in those frames unchanged) and
measures its enclosures again; its classification is the set function of
both inputs' classifications within its bounds.

S9b.2: an input other than a line-profile prism (a Boolean's stack or
polyhedral result, a plane's piece, any solid of planar faces and straight
edges) is decided on its stored geometry: its vertices as rationals, each
face cut into triangles of its own vertices in its projection on the
normal's largest coordinate plane (ears clipped, the fattest first, holes
bridged; so each triangle is exactly planar and neighbouring faces share
their stored edges exactly), or, when that fails the exact area check, its
trapezoids zipped through every corner on their sides; its side by exact
ray parity. Fragments of one stored face join back into one face (their
triangles' planes differ within rounding); a prism's fragments join by
their exact planes as before. Its entities are named by their ids.

Cost: a part is split only by the other's planes whose faces' bounding box
meets its own (the other's surface near it is all that can cross it); a
fragment is classified at a point of short dyadic coordinates inside it
(the centroid of its shortest fan triangle for a sliver), in one ray cast
when it lies on none of the other's planes; the ray casts and the
conforming step's collinearity tests take their signs from integer forms
over common denominators, after binary64 filters.

History: a face continues the input faces its fragments come from facing
their way (a cut's tool's, or one facing the other way, it touches); an
edge along an input edge (exactly) continues it, one elsewhere lies on the
faces meeting there; a vertex at an input vertex continues it, one
elsewhere lies on the input edges and faces through it; a solid continues
the regions of the inputs whose faces it continues (a cut's: the object's),
a cavity's void is generated from the tool's region.

### Prisms with arcs in any position (S9c.1)

Two prisms of line, arc and circle profiles in any frames, at least one
with an arc, are decided on their exact models (`curved/model.rs`): a
point of a prism is `o + u x + v y + w n` on its stored axes as rationals
(a cap's plane holds `x` and `y`, its normal `x * y`), an arc wall is the
prism's cylinder over its profile circle (elliptic in the world when the
stored axes are not orthonormal), a full circle is split into two arcs at a
rational point `(1 - s^2, 2 s) / (1 + s^2)` of it (each input at another
`s`, tried again at others when a meeting falls on a seam). Arcs must end
on their circles exactly.

Every vertex is exact (`curved/num.rs`, `curved/meet.rs`): an edge of one
input meets a face surface of the other in a quadratic surd `a + b sqrt(d)`
(a line against a plane or a cylinder, an arc against a plane, two
circles), kept where it lies strictly inside both the edge and the face's
region (a vertex of one on the other's surface, an edge meeting an edge,
or a tangency is `Degenerate`); the two ellipses of equal circular
cylinders whose axes cross add their two crossing points. Two faces'
surfaces meet in lines, generatrices, plane sections of cylinders (the
cylinder's `w = a0 + a1 cos t + a2 sin t` as a conic `c + a cos t + b sin t`)
or those two ellipses; each branch is split at the vertices on it and its
pieces kept where a rational point strictly between lies inside both faces.
Two cylinders meet in such curves only when both are circular in a common
measure (their frames' axes equal, or both exactly orthonormal): parallel
(two circles in a section), coaxial (one surface, or none), or equal with
crossing axes; others are S9c.2's unless their faces' bounds or their
sections' reach are certainly apart. A plane within rounding of a
cylinder's axis direction (its section's axis past `10^12` radii) is
`Degenerate`. Orders along a curve and around a vertex are exact signs of
one surd or of two (`x + y sqrt(e)`, `x` and `y` in `Q(sqrt(d))`).

Each face's pieces (`curved/graph.rs`) are traced from its edges split at
their vertices and its sections: at a vertex the next edge is the first
clockwise from the way back (exact angles about the face's normal); loops'
orientation and nesting come from their binary64 image in the face's
parameters, a hole touching another loop within `1e-9` of the face's size
refused. A piece is classified at a rational point of one of its edges
pushed into it and then off the face either way, against the other
input's exact membership (the push's first-order sign at each boundary it
lies on, a push along a cylinder's circle keeping to it), and kept as
S9b.1 keeps fragments. Faces of both inputs on one surface (coplanar caps
or walls, one cylinder) hold each other's edges within them, the edges
crossing on it adding vertices; their pieces facing one way join.

The result (`curved/assemble.rs`): pieces of one input face kept the same
way join across the edges they share (a circle's halves among them), with
coincident pieces of the other input facing the same way; edges join where
they run on along one exact curve between the same faces (a seam's
vertices dropped, a whole circle or ellipse left as a ring edge); curves
are rounded once (lines, circular arcs on the cap's arc frame, ellipse arcs
on principal axes); pcurves lie on the input faces' stored surfaces
(lines and sinusoids on cylinders where they fit the edge at the same
fractions, exact projections otherwise, lifted continuously along each
loop and each hole lifted to its outer loop's turn); solids are the shells
joined by edges, a cut's shell of the tool alone its cavity (one among
several solids `OutOfDomain`), shells meeting at a vertex `Degenerate`.
Names follow S9b.1's rules; a result is a `Polyhedron` (both inputs, the
operation, its index), classified by the set function and moved by its
stored geometry. A result with arcs given to another Boolean is still
`OutOfDomain` (S9c).

### Cylinders meeting in quartics (S9c.2a)

Two circular cylinders in exact frames (the world's axes permuted or
reversed) whose axes cross are perpendicular; unless they are S9c.1's equal
cylinders with meeting axes, they meet in a quartic (`curved/procedural.rs`).
Across the common perpendicular `e = nA x nB` each cylinder reaches
`[e_k - r_k, e_k + r_k]`; exactly: apart, touching (outside, or inside
with a node) `Degenerate`, one extent strictly inside the other's (two
rings about the inner cylinder), or overlapping in part (one loop). Each
piece is a graph over one cylinder's angle, the ruling at `(cos, sin)`
meeting the other where `A w^2 + 2 B w + C = 0` and the branch the sign of
`A w + B`: a ring over its cylinder's whole turn, a loop in four graphs,
over each cylinder's angle about that cylinder's extreme (where the other's
rulings are tangent to the curve), switched at rational points of the
thinner cylinder's angle between the two kinds of turning point, which
become vertices (`w` a surd of one field). A piece's points are placed by
its carrier's angle, its tangent is the two gradients' cross product turned
to the carrier's run, and a rational point strictly between two places is
a point of the piece; everything else is S9c.1's arrangement. Vertices on
a quartic are where an input's generatrix (an edge or a seam) meets the
other cylinder (quadratic surds); where a cap's circle meets the other
cylinder the vertex is algebraic (S9c.2b.2, below). Edges are
`Curve3::Meet` (the carrier's and the other's stored cylinders, the branch,
the angle's range: `TOPOLOGY_MODEL.md`), their pcurves exact projections
(a meeting on its own carrier evaluated by its angle). The validator
certifies a band face's holes bounded by projections and sinusoids (the
signed `+v` ray over every `u` alias, projections by certified pieces:
`BREP_VALIDATION.md`).

### Cylinders in turned frames (S9c.2b.1)

Two cylinders whose frames' stored axes are not exactly orthonormal (their
affine models elliptic in the world) and whose axes cross meet in a
quartic taken the same way (`curved/turned.rs`), its structure decided
from each cylinder's discriminant: over a carrier's ruling at `(cos, sin)`
the other cylinder gives `A w^2 + 2 B w + C = 0`, and `D = B^2 - A C` is a
quadratic form in `(cos, sin)`, a quartic in the half-angle tangent `t` of
a chart (the base turned by `2 atan t`, its antipode at infinity) whose
real roots, the carrier's turning points, are isolated exactly
(`polynomial/real.rs`'s Sturm sequences). A chart is chosen whose antipode
has `D < 0`; a cylinder with no such point has `D >= 0` all round and
carries two rings, one with `D < 0` all round meets nothing. Otherwise
every root pair of the first cylinder's `D` bounds a loop: the second
cylinder's turning points (where its ruling touches the first: `w = -B /
A`) are placed on the loops by their angle and branch on the first (binary64
views of isolators narrowed 160 bisections), the events ordered along each
loop, and a switch put at a rational angle of the first cylinder midway
between each adjacent pair of different kinds. Each run between switches
is a graph over the second cylinder's angle (about the first's turning
points; its branch and range read at the switches and a point inside) or
over the first's (about the second's, on one branch), and is verified
exactly: no root of its carrier's discriminant within its range (the
roots against rational ends; Sturm counts at surd ends, the range clear
of the chart's antipode). A repeated root, or a critical point of the
chart's quartic where `|D| < (A res / 2)^2` (two branches within the
resolution: equal cylinders with meeting axes in stored turned frames
among them, whose extents across the common perpendicular are equal
exactly), is `Degenerate`.

### Algebraic vertices (S9c.2b.2)

A cap's circle `c + a cos + b sin` meets the other cylinder where a
quartic in its half-angle tangent `t` vanishes (`curved/algebraic.rs`):
each real root `alpha` (isolated exactly) gives a vertex at `(cos, sin) =
((1 - alpha^2), 2 alpha) / (1 + alpha^2)` with coordinates in `Q(alpha)`
(the antipode of `(1, 0)`, where the quartic drops a degree, a rational
vertex), in exact frames (S9c.2a's former nested surds) and turned ones
alike. The surds `a + b sqrt(d)` of `curved/num.rs` take `a` and `b` in
`Q` or one `Q(alpha)`: elements are polynomials in `alpha` reduced by its
polynomial, inverses by the extended Euclidean algorithm (a common factor
that vanishes at `alpha` is a zero), signs exact by Sturm-Tarski at
`alpha` (`polynomial/real.rs`), a surd's sign by S9c.1's tower rule.
Numbers of two different fields (two such vertices compared) are ordered
by enclosures of their generators' isolators refined up to 480 bisections
(`2^-192` interval grid); unseparated below `1e-40` of their magnitude
they count as equal, which the arrangement refuses. Everything else
(places, orders, membership, tangents, pushes) is the arrangement's own,
now over these fields. Parallel cylinders not circular in a common measure
meet in the generatrices through the first cylinder's circle's crossings
with the other (`CylPair::Lines`): both cylinders hold them, their axes
being parallel exactly.

### Spline profiles (S9a.2)

Either profile may hold S8b's nonrational spline segments
(`profile/boolean/splines.rs`). A spline meets a line where the line's
equation has a root on one of its Bézier arcs (S8b.3's `meets`), a circle
where `(x - c_x)^2 + (y - c_y)^2 - r^2` does (degree `2p`, the same
isolation), and another spline where each arc's parameter is a root of the
resultant of the other arc's implicit equation on it: `Res_t(x(t) - X,
y(t) - Y)` (a Sylvester determinant, exact) evaluated on the other arc's
points at `p q + 1` rational parameters and interpolated exactly, both ways,
each root paired with the one root whose certified point box meets it
(refined until one partner or none remains; two is `ComputationLimit`). A
crossing is at the spline's parameter rounded to binary64, the vertex the
curve's exact point there rounded; a root of even multiplicity is a
tangency and cuts nothing, one of two splines is refused (`Degenerate`),
and an identically vanishing resultant or equation (a spline along another
curve) is `OutOfDomain` unless the two segments are one curve, equal or
with reversed poles and mirrored knots, whose pieces are then shared like
one circle's. A vertex within the resolution of a spline cuts it at its
nearest parameter (the distance's derivative's roots, decided exactly at
the rounded parameter). Pieces are ordered by parameter, classified at the
fraction `0.4453125` of their parameter range, and a traced result joins
consecutive pieces of one segment into its exact restriction between their
outer ends (S8b.3's knot insertion), its end poles set to the result's
vertices (within rounding of the restriction's own ends). A stack's walls
on one spline segment are one face on that segment's whole degree-`(p, 1)`
wall (S8b's) where they join, their pcurves lines in (curve parameter,
height); its horizontal edges on a spline are its lifted restrictions.

## Evidence

* **Case protocol.** A Boolean case (`identity_reference.
  encode_boolean_case`) is the object's prism rows (`case`, `op 91`,
  `frame`, `offsets`, `boundary`...), a row `boolean fuse|cut|common 93`
  (the operation and its id), then the tool's rows without their `case`
  row (`op 92`, `frame`, `offsets`, `boundary`...) and `end`: a reader
  splits the block at the `boolean` row and parses each side as an identity
  case, the tool's with the object's `case` row. Natively
  (`native_boolean_case`) the two prisms' construction rows are joined by a
  `boolean OP` row. Every existing generator still passes `--check`.
* **Independent reference.** `boolean_reference.py` (mpmath, 40 digits)
  puts the tool in the object's frame (the origins' offset solved exactly in
  the stored axes and required binary64) and reduces every slab's region to
  the three atoms inside both profiles, the object's only and the tool's
  only. The atoms' areas and first moments come from slicing in `v`: breaks
  at every vertex, every circle's extremes and every meeting of the two
  boundaries (exact discriminants for lines and circles), and in each band
  the elementary intervals' ends integrated in closed form (lines as
  polynomials of canonical Fractions, so coincident edges are one curve;
  circles' branches by `asin` and `(r^2 - t^2)^(3/2)`). Each boundary
  element is cut at its meetings with the other boundary and at the other's
  vertices on it, and its pieces classified at their midpoints (on the
  other boundary with the same or the opposite direction, else inside or
  outside by a ray's parity); a region's perimeter sums the classes whose
  two sides it holds one of. The result: volume and moments from the slabs,
  area as walls (perimeter times height) and caps (the symmetric difference
  of the regions below and above each slab boundary), the centre in world
  coordinates, solids by union-find over the slabs' intervals (positive
  overlaps only: regions meeting at a point or along an edge are separate
  solids, the regularized Boolean), empty for zero volume.
* **Fixtures.** `generate_boolean_fixtures.py --check` writes 45 cases
  (`boolean-cases.txt`; 15 fuses, 17 cuts, 13 commons): overlapping
  rectangles, a bar cutting a square in two, a rectangle and a circle or a
  stadium, lenses, one inside the other (a hole cut, a common), identical
  profiles, disjoint ones, touching along a whole and a partial edge (a
  shared wall), at a vertex, collinear overlapping edges, arcs of one circle,
  tangent circles outside and inside, a vertex on the other's edge, holes (a
  pin through a hole, a hole filled, a common with a hole), a tool missing,
  spanning and partly covering the object's heights, caps at equal and
  touching heights, pockets and a slot (S9a.2), origins offset in the plane
  and along the axis, in the `XY` and tilted `TILT` frames (13 tilted).
  Declared outcomes: 31 `prisms`, 6 `stack` (S9a.2), 5 `empty`, 3
  `degenerate`.
  `boolean-expected.tsv` gives per case `expect KIND STEP` (the declared
  outcome: `prisms`, `stack`, `empty` or `degenerate`, and the computed
  sub-step `S9a.1` or `S9a.2`), `result N volume area cx cy cz` (totals over
  the N solids) or `empty`, and `slab w0 w1 area perimeter` per run of equal
  slab regions. Before writing, the reference is checked against each
  profile's Green's-theorem area and moments and perimeter, `fuse = A + B -
  common` and `cut = A - common` for all three operations, exact Fraction
  clipping (Sutherland-Hodgman by an ear-clipping's triangles, holes
  subtracted) and exact boundary classes on the 21 polygon pairs, the
  closed-form lens (area, centre, arc lengths) on the 8 pairs of circles, and
  7 hand-computed results, all within 1.2e-38 relative; each case's
  declared solid count must equal the reference's.
* **Native.** `occt_boolean_oracle.cpp` builds both prisms as the split
  oracle does and runs `BRepAlgoAPI_Fuse`, `Cut` or `Common`: validity, the
  solids, each solid's volume, area and centre (`BRepGProp`) and its counts
  as built and after `ShapeUpgrade_UnifySameDomain`. Its capture came before
  any kernel code (`fixtures/occt-boolean-preimplementation`): every result
  valid with the reference's solid count, all 45 within the 2e-8 allowance
  (6.1e-15 at worst), no review. Where the decisions differ: OCCT returns two
  valid solids for prisms touching along a vertical edge and one valid solid
  for a hole tangent to the outer circle (the decisions: `Degenerate`), and
  keeps each input's cap and wall images and its cylinders' seams (counts
  compared after unifying). `compare_boolean.py` runs the probe
  (`examples/boolean_probe.rs`) on every case.
* **Spline profiles (S9a.2), before their kernel code.** The reference's
  `SplinePair` (`boolean_reference.py`, chosen by `make_pair` when either
  profile holds S8b's spline segments; S9a's pairs are unchanged) takes each
  spline as its Bezier spans (blossoms, exact). Meetings without resultants:
  with a line or a circle, the exact polynomial of its equation on the span
  split into square-free factors (Yun's algorithm in Fractions: a root's
  multiplicity is exact, even is a tangency), roots isolated by Bernstein
  subdivision and Descartes' rule, refined by bisection to 40 digits and
  polished by `mp.findroot`, checked against `mp.polyroots`; with another
  span, both Bezier forms subdivided while their control boxes meet, then
  Newton on `C1(s) = C2(u)` (residual below 1e-35, transversal); segments
  of one curve share their pieces. Pieces are classified at their midpoints
  (ray parity, rays meeting a span by the same isolation) and must agree at
  their quarter points; the atoms are Green's theorem over the classified
  pieces (spline pieces as exact antiderivatives of polynomials in the span
  parameter), the lengths `mp.quad` of the speed, and S9a's slicing (breaks
  at every span end, `y` extreme and meeting) gives the check and the
  solids. `generate_boolean_fixtures.py --check` writes 46 cases
  (`boolean-spline-cases.txt`, `boolean-spline-expected.tsv`; 14 fuses, 16
  cuts, 16 commons, 18 tilted): a spline crossing lines, arcs, circles and
  another spline (four times, and at both splines' interior knots),
  tangencies of a spline and a line (inside and outside) and a circle,
  identical profiles and one spline shared in the same and the opposite
  direction, a spline hole (cut through, filled, a common), a vertex on a
  spline and a spline's end on an edge, profiles inside and apart, and
  stacks (a step, a pocket, a slot through a spline wall, a plug in a spline
  hole, crossing waves); declared 36 `prisms`, 5 `stack`, 4 `empty`, 1
  `degenerate` (a disc's hole tangent to the dome's apex). A tool holding a
  spline is given in the object's frame coordinates (offset along the axis
  only). Checks before writing: Green over pieces against the slicing
  (8.3e-40), each profile against exact Green's theorem (5.2e-41) and
  Bernstein products (3.9e-41), Gauss-Legendre quadrature of every spline
  piece (3.7e-39), lengths (2.8e-41), the identities (1.6e-40), straight
  splines against their polygon (exact), 8 hand results (dome cut by lines:
  parabolic segments and `asinh` lengths; rectangles, squares, a disc;
  8.6e-41), roots against `polyroots` (4.8e-41), Newton residuals
  (1.8e-40), and the atoms against 16 and 32 chords per span extrapolated
  (2.1e-4 relative, within a quarter of the chords' own difference).
  Natively (`compare_boolean.py --splines`, capture
  `occt-boolean-spline-preimplementation`, `rust_spline_boolean_exists`
  false): every result valid with the reference's solid count, 37 within
  2e-8, 9 reviewed (`occt-boolean-spline-divergences.json`): BRepGProp's
  integration of faces bounded by B-spline edges errs by 2.6e-8 to 8.0e-4
  (the three-span wave S8b's split capture reviewed), while Green's theorem
  over OCCT's own cap edges agrees with the reference within 4.5e-8. OCCT
  keeps tangency points as vertices and edges (a tangency cuts nothing in
  the decisions) and returns one valid solid for the degenerate case.
  The kernel (S9a.2's splines): 45 results inside the reference, the
  degenerate one refused, every count OCCT's after unifying but four
  reviewed tangencies (OCCT keeps the touching point as vertices and edges
  on the dome's wall); `tests/booleans.rs` checks the same and every
  fixture's history.
* **Kernel (S9a).** All 45 cases: 42 results inside the reference with
  the reference's solid count (each solid's volume, area and centre
  enclosed), 5 of them empty and 6 S9a.2 stacks among them, and the 3
  degenerate ones refused; every kernel count equals OCCT's after
  unifying (the kernel's joined collinear lines and cocircular arcs are
  OCCT's unified faces and edges; its seamless circles counted as OCCT's
  seams by `Topology::occt_counts`), no review. `tests/booleans.rs` checks
  the same against the reference, every fixture's history (independent
  check, every input entity covered, repeated exactly) and hand cases of
  every class (overlaps, a lens, a hole, identical, disjoint, inside, a
  shared wall, tangent cylinders inside and outside, arcs of one circle, a
  vertex on an edge, a holed box, frames with equal axes and offset
  origins, a result taken as an input again, inputs sharing ids refused;
  S9a.2: a tower and a pocket with their counts, classification and rigid
  motion, a plug filling a hole over part of its height, a cavity, a tool
  through a round wall).
* **S9b evidence (polyhedra in any position).** `polyhedral_reference.py`
  decides each prism on its exact model (the stored axes as rationals),
  cuts each profile into trapezoids so each prism is a union of convex
  cells, and clips every pair of cells by exact half-spaces in Fractions:
  the common's volume and moments exactly, the fuse and cut by inclusion
  and exclusion (checked against the result's own convex cells), areas by
  splitting each boundary face by the other prism's cells' planes and
  classifying each piece on both sides by an infinitesimal push, solids by
  convex cells sharing positive area. `generate_polyhedral_fixtures.py
  --check` writes `boolean-polyhedra-cases.txt` and
  `boolean-polyhedra-expected.tsv` (45 cases: a turned box, a tilted bar
  cutting a box in two, coplanar caps and walls of either orientation, an
  edge and a vertex on a face, a tilted corner, a box inside another, apart,
  an L profile, a bar through a hole and across it, both inputs turned),
  after checking the reference against S9a's slicing on its 21 polygon
  cases (within 2.3e-17: the stored axes' departure from orthonormal),
  `area(A u B) + area(A n B) = area(A) + area(B)` (3e-41) and the closed
  forms of axis-aligned boxes. `compare_polyhedral.py` reproduces the
  `BRepAlgoAPI` capture `occt-boolean-polyhedra-preimplementation`, taken
  before S9b's kernel module: every result valid with the reference's solid
  count, all 45 within 7.9e-16, no review. The kernel (S9b.1): 41 results
  inside the reference, the 4 degenerate ones refused (two solids touching
  along an edge or at a vertex, two walls of a tool's edge passing within
  rounding of the object's edge: a neck thinner than the resolution, the
  last declared after the kernel met it), every count OCCT's after
  unifying but two reviewed (OCCT keeps a tool's touching edge or vertex as
  an imprint on the object's face; the kernel's regularized cut is the
  object). `tests/polyhedral_booleans.rs` checks the same, every fixture's
  history (independent check, every input entity covered, repeated
  exactly) and hand cases (a turned box's quarter, its rigid motion and
  classification, a cavity, a cut in two, touching solids refused).
* **S9c.1 evidence (arcs in any position), before its kernel code.**
  `curved_boolean_reference.py` (mpmath, 40 digits) takes each prism on its
  exact model (the stored axes as rationals) and slices both solids by the
  planes parallel to both axes (`d = n_A x n_B`, or `n x e` for parallel
  axes): each section is a union of parallelograms, one per chord of the
  profile (its ends `L(s) + k sqrt(Q(s))`, exact), so a slice's common, cut
  and fuse are convex clippings whose areas and moments are Green's theorem
  over segments (the decisions' ellipse arcs appear only in slices of
  another direction; here they bound the planar faces' regions of the area
  part). Breakpoints are the real roots of exact polynomials (three lines
  of the slice concurrent, two parallel ones coinciding, a trace through a
  vertex or tangent to a circle: surds squared away, spurious roots
  filtered by the meeting lying on both sections' boundaries), and each
  interval is integrated by Gauss-Legendre after `s = a + (b - a)(1 - cos
  t)/2` (end-point square roots analytic), refined to 1e-33 of the case's
  size to the fourth. Areas: every input face is swept by lines of its own
  parameters (`v` on caps, the height on flat and cylindrical walls), each
  line cut at its crossings of the other solid's boundary and its pieces
  classified at their midpoints (inside, outside, or on a coplanar or
  coincident face of the same or opposite orientation, decided in
  rationals), the class lengths integrated with the surface's own element
  between breakpoints found the same way (in `tan(theta / 2)` on
  cylinders). Solids: the slices' arrangements of lines, faces joined
  across one line within an interval and by overlapping limits across a
  breakpoint. `generate_curved_boolean_fixtures.py --check` writes
  `boolean-curved-cases.txt`, `boolean-curved-expected.tsv` and
  `boolean-curved-frames.tsv` (the stored axes' bits) with 44 cases (14
  fuses, 13 cuts, 17 commons): a tilted cylinder through a box, a box
  corner in a cylinder (exact and leaning), Steinmetz solids (perpendicular,
  oblique, both axes tilted), parallel cylinders with coplanar caps, coaxial
  cylinders, a tilted pin through a square hole, a coaxial pin filling a
  round hole (coincident cylinders of opposite orientations), a tilted
  cylinder across a hole's wall, a stadium in a turned and in a tilted
  frame, quarter cylinders, and 7 `degenerate` with reasons (a plane
  tangent along a generatrix, cylinders tangent outside and inside, the
  Steinmetz cut touching itself at two points); 35 `solid`, 2 `empty`, 16
  in exact frames only. Frames are those whose stored axes the kernel gives
  bit for bit: not `ROT`, whose `x` differs from `stored_axes` in its last
  bit on macOS arm64 (the platform `hypot`); `R125` (`x` along `(12, 5,
  0)`) turns instead. Checks before writing: closed forms (9.2e-41 in exact
  frames, 1.9e-16 in turned ones: the stored axes' departure from
  orthonormal), `fuse = A + B - common` and `cut = A - common` with each
  operation sliced apart (2e-41), the area identity and every face's
  classes against its closed-form area (2e-40), a second slicing direction
  for parallel axes (1.4e-41), Monte-Carlo volumes and centres (2.7
  standard errors at worst), and S9a's and S9b's references on their 90
  fixtures (every count; S9b's 25 printed digits exactly, S9a's within
  2.4e-17, its frame coordinates taking the axes as orthonormal). A scan
  for near coincidences moved two fixtures off coincidences of the slicing.
  `compare_curved_boolean.py` (`compare_boolean.py` through its `make_set`
  hook, which also repaired `compare_polyhedral.py`: since the spline set it
  had checked S9a's capture) reproduces `occt-boolean-curved-preimplementation`
  (`rust_curved_boolean_exists` false): every result valid with the
  reference's solid count, the degenerate ones included, all 44 within
  8.6e-9 (`across_hole_common`, BRepGProp on faces bounded by ellipses), no
  review; four solids' counts change when unified (coplanar caps merged, a
  tangency's imprint); the kernel's probe `unsupported` on all 44.
* **Kernel (S9c.1).** `tests/curved_booleans.rs`: the 44 fixtures as the
  reference (31 results and empties inside its measures, the 7 declared
  degenerate refused, the 6 of cylinders not circular in a common measure
  `OutOfDomain`, S9c.2's), every history checked, results deterministic and
  moved rigidly with their ids, volumes and boundaries.
  `compare_curved_boolean.py`: 42 matches, 2 reviewed (the Steinmetz fuse's
  and common's counts, OCCT's seam edges on their faces;
  `occt-boolean-curved-divergences.json`), no failure. The `boolean`
  target's arcs in turned, leaning and tilted frames now reach S9c.1; its
  corpus (1,080 inputs), the spline variants and the regressions replay
  clean with the history check.
* **S9c.2 evidence (two cylinders meeting in quartics), before its kernel
  code.** `generate_procedural_boolean_fixtures.py --check` writes
  `boolean-procedural-cases.txt`, `boolean-procedural-expected.tsv` and
  `boolean-procedural-frames.tsv` with 28 cases (10 fuses, 9 cuts, 9
  commons; 22 solid, 6 degenerate) from S9c.1's reference, frames, checks
  and near-coincidence scan, imported unchanged (the slicing by planes
  parallel to both axes takes cylinders of any radii and offset: each
  section is still a union of parallelograms). S9c.2a, 22 in exact frames
  (a thick cylinder along `z` in `XY`, `TURN` or `DOWN` against a pipe along
  `x` in `SIDE`, offset along `y`): a pipe through with offset and with
  crossing axes (two rings; the pipe the object in the second, its cut two
  stubs), a partial bite and equal radii offset (one loop), a pipe ending
  inside (its cap inside), a box's round hole crossed by a pipe (rings on
  the hole's wall, a common of two pieces), and internal and external
  tangency declared `degenerate` for every operation (S9c.1's evidence
  amendment (c)). S9c.2b, 6 in turned frames: oblique crossing axes
  (`TILT`), skew unequal cylinders (`LEAN`), a pin in `R125` parallel to a
  box's round hole cutting its wall. Checks before writing, besides S9c.1's
  per pair: the common of perpendicular cylinders by one tanh-sinh
  quadrature in `eta` between its kinks (`int len([xA - wA, xA + wA] & [bLo,
  bHi]) len([zB - wB, zB + wB] & [aLo, aHi]) d eta`, `4 int wA wB d eta` for
  cylinders crossing whole), its moments, and its area from the walls' own
  angles and the caps' parts inside the other, so all three operations'
  volumes, areas and centres (a round hole as the box less its cylinder);
  Legendre's `8 rA / 3 ((rA^2 + rB^2) E(k) - (rA^2 - rB^2) K(k))`, `k = rB /
  rA` (mpmath's `ellipe`, `ellipk` of `m = k^2`; at equal radii `16 r^3 /
  3` exactly, `K(1)`'s zero factor dropped), against the quadrature
  (2.2e-41) and the reference, the common's centre at the axes' crossing;
  in turned frames the perpendicular common over `sin phi` and parallel
  cylinders' lenses. Closed forms 7.9e-41 in exact frames and 5.0e-16 in
  turned ones, inclusion and exclusion 1.2e-41, the area identity 1.3e-40,
  every face's classes 9.9e-41, a second slicing direction 6.4e-42, Monte
  Carlo 3.2 standard errors; no near coincidence, no fixture moved, and the
  curved fixtures regenerate byte for byte (their `--check` passes).
  `test_procedural_boolean_reference.py` checks the closed forms alone.
  `compare_procedural_boolean.py` (`compare_boolean.make_set`) reproduces
  `occt-boolean-procedural-preimplementation`
  (`rust_procedural_boolean_exists` false; the kernel's probe `unsupported`
  on all 28, `OutOfDomain(... (S9c.2))`): every result valid with the
  reference's solid count, the degenerate ones included (the internal
  tangency's node unmarked, the external one's fuse two solids touching at
  a vertex); 5 matches without a quartic edge (parallel cylinders, the
  external tangency) and 23 reviewed (`occt-boolean-procedural-
  divergences.json`): BRepGProp's default integration on faces bounded by
  the intersection's approximated curves misses by up to 3.0e-4 (`skew`,
  7.8e-5 in exact frames), the same results measured adaptively (`Eps =
  1e-10`) within 2.7e-9; six solids' counts change when unified.
* **Kernel (S9c.2a).** `tests/procedural_booleans.rs`: the 22 exact-frame
  fixtures as the reference (16 results inside its measures, the 6
  declared tangencies refused), the 6 turned-frame ones `OutOfDomain`
  (S9c.2b's), every history checked, results deterministic and moved
  rigidly. `compare_procedural_boolean.py`: 5 matches, 23 reviewed (the 23
  measures, now with the 11 kernel results whose counts differ: the loops'
  switch vertices, OCCT's own splits of its approximated rings), no failure;
  S9a's, S9b's and S9c.1's comparisons unchanged. The DRAW survey's two
  S9c.1 failures are fixed: a hole no piece holds (a sliver whose binary64
  image turned the wrong way) is `Degenerate`, and a band with holes
  bounded by ellipses' projections now validates.
* **Kernel (S9c.2b.1).** `tests/turned_booleans.rs`: the 15 fixtures as the
  reference (13 results inside its measures, the 2 near nodes refused),
  every history checked, results deterministic and moved rigidly; the
  procedural set's `oblique` and `skew` now inside the reference too (its
  `parallel_hole` S9c.2b.2's), and S9c.1's `steinmetz_oblique` and
  `steinmetz_tilted` declared `degenerate` (their stored models touch at
  the ends of their common extent). `compare_turned_boolean.py`: 2
  matches, 13 reviewed (the measures; the loops' counts with their switch
  vertices added), no failure; `compare_procedural_boolean.py` 5 and 23,
  `compare_curved_boolean.py` 42 and 2, no failure.
* **S9c.2b.1 evidence (cylinders in turned frames), before its kernel
  code.** `generate_turned_boolean_fixtures.py --check` writes
  `boolean-turned-cases.txt`, `boolean-turned-expected.tsv` and
  `boolean-turned-frames.tsv` with 15 cases (5 fuses, 5 cuts, 5 commons; 13
  solid, 2 degenerate), every pair two cylinders with axes not parallel and
  at least one frame turned (affine models on the stored axes, elliptic in
  the world), from S9c.1's reference, frames and per-pair checks and
  S9c.2's closed-form machinery, imported unchanged: a partial bite, one
  loop (`bite_tiltx`, `XY` against `TILTX`, all three operations); a thin
  pipe through a thick cylinder as the object, two rings (`rings_pipe`,
  `TILT2` against `XY`; its cut two stubs); a pipe ending inside, its cap's
  disc inside (`blind_tilt`, `TILT` against `SIDE`, perpendicular); a box's
  round hole crossed by a tilted pipe, two rings on the hole's wall
  (`hole_tiltx`, its common two pieces); unequal skew cylinders offset, one
  loop (`skew_lean`, `LEAN` against `TILTX`); and equal radii with meeting
  axes declared `degenerate` (`node_lean`, `LEAN` against `TILTX`;
  `node_r125`, `R125` against `SIDE`: the stored models' intervals' ends
  2.4e-17 and 8.3e-17 apart, two rings about 1e-8 apart at the nodes).
  Checks before writing, besides S9c.1's per pair: closed forms in the
  ideal frames (the common of cylinders crossing whole as the perpendicular
  one over `sin phi`, centred on the axes' common perpendicular; the
  perpendicular pair by the quadrature with its cap, areas included; the
  hole as the box less its cylinder and the pipe between the box's walls
  by Cavalieri; the nodes by Legendre's form at `k = 1`) within 4.2e-16;
  the curve's class from the ideal and the models' intervals (as declared,
  ends at least 0.175 apart, the nodes' within 1e-15 but apart); every
  cap's circle against the other input's cylinders (crossings of the model
  found by sign changes along the circle and bisection, required outside
  the face along its axis; none crosses, the least clearance 0.45 of the
  radius squared), so the quartic meets no cap's circle and every vertex is
  a quadratic surd. Inclusion and exclusion 1.2e-41, the area identity
  6.0e-41, every face's classes 6.9e-41, Monte Carlo 2.4 standard errors;
  the scan finds only the nodes' breakpoints (about 1e-17 apart).
  `test_turned_boolean_reference.py` checks the closed forms, the node
  models and the cap-circle test alone. Where one frame's stored `x` or `y`
  lies along the axes' common perpendicular exactly (`XY` against `TILT`
  or `LEAN`, `TILT` against `TILTX`: S9c.1's `steinmetz_oblique` and
  `steinmetz_tilted`) the models' intervals are equal exactly, a double
  tangency whose section is two conics, not two loops.
  `compare_turned_boolean.py` (`compare_boolean.make_set`) reproduces
  `occt-boolean-turned-preimplementation` (`rust_turned_boolean_exists`
  false; the kernel's probe `unsupported` on all 15, `OutOfDomain(...
  (S9c.2))`): every result valid with the reference's solid count; the two
  nodes match within 1.5e-10 (OCCT's section two ellipses within its
  tolerance) and 13 are reviewed (`occt-boolean-turned-divergences.json`):
  BRepGProp's default integration on faces bounded by approximated
  intersection curves misses by up to 3.3e-4 (`skew_lean`), the same
  results measured adaptively within 9.1e-10; four solids' counts change
  when unified.
* **S9c.2b.2 evidence (sections crossing caps' circles), before its kernel
  code.** `generate_capped_boolean_fixtures.py --check` writes
  `boolean-capped-cases.txt`, `boolean-capped-expected.tsv` and
  `boolean-capped-frames.tsv` with 18 cases (5 fuses, 6 cuts, 7 commons;
  all solid), every pair two cylinders whose quartic section crosses a
  cap's circle within both faces (the vertex there a nested surd in exact
  frames, a root of the circle's quartic in turned ones), from S9c.1's
  reference, frames and per-pair checks, S9c.2's `Perpendicular` and
  S9c.2b.1's classes, imported unchanged. Exact frames (`XY` against
  `SIDE`): a pipe across the top cap's rim, its axis 0.1 below the cap's
  plane, two rings each crossing the cap's circle twice (`rim_pipe`, all
  three operations); a bite at the rim, one loop crossing it twice
  (`rim_bite`); a pipe ending partway through the wall, its end cap's
  circle crossing it (`blind_wall`); a box's round hole crossed at its top
  rim by a pipe along the top face (`hole_rim`, its common two pieces).
  Turned frames, oblique: a `TILTX` pipe across the rim (`rim_tiltx`), a
  `LEAN` pipe entering through the rim and ending inside (`enter_lean`, a
  blind hole drilled at the rim), a `TILT` pipe ending partway through the
  wall (`blind_tilt`), a `TILTX` pipe crossing a round hole's top rim and
  leaving through the box's bottom (`hole_rim_tiltx`, two pieces). Checks
  before writing, besides S9c.1's per pair: the exact-frame pairs' closed
  forms (the common of perpendicular cylinders clipped by both inputs' caps
  by one quadrature in `eta`, moments and areas included; the hole as the
  pipe's circular segment below the top face less its common with the
  hole's cylinder, the top face's strip less the hole's disc between two
  chords) within 6.4e-41; the classes from the ideal and the models'
  intervals (ends at least 0.6 apart); every cap's circle against the
  other input's cylinders, the opposite of S9c.2b.1's check: crossings of
  the model by sign changes along the circle and bisection, each at least
  0.1 from the other face's ends (2.0 achieved) at a slope of at least
  0.05 (0.25), at least one per pair on the other face (20 of 22), a
  circle without crossings clear by 1e-3 of the radius squared (0.66).
  Inclusion and exclusion 1.5e-41, the area identity 5.9e-41, every face's
  classes 3.3e-41, Monte Carlo 2.9 standard errors; no near coincidence
  (perpendicular pairs in stored turned frames were dropped: a cap plane
  parallel to the other's generatrices within rounding puts two slicing
  breakpoints about 1e-17 apart). `test_capped_boolean_reference.py`
  checks the closed forms and the cap-circle test alone.
  `compare_capped_boolean.py` (`compare_boolean.make_set`) reproduces
  `occt-boolean-capped-preimplementation` (`rust_capped_boolean_exists`
  false; the kernel's probe `unsupported` on all 18, `OutOfDomain(...
  crossing a cap's circle (S9c.2b))` in exact frames and `(... in turned
  frames (S9c.2b.2))` in turned ones): every result valid with the
  reference's solid count; 2 match and 16 are reviewed
  (`occt-boolean-capped-divergences.json`): BRepGProp's default
  integration on faces bounded by approximated intersection curves misses
  by up to 1.4e-5 in exact frames and 1.9e-5 in turned ones, the same
  results measured adaptively within 1.0e-9 and 1.9e-8 (a small result's
  approximated section, unchanged at `Eps = 1e-12`); eight solids' counts
  change when unified.
* **Kernel (S9c.2b.2).** `tests/capped_booleans.rs`: the 18 fixtures as the
  reference, every history checked, results deterministic and moved
  rigidly; the procedural set's `parallel_hole` and S9c.1's
  `parallel_cylinders` now inside their references too, no case of any
  cylinder set left `OutOfDomain`. `compare_capped_boolean.py`: 18
  reviewed (the measures; ten results' counts with the same faces as OCCT's
  unified result, edges and vertices where each splits its section
  curves), no failure; `compare_procedural_boolean.py` 4 matches and 24
  reviewed (`parallel_hole_fuse`'s counts), `compare_turned_boolean.py`,
  `compare_curved_boolean.py`, `compare_polyhedral.py` and
  `compare_boolean.py` without failure. The bridge self-test's three
  cylinder gaps (a turned cylinder's quartic across its cap, parallel
  turned cylinders, a pipe through a cap's rim) are decided now.
* **S9d.1 evidence (spheres against polyhedral prisms), before its kernel
  code.** `sphere_boolean_reference.py` (mpmath, 40 digits) takes the
  sphere, cap or zone on its exact model (`|X - o|^2 <= R^2`, the end
  planes `w = R sin(latitude)` in the frame's affine coordinates, the
  heights as stored) and the prism on S9c.1's, a nonconvex profile
  ear-clipped into convex pieces; it slices both by the planes normal to
  the zone's axis `m = x * y` (any rational direction accepted): each
  slice's disc, cut by the zone's lines, against the pieces' convex
  sections, the four operations apart by the boundaries cut and classified
  (Green's theorem over segments and arcs of the one circle in closed
  form), the common again by clipping each piece; breakpoints the roots of
  exact quadratics (a slice tangent to the sphere, through a vertex, an
  edge's meeting with the sphere or a zone plane, tangent to a face's or a
  zone plane's circle, through a zone plane's and a face's meeting with
  the sphere), Gauss-Legendre between them refined to 1e-33 of the case's
  size to the fourth. Areas: the sphere's face by Archimedes (`R / |d|`
  times the integral of its section's angle inside the pieces), every
  planar face and end disc in closed form in its own plane (inside,
  outside, or on a coplanar face of the same or the opposite orientation).
  Solids by convexity: pieces meeting the sphere (an exact squared distance
  from the centre to a polyhedron against `R^2`), joined across shared
  faces; `S - K` by the overlaps of the sphere's parts beyond `K`'s faces,
  `K - S` by the runs of `K`'s faces' boundaries outside the sphere.
  `generate_sphere_boolean_fixtures.py --check` writes
  `boolean-sphere-cases.txt`, `boolean-sphere-expected.tsv` and
  `boolean-sphere-frames.tsv` with 30 cases (4 fuses, 12 cuts, 14 commons;
  25 solid, 1 empty, 4 degenerate; 25 in exact frames): a sphere through a
  box's face, a box's corner in a sphere (exact, a `TILT` box, a `LEAN` box
  against a `TILTX` sphere), a box's edge through a sphere, a square post
  through a zone's flat end, a box through a zone's flat end and band (the
  zone along `x`), a hemisphere crossed through its disc and dome, a dome
  on a box's face (coplanar discs, opposite; the common empty) and in a box
  on its bottom face (the same orientation), a sphere in a box and a box in
  a sphere (cavities), a `TILT` bar through a sphere (the bar's cut two
  solids), a slab cutting a sphere in two, an L prism severed at its
  corner (two solids), an octant, and `degenerate` a face tangent to the
  sphere, a box's eight vertices on it, an edge tangent to it. Caps and
  zones only in exact frames centred at the origin (every reading of their
  end planes the same plane). Checks before writing: closed forms of every
  pair (a rectangle's part of each section integrated along the box's
  axis, both surfaces' parts inside the other) within 9.3e-40 in exact
  frames and 2.2e-16 in turned ones, both inputs' slicings against their
  closed forms and inclusion and exclusion 9.2e-41, the common two ways
  1.4e-42, the area identity 3.3e-40, every face's classes 1.8e-40, a
  second slicing direction 1.6e-40, Monte Carlo 2.7 standard errors, no
  near coincidence. `test_sphere_boolean_reference.py` checks the closed
  forms (a cap `pi h^2 (3R - h)/3` and `2 pi R h`, an octant `pi R^3 / 6`
  and its centre `3R/8`, a zone, a sphere and a box inside) and the
  reference on the cap. The Boolean protocol takes a sphere on either side
  (its identity rows `frame` and `sphere R LOW HIGH`, which the kernel's
  test support already reads); `native_case` gives it one `sphere` row,
  built by `BRepPrimAPI_MakeSphere` in `occt_boolean_oracle.cpp`.
  `compare_sphere_boolean.py` (`compare_boolean.make_set`) reproduces
  `occt-boolean-sphere-preimplementation` (`rust_sphere_boolean_exists`
  false; the kernel's probe `unsupported` on all 30, `OutOfDomain("a
  Boolean of a solid with curved faces or edges in any position (S9c)")`):
  every result valid with the reference's solid count, all 30 within 3.3e-9
  (BRepGProp on spherical faces bounded by circles off their parallels;
  sections along parallels within 4e-15), no review; three solids' counts
  change when unified.
* **Fuzzing.** The `boolean` target (`FUZZING.md`): the split target's line
  and arc profiles, the tool offset exactly in the axis-aligned frame or
  sharing the tilted one's origin, heights equal, spanning, overlapping,
  disjoint, inside (pockets, cavities) or on top (touching stacks); fuse,
  cut and common with the volume identities, rigid motions of the results
  and every result vertex classified on its boundary.

## DRAW commands

The DRAW adapter (`examples/draw_worker.rs`, `UPSTREAM_TESTS.md`) runs
OCCT's Boolean commands of one object and one tool through `Solid::fuse`,
`cut` and `common`: `bfuse`, `bcut`, `bcommon`, `btuc`, `bop` with
`bopfuse`, `bopcut`, `boptuc` and `bopcommon`, and `bbop`/`bapibop` 0 to 3
on the General Fuse arguments, returning OCCT's compound of the result's
solids and reporting every refusal (`OutOfDomain`, `Degenerate`, ...)
unsupported. S9a.2's stacks are Boolean results like the prisms:
`checkshape`, `nbshapes`, `vprops`, `sprops` and `lprops` read their
topology (a closed cavity is a solid of two shells), and `unifysamedom`
returns them unchanged, since the kernel builds them unified; native
DRAW's unified counts agree on every stack checked. S9b.1's polyhedra are
handled alike: their measures and counts are read from their topology,
`unifysamedom` returns them unchanged (coplanar fragments are joined into
maximal faces and collinear edges joined as they are built), and native
DRAW's unified counts agree on every polyhedron checked. A stack or a
polyhedron given to another Boolean is taken as it is (S9b.2, on its
stored geometry; the prism argument built again when they share ids).
`ttranslate` and `trotate` move a prism by
the kernel's rigid motion, a prism in the moved frame; `tcopy`, like
`copy`, gives the same shape. A `trotate` by whole quarter turns about a
coordinate axis turns a prism's frame exactly (a signed permutation of
coordinates, the origin by DRAW's location arithmetic): the kernel's
rotation rounds the cosine of a quarter turn to 6.1e-17, as OCCT's
`gp_Trsf` does, and OCCT's tolerances absorb it, while the kernel's exact
decisions would find a wall turned onto another's plane tilted off it (a
box and its quarter-turned copy fused into an L would keep a crease: 9
unified faces where OCCT has 8). The derived cases `boolean_prisms`,
`boolean_stacks` (a step, a pocket, a box cut in two by a slab, a closed
cavity, a tool through a round wall) and `boolean_polyhedra` (quarter
turns, a bar turned 45 degrees through a box, a tilted bar cutting a box
in two, a turned box inside another) and 822 cases of upstream's
`boolean` group (86 of them stacks, 300 polyhedra of two boxes, one
turned, 35 pockets cut one after another, a stack or a polyhedron given
to the next `bcut`, 65 of S9c.1's prisms with arcs in any position: a
cylinder and a turned box, a cylinder turned about its axis, equal
cylinders crossed at right angles, and 16 of S9c.2a's: a cylinder of
radius 0.5 through one of radius 1 at right angles, two quartic rings)
evaluate on both backends; none fails on the kernel, and S9c.2b.1 adds
none. Of the upstream cases in frames with different axes the rest are
refused: parallel cylinders in turned frames S9c.2b.2 takes, solids
other than prisms, S9b.1's, S9c.1's, S9c.2a's and S9c.2b.1's
`Degenerate` (a turned box's corner on another's wall, edge or corner, or
on a cylinder, within rounding; a wall tangent to a cylinder; two
cylinders touching at a point; equal cylinders in a turned frame whose
axes meet, their section within the resolution of a node), and an arc
ending off its circle; of the stacks given to another Boolean, those
with cylindrical walls (S9c; `UPSTREAM_TESTS.md`).
