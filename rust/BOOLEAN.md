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
parallel cylinders in turned frames (`curved/algebraic.rs`), and S9d.1's
spheres against polyhedral prisms (`curved/sphere.rs`), and S9d.2a's
spheres against prisms with arcs and two spheres (`curved/spheres.rs`,
S9d.2b's loops in exact frames among them, S9d.2c's turned caps' circles
and loops in turned frames in `curved/spheres_turned.rs`), and S9d.3a's
cones and frusta against polyhedral prisms (`curved/cone.rs`), and S9d.3b.1's
cones against cylinders, spheres and cones meeting in rings or on a plane
(`curved/cones.rs`), S9d.3b.2's cones and spheres in loops, S9d.4a's
whole tori against polyhedral prisms (`curved/torus.rs`), S9d.4b.1's
torus v-segments and wedges against them (`curved/torus_segment.rs`), and
S9d.4b.2a's whole tori against prisms with arcs, spheres and cones and
S9d.4b.2b's two whole tori (`curved/torus_curved.rs`); parts of a torus
against curved faces are not.

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
  end, whose order its enclosures leave undecided, or a result's cavity
  whose containment the validator's rays leave undecided (a ray meeting a
  sphere or torus face with loops).
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

### Spheres against polyhedral prisms (S9d.1)

A sphere, a cap or a zone (S3's `sphere_*`) against a prism takes S9c's
arrangement with a second kind of model (`curved/sphere.rs`): the set
`|p - c|^2 <= r^2` of its stored centre and radius, a zone's end planes
through `c + h n` normal to its stored axis at the stored heights, their
discs as faces, and its wall split into two hemispheres by a plane through
an axis and a rational point of the equator (a cap's or zone's axis is its
frame's; a whole sphere's a generic rational direction that moves with the
seam, so its poles leave special points; each input its own seam, tried
again at another when a meeting falls on it). The split's great circle and
the rims are the model's edges, cut at the rims' points on the split (quadratic
surds) and at the poles (virtual vertices). Circles of a surd radius
(sections, rims, the split) are `Crv::Circle`: points `c + dx x + dy y`
over a rational orthogonal basis of their plane, placed by `(dx, dy)` (the
arrangement's angle orders are sign tests, unchanged), a rational direction
scaled onto the circle giving a point of one quadratic field. A plane's
section of the sphere is such a circle (centre the centre's projection,
radius squared `r^2 - d^2`); a line meets the sphere in a quadratic, a
circle meets a plane on its line of the circle's plane. Membership is exact
at first order (a push along the sphere's tangent plane keeps to it), a
hemisphere's region the split plane's side and the ends' heights, a disc's
its rim's inside. Each hemisphere is traced in its stereographic projection
from the opposite pole onto the split plane (conformal: the rim is not
compressed as an orthographic projection compresses it). A circle meeting
no seam and lying inside both faces is a ring, given one vertex at a
rational point and joined back into a ring edge. In the result a circle is
`Curve3::Circle` or `CircularArc` on its basis's frame, its pcurve on a
sphere a line where it is a meridian (its `u` taken inside it, so an end at
a pole needs none) or a parallel, an exact projection otherwise; a sphere
face whose loops wind once in all holds its pole as a vertex loop (S3's
caps). A section through the stored sphere's pole, or within the
resolution of it, has a vertex there, kept in the result, its meridian
arcs' pcurves lines on either side; a section through a pole off the
axis's planes is `OutOfDomain`. The validator decides a ray against a whole sphere (a cavity in a
sphere); a result's bounds hold its spheres' boxes. A sphere against a
cylinder or another sphere is S9d.2's.

### Spheres against cylinders and spheres (S9d.2a)

Two spheres meet in their radical plane's section, a sphere's own circle
meets another sphere on that plane (`curved/spheres.rs`). A cylinder and a
sphere meet where the cylinder's rulings meet the sphere: the ruling's
quadratic has the discriminant of S9c.2b.1 with the sphere as the other
quadric (three world rows instead of a cylinder's two), positive all round
(two rings over the cylinder's angle, `Curve3::Meet` with `other_sphere`:
a pipe through a ball, centred on the axis or not), negative all round
(apart), or changing sign (a loop: S9d.2b's, `OutOfDomain`); only a real
repeated root is a tangency. A prism's cap circle meets a sphere at the
roots of its quartic (`algebraic.rs`'s circle points against any quadric).
A sphere's own circle (a rim, the split) meets a cylinder where `F0 + s F1
= 0`, `s` its radius over its basis's length, a basis of equal axes: a
quartic when `s` is rational, else `F0^2 - s^2 F1^2` with its roots of the
right signs kept, the points surds over `Q(alpha)`. A whole sphere's split
therefore uses two rows of a rational rotation (world axes of one rational
length, picked by the seam); a zone's in an exact frame has equal axes too;
a zone's circle of unequal axes (a turned cap's) against a cylinder is
S9d.2c's (below), against a cone S9d.3c's (`OutOfDomain`). A meeting
on a hemisphere is a seam's only on the split (a prism's vertex on the
sphere elsewhere is `Degenerate` at once).

A loop of a sphere and a cylinder (S9d.2b) in an exact frame: at height
`w` the cylinder's circle meets the sphere where `alpha cos u + beta sin u
= g(w)`, `g` quadratic, so near a ruling's tangency (a root of the
cylinder's discriminant, the curve vertical) the curve is a graph over the
height, `u = phi +- acos(g(w) / rho)`, whose points at rational heights
are quadratic surds (`trig`); near a circle's tangency (a root of `rho^2 -
g(w)^2`, the curve horizontal) it is S9c.2b.1's graph over the angle. The
two kinds of turning point are ordered along each loop (binary64 views of
isolators narrowed 160 bisections) and a rational switch of the cylinder's
angle placed between adjacent ones of different kinds; each piece over the
height is verified exactly (no root of `rho^2 - g^2` between its ends'
heights, Sturm counts at surds), each over the angle as S9c.2b.1's. A
piece over the height is `Curve3::Rise` (`TOPOLOGY_MODEL.md`), its range in
the stored cylinder's heights.

### Turned caps' circles and loops in turned frames (S9d.2c)

Where a frame is turned (`curved/spheres_turned.rs`). A cap's or zone's
own circles on a turned frame's stored axis lie on bases of unequal
lengths, `c + X x + Y y` with `xx X^2 + 2 xy X Y + yy Y^2 = r2`: against a
cylinder, in coordinates turned by a rational rotation (the identity, then
`(3, 4, 5)`, `(5, 12, 13)`, `(8, 15, 17)`), the circle's and the
cylinder's functions are quadratics in `t` whose resultant is a quartic in
`s`; each simple real root's `t` is the two quadratics' common root, so the
meeting is a point of `Q(alpha)` (degree four), placed on the circle by
its `(X, Y)`. A repeated root or a vanishing denominator tries the next
rotation, a tangency failing all four; a crossing within the resolution of
a tangency (the meetings with the cylinder's function offset by `+-2 rho
res` counted differently) is `Degenerate`, since a turned frame's rounding
leaves no exact one. On a turned cylinder the model's circle at height `w`
is an ellipse in the world and the sphere's function on it, `F(u, w)`, is
of degree two in `(cos u, sin u)`: a loop's pieces over the angle stay
S9d.2b's, and a piece over the height is the one root of `F(., w)` on its
branch's half-turn (the side of `(alpha, beta)` in the cylinder's local
coordinates), its point at a rational height a root of a quartic. Each is
verified exactly: no root on the half-turn's boundary over its heights (a
quartic `E(w)` in the height), no double root (the discriminant of `F`'s
quartic in the half-angle tangent, degree twelve in `w`, whose real roots
are the height graph's turning points and order the loop's events), one
root on the branch at a rational height inside. Heights, branches and the
cylinder's membership are read in its local coordinates (the inverse
frame's rows, its axes in an exact frame, where every result is
unchanged). `Curve3::Rise` keeps its closed form, the binary64 reading on
the stored frame. A cone's loops in a turned frame and a turned cap's
circle against a cone are S9d.3c's (`OutOfDomain`).

### Cones against polyhedral prisms (S9d.3a)

A cone or frustum (S3's radii `b` and `t` over `0..h` on its stored frame)
is `u^2 + v^2 <= (b + k w)^2`, `0 <= w <= h`, `k = (t - b) / h`, in the
frame's coordinates as rationals (affine where the stored axes are not
orthonormal), its membership pushed at first order as a sphere's
(`curved/cone.rs`). Its wall is a graph over the plane of `(u, v)` (`k` is
never zero), so it is one face traced in that projection with no seam: a
rim is one closed edge with a vertex at a rational point (a meeting there
is tried again at another seam), the apex a point inside the face, and a
result's wall whose loops wind once round the axis closes at the apex, a
vertex loop continuing the input's. A line meets the wall at a quadratic's
roots (a double root is a tangency where the wall is, nothing beyond its
ends: a line through a frustum's virtual apex outside it). A plane
`alpha u + beta v + mu w + kappa = 0` meets it in a graph over the cone's
angle: the ruling at `(cos, sin)` meets the plane at `rho = G / (mu + k
(alpha cos + beta sin))`, `G = b mu - k kappa`, so a rational angle gives a
rational point and `rho > 0` picks the solid's nappe; the directions
where the denominator vanishes (the plane parallel to a ruling) bound a
hyperbola's branch or a parabola's, whose ends at infinity lie beyond the
end planes, an ellipse is closed, a plane normal to the axis a circle
about it; `G = 0` (a plane through the apex, virtual for a frustum, or
tangent along a ruling) is `Degenerate`. Each section edge is rounded once
by S8d.2's conic (`solid/split/conic.rs`'s `cone_conic`, the plane in the
cone's frame): `EllipseArc`, `HyperbolaArc`, `ParabolaArc`, or a
`CircularArc` about the axis, its parameter running as the chain runs;
its pcurves on the cone are `Projection`s, a circle about the axis a line
of constant `v`. A cone against a prism with arcs, a sphere or a cone is
S9d.3b's (`OutOfDomain`). Two parallel faces of the inputs apart by no
more than the resolution (not on one plane exactly) are `Degenerate`: the
sliver between them is thinner than any loop's binary64 image can nest.

### Cones against cylinders, spheres and cones (S9d.3b.1)

A cone's ruling is a line: `o + b e + w (n + k e)`, `e = cos x + sin y` in
its frame, so a cone carries S9c.2a's graphs as a cylinder does
(`curved/cones.rs`). Against the other quadric `sum (g_i . p - e_i)^2 = (r
+ t (h . p - e_h))^2` (a cone's slope `t`, zero for a cylinder or a
sphere) the ruling's quadratic `A w^2 + 2 B w + C` has coefficients that
are quadratic forms in `(cos, sin)`, and `D = B^2 - A C` a quartic form
(`turned.rs`'s forms of any degree, their chart polynomials times `(1 +
t^2)^deg`). The carrier is the first input, a cylinder before a cone,
whose `A` has no real root: `D` positive all round gives two rings over
its angle (coaxial pairs among them, `D` constant), negative all round
the surfaces apart; `A` within rounding of zero (a cylinder along a
cone's ruling) is `Degenerate`, as is `D` vanishing identically (a sphere
inscribed in a cone) or a carrier's apex on the other. Two cones whose
quadrics differ by an affine function (every ruling asymptotic to the
other: parallel axes, or one axis, with equal slopes) meet on that plane:
S9d.3a's plane section of one. Rings on a cone's wall have no seam to
cross, so one gets a vertex at a rational angle. Circles against a cone
(a cap's arc, a rim, a sphere's circle) solve the other quadric along
them with its cone term. `Curve3::Meet` takes the cone carrier and the
cone as the other quadric (`TOPOLOGY_MODEL.md`). A carrier-less pair (no
input's rulings meet the other all round: spheres off a cone's axis) or
one whose only carrier has rulings along the other's asymptotes (their
branches to infinity) is S9d.3b.2's. A cone and a sphere in loops take
S9d.2b's graphs over the height with the cone's circle of radius `b + k
w` (`Curve3::Rise` on a cone), in exact frames; a carrier whose `A` has
simple real roots is split at those directions (algebraic over
`Q(alpha)`), each branch an open piece between them, the finite branch
switched at `w = -C / 2B`; the ruling's roots taken as `C / (-B - s
sqrt(D))` where `(-B + s sqrt(D)) / A` cancels. A loop in a turned frame is
`OutOfDomain`.

### Tori against polyhedral prisms (S9d.4a)

A whole torus is `(|l|^2 + R^2 - r^2)^2 <= 4 R^2 (l_u^2 + l_v^2)` in its
frame's coordinates as rationals (`curved/torus.rs`). On it the distance
from the axis is `rho = (|l|^2 + R^2 - r^2) / 2R`, rational in a point's
coordinates, so both angles have exact places, `u` along `(l_u, l_v) /
rho` and `v` along `(rho - R, l_w) / r`. Its wall is traced in four
patches of `(u, v)`, cut at the meridians through a rational direction and
its opposite and at the parallels of a rational angle `v0` and `v0 + pi`,
every seam a circle with rational centre and axes. A line meets the torus
at its quartic's roots (`Q(alpha)`); a plane in a spiric section, a graph
over `u` (the tube's circle at `u` meeting the plane where `r A cos v + r
mu sin v = -(R A + kappa)`, one quadratic surd at a rational `u`) or over
`v` symmetrically: rings over `u` where `D_u` is positive all round (loops
about the axis), over `v` where `D_v` is (loops about the tube), else loops
of graphs over `v` about their `u` turning points and over `u` about their
`v` ones, switched at rational `u` between turning points of different
kinds and verified exactly (S9d.2b's rule). Edges are S8d.3's
`Curve3::Section`, their pcurves `Projection`s on the torus; a face's loops
lie on one sheet of the torus's cover by their material's side. A plane
tangent to the torus or within the resolution of it is `Degenerate`.

### Torus segments and wedges against polyhedral prisms (S9d.4b.1)

A v-segment or a wedge is S9d.4a's torus with its ends
(`curved/torus_segment.rs`). A v-segment (latitudes `low < high`, a full
turn) is the region between the tube's arc and the axis, revolved: its end
discs lie in the planes `w = z` at the stored heights `r sin(latitude)`,
from the axis to the arc's ends, where the planes meet the tube on the side
of the latitude's cosine (`rho = R +- sqrt(r^2 - z^2)`). Between its
critical heights (the ends and `+-r`) the tube's outer and inner points at a
height lie on the arc or not all along, decided once, so a point's side is
exact in `t^2 - R^2`, the torus's quartic and its height (the tube's disc,
`t < R + q`, `t < R - q` or nothing), pushed at first order as the whole
torus's. A wedge (the whole tube over `0 < angle < 2 pi`) is the torus in
the sector from the half-plane of `x` to that of its chart direction `(cos
angle, sin angle)` rounded: two half-planes' common within a half turn,
their union beyond. The wall is traced in two patches, cut at the ends
instead of a seam where an end lies (a segment's at its meridian seam, a
wedge's at its parallel seam). A segment's rims are its end planes' rings
over `u` (S9d.4a's spiric sections; where the plane is tangent to the torus
along the rim, the outer and inner halves' `w = +-r`, a circle of radius
`R`), a wedge's its end half-planes' rings over `v` (surds of `cos^2 + sin^2`
of the end direction); a rim meets a plane where the two planes' line meets
the torus on the rim's branch (a repeated root on the other ring or circle is
no contact), and a result's rim is the input's circle. An inside-out
segment's wall (the inner half) has its material outside the tube. A plane
tangent to the torus, or within the resolution of it, off the part's wall is
no contact: the part's section is then the two branches of the graph over
the wall's range (`v` for a segment, `u` for a wedge) where that graph's
discriminant is positive all over a range holding the wall's exactly (Sturm
counts, no extremum within the resolution of zero), or nothing where it is
negative; on the wall it is `Degenerate`. A face on a torus none of whose
loops winds holds its holes inside its outer loop on the cover (a wedge's
over more than half a turn). `Solid::classify` decides a segment or a wedge
(`decide::torus_part_location`: membership exactly, the distances from the
wall within its range and from the end discs in rational intervals).

### A whole torus against prisms with arcs, spheres and cones (S9d.4b.2a)

The other input's curved face is a quadric (`procedural::Other`: a
cylinder or cone on its affine frame, a sphere in the world), its function
`G(u, v)` on the torus's angles a trigonometric polynomial of degree two in
each (`curved/torus_curved.rs`, exact, reduced by `sin^2 = 1 - cos^2`). At
a rational `u` the roots of `G` in `v` are a quartic's in a chart's
half-angle tangent, so a point there lies in one `Q(alpha)`, and
symmetrically at a rational `v`. The critical values of `u` are the real
roots of the quartic's discriminant in `v` (degree 24 in a chart of `u`
whose antipode is none); a tube's circle wholly on the quadric (every
coefficient zero at one `u`) is `OutOfDomain`. On the line of each critical
value every box is certified clear of `G`, of `G_v` or of `G_u` (binary64
intervals, mean-value forms, halved to `1e-10` radians): a point where none
clears is a singular point of the meeting, a tangency of the surfaces,
`Degenerate`; a turning point is not. A rational `u` between each two
critical values (and the chart's antipode) is a line whose roots, exact,
seed a numerical trace of each component; runs where its slope in the
angles is at most one (`|G_u| <= |G_v|`, turning points in `v` among them)
are graphs over `u`, the others (turning points in `u`, `G_v = 0`) over
`v`, switched where the slope is one (each graph's series as far from its
own turning points as the other's), a component of one kind all round a ring
over its parameter where it winds once in it alone. A piece's
branch is a window of the other angle (rational directions) between its
own values and the other roots at samples along its run; its switches are
its points at rational parameters (on the graph over `u` beside them), its
range their directions. Each piece is verified exactly: `G` has no zero on
the window's ends over the range (Sturm counts, the range's algebraic ends
widened to rational ones outside them), one root inside the window at a
rational parameter in the range, and no double root in the rectangle
(certified boxes clear of `G` or of its derivative in the other angle); a
piece that fails is split and tried again. Every root on every line must
lie on a verified piece: no component is missed. A coaxial quadric (`G`
independent of `u`) meets the torus in circles: rings over `u` at the roots
in `v`, stored as `Curve3::Circle` about the axis. Other sections are
`Curve3::Toric` (`TOPOLOGY_MODEL.md`), their pcurves `Projection`s. A cap's,
a rim's or a whole sphere's great circle meets the torus where the torus's
quartic along it vanishes, a polynomial of degree eight in its half-angle
tangent (a sphere's circle of a surd radius against a torus is
`OutOfDomain`); a seam's tangency with a torus is retried at another seam.
A v-segment or wedge against a curved face is `OutOfDomain("a torus
segment or wedge against a curved face (S9d.4b)")`.

### Two whole tori (S9d.4b.2b)

The other surface may be a whole torus (`torus_curved::Far::Torus`: its
model's affine frame and radii), the carrier the object's torus. On the
carrier's angles its function `(|l|^2 + R2^2 - r2^2)^2 - 4 R2^2 (l_u^2 +
l_v^2)`, `l` its local coordinates, is of degree two in each angle where
both stored frames are exactly orthonormal (along a round circle `|l|^2` is
affine in its angle's cosine and sine), a quadric's degree: S9d.4b.2a's
discriminant, lines, traces, pieces and checks apply unchanged (every
fixture). In frames not exactly orthonormal (axes rounded to binary64 by a
turn: a circle of one is an ellipse in the other's coordinates) it is of
degree four in each, a point at a rational parameter algebraic of degree
eight, and the discriminant (degree 112 in a chart of `u`) out of reach.
The turning points in `u` (`G = G_v = 0`) are enclosed instead by
subdividing `[-pi, pi]^2` into boxes clear of `G` or of `G_v` (mean-value
forms in binary64 intervals) and boxes under `1e-6` radians clear of `G_u`
(the meeting regular there: a turning point), merged where they touch;
a box clear of none at `1e-10` radians is a tangency, `Degenerate`. Lines
of `u` in the gaps between the boxes seed the traces, and besides every
line's roots every turning point's box must lie inside a verified piece
over `v` (its window of `u` and its range of `v`, by a margin of `1e-9`):
the piece's one root there is the turning point, so every component (one
that does not turn in `u` winds about the axis and crosses every line) is
on the pieces. Pieces, their verification and `Curve3::Toric` (its
`other_minor` the other torus's minor radius, `TOPOLOGY_MODEL.md`) are
S9d.4b.2a's with forms of either degree; the tori's seams are circles
(`conic_torus`). Coaxial tori meet in circles, tori tangent along one
refused (`tori_kiss`). Parallel tori whose top or bottom circles lie at
one height and cross are tangent there (both normals along the axes), a
singular point of the meeting, `Degenerate`; with their equators at one
height they cross transversally (the normals along the two radii) and the
meeting only turns in `u` there. Linked tori apart fuse as two solids,
their common empty. Fields of degree eight with coefficients of a few
thousand bits (the rounded frames' inverses) make a turned pair's Boolean
take 15 to 50 s: products in `Q(alpha)` reduce by `x^j mod p` kept over
one denominator (`num::Gen::mul_mod`, an integer product and one
reduction per coefficient), a sign is tried by a binary64 enclosure over
the generator's isolator before Sturm-Tarski, and a meeting keeps its
tangents by point; results are unchanged. The parallel fixtures' Booleans
take about 20 s, most of it the certified integrals (the validator's and
the mass's) along their meetings.

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
* **Kernel (S9d.1).** `tests/sphere_booleans.rs`: the 30 fixtures as the
  reference (25 results and an empty common inside its measures, the 4
  declared degenerate refused), every history checked, results
  deterministic and moved rigidly with every vertex on their boundary.
  `compare_sphere_boolean.py`: 30 matches, counts included (OCCT's seams
  counted), no review; every other Boolean comparison unchanged.
* **S9d.2 evidence (spheres against cylinders and spheres), before its
  kernel code.** `spheres_boolean_reference.py` (mpmath, 40 digits) takes a
  sphere or cap on S9d.1's model, a prism of lines, arcs and circles
  (holes reversed) on S9c.1's, and two spheres. It slices both inputs by
  the planes normal to the prism's axis (two spheres: the first one's) in
  an affine rational chart of the slices (the prism's own `(u, v)`), where
  a sphere's section is an ellipse `q^T G q + 2 q^T g(s) + k(s) <= 0` (a
  circle in an exact frame) and a prism's the preimage of its profile, an
  arc `C + A (cos t, sin t)` of the same parameter, cut by the caps' strip
  when they are not slices: every boundary a segment or such an arc, a
  sphere's taken with `A = sqrt(rho^2) L^-1` (`G = L^T L`) so that `t` is
  its true angle. The operations apart as S9d.1's (boundaries cut where
  they cross, a segment against an arc a quadratic, two arcs the quartic
  `z^2 f(z)`, `z = e^(it)`, its unit roots refined by Newton; pieces
  classified at their midpoints), areas and moments by Green's theorem
  over segments and arcs exactly (trigonometric polynomials of degree
  three), the sphere's face by Archimedes. Breakpoints, roots of exact
  polynomials: planes tangent to a sphere, prism vertices, a zone's
  planes, edges meeting a sphere, faces' planes tangent to its sections, a
  cap circle's meetings with it (a quartic in `tan(theta / 2)`), a cap's
  trace tangent to its circle, and two sections tangent: the discriminant
  in `lambda` of `det(lambda Q1(s) + Q2(s))`, or where concentric sections
  of proportional Gram matrices coincide (a sphere centred on a cylinder's
  axis, its parallels). Other faces in closed form in their planes (caps,
  a zone's discs, a flat wall's parallelogram against the sphere's
  ellipse), a cylindrical wall by its angle (each generatrix inside the
  sphere between a quadratic's roots, breakpoints quartics in `tan(theta
  / 2)`). Solids: the result's section in each interval chained into
  loops and components keyed by their faces, followed through an interval
  by key and joined across a breakpoint when most of 64 points of the
  smaller section lie in the other. `generate_spheres_boolean_fixtures.py
  --check` writes `boolean-spheres-cases.txt`, `boolean-spheres-
  expected.tsv` and `boolean-spheres-frames.tsv` with 33 cases (4 fuses,
  15 cuts, 14 commons; 28 solid, 1 empty, 4 degenerate; 31 in exact frames,
  whole spheres in any): two spheres crossing, nested (a cavity), apart
  (fused two solids, the common empty), crossing in `TILTX` and `LEAN`; a
  sphere and a coaxial cylinder (parallels at rational heights): a pipe
  through a sphere (the spherical ring; the pipe's two ends), in `SIDE`, a
  box with a coaxial hole whose walls cut the sphere (the sphere less it
  five solids), a hemisphere against a pipe through its disc; an off-axis
  cylinder (the quartic): a rod through the sphere, two rings (in exact
  and `TILT` frames; the rod less the sphere two solids), a bite, one loop
  (exact and `TILT`), a cylinder ending inside the sphere, its cap circle
  crossing it; a stadium through a sphere (the sphere less it two solids);
  and `degenerate` two spheres tangent, a cylinder tangent outside and
  inside. Checks before writing: closed forms (two spheres' lens as two
  caps; a coaxial cylinder's axial integrals of `pi min(a^2, R^2 - z^2)`,
  the ring `pi h^3 / 6` whatever `R`; the holed box as S9d.1's box less
  the core; an off-axis cylinder by the lens of two discs, its moments by
  circular segments and each circle's angle inside the other, along the
  axis) within 9.2e-41 in exact frames and 4.2e-17 in turned ones, both
  inputs' slicings and inclusion and exclusion 9.2e-41, the area identity
  2.3e-40, every face's classes 1.8e-40 (a cylindrical wall's area by
  mpmath's quadrature of its element), a second slicing direction `(2, -3,
  5)` (the prism cut obliquely, its circles ellipses) 1.8e-40, Monte Carlo
  3.1 standard errors, no near coincidence, every cap circle at least
  0.48 of the radius squared from tangency with the sphere.
  `test_spheres_boolean_reference.py` checks the closed forms (the ring,
  a cap inside a coaxial cylinder, the lens of two spheres and of two
  discs), an ellipse's Green integrals, the tangency polynomial and the
  reference on the ring. No protocol or oracle change: either input or both
  may already be a `sphere` row. `compare_spheres_boolean.py`
  (`compare_boolean.make_set`) reproduces
  `occt-boolean-spheres-preimplementation` (`rust_spheres_boolean_exists`
  false; the kernel's probe `unsupported` on all 33, `OutOfDomain("a
  sphere against a cylinder or a sphere (S9d.2)")` against a prism,
  `OutOfDomain("a Boolean of a solid with curved faces or edges in any
  position (S9c)")` for two spheres): every result valid with the
  reference's solid count; 18 match (parallels of the sphere's own frame
  within 2.3e-14, other circles within 1.7e-8) and 15 are reviewed
  (`occt-boolean-spheres-divergences.json`): BRepGProp's default
  integration on faces bounded by approximated quartics (and on the holed
  box's circles off the sphere's parallels) misses by up to 2.4e-6, the
  same results measured adaptively within 3.4e-9; seven solids' counts
  change when unified.
* **Kernel (S9d.2b).** The loops: `tests/spheres_booleans.rs` now takes 32
  of the 33 fixtures as the reference (the bites, the cap's circle crossing
  the sphere, the stadiums), the bite in `TILT` `OutOfDomain`;
  `compare_spheres_boolean.py` 12 matches and 21 reviewed (seven more
  results whose counts differ at the loops' switch points), no failure.
* **Kernel (S9d.2a).** `tests/spheres_booleans.rs`: 24 of the 33 fixtures
  as the reference (two spheres; coaxial pipes and holes, parallels as rings
  of `Curve3::Meet` with the sphere as its other quadric; rods through
  spheres, two rings; the 4 tangencies refused), the 9 loops (a bite, a
  cap's circle crossing the sphere, the stadiums) `OutOfDomain` (S9d.2b's),
  every history checked, results deterministic and moved rigidly.
  `compare_spheres_boolean.py`: 12 matches, 21 reviewed (the measures, now
  with ten kernel results whose faces are OCCT's unified ones and whose
  edges and vertices differ where each splits its circles and quartics), no
  failure; every other Boolean comparison unchanged.
* **S9d.2c evidence (turned caps' circles against cylinders, loops in
  turned frames), before its kernel code.** `spheres_boolean_reference.py`
  reads a cap's end planes as the kernel's `Ball` does (`AxisSphere`:
  through `o + h n` normal to the stored axis, S9d.1's affine plane in an
  exact frame) and slices a turned cap along its own axis, the prism cut
  obliquely: its caps classified in their own `(u, v)` against the ball's
  ellipse and zone half-planes, a cylindrical wall's zone bound varying
  along its generatrices (events where the rim crosses the wall, a quartic
  in `tan(theta / 2)`, or the bound meets the prism's ends), components of
  one key numbered by their centres. `generate_spheres_turned_boolean_
  fixtures.py --check` writes `boolean-spheres-turned-cases.txt`,
  `-expected.tsv` and `-frames.tsv` with 18 cases (3 fuses, 7 cuts, 8
  commons; 17 solid, 1 degenerate): hemispheres in `TILT`, `LEAN` and
  `TILTX` against a coaxial pipe, S9d.2's bite and S9d.2's rod in `XY`;
  loops of the bite in `LEAN` and `TILT`, of a thick cylinder in `TILTX`
  and a thin one in `R125`; a hemisphere in `TILT` against the bite in
  `LEAN`; a hemisphere's rim tangent to a cylinder within rounding,
  `degenerate`. Checks before writing: closed forms (S9d.2's lens along a
  turned cylinder's axis; the hemisphere against the coaxial pipe by
  circular segments along the axis, its faces by central symmetry and a
  circle and an ellipse) within 8.1e-17, a hemisphere and its complement
  against the whole sphere sliced along the prism's axis within 1.7e-41,
  inclusion and exclusion, the area identity and every face's classes
  within 1.4e-40, a second direction for whole spheres, Monte Carlo 2.2
  standard errors, no near coincidence but the declared pair's.
  `test_spheres_turned_boolean_reference.py` checks the closed forms by
  quadrature, the kernel's cap reading, the oblique wall against the axial
  one and the halves. `compare_spheres_turned_boolean.py` reproduces
  `occt-boolean-spheres-turned-preimplementation`
  (`rust_spheres_turned_boolean_exists` false; the kernel `unsupported` on
  all 18): every result valid with the reference's solid count, 1 match,
  17 reviewed (`occt-boolean-spheres-turned-divergences.json`: BRepGProp's
  default integration on approximated quartics, up to 1.4e-5; adaptively
  converged within 1.5e-8, or 1.0e-8 absolute on a small result within its
  edges' tolerance), five solids' counts change when unified.
* **Kernel (S9d.2c).** `tests/spheres_turned_booleans.rs`: all 18 fixtures
  as the reference (17 within the kernel's enclosures, the rim tangent
  within rounding refused), every history checked, results deterministic
  and moved rigidly, the turned loops' `Curve3::Rise` pieces within `1e-12`
  of both surfaces, the turned caps' rim vertices on the rim's plane, a
  turned cone's loop and a turned cap's circle against a cone `OutOfDomain`
  (S9d.3c). `compare_spheres_turned_boolean.py`: 0 matches, 18 reviewed
  (the measures, and fifteen kernel results whose faces are OCCT's unified
  ones and whose edges and vertices differ where each splits its loops),
  no failure; `compare_spheres_boolean.py` unchanged in its counts with
  `bite_tilt_cut` now within the reference; every other Boolean comparison
  unchanged. `tests/spheres_booleans.rs` takes all 33 of S9d.2's fixtures.
* **S9d.3a evidence (cones against polyhedral prisms), before its kernel
  code.** `cone_boolean_reference.py` (mpmath, 40 digits) takes the cone or
  frustum on its exact model (the stored frame's axes as rationals, `0 <= w
  <= h` and `u^2 + v^2 <= r(w)^2`, `r` linear from `bottom` to `top`, a
  zero radius an apex) and the prism on S9d.1's, and works in the cone's
  chart `(u, v, w)` (volumes and moments carried back by `det(x, y, n)`,
  planar areas by Nanson's formula). It slices both by planes normal to the
  axis (each slice of the cone a disc of radius `r(s)`, the end planes
  slices) or tilted from it by less than the cone's complement (every
  section an ellipse cut by the end planes' lines, made a disc by
  coordinates orthonormal for the quadric's form there, a Cholesky factor):
  S9d.1's classification of the disc's and the pieces' boundaries, the
  four operations apart and the common again by clipping each piece, a
  line within 1e-30 of tangency taken as tangent (a face tangent along a
  ruling touches every slice's circle); the wall by its area element `r(w)
  N(theta) dtheta dw` (`N = sqrt(1 + k^2)` in an orthonormal frame, else
  `|cos theta (y x n) - sin theta (x x n) - k (x x y)|` by Gauss-Legendre
  along each arc) over the arcs inside; breakpoints the roots of exact
  quadratics (vertices, edges meeting the quadric and the end planes, a
  face's line tangent to the section as the discriminant along the face's
  and the slice's line, the end circles' extreme levels, the apex, a face's
  and an end plane's line meeting the quadric), Gauss-Legendre between them
  refined to 1e-33 of the case's size to the fourth. Planar faces by their
  own parameters: normal to the axis against the level's disc in closed
  form (or on an end plane, the same or the opposite orientation), any
  other sliced by lines of constant `w` (the chord inside the disc a
  quadratic's roots, the spacing `|a| / |a_uv|`); the end discs against
  the prism's sections in closed form. Solids by convexity: the common's
  pieces of positive volume joined across internal faces inside the cone,
  `K - P` as the components of the open sets `int K n H_f` (joined when
  `int K n H_f n H_g` has volume), `P - K` by runs of the faces' boundaries
  outside the cone. `generate_cone_boolean_fixtures.py --check` writes
  `boolean-cone-cases.txt`, `boolean-cone-expected.tsv` and
  `boolean-cone-frames.tsv` with 30 cases (4 fuses, 14 cuts, 12 commons;
  25 solid, 1 empty, 4 degenerate; 24 in exact frames): a frustum cut by a
  face oblique to its axis (an ellipse, all three operations), by a face
  parallel to a ruling (a parabola clipped by both end planes), by a box's
  wall parallel to the axis (a hyperbola); a box through a cone's apex
  region, its bottom a circle and a wall a hyperbola, and a box's corner
  below the apex; a cone inside a box and a box inside a frustum
  (cavities); a half-space normal to the axis (a frustum, all three); a
  slab through the axis and a bar across it (two solids); a box on the top
  disc (coplanar, opposite; the common empty) and on the base plane (the
  same orientation); in turned frames a frustum in `TILT` below a box's
  face (the ellipse crossing the base), a cone in `LEAN` with its apex in a
  box (an oblique cone on an ellipse), a frustum in `TILTX` inside a box and
  across a slab (two solids); and `degenerate` a box's wall through the
  apex (lines), a face tangent along a ruling, a box's vertex on the wall.
  Checks before writing: closed forms of every pair in the cone's ideal
  frame (an aligned box by S9d.1's rectangle-in-disc antiderivatives and
  the walls' hyperbolic chords `|[v0, v1] n [-q, q]|`, `q^2 = r^2 - u^2`; a
  half-space by circular segments `r^2 acos(d/r) - d sqrt(r^2 - d^2)`, its
  moment `2/3 (r^2 - d^2)^(3/2)`, the wall's angle `2 acos(d/r)`, the
  face's chord; a slab as two half-spaces; normal to the axis a frustum)
  within 1.2e-40 in exact frames and 1.2e-16 in turned ones, both inputs'
  slicings and inclusion and exclusion 9.2e-41, the common two ways
  2.2e-42, the area identity 1.8e-40, every face's classes 9.2e-41, a
  second slicing direction `(1, -2, m)` (every section an ellipse; pairs
  without coplanar faces) 1.8e-40, Monte Carlo 2.7 standard errors, no near
  coincidence. `test_cone_boolean_reference.py` checks the closed forms
  (a cone's and a frustum's measures, the frustum above a normal plane, an
  oblique cone as a third of its ellipse times the apex's height, its wall
  as the projected ellipse over `sin(half angle)`, the hyperbolic segment
  by its antiderivative `r^3 acos(d/r)/3 - 2 d r q / 3 + d^3 ln(r + q)/3`, a
  half cone's centroid `R / pi`) and the reference on a frustum. The
  protocol takes a cone on either side (`encode_boolean_case` accepts a
  tool's `frame` and `cone BOTTOM TOP HEIGHT` rows, which the kernel's test
  support already reads); `native_case` gives it one `cone` row, built by
  `BRepPrimAPI_MakeCone` in `occt_boolean_oracle.cpp`.
  `compare_cone_boolean.py` (`compare_boolean.make_set`) reproduces
  `occt-boolean-cone-preimplementation` (`rust_cone_boolean_exists` false;
  the kernel's probe `unsupported` on all 30, `OutOfDomain("a Boolean of a
  solid with curved faces or edges in any position (S9c)")`): every result
  valid with the reference's solid count, all 30 within 1.4e-8 (BRepGProp's
  default integration; adaptively within 5.0e-9, every section an exact
  conic), no review; thirteen solids' counts change when unified (the
  wall's faces split at OCCT's seam).
* **S9d.3b evidence (cones against curved faces), before its kernel
  code.** `cones_boolean_reference.py` (mpmath, 40 digits) takes a cone or
  frustum on S9d.3a's model, a prism of lines, arcs and circles on S9d.2's,
  a sphere, cap or zone on S9d.1's, or two cones, each input a set of exact
  surfaces `X^T M X + 2 m . X + c` in world coordinates (the cone's quadric
  from the rows of its frame's inverse, its end planes; the sphere and its
  zone planes; the prism's walls, its arcs' cylinders and its caps). It
  slices both inputs in S9d.2's affine rational chart along a direction in
  which every quadric's section is an ellipse (a cone's own axis where the
  other allows it, else a prism's, else a rational combination of the
  axes, checked exactly: the cone's form positive definite on the slices,
  no slice parallel to a prism's axis, a cap's planes slices): a cone's
  section `C + sqrt(rho) L^-1 (cos t, sin t)` (`G = L^T L`) cut by its end
  planes' lines, the others S9d.2's. The operations apart and Green's
  integrals as S9d.2's; the sphere's face by Archimedes; the cone's wall a
  second way where the slices are normal to its axis (its element over its
  circle's arcs inside). Breakpoints: vertices, every edge of one input (a
  line; a circle as a conic in its plane) meeting every surface of the
  other (a quadratic; the resultant of the two traces in one coordinate, a
  quartic), every plane tangent to every quadric's section, two sections
  tangent (S9d.2's pencil, square-free), a section shrinking to a point, the
  planes that are slices. Every face but the sphere's is swept by lines on
  it (a cone's rulings, with their element `r(w) |U' x (n + k U)|`, a
  cylinder's and a flat wall's generatrices, parallel lines of a skew
  direction on a planar face), each cut where every surface's quadratic in
  the line's parameter has a root and its pieces classified (a coplanar face
  of the other apart as the same or the opposite orientation); breakpoints
  in the sweep's parameter the roots, at 100 digits, of exact polynomials:
  each surface's leading coefficient and discriminant along the lines and
  every two surfaces' resultant (trigonometric along rulings, through
  `tan(theta / 2)`). Solids as S9d.2's, components of one key (both ends of
  a rod in one tilted slice) followed by their centroids, sections vanishing
  on both sides of a breakpoint (two apexes, a pinch along a circle) not
  joined; equal cones as the cone. `generate_cones_boolean_fixtures.py
  --check` writes `boolean-cones-cases.txt`, `boolean-cones-expected.tsv`
  and `boolean-cones-frames.tsv` with 40 cases (7 fuses, 14 cuts, 19
  commons; 35 solid, 1 empty, 4 degenerate; 36 in exact frames): a coaxial
  pipe through a frustum (circles, all three; the pipe less it two
  solids), a cylinder through its wall and base off the axis, a rod across
  it (two rings; the rod less it two), a rod across a cone's tip (the cone
  less it two) and around its apex, a box with a coaxial hole, a stadium
  (hyperbolas and quartics), a pipe in `TILT`; coaxial spheres through a
  frustum's top disc (all three) and on a cone's apex, a sphere off the
  axis (a loop), a dome through which a cone passes, a sphere against a
  frustum in `TILT`; coaxial cones tip in tip, frusta of opposite slopes, a
  frustum standing on another (coplanar discs, the common empty), a cone
  across a frustum (two rings; its ends two solids), a cone in `LEAN`, two
  cones of parallel axes; and `degenerate` a cylinder tangent along a
  ruling, a sphere inscribed along a circle, two cones apex to apex, equal
  cones. Each expected row names the sub-step that would take it if
  S9d.3b were split as S9d.2 was: `S9d.3b.1` (36) where some input's
  rulings meet the other's quadrics transversally wherever their curve runs
  (circles, and graphs over a cylinder's or a cone's angle), `S9d.3b.2` (4:
  the spheres off the axis, the inscribed sphere) where both inputs'
  rulings, or the cone's against a sphere, are tangent to the other on
  their faces, decided by the rulings' exact discriminants. Checks before
  writing: closed forms of 22 of the 27 pairs along the cone's axis in its
  ideal frame (a cylinder, a sphere or cap, or a cone of parallel axis,
  coaxial or not, by the lens of two discs: areas, first moments by
  circular segments, each circle's angle inside the other with each wall's
  element, a sphere's by Archimedes, end discs by the lens, coplanar ones
  apart; a rod across the axis by the strip `|y - y0| <= sqrt(a^2 - (z -
  z0)^2)` against the disc and the rod's generatrices' chords; the holed
  box as S9d.3a's box less the coaxial core; equal cones) within 1.7e-40 in
  exact frames and 4.4e-17 in turned ones, both inputs' slicings 8.6e-42,
  inclusion and exclusion 1.2e-41, the area identity 1.1e-40, every face's
  classes 5.0e-41, both sides' shared areas 5.7e-42, the cone's wall two
  ways 3.6e-41, a second slicing direction (19 pairs: none with a cap,
  coplanar faces or a declared degeneracy) 4.6e-41, Monte Carlo 2.9
  standard errors, no near coincidence, every edge and vertex of one input
  at least 0.0126 (relative) from tangency with or incidence on a surface
  of the other. `test_cones_boolean_reference.py` checks the closed forms
  (a coaxial cylinder and an inscribed sphere by hand, two cones of
  parallel axes by the lens of equal discs, the strip by a double
  quadrature), the resultant and the roots of trigonometric polynomials,
  the rulings' quadratics and the reference on the coaxial pipe. No protocol
  or oracle change: either input or both may already be a `cone` or
  `sphere` row. `compare_cones_boolean.py` (`compare_boolean.make_set`)
  reproduces `occt-boolean-cones-preimplementation`
  (`rust_cones_boolean_exists` false; the kernel's probe `unsupported` on
  all 40, `OutOfDomain("a cone against a prism with arcs, a sphere or a
  cone (S9d.3b)")`): every result valid, the reference's solid count in 39
  (`ball_inscribed_cut`, declared degenerate, one solid touching itself
  along the circle natively); 23 match (sections circles about the axis
  within 4.1e-14; two cones of equal half-angle and parallel axes meet in
  a plane, OCCT's hyperbolas within 6.1e-10) and 17 are reviewed
  (`occt-boolean-cones-divergences.json`): BRepGProp's default integration
  on faces bounded by approximated quartics misses by up to 4.1e-6, the
  same results measured adaptively within 6.0e-9, and the inscribed
  sphere's count; eleven solids' counts change when unified.
* **S9d.4a evidence (tori against polyhedral prisms), before its kernel
  code.** `torus_boolean_reference.py` (mpmath, 40 digits) takes a whole
  torus on its exact model (the stored frame's axes as rationals, `(|p|^2 +
  R^2 - r^2)^2 <= 4 R^2 (u^2 + v^2)` in its chart, `R > r`) and the prism on
  S9d.1's, and works in the torus's chart (volumes and moments carried back
  by `det(x, y, n)`, planar areas by Nanson's formula, the wall by its
  element `r (R + r cos phi) |cof(A) N|` in its own `(theta, phi)`). It
  slices both inputs two ways. Normal to the axis, each slice of the torus
  is the annulus between the circles of radii `R -+ sqrt(r^2 - s^2)`:
  S9d.3a's classification run on both discs against the pieces' sections,
  the operations their combinations, the common again by clipping each piece
  by both discs, the wall by the circles' angles inside with the element `r
  rho / q`; breakpoints the pieces' levels, `+-r`, an edge's meetings with
  the torus and a face's line tangent to either circle (quartics, their real
  roots through Yun's square-free factorization and `polyroots`). By the
  meridian half-planes about the axis, each section of the torus is the
  tube's disc about `(R, 0)`, the prism's a convex polygon clipped at the
  axis: volumes and moments by Green's integrals of `t`, `t^2` and `t w`
  (the cylindrical element), the wall over the disc's arcs inside by its own
  element; breakpoints the vertices' angles, edges' and faces' lines meeting
  the torus, a face's line tangent to the tube's circle (the decisions'
  spiric quadratic in `alpha = a_u cos theta + a_v sin theta`) or parallel
  to the axis. Planar faces by their own parameters as S9d.3a's, the chord's
  part inside the annulus. Solids (a convex prism) followed through the
  meridians: the common's convex sections, `D - C` as the union of the
  discs' parts beyond each face (joined where two meet), `C - D` by runs of
  the section's boundary outside the disc, each holding a vertex; joined
  across a breakpoint where a face's part of the disc, the common, or a
  vertex outside the disc persists there. Fuse: one solid when the inputs
  overlap. `generate_torus_boolean_fixtures.py --check` writes
  `boolean-torus-cases.txt`, `boolean-torus-expected.tsv` and
  `boolean-torus-frames.tsv` with 35 cases (4 fuses, 15 cuts, 16 commons; 30
  solid, 1 empty, 4 degenerate; 29 in exact frames), the torus of radii 5/2
  and 3/2: a bar across the whole torus through the hole, its walls cutting
  loops about the tube (all three; two chunks, two halves), a strip along
  the equator (less the torus three solids); a slab normal to the axis
  (loops about the axis; two solids both ways) and a half-space above the
  equator; a cap cut from the tube's outside by a wall, a box through the
  hole whose four walls cut such caps (four solids), caps cut from the
  tube's inside by a box's vertical edges (four) and by a wedge's edge; a
  half-space through the axis (loops about the tube, all three); a box
  inside the tube, the torus inside a box (cavities), a box in the hole
  (fused two solids, the common empty); in turned frames a torus in `TILT`
  above a plane of a Villarceau plane's inclination 3/4 from the centre (one
  contractible loop), a torus in `TILTX` across a slab between its saddle
  levels (two loops about the tube in each plane: two solids both ways), a
  bar in `LEAN` across the torus; and `degenerate` a face tangent along the
  top circle, a wall tangent to the inner equator (a figure eight), a face
  on a Villarceau plane (tangent at two points) and a box's vertex on the
  torus. Checks before writing: closed forms of 19 of the 20 pairs along the
  axis in the torus's ideal frame (an aligned box by S9d.1's
  rectangle-in-disc antiderivatives on both circles, the wall by its
  latitude, the box's walls by their chords inside the annulus; the torus
  less half-spaces by circular segments of both circles, normal to the axis
  `4 pi R q` per slice) within 3.3e-40 in exact frames and 1.5e-16 in turned
  ones (the degenerate pairs within 2.3e-29), both inputs' slicings 1.4e-40,
  inclusion and exclusion 1.1e-40, the common two ways 4.6e-41, the area
  identity and every face's classes 3.0e-35 (1.3e-39 but for the vertex on
  the torus), the meridians as a second direction 1.6e-40, the wall two ways
  2.1e-37, Monte Carlo 2.8 standard errors, no near coincidence outside the
  degenerate pairs, every other pair's vertices at least 0.18, faces' planes
  0.031 and edges' crossings (a sine) 0.54 from tangency with the torus.
  `test_torus_boolean_reference.py` checks the closed forms by Pappus's
  theorems (the torus, the half above the equator and the half beyond a
  plane through the axis with its centroid `(4 R^2 + r^2) / (2 pi R)`, a
  band) and the reference on a half-space through the axis, a band, a slab
  and a box in the hole. The protocol takes a torus on either side
  (`encode_boolean_case` accepts a tool's `frame` and `torus MAJOR MINOR LOW
  HIGH ANGLE` rows, which the kernel's test support already reads;
  `native_case` gives a whole torus one `torus` row, built by
  `BRepPrimAPI_MakeTorus(gp_Ax2, R, r)` in `occt_boolean_oracle.cpp`; every
  older capture reproduces). `compare_torus_boolean.py`
  (`compare_boolean.make_set`) reproduces
  `occt-boolean-torus-preimplementation` (`rust_torus_boolean_exists` false;
  the kernel's probe `unsupported` on all 35, `OutOfDomain("a Boolean of a
  solid with curved faces or edges in any position (S9c)")`): every result
  valid with the reference's solid count; 16 match (circles within 4.7e-16)
  and 19 are reviewed (`occt-boolean-torus-divergences.json`): BRepGProp's
  default integration on faces bounded by approximated spiric sections
  (B-splines, every section but circles) misses by up to 6.3e-6, the same
  results measured adaptively within 1.8e-9 but for a small cap (1.2e-8) and
  the Villarceau plane (2.0e-8); nine solids' counts change when unified.
* **S9d.4b.1 evidence (torus segments and wedges against polyhedral
  prisms), before its kernel code.** `torus_segment_boolean_reference.py`
  (mpmath, 40 digits) takes a torus other than a whole one on the model S3
  builds (`Topology::torus`, as `BRepPrimAPI_MakeTorus(gp_Ax2, R, r, low,
  high, angle)`): a v-segment is the region between the tube's arc from
  `low` to `high` and the axis, revolved, its end faces planar discs normal
  to the axis at the stored heights `r sin(latitude)`, from the axis to the
  arc's end (not a cone, cylinder or plane swept by the tube's point); in
  the meridian half-plane its section `M` has an odd number of the arc's
  points at its height to its right (the disc `[0, c]` between the end
  heights, the tube's cap `[R - q, R + q]` beyond); a wedge is the whole
  tube between the half-plane of `x` and that of the rounded `(cos angle,
  sin angle)` of its chart. Both slicings of S9d.4a with `M`: normal to the
  axis, each slice a disc or an annulus (a wedge's cut to its sector, two
  convex sectors beyond a half turn), the common by clipping each polygon
  to the sector and the signed discs, the rest by inclusion and exclusion;
  by the meridian half-planes, `M` against the polygon with its boundaries
  cut and classified (the arc by the polygon's lines, the polygon's edges
  by the tube's circle and the end heights, `M`'s arc signed by its
  orientation carrying all of Green's integrals: the axis has `t = 0`, the
  end segments `dw = 0`), a wedge's section empty outside its turn.
  Breakpoints S9d.4a's with the end heights, the rings and the end
  half-planes and circles against the prism's edges and faces. Every face
  in its own parameters: the wall two ways, a segment's end discs exactly
  at their planes (`same` or `opp` where a prism face lies there) and by the
  meridians' radii, a wedge's exactly in their half-planes and by the
  normal slices' chords, the prism's faces by `K`'s regions just above and
  below a level (normal to the axis) or by chords (a face on a wedge's end
  half-plane coplanar on the end's ray). Solids by sweeping the meridian
  half-planes: each section's components in bands of height, sections
  sampled between breakpoints joined one to one (refined otherwise), across
  a breakpoint only through its own section, around the turn and through
  the axis; the fuse by overlap or a shared face, the sweep agreeing.
  `generate_torus_segment_boolean_fixtures.py --check` writes
  `boolean-torus-segment-cases.txt`, `-expected.tsv` and `-frames.tsv` with
  29 cases (9 fuses, 9 cuts, 11 commons; 25 solid, 1 empty, 3 degenerate;
  26 in exact frames) on radii 5/2 and 3/2: the outer half (a barrel, its
  end discs tangent to the wall along their rings) with a bar through the
  axis (all three) and a box across its upper end disc; the inner half (a
  spool) with a slab through its waist (two solids) and a column through it
  (the column less the spool six solids); an upper band (latitudes
  0.5..2.25) with a box across its lower end disc (all three), a box
  standing on its upper end disc (coplanar, opposite) and a box through its
  axis; a half turn with a bar along both end discs (the common two chunks)
  and a slab across it; a quarter turn with a box across its start disc, a
  box in the hole cutting a cap over both ends, a box beyond its turn and a
  box against its start disc (coplanar, opposite: fused one solid, the
  common empty); in turned frames the outer half in `TILT` above a plane
  and the half turn in `TILT2` across a slab; `degenerate` a face on the
  outer half's end plane (tangent to the wall along the ring), a face
  tangent to the band along its top circle and a face tangent to the
  quarter turn at its start circle's outermost point. Checks: closed forms
  of 17 of the 19 pairs along the axis in the part's ideal frame, its
  sections given by hand (an aligned box by S9d.1's rectangle-in-disc
  antiderivatives, a wedge's box clipped to its sector, a half-space by
  circular segments, inputs meeting on a face or at a point by their sums)
  within 8.2e-41 in exact frames, 1.3e-16 in turned ones and 4.3e-16 where a
  wedge's rounded end (1.2e-16 or 6.1e-17 off its ideal plane) cuts the box;
  every operation two ways 6.4e-41, inclusion and exclusion 5.8e-41, the
  wall two ways 7.0e-37, the end discs two ways 2.9e-41, the area identity
  (`+ 2 opp`) and every face's classes 2.8e-40, Monte Carlo 3.3 standard
  errors, no near coincidence outside the degenerate pairs, every other
  pair's vertices 0.072, faces' planes 0.063 and edges' crossings (a sine)
  0.28 from tangency and incidence, rims crossing faces at a sine of 0.17.
  `test_torus_segment_boolean_reference.py` checks the closed forms by
  Pappus's theorems (both halves, a band, a wedge's centroid) and the
  reference on a half-space through the axis, a slab through the spool and
  a box beyond a wedge. The protocol already carried a part's latitudes and
  turn (`torus MAJOR MINOR LOW HIGH ANGLE`, which the kernel's test support
  reads and `Solid::torus_with` builds for every fixture); `native_case`
  gives a segment or wedge the long row `torus ... R r LOW HIGH ANGLE`, built
  by `BRepPrimAPI_MakeTorus(gp_Ax2, R, r, a1, a2, angle)` and reversed when
  its volume is negative (OCCT's inside-out inner half), a whole torus
  keeping the short row (every older capture reproduces).
  `compare_torus_segment_boolean.py` reproduces
  `occt-boolean-torus-segment-preimplementation`
  (`rust_torus_segment_boolean_exists` false, keyed on
  `solid/boolean/curved/torus_segment.rs`; the kernel's probe `unsupported`
  on all 29, `OutOfDomain("a Boolean of a torus segment or wedge
  (S9d.4b)")`): every result valid with the reference's solid count; 16
  match (lines and circles within 2.6e-16) and 13 are reviewed
  (`occt-boolean-torus-segment-divergences.json`): BRepGProp's default
  integration on faces bounded by B-spline sections misses by up to 1.7e-5,
  adaptively within 5.7e-9 but for two results whose approximated sections
  bound a region 3.7e-8 and 3.2e-8 off (unchanged at an accuracy of
  1e-12); nine solids' counts change when unified. With the kernel's module
  every result lies within the reference (its enclosures within 1e-9) and
  the three degenerate cases are refused: 15 match and 14 are reviewed,
  seven of them for entity counts too (OCCT's B-spline sections split at
  its own points, the kernel's in exact pieces over `u` and `v` switched at
  rational points; `ub_hole_cut`, a native match, by its counts alone).
* **S9d.4b.2 evidence (a whole torus against prisms with arcs, spheres,
  cones and tori), before its kernel code.**
  `torus_curved_boolean_reference.py` (mpmath, 40 digits) takes each input's
  exact model (the torus S9d.4a's, the cone S9d.3a's, the sphere S9d.1's, a
  convex prism of segments and arcs S9d.2's), each a set of exact surfaces
  of degree one, two or four with a membership test, and sweeps every face
  of both inputs by two families of circles and lines: the tori's tube
  circles and parallels, the sphere's meridians and parallels, a cylinder's
  or cone's rulings and sections, a flat wall's generatrices and rows, a
  planar face's chords in two directions. Along a curve each surface of the
  other input is a trigonometric polynomial of its degree (from `2d + 1`
  samples; real roots the unit-circle roots of `z^d f`) or a line's
  polynomial; the curve is cut at the roots (a tangency's double root split
  by rounding dropped), each piece classified inside it, and the pieces'
  areas, volumes and first moments integrated by the divergence theorem
  (closed forms but the area element, by 24-point Gauss-Legendre). The
  outer parameter is integrated between its events, where the curves'
  structure (their classes with the surfaces bounding them) changes: a scan
  of 360 curves, bisection to 1e-35 (1e-20 on a stadium's tangent edges,
  where two surfaces' roots coincide to second order), every quadrature
  node's structure checked against its interval's. Results: common the
  `in` pieces of both, fuse the `out` pieces, cut the object's `out` and the
  tool's `in` reversed. Solids by sweeping the torus's normal slices (its
  annuli) along 512 rays from the axis in 200 slices, intervals joined where
  they overlap. `generate_torus_curved_boolean_fixtures.py --check` writes
  `boolean-torus-curved-cases.txt`, `-expected.tsv` and `-frames.tsv` with
  44 cases (10 fuses, 16 cuts, 18 commons; 39 solid, 2 empty, 3
  degenerate; 41 in exact frames) on the torus of radii 5/2 and 1: coaxial
  cylinders through the hole (all three; the pipe first) and around the
  tube, coaxial spheres about the centre and on the axis, a coaxial cone and
  a coaxial torus, all meeting in circles; a rod across the torus (the
  common two solids, the rod less the torus three), a vertical cylinder and
  a frustum through the tube, a stadium prism across it, a sphere on the
  tube's top (all three) and one swallowing the tube (either first), a cone
  in `LEAN` whose apex lies in the tube, a torus linked with it apart (the
  fuse two solids), one ringing the tube, one of a parallel axis (the common
  two solids) and one linked through the tube; a pin and a ball in the hole
  apart; a coaxial pipe and a sphere with the torus in `TILT`; `degenerate`
  a cylinder tangent along the outer equator, a sphere touching the tube at
  a point and a coaxial torus touching along a circle. Checks: closed forms
  of 23 of the 27 pairs by one quadrature along the axis of the annulus
  against the other's sections (lenses of signed discs about a parallel
  axis with each circle's angle inside the other, walls by their own
  elements and tori by their latitudes; a rod's strip by rectangle-in-disc
  antiderivatives; sums apart) within 7.5e-39 in exact frames and 4.6e-17
  in turned ones; every operation two ways (each face's two families)
  6.6e-35, both inputs from their faces 1.4e-39 of their closed forms,
  inclusion and exclusion 1.1e-40, the area identity 1.4e-39, Monte Carlo
  2.8 standard errors, the degenerate pairs 3.8e-21; no near coincidence
  outside them, every other pair's surfaces meeting at a sine of at least
  0.42, 0.066 apart where they do not, its edges crossing the torus at a
  sine of 0.62 and 0.031 from tangency, the apex 0.13 inside.
  `test_torus_curved_boolean_reference.py` checks the closed forms against
  Pappus's theorems (the torus, its part inside a coaxial pipe through the
  hole) and the reference on that pipe, a ball holding the torus, a ball
  in the hole and the rod's pieces. No protocol row is new: the kernel's
  test support and `occt_boolean_oracle.cpp` already built a torus, sphere,
  cone or prism on either side (every older capture reproduces).
  `compare_torus_curved_boolean.py` reproduces
  `occt-boolean-torus-curved-preimplementation`
  (`rust_torus_curved_boolean_exists` false, keyed on
  `solid/boolean/curved/torus_curved.rs`; the kernel's probe `unsupported`
  on all 44, `OutOfDomain("a torus against a curved face (S9d.4b)")`):
  every result valid with the reference's solid count; 24 match (circles
  within 1.6e-13, three with B-spline sections within 1.9e-8) and 20 are
  reviewed (`occt-boolean-torus-curved-divergences.json`): BRepGProp's
  default integration on faces bounded by B-spline sections misses by up
  to 9.8e-6, adaptively within 4.8e-9 and unchanged at an accuracy of
  1e-12; thirteen solids' counts change when unified. With the kernel's
  module (S9d.4b.2a) every case but the tori lies within the reference (its
  enclosures within 1e-9), the two degenerate ones refused and the 12 tori
  `unsupported` (S9d.4b.2b's): 19 match and 25 are reviewed, 20 of them
  for the kernel's entity counts (OCCT's B-spline sections split at its own
  points and faces split at the torus's seams, the kernel's sections in
  exact pieces switched at rational points; five, native matches, by their
  counts alone: four coaxial pairs, whose bands cross OCCT's seam parallel,
  and `bore_cut`). With S9d.4b.2b every case lies within the reference, the
  three degenerate ones refused: 15 match and 29 are reviewed, the eight
  tori with unified counts other than the kernel's by them (coaxial and
  ringing pairs whose bands cross OCCT's seams, the others' B-spline
  sections split at OCCT's points; four of them native matches, reviewed by
  their counts alone).
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
`copy`, gives the same shape. Since S9d.4b's survey each `pcone`,
`psphere` and `ptorus` is built under an operation of its own, so two of
them share no ids, and a copy of one sharing ids with the other argument
is built again from its constructor's numbers in its own frame (before,
every one was built under the unspecified operation, and two in one
Boolean were refused as solids other than prisms sharing ids). A `trotate` by whole quarter turns about a
coordinate axis turns a prism's frame exactly (a signed permutation of
coordinates, the origin by DRAW's location arithmetic): the kernel's
rotation rounds the cosine of a quarter turn to 6.1e-17, as OCCT's
`gp_Trsf` does, and OCCT's tolerances absorb it, while the kernel's exact
decisions would find a wall turned onto another's plane tilted off it (a
box and its quarter-turned copy fused into an L would keep a crease: 9
unified faces where OCCT has 8). Since S9d.3b.1's survey other solids
(cones, spheres, tori) are turned by the quarter turn's exact matrix
too: a cone turned by the rounded rotation crossed a cylinder's cap
within rounding of parallel. The derived cases `boolean_prisms`,
`boolean_stacks` (a step, a pocket, a box cut in two by a slab, a closed
cavity, a tool through a round wall) and `boolean_polyhedra` (quarter
turns, a bar turned 45 degrees through a box, a tilted bar cutting a box
in two, a turned box inside another) and 975 cases of upstream's
`boolean` group (86 of them stacks, 300 polyhedra of two boxes, one
turned, 35 pockets cut one after another, a stack or a polyhedron given
to the next `bcut`, 65 of S9c.1's prisms with arcs in any position: a
cylinder and a turned box, a cylinder turned about its axis, equal
cylinders crossed at right angles, 16 of S9c.2a's: a cylinder of radius
0.5 through one of radius 1 at right angles, two quartic rings, and 8 of
S9c.2b.2's: equal parallel cylinders, one moved and turned off whole
quarter turns, meeting in generatrices, and 1 of S9d.1's: a sphere placed
on a DRAW `plane`, which the adapter takes since S9d.1's survey, and a
box whose wall cuts a cap off it, and 20 of S9d.2's: a cylinder and a
sphere centred on its top cap holding the cap, meeting the wall in a
parallel, the sphere turned by quarter turns, and 36 of S9d.3a's: a box
and a frustum inside it, standing on it or through a face, meeting its
faces in circles, or with its axis outside a wall, meeting the walls
across it in hyperbolas, in `ZH3` and `ZH4` turned about its axis, and
16 of S9d.1's pole follow-up: the cylinder and a sphere of radius 2 on
its cap turned so the cap's plane holds its axis, the section through
its poles, and 36 of S9d.3b.1's: the cylinder and a coaxial frustum on
its cap, inside it or through its caps, meeting it in circles, or a
frustum across it, its axis crossing the cylinder's at right angles,
meeting the wall in quartics, its wide end through the caps in `ZK9`
and `ZL1`, and 20 more of S9d.3b.1's, two `pcone`s: a frustum and a
narrower coaxial one on its top disc, inside it or through its discs, and 16 of
S9d.4b.2a's: the cylinder and a coaxial torus whose tube its wall cuts in
two circles) evaluate on both backends, and
S9c.2b.1 adds none. One of that
sphere's turns, `ZI5`, was wrong until a sphere face's closing chord at
a pole was enclosed narrowly (its volumes off by `32 pi / 9`, its `btuc`
refused by the kernel's validation), and is registered since. A volume
audit of the registered cases found Rust's volumes and centres of
gravity native DRAW's to its printed digits, or where they differ (eight
cases) nearer the closed forms; run again in S9d.4b's survey, Rust's
values were unchanged bit for bit, the new cases' the closed forms'
within 2.6e-15 relative. Eight more of S9d.3b.1's, a wider
frustum across the cylinder, evaluate right but take 23 to 120 seconds
on the debug worker, past the contract's 30, and are not registered. Of
the upstream cases in frames with different axes the rest are refused:
S9b.1's, S9c.1's, S9c.2a's, S9c.2b.1's,
S9d.1's, S9d.2a's, S9d.3a's, S9d.3b.1's and S9d.4b.2's `Degenerate` (a turned box's corner on
another's wall, edge or corner, or on a cylinder, within rounding; a wall
tangent to a cylinder; two cylinders touching at a point; equal cylinders
in a turned frame whose axes meet, their section within the resolution
of a node; a box's corners on a sphere, reported as a vertex of one input
on the other's face since S9d.2a, and its wall tangent to one; a
cylinder's rim on a sphere's equator, the wall tangent to the sphere
along it, its discriminant vanishing identically, a tangency since
S9d.2's survey; a box's wall through a frustum's axis, its corners on the
rim, and a frustum's base circle tangent to a face's edges; a cylinder's
cap through a frustum's virtual apex, a frustum's base rim on a
cylinder's rim, turned about the axis or not; two frusta on one cone or
with one virtual apex, coincident rims or bases on one plane; three
copies of a torus about perpendicular axes, their tubes touching), and an
arc ending off its circle; of the stacks
given to another Boolean, those with cylindrical or conical walls (S9c;
`UPSTREAM_TESTS.md`).
