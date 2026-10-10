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
(`curved/cones.rs`), S9d.3b.2's cones and spheres in loops, S9d.3c's
cones' loops against cylinders and cones in any frames
(`curved/cones_loops.rs`) with turned cones' loops against spheres and
turned caps against cones, S9d.4a's
whole tori against polyhedral prisms (`curved/torus.rs`), S9d.4b.1's
torus v-segments and wedges against them (`curved/torus_segment.rs`), and
S9d.4b.2a's whole tori against prisms with arcs, spheres and cones and
S9d.4b.2b's two whole tori (`curved/torus_curved.rs`), and S9d.4c's caps
and zones against tori and torus parts against curved solids
(`curved/torus_parts.rs`). Since S9e a Boolean's result is an input again
(below), and S9e.4a's imported solids (a `.brep` or STEP body without a
construction, `solid/imported.rs`) are decided on the construction their
stored surfaces give, S9e.4b.1's imported prisms with their arcs' ends
taken onto their circles (`curved/snapped.rs`), S9e.4b.2's imported
polyhedra on their stored vertices (`polyhedra/imported.rs`), S9e.4b.3a's
imported plane pieces of a sphere, cylinder or cone as their primitive
common their planes' half-spaces (`curved/pieces.rs`), and S9e.4b.3b's
split pieces (S8's `Clipped` and `Half`) against curved faces on the same
model (`curved/splits.rs`), S9e.4b.3c.1's two inputs on one sphere
in general position and sphere pieces with their rims split by stored
vertices (`curved/graph.rs`, `curved/pieces.rs`), and S9e.4b.3c.2's exact
incidences of two inputs on one sphere (a vertex, a circle or a line of
both, plane faces on one plane with overlapping edges: `curved/graph.rs`),
and S9e.4b.3c.3a's imported bodies of one curved face and planes that are
another Boolean of their primitive and the hull of their other planes (a
groove, a bite, a boss: `solid/imported.rs`'s `Form`, `curved/pieces.rs`),
and S9e.4b.3c.3b's that are a Boolean tree of their primitive and several
hulls (a groove or boss on a U, a primitive bitten twice, a sphere's disc
or a cylinder's flat with a box: `solid/imported.rs`'s `Tree`,
`curved/pieces.rs`), and S9e.4b.4a's imported prisms whose arcs of two
circles meet at a joint, each such arc taken through its two ends
(`curved/snapped.rs`), and S9e.4b.4b.1's imported bodies of several
sphere, cylinder and cone faces whose plane faces are their primitives'
ends, a Boolean chain of those primitives (a stepped shaft, a cup, a dome
and a pin on one ball: `solid/imported.rs`'s `primitives_piece`,
`curved/pieces.rs`), and S9e.4b.4b.2a's such bodies led by a prism leaf
(`bfuse_complex/K1`'s rounded box less a bore, a plate with a boss, a dimple,
a dome or a conical pocket: `solid/imported.rs`'s `leaves`), and
S9e.4b.4b.2b.1's with a primitive's flat, a prism of one cap or a pocket (a
ball's half with a frustum on its disc, a prism under or on another, a plate
less a slot holding a crescent: `solid/imported.rs`'s `bosses` and
`chain_piece`), and S9e.4b.4c.1's imported polyhedra against curved faces
and with cavities, their stored triangles a leaf of the curved engine
(`curved/meshes.rs`).

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
cylinders whose axes cross add their two crossing points, where the
surfaces are tangent: inside both faces with the faces' outward normals
opposite there (the inputs touching each other, a bore against a rod of its
radius across it) they are a tangency between the inputs, `Degenerate` for
every operation; with the normals alike (whole rods, overlapping there)
they are vertices. Two faces'
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
`Degenerate`, and so are two cylinders whose axes are within rounding of
parallel without being parallel (the sine of their angle at most
`10^-12`, as a frame's normal normalized again is, its meeting a sliver
of an ulp of the carrier's angle) unless they are certainly apart within
their faces' bounds: the other's section beyond the first circle, within
it or holding it by a certified margin over its axis's drift there. Orders along a curve and around a vertex are exact signs of
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
S9b.1 keeps fragments. The arrangement is the operation's own only in that
keeping (which side of a piece the other input holds is decided without
it), so the last two arrangements are kept by their inputs' content (the
inputs' and the seams' `Debug` text): fuse, cut and common of one pair
share one. Faces of both inputs on one surface (coplanar caps
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
several solids `OutOfDomain`), shells meeting at a vertex `Degenerate`,
and so is a shell touching itself at a vertex (the edges there, linked
where a loop runs from one to the next, in more than one fan; a hole's
wall and a rod of its radius crossing its axis, which first showed it, is
now refused before as a tangency between the inputs, below).
Names follow S9b.1's rules; a result is a `Polyhedron` (both inputs, the
operation, its index), classified by the set function and moved by its
stored geometry. A result with arcs given to another Boolean is decided on
its construction's exact model (S9e.1, below).

### Cylinders meeting in quartics (S9c.2a)

Two circular cylinders in exact frames (the world's axes permuted or
reversed) whose axes cross are perpendicular; unless they are S9c.1's equal
cylinders with meeting axes, they meet in a quartic (`curved/procedural.rs`).
Across the common perpendicular `e = nA x nB` each cylinder reaches
`[e_k - r_k, e_k + r_k]`; exactly: apart, touching (outside, or inside
with a node) `Degenerate`, one extent strictly inside the other's (two
rings about the inner cylinder), or overlapping in part (one loop); ends
of the extents near each other are a near node, `Degenerate` by S9c.2b.1's
margins (below). Each
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
exactly), is `Degenerate`. So is an extremum of either cylinder's `D`
nearer zero than the meeting's stored binary64 image holds, in exact
frames too: a stored meeting is the root `(-B + s sqrt(D)) / A` with its
coefficients rounded from the stored frames, `D` within about `eps A L^2`
(`L^2 = 2 rho^2 + r^2`, `rho` the reach from the carrier's stored origin
to the other's axis plus its radius, `r` the other's radius), so its
height near the extremum, and a vertex or a ring's closing point there,
is uncertain by about `eps L^2 / sqrt(D)`; `|D| < (8 eps L^2 / res)^2` is
refused (`turned::conditioned_node`, each chart at an axis point read
over the quarter turns either side of it).

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
axis's planes is `OutOfDomain`, and one other than a circle passing within
rounding of it, or so near it (within about `1e-4` of a ball of radius
`1.25`) that 256 anchors leave its projection's lift unpinned, is
`ComputationLimit` (a ball centred within rounding on a parallel
cylinder). A plane within the resolution of tangency to a sphere is
`Degenerate` ("a plane crossing a sphere within the resolution of
tangency"): crossing it wherever the faces' boxes meet (a cap no higher
than the resolution), and missing it where the plane's point nearest the
centre lies in the plane face, the sphere's point nearest the plane in the
sphere's, and the gap between them outside either input by the faces'
outward normals (the inputs touching across it, or one inside the other
behind a wall thinner than the resolution; their faces' boxes widened by
it), as are two spheres apart or nested by at most the resolution with
each one's point nearest the other in its face (`graph::near_miss`). A
gap inside both inputs (a slab's floor beneath a dimple's sphere) is no
contact, as equal rods overlapping at their tangent points are not. A
sphere whose distance from an input edge's line or circle (its nearest
point strictly inside the edge: a line's foot, a circle's point exact in
the quadratic field of the centre's projection, a conic that is no circle
exactly at a rational point of it at its least distance to rounding) or
from an input vertex lies within the resolution of its radius, the
sphere's point nearest it in its face, is `Degenerate` too ("a sphere
within the resolution of tangency to an edge", "a sphere within the
resolution of a vertex"):
crossing it wherever so, missing it where the gap lies outside either
input by each one's membership at its point pushed across it
(`graph::edge_near_misses`, after `near_miss`; edges between faces on one
surface and seams' edges and vertices are none; on a given result's
meetings of curved faces, cone, torus and spline curves the nearest point
at a rational parameter at each least distance found in binary64,
`near::nearest_on_run`). A sphere within the resolution of tangency to a
cylinder or a cone face ("a sphere within the resolution of tangency to a
cylinder or cone"), its point nearest the centre exact in the quadratic
field of the centre's distance from the axis (on the local axes, a turned
frame's within rounding) and in the face, the sphere's toward it in its
face, is `Degenerate` by the same rule, before the meetings are found
(`near::sphere_quadrics`); and so is an input edge's line, conic or circle,
or a given result's other curve, within the resolution of tangency to a
plane, cylinder or cone face of the other input ("an edge within the
resolution of tangency to a face"), at a point strictly inside the edge
where the surface's function along it is least or greatest (exact where
it is linear or quadratic along the curve, else a rational point at an
extremum found in binary64) and the surface's point nearest it in the face,
where both of the edge's faces leave the point away from the surface (a
face heading toward it crosses it there, and its section near the edge is
the arrangement's own: S9f.2b.2's `lens_tilt_loop`, a rod's top rim on a
lens's top plane but for rounding; a face running along it, its direction
into the face in the surface's tangent plane exactly or within the
`10^-12` the faces' own rules take as parallel, decided exactly, is a
tangency or incidence of the two faces, theirs: S9c.1's rod and a turned
box on its tangent plane; `near::edge_faces`). The
validator decides a ray against a whole sphere (a cavity in a sphere); a
result's bounds hold its spheres' boxes. A sphere against a cylinder or
another sphere is S9d.2's.

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
S9d.2c's (below), against a cone likewise (S9d.3c). A meeting
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
circle against a cone are S9d.3c's (below).

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
S9d.3c's (below).

### Cones' loops against cylinders and cones, turned cones (S9d.3c)

Two ruled faces (a cone and a cylinder, two cones) whose `D` changes sign
over both carriers' angles meet in loops, met as S9c.2b.1 meets two turned
cylinders, generalised to ruled carriers on any affine frames
(`curved/cones_loops.rs`): carrier 0 the cylinder (its `A` constant), else
the object's cone; each interval of `D_0 > 0` in a chart based at a
negative point is one component, its plus branch up and minus branch back;
carrier 1's turning points (`D_1`'s real roots, its ruling touching the
first quadric at `w = -B_1 / A_1`) are placed along it by binary64 views,
and a switch at a rational angle of carrier 0 between adjacent turning
points of different kinds (a point with one surd). The runs between
consecutive switches are graphs over carrier 0's angle about carrier 1's
turning points and over carrier 1's about carrier 0's, `Curve3::Meet` with
either carrier; each is verified exactly: no root of its carrier's `D` or
`A` strictly inside its range (Sturm counts at surd or algebraic ends, a
root at an end left out), its ends and, for carrier 1's, an exact interior
point of the run on its branch, the switches' order along the component
exact (the plus branch by ascending chart `t`, the minus by descending).
Where two cones' direction cones cross (`A`'s simple real roots) a
component runs through infinity: at a root of `A_0` the branch `-sign B_0`
runs off (a cut, no vertex) and the finite branch's point `-C / 2B` (in
`Q(alpha)`) is a switch, carrier 1's roots alike; a graph over carrier 1
reaching a cut ends at one of carrier 1's roots on its running-off branch.
A cylinder exactly along a cone's ruling (its `A` zero; S9d.3b.1's rounding
case in the limit) and a repeated real root of `A` are `Degenerate`, with
S9d.3b's tangencies, near nodes, a carrier's apex on the other quadric and
two turning points (or a turning point and a cut) within rounding along a
component. A turned cone's loops with a sphere take S9d.2c's height graph
with the cone's radius `rho(w) = r + k w` (`F`'s coefficients of degree two
in `w`; the apex's height, where the circle is a point and the
discriminant vanishes without a turning point, left out), `Curve3::Rise`
with the cone's half angle keeping its closed form; a turned cap's circle
against a cone S9d.2c's resultant with the cone's radius term, its
tangency band from that term's bound over the circle. A point at a cone
carrier's apex height is on none of its `Meet` pieces (the apex itself is
refused before).

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
tangent (a sphere's circle of a surd radius, S9d.4c's below); a seam's
tangency with a torus is retried at another seam. A v-segment or wedge
against a curved face is S9d.4c's (below).

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
refused (`tori_kiss`). Two tori of equal radii whose centres lie within
`10^-12` of their size of each other and whose axes are within rounding of
parallel (a frame's normal normalized again, or one normal with the axes
turned and rounded) are one surface within rounding, `Degenerate`, as one
surface exactly is (their meeting's projections were left unpinned,
`PrecisionLoss`, before). Parallel tori whose top or bottom circles lie at
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

### Caps and zones against tori, torus parts against curved solids (S9d.4c)

A sphere's own circle (a cap's or zone's rim of radius `sqrt(R^2 - h^2)`,
or a circle of a turned frame's cap on a basis of unequal lengths) meets a
torus where the torus's quartic on the circle's plane vanishes on it
(`curved/torus_parts.rs`, `circ_torus`): in coordinates `(s, t)` turned by
a rational rotation the circle is a quadratic in `t` of constant leading
coefficient, the quartic reduced modulo it `r1(s) t + r0(s)`, and their
resultant, of degree at most eight in `s`, has one simple real root per
crossing, `t = -r0 / r1` there (`Q(alpha)`); a repeated root or `r1`
vanishing tries the next of eight rotations, a resultant vanishing
identically is the circle on the torus, and a crossing within the
resolution of a tangency (the count changing with the torus's function
offset by `8 R r (R + r) res`, its gradient's bound on the torus times the
resolution) is `Degenerate`. A round circle of a rational scale keeps
S9d.4b.2a's octic. The cap's discs meet the torus in S9d.4a's spiric
sections, its wall in S9d.4b.2a's meeting with the whole sphere.

A torus v-segment or wedge (S9d.4b.1's model) against a prism with arcs, a
sphere, a cone or a whole torus: its wall's meeting is S9d.4b.2's with the
whole torus's surface (traced and verified over the whole torus, a
tangency anywhere on it `Degenerate`), its pieces kept where the part's
patches hold them; its end discs meet quadrics in plane sections and
another torus in S9d.4a's spiric sections. A rim against a curved surface
(`rim_far`) fixes one of the torus's angles at a quadratic surd (a
segment's `v` at `(s sqrt(r^2 - z^2), z) / r`, a wedge's end `u` at its
rounded direction over its length): the other surface's function `G` there
is `A + sqrt(d) B` in the free angle, and the rim's crossings are the real
roots of the norm `A^2 - d B^2` (degree eight, sixteen for two tori in
rounded frames) at which `A + sqrt(d) B` vanishes, decided exactly in
`Q(alpha)(sqrt(d))` (the conjugate's belong to the other ring or the
opposite half-plane's circle); a root where its derivative vanishes too is
a tangency. The halves' rims (circles of radius `R`) keep the conics'
crossings. Signs of surds over algebraic fields are tried by a binary64
enclosure before their exact products.

Two corrections to the arrangement's assembly, both reachable only here:
a shell of both inputs' faces enclosing a void (a band's end disc and inner
wall under another input's face across its hole, in a fuse) is found by its
certified flux built alone and made a cavity of the result (the
validator's rays then leave its containment undecided: `ComputationLimit`,
as S9d.4b.1's cavities), and a loop through a sphere's pole (a wedge's
rounded end half-plane within rounding of a sphere's axis, a meridian's
section through the pole) winds by its pcurves' lifted ends, the half turn
at the pole included.

### A Boolean's result given to another Boolean (S9e.1)

What a result stores is rounded: its faces keep their inputs' stored
surfaces (a wall's cylinder on a frame at its circle's rounded centre, not
the prism's affine cylinder), its vertices and edges' curves are rounded
once, and what the inputs' exact models decide (a profile's line tangent to
its arc, a shared cap, a section through a vertex) does not survive the
rounding. A result of the curved arrangement given to another Boolean is
therefore decided on its construction's exact model, its *given model*
(`curved/given.rs`): the arrangement of its inputs' exact models that built
it, run again (the seams tried in `build`'s order, the thread's cache of
arrangements), its kept pieces the body. The model's faces are the input
faces holding kept pieces, each on its input's exact surface in its input's
frame, reached through the input's model (`Prism::view`: the second
arrangement's cylinder pairs, sections, parameters and seams take a face's
input model and face; the faces of one given model lie in two frames); a
face's region is where its input's `in_face` holds and the other input's
sides at the point (pushed off the face both ways where it lies on a face
of the other on its surface) are those the first Boolean keeps, `On` on
the other's boundary elsewhere (the first arrangement's sections); its
orientation is the kept pieces' (a cut's tool's faces reversed). The
model's edges are the arrangement's edges on the kept pieces' boundaries
(exact lines and conics, each conic with the model face whose cylinder's
angle is its parameter), its vertices the arrangement's exact vertices;
its membership is the first Boolean's set function over both inputs'
exact membership, pushed alike (a boundary point `On` unless the other
input decides it). The stored topology only names the model: the re-run's
assembly gives the result's slots in the stored order (its provenance: each
face slot's pieces, each edge slot's arrangement edges, each vertex slot's
arrangement vertex), checked against the stored vertices (bit for bit, or
within the resolution after a rigid motion, whose moved inputs' models are
exact in the moved frames), and every model face, edge and vertex carries
the ids of the result's entities it lies in (none for an edge inside a
result face, a full circle's seam, or a vertex inside a result edge).

The second Boolean then arranges the given model as either input: its
edges are pierced by the other input's faces and the other's edges by its
faces (lines and conics against planes and cylinders, S9c.1's meetings: a
new vertex on a given section is three surfaces' point without a new
computation), its faces meet the other's in S9c.1's sections restricted to
both regions. Names follow S9a's rules over the given result's ids, so the
second history is over the first result's entities: a model face holding
several result faces (a slot through an input face left it in parts) names
each second-arrangement piece by the result face it lies in (its pieces
joined across the edges they share, each group by a model edge on its
boundary and the result face on its side); an edge along a given edge keeps
that edge's stored circle or ellipse frame, its new ends' angles on it. In
scope: a result of prisms of lines, arcs and circles in any frames (S9c.1's
pairs) that is the only solid of its Boolean, as object or tool, with a
prism of lines, arcs and circles or another such result (S9e.2 adds
stacks, S9b.1 results and one solid of several, below). Refused: a plane's
piece (S9e.4); results of spheres, cones and tori, with procedural
edges, against a sphere, cone or torus, and deeper chains are S9e.3a's
(below), a meeting of two curved faces met by the partner S9e.3b's;
a given edge through an irrational point against a cylinder
(`ComputationLimit`, as a section's line in S9c.1); every S9c.1
degeneracy, including where the other input meets a given model's section
at a point where the first result has no face.

Two corrections found here: parallel circular cylinders' relation holds the
second cylinder's circle in the first's frame, so an arc of the second
input's cylinder is met against the first's with the relation taken from
its own side (read in the first's frame its crossings were missed and the
result left open: a cylinder's cap circle across a parallel one's wall, in
exact frames), and the independent history check takes an ellipse split
into arcs (a given result's section edge cut by the second Boolean) as one
curve, as it takes a circle.

### Stacks, S9b.1 results and one solid of several given (S9e.2)

A stack (S9a.2) with a cylindrical wall, or of planes against a partner
with an arc or a cylindrical face, and an S9b.1 result of two line prisms
under the same condition, are decided on the same given model
(`curved/given.rs`): their construction is `A op B` of two prisms (a
stack's two profiles as prisms on its frame over their heights), whose
point set is the inputs' exact set function whichever arrangement built
it, so S9c.1's arrangement of the two prisms' exact models is run again
(same-frame prisms included: walls parallel, caps coplanar) and its kept
pieces are the body. Their stored slots come from another assembly
(stack.rs's slabs, polyhedra.rs's fragments), so the re-run's assembly is
matched to the stored topology geometrically (`curved/matched.rs`): each
stored vertex the one re-run vertex within the resolution of it, one to
one; each stored edge the re-run edge between the matched ends through
the stored edge's points at a quarter, a half and three quarters of its
range (a ring through the stored ring's curve); each stored face the
re-run face bounded by the matched edges on a surface of the same kind;
the given solid the one re-run solid so matched. Anything unmatched is
`ComputationLimit("a given result rebuilt differently")`. The model's
faces, edges and vertices carry the matched ids, its edges the matched
stored curves (a stored circle or ellipse turning against the conic's
parameter read the other way). With a polyhedral partner a stack of planes
or an S9b.1 result stays S9b.2's stored model.

One solid of a result of several (of S9e.1's class, a stack or an S9b.1
result) holds part of its construction's region only. Its given model is
the whole construction (every solid's faces and edges, the construction's
regions and membership), and once the second arrangement is built its
pieces are sorted by solid, exactly and by adjacency alone
(`matched::keep_solid`): a given face's pieces joined across the edges
they share take the solid of a given edge on the group; a piece of the
other input with a side inside the construction takes the solid of the
given pieces sharing an edge with it, else of its neighbours across the
other input's own edges inside the construction. Pieces of the other
solids' faces are dropped, and the other input's pieces inside another
solid are outside the given one on both sides. The other input strictly
inside one solid meets no face to go by: `ComputationLimit("a solid
inside one of a given result's several solids")`. The arrangement holding
every solid, a degeneracy between the other input and a solid not given
is S9c.1's `Degenerate` too (a conservative refusal).

Two corrections found here: pieces of both inputs on one surface joined
into one face (the other input's piece reversed beside a kept one, the
rollex's pocket floor cut by a cylinder standing on it) ran their loops
about their own faces' opposite normals and did not close; each piece's
loops now run about the face's normal, every vertex's way on chosen about
it. And a stack moved rigidly kept no height range (every point off its
heights), so a result with a moved stack among its inputs classified its
own vertices outside; the moved stack keeps its range.

### Given results of spheres, cones and tori, deeper chains (S9e.3a)

A Boolean's result whose construction holds spheres, caps, zones, cones,
frusta, tori or their parts, or given results themselves, is decided on the
same given model (`curved/given.rs`, `curved/chain.rs`): its construction
is re-run whichever arrangement decided it (S9c's or S9d's), its faces are
its leaves' faces holding kept pieces on their exact surfaces (a view goes
down every given level to the primitive model; a face's outward normal and
its parameters' sign compose each level's orientation, so a sphere's,
cone's or torus's face kept as a cut's tool is reversed there alone), and
every first-arrangement edge on a kept piece is a model edge with its exact
curve and places: lines, conics, a sphere's circles, a cone's sections (one
normal to the axis given as the circle it is), a torus's sections and
rims, and the meetings of two curved faces (`Meet`, `Rise`, `Toric`). The
construction tree is evaluated level by level, each level's arrangement
over its two leaves' models, to three Booleans (deeper,
`ComputationLimit`). A given result may meet a sphere, cone or torus: the
leaves' surfaces meet the partner's in S9d's pairs.

The second arrangement meets a given edge by S9c.1's and S9d's edge
meetings through the edge's own surface (a sphere's circle its sphere,
reached through its faces' views): lines and conics against every
surface, a sphere's circle against planes, cylinders, spheres, cones and
tori, a torus's section against planes and, where it is a circle of the
torus at a fixed angle, against quadrics and tori. One meeting is new: a
cone's section against a plane, on the line of the two planes
(`ConeSec::meet_plane`). A meeting of two curved faces, a cone's section
against a curved face or a torus's general section against a curved face,
met by the partner, is three surfaces two of them curved: S9e.3b's (below);
where one of the edge's faces' surfaces is apart from the partner's face the
edge is not met and stays given. A result edge over the
whole of a given edge keeps its stored curve; a piece of a cone's or a
torus's section is rounded from its exact curve on its primitive model.

In scope besides S9e.1's and S9e.2's: results of spheres, cones and tori
with prisms, spheres, cones and tori; results with procedural edges the
partner does not reach; deeper chains to three Booleans; given results
against spheres, cones and tori. Refused until S9e.3b: its meetings.
Refused: deeper trees; splines (S9f); every S9d degeneracy in the second
arrangement (a plane's piece, refused here until S9e.4, is given on since
S9e.4b.3b: below). DRAW's `bcut_simple/G9` and `H3` are such a degeneracy: the
rod of radius 1 about `(5, 0)` touches the frustum's top circle of radius 6
at `(6, 0, 4)`, a tangency between the inputs.

Corrections found here: S9d.4c's meeting of a part's rim with a curved
face assumes a circle at a fixed angle of the torus, and a given result's
general section was met by it and the result left open (the fuzz replay's
sphere partner): such a section met by a curved face is now S9e.3b's; and
the independent history check takes a torus's plane section, a hyperbola
and a parabola split by another Boolean as one curve.

### Given meetings of curved faces met by the partner (S9e.3b)

A given result's edge on a meeting of two curved faces (`Meet`, `Rise`,
`Toric`), or on a plane's section of a cone or a torus other than a circle
of the torus, is met by the partner's faces where three surfaces meet
(`curved/triple.rs`). With a plane among the three, the other two are
restricted to the plane's rational affine coordinates (a quadric's conic, a
torus's spiric quartic); with three quadrics, the two others are restricted
to a ruled one's rulings over its rational chart (the given curve's
carrier, else the partner's). The resultant in the second parameter (degree
4 for two conics, 8 for a conic and a quartic or for three quadrics) has
its real roots isolated exactly; at each the second parameter is the
fibre's gcd over `Q(alpha)`, so every vertex lies in one `Q(alpha)` with no
tower. A fibre holding two points (a plane along a cylinder's rulings, a
partner parallel to the carrier) retries the projection (the parameters
exchanged, sheared, another ruled carrier). Each point is verified on all
three surfaces exactly and kept where the given curve's own test holds it
(branch, window, range); the second arrangement then takes it as any
vertex of an edge met by a face. A result edge over part of a given
meeting keeps its stored curve between its ends' parameters.

`Degenerate`: the partner tangent to the given curve at a meeting (the
three gradients dependent there, exactly), two meetings of one edge within
the resolution, and S9's rules. Refused by cost: a torus among three
curved surfaces (a `Toric` met by a curved face, a `Meet` or `Rise` met by
a torus) and two tori in one plane, eliminants of degree 16 or more
(`OutOfDomain`). Found with it: a torus face wound in `v` (a plane's cut
across the tube) with a hole is now decided by the `+u` ray's signed
crossings (before, `UncertifiedContainment` always), a result's bounds take
its torus faces' whole box and its cylinders' and cones' circles between
their edges' heights (a given result's faces bulge past its edges, and the
next Boolean's vertices classify against them), a vertex at a pole of a
sphere face it bounds stays in the result (a meridian's pcurve is undefined
across it), and the history check takes a split meeting's pieces as one
curve. An eliminant whose coefficient bound passes 4,096 bits is a
`ComputationLimit` (torus sections in turned frames against a turned
cylinder: 8,000 to 11,000 bits, minutes an operation).

### Imported solids (S9e.4a)

A body without a construction (a `.brep` or STEP converter's cell topology
and resolution) is a solid through `Solid::imported_with(operation,
topology, resolution)` (`solid/imported.rs`): its entities renamed under the
operation (each slot its ordinal, so two imports never share ids), each
generated in the history from a label of its slot. What it stores is
rounded (OCCT writes surfaces and curves with 17 significant digits,
vertices with 15; the converter normalizes every frame again), so it is
decided on the construction its stored surfaces give, read off them once: a
prism of lines, arcs and circles (two planar caps facing apart, every other
face a plane or cylinder along their normal; the bottom cap's stored frame,
turned exactly to point into the material, its loops' vertices and arc
centres in that frame's exact affine map rounded once, the stored radii,
the top cap's height rounded once); a sphere, cap or zone (the stored frame
and radius, each disc's latitude from its ring's height); a cone or frustum
(the stored frame moved along its axis to the lower end, the rings' radii,
zero at the apex); a whole torus (the stored frame and radii). The
construction's topology is matched to the stored one by S9e.2's geometric
match (vertices within the resolution, edges through their points, faces by
their edges and surface kinds), so with the converter's own validation
every stored surface lies within the resolution of the construction's. In a
Boolean the construction stands in the imported solid's place
(`polyhedra::build`), every pair's engine and rule its own, and the
result's plans are renamed through the match: the history is over the
stored ids. The imported solid classifies by its construction and moves
with it (its stored geometry moved, its construction rebuilt in the moved
frame); its mass is the stored topology's certified enclosure. A result of
an imported solid given to another Boolean re-runs with the same
construction. No stored edge is trusted as an exact curve, so a line
tangent to its arc and a periodic face split at a seam (two arcs of one
circle) are the profile's data, decided as S9c decides them.

An imported prism's arc whose ends, rounded into its cap's frame, lie off
its circle (a cap turned or tilted from the world's axes, a vertex at an
angle whose cosine is irrational: OCCT's 15 digits) is S9e.4b.1's: the
construction keeps the rounded profile (its topology and the match are
S9e.4a's), and its exact model (`model::Prism::new`, for a profile the
import flags, `Profile::rounded_arcs`) takes each arc's end onto the
arc's circle (`curved/snapped.rs`): the circle's rational point at the
binary64 rounding of the end's half-angle tangent `dy / (r + |dx|)`,
computed exactly, reflected where the end lies left of the centre; an end
on its circle is kept, a joint of a line and an arc takes the arc's point
(the line's end with it), a joint of two arcs of one circle that circle's
point. An end moves by at most its distance from the circle plus `r
2^-52`, within the resolution (the profile's validation), so the model's
vertices stay within the resolution of the stored ones. A joint tangent
in the body OCCT was given is tangent in no binary64 data; taken onto the
circle, the line meets it at the joint at an angle within rounding of
tangency, which the arrangement decides by its exact turn as any joint. A
kernel profile's arc ending off its circle stays refused (S9c: the caller
placed it).

Where arcs of two different circles meet at a joint whose rounded point
lies off either (a crossing joint: four discs' common, a fillet chain, an
arc tangent inside another, circles crossing at a small angle), their
common point is a quadratic surd, which a profile segment's rational ends
do not hold: S9e.4b.4a takes each arc ending there whose circle no other
arc of its path shares through its two ends instead (`snapped::path`). Its
circle in the exact model is the one through both ends `a` (its start) and
`b`, of centre `a + rho e` and radius `rho`, `e` the rational unit vector
nearest the stored centre's direction from `a` (the same half-angle
tangent, rounded once) and `rho = |b - a|^2 / (2 (b - a) . e)` exactly;
its centre and radius must lie within the resolution of the stored ones
(else `Degenerate("an imported prism's arc too short to take through its
ends")`). The joint is its rounded point, or the point of an arc's circle
kept (one several arcs share, or one meeting only lines). A joint tangent in
the body OCCT was given is decided by its exact turn, as S9e.4b.1's. A
partner's cylinder on such an arc's stored circle lies within rounding of
the circle the arc is taken through: two parallel circular cylinders
within the resolution of one and not one are
`Degenerate("two cylinders within the resolution of one cylinder")`
(`curved/meet.rs`'s `cyl_pair`), the sliver between their walls refused as
S9a refuses two boundaries within the resolution in one frame. A joint of
two circles each holding other arcs of the path stays refused
(`OutOfDomain`, S9e.4b.4).

An imported body of plane faces and line edges that is no S9e.4a prism (its
recognition, construction or match failing: a pyramid, a frustum of one, a
wedge with slanted faces, a result of boxes in different frames) is
S9e.4b.2's: a polyhedron decided on its stored vertices
(`polyhedra/imported.rs`, S9b.2's stored model in `polyhedra.rs`): each
vertex its stored binary64 point, each edge the segment between its stored
vertices, each face the polygon of its stored vertices cut into exactly
planar triangles in its projection (ears clipped, holes bridged), its
outward normal from the region behind it, membership by an exact ray's
parity. Its planes' common points are not taken: a vertex of four or more
faces has none once they are rounded, and two bodies sharing a face would
meet in slivers where their stored points agree (the DRAW survey's
`buc60803a` and `b`: their planes' points 5.3e-15 to 7.2e-15 apart). A face
whose stored vertices are coplanar exactly is that plane, its fragments
joined with any face's on it; another is joined by the face (S9b.2's
rule). Its entities are its stored ones, so the history is over them
directly; it classifies `Boundary` within the resolution of a face's
triangle, else by parity, and moves with its stored topology. A face its
stored vertices fold is `Degenerate`. An imported prism of lines against an
imported polyhedron is decided on its stored vertices too, not S9e.4a's
construction, whose corners re-derived from rounded local coordinates can
miss the vertices two files share by an ulp (`bopfuse_complex/K5`).
S9e.4b.4c.1: against a solid with curved faces or edges, and in a result of
one given to a Boolean of curved faces, the polyhedron is a leaf of the
curved engine (`curved/meshes.rs`): each triangle of its stored model a
model face on its exact plane, named by its stored face (so a result face
of one stored face joins its triangles' pieces across their diagonals), each
triangle's edge a line between stored vertices named by the stored edge with
those ends (a diagonal by none), its membership the parity of an exact ray's
crossings, a point on its surface pushed along directions decided by the
wedges of the triangles holding it. A body of one solid region of several
shells (a cavity) is such a polyhedron too, in both engines; each cavity of
a result is a shell and a void region of its own, and where an input holds
a cavity the curved assembly tries every shell's orientation.

An imported body of one sphere, cylinder or cone face and plane faces that
is none of S9e.4a's constructions is S9e.4b.3a's plane piece
(`curved/pieces.rs`): its primitive (the whole sphere on the stored frame
and radius, or a cylinder or a cone on the curved face's stored frame
reaching past the body's ends) common the half-spaces of its plane faces'
stored planes (`o + u x + v y`, normal `x * y`, the side by the face's
region), decided as the given model (S9e.1) of that first Boolean. Its
second input is a hull leaf model (`Hull`): the planes' convex polygons of
three planes' points (rationals, bounded by a cube about the primitive),
lines between them, membership every half-space. The first arrangement's
assembly is matched to the stored topology on import (S9e.2's match), so
the model's faces, edges and vertices carry the stored ids and the history
is over them directly; it classifies by the primitive and the planes'
sides within the resolution and moves with its stored topology, read off it
again. Two rules of the engine were widened by it: a section within the
resolution of a stored sphere's pole takes a vertex of both faces already
within the resolution of the pole as its pole vertex (two planes' line
through a turned frame's pole, rounded); and a loop through the pole its
face would close at, whose pcurves turn half a turn there (a meridian
circle through both poles, or a wall through the axis), winds none, its
pcurves from that pole on lifted by the turn, the face closing on it (no
pole vertex loop beside it). A given
model's circle matched to a stored circle whose frame turns against it is
read the other way (`Given::flip`, S9e.2's rule, for spheres' circles too).
Only spheres' pieces come from `.brep` files today: a cylinder's or cone's
oblique section is an ellipse record the reader does not read, and a
sphere's section other than a meridian or a parallel of its stored frame
carries an approximated pcurve the converter does not certify; the
kernel's own split pieces' topologies are pieces too.

The kernel's own plane pieces, S8's split pieces (`Clipped`: a prism of
lines, arcs and circles split by a plane oblique to its axis; `Half`: a
cone's, frustum's, zone's or cap's piece by a plane through its axis or
across it, a whole torus's band by a plane normal to its axis), are
S9e.4b.3b's (`curved/splits.rs`) where either input has a curved face: the
same model read off the split itself, its primitive the split's own (the
prism on its frame between its heights, the cone, the whole torus; a
zone's or cap's whole sphere with its ends' parallels' planes) common the
half-space of the plane the piece was built on, taken exactly into the
world through the frame's axes (a plane within the resolution of the axis
through it, a torus's band's at its rounded height), matched to the piece's
stored topology so the history is over its ids; a profile not convex
across the plane leaves several pieces on a side, each one solid of its
construction (S9e.2's sorting by solid). A `Clipped` prism of lines against
a polyhedral partner stays S9b.2's stored model. A cone's half by a plane
through its axis is `Degenerate` (S9's apex rule); a torus's spiric piece
(S8d.3) and a spline prism's piece are refused (`OutOfDomain`). Three rules
of the assembly were widened by it, each an imported piece's too: a loop
through a stored sphere's pole that winds none has its fins from the pole
on lifted by the turn it made and the face's holes lifted with it; a stored
circle matched as a ring runs the result's way (turned over where it runs
the other); a cone's section stored as a ring likewise about the cone's
axis.

Two inputs whose spheres are one (equal centres and radii as rationals: a
sphere, cap or zone, an imported piece's primitive, a split zone) are
S9e.4b.3c.1's: their sphere faces are faces on one surface, as equal
cylinders' and coplanar planes' are (each holds the other's edges inside
it, its pieces classified by the other's membership pushed off the sphere
both ways, the kept pieces of both facing one way joined); where a circle
of each input on the sphere crosses one of the other's, the first's point
in the second's plane is a vertex of both (`VKey::Circles`, a whole
sphere's split great circle among them, which bounds no plane face), and a
piece's sphere is split at the second arrangement's seam so the two
pieces' great circles differ and move with the seams tried. An imported
sphere piece whose rim OCCT stored as two arcs (the survey's `so1` and
`so4`) has its own arrangement's ring split at the circle's point in the
direction of each stored vertex there (`VKey::Stored`, a vertex the
assembly keeps), so S9e.2's match takes the stored arcs. Spheres within the
resolution of one and not one, and circles of the inputs tangent on the
sphere, are `Degenerate`; concentric spheres of different radii stay
`Degenerate`. The exact incidences of two inputs on one sphere are
S9e.4b.3c.2's, merged in the arrangement: model vertices of both at one
exact point one vertex (`VKey::Both`, continuing both), a model vertex of
either strictly inside a line or circle edge of the other splitting it
(and lying on its faces), and, each input's edges split at every vertex of
both on them, a part of B's edge with the ends, the line or circle and the
arc of a part of A's one arrangement edge, its half-edges in B's faces too
(`Arr::shared`, continuing both edges); a face on one surface with the
other's has it as its boundary already, a section of two faces along it
(two wedges' half-planes crossing along their axis) is taken by it, and an
edge of one along the other's face on one line or circle with an edge of
that face is an edge of both where they overlap and outside the face
elsewhere (else `Degenerate`). The splits' own great circles take part
alike (a split's vertex on the other's edge is no longer a seam's
conflict). Incidences within the resolution and not exact stay S9's
`Degenerate` (the survey's `so6` and `so7`, any exact incidence turned by
a rotation that rounds its frames), as do a result touching itself along
an edge of both and faces meeting only along a line (S9's near-plane
guard: two octants about one axis).
A body of one curved face and planes that is no common of its primitive
and its planes is S9e.4b.3c.3a's: one Boolean of its primitive and the
convex hull of its other planes (`imported::Form`), the hull less the
primitive where its curved face's material lies outside its quadric (a
hole, groove, slot or dimple; the hull the first input), else the common,
the primitive less the hull of its planes turned over (a bite) and their
fuse (a boss) in turn, the first whose given model matches the stored
topology the body. A cylinder's or cone's primitive spans its curved face's
own axial range, past each end by a quarter of it but ending at a cap (a
plane face normal to its axis at an end of that range, its outward normal
away from the range, or into it for the hull less the primitive), the caps'
planes not the hull's; two faces on one plane facing one way are one plane
of the hull (facing apart `Degenerate`); a stored vertex loop at a pole of
its sphere face (OCCT writes one at every pole inside a face) the re-run
lacks stays unmatched, deleted in a Boolean's history. A tangency in the
piece's own arrangement is its curved face's with its own planes:
`Degenerate("an imported plane piece whose curved face is tangent to its
plane faces")` (the DRAW survey's `bcut_complex/I6` tool). A body no form
matches is S9e.4b.3c.3b's: a Boolean tree of its primitive and the convex
hulls of its planes (`imported::Tree`), read off how its edges bend (convex
or concave at their middle points). Its plane faces fall in two groups: the
primitive's (sharing an edge with the curved face convex inside the
quadric, concave outside) and the other's, every other face joining
through edges that keep a group (convex, but concave for the primitive's
group outside the quadric), a component met only across the other kind
tried in the other group and in the same, in turn. A group's region is the
hull of its faces' planes less its pockets (faces joined by concave edges
within it, their planes turned over); the body the primitive common its
group's region fused with the other's (or, where no face of the other
meets the curved face, the primitive common both regions' union), or
outside the quadric the other's less that; inner Booleans' results their
given models, the root's matched to the stored topology. Planes of faces
within the resolution of one plane facing one way are one plane; an edge of
the curved face along which a plane face is tangent to it is the body's own
`Degenerate`. Bodies no tree of at most four Booleans matches (a pocket
within a pocket) are `OutOfDomain("an imported plane piece other than a
Boolean tree of its primitive and its planes' hulls (S9e.4b.4)")`. A given
piece's edge along its cylinder's section by a plane along its axis (a
slot's rim, its base a quadratic surd) meets another cylinder only where it is apart from it (no
real root); elsewhere it stays `ComputationLimit`.

A body of several sphere, cylinder and cone faces whose plane faces are all
ends of their primitives is S9e.4b.4b.1's: a Boolean chain of its
primitives (`imported.rs`'s `primitives_piece`: `bfuse_complex/E5`'s
stepped shaft, `bcut_complex/G9`'s dome and pin on one ball). The curved
faces on one stored surface are one primitive's: a whole ball on its stored
frame, or a cylinder or a cone over its faces' axial range, past each end
by a quarter of it but at a cap (a plane face normal to its axis at that end
facing away from the range, or into it where the primitive's material lies
outside its quadric: a blind bore's floor), and past the other primitives'
bounds too where it is the first, cut or in common; a cylinder or a cone
whose stored axis lies within the resolution of an earlier primitive's is
built on that one's frame, so coaxial primitives share one exact axis and
their caps on one plane one height. Widest first (a sphere's radius, a
cylinder's, a cone's widest end), each next one is cut where its material
lies outside its quadric, in common where its faces meet the earlier ones'
curved faces along convex edges, fused where along concave edges or none;
each Boolean is S9e.4b.3c.3b's tree's, the last one's given model matched to
the stored topology. The arrangement splits a cylinder's or a cone's meeting
with a sphere at a stored vertex within the resolution of it as a sphere's
circles (`Arr::split_at`, at the carrier's rational unit direction nearest
the vertex's), turns it with OCCT's stored circle the way its carrier's
angle does (`flip`, a ring by `ring_about`), and takes a given meeting of
two curved faces against a face of the partner on one surface with one of
its faces as lying on it (no three surfaces' points there). S9e.4b.4b.2a:
where some plane faces are no primitive's ends, or its faces are tangent
along an edge, the chain is led by a *prism leaf*: two plane
faces facing apart, the bottom's loops each of whose edges' other faces is a
wall along its normal (a plane parallel to it, a cylinder along it) sharing
an edge with the top cap, its outer loop among them; S9e.4a's prism read off
those loops on the bottom's stored frame (its arcs onto their circles, its
walls' tangent joints the construction's own), its faces the caps, the walls
and every face on their surfaces (`imported.rs`'s `leaves`, each tried in
turn, the match the arbiter); the other surfaces' primitives follow as
above, one along the leaf's axis on its axes bit for bit at its own stored
origin (`bfuse_complex/K1`'s rounded box less its bore, a stadium plate with
a boss, a dimple, a dome, a conical pocket). S9e.4b.4b.2b.1: then each such
chain with the other plane faces (no end, no leaf's): a *prism of one cap*, a
plane face whose outer loop's every edge's other face is a wall along its
normal (a prism boss read off its one cap, its other end hidden in the part
it stands on), S9e.4a's prism on that loop past its walls' far end by a
quarter, built on the side its stored walls' axes point to and, parallel
to the leaf, on its axes bit for bit, fused after the leaf; the remaining
plane faces in components through their edges with each other, one convex
within itself a *flat* of the one part it meets along convex edges (that
part common the hull of its planes: a ball's half under a frustum), one
concave within itself a *pocket* (the hull of its planes turned over, less
the primitives it meets along concave edges alone, its teeth) cut from the
chain (a slot across a plate leaving a crescent of a cylinder), each plane
of a hull parallel to an axis of the chain on that axis's axes
(`imported.rs`'s `bosses`, `plane_along`), the match the arbiter. A body's
own faces tangent along an edge other than a prism's walls' joints (two
walls along one axis tangent along a line along it; two faces of one stored
surface, a face OCCT split at its seam, none) are `Degenerate("an imported
body of several primitives whose faces are tangent along an edge")`; plane
faces other than
the primitives' ends, a prism's, a flat's or a pocket's `OutOfDomain("an
imported body of several primitives with plane faces other than their ends,
a prism's, a flat or a pocket (S9e.4b.4b.2b.2)")` (a hull across several
parts, a polyhedral boss, a prism cut at both caps' rims, a pocket holding a
prism or meeting a part both ways, two leaves of two caps); faces of one
surface on different sides, a widest primitive whose material lies outside
it, a primitive meeting the earlier ones along both kinds of edge and a
chain the match rejects `OutOfDomain("an imported body of several primitives other than a Boolean
chain of them (S9e.4b.4c)")`.

Refused: a joint of two arcs whose circles each hold other arcs of the
path, off either (`OutOfDomain`, S9e.4b.4: neither circle can move); a
polyhedron of several solids (S9e.4b), a cavity among several solids (S9c,
by design), a kept cavity whose containment the validator's rays leave
undecided (`ComputationLimit`, as the kernel's own); a plane piece no Boolean tree of its primitive
and its planes' hulls matches (a pocket within a pocket: `OutOfDomain`,
S9e.4b.4); a piece whose planes pass through its cone's
apex (`Degenerate`, S9's rule: `shading_132`); a body of several
primitives with plane faces other than their ends, a prism's, a flat's or
a pocket's (`OutOfDomain`, S9e.4b.4b.2b.2: a hull across several parts, a
polyhedral boss, a prism cut at both caps' rims, a mixed pocket, two
leaves of two caps) or no Boolean chain of them (`OutOfDomain`,
S9e.4b.4c); a
torus's v-segment or wedge and every other body of curved faces, a torus
among several among them (`OutOfDomain`, S9e.4b: a general body on its
stored surfaces); spline faces or edges (S9f); S9's
degeneracies unchanged (an
imported cylinder tangent to the partner's wall, a box on the plane of an
imported prism's wall, a corner of one input on the other's face, a slab
within the resolution of an imported polyhedron's face). OCCT 8.1 writes `.brep` version 3 under a copyright line
the reader refuses, so the fixtures are version 1. Found with it: a point
on a profile arc's chord is displaced alike for every chord (two arcs of
one circle share their chord, run either way, and the point counted inside
both circular segments, so a box crossing such a prism in the same axes
was left open, latent since S9c.1).

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

### Spline walls in any position (S9f.1)

A prism whose profile holds S8b's spline segments (with lines and arcs)
meets a prism of lines in any relative position, a same-axis pair whose
offset or heights round included, in the curved engine
(`curved/spline_walls.rs`): a spline segment is its exact Bézier arcs run
as the profile runs (`Seg::Spline`), its wall `o + S_x(tau) x + S_y(tau) y
+ w n` a face (`Surf::Spline`), its cap edges curves over it at a constant
height (`Crv::Spline`). A plane `m . (P - p) = 0` meets a wall where
`a(tau) + b w = 0` (`a` affine in `S`, `b = m . n`): for `b != 0` in the
crease `w = -a(tau) / b` (`Crv::Spline` with that height, its rounded
curve the affine image of the profile's piece, S8b.3's `plane_image`, and
its pcurve on the wall S8b.3's `wall_pcurve`), for `b = 0` in generatrices
at the roots of `a` on each arc (degree `p`); a line meets a wall at the
roots of its trace's equation on each arc, a cap edge or a crease meets a
plane at the roots of the plane's function along it. Every root is the
arc's own parameter, rational or the generator `alpha` of `Q(alpha)` (S9c.2b.2's
fields, `num.rs`), one generator for every root of one polynomial
wherever it is found (its primitive coefficients and root index key it), so
every vertex on a wall is a point of `Q(alpha)`, degree at most `p`, at a
known parameter. Membership is exact: a point of a spline segment (found
again from its coordinates at its generator, or for a rational point as the
common roots of `S_x - u` and `S_y - v`) is on the profile's boundary, its
side from the run's tangent there; any other rational point counts the
`+u` ray's crossings with each arc at the roots of `S_y - v`, by the
chords' predicate (an end strictly above the point), right of it by the
exact sign of `S_x - u` there. A spline joint's vertical edge (between two
spline walls) lying in a plane of the other prism, across the plane's
face, is part of that face when the segments arriving and leaving lie on
the plane's two sides there (the plane crosses the solid along it), a
contact otherwise.

`Degenerate`: a plane tangent to a wall along a generatrix (a root of
multiplicity above one), at a knot (above one on either side), or along a
joint's edge; a line or a curve over a spline tangent to the other's face;
a plane within rounding of a wall's axis (the crease's slope past `10^12`,
S9c.1's (c) on the exact models: `TILT`'s axis against `TILTX`'s caps is
`-8.9e-17` off, against `SIDE`'s and `XY`'s `x` planes exactly parallel,
generatrices). An edge of one input in the plane of a face of the other,
even off that face, is refused as S9c.1 refuses it. Refused, with S9f's
labels: a spline prism against a prism with arcs or another spline prism in
any position (S9f.2), a sphere or a cone (S9f.3), a torus, a given result
or a plane's piece (S9f), and a Boolean's result with spline walls given
to another Boolean (S9f).

R4 under rounding (found by the evidence before any Boolean): an interior
knot of multiplicity `p` where the profile is exactly C1 is removed once,
exactly, before every lift and placement of a spline (`topology::c1_reduced`:
the profile's own poles less the knot's, C1 by construction), in prisms'
walls and cap edges (`lifted_spline`, `spline_wall`), in creases and
pieces (`placed_spline`, `spline::piece`) and in a general body's rigid
motion where the motion's rounding breaks C1; a lift still not C1 is
`PrecisionLoss`. The validator's exact Green path and its certified
quadrature take a spline pcurve piece whose control polygon leaves its
patch: by the surface's domain while the curve keeps to it (a crease
nearly touching a cap), boxed by that bound; across a knot line by a sliver
at most `2^-20` of the patch (a crease's rounded identity in `u` a step
past the line), integrated on the patch's polynomial with the slivers'
error bounded and added (`BREP_VALIDATION.md`).

### Spline walls against parallel curved walls (S9f.2a)

A spline prism meets a prism with arcs, circles or splines whose stored
normal is exactly parallel to its own (bitwise equal or opposite) in any
other frame, turned about the axis or offset by an amount that rounds
(`curved/spline_parallel.rs`; S9a.2 keeps one exact frame). The other
frame's `(u', v')` are rational affine functions of this one's `(u, v)`
(`map2`), so every wall against wall section is generatrices over the 2D
crossings of the profiles' curves: a spline arc against a cylinder at the
roots of `|m(S(s)) - c|^2 - r^2` (degree `2 p`); against another spline
arc at the roots of `f(m(S(s)))`, `f` the other arc's implicit equation in
its own frame (the Sylvester resultant of `B_x(sigma) - X` and `B_y(sigma)
- Y`, interpolated exactly on an integer grid, once per arc), degree at
most `p q`, the other's parameter at each root `-b / a` of their first
subresultant (an element of the root's field), the crossing kept when it
lies in `[0, 1]`. Every root is the spline arc's own parameter, S9f.1's
generator: against a cylinder the spline's, between two splines always the
object's, so a crossing found from either input's edges or faces is one
number. Cap edges and creases over a spline meet the other's curved wall at
the same roots, a cylinder's cap edge (a conic) the spline's wall at its
angle there. A point of another field is found on a spline segment by the
arc's implicit equation's exact sign and its inversion; an irrational
point off every spline segment is classified at a rational point of a box
about it that no element of the profile meets (lines by sides, circles by
distance, spline arcs by their control boxes under exact subdivision).

`Degenerate`: walls tangent along a generatrix (a root of multiplicity
above one), a crossing at a knot of either curve where the legs beside it
do not lie on the other's two sides, at knots of both. `ComputationLimit`:
arcs whose degrees' product exceeds 16, a crossing at a node or cusp of the
other arc's curve. `OutOfDomain`: coincident spline walls (a resultant
identically zero on an arc, refused, S9f); crossing axes (against a
cylinder S9f.2b's, against a spline wall refused, S9f).

### Spline walls against crossing cylinders (S9f.2b)

A spline prism meets a prism with arcs or circles whose stored normal
crosses its own (`n_A x n_B != 0` exactly), a leaning, tilted or side
rod, ring or stadium (`curved/spline_crossing.rs`). Along the spline
wall's ruling at the run parameter `tau` the cylinder's function in its
own frame's exact rows is `F = A w^2 + 2 B(tau) w + C(tau)` (`A` a
positive constant, `B` of degree `p`, `C` of `2 p` on each Bézier arc),
so the section is each branch `w = (-B +- sqrt(D)) / A` over each maximal
range of the segment's run where `D = B^2 - A C > 0` (`Crv::WallMeet`, a
graph over `tau`, placed by `tau`): a turning point (a root of `D`) ends a
range where it lies outside either face, and the pieces beside it are
outside too. A point is on a branch when its profile point lies on the
segment within the range, it lies on the cylinder and `A w + B` has the
branch's sign, exactly; a piece's midpoint at a rational `tau` lies in
`Q(sqrt(D(tau)))`. Vertices: a cap edge or a crease over the segment meets
the cylinder at the roots of `F` along it (degree `2 p`, `wallcrv_cyl`); a
cylinder's cap circle meets the wall at the roots of `F` on its cap plane's
crease, at its angle there (`conic_wall`); vertical edges by S9f.1's
`line_wall` and quadratic surds. Every root is the arc's own parameter, so
every vertex lies in `Q(alpha)` of degree at most `2 p` or in
`Q(sqrt(d))`.

S9f.2b.2 adds the loops and the towers. A turning point strictly inside
both faces gets a graph over the spline prism's height about it: the run
parameter the one root of `F(., w)` in a rational window of one Bézier arc,
over the heights `[w_-(tau_s), w_+(tau_s)]` of a rational switch `tau_s` on
the side where `D > 0` (a quarter of the way to where the gentler branch's
slope over the profile has fallen to one), placed by the height; the
graphs over `tau` on both branches end at `tau_s`, and the two switch
points are vertices of the arrangement (`CylPair::Mixed`'s switches: the
meeting of a spline wall with a crossing cylinder is now found once per
pair of faces in the pairs' pass). Each piece is verified exactly before it
is kept: `F`'s roots at the window's near end outside the range and `D < 0`
at its far end, one root at the rational height `-B(tau_s) / A`, and `H = A
C'^2 - 4 B B' C' + 4 C B'^2` without a root in the closed window (no double
root at any height), the distances halved otherwise (at most twenty times,
then `ComputationLimit`). A cylinder's cap circle in a plane holding the
wall's axis direction meets the wall where each arc's implicit equation
vanishes along the circle's half-angle tangent `t` (degree at most `2 p`),
its points and angles in one field `Q(t)` instead of the tower
`Q(alpha)(sqrt(delta))` (`tower_points`). In the topology the graph over
the height is `Curve3::WallMeet` with a `window` of the wall's `u`
(`TOPOLOGY_MODEL.md`), certified by interval Newton and the implicit
function theorem on the window's span (`MATHEMATICS.md`).

In the topology the meeting is a new procedural curve,
`Curve3::WallMeet` (`TOPOLOGY_MODEL.md`): the wall face's stored surface,
the other cylinder's stored frame and radius, the branch's sign and the
`u` range; on its own wall its pcurve reads its own `(u, v)`, on the
cylinder it is the cylinder's inverse. Its certified jets come from one
knot span's polynomial (across a knot the union of both spans' jets to
the order the wall's continuity allows), the integrals along it are split
at the knots exactly, and `a`, `b`, `c` and the discriminant are exact
Bernstein polynomials per wall and cylinder (`MATHEMATICS.md`).

`Degenerate`: the cylinder tangent to the wall (a root of `D` of
multiplicity above one near the faces); a turning point at an interior
knot; a turning point on a face's boundary, at a segment's end, or outside
the faces within the resolution of both (a rounding away from turning back
on a cap's edge, or a loop's on a cap's rim); a vertex's polynomial with a
multiple root (an edge tangent to the other's face; a cap circle tangent to
the wall's generatrix, the meeting turning back on the rim); a cylinder
whose axis is within rounding of the wall's without being parallel (the
sine of their angle at most `10^-12`, as a frame's normal normalized again
is: `A` within rounding of zero, the meeting within the faces a sliver of
the run no binary64 edge holds), as two such cylinders are, unless they
are certainly apart within their faces' bounds: on every Bézier arc the
wall's points at the overlap's middle height lie beyond the cylinder's
radius, or within it, by more than the ruling's drift over half the
overlap's heights (`|X|^2 - (r +- m)^2` of one sign on the arc, exactly).
`ComputationLimit`: a graph over the height whose window does not verify
after twenty halvings. Spline walls against spline walls on crossing axes
stay refused (S9f).

### Spline walls against spheres (S9f.3a)

A spline prism meets a sphere, a cap or a zone (`Solid::sphere_with`) in
any position, either the object. Along the wall's ruling at the run
parameter `tau` the sphere's function `|X - c|^2 - r^2` (its stored centre
and radius, the world's three rows: `procedural::other_sphere`) is `A w^2
+ 2 B(tau) w + C(tau)` with `A = n . n` positive for every ruling, `B` of
degree `p`, `C` of degree `2 p` and `D = B^2 - A C = A r^2 - |(P - c) x
n|^2` of degree `2 p`: S9f.2b's quadratic with three rows instead of a
cylinder's two. The meeting is S9f.2b's (`spline_crossing::meeting_with`,
the partner's rows and labels): branches over the run between the roots of
`D`, graphs over the height about each turning point inside both faces,
their switches vertices, found once per pair of a spline wall and a
hemisphere (S9d.1's split; a turning point inside one hemisphere lies
outside the other, whose branches end there). A sphere smaller than the
wall's height straddling it meets it in a loop; a larger one, or one
crossing a cap, in branches whose turning points lie outside a face. A
loop's switch lies a sixty-fourth of the way from its turning point to
where the gentler branch's slope has fallen to one (a cylinder's a
quarter): its graph over the height spans the sphere's section there, and
at a quarter the slowest fuzz variants spent ten times longer integrating
it. A turning point on the hemispheres' split, which is no edge of the
input, tries the split at another seam; one inside both faces but within
the resolution of a cap's or a rim's plane is `Degenerate` like one on it
(a sphere about a turned prism's frame origin).

Vertices (`curved/spline_sphere.rs`): a spline cap edge or crease against
the sphere at the roots of the sphere's function along it, degree `2 p`
(`wallcrv_sphere`); the spline prism's vertical edges at quadratic surds
(S9d.1's `line_sphere`); a sphere's circle (a rim of a cap or zone, the
split's great circle) against a spline wall along its plane's crease on
the wall, degree `2 p`, its place on the circle read off rationally
(`circ_wall`), or, in a plane holding the wall's axis direction (a cap's
split great circle on the prism's axis, a rim on its side), by a primitive
element of the tower: each arc's implicit equation at the circle's
projection `l0 + dx lx + dy ly`, reduced by `dx^2 |x|^2 + dy^2 |y|^2 = r2`
to `E(dx) + dy O(dx)`, vanishes where `E^2 - (r2 - |x|^2 dx^2) / |y|^2 O^2`
does (degree at most `2 p`), and there `dy = -E / O` in the same field;
where `O` vanishes at a root (a great circle whose `x` is the wall's axis:
its points on one generatrix share `dx`) the roles of `dx` and `dy` are
swapped (`tower_points`). In the topology the meeting is `Curve3::WallMeet`
with `other_sphere` (the sphere's stored frame and radius;
`TOPOLOGY_MODEL.md`), certified over the sphere's three rows
(`MATHEMATICS.md`); on the sphere its pcurve is the sphere's inverse.

`Degenerate`: the sphere tangent to the wall ("a sphere tangent to a spline
wall"); a turning point at an interior knot; a turning point on a face's
boundary (a cap's edge, a rim, a segment's end); a vertex's polynomial with
a multiple root (an edge tangent to the other's face, a circle tangent to a
generatrix). A spline prism against a cone is S9f.3b's (below).

### Spline walls against cones (S9f.3b)

A spline prism meets a cone or a frustum (`Solid::cone_with`) in any
position, either the object. Along the wall's ruling the cone's function on
its exact model (`procedural::other_cone`: `u^2 + v^2 - (b + k w)^2` in its
frame's exact rows, `k` its rational slope) is `A w^2 + 2 B(tau) w +
C(tau)`: S9f.2b's quadratic with a radius row of negative sign
(`spline_crossing::terms`: every row carries its sign, so the coefficients,
a point's `F_w`, the tangent, the binary64 views and the graphs over the
height are one code for cylinders, spheres and cones). `A = q_u^2 + q_v^2 -
k^2 q_w^2` (`q` the prism's axis in the cone's rows) is one constant for
the pair, and its sign is the case. `A > 0` (the prism's axis farther from
the cone's than its half angle): each ruling meets the quadric twice on one
nappe or not at all, and the meeting is S9f.3a's (branches over the run
between the turning points, loops' graphs over the height about turning
points inside both faces, switched a sixty-fourth of the way as a sphere's;
a turning point within the resolution of a rim's plane `Degenerate`). `A <
0` (within the half angle, the prism's axis along the cone's among them):
each ruling meets the double cone once on each nappe, `D > 0` but where a
ruling passes the apex, and the plus and minus branches run over the whole
run, one on each nappe; the other nappe's lies beyond the apex, outside the
cone's face (a cone's radii are nonnegative, so its apex, real or virtual,
lies at or beyond an end), and the arrangement drops it as any piece outside
a face. `A = 0` (the axis exactly along a generatrix direction: one finite
root per ruling, running to infinity where `B` vanishes) is
`Degenerate("a spline wall along a cone's ruling")`, as S9d.3c's cylinder
along a cone's ruling; the apex, a cone's or a frustum's virtual one, on a
spline wall's surface is `Degenerate("a cone's apex on the other input's
surface")`, `cone_pair`'s rule; both are checked first
(`curved/spline_cone.rs`'s `meeting`). The cone pairs' pass leaves spline
walls to it.

Vertices (`curved/spline_cone.rs`): a spline cap edge or crease against the
cone at the roots of its function along it, degree `2 p` (`wallcrv_cone`);
the spline prism's vertical edges by S9d.3a's `line_cone`; a rim (a circle
of rational radius on the cone's stored axes) against a spline wall along
its plane's crease in the rim's own elliptic cylinder (the rows dual to the
conic's axes, `alpha^2 + beta^2 - 1`, degree `2 p`, placed by its angle:
`spline_crossing::conic_crease`), or, in a plane holding the wall's axis
direction (the cone's axis across the prism's), at S9f.2b.2's points in the
rim's half-angle chart (`tower_points`: its radius rational, no elimination
needed). In the topology the meeting is `Curve3::WallMeet` with
`other_half_angle` (the cone's stored base frame, bottom radius and half
angle; `TOPOLOGY_MODEL.md`), certified with its polynomials kept per power
of the half angle's tangent (`MATHEMATICS.md`); on the cone its pcurve is
the cone's inverse.

`Degenerate`: the two above; the cone tangent to the wall ("a cone tangent
to a spline wall"); a turning point at an interior knot or on a face's
boundary, as S9f.3a's with the cone's reasons; a vertex's polynomial with a
multiple root; and the engine's rules.

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
* **S9d.3c evidence (cones' loops against cylinders, cones and spheres,
  turned caps against cones), before its kernel code.**
  `cones_boolean_reference.py` takes a turned cap as the kernel's `Ball`
  reads it (S9d.2c's `AxisSphere`: its end plane through `o + h n` normal
  to the stored axis), its circle's edge and its disc swept on vectors
  exactly in that plane (`y x a`, `a x (y x a)`; the frame's `x` and `y` in
  an exact frame, S9d.3b's rows unchanged). `generate_cones_loops_boolean_
  fixtures.py --check` writes `boolean-cones-loops-cases.txt`,
  `-expected.tsv` and `-frames.tsv` with 31 cases (3 fuses, 14 cuts, 14
  commons; 29 solid, 2 degenerate; 9 in exact frames): cone-cylinder loops
  (a rod across a frustum's axis grazing its wall, all three; the rod less
  it; a rod beside a cone's apex; a rod in `LEAN`; a vertical rod against a
  frustum in `TILT`), cone-cone loops (a thin frustum across the axis
  grazing the wall, all three; a frustum in `TILT` against one in `LEAN`; a
  cone in `LEAN` whose direction cone crosses the frustum's, its curve
  through infinity), turned cones' loops against spheres (in `TILT`, all
  three; `LEAN`; `R125`), turned hemispheres against cones (a coaxial cone,
  a frustum off the axis under the lower hemisphere in `LEAN`, a frustum in
  `LEAN`), and `degenerate` a cylinder exactly along a cone's ruling (the
  cone's slope `fl(0.6) / fl(0.8)`, `LEAN`'s stored axis: the cylinder's
  `A` zero) and a turned rim tangent to a cone within rounding. Checks
  before writing: closed forms (a rod across the axis by S9d.3b's strip, a
  whole sphere against a turned cone by the lens along its ideal axis, a
  hemisphere against a coaxial cone by circular segments along the axis,
  its sphere face and the wall by their arcs above the chord, its disc by
  its chords in its plane) within 2.2e-40 in exact frames and 1.3e-16 in
  turned ones, an input cut in two along its axis (a prism in its frame, a
  frustum in an exact frame at a dyadic height) against the whole within
  5.3e-42, a hemisphere and its complement against the whole sphere within
  1.4e-41, inclusion and exclusion 1.2e-41, the area identity and every
  face's classes 6.9e-41, the cone's wall two ways 2.7e-41, a second
  direction 1.7e-41, Monte Carlo 2.5 standard errors, no near coincidence
  but the declared pairs', every other edge and vertex at least 0.0156
  from tangency with or incidence on the other's surfaces, every plane at
  least 0.05 of the radius from tangency to the other's spheres and
  parallel cylinders; Python 3.9 and
  3.12 write the same files. `test_cones_loops_boolean_reference.py` checks
  a circular segment by polar quadrature, a cone wider than the ball giving
  the hemisphere, the split parts on the input's exact model, a split on a
  coaxial pipe and the tangent cone's construction.
  `compare_cones_loops_boolean.py` reproduces
  `occt-boolean-cones-loops-preimplementation`
  (`rust_cones_loops_boolean_exists` false; the kernel `unsupported` on all
  31): every result valid with the reference's solid count, 7 match (the
  turned hemispheres within 1.8e-8), 24 reviewed
  (`occt-boolean-cones-loops-divergences.json`: BRepGProp's default
  integration on approximated loops, up to 5.8e-6; adaptively converged
  within 1.7e-8), five solids' counts change when unified.
* **Kernel (S9d.3c).** `tests/cones_loops_booleans.rs`: all 31 fixtures as
  the reference (29 within the kernel's enclosures, the cylinder along a
  cone's ruling and the turned rim tangent within rounding refused), every
  history checked, results deterministic and moved rigidly, the loops'
  `Curve3::Meet` pieces over both carriers within `1e-9` of both surfaces,
  the turned cones' `Curve3::Rise` pieces within `1e-12`, the turned rims'
  crossings on the rim's plane. `compare_cones_loops_boolean.py`: 5 matches,
  26 reviewed (the measures, and the kernel's edges and vertices where it
  splits its loops at its switches), no failure; `compare_cones_boolean.py`
  unchanged in its counts with `ball_tilt_common` now within the reference;
  every other Boolean comparison unchanged.
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
* **S9d.4c evidence (a sphere's cap or zone against a torus, torus
  v-segments and wedges against curved solids), before its kernel code.**
  `torus_parts_boolean_reference.py` (mpmath, 40 digits) is S9d.4b.2's
  reference (every face of both inputs swept by two families of curves, the
  divergence theorem on the classified pieces) with S9d.4b.1's part model
  (a segment's end planes at the stored heights and its wall over its arc,
  oriented by its `sigma`: the inner half inside out; a wedge's half-planes
  on the rounded end direction, its wall over the turn and its meridian
  discs) and a cap or zone as S9d.2c's `AxisSphere` reads it (its wall
  about the unit axis between its ends' latitudes, its discs by chords);
  solids by the part's normal slices, a region sheared between two slices
  joined through intervals at levels between them. `generate_torus_parts_
  boolean_fixtures.py --check` writes `boolean-torus-parts-cases.txt`,
  `-expected.tsv` and `-frames.tsv` with 53 cases (7 fuses, 23 cuts, 23
  commons; 51 solid, 2 degenerate; 45 in exact frames) on the torus of
  radii 5/2 and 1: caps and zones against it (coaxial ones in circles, a cap
  on the tube's top with a rim of surd radius across it, a zone about `x`
  across the outer side, a hemisphere in `LEAN`, a zone against the torus
  in `TILT`), the outer and inner halves against coaxial pipes, a cone, a
  sphere and a lower hemisphere, a sphere about a point of the outer half's
  rim and a torus ringing its outer side, the band `0.5..2.25` against a
  sphere, a frustum and a ringing torus across its rims and in `TILT`
  against a coaxial pipe, wedges of a quarter, a half and three quarters
  against a pipe, a small torus, a sphere and a frustum across their end
  discs and rims (the quarter in `LEAN` too); `degenerate` a cap's rim
  within the resolution of tangency to the outer equator and a quarter
  wedge's start rim touching a sphere. Checks: closed forms of the 9
  coaxial pairs (both inputs' radial intervals along the axis, the walls by
  their own elements) within 6.9e-40 in exact frames and 3.8e-17 in turned
  ones; a half and the other half against the tool as the whole torus, a
  cap or zone and the rest of its sphere against the torus as the whole
  sphere, within 4.3e-42; every operation two ways 1.5e-31, both inputs
  from their faces 8.3e-35, inclusion and exclusion 1.8e-35, the area
  identity 8.3e-35, Monte Carlo 2.5 standard errors; no near coincidence
  outside the declared pairs, surfaces meeting at a sine of at least 0.38,
  edges crossing at 0.29, every plane at least 0.0059 of the case's size
  from tangency to the other's spheres and tori.
  `test_torus_parts_boolean_reference.py` checks the halves' and a
  wedge's closed forms against their volumes of revolution, a half's
  complement, and the reference on the outer half about a pipe, a
  hemisphere in the tube and a quarter wedge in a ball. No protocol row is
  new. `compare_torus_parts_boolean.py` reproduces
  `occt-boolean-torus-parts-preimplementation`
  (`rust_torus_parts_boolean_exists` false, keyed on
  `solid/boolean/curved/torus_parts.rs`; the kernel's probe `unsupported`
  on all 53): every solid count the reference's; 26 match (the coaxial
  pairs within 2.4e-14) and 27 are reviewed
  (`occt-boolean-torus-parts-divergences.json`): BRepGProp's default
  integration on faces bounded by B-spline sections misses by up to 7.4e-6,
  adaptively within 2.3e-8 and unchanged at 1e-12, and the quarter wedge in
  `LEAN` against a sphere is wrong natively (its common 1.6e-4 off and
  invalid under `BRepCheck_Analyzer`); fourteen solids' counts change when
  unified.
* **Kernel (S9d.4c).** `tests/torus_parts_booleans.rs`: all 53 fixtures as
  the reference (51 within the kernel's enclosures, each at most `1e-9`
  wide; the cap's rim within the resolution of tangency and the wedge's
  rim touching a sphere refused), every history checked, results
  deterministic and moved rigidly, the coaxial pairs meeting in circles and
  the others in `Curve3::Toric` pieces, a band, a quarter wedge and a cap
  turned by a rigid motion with their tools keeping the reference's
  volumes, a void under a band's hole a cavity (`ComputationLimit`) and a
  loop through a sphere's pole winding by its pcurves' ends.
  `compare_torus_parts_boolean.py`: 16 matches, 37 reviewed (the capture's
  27, 23 of them with the kernel's entity counts, and ten native matches by
  their counts alone: coaxial bands crossing OCCT's seams, and sections
  split at OCCT's points against the kernel's exact pieces), no failure;
  every other Boolean comparison unchanged in its counts.
* **S9e.1 evidence (a Boolean's result given to another Boolean), before
  its kernel code.** The case protocol chains a second Boolean on the
  first's one solid: a `then OP ID` row (`then OP ID swapped` when the
  first result is the tool) and a third prism's rows after the tool's
  (`identity_reference.encode_chained_case`, `native_chained_case`; the
  native oracle and `tests/support/boolean_protocol.rs` read it, the rows
  the second Boolean's). `chained_boolean_reference.py` (mpmath, 40 digits)
  is S9c.1's reference generalized to three prisms and the chained set
  function `op2(op1(a, b), c)` (swapped `op2(c, op1(a, b))`): slices in
  planes holding every axis (the three axes take at most two directions),
  each slice the first operation's convex pieces of the first two
  sections' parallelograms and then the second's against the third's,
  breakpoints over the three sections' lines; every face of each prism
  swept against both others at once, a piece bounding the chain where the
  chained function differs across the face, one on faces of several prisms
  counted once. `generate_chained_boolean_fixtures.py --check` writes
  `boolean-chained-cases.txt`, `-expected.tsv` and `-frames.tsv` with 30
  cases (10 chains, each operation; 23 solid, 1 empty, 6 degenerate) in
  `XY`, `SIDE` and `TILT`: a box less a tilted hole halved by a box whose
  wall holds the hole's axis, and as the tool of a slab crossing the hole;
  a box fused with a tilted pin, stepped, and as the tool of a slab above
  it; a box on the holed top face (coplanar, edges crossing there); a box
  less a groove across its top (its top face in two result faces) drilled
  through one of them; a column cutting the hole in a circle; a quarter
  cylinder bored coaxially; `degenerate` a third cylinder tangent to the
  groove and a third box whose edge lies on the groove's wall. Checks:
  closed forms within 1.3e-40; S9c.1's pair reference for the first result
  in the volume identities within 1.3e-41 and the area identity (no face
  of the third on another's) within 7.0e-41; face classes their closed
  forms within 2.6e-41; Monte Carlo 2.7 standard errors; no near
  coincidence. `test_chained_boolean_reference.py` checks three boxes by
  their grid cells exactly, a far third prism and the common's symmetry.
  `compare_chained_boolean.py` reproduces
  `occt-boolean-chained-preimplementation`
  (`rust_chained_boolean_exists` false, keyed on
  `solid/boolean/curved/given.rs`; the kernel's probe `unsupported` on all
  30): every result valid with the reference's solid count, all 30 match
  (within 8.0e-10; exact frames 2.9e-16), no review; seven solids' counts
  change when unified.
* **Kernel (S9e.1).** `tests/chained_booleans.rs`: all 30 fixtures as the
  reference (24 within the kernel's enclosures, each at most `1e-9` wide;
  the 6 declared degenerate refused), every second history complete over
  the first result's ids, results deterministic and moved rigidly, first
  results translated with their third prisms keeping the reference's
  volumes, the grooved top's two result faces each named by its own part,
  the stored frames the reference's bit for bit; `tests/curved_booleans.rs`
  a cylinder's cap circle across a parallel cylinder's wall, either input
  first, by its lens's closed form. `compare_chained_boolean.py`: 24
  matches, 6 reviewed for their entity counts (OCCT's vertices at the
  hole's and pin's stored cylinders' seam angles and a touching box's
  imprints kept after unifying, found by a diagnostic build printing the
  unified vertices), no failure; every other Boolean comparison unchanged
  in its counts.
* **S9f.1 evidence (spline prisms against polyhedral prisms in any
  position), before its kernel code.** `curved_boolean_reference.py`
  (mpmath, 40 digits) takes spline walls: each of S8b's profile splines is
  its exact Bezier spans (blossoms), one element per span, its joins the
  profile's vertices; a line meets a span at the real roots of the
  degree-`p` polynomial `(S(tau) - b) x d` (monotone runs, Newton steps; no
  resultant); S9c.1's slices and face sweeps keep their structure, a spline
  crossing implicit: its events with a line's crossing moving linearly are
  the exact degree-`p` roots of `(S(tau) - Q0) x Q1`, tangencies those of
  `S'(tau) x d`; a span's wall is swept by its generatrices, its crossings
  of the other's planes exact polynomials in `tau`; a profile with splines
  takes the symmetric Green forms on every element; and each face sweep
  also integrates `X . N` and `x_i^2 N_i`, so every operation's volume and
  moments come a second way, by the divergence theorem over the kept
  pieces. `generate_spline_any_boolean_fixtures.py --check` writes
  `boolean-spline-any-cases.txt`, `-expected.tsv` and `-frames.tsv` with 38
  cases (12 fuses, 13 cuts, 13 commons; 32 solid, 1 empty, 5 degenerate):
  the bulge under a tilted box (an oblique crease), the dome between two
  `SIDE` planes (generatrices; the common `44/3` exactly), a leaning box
  whose edges pierce the blob's wall, the wave in a `TILTX` slab, a turned
  (`R125`) box on the capsule's base (coplanar caps) and one standing on
  its top (a face shared with the opposite orientation), a tilted pin
  across the lens hole's walls (a common of two solids), the bulge against
  a box in the same `TILT` axes at an offset that rounds (`bulge_offset`), the
  dome in `TILT` under an `XY` box, R4's `kink` (a quadratic whose apex is
  a knot of multiplicity two, C1) in `TILT` against a leaning box (its cut
  two solids), the blob as a leaning tool through a box; `degenerate` a
  `TURN` wall tangent along the dome's apex generatrix, one tangent to the
  kink at its knot, and a `TILTX` cap plane 8.9e-17 from a `TILT` blob's
  axis. Checks: the divergence theorem 4.2e-41, inclusion and exclusion
  4.5e-41, the area identity 1.2e-40, every face's classes 3.5e-41, a
  second slicing direction for parallel axes 2.1e-41; S8b's split reference
  on each spline profile's prism against a half-space box (oblique and
  parallel planes, the prism in `XY` and `TILT`) 1.4e-41 exact and 1.9e-18
  turned; S9a.2's `SplinePair` on `bulge_offset` (its offset taken exactly)
  and on its own 17 one-spline fixtures 2.2e-17 (the 10 in `XY` to every
  printed digit); solid counts as declared; margins outside the declared
  pairs (vertices 0.019 from the other's faces, creases' extremes 0.12
  from the caps, parallel planes 0.16 from tangency, edges crossing spline
  walls at a sine of 0.37, oblique planes at 0.48 to the axis), the
  declared pairs' below 1e-12. `test_spline_any_boolean_reference.py`
  checks the root finder, the implicit events, every spline profile's
  moments against S9a.2's Green, the dome between side planes both ways,
  the degenerate margins and the kink's liftings. No protocol or oracle row
  is new. `compare_spline_any_boolean.py` reproduces
  `occt-boolean-spline-any-preimplementation`
  (`rust_spline_any_boolean_exists` false, keyed on
  `solid/boolean/curved/spline_walls.rs`; the probe `unsupported` on all 38,
  32 from the polyhedral engine and the capsule's 6 from the curved
  engine's model): every result valid with the reference's solids; 23 match
  and 15 are reviewed (`occt-boolean-spline-any-divergences.json`):
  BRepGProp's default integration on faces bounded by B-spline edges
  misses by up to 6.7e-7 (the wave 1.0e-3), adaptively within 6.0e-9 in
  volume but for the wave, while Green's theorem over OCCT's own faces and
  pcurves (a diagnostic build) gives the reference within 1.3e-8 on all
  38; two solids' counts change when unified. With S9f.1's kernel: 22
  matches and 16 reviewed, the kernel within the reference on the 33
  results and refusing the 5 degenerate cases; the new review is
  `capsule_stand_cut`'s counts (OCCT's unified cut keeps the capsule's top
  spline edge split in four where the standing box's footprint crosses it,
  three vertices a diagnostic build printing the unified solid's vertices
  placed). `boolean-spline-any-r4-*` (written by the same generator after
  the capture, 6 cases, no native rows) adds R4's rounded knot in `TILT`
  and `TILT`'s walls against planes exactly parallel to its axis (`XY`'s
  `x` planes, `SIDE`'s caps), checked by the kernel's tests.
* **S9f.2b evidence (spline walls against cylinder walls on crossing
  axes), before its kernel code.** `curved_boolean_reference.py` takes the
  other prism's arcs and circles when the axes cross: the slicing's planes
  hold both axes, so a spline wall's section is still its chords'
  parallelograms and a cylinder's its chords' (surds), and an event where a
  spline chord, a cylinder's chord and a height line are concurrent is a
  point of the walls' meeting at a cap's height: along the spline wall's
  ruling `X = o + S(tau) + w n` the cylinder's equation in its own frame is
  `A w^2 + 2 B(tau) w + C(tau)` (degrees 0, `p`, `2 p`), so at the spline
  prism's cap (`w` fixed) or the cylinder prism's (`w` linear in `S(tau)`)
  the points are the real roots of an exact polynomial of degree `2 p`
  (`wall_cylinder_events`). A cylinder wall's generatrices swept against a
  spline wall cross it implicitly, their events the same meetings and the
  generatrices whose projection along the spline's axis touches a span;
  the caps' sweep lines meet the other's walls at the same roots. The
  meeting's turning points (`wall_cylinder_turns`, the roots of `B^2 - A
  C`) are the trace's tangencies, already events. `generate_spline_crossing_boolean_fixtures.py
  --check` writes `boolean-spline-crossing-cases.txt`, `-expected.tsv`
  (its `expect` rows name the sub-step: `S9f.2b.1`, or `S9f.2b.2` for a
  meeting turning back inside the faces) and `-frames.tsv` with 34 cases
  (11 fuses, 11 cuts, 12 commons; 30 solid, 4 degenerate; 6 of S9f.2b.2):
  tilted rods along the bulge's and the capsule's walls (the capsule's
  crossing its double knot), a leaning rod along the wave's top across its
  knot, a perpendicular rod covering the dome, a steep cylinder holding
  most of the blob as the object (`STEEP`, the normal (0, 5, 12) / 13), a
  tilted rod ending inside the blob, a stadium's arc and edges against the
  bulge, a tilted ring's hole around the dome; loops (a perpendicular rod
  through the bulge's wall, a tilted rod through the lens's); and
  `degenerate` a rod touching the dome's apex and one whose meeting with
  `knot`'s wall turns back at its knot. Checks: the divergence theorem,
  inclusion and exclusion and the area identity within 9.4e-41, every
  face's classes 4.5e-41, the perpendicular pairs' commons as the product
  of the profile's chord along `x` and the disc's height chord within
  7.1e-44; solid counts as declared; margins outside the declared pairs at
  least 0.011 (vertical edges piercing curved walls, vertices 0.037 from
  the other's faces, meetings crossing caps at 0.11, turning points 0.56
  outside a face, the loops' 1.5 inside both faces), the declared pairs'
  below 1e-12. Python 3.9 and 3.12 write the same files.
  `test_spline_crossing_boolean_reference.py` checks the meeting at a cap
  and its turning points against closed forms, the tower field's points on
  a cap holding the wall's axis direction, the product of chords and the
  degenerate margins. No protocol or oracle row is new.
  `compare_spline_crossing_boolean.py` reproduces
  `occt-boolean-spline-crossing-preimplementation`
  (`rust_spline_crossing_boolean_exists` false, keyed on
  `solid/boolean/curved/spline_crossing.rs`; the probe `unsupported` on all
  34, refused by the curved engine's `spline_pairs`): every result valid
  with the reference's solids; 13 match and 21 are reviewed
  (`occt-boolean-spline-crossing-divergences.json`): BRepGProp's default
  integration on faces bounded by B-spline intersection edges misses by up
  to 1.5e-6 (`knot_turn_cut` 1.3e-5, the wave 2.9e-4), while a diagnostic
  build's adaptive BRepGProp or Green's theorem over OCCT's own faces and
  pcurves gives the reference within 4.0e-9 in volume, 1.4e-9 in area and
  8.6e-9 in the centre on all 21; five solids' counts change when unified.
  With S9f.2b.1's kernel: 13 matches and 21 reviewed, the kernel within
  the reference on the 24 results of S9f.2b.1 with OCCT's unified counts
  on every solid but `capsule_tilt`'s three (OCCT splits the meeting's
  ellipse arcs on the capsule's caps at the rod's seam generatrix, and in
  the fuse its seam at the hole: those three reviews now name the counts
  too), refusing the 4 degenerate cases and S9f.2b.2's 6 (`unsupported`).
* **S9f.3a evidence (spline walls against spheres), before its kernel
  code.** An independent reference of its own,
  `spline_sphere_boolean_reference.py` (the spline prism's exact model and
  profile parsing shared with `curved_boolean_reference.py`): the pair
  sliced by planes, each slice's sections exact regions of the profile's
  plane (the profile cut by the caps' lines, the sphere's circle projected
  along the prism's axis and cut by a hemisphere's line), their Boolean
  pieces by Green's theorem in closed form, breakpoints at every vertex's
  slice and every edge's extremes (the meeting's from `F` and its tangent
  condition, `w` eliminated exactly), volumes and moments two ways (along
  the caps' normal and along `(2, -3, 5)`), the walls swept along their
  generatrices, the sphere's face by Archimedes' area element, a
  hemisphere's disc in its plane, solids by the slices' union-find.
  `generate_spline_sphere_boolean_fixtures.py --check` writes 33 cases (27
  solid, 6 degenerate): loops on the bulge, the lens, the blob under its
  top cap, `knot`'s span in `TILT` and a sphere object across the wave's
  knot, branches over the dome and under the blob in `TILT`, hemispheres on
  the bulge's axis and on its side (towers on the split's great circle and
  on the rim), and declared degenerate a touch at the dome's apex, a turning
  point at `knot`'s knot and a loop turning back on a cap's edge. The two
  slicings within 1.8e-41 of the size, the inputs' closed forms 5.7e-42,
  faces' classes 2.8e-41, the area identity 8.0e-41, Monte Carlo within 2.8
  standard errors, margins outside the declared pairs at least 0.023,
  Python 3.9 and 3.12 the same files. The capture
  `occt-boolean-spline-sphere-preimplementation` (keyed on
  `curved/spline_sphere.rs`, the probe `unsupported` on all 33): every
  result valid with the reference's counts, 4 matching and 29 reviewed
  (BRepGProp's default measure up to 4.1e-5, the wave's 1.0e-3; adaptive
  BRepGProp or Green's theorem over OCCT's faces within 6.0e-9 in volume,
  7.8e-9 in area and 8.5e-9 in the centre, the wave's cut and common
  within 9.3e-7 inside their edges' tolerance of 3.7e-5). With S9f.3a's
  kernel: 4 matches and 29 reviewed, the kernel within the reference on all
  27 results and refusing the 6 degenerate cases; the 18 loop cases'
  reviews name their counts (the kernel's loops cut at their switches and
  at the hemispheres' split, OCCT's edges at its own seams).
* **S9f.3b evidence (spline walls against cones), before its kernel
  code.** `spline_cone_boolean_reference.py`, an extension of S9f.3a's (its
  prism, profile elements, Green's integrals of the pieces reused): the cone
  on an exact frame, the pair sliced along two directions cutting it in
  ellipses (its axis leaned toward the prism's), each slice's cone section
  its quadric's restriction on its principal axes projected along the
  prism's axis and cut by the end planes' lines; the cone's wall by its area
  element `sqrt(1 + k^2) |G(s)| / D(theta)^2` in closed form over the
  section's arcs inside the prism (`r = G / D` along the slice); the prism's
  walls along their rulings with `A` of either sign; the caps by chords
  against the cone's conic of any type; the end discs as the hemisphere's;
  solids by intervals on chords of slices across the prism's axis.
  `generate_spline_cone_boolean_fixtures.py --check` writes 27 cases (21
  solid, 6 degenerate): `A < 0` a frustum on the bulge's axis across both
  caps, a cone with its apex inside the lens prism (the other nappe meeting
  the walls outside the cone's face), a cone hanging over the blob with its
  apex below, a frustum object across the wave; `A > 0` a thin frustum
  piercing the dome's arch (a loop), one whose rim's plane holds the dome's
  axis (the rim's tower points, the top cap cutting the cone in a
  hyperbola), a frustum on `z` against `knot` in `TILT` (a loop); declared
  degenerate the apex on the dome's arch, `A = 0` exactly (the bulge in
  `TILT` against a cone of the stored axis's slope) and a cone touching the
  arch. The two slicings within 1.1e-39 of the size, the inputs' closed
  forms 1.0e-39, faces' classes 1.4e-40, the area identity 9.2e-41, Monte
  Carlo within 2.1 standard errors, margins outside the declared pairs at
  least 0.0156 (`|A|`), Python 3.9 and 3.12 the same files. The capture
  `occt-boolean-spline-cone-preimplementation` (keyed on
  `curved/spline_cone.rs`, the probe `unsupported` on all 27): every result
  valid with the reference's counts, 4 matching and 23 reviewed (BRepGProp's
  default measure up to 5.5e-6 on the quadratic and cubic walls, 7.4e-5 on
  `knot` in `TILT`, the wave's 1.0e-3; adaptive BRepGProp or Green's
  theorem over OCCT's faces within 3.2e-8 in volume, 7.0e-9 in area and
  9.0e-9 in the centre). With S9f.3b's kernel: 4 matches and 23 reviewed,
  the kernel within the reference on all 21 results and refusing the 6
  degenerate cases; the 6 loop cases' reviews name their counts (the
  kernel's loops cut at their switches; in `knot_tilt_cone_fuse` OCCT's
  cone face outside the prism two faces, the kernel's one).
* **S9f.2b.2 evidence (loops and towers), before its kernel code.** The
  same generator writes 51 cases, 17 new (`step` `S9f.2b.2` exactly when a
  turning point lies inside both faces or a cap circle meets a wall in a
  plane holding its axis inside its heights, `towers`): a rod ending
  inside the dome on its arch and a wide rod ending inside the lens across
  its two cubics (towers alone), loops cut by a tower's cap circle and by
  the spline prism's top cap (the bulge) and by a tilted rod's cap circle
  (the lens), and `cap_turn`, declared degenerate (a loop's turning point
  on the rod's cap rim). Checks as S9f.2b's within 9.4e-41, margins outside
  the declared pairs at least 0.006, Python 3.9 and 3.12 the same files;
  the reference test checks the dome's tower points in closed form. The
  capture was taken again, whole (`rust_spline_loops_boolean_exists`
  false, keyed on S9f.2b.2's two refusals in `spline_crossing.rs`): the 34
  earlier rows to the bit, the 17 new done and valid with the reference's
  counts, 4 matching and 13 reviewed (BRepGProp's default measure; adaptive
  BRepGProp or Green's theorem over OCCT's faces within 2.1e-9 in volume,
  3.7e-9 in area, 1.8e-10 in the centre). With S9f.2b.2's kernel: 14
  matches and 37 reviewed, the kernel within the reference on all 45
  results and refusing the 6 degenerate cases; the 15 loop cases' reviews
  name their counts (the kernel's loops cut at its switches and at its
  rod's faces' boundaries, OCCT's at its rod's seam touching a turning
  point and at its intersection edges' own splits), `bulge_cap_loop`'s
  three for counts alone.
* **S9f.2a evidence (spline walls against arc, circle and spline walls on
  exactly parallel axes), before its kernel code.**
  `curved_boolean_reference.py` takes the other prism's arcs, circles and
  splines when the axes are exactly parallel: every wall is swept along the
  common axis, so a spline's chord meets the other's arc or spline chord,
  and a curved wall's generatrices change class, only where the slice's
  trace or the generatrix passes a 2D crossing of the profiles' curves
  projected along the axis (`parallel_crossings`, in the exact affine map
  between the frames' `(u, v)`): a span against a circle at the real roots
  of the exact degree-`2 p` polynomial `|S(tau) - c|^2 - r^2` of the span
  mapped into the circle's frame (Yun's factors: a tangency a root of even
  multiplicity), two spans by subdivision of their Bezier forms in
  fractions while their control boxes meet, then Newton's method at 40
  digits (on the distance's gradient where tangent), never a resultant;
  each crossing's slice parameter and face parameters are breakpoints.
  `generate_spline_parallel_boolean_fixtures.py --check` writes
  `boolean-spline-parallel-cases.txt`, `-expected.tsv` and `-frames.tsv`
  with 43 cases (15 fuses, 14 cuts, 14 commons; 39 solid, 4 degenerate):
  spline walls against a disc (the bulge, `R125`), a stadium's arc and
  line walls (the dome in `TILT`, the stadium in `TILT2`: a quarter turn
  about the tilted axis), a holed square's round hole and lines (the blob,
  `TURN`), a disc across the lens hole (`FLIP`, a half turn), a disc across
  the capsule's spline; spline walls against spline walls (the bulge and a
  turned dome, the blob and an `R125` lens, the wave and its half-turned
  mirror crossing it four times (a common of two solids), two domes about
  the tilted axis, a quartic hump against a cubic lens, two capsules);
  rounding offsets in `TILT` (the bulge against a disc and a dome); and
  `degenerate` a disc's cylinder and a half-turned dome each tangent to the
  dome's wall along its apex generatrix. Checks: the divergence theorem
  1.3e-40, inclusion and exclusion 1.3e-40, the area identity 9.8e-41,
  every face's classes 3.0e-41, a second slicing direction 1.3e-40; S9a.2's
  `SplinePair` on the 33 cases whose map is an exact turn or an offset in
  one frame (the tool's profile turned exactly into the object's frame)
  2.4e-39 in `XY` and 2.2e-17 in `TILT`; solid counts as declared; margins
  outside the declared pairs (crossings at a sine of 0.47, near misses
  0.37 apart, vertices 0.13 from the other's faces, creases 0.06 from the
  caps), the declared pairs' crossing sines below 5.7e-21. Python 3.9 and
  3.12 write the same files. `test_spline_parallel_boolean_reference.py`
  checks the exact map between parallel frames, a span's crossings with a
  circle (closed forms, a tangency's multiplicity, a near miss) and with
  another span (crossing, touching and missing parabolas), containment in
  closed form both ways, S9a.2's `SplinePair` on a quarter-turned pair and
  the degenerate margins. No protocol or oracle row is new.
  `compare_spline_parallel_boolean.py` reproduces
  `occt-boolean-spline-parallel-preimplementation`
  (`rust_spline_parallel_boolean_exists` false, keyed on
  `solid/boolean/curved/spline_parallel.rs`; the probe `unsupported` on all
  43, refused by the curved engine's `spline_pairs`): every result valid
  with the reference's solids; 31 match and 12 are reviewed
  (`occt-boolean-spline-parallel-divergences.json`): BRepGProp's default
  integration on faces bounded by B-spline edges misses by up to 1.9e-7
  (the waves 8.6e-4), adaptively within 2.3e-9 in volume but for the
  waves, while Green's theorem over OCCT's own faces and pcurves (a
  diagnostic build) gives the reference within 1.1e-8 on all 39 results;
  four solids' counts change when unified. With S9f.2a's kernel: 31
  matches and 12 reviewed, unchanged, the kernel within the reference on
  the 39 results with OCCT's unified counts on every solid, and refusing
  the 4 degenerate cases.
* **S9e.2 evidence (a stack, an S9b.1 result or one solid of several given
  to another Boolean), before its kernel code.** The protocol's `then` row
  may end `solid X Y Z`: the first result's solid holding the point
  strictly inside is the second's argument (natively
  `BRepClass3d_SolidClassifier`, in the kernel `Solid::classify`).
  `chained_boolean_reference.py` takes a fourth prism, a selector box `D`
  holding that solid and no part of another: the chain `op2(op1(a, b) and
  d, c)` (swapped `op2(c, op1(a, b) and d)`), each slice's first pieces
  clipped to `D`, every face swept against the three others, and the area
  of `D`'s faces with the first result on either side (`separation`, zero
  when `D` picks whole solids). `generate_given_boolean_fixtures.py
  --check` writes `boolean-given-cases.txt`, `-expected.tsv` (its `expect`
  rows naming the class: `stack`, `polyhedral`, `several`) and
  `-frames.tsv`: 36 cases (12 chains, each operation; 30 solid, 6
  degenerate) in `XY`, `SIDE`, `DOWN`, `TILT` and `R125`: DRAW's rollex of
  `bcut_simple/L3` to `L6` (a disc less a pocket across its rim, a stack,
  then a cylinder on the pocket's floor, turned and in the stack's frame);
  a box with a boss bored coaxially, sliced by a tilted slab and as a
  tilted box's tool; a box fused with a turned box and a box with a
  slanted top (S9b.1), each bored; a box severed by a cylinder, its lower
  solid drilled and as a box's tool; a box severed by a tilted slab, its
  lower piece bored; `degenerate` a cylinder tangent to the rollex's rim
  and a box tangent to the severed box's other solid. Checks: closed forms
  within 3.3e-41; the pair identities for the given solid (a picked one's
  measures from the three-prism chain `op1(A, B) and D`) within 7.7e-42 and
  the area identity 8.2e-41; separation exactly zero; face classes
  2.2e-41; Monte Carlo 3.1 standard errors; no near coincidence.
  `test_given_boolean_reference.py` checks four boxes by their grid cells
  exactly, a selector across the first result and one around all of it.
  `compare_given_boolean.py` reproduces
  `occt-boolean-given-preimplementation` (`rust_given_boolean_exists`
  false, keyed on `solid/boolean/curved/matched.rs`; the kernel's probe
  `unsupported` on all 36): every result valid with the reference's solid
  count, all 36 match (within 5.7e-9; exact frames 3.3e-15), no review; the
  rollex's cut has DRAW's area, 30152.95.
* **Kernel (S9e.2).** `tests/given_booleans.rs`: all 36 fixtures as the
  reference (30 within the kernel's enclosures, each at most `1e-9` wide;
  the 6 declared degenerate refused as tangencies), every second history
  complete over the given solid's ids and none of the other solid's, the
  first result's solids and the pick as the protocol's, results
  deterministic and moved rigidly, given solids translated with their third
  prisms keeping the reference's volumes, a partner strictly inside one
  solid of several `ComputationLimit`, the stored frames the reference's
  bit for bit. `compare_given_boolean.py`: 36 matches, every evaluated
  case with OCCT's unified entity counts, no review; every other Boolean
  comparison unchanged in its counts. DRAW's `bcut_simple/L3` to `L6`
  evaluate and are registered (`UPSTREAM_TESTS.md`).
* **S9e.3a evidence (given results of spheres, cones and tori, deeper
  chains, given results against a sphere, cone or torus), before its kernel
  code.** The protocol chains further Booleans (several `then` rows, each
  with its solid's rows). `chained_curved_boolean_reference.py` sweeps every
  face of every solid against all the others (S9d.4b.2's face sweeps,
  several solids and a set function; pushes off a face only for a solid
  with a surface on it), solids counted by rays at two resolutions.
  `generate_given_curved_boolean_fixtures.py --check` writes
  `boolean-given-curved-cases.txt`, `-expected.tsv` (its `expect` rows
  naming the class: `sphere`, `cone`, `torus`, `procedural`, `deep`,
  `partner`) and `-frames.tsv`: 48 cases (16 chains; 42 solid, 6
  degenerate): a domed box against a tilted slab and a cylinder through
  the dome's circle; DRAW's `G9` body with its rod moved clear, crossed by
  a box and bored coaxially; a countersunk box drilled; a box grooved by a
  torus against a tilted slab and a cylinder; a sphere and peg cut clear of
  their meeting; two deeper chains; a holed box against a sphere, a cone
  and a torus; `degenerate` `G9` itself and a box touching the dome.
  Checks: coaxial closed forms within 5.9e-41, two families of curves
  within 3.3e-35, the pair identities 4.1e-41, the area identity 1.5e-40,
  Monte Carlo 2.4 standard errors, no near coincidence.
  `compare_given_curved_boolean.py` reproduces
  `occt-boolean-given-curved-preimplementation`
  (`rust_given_curved_boolean_exists` false, keyed on
  `solid/boolean/curved/chain.rs`; the kernel's probe `unsupported` on all
  48): every result valid with the reference's solid count, 33 match, 15
  reviewed (BRepGProp's default integration on faces bounded by
  approximated meetings of two curved faces, up to 6.4e-6; adaptively
  within 2.4e-9); `G9`'s cut has DRAW's area 727.481.
* **Kernel (S9e.3a).** `tests/given_curved_booleans.rs`: all 48 fixtures
  as the reference (42 within the kernel's enclosures, each at most `1e-9`
  wide; the 6 declared degenerate refused as tangencies, `G9` among them),
  every stage's history complete and chaining, results deterministic and
  moved rigidly, given results translated with their last solids keeping
  the reference's volumes, a tree of four Booleans `ComputationLimit`, a
  given `Rise` met by a box S9e.3b's, the stored frames the reference's bit
  for bit. `compare_given_curved_boolean.py`: 25 matches and 23 reviewed
  (the native measures, and entity counts: the kernel's exact meeting
  pieces and whole periodic faces against OCCT's split approximations and
  seams); every other Boolean comparison unchanged in its counts. DRAW's
  `G9` and `H3` refused as the tangency they are (`UPSTREAM_TESTS.md`).
* **S9e.3b evidence (given meetings of two curved faces, and cones' and
  tori's general sections, met by the partner), before its kernel code.**
  `generate_given_met_boolean_fixtures.py --check` writes
  `boolean-given-met-cases.txt`, `-expected.tsv` (its `expect` rows naming
  the given edges' class, `rise`, `meet`, `toric`, `cone` or `spiric`, and
  the partner's surface; `meet N SINE` rows the triple points found) and
  `-frames.tsv`: 50 cases (17 chains; 45 solid, 5 degenerate): a sphere and
  a peg met by a tilted slab, a wall through the peg's axis, a pipe, a ball
  and a frustum; two crossing rods met by a slab, a wall and a ball; a
  frustum pierced by a pipe met by a wall; a torus with a rod through its
  tube met by a slab and a wall; a frustum cut obliquely and a torus cut off
  its axis met by rods, pipes and balls; `degenerate` a box touching the
  peg's meeting at its lowest point and a ball there with dependent
  normals. Rows from S9e.3a's chained reference (two families within
  1.1e-35, closed forms 1.6e-36, the pair identities 9.2e-41, the area
  identity 1.8e-40, Monte Carlo 3.2 standard errors);
  `given_met_reference.py` traces the meetings independently (34 triple
  points on the given edges inside the partner's faces, every sine at least
  0.38; the declared tangencies' normals dependent exactly).
  `compare_given_met_boolean.py` reproduces
  `occt-boolean-given-met-preimplementation`
  (`rust_given_met_boolean_exists` false, keyed on
  `solid/boolean/curved/triple.rs`; the kernel's probe `unsupported` on all
  50): every result valid with the reference's solid count, 12 match, 38
  reviewed (BRepGProp's default integration on approximated sections, up to
  4.8e-6; adaptively within 3.2e-8).
* **S9e.4 evidence (imported solids), before its kernel code.** 13 bodies
  written by OCCT itself (`occt_boolean_oracle.cpp`'s `write` blocks,
  `BRepPrimAPI` makers and `MakePrism`, `BRepTools::Write` in version 1)
  under `rust/fixtures/imported/`; the case protocol's `brep PATH` row.
  `generate_imported_boolean_fixtures.py --check` writes
  `boolean-imported-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 69 cases (23 groups; 57 solid, 9 degenerate, 3
  unsupported): boxes (one tilted), cylinders (one along `x`), a profile of
  lines and a tangent arc, a cylinder of two half faces, a plate with a
  hole, a sphere, a hemisphere, a frustum, a cone with its apex and a torus
  against slabs, rods, boxes, a sphere and another imported body, as object
  and tool, and one chain. Rows from the constructions OCCT was given:
  S9e.3a's chained reference, S9d.1's sphere reference and coaxial closed
  forms (closed forms within 4.3e-41, two families 8.6e-32, pair identities
  2.3e-41, Monte Carlo 2.4 standard errors); every file read independently
  with its locations, each stored vertex on its construction's surfaces.
  `compare_imported_boolean.py` reproduces
  `occt-boolean-imported-preimplementation`
  (`rust_imported_boolean_exists` false, keyed on `solid/imported.rs`; the
  probe reading every file and reporting `unsupported`): every result
  valid, 61 match, 8 reviewed (BRepGProp's default integration on the
  crossing cylinders' and the torus's approximated sections, up to 2.8e-5,
  adaptively within 1.7e-9; the declared tangent fuses OCCT keeps as two
  solids touching along a line).
* **Kernel (S9e.4a).** `tests/imported_booleans.rs`: all 69 fixtures as
  the reference (57 within the kernel's enclosures, each at most `1e-9`
  wide; the 9 declared degenerate refused as tangencies, the 3 turned
  profiles `OutOfDomain` as S9e.4b's), every history complete over the
  imported bodies' stored ids (no relation from a construction's), results
  deterministic and moved rigidly, both inputs translated keeping the
  reference's volumes, every body its construction (its closed-form volume,
  points classified, its stored vertices on its boundary, two imports'
  ids apart), STEP solids imported and cut alike (a box, cylinders, a
  sphere, a hemisphere, cones, a torus, an L prism, a plate with a hole; a
  cavity and a spline prism refused), the stored frames the reference's bit
  for bit. `compare_imported_boolean.py`: 54 matches and 15 reviewed (the
  native measures, and entity counts: the kernel's whole periodic faces and
  exact meeting pieces against OCCT's seams and split approximations, the
  imported seam-split cylinder's two faces kept).
* **S9e.4b.1 evidence (imported prisms whose arcs round off their
  circles), before its kernel code.** 5 bodies written by OCCT
  (`MakePrism` of profiles in turned frames: a stadium in `TILT`, a
  rectangle with four fillets in `R125`, a circle as two arcs, an arc
  crossing its lines and a lens of two circles in `TURN30`) under
  `rust/fixtures/imported/`. `generate_imported_arcs_boolean_fixtures.py
  --check` writes `boolean-imported-arcs-cases.txt`, `-expected.tsv`,
  `-frames.tsv` and `-bodies.txt`: 36 cases (12 groups; 27 solid, 6
  degenerate, 3 unsupported) against boxes, a rod, a ball, a slab and
  another imported body, as object and tool, and one chain, from S9e.3a's
  chained reference on the constructions OCCT was given (two families
  4.1e-32, pair identities 4.4e-38, Monte Carlo 3.2 standard errors),
  every file read independently (stored vertices on the constructions'
  surfaces within 3.9e-17 of the size, the arcs' ends off their circles
  once rounded into the construction's frame).
  `compare_imported_arcs_boolean.py` reproduces
  `occt-boolean-imported-arcs-preimplementation` (keyed on
  `curved/snapped.rs`, the probe `unsupported` on all 36 before it):
  every result valid, 30 match, 6 reviewed (BRepGProp's default
  integration on the ball's and crossing cylinders' approximated sections,
  up to 1.4e-5, adaptively within 3.2e-9; the declared tangent fuse OCCT
  keeps as two solids touching along a line).
* **Kernel (S9e.4b.1).** `tests/imported_arc_booleans.rs`: all 36
  fixtures as the reference (27 within the kernel's enclosures, each at
  most `1e-9` wide; the 6 declared degenerate refused, the coplanar box as
  two faces within the resolution of one plane and the touching box as a
  tangency; the lens `OutOfDomain` as S9e.4b.4's), every history complete
  over the imported bodies' stored ids, results deterministic and moved
  rigidly, both inputs translated (and turned where no face of one is
  exactly parallel to the other's cylinder) keeping the reference's
  volumes, every body its construction. `compare_imported_arcs_boolean.py`:
  27 matches and 9 reviewed (the native measures, the tangent fuse, and
  entity counts: the ball's meeting with the slot's wall in exact pieces,
  the split circle's two faces kept). S9e.4a's `dee_turn` cases, refused
  before as S9e.4b's, are degenerate under S9's rules (the turned
  profile's corner on the box's), declared so.
* **S9e.4b.2 evidence (imported polyhedra other than prisms), before its
  kernel code.** 15 bodies written by OCCT under `rust/fixtures/imported/`
  (the oracle's `wedge` and `polyhedron` rows and a `write` block of a
  Boolean): a tetrahedron, an octahedron, a pyramid, a pyramid's frustum
  and a slanted wedge in turned frames, a box less a skew box, a box fused
  with a turned box, a box with a cavity, and the DRAW survey's shapes from
  exact points (two frustums sharing a face, a frustum on a box, a partly
  drafted prism with a prism on its cap, two frustums at right angles).
  `generate_imported_polyhedra_boolean_fixtures.py --check` writes
  `boolean-imported-polyhedra-cases.txt`, `-expected.tsv`, `-frames.tsv`
  and `-bodies.txt`: 48 cases (16 groups; 40 solid, 4 empty, 1 degenerate,
  3 unsupported) against boxes, slabs, rods and another imported body, as
  object and tool, and two chains, from an independent exact reference of
  convex cells (`imported_polyhedra_boolean_reference.py`: volumes three
  ways and the pair identities exactly, areas two ways within 1.8e-40,
  closed forms exactly, Monte Carlo 3.3 standard errors, vertex and edge
  clearances of at least 4.0e-3 of the size outside 196 exact contacts),
  every file read independently (stored vertices on the construction's
  boundary within 6.6e-16 of the size, no body a prism but the ridge).
  `compare_imported_polyhedra_boolean.py` reproduces
  `occt-boolean-imported-polyhedra-preimplementation` (keyed on
  `polyhedra/imported.rs`, the probe `unsupported` on all 48 before it):
  every result valid, 47 match (within 8.9e-16), 1 reviewed (the declared
  flush fuse, one solid sharing the pyramid's base within OCCT's
  tolerance).
* **Kernel (S9e.4b.2).** `tests/imported_polyhedra_booleans.rs`: all 48
  fixtures as the reference (40 within the kernel's enclosures, each at
  most `1e-9` wide, and 4 empty; the flush fuse refused as a face using an
  edge both ways; the cavity `OutOfDomain` as S9e.4b.4's), every history
  complete over the imported bodies' stored ids, results deterministic and
  moved rigidly, both inputs translated and turned keeping the reference's
  volumes, every body its construction, curved partners refused as
  S9e.4b.4's. `compare_imported_polyhedra_boolean.py`: 47 matches and 1
  reviewed, every entity count OCCT's after unifying (an exactly planar
  imported face joined with a partner's on its plane). A trial of the DRAW
  survey's 7 restore cases with such polyhedra: 4 evaluate on both
  backends with every check, 3 refused by S9's rules (a fuse of frustums
  apart by more than the resolution, two solids; a shared corner stored
  1e-13 apart).
* **S9e.4b.3a evidence (imported plane pieces), before its kernel code.**
  8 bodies written by OCCT under `rust/fixtures/imported/` (a `write`
  block of a sphere common or less a box or prism, on turned rational
  frames): two octants, a wedge above a parallel's plane, a lune, a half,
  a zone's wedge, and S9e.4b.3c's octant below a parallel's plane and
  bitten ball. `generate_imported_pieces_boolean_fixtures.py --check`
  writes `boolean-imported-pieces-cases.txt`, `-expected.tsv`,
  `-frames.tsv` and `-bodies.txt`: 45 cases (15 groups; 33 solid, 6
  degenerate, 6 unsupported) against boxes, a slab, rods, a ball and a
  cone, as object and tool, two pieces, a piece and an imported box, and a
  chain, from S9e.3a's chained reference on the constructions OCCT was
  given (two families within 7.5e-37, closed forms of each solid and each
  piece within 2.0e-40, the pair and area identities, Monte Carlo 3.5
  standard errors), every file read independently (stored vertices on the
  construction's surfaces within 7.6e-16 of the size).
  `compare_imported_pieces_boolean.py` reproduces
  `occt-boolean-imported-pieces-preimplementation` (keyed on
  `curved/pieces.rs`, the probe `unsupported` on all 45 before it): every
  result valid, 35 match, 10 reviewed (BRepGProp's default integration on
  the rods' and the cone's meetings with a sphere, adaptively within
  3.2e-9; the declared touching fuse, two solids sharing a point).
* **Kernel (S9e.4b.3a).** `tests/imported_piece_booleans.rs`: all 45
  fixtures as the reference (33 within the kernel's enclosures, each at
  most `1e-9` wide; the 6 declared degenerate refused, the turned box on
  the octant's base plane as two faces within the resolution of one plane
  and the touching box as a tangency; the 6 one-sphere and bitten cases
  `OutOfDomain` as S9e.4b.3c's), every history complete over the imported
  bodies' stored ids (none from a primitive or a hull), results
  deterministic and moved rigidly, both inputs translated and turned
  keeping the reference's volumes (the piece read again off its moved
  stored topology), every body a piece (its volume its closed form, a point
  inside it and its mirror in the centre outside, its stored vertices on
  its boundary), the kernel's own split pieces of a cylinder, a frustum and
  a sphere's zone imported as pieces, their Booleans with a box and a ball
  obeying the pair identities, the octant against a torus, the stored
  frames the reference's bit for bit.
  `compare_imported_pieces_boolean.py`: 32 matches and 13 reviewed (the
  native measures, the touching fuse, and entity counts: each splits its
  sections at its own points and seams). A trial of the DRAW survey's 39
  S9e.4b.3 restore cases: none evaluates, 38 pairs of pieces of one sphere
  S9e.4b.3c's and `buc60926` a plane through a cone's apex.
* **S9e.4b.3b evidence (split pieces), before its kernel code.** A
  `split` row in both protocols (the piece `Solid::split_by_plane` leaves
  on one side; natively the solid common `BRepPrimAPI_MakeHalfSpace`).
  `generate_split_pieces_boolean_fixtures.py --check` writes
  `boolean-split-pieces-cases.txt`, `-expected.tsv` and `-frames.tsv`: 48
  cases (16 groups; 39 solid, 6 degenerate, 3 unsupported) on 8 split
  pieces (a cylinder's, a dee's and a box's oblique pieces, a frustum's and
  a zone's across their walls, a zone's half, a torus's band, a frustum's
  half through its axis) against boxes, balls, rods, a cone, each other and
  in a chain, from S9e.3a's chained reference on the constructions OCCT is
  given (each piece its solid common a half-space box, a zone's its whole
  sphere and the slab between its parallels; two families within 6.2e-36,
  closed forms of each solid and the cylinder's piece, the pair and area
  identities, Monte Carlo 2.8 standard errors).
  `compare_split_pieces_boolean.py` reproduces
  `occt-boolean-split-pieces-preimplementation` (keyed on
  `curved/splits.rs`, the probe `unsupported` on all 48 before it): every
  result valid, 33 match, 15 reviewed (BRepGProp's default integration on
  the balls' and rods' meetings with the pieces, adaptively within 1.4e-8).
* **Kernel (S9e.4b.3b).** `tests/split_piece_booleans.rs`: all 48 fixtures
  as the reference (39 within the kernel's enclosures, each at most `1e-9`
  wide; the frustum's half through its axis refused at its apex, the box on
  the cut plane as two faces within the resolution of one plane, the zone's
  half against its own ball `OutOfDomain` as S9e.4b.3c's), every history
  complete over the pieces' ids, results deterministic and moved rigidly,
  moved inputs keeping the reference's volumes, the stored frames and
  planes the reference's bit for bit, a U profile's two pieces above one
  plane, a hemisphere's halves through its pole, pieces against an imported
  piece and a given result, a torus's spiric piece refused, and three
  pieces the split target found (a zone's ring, a cap's half's hole, a
  frustum's ring), split and imported. `compare_split_pieces_boolean.py`:
  26 matches and 22 reviewed (the native measures, and entity counts: each
  splits its sections at its own points and seams).
* **S9e.4b.3c.1 evidence (two pieces of one sphere), before its kernel
  code.** A `write` block's row `divide` in the native oracle
  (`ShapeUpgrade_ShapeDivideClosedEdges`: a rim two arcs, as the DRAW
  survey's `so1` and `so4`). `generate_one_sphere_boolean_fixtures.py
  --check` writes `boolean-one-sphere-cases.txt`, `-expected.tsv`,
  `-frames.tsv` and `-bodies.txt`: 22 cases (9 groups; 16 solid, 1 empty,
  2 degenerate, 3 unsupported) on five pieces of one sphere OCCT wrote (a
  hemisphere and a cap with divided rims, two octants, a wedge above a
  parallel's plane) against each other, a ball of their sphere and in a
  chain, from S9e.3a's chained reference on the constructions OCCT was
  given (two families within 3.4e-36, each piece's closed form, the pair
  identities, Monte Carlo 2.4 standard errors, solid counts by rays joined
  within two grid spacings).
  `compare_one_sphere_boolean.py` reproduces
  `occt-boolean-one-sphere-preimplementation` (keyed on the refusal of two
  inputs on one sphere in `curved/pieces.rs`, the probe `unsupported` on all
  22 before it): 8 match, 14 reviewed (10 results OCCT leaves invalid on
  same-domain spherical faces, their measures adaptively within 7.1e-12; a
  common by BRepGProp's default integration; the chain OCCT gets wrong).
* **Kernel (S9e.4b.3c.1).** `tests/one_sphere_booleans.rs`: all 22
  fixtures as the reference (17 within the kernel's enclosures, each at
  most `1e-9` wide; the ball within the resolution of the cap's sphere
  refused as `Degenerate`, the hemisphere and the octant on one frame
  `OutOfDomain` as S9e.4b.3c.2's), every history complete over the imported
  bodies' stored ids, results deterministic and moved rigidly, moved
  inputs keeping the reference's volumes, the divided rims' pieces of two
  arcs, two whole balls of one sphere, and a groove refused on import as
  S9e.4b.3c's. S9e.4b.3a's and S9e.4b.3b's one-sphere cases are solid
  within the reference. `compare_one_sphere_boolean.py`: 8 matches and 14
  reviewed (the native measures and validity, OCCT's wrong chain, and
  entity counts: the kernel joins the pieces of both inputs on the sphere
  facing one way). A trial of the DRAW survey's 38 restore cases of `so1` to
  `so7`: 14 evaluate on both backends, 13 S9e.4b.3c.2's, 11 `Degenerate`
  (the turned copies' corners within rounding of the partner's).
* **S9e.4b.3c.2 evidence (exact incidences on one sphere), before its
  kernel code.** `generate_one_sphere_incidence_boolean_fixtures.py
  --check` writes `boolean-one-sphere-incidence-cases.txt`,
  `-expected.tsv`, `-frames.tsv` and `-bodies.txt`: 45 cases (19 groups;
  29 solid, 4 empty, 12 degenerate) on eleven pieces of one ball OCCT wrote
  on the world's axes (`MakeBox` on frames whose normal is the axis and
  whose origin lies on it: the hemisphere with its rim divided, wedges as
  the DRAW survey's `so2`, `so3` and `so5`, a half, two octants, wedges
  about another axis, two near copies) against each other, caps of the
  ball on one frame and in a chain, from an exact reference of its own
  (the ball cut into cells by every input's heights along the axis and
  half-planes about it, closed forms per cell; the volume again by the
  divergence theorem within 4.7e-41, the pair identities, every
  Monte-Carlo point's membership that of its cell, Monte Carlo 2.5 standard
  errors), every stored plane read back holding the axis, normal to it or
  apart from the ball exactly.
  `compare_one_sphere_incidence_boolean.py` reproduces
  `occt-boolean-one-sphere-incidence-preimplementation` (keyed on the
  refusal of exact incidences on one sphere in `curved/graph.rs`, the probe
  `unsupported` on the 36 solid, empty and `half_wedge` cases before it,
  the 9 near and touching ones refused by S9's near-plane guard): 42 match,
  3 reviewed (the hemisphere less the wedge invalid with the reference's
  measures; the near cases' slivers below OCCT's tolerance, dropped).
* **Kernel (S9e.4b.3c.2).** `tests/one_sphere_incidence_booleans.rs`: all
  45 fixtures as the reference (33 solid and empty within the kernel's
  enclosures, each at most `1e-9` wide; the 12 degenerate refused with S9's
  reasons, none as an incidence), every history complete over the imported
  bodies' stored ids, an edge or vertex of both continuing both inputs'
  (the wedges' axis edges, corner and pole), results deterministic and
  moved rigidly, both inputs moved by exact motions (a dyadic translation,
  a quarter turn) keeping the reference's volumes and turned by a rounding
  rotation refused or within it. S9e.4b.3c.1's `hemi_octant` is solid
  within its reference. `compare_one_sphere_incidence_boolean.py`: 42
  matches and 3 reviewed; `compare_one_sphere_boolean.py` 8 and 14
  unchanged. A trial of the DRAW survey's restore cases: the 13 of `so1`
  and `so2`, `so2` and `so3`, `so5` and `so2` evaluate on both backends
  with every check, 27 of the 38 `so` cases in all, the other 11 (`so6`
  and `so7`) `Degenerate`.
* **S9e.4b.3c.3a evidence (plane pieces of another Boolean), before its
  kernel code.** `generate_piece_forms_boolean_fixtures.py --check` writes
  `boolean-piece-forms-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 39 cases (13 groups; 33 solid, 6 degenerate) on nine
  bodies OCCT wrote, every section a circle or a line (`form_boss`, the
  survey's `bcut_complex/G4` part on the world's axes: a box with a
  cylindrical boss; a scoop with its top face in two faces on one plane, a
  slot ending at its floor, a dimple, a ball's boss, a conical hole, a
  cylinder bitten between its caps; declared `degenerate` `I6`'s notched
  tool and `shading_132`'s three-quarter frustum) against rods, slabs,
  balls and boxes, as the tool, two imported, in a chain and in `G4`'s
  configuration, from S9e.3a's chained reference on the constructions OCCT
  was given (two families within 3.4e-36, each body's closed form 3.7e-38,
  the pair identities, Monte Carlo 2.6 standard errors, solid counts by
  rays joined within two grid spacings).
  `compare_piece_forms_boolean.py` reproduces
  `occt-boolean-piece-forms-preimplementation` (keyed on `imported.rs`'s
  refusal of a piece other than its primitive common its planes, the probe
  `unsupported` on the 36 solid and notch cases before it, the frustum's 3
  refused at its apex): 26 match, 13 reviewed (BRepGProp's default
  integration, within 1.4e-9 adaptively).
* **Kernel (S9e.4b.3c.3a).** `tests/piece_form_booleans.rs`: all 39
  fixtures as the reference (33 within the kernel's enclosures, each at
  most `1e-9` wide; the notch refused for its wall tangent to its own
  faces, the frustum for its planes through its apex), every history
  complete over the imported bodies' stored ids, faces on one plane keeping
  their ids, results deterministic and moved rigidly, moved inputs keeping
  the reference's volumes, every body in its form, and the kernel's own
  grooves, bosses and bites written, read back and imported. S9e.4b.3a's
  bitten ball is a bite (`bitten_box` solid within its reference).
  `compare_piece_forms_boolean.py`: 20 matches and 19 reviewed (the native
  measures and, with the kernel, entity counts: OCCT splits its cylinders'
  and spheres' faces at its seams); `compare_imported_pieces_boolean.py` 30
  and 15. A trial of the DRAW survey's restore cases: `bcut_complex/G4`
  evaluates on both backends with every check, `I6` is refused for its
  tool's own tangency (its part tangent to the tool besides),
  `buc60926` at its frustum's apex as before.
* **S9e.4b.3c.3b evidence (plane pieces of Boolean trees), before its
  kernel code.** `generate_piece_trees_boolean_fixtures.py --check` writes
  `boolean-piece-trees-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 30 cases (10 groups; 24 solid, 3 degenerate, 3
  unsupported) on eight bodies OCCT wrote, every section a circle or a line
  (a U prism's groove by a ball, the fuzz target's; a U's cylindrical boss;
  a cylinder bitten twice by a U; a box fused with or less an upper
  hemisphere whose disc lies below or inside it, the fuzz target's ball's
  half; a box with a cylinder flattened along its axis; declared
  `degenerate` `I6`'s notch on a U and `unsupported` a tooth left in a
  cylinder's U-shaped bite) against rods, slabs, a ball and a box, as the
  tool, two imported and in a chain, from S9e.3a's chained reference on the
  constructions OCCT was given (a U as a box less its slot's box, a
  hemisphere as its ball common a cylinder above its equator; two families
  within 1.6e-35, each body's closed form 4.5e-37, the pair identities,
  Monte Carlo 2.5 standard errors, solid counts by rays joined within two
  grid spacings). `compare_piece_trees_boolean.py` reproduces
  `occt-boolean-piece-trees-preimplementation` (keyed on `imported.rs`'s
  refusal of a piece other than one Boolean of its primitive and its
  planes' hull, the probe `unsupported` on all 30 before it): 27 match, 3
  reviewed (BRepGProp's default integration, within 4.8e-9 adaptively).
* **Kernel (S9e.4b.3c.3b).** `tests/piece_tree_booleans.rs`: all 30 fixtures
  as the reference (24 within the kernel's enclosures, each at most `1e-9`
  wide; the U's notch refused for its wall tangent to its own faces, the tooth
  as S9e.4b.4's), every history complete over the imported bodies' stored ids,
  the U's arm ends on one plane keeping their ids, results deterministic and
  moved rigidly, moved inputs keeping the reference's volumes, every body its
  tree, and the kernel's own bodies of trees imported (the fuzz target's two,
  a rod bitten twice, a cylinder's oblique end fused with a box).
  `compare_piece_trees_boolean.py`: 22 matches and 8 reviewed (the native
  measures and, with the kernel, entity counts: its meetings of a cylinder
  with a ball split at their turning points, OCCT's unified result keeping a
  boss's seam edge). A trial of the DRAW survey's 171 restore cases the import
  reaches: no case moves (55 evaluate on both backends), none refused as
  S9e.4b.3c's.
* **S9e.4b.4a evidence (imported prisms whose arcs of two circles meet at a
  joint), before its kernel code.** 5 bodies written by OCCT (`MakePrism` of
  profiles in turned frames, every joint of two circles a rational common
  point of both in the profile's frame: four discs' common, a pointed arch,
  an arc tangent inside another, two circles crossing at 0.11 rad, and a
  lens whose arcs are split) under `rust/fixtures/imported/`.
  `generate_imported_joints_boolean_fixtures.py --check` writes
  `boolean-imported-joints-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 39 cases (13 groups; 33 solid, 3 degenerate, 3
  unsupported) against boxes, rods, a ball and a `TILT` slab, as object and
  tool, two imported and a chain, from S9e.3a's chained reference on the
  constructions OCCT was given (two families 2.1e-36, pair identities
  3.4e-42, Monte Carlo 3.3 standard errors), every file read independently
  (stored vertices on a wall and a cap within 1.9e-15 of the size, all 18
  stored joints off their circles once rounded into the construction's
  frame). The S curve of two arcs touching from either side (`P4`'s) is not
  convex, which the reference does not take: the kernel's tests write one.
  `compare_imported_joints_boolean.py` reproduces
  `occt-boolean-imported-joints-preimplementation` (keyed on `snapped.rs`'s
  refusal of arcs of two circles meeting at a joint, the probe
  `unsupported` on all 39 before it): every result valid, 36 match, 3
  reviewed (BRepGProp's default integration on a ball's faces met by the
  quad's walls, up to 7.0e-7, adaptively within 2.3e-9).
* **Kernel (S9e.4b.4a).** `tests/imported_joint_booleans.rs`: all 39
  fixtures as the reference (33 within the kernel's enclosures, each at most
  `1e-9` wide; `quad_seat`, a rod on one of the quad's stored circles below
  it, refused as two cylinders within the resolution of one cylinder; the
  split lens `OutOfDomain` as S9e.4b.4's), every history complete over the
  imported bodies' stored ids, results deterministic and moved rigidly,
  moved inputs keeping the reference's volumes, every body its construction,
  and the kernel's own four discs' common, lens, S curve and internal
  fillet written, read back and imported giving its own Booleans' volumes.
  `compare_imported_joints_boolean.py`: 36 matches and 3 reviewed (the
  native measures and, with the kernel, entity counts: a ball's meetings
  and sphere face split at the kernel's own points and at OCCT's).
  S9e.4b.1's lens cases are solid within their reference
  (`compare_imported_arcs_boolean.py` 27 and 9, all 30 solid cases within
  the reference). A trial of the DRAW survey's 171 restore cases the import
  reaches: only the 6 of two circles at a joint move, `bcut_complex/P4`
  evaluating on both backends with every check, `E8`, `D5` and `E1`
  refused as two cylinders within the resolution of one cylinder (each
  tool on a stored circle of the part's arcs), `bug4993_1` and `_2` as two
  faces within the resolution of one plane.
* **S9e.4b.4b.1 evidence (imported bodies of several primitives), before
  its kernel code.** 8 bodies written by OCCT, each one Boolean of two
  coaxial primitives, every section a circle (E5's stepped shaft, a cup with
  a blind bore, G9's dome and pin on one ball, a bead, a knob; declared
  `degenerate` a capsule, tangent along its rim, given to no case for the
  reference's cost there, and declared `unsupported` K1's rounded box less a
  cylinder) under `rust/fixtures/imported/`.
  `generate_primitive_chains_boolean_fixtures.py --check` writes
  `boolean-primitive-chains-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 24 cases (8 groups; 21 solid, 3 unsupported) against a
  ball, `TILT` boxes and slabs, as object and tool, the dome and the pin
  both imported and a chain, from S9e.3a's chained reference on the
  constructions OCCT was given (two families 1.8e-35, each body's closed
  form 4.3e-42, Monte Carlo 2.7 standard errors), every stored vertex
  within 2.3e-14 of the size on its construction.
  `compare_primitive_chains_boolean.py` reproduces
  `occt-boolean-primitive-chains-preimplementation` (keyed on
  `imported.rs`'s general refusal's text before the step, the probe
  `unsupported` on all 24 before it): every result valid, 20 match, 4
  reviewed (the cup's ball cases, OCCT's faces bounded by its approximated
  quartics of tolerance up to 8.0e-6; the knob's common, BRepGProp's
  default integration).
* **Kernel (S9e.4b.4b.1).** `tests/primitive_chain_booleans.rs`: all 24
  fixtures as the reference (21 within the kernel's enclosures, each at most
  `1e-9` wide; the rounded box `OutOfDomain` as S9e.4b.4b.2's), every
  history complete over the imported bodies' stored ids, results
  deterministic and moved rigidly, moved inputs keeping the reference's
  volumes, every body its chain, the capsule refused for its own tangency,
  and the kernel's own stepped shaft, cup and comb written, read back and
  imported giving its own Booleans' volumes.
  `compare_primitive_chains_boolean.py`: 14 matches and 10 reviewed (the
  native measures and, with the kernel, entity counts: its meetings of a
  cylinder with a ball split at their turning points, OCCT's seams). A
  trial of the targeted DRAW restore cases: `bfuse_complex/E5`,
  `bcut_complex/G9` and `bugs/modalg_2/bug417` evaluate on both backends
  with every check; `bugs/modalg_6/bug28773` is refused as two cylinders'
  axes within rounding of parallel (its tube's disc frame leaning 2.2e-33
  off its walls), `bfuse_complex/K1` as S9e.4b.4b.2's.
* **S9e.4b.4b.2a evidence (imported bodies led by a prism leaf), before its
  kernel code.** 7 bodies written by OCCT, each one Boolean of a prism of
  lines and arcs (a square its corners rounded, or a stadium) and a
  primitive, every section a circle or a line (K1's rounded box less a bore
  on the world's axes, a stadium plate with a boss, a rounded plate with a
  dimple, a stadium plate with a conical pocket, a rounded plate with a
  dome; declared `degenerate` a post on a fillet's circle and declared
  `unsupported` a notch through both caps at the rim, both given to no case)
  under `rust/fixtures/imported/`. `generate_prism_leaves_boolean_fixtures.py
  --check` writes `boolean-prism-leaves-cases.txt`, `-expected.tsv`,
  `-frames.tsv` and `-bodies.txt`: 24 cases (8 groups; 21 solid, 3
  degenerate: a rod tangent to the bore) against `TILT` rods and boxes, a
  ball, as object and tool, the boss and the pocket both imported and the
  dimple's chain, from S9e.3a's chained reference on the constructions OCCT
  was given (two families 3.8e-36, each body's closed form 4.0e-42, Monte
  Carlo 2.5 standard errors), every stored vertex within 6.7e-16 of the size
  on its construction. `compare_prism_leaves_boolean.py` reproduces
  `occt-boolean-prism-leaves-preimplementation` (keyed on S9e.4b.4b.1's
  refusal of plane faces other than the primitives' ends, the probe
  `unsupported` on all 24 before the step): every result valid, 19 match, 5
  reviewed (BRepGProp's default integration).
* **Kernel (S9e.4b.4b.2a).** `tests/prism_leaf_booleans.rs`: all 24
  fixtures as the reference (21 within the kernel's enclosures, each at most
  `1e-9` wide; the tangent rod `Degenerate`), every history complete over
  the imported bodies' stored ids, results deterministic and moved rigidly,
  moved inputs keeping the reference's volumes, every body its chain led by
  its prism, the post refused for its own tangency and the notch as
  S9e.4b.4b.2b's, and the kernel's own boss, dimple and dome plates written,
  read back and imported giving its own Booleans' volumes (K1's shape the
  writer refuses: its bore's rims stored as ellipses of equal axes).
  `compare_prism_leaves_boolean.py`: 15 matches and 9 reviewed (the native
  measures and, with the kernel, entity counts: its meetings of two
  cylinders split at their turning points, OCCT's seams and its closed
  edges' vertices); S9e.4b.4b.1's rounded box solid, its
  `compare_primitive_chains_boolean.py` 14 and 10 with the kernel within
  the reference on all 24. A trial of the targeted DRAW restore cases:
  `bfuse_complex/K1`'s part imports as its chain and the fuse is refused as
  a tangency between the inputs (its tool a rod of the bore's radius whose
  axis crosses it); `bugs/modalg_6/bug28773` stays two cylinders' axes
  within rounding of parallel; `E5`, `G9` and `bug417` evaluate as before.
* **S9e.4b.4b.2b.1 evidence (imported bodies with a primitive's flat, a
  prism of one cap or a pocket), before its kernel code.** 4 bodies written
  by OCCT, each one Boolean, every section a circle or a line (a ball's half
  fused with a frustum standing on its disc off its axis, a hexagon under a
  stadium, flush, a stadium boss on a box read off its top, a plate with a
  hole less its copy's band turned half a turn, a slot holding a crescent)
  under `rust/fixtures/imported/part_*.brep`, the `boolean` fuzz target's
  three refused first results' kinds. `generate_plane_parts_boolean_fixtures.py
  --check` writes `boolean-plane-parts-cases.txt`, `-expected.tsv`,
  `-frames.tsv` and `-bodies.txt`: 24 cases (8 groups; 21 solid, 3
  degenerate: a ball resting on the frustum's top) against rods, a box and a
  ball, as object and tool, the cake and the stack both imported and the
  cake's chain, from S9e.3a's chained reference on the constructions OCCT
  was given (two families 1.4e-31, each body's closed form 6.2e-33, Monte
  Carlo 2.2 standard errors), every stored vertex within 4.7e-16 of the size
  on its construction. `compare_plane_parts_boolean.py` reproduces
  `occt-boolean-plane-parts-preimplementation` (keyed on S9e.4b.4b.2a's
  refusal of plane faces other than the primitives' ends or a prism's, the
  probe `unsupported` on all 24 before the step, the cake's refusal as its
  own tangency counted so): every result valid, 19 match, 5 reviewed
  (BRepGProp's default integration; a touching ball's fuse OCCT keeps as two
  solids).
* **Kernel (S9e.4b.4b.2b.1).** `tests/plane_part_booleans.rs`: all 24
  fixtures as the reference (21 within the kernel's enclosures, each at most
  `1e-9` wide; the touching ball `Degenerate`), every history complete over
  the imported bodies' stored ids, results deterministic and moved rigidly,
  moved inputs keeping the reference's volumes, every body its chain with its
  flat, prism of one cap or pocket, S9e.4b.4b.2a's notch refused as
  S9e.4b.4b.2b.2's, and the kernel's own cake and stack and the fuzz
  target's three bodies written, read back and imported (their Booleans with
  a turned box the kernel's own results', but the kernel's stack, which
  given to a Boolean it refuses: `an edge of one input on a face of the
  other`). `compare_plane_parts_boolean.py`: 15 matches and 9 reviewed (the
  native measures and, with the kernel, entity counts: OCCT's seams and edge
  divisions, its meetings of a ball with a cylinder approximated where the
  kernel's split at their turning points). A trial of the 171 DRAW restore
  cases the import reaches: no case moves (59 evaluating on both backends,
  `bfuse_complex/K1` a tangency between the inputs as S9e.4b.4b.2a's trial
  found).
* **Equal cylinders with crossing axes (S9c.1's tangent points).** A
  plate's round hole against a rod of its radius whose axis crosses the
  hole's at a right angle (`K1`'s part bore and tool rod, built or as a box
  already cut by a rod) gave an `InvalidTopology` fuse, a cut and a
  `Degenerate` common. The two points where the equal cylinders' ellipses
  cross, the surfaces tangent there, are now a tangency between the inputs
  where they lie inside both faces with the faces' outward normals opposite
  (the inputs touching each other there), every operation `Degenerate` as
  a rod tangent to the hole at one point is; whole rods (normals alike)
  keep their Steinmetz fuse and common, their cut touching itself. A result
  found touching itself at a vertex is `Degenerate` on every path
  (REVIEW_NOTES.md's "Equal cylinders with crossing axes").
  `tests/curved_booleans.rs`'s `equal_crossing_cylinders_are_a_tangency`;
  S9c.1's fixtures and `compare_curved_boolean.py` (42 matches, 2
  reviewed) unchanged, and the upstream whole-rod cases
  `bfuse_complex/J5`, `bopfuse_simple/ZD8` and `ZE1` and
  `bopcommon_simple/ZD8` and `ZE1` still evaluating.
* **S9e.4b.4c.1 evidence (imported polyhedra against curved faces and with
  cavities), before its kernel code.** S9e.4b.2's bodies (a pyramid, a
  frustum of one, a slanted wedge, a tetrahedron, an octahedron, a notched
  box, a hollow box) and S9e.4a's ball read again, and one body OCCT wrote (a
  box in a turned frame less a box in another inside it: a cavity whose faces
  are two triangles each) under `rust/fixtures/imported/poly_cavity.brep`.
  `generate_polyhedra_curved_boolean_fixtures.py --check` writes
  `boolean-polyhedra-curved-cases.txt`, `-expected.tsv`, `-frames.tsv` and
  `-bodies.txt`: 34 cases (12 groups; 30 solid, 3 degenerate: a ball
  tangent to the tetrahedron's base, 1 unsupported: a ball inside the hollow
  box's cavity fused, a cavity among several solids) against rods, balls
  and a frustum, as object and tool, the tetrahedron and the imported ball,
  two chains and the cavities' ball and slab, from
  `polyhedra_curved_boolean_reference.py` (S9e.3a's chained reference with
  convex hulls of exact points; two families 4.2e-36, each body's closed
  form 1.0e-41, Monte Carlo 2.7 standard errors), every stored vertex within
  2.2e-16 of the size on its construction's planes.
  `compare_polyhedra_curved_boolean.py` reproduces
  `occt-boolean-polyhedra-curved-preimplementation` (keyed on
  `polyhedra/imported.rs`'s refusal of an imported polyhedron against curved
  faces, the probe `unsupported` on all 34 before the step): every result
  valid, 32 match, 2 reviewed (BRepGProp's default integration; the
  tetrahedron's common with the ball no closer than its edges' tolerance).
* **Kernel (S9e.4b.4c.1).** `tests/polyhedra_curved_booleans.rs`: all 34
  fixtures as the reference (30 within the kernel's enclosures, each at most
  `1e-9` wide; the touching ball `Degenerate`, the inner ball's fuse a cavity
  among several solids), every history complete over the imported bodies'
  stored ids, results deterministic and moved rigidly, moved inputs keeping
  the reference's volumes, every body its stored triangles (its cavity
  outside it), the kernel's own notched and hollow boxes written, read back
  and imported (their Booleans with a rod and a ball the kernel's own
  results'), and a cavity split in two by a slab given to a ball.
  `compare_polyhedra_curved_boolean.py`: 18 matches and 16 reviewed (the
  native measures and, with the kernel, entity counts: its meetings across a
  stored face's diagonals, OCCT's seams). S9e.4b.2's `hollow_slab` solid
  with it, STEP-b's `box_void` imported. A trial of the 171 DRAW restore
  cases the import reaches: no case moves (59 evaluating on both backends, `bfuse_complex/K1` a tangency between the inputs and `bugs/modalg_6/bug28773` axes within rounding of parallel as before, none reaching an imported polyhedron against curved faces or a cavity; native DRAW's statuses as before).
* **Fuzzing.** The `boolean` target (`FUZZING.md`): the split target's line
  and arc profiles, the tool offset exactly in the axis-aligned frame or
  sharing the tilted one's origin, heights equal, spanning, overlapping,
  disjoint, inside (pockets, cavities) or on top (touching stacks); fuse,
  cut and common with the volume identities, rigid motions of the results
  and every result vertex classified on its boundary; one operation's first
  result (of at most 12 faces; curved ones too since S9e.1, `GIVEN_CURVED`)
  cut by a turned box and in common with it (S9e.2: or a turned cylinder,
  `GIVEN_ROUND`; S9e.3a: results of spheres, cones and tori too, or a
  sphere, `GIVEN_BALL`); S9e.4a: the object written by the kernel's
  `.brep` writer, read back and imported, given the chosen operation again,
  its volume the object's own result's (`IMPORTED`; since S9e.4b.1 its
  arcs' ends rounded off their circles taken onto them; since S9e.4b.2 the
  chained stage's first result of plane faces too, imported and cut by the
  turned box again, its volume the chained cut's; since S9e.4b.3a a first
  result of one sphere, cylinder or cone face and plane faces too, a plane
  piece; since S9e.4b.4b.1 one of several such faces too, a Boolean chain
  of their primitives, since S9e.4b.4b.2a one led by a prism leaf too,
  since S9e.4b.4b.2b.1 one with a flat, a prism of one cap or a pocket);
  S9e.4b.4a: a lens in the object's frame given the chosen
  operation, then written, read back, imported and given it again, both
  volumes equal (`JOINTS`, by the chained byte's next bit; in the tilted
  frame its joints round off both circles once read back); S9e.4b.4c.1: the
  chained cut's first solid of plane faces written, read back and imported
  and the kernel's own cut each less a ball about the turned box's axis, one
  volume (`MESHES`, by the chained byte's next bit).

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
in two, a turned box inside another) and 987 cases of upstream's
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
its poles, and 44 of S9d.3b.1's: the cylinder and a coaxial frustum on
its cap, inside it or through its caps, meeting it in circles, or a
frustum across it, its axis crossing the cylinder's at right angles,
meeting the wall in quartics, its wide end through the caps in `ZK7` to
`ZK9` and `ZL1`, and 20 more of S9d.3b.1's, two `pcone`s: a frustum and a
narrower coaxial one on its top disc, inside it or through its discs, and 16 of
S9d.4b.2a's: the cylinder and a coaxial torus whose tube its wall cuts in
two circles, and 4 of S9e.2's: DRAW's rollex, a disc less a pocket across
its rim, a stack with a cylindrical wall, given to the next `bcut` with a
cylinder standing on the pocket's floor) evaluate on both backends, and
S9c.2b.1 adds none. One of that
sphere's turns, `ZI5`, was wrong until a sphere face's closing chord at
a pole was enclosed narrowly (its volumes off by `32 pi / 9`, its `btuc`
refused by the kernel's validation), and is registered since. A volume
audit of the registered cases found Rust's volumes and centres of
gravity native DRAW's to its printed digits, or where they differ (eight
cases) nearer the closed forms; run again in S9d.4b's survey, Rust's
values were unchanged bit for bit, the new cases' the closed forms'
within 2.6e-15 relative; in the survey of S9d.2c, S9d.3c and S9d.4c 129
differ from it within rounding (2.2e-15 relative at most: the certified
integrals' speed-up 118, S9d.4c 11 centres), the closed forms' and the
references' within 2.8e-15. Eight more of S9d.3b.1's, a wider
frustum across the cylinder (`ZK7`, `ZK8`), evaluated right but took 23
to 120 seconds on the debug worker, past the contract's 30, until the
certified integrals' speed-up (9 to 13 since), and are registered since
that survey, in which no upstream case's status changed with S9d.2c,
S9d.3c or S9d.4c (the group's spheres and tori are whole, no cap, zone
or torus part, and its cones meet the cylinder in circles or rings, not
loops). The survey of S9e.1 and S9e.2 found no other case newly
evaluating (the rollex registered in S9e.2's run of the chained cases),
the volume audit's values of the 983 before bit for bit and the rollex's
`generate_given_boolean_fixtures.py`'s within 3.2e-16 relative. The
survey of S9e.3a and S9f.1 (with `b2765f20`'s plane crossing a sphere
within the resolution of tangency) found no status changed and no case
newly evaluating (S9f.1's spline solids in the group come from
`nurbsconvert`, which the adapter does not run; S9e.3a's given results
reach only `G9` and `H3`), the volume audit's values of the 987 bit for
bit, and no registered sphere near a tangency. The survey of S9e.3b,
S9f.2a and S9f.2b.1 found every status and refusal unchanged and no case
newly evaluating (no spline wall of the group meets another prism's wall,
and no given result the adapter makes holds a meeting of two curved
faces met by the partner: `G9` and `H3` stay a tangency), the volume
audit's values of the 987 bit for bit, the tori of `ZL2` to `ZL5` and the
spheres through their poles of `ZI4` to `ZI7` among them, so S9e.3b's
validator for torus bands, bounds of curved faces and pole vertices move
none. Of
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
arc ending off its circle; of the results
given to another Boolean, a frustum fused onto a cylinder (`bcut_simple/
G9`, `H3`: a result of solids other than prisms, S9e.3, until S9e.3a; since
then the rod's circle touching the frustum's top circle, a tangency
between the inputs; `UPSTREAM_TESTS.md`).
