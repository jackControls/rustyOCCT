# Booleans of solids

S9 of `REVIEW_NOTES.md` fuses, cuts and intersects solids (the Combine
job): the faces intersected (S7), split (S8), classified by regions and
assembled into shells and regions, with complete histories. This document
describes what is implemented; the decisions are in `REVIEW_NOTES.md` (S9).
S9a is implemented (`profile/boolean.rs`, `solid/boolean.rs`, S9a.2's
stacks in `solid/boolean/stack.rs`); S9b on are not.

## Contract

`Solid::fuse(operation, other)`, `Solid::cut(operation, tool)` and
`Solid::common(operation, other)` return the result's solids (each a
maximal connected solid region, ordered by the lowest input face id it
keeps, none for an empty result) and the operation's history (`Fuse`,
`Cut`, `Common`). A result that is itself a prism of one of the inputs'
frames (one profile over one height range) is built as one, keeping the
prism's exact queries; any other is a general body built through
`TopologyParts` and validated before it is returned. An error is one of:

* `OutOfDomain`: a pair of a later sub-step (frames whose axes differ, or
  an offset or a profile that rounds, before S9b; spline profiles), or a
  cavity in a result of several solids.
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
  count, all 45 within 7.9e-16, no review.
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
unsupported; the derived case `boolean_prisms` and 320 cases of upstream's
`boolean` group evaluate on both backends.
