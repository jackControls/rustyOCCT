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

* `OutOfDomain`: a solid or plane of a later sub-step (S8a.2: planes oblique
  to a prism's axis; S8c: primitives).
* `Degenerate`: the split would leave an edge or a piece thinner than the
  resolution (a crossing within the resolution of a vertex, a plane within
  binary64 of a cap), or pinch a piece (a plane tangent to a hole, or
  touching the profile at a vertex, inside the solid): a profile cannot hold
  a hole touching its boundary.
* `ComputationLimit`: a certified comparison it could not decide (an arc's
  crossing at its end).

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
* **Oblique** planes (S8a.2) build general bodies with ellipse-arc edges.

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
  kernel two reviewed count differences: OCCT splits an arc wall along a
  tangent ruling, and keeps a seam in a cut cylinder.
* **Kernel.** `tests/split.rs`: every side's sums of the kernel's enclosures
  contain the reference's volume, area and moments; histories pass the
  independent check, cover every input entity and repeat exactly;
  `compare_split.py`: the 13 cases of S8a.1 inside the reference, the 13
  oblique ones pending S8a.2.
* **Fuzzing.** The `split` target cuts rectangles, regular polygons,
  stadiums, U shapes and holed squares in two frames with planes chosen
  degenerate on purpose; volumes add up, pieces lie on their sides, and the
  split's own history check runs in its debug build.
