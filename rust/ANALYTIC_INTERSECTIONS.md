# Intersections of analytic surfaces

S7 of `REVIEW_NOTES.md` intersects planes, cylinders, cones, spheres and tori:
the mathematical core the split (S8) and the Booleans (S9) stand on. S7a,
described here, covers every pair whose intersection is empty, the same
surface, points, lines or a conic; S7b adds the procedural curves of the other
quadric pairs and of tori (D13 of `TOPOLOGY_MODEL.md`), S7c curve/surface and
S7d curve/curve intersections.

## Contract

`intersection::surface_surface(a, b)` takes two topology `Surface` values and
returns `Empty`, `Same`, `Items` (points, lines, circles, ellipses and
hyperbolas) or `NotConic` (S7b); a spline surface is out of domain. It is
symmetric in its arguments.

* **Exact surfaces.** A surface is the exact point set its stored binary64
  data define: a plane through the stored origin with the stored normal; a
  cylinder, cone or sphere about the line through the stored origin along
  the stored normal, normalised exactly (`n n^T / |n|^2` is rational), with
  the stored radius and half-angle. A cone is both nappes.
* **Exact degeneracy.** Parallel, perpendicular, coaxial, coplanar,
  tangent, coincident and through-an-apex configurations are decided by
  exact rational predicates on the stored data. A near-degenerate input gets
  the generic result: two planes `2^-30` rad apart meet in a line `2^30`
  away, a cylinder nearly parallel to a plane gives an ellipse `2^30` long.
  Merging within a tolerance belongs to the face algorithms (S8, S9), which
  have bounded faces and a resolution.
* **Cones.** `cos a` and `sin a` of a nonzero binary64 angle are
  transcendental, so a plane is never exactly parallel to a generatrix (no
  parabola), and a cone's apex is rational only when its radius at the
  origin is zero. Comparisons involving them are certified in rational
  intervals; one that cannot be separated is
  `Error::ComputationLimit`.
* **Enclosures.** Every number returned is a canonical parameter with
  binary64 bounds `[lo, hi]` containing its exact value, typically adjacent
  binary64 values.
* **Canonical forms.** A line by its point nearest the global origin and a
  unit direction; a circle by its centre, unit normal and radius; an ellipse
  by its centre, unit normal, unit major direction and semi-axes (major
  first); a hyperbola (both branches) by its centre, unit normal, unit
  transverse direction, semi-transverse and semi-conjugate axes. A unit
  direction's first nonzero coordinate is positive; items are sorted by
  kind, then by their numbers.

## Pairs

| Pair | Results and exact degeneracies |
| --- | --- |
| Plane/plane | a line; parallel: empty or the same plane |
| Plane/sphere | a circle; tangent: a point; else empty |
| Plane/cylinder | an ellipse; axis normal to the plane: a circle; axis parallel to it: two lines, one (tangent) or none |
| Plane/cone | an ellipse or a hyperbola; axis normal: a circle, or the apex; through the apex (a rational apex, or a plane containing the axis): two lines or the apex |
| Sphere/sphere | a circle; tangent: a point; concentric: empty or the same; else empty |
| Cylinder/cylinder | parallel axes: two lines, one or none; coaxial: the same or empty; equal radii with crossing axes: two ellipses in the bisecting planes; otherwise not a conic |
| Coaxial pairs | a cylinder, cone or sphere centred on the other's axis: circles where the radius functions of the axial coordinate agree, the apex where both pass through it, or the same surface |
| Other quadric pairs, tori | not a conic (S7b) |

The plane/cone section uses a closed form: with the apex `V`, the plane's
unit normal `n`, `cos b = n . a` (the unit axis), `D = (V - o) . n` and
`e1` the unit projection of the axis on the plane,
`A (s - s0)^2 + cos^2(a) t^2 - D^2 cos^2(a) sin^2(a) / A = 0` in the plane
coordinates along `e1` and `n x e1` about the foot of `V`, with
`A = cos^2(b) - sin^2(a)` and `s0 = -D cos(b) sin(b) / A` (`MATHEMATICS.md`).

## Procedural curves (S7b.1)

Two cylinders with crossing axes (other than S7a's equal radii through a
common point) and a cylinder and a sphere off its axis meet in a curve that is
not a conic: `SurfaceIntersection::Procedural(ProceduralCurve)`, D13's
procedural curve. It keeps both surfaces and is parameterised on the ruled
one (the thinner cylinder; of equal ones the first by stored normal, then
origin; the cylinder of a cylinder/sphere pair) by the angle `u` of its
ruling in an exactly orthonormal frame whose `x` points along the common
normal towards the other axis (towards the sphere's centre). A ruling meets
the other quadric where `A v^2 + 2 B(u) v + C(u) = 0`, so each `u` gives the
points of two branches, `v = (-B +- sqrt(D)) / A`.

* **Classes.** Exact rational predicates on the stored data: for two
  cylinders the squared axes' distance `d^2` against `(r1 + r2)^2` (empty,
  a tangent point) and `(r1 - r2)^2` (a figure-eight when equal, two rings
  around the thinner cylinder below); for a sphere its radius against
  `|e - r|` and `e + r`, `e` the centre's distance from the axis (empty, a
  tangent point, a loop, Viviani's figure-eight, two rings).
* **Components.** A `Loop` over `[-t, t]` (both branches, joined at its ends),
  two `Ring`s (one branch each over a whole turn) or a `FigureEight` (both
  branches touching at the node, `u = pi`). In the canonical frame the
  discriminant has a closed form (`MATHEMATICS.md`), so `t = arccos c` is one
  certified arctangent.
* **Points.** `ProceduralCurve::point_at(u, branch)` encloses the points for
  every parameter in an enclosure of `u`: binary64 intervals, rational ones
  when those cannot decide or are wider than `1e-12` relative (near a loop's
  end, where the square root of a small discriminant widens them). An error
  where the ruling certainly misses the other surface. The pcurve on the
  ruled surface is `(u, v)` exactly.

Evidence: `procedural_intersection_reference.py` (exact classes, 80-digit
roots by dense sampling and `findroot`, points) and 20 fixture cases with
every class exact and near (`procedural-intersection-*.txt|tsv`);
`procedural_intersections.rs` requires the class, loop ranges containing the
reference's, the curve's points at the loop's ends, the rings' and the
figure-eight's points containing the reference's, the loop's middle within
`1e-12`, order independence and sampled points on both surfaces;
`GeomInt_IntSS` was captured before any kernel code
(`fixtures/occt-procedural-intersection-preimplementation`) and
`compare_procedural_intersections.py` finds every native sample within
`1e-6` of the exact curve and every component covered in 17 cases, with 3
reviewed differences (two missed tangent points, one sample `2.3e-5` off near
Viviani's node), and the kernel inside the reference on all 20.

## Evidence

* **Independent reference.** `analytic_intersection_reference.py` computes
  every case from the same stored data by another method: it substitutes a
  rational basis of the plane into the quadric's implicit equation (exact
  rationals for spheres and cylinders, 80-digit mpmath for cones) and
  classifies the quadratic, and it decides every degeneracy by exact
  predicates. `generate_analytic_intersection_fixtures.py` writes 67 cases
  (`analytic-intersection-cases.txt`, `-expected.tsv`) with an exact and a
  near case of every degeneracy class, and each surface's stored normal
  (`-frames.tsv`): inputs are chosen so that `Frame3::new` stores the same
  bits as the reference's normalisation (a normal `(0, 3, 4)`, not
  `(0, 0.6, 0.8)`, whose normalisation depends on the platform's `hypot`).
* **Kernel tests.** `analytic_intersections.rs` requires the stored normals
  to match, every reference number to lie inside the kernel's enclosures
  (up to `1e-25` relative, the reference's printing), results symmetric in
  their arguments, enclosures within four units in the last place, every
  returned curve on both surfaces and translations to move the items.
* **Native.** `occt_analytic_intersection_oracle.cpp` runs
  `IntAna_QuadQuadGeo` with `Precision::Angular()` and
  `Precision::Confusion()`; its observations were captured before any kernel
  intersection code (`fixtures/occt-analytic-intersection-preimplementation`).
  `compare_analytic_intersections.py` reproduces the capture, finds 63
  cases where OCCT's canonical results equal the reference's within `1e-9`
  and four reviewed differences, each a near-degenerate configuration IntAna
  snaps with its tolerances (`occt-analytic-intersection-divergences.json`),
  and requires the kernel inside the reference on all 67 and within `1e-9`
  of OCCT wherever OCCT agrees (43 cases with items).
* **Fuzzing.** The `analytic_intersections` target builds pairs with exact
  degeneracies on purpose (a shared axis, a shared origin, a tangent offset)
  and checks symmetry, ordered enclosures, every returned curve on both
  surfaces and exact translations.

## Acceptance

S7a's acceptance needs kernel CI green at the accepted revision, the
schedule run's full fuzz replay green there, and a clean local 600-second
`analytic_intersections` campaign. Recorded here once met.
