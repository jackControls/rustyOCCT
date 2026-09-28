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
hyperbolas), `Procedural` (S7b), `Traced` (S7b.3b) or `NotConic` (a pair not
yet parameterised); a spline surface is out of domain. It is symmetric in its
arguments.

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
| Other quadric pairs | procedural curves (S7b.1, S7b.2); two cones and a cone's rational apex on a sphere not yet |
| Tori | a plane or a sphere: procedural curves, or circles in the special cases; coaxial pairs: circles (S7b.3a); a cylinder, a cone or another torus off the axis: traced curves (S7b.3b) |

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

## Curves with a cone (S7b.2)

A sphere and a cone (the centre off the axis) and a cylinder and a cone
(axes not coaxial) are procedural curves too; two cones and a cone whose
rational apex lies on a sphere stay `NotConic` (a later part of S7b).

* **Sphere and cone.** Parameterised on the cone's rulings through its apex
  `V`: `V + v d(u)`, `d(u) = cos h a + sin h (cos u x + sin u y)`, `x`
  towards the centre. Then `A = 1`, `C = |V - c|^2 - R^2` is constant and
  `B(u) = b0 + b1 cos u` a sinusoid: the apex inside the sphere gives two
  rings (one on each nappe); outside, loops where `B > sqrt(C)` or
  `B < -sqrt(C)`, a `cos u` threshold each, so their ends are `arccos` in
  closed form (a loop around `u = 0` is `[-t, t]`, around `u = pi`
  `[t, 2 pi - t]`), or two rings when `B` stays beyond `sqrt(C)`.
* **Cylinder and cone.** Parameterised on the cylinder, `A` constant.
  `D(u)` has no convenient closed form: its roots are isolated by certified
  subdivision (binary64 intervals, a piece they cannot settle again in
  rational intervals) with the mean-value enclosure `D(m) + D'(piece)
  (piece - m)` and a work budget, each certified root then narrowed by
  bisection on certain signs to a few units in the last place. Loops lie
  between consecutive roots where `D > 0`, two rings where `D > 0`
  throughout. A double root (a tangency) cannot be certified and is
  `ComputationLimit`; it cannot occur exactly with a binary64 half-angle.

Evidence: 15 more fixture cases (`ck_`, `ks_`: crossing, parallel, skew and
tilted axes, a miss, rings, two loops, an irrational apex, the apex inside
the sphere, two loops on both nappes, a near-tangent sphere), the reference
extended with the cone's parameterisation, and a second `GeomInt_IntSS`
capture taken before any kernel cone code
(`fixtures/occt-procedural-cone-preimplementation`): 14 native matches and
one reviewed difference (samples of the near-tangent loops, ill-conditioned,
`1.55e-6` from the exact curve), and the kernel inside the reference on all
15. Native samples are compared by closest-point distance: near a loop's end
the branches are vertical in the ruled parameterisation, so a sample's own
angle is not its nearest parameter.

## Tori (S7b.3a)

A torus (major `R`, minor `r < R`) is parameterised by the angle `phi` of
its meridian in an exactly orthonormal frame whose `x` is the component
normal to the axis of the plane's stored normal, or of the direction from
the torus's origin to the sphere's centre:
`C(phi) + r (cos t e + sin t a)`, `C = o + R e(phi)`. On a meridian circle a
plane or a sphere is `f0 + alpha cos t + beta sin t`, so each `phi` gives two
points, `t = atan2(beta, alpha) +- arccos(-f0 / sqrt(alpha^2 + beta^2))`,
evaluated algebraically (`MATHEMATICS.md`); the branches are `+` and `-`.

* **Classes.** `D = alpha^2 + beta^2 - f0^2` is a quadratic `P(c)` in
  `c = m cos phi` with rational coefficients and a negative leading one,
  `m` the length of that normal component. The exact signs of `P(m)`
  (`phi = 0`) and `P(-m)` (`phi = pi`), numbers `u + v sqrt(q)`, and of the
  vertex against `+-m` give every class: empty, a tangent point at `0` or
  `pi`, one loop around `0` or `pi`, two mirror loops, two loops sharing an
  end at `0` or `pi` (a node), two rings, a figure-eight with its node at
  `0` or `pi`, and Villarceau-like pairs of loops `[0, pi]`, `[pi, 2 pi]`
  sharing both ends. A loop's ends are `arccos` of a root of `P` over `m`.
  `P`'s discriminant is a sum of squares, zero only for a plane containing
  the axis or a sphere centred in the equatorial plane with
  `rho^2 = r^2 + |w|^2 - R^2`: empty when the centre is nearer the axis than
  `R`, else the sphere contains a meridian circle (`NotConic`).
* **Special cases.** A plane containing the axis: two meridian circles. A
  plane normal to the axis: circles of radius `R +- sqrt(r^2 - z^2)` (one of
  radius `R` when tangent). Coaxial pairs meet where their meridians meet in
  a half-plane: a sphere or another torus (two circles, exact classes, one
  when tangent; the same torus is `Same`, a concentric one of another minor
  radius `Empty`), a cylinder (a vertical line), a cone (two lines through
  the apex, a quadratic per nappe with a certified discriminant, never
  exactly tangent). All are circles about the axis.
* **Points.** `ProceduralCurve::point_at` on a torus curve evaluates the
  meridian's point in binary64 intervals, rational ones when those cannot
  decide or are too wide; `carrier()` is the torus.

Evidence: 58 fixture cases (`tp_`, `ts_`, `tx_`: every class exactly, with
spheres on rational centres for the tangencies a binary64 plane normal
cannot make exactly, near cases beside them, tilted frames), the reference
extended with the meridian parameterisation (`D` by evaluating the other
surface on the meridian circle at three angles, its own exact predicates on
`P` found by exact interpolation, every simple loop end checked as a sign
change of that `D`), and a `GeomInt_IntSS` capture taken before any kernel
torus code (`fixtures/occt-procedural-torus-preimplementation`): 55 native
matches and 3 reviewed differences (tangent points `GeomInt_IntSS` misses),
the kernel inside the reference on all 58. The fixture's stored normals are
checked bit for bit (`procedural-intersection-frames.tsv`): a near-Villarceau
case is sensitive to one unit in the last place of the normal, and the
platform's `hypot` can round a non-Pythagorean normal differently.

## A torus with a cylinder, a cone or another torus (S7b.3b)

A cylinder or a cone off the torus's axis meets a meridian circle in up to
four points and no closed form separates them: the result is a
`TracedCurve` (`SurfaceIntersection::Traced`), the zero set of
`G(phi, t) = f(p(phi, t))` on the flat parameter torus (the meridian
parameterisation of S7b.3a; `x` along the other axis's component normal to
the torus's axis), as a graph. For two tori (S7b.3b.2) the carrier is the
first by stored data and `f` the other's quartic; along a meridian circle
`|p - o2|^2` is affine in `cos t`, `sin t`, so `G` is of degree two in them
for every pair.

* **Tangencies** (the curve's singular points) are decided exactly: a torus
  and a cylinder are pipes about the spine circle and the axis, and touch
  where a critical distance between them is `r + r_c` or `|r - r_c|`, found
  from exact resultants of three conics in the spine's plane and decided on
  the algebraic roots (`tangency.rs`). Each is a `Node`, crossing or isolated
  by the certified sign of the Hessian of `G`. Two tori touch where their
  spines have a critical pair at the distance `r1 + r2` or `|r1 - r2|`: the
  criticality along the first spine (a cubic once the distance condition is
  used) and the distance (a quartic) are reduced modulo the first spine's
  circle, and the critical points are the real roots of the resulting
  resultant, decided on the algebraic roots. A cone is never exactly tangent
  to a torus for a binary64 half-angle.
* **Folds**, where a component turns in `phi` (`G = G_t = 0`), are found by
  subdivision of the parameter torus with mean-value exclusion and
  certified by the Krawczyk operator; each is a `Fold` with enclosed angles
  and point.
* **Boxes and tracks.** Every fold and tangency gets a box with no zero of
  `G` on its top and bottom and the certified simple roots of its local
  picture on its sides. Outside them each branch is a `Track`, a chain of
  certified windows (`G` changes sign across the window for every `phi` of
  the step, `G_t` keeps a sign), followed once round the torus from a
  meridian outside every box; the windows fix how branches continue (through
  `t = pi` too) and bracket every point evaluation.
* **Components.** `Smooth` (a closed curve, its folds and winding numbers
  `[w_phi, w_t]`), `Crossing` (branches through crossings) and `Isolated`
  (a tangency point). `TracedCurve::point_at(track, phi)` and `t_at` give
  enclosed points: binary64 intervals, then rational ones narrowed by the
  interval Newton operator from the binary64 bracket when those are wider
  than `1e-12` relative. A budget exhausted anywhere is `ComputationLimit`;
  a higher contact, or an axis meeting the spine with `r = r_c`, too.

Evidence: `torus_curve_reference.py` (80-digit polynomial roots of every
meridian in `z = e^{it}`, the critical meridians as the real roots of a
resultant of `z^2 G` and `z^2 G_t`, components by continuity, tangencies by a
Groebner basis in sympy: of three conics for a cylinder, of both spines'
points for two tori) and 38 fixture cases (24 with a cylinder or a cone, 14
pairs of tori: side by side, linked, a ring round the tube, through the
hole, over it, tilted, touching at outer equators, at a saddle from inside,
at two saddles crossing, at right angles, near cases); for the first 24 (`torus-curve-*`: loops through
the tube, round the hole, rings, curves winding round the tube, tilted
frames, isolated tangencies from outside, inside the hole and over the tube,
a crossing at a saddle, near cases of each); a `GeomInt_IntSS` capture before
any kernel code (`fixtures/occt-torus-curve-preimplementation`) and for the
14 pairs another before the kernel intersected two tori
(`fixtures/occt-torus-pair-preimplementation`): 28 native matches and 10
reviewed differences (six missed isolated tangencies, a missed tiny loop,
three approximated lines off the surfaces), and the kernel inside the
reference on all 38 (`compare_torus_curves.py`, `torus_curves.rs`: folds,
tangencies and rings' points enclosed, components and winding numbers equal,
order independence, points along every track on both surfaces).

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
