# Sustained geometry fuzzing

`rust/fuzz` is a separate cargo-fuzz workspace. Its fuzz-specific dependencies and C++
libFuzzer runtime are test-only. The production kernel remains native Rust.
Both workspaces have checked-in dependency locks. The campaigns run optimized
code with AddressSanitizer, debug assertions, overflow checks, and libFuzzer
coverage feedback. They are distinct from the deterministic invariant tests.

| Target | Generated inputs | Assertions |
| --- | --- | --- |
| `predicates` | Full binary64 bit patterns, scaled integers, coplanarity, sphere-boundary points, non-finite coordinates | Exact 2D/3D orientation and insphere against independent rational matrix elimination; sphere ordering; typed rejection |
| `intersections` | Lines/segments, three-point planes/triangles, full binary64 exponents, coplanarity and degeneracies | All three intersection APIs against a rational barycentric oracle; exact classification; minimal finite coordinate/parameter enclosures; explicit unrepresentable results |
| `linear_sets` | All 25 linear primitive pairings, raw binary64, coplanar/shared-vertex modes, full exponent range, invalid definitions and collapsed segments | Complete exact sets from independent boundary crossings and a gift-wrapping hull; canonical operand/winding symmetry; minimal construction bounds; zero-distance consistency on scaled-integer modes |
| `proximity` | All 25 point/line/segment/plane/triangle pairings; arbitrary binary64 coordinates, scaled/coplanar/shared-vertex cases, collapsed segments and invalid inputs | Independent exact membership and supporting-plane certificates prove each returned pair globally minimal; minimal output enclosures, threshold comparisons, operand symmetry, coordinate permutation/vertex reversal and deterministic witnesses |
| `spline_proximity` | Positive-weight rational space curves through degree 25 with known zero and nonzero minima, nonplanar cylinder curves, rational/irrational factor roots, composed exact knot/degree edits, trimmed/singleton ranges, near-tie parabolas, rational circles and periodic polylines | Complete minimum parameter sets and whole intervals from independent factor and radical signs, a monotone cubic, circle identities and exhaustive rational segment projections; tight parameter/coordinate/distance bounds; exact ties and periodic aliases |
| `spline_linear` | Known-factor rational Bézier curves through degree 25 with sub-float root pairs, rational polylines with periodic turns, constant spans and retracing, a rational circle against rational slopes, a quadratic retrace, lines/collapsed/reversed segments, invertible integer affine maps, parameter domains through 2^±2048 and exact knot/degree edits | Complete point/interval sets from constructed factors, per-span rational linear inequalities and closed-form quadratic sign functions; exact parameter, coordinate and line-parameter identity within 2^-96 independent brackets; tight or typed unrepresentable views; reversal symmetry, malformed-rational and work-limit rejection |
| `brep_validation` | Valid cell complexes: star-outline prisms with no hole, a seamless round hole, a square hole or an inverted box cavity whose shells and regions are merged by fuzz-crate code, at scales through 2^±10 with random frames, plus the stadium fixtures; one of 33 mutations: extra vertex/edge, empty shell/loop, bad references, faces leaving or repeated in shells, flipped fins and faces, moved vertices, shifted pcurves, inverted shells, swapped loops, and the cell model's failure modes (a fin shifted by whole periods, a flipped winding, fin lists exchanged between edges, a side at the wrong shell, a vertex loop off its surface, an open face's only loop removed), and enclosures (M5): a bound removed, set outside `[0, resolution]` (twice it, negative, NaN, infinite), a vertex moved half the resolution, circle prisms up to `2^48` from the origin at a tolerance of a few ulps there, and cones or frusta (S3) with the pole moved along a ruling or the axis, the pole removed, a right semi-angle or a ring pcurve shifted, and whole spheres, hemispheres and zones with the pole moved over the sphere or along the axis, removed, a zero radius, a ring pcurve shifted or a whole sphere turned inside out, and whole tori, v-segments and wedges with a loop winding twice, a pcurve shifted off the tube, the tube reaching the axis or a whole torus turned inside out, and (R4) a prism's line edge, line pcurve or plane made a spline of degree 2 or 3 with a knot repeated to the degree, exactly C1, then broken there, and (S4) the valid spline fixtures (spline prisms, and stadiums with spline geometry on their cylinder) moved by exact similarities (the rational corner included, and since F8 the knotted walls and the stadium's spline parallel), then a spline use's pcurve shifted or a ruled wall reversed | Clean report for every base; exact issue lists for local mutations and a required issue on the mutated entity otherwise; a `10·tol` pcurve shift clean at `100·tol`; a vertex loop on a cap and a reversed two-fin radial order stay valid; exactly `enclosure_missing`, `enclosure_exceeds_resolution` or `enclosure_unsound` for the enclosure mutations; a far prism fails to construct or encloses every entity within its resolution; a cone constructs unless its radii are equal, encloses within its resolution, has a certified volume within `1e-9` of the frustum formula and its synthesized counts, and reports exactly `pole_off_apex`, or `vertex_loop_off_surface`, `loop_without_face` with `winding_mismatch`, `degenerate_surface` or `pcurve_off_edge`; an exact repeat reports no continuity issue and the broken knot exactly one more (`edge_not_c1`, `pcurve_not_c1` or `face_not_c1`); a moved spline prism is valid with a certified volume, area (F8) and centroid enclosing the reference's, moved, a shift far beyond the tolerance is `pcurve_off_edge`, one within it at most an unsound enclosure, and a reversed wall `loop_winding`; deterministic, duplicate-free reports identical to `from_parts` |
| `identity` | Structure-aware (`arbitrary::Unstructured`) profiles: 3–12-sided or circular outlines, clockwise or counter-clockwise input, up to three square or round holes, labels absent or present, random operation ids, frames, directions, scales through 2^±8 and up to two rigid motions; the same bytes also make a cone or frustum with apices at either end, a whole sphere, hemisphere or zone, and a whole torus, v-segment or wedge; and (S5) a filleted rectangle with an optional concave notch at dyadic sizes, entered either way; (S6) the prism's profile as a face body and its outer boundary as a wire body | An independent version-1 encoder and FNV-1a-128 recompute every id from its derivation; ids unique and slot maps inverse; rebuilding, rigid motion, reversed direction and moved labelled points keep ids; permuting label values permutes parents bijectively; unlabelled builds share no id; entity counts and roles follow the profile (a circle sweeps two ring edges and a wall, no vertices), one region generated from every boundary, and the synthesized OCCT counts add a seam and two seam vertices per circle; a cone's entities, parents (its meridian elements) and counts follow its radii, and rebuilding, stretching and rigid motion keep its ids; likewise a sphere's, whose pole is a vertex only when it closes a band, and a torus's, whose wedge sides come from the boundary and the arc; an arc path validates with a polygon's counts, clockwise input and reversal keep its ids, its closed-form area and volume lie in the certified enclosure, and its filleted corner classifies outside; the face body validates as a sheet with a free face's counts (a ring edge closed at one vertex), its face, edges and vertices derive (kind `MakeFace`) from exactly the prism's start cap, bottom edges and bottom vertices' parents, rigid motion keeps its ids and its certified area contains the profile's; the wire body validates as a wire (kind `MakeWire`) with one wire when it has several edges |
| `history` | The `identity` structure (profiles, circles included, holes, labels, frames, scales and up to two rigid motions) extruded in either direction, plus one of seven mutations of the produced history, or a synthetic split of a real edge; and a cone, a sphere and a torus from the same bytes | Every construction and transform history checks clean, resolves every input to itself, composes associatively and replays identically at its recorded algorithm level; the construction generates exactly one region; a dropped relation reports its missing target and source, a duplicate its order and source, a changed kind a dangling id or differing geometry, a dangling target both the phantom and the lost target, swapped split children both ordinals, and a construction after its own output an invalid composition; a cone's construction checks clean with kind `Revolve`, replays at its level and at no other, generates every entity once, composes with rigid motions and reports a dropped relation's target |
| `split_merge` | The `identity` structure (profiles, circles, holes, labels, frames, scales and rigid motions), a split height anywhere in the prism, one of six history mutations, and a separate construction stacked on the prism's end | A split fails exactly when a piece would be thinner than the tolerance; split, fuse and their composition check clean; the composition resolves every input to `Same`; every input id resolves and no new id does; volumes add up within the rounding of area x height; fusing the pieces back gives the original properties; a fuse with itself or with a piece twice is refused; the stacked fuse is independent of the call order; swapped split children report both ordinals, a dropped deletion its missing source, a deleted entity called unchanged a dangling id, a split given another parent's children a differing support or dimension, a merge naming one parent twice the lost and duplicated sources, and a fuse before its split an invalid composition; no continued entity's enclosure falls below its parents' through a transform, split or fuse, and every enclosure lies within the resolution |
| `attributes` | The `identity` structure, six keys with random policies, random keys and small values on every entity, a transform, a split at a random height, changes on one piece and the fuse of the pieces | `attributes::check` clean on every step; a merge conflicts exactly when its parents' values differ or one lacks the key; no attribute on a deleted id; re-running a step gives the same bodies and history; a key with no policy is refused |
| `brep_io` | A prism of the `identity` structure as the kernel writes it, or one of three upstream `data/occ` files, a cone and a sphere written with OCCT's seam and degenerated pole edges, a torus with its seams in `u` and `v`, (S4) six valid spline fixtures (spline prisms, a rational corner, stadiums with spline geometry on their cylinder, a spline face) written with B-spline records, and (S6) a valid sheet, shell, wire or acorn fixture or a profile's face or wire body written as a free shape, then up to six structure-aware token and line mutations: numbers moved, negated, truncated or replaced by record words, flags, references and special values; tokens deleted; lines duplicated, deleted, swapped or truncated; shape orientations flipped or made internal or external | Reading never panics and malformed text is a typed `BrepError` (every table and shape reference is range-checked); every imported solid and free shape validates; a body the writer expresses reads back to one solid or free shape of the same class with the same synthesized counts, cell counts and (without ring edges) bit-identical vertices; an unmutated prism, cone, sphere, torus, spline fixture or free shape always writes and round-trips |
| `analytic_intersections` | Two planes, cylinders, cones, spheres or (a kind byte of 224 or more, S7b.3a) tori on dyadic origins, small integer axes, dyadic radii and half-angles, independent or with an exact degeneracy made on purpose: a shared axis (parallel; a plane normal to a torus's axis), a shared origin and axis (coaxial, concentric), a tangent offset, a shared origin (S7a) | No panic, only a certified comparison it cannot decide (`ComputationLimit`); symmetric in its arguments (the same items, overlapping enclosures); ordered enclosures; every returned point, line, circle, ellipse and hyperbola on both surfaces within `1e-9` of its scale; an exact dyadic translation keeps the items; a procedural curve (S7b.1: two cylinders, a cylinder and a sphere; S7b.2: a cylinder and a cone, a sphere and a cone; S7b.3a: a torus and a plane or a sphere) is the same whatever the argument order and points along every loop, ring and figure-eight, and every finite point of every track of a traced curve (S7b.3b: a torus and a cylinder, a cone or another torus; S7b.4: two cones, a cone's apex on a sphere or a cylinder), lie on both surfaces |
| `tessellation` | Structure-aware bodies (T-a of `REVIEW_NOTES.md`): the `identity` structure's polygon and circle prisms with square and round holes and its cones, spheres and zones and tori, v-segments and wedges, each possibly moved rigidly; filleted and notched arc paths as prisms and as face bodies; a prism's profile as a face body; (T-b) prisms of a rectangle whose right side and maybe top are clamped quadratic or cubic B-splines of up to two spans, rational or not (weights 0.7 to 1.5), on axis-aligned frames at dyadic offsets, and face bodies on whole spline surfaces of degrees 2 and 3 bounded by their boundary rows, axis-aligned or tilted, built as validated parts (an invalid one skipped); a deflection from an eighth to a 2048th of the body's size (a 64th for splines) and an angle from 0.2 rad (0.35 for splines) to `π/2` | Every buildable body tessellates (no error); distinct nodes per triangle; every mesh edge in two triangles traversed once each way (a face body's boundary edges once, at least three per loop); the Euler characteristic of the boundary (`2 - 2 holes`, a whole torus 0, a face `1 - holes`); every reported bound within the request; 12 barycentric samples of every triangle within its bound of its face's whole surface by closed-form distances independent of the kernel's evaluation (on a spline surface the centroids of up to 64 triangles per face, by a projection with the target's own binary64 Cox-de Boor evaluation); a solid's enclosed volume within the deflection times the areas of its mass properties (a spline prism's by Gauss-Legendre quadrature of its profile); the same mesh again, bit for bit (a spline body's meshes of at most 500 triangles) |
| `curve_surface` | A plane, cylinder, cone, sphere or torus (the bytes of `analytic_intersections`) and a line or a circle on dyadic data, independent or with an exact degeneracy made on purpose: a line along the axis at the surface's radius (a cylinder's ruling), across the axis at that offset (tangent to a cylinder or a sphere) or through the origin; a circle coaxial at the radius, beside the surface at the sum of the radii, or in a plane through the axis (S7c.1); a shape byte of 128 or more: an ellipse or a hyperbola in the same modes, or a rational B-spline against a torus or a cone, with independent poles or the exact quarter circle of the surface's reference circle (S7c.2) | No panic, only a certified comparison it cannot decide (`ComputationLimit`); points sorted with ordered enclosures, each the curve's point at its parameter and on the surface within `1e-9` of its scale; a contained curve and every overlap on the surface; an arc gives its whole circle's result; an exact dyadic translation that keeps the stored frames keeps the classes, contacts and overlapping parameters |
| `curve_curve` | Two lines, circles, ellipses or hyperbolas on dyadic data, independent or with an exact degeneracy made on purpose: the second in the first's plane, the same curve in another representation (a line through two other of its points, a circle with its normal reversed, a conic with its axis reversed: an ellipse's same set, a hyperbola's other branch), two circles touching in their plane, a line along the first's axis (S7d.1); a first kind byte of 128 or more: a rational B-spline against a circle, an ellipse or a hyperbola, with independent poles, poles in the conic's plane, or the exact quarter circle on a circle (S7d.2) | No panic, only a certified comparison it cannot decide (`ComputationLimit`); swapping the curves swaps the parameters (overlapping enclosures) and keeps the contacts; points sorted, each both curves' point at its parameters within `1e-9` of its scale; a coincident pair's points on both curves; a spline's points on both curves at their parameters and its overlaps on the conic |
| `split` | A rectangle, a regular polygon, a stadium, a U or a square with a round or square hole on dyadic sizes, in an axis-aligned or a tilted frame, split by a plane chosen degenerate on purpose: parallel to the axis at a dyadic offset, through a profile vertex, along a profile edge, tangent to an arc or the hole, normal to the axis at a dyadic height or in a cap, or oblique (S8a.2: at a dyadic point, through a cap's vertex, or touching a cap's circle or arc) | No panic; only a sub-resolution or pinching split `Degenerate` (an oblique plane parallel to the tilted frame's axis to binary64 among them) or an undecided comparison `ComputationLimit`; every piece validates as it is built and the split's history check passes (debug build); the pieces' volumes add up to the solid's within `1e-9`; each piece's centre on its side; each piece moves rigidly with its ids and volume; a missing or touching plane returns the solid. A first byte of 224 on selects S8c and S8d: a cone, apex cone, whole sphere, zone, cap or whole torus in either frame, split normal to its axis (at a dyadic height or an end), through a dyadic point at an angle, containing the axis, oblique, through a cone's end-circle point containing its tangent (touching the rim), parallel to a ruling (a parabola) or through a frustum's virtual apex (its rulings; S8d.2); only a torus cut in two caps or a plane through an apex or pole off the axis may be `OutOfDomain` (S8d.3 splits tori by any other plane). A first byte of 192 to 223 selects S8b.3's spline prisms (a quadratic bulge, a three-span wave, a lens hole of two cubics given either way round) under the same plane modes, the tangent ones touching a spline at its apex, with the same oracle. The target's per-input limit is 60 s (a cap cut near its pole takes 8.8 s under ASan) |
| `step` | A first byte picks the input (the STEP import track of `REVIEW_NOTES.md`): the rest as a whole Part 21 file, the rest as the data section of a file with a valid header and a length unit, or one of the 22 authored fixtures of `fixtures/step` unmutated or with up to six mutations of its instances (numbers moved, flipped, scaled or zeroed; `.T.` and `.F.` flipped; references redirected, dangling included; entity names and tokens replaced by other entities, loop and shell kinds, flags, special numbers and punctuation; instances deleted, duplicated or swapped; the text truncated) | No panic; malformed text is a `StepError`; importing twice gives the same result; every imported body passes `Topology::check` at its resolution and, when the `.brep` writer expresses it, reads back with the same synthesized counts; every rejected body names its constructs or lists its issues; an unmutated fixture imports every body |
| `boolean` | Two prisms in one frame (S9a): the split target's line and arc profiles on dyadic sizes, the tool offset by a dyadic vector in the axis-aligned frame or sharing the tilted frame's origin, its heights equal, spanning, overlapping, disjoint, inside the object's (pockets, cavities; S9a.2) or on its top (touching stacks); fuse, cut and common | No panic; only a tool whose offset or profile rounds, or a cavity among several solids, `OutOfDomain`, a sub-resolution or self-touching result `Degenerate`, or an undecided comparison `ComputationLimit`; every result, S9a.2's stacks among them, validates as it is built and its history passes the independent check (debug build); when all three succeed, `V(A ∪ B) = V(A) + V(B) - V(A ∩ B)` and `V(A - B) = V(A) - V(A ∩ B)` within `1e-9`; every result moves rigidly with its ids and each of its vertices classifies on its boundary |
| `modeling` | Valid radial polygons/circles, optional holes, 49 scales, up to eight operations, plus raw invalid input | Repeated rigid transforms, reversed winding/offsets and planar splits; mass/first-moment conservation, topology, classification, bounds and finite positive properties |
| `curved` | Full binary64 coefficients, centers/radii/axes and line endpoints; scaled integers; exact/neighboring tangencies; generator and point segments | Quadratic root count, multiplicity, exact comparisons, minimal enclosures; circle/sphere/cylinder hits against independent polynomial-sign and axial-projection oracles; endpoint clipping and typed failures |
| `splines` | Raw binary64 poles/weights/knots, scaled geometry, degree 1..25, repeated knots, periodic seams, full-range wrapped parameters, explicit sides and derivative requests | Exact basis-function derivatives and closed quotient formulas independently check homogeneous pole interpolation; minimal position/derivative bounds; discontinuity, domain, nonfinite and overflow errors |
| `bezier_editing` | Raw binary64 and scaled rational controls, degree 1..25, clamped/unclamped/periodic splines, subnormal spans, multi-period extraction and composed edits with rational cuts | Complete homogeneous polynomial identities from independent Cox basis coefficients and affine substitution; exact jets/minimal bounds, positive weights, shared endpoints, reversal, elevation/split commutation and preflight limits |
| `knot_editing` | Degrees 1..25, raw binary64 and scaled controls, rational cuts, unclamped inactive controls, repeated knots, periodic seams/origin changes, refinement/removal sequences and invalid data | Entire raw-support homogeneous identities; independent coefficient-equation removal feasibility, exact complete controls, rational jets/extraction, minimal enclosures, batch order/duplicate behavior, round trips and unchanged rejection limits |
| `surface_knots` | Tensor degrees 1..25 in both axes, independent periodicity and raw unclamped support, arbitrary binary64 and scaled controls, rational insertions, batch/transpose identities, inverse edits and periodic origin deletion | Every transverse homogeneous Cox coefficient; shared independent fraction-free removal equations; full extracted patch coefficients and quotient jets; retained isocurves, parameter wrapping, invalid rational rejection and atomicity |
| `degree_elevation` | Curve/tensor degrees 1..25; clamped, unclamped, inactive and periodic axes; raw binary64 and extreme rational atoms; invalid targets and full 4096-control boundaries | Complete independent Cox-equation control reconstruction; original domains, positive weights, staged edits, refinement composition, tensor order/transposition, exact jets and typed degree/count rejection |
| `exact_spline_intersections` | Rational controls/weights and knots, degree 1..25, rational trims, close roots through 2^-2048 spacing, huge/subnormal domains, periodic seams and edit sequences | Complete known-factor contacts, rational secants and exact circle overlaps; independent full polynomial edit identity; exact parameter/coordinate comparisons, contact orders through 50, minimal finite bounds or typed conversion failure, traversal limits and malformed-rational rejection |
| `surface_editing` | Tensor degrees 1..25 in both directions, raw binary64 and scaled controls, independent periodicity, unclamped knots, subnormal domains, full low-degree multi-period queries and high-degree selected spans, composed edits and rational cuts | Complete tensor polynomial identities; all exact partials through order two and minimal bounds; exact isocurves/shared boundaries, reversals, transposition, elevation/split commutation and Cartesian output limits |
| `surfaces` | Rational tensor grids, independently periodic U/V, repeated knots, high degree in either direction, full-range wrapped parameters and malformed data | Independent tensor basis plus closed bivariate quotient formulas; exact mixed-partial and quadrant continuity decisions; minimal enclosures and typed failures |
| `roots` | Products of rational/irrational/complex factors through degree 25, repeated roots, closed-domain clipping, power-of-two coefficient scaling and arbitrary binary64 query polynomials | Complete expected root list from known factors; algebraic signs reduce independently in Q(sqrt(d)); multiplicities, cross-equation root ordering and equality against independently constructed factors, minimal enclosures and nonfinite rejection |
| `spline_intersections` | Rational Bézier curves with known factored plane numerators and squared sphere/cylinder contact equations; weighted rational lines against quadrics; nonperiodic/periodic rational polylines; varying weights, parameter/space scales, oblique planes, tangencies, knots and zero spans; explicit trim bounds, neighboring floats and large periodic offsets | Complete parameter/position bounds from rational Bernstein evaluation, affine span equations, or independent geometric quadratic roots with rational weight-parameter conversion; one-sided orders, crossing/tangent/boundary classification, maximal clipped overlaps, closed seam events and repeated turns |

Planned targets and mutations are specified with their milestones:
T1 added the six cell-complex mutations to `brep_validation`, M3
`split_merge`, M4 `attributes`, T2 `brep_io`, and M5 the enclosure mutations
of `brep_validation` and the enclosure law of `split_merge`; S3 added cones
spheres and tori to `brep_validation`, `identity`, `history` and `brep_io`. R4 added spline continuity to `brep_validation`; S4 added the valid spline fixtures to `brep_validation` (mutation 32) and `brep_io`; F8 put the rational corner back into mutation 32 (its certified quadrature takes the bulge's time; its first-order enclosures had taken seconds per input under AddressSanitizer) with the knotted walls, the stadium's spline parallel and an area check; S6 added face and wire bodies to `identity`, free shapes (seeds starting `F`) to `brep_io` and sheets and wires (seeds starting `W`) to `brep_validation`; S7a added `analytic_intersections`; T-a of the tessellation track added
`tessellation`, and T-b its spline prisms and sheets; S7c.1 added `curve_surface`, S7c.2 its conics and splines, S7d.1 `curve_curve`, S8a `split`, STEP-a of the STEP import track `step`. A planned target is not evidence until its
clean campaign is recorded there.

The rational oracle uses `num-rational` with Gaussian elimination and
barycentric coordinates, plus polynomial-sign/vertex comparisons and cylinder
axial projection. Production uses a fixed binary64 integer lattice, determinant
expansion, edge half-planes, exact radical comparisons and, for analytic lines,
cylinder cross products. Spline production uses differentiated de Boor pole interpolation;
its fuzz oracle instead evaluates basis functions and their derivative identity.
Spline/quadric line fuzzing independently solves the physical line's quadratic
and converts its parameter through the rational weight map; production instead
isolates the spline's homogeneous implicit polynomial. Separate Python quadric
fixtures use cylinder cross products to check production's dot-product form.
The proximity checker uses convex supporting-plane inequalities, independently
of the production face enumeration and exact normal-system solve. Python
fixtures additionally use closed analytic projection/cross-product formulas.
Complete linear-set fuzzing uses cross-product line/plane formulas and boundary
edge candidates with a gift-wrapping hull; production solves affine equalities
and enumerates feasible halfspace vertices before a monotone-chain hull.
Editing production uses spline blossoms (curves), exact boundary knot insertion
(tensor extraction) and homogeneous de Casteljau;
its checker uses basis polynomials, binomial affine substitution and Bernstein
coefficient expansion. Clearing common denominators before those linear
transforms reduces repeated GCD work without changing the assertions.
Knot editing uses exact Boehm insertion and inverse insertion in production;
its independent checker solves Cox power coefficient equations with integer
elimination, including inactive unclamped support and failed removals.
The Rust oracles and production share `num-bigint`/`num-rational`;
the checked-in Python `Fraction` fixtures provide a separate integer runtime.
Coverage counts include the oracle and dependencies: they are not kernel-only
coverage percentages or evidence of exhaustive input coverage.

## Continuing campaigns

[Rust geometry fuzzing](../.github/workflows/rust-fuzz.yml) runs all twenty-nine targets
for 60 seconds of mutation each on relevant pushes/PRs, and 600 seconds each every day at
06:23 UTC on the default branch. Manual runs accept 1–3,600 seconds per target.
GitHub can delay scheduled jobs. The schedule must remain enabled on the fork.

Every run restores the previous evolving corpus, adds seeds derived from the
exact fixtures, and retains newly discovered inputs for the next run. PRs read
the corpus but do not publish corpus caches. All runs upload logs, JSON metadata,
corpora, and crash/timeout/OOM artifacts for 30 days, including failed runs.
Cache eviction does not remove the checked-in fixtures or regressions.

Each input has a 20-second limit, except complete surface-knot and degree-elevation verification
and surface editing, which have 60 seconds. All targets retain the 2 GiB process RSS limit. The new
tensor target checks both axes at degree 25 and every raw-support homogeneous
equation; the densest retained input takes approximately 23 seconds with
instrumentation on the development machine. Its larger verification budget is
explicit in each campaign report.

A timeout whose input completes is triaged by measurement, never by
re-running to green. The input is timed locally with AddressSanitizer. The
Linux runner is about 2.6 times slower than the development Mac, so if it
takes more than 20 s ÷ 2.6 there, the target moves to the 60-second
budget; otherwise the runner was slow and the budget stays. Either way the
input is retained under `fuzz/regressions` with its measured times.
`surface_editing` moved on 2026-09-26 (8.2 s locally under the sanitizer).

Surface editing
and surface knot/degree editing permit 4096 bytes to populate full tensor grids; curve knot editing permits 512 bytes;
other targets permit 256
bytes. The modeling harness bounds geometry to 24 vertices and eight operations;
curve spline/Bézier-editing inputs have at most 51 poles; knot editing starts with
up to 101 poles and can insert up to 50 more; surfaces have up to 108 poles with
degree 25 in either direction. The surface-editing target allows up to 784
poles, including degree 25 in both axes (the kernel allows 4096 total poles).
It queries one selected knot rectangle for high-degree grids and also full
domains/multiple turns for degrees whose product is at most 25. Full high-degree
decompositions and all edits remain in the complete native coefficient bridge;
ordinary fixtures retain extraction and composed high-degree edits for all four
periodicity combinations. A Bézier
known-factor intersection input has degree at most eight, and a rational polyline
has up to eight spans. Power-curve quadric inputs reach degree 25 (degree-50
equations), checked by independent monotonic rational powers. Independent exact
fixtures additionally cover dense degree-25 intersections.
The exact edited-curve target starts with at most 26 poles and refines to at
most 52. It checks degree-25 known-factor plane numerators and squared quadric
contacts through order 50, rational secants, periodic line crossings and
exact periodic circle overlaps. Parameter domains include rational offsets and
widths through 2^2048 and 2^-2048; positive rational weights, common homogeneous
scaling and coordinate permutations exercise conditioning independently of
root identity. These cases retain exact results when floating conversion fails.
Spline minimum-distance inputs additionally offset known-factor curves so the
positive squared minimum is exactly `h²`. Their stationary equations can reach
degree 73, while an attained affine-hull bound can certify minima directly. Rational
parameter ties reach degree 25, while positive-distance irrational ties and
composed edits start at degree at most eight, with optional elevation to nine.
The independent factor identities still certify every minimizing parameter;
all original zero-distance degree-25 cases remain. Nonplanar rational cylinder
curves additionally have `D=h²+scale²*s²*(s²-a)²`, with `0<a<1`. Their full 3D
affine hull has an unattainable zero lower bound, so non-singleton ranges retain
the general stationary solver. Independently known roots and coordinates
certify all minima through parameter changes, translations, axis permutations,
degree-five-to-six elevation, knot insertion and closed trims. The input and
RSS budgets are unchanged; the original degree-24 timeout seed remains included.
Spline/linear inputs use at most 128 bytes and the standard 20-second/2 GiB
limits. Their known-factor curves reach degree 25; polylines can be elevated to
degree 25 after optional knot insertion. Each input also runs the reversed
query and typed rejection checks. Early development replay of 3,000 random
inputs found known-factor sub-float root pairs taking up to 65 seconds. Profiling
traced this to schoolbook integer gcds during polynomial content removal and to
late rational-root recognition. The binary gcd and rational-root-theorem
candidate (see [the capability](SPLINE_LINEAR_INTERSECTIONS.md)) reduced that
input to 0.79 seconds in release, without changing any limit or assertion.
General root products reach degree 25. An
outer deadline kills the build/fuzzer process group. Corpus replay without any
subsequent mutation is an incomplete run. An incomplete run, crash,
timeout, OOM, changed dependency lock or mathematical disagreement fails CI.
**Per-push and scheduled runs (U6 of `REVIEW_NOTES.md`).** The daily
schedule replays every retained input, then mutates for 600 seconds, and on
success records the replayed corpus in a manifest (`corpus/<target>.replayed`,
cached beside the corpus). A push or pull request instead runs
`run_fuzz.py --per-push --sample-seed <run id>`: it replays every checked-in
regression, every input absent from the manifest (seeds added since the last
full replay) and `SAMPLE_SIZE` (64) others drawn by `random.Random(seed)` from
the sorted remaining names, then mutates for 60 seconds. The seed is the
workflow's run id and every report records it with the counts of each group,
so `run_fuzz.py --replay sample --sample-seed <id>` on the same corpus and
manifest replays the same inputs. Without a manifest the whole corpus
replays. `analytic_intersections` also keeps allocation stack traces to five frames
(`malloc_context_size=5`): with thirty, AddressSanitizer's stack depot grew
it to 1,489 MB in 120 s against 33 MB without a sanitizer; five frames keep
227 MB and still name each allocation's site. Campaigns pass `-reload=0`: libFuzzer otherwise rereads the corpus
directory every second and reruns every file newer than its first read that
added no coverage, outside the stop-file check. The first per-push runs after
U6 (`88ae24df`) copied a whole corpus into the sample directory with fresh
modification times and overran the shutdown grace on thirteen targets
(exit 124, no findings); the sample now keeps each seed's time as well. A sampled run's new inputs are not retained (they are new only
against the sample); the corpus grows on full replays. `surface_knots`,
`degree_elevation` and `surface_editing` replay only their regressions per
push (`--regressions-only`) and fuzz on the schedule. Acceptance still needs
the full replay green on the schedule run of the accepted revision, plus the
clean local 600-second campaign. Locally, a sampled `brep_io` run replayed 65
inputs in 5.6 s where the full corpus took 126 s.

Retained corpora are minimised weekly (R12 of `REVIEW_NOTES.md`): on Sundays
the fuzzing workflow runs `run_fuzz.py --minimize` instead of a campaign. It
seeds each corpus as usual, merges it into a fresh directory with libFuzzer's
`-merge=1` under the target's own input and 2 GiB limits and startup budget,
and replaces the corpus only when the merge completes; the cache then keeps
the smaller corpus. Inputs that add no coverage are dropped; retained
regressions are re-seeded on every run. Locally the `brep_io` corpus went
from 1,872 inputs to 1,108 in 190 s.

These are harness limits, not kernel production performance guarantees.

B-rep validation inputs use at most 256 bytes and the standard 20-second/2 GiB
limits. Every prism is built by `Solid::extrude` and validated before mutation,
so the base's validity is also re-checked on every input. Since S4d (exact spline
composition, Bernstein fluxes and mass integrals) its replay and mutation churn
BigRationals like the exact-oracle targets: a 180-second local campaign reached
2,064 MB of RSS while the input it stopped on runs alone in 66 ms, so it
joins them below (allocator purge between inputs, 64 MiB quarantine).

The surface-knot runner enables the test-only `asan-allocator` feature together
with AddressSanitizer. At most once per second, after a complete input and all
its mathematical checks have returned, the harness invokes
`__sanitizer_purge_allocator`. The pinned libFuzzer already invokes this API
during mutation, but omits it during seed replay. The additional call covers
replay as well. It drains freed-allocation quarantine and releases unused pages;
it does not free live geometry or change the RSS limit. This does shorten the
quarantine across independent inputs, as the runtime's existing purge does;
the normal sanitizer checks remain active throughout each input. Ordinary
replay binaries and the production kernel have no sanitizer FFI dependency.
See the pinned runtime's `FuzzerLoop.cpp::{ReadAndExecuteSeedCorpora,PurgeAllocator}`
and [LLVM's allocator implementation](https://github.com/llvm/llvm-project/blob/main/compiler-rt/lib/asan/asan_allocator.cpp).

This target also sets `ASAN_OPTIONS=quarantine_size_mb=64`, retaining a 64 MiB
freed-block quarantine while keeping the same 2 GiB process gate. Other targets
use their existing sanitizer settings. Default quarantine plus allocator
cleanup still exhausted macOS RSS during retained-corpus replay despite roughly
38 MB of live allocations at the failure. Ordinary Rust replay of the same 197
inputs twice peaked at approximately 62 MiB. The smaller quarantine is a
documented instrumentation tradeoff: it can miss a stale-pointer access after
its freed allocation leaves quarantine sooner. It does not suppress a reported
error, remove a mathematical assertion, or raise the process memory cap.
Campaign manifests record the effective sanitizer options. These measurements
diagnose this corpus and runtime; they do not establish a general kernel memory
bound.

`spline_linear` uses the same two settings, for the same measured reason. Its
second clean local 600-second campaign at `809f8b4b` stopped with an OOM at
2,055 MB. RSS was already about 1.9 GB when the 535-input corpus replay ended.
The saved OOM input alone peaks at 65 MB under AddressSanitizer and 3.6 MB
without it. The campaign-wide growth is therefore freed-block quarantine and
allocator retention from exact BigInt churn, not live data. The input limit,
2 GiB gate and all assertions are unchanged.

`boolean` uses them too (S9a.1): a 600-second campaign at `d7e515d9` left
its 348-input replay at 1,443 MB and stopped at the 2 GiB gate after 1,125
mutations, on an input that peaks at 6 MB alone; its exact arrangements churn
BigRationals as the others do. `split`, listed since S8d.2, only now calls
the purge between inputs from its target (it had the quarantine setting
alone). With S9a.2's stacks the purge was not enough: a campaign at
`f510f3a6` left its 672-input replay at 1,849 MB and stopped at the gate
with 28 MB live and 53 MB quarantined, the rest AddressSanitizer's stack
depot, so `boolean` keeps five frames per allocation too, as
`analytic_intersections` does.

The mutation timer starts when the pinned libFuzzer reports `INITED`, after
corpus replay. Its `max_total_time` flag includes initialization and previously
allowed a growing corpus to consume the entire short campaign; this was caught
as an incomplete CI run. The runner now creates a nonempty `stop_file` after
the full requested mutation budget. LibFuzzer stops normally and emits final
statistics. Build and corpus replay have a separate deadline of
`min(3600, 600 + input_limit_seconds * (initial_corpus_files + 1))` seconds.
This reserves build time and each input's configured allowance, including the
initial empty input, subject to a one-hour cap. An earlier average-cost estimate
exhausted startup on a 385-input surface corpus. Every saved input remains in
replay, with its own timeout, and mutation still receives its full separate
budget. A late `INITED` marker cannot borrow mutation or
shutdown time. The pinned libFuzzer checks `stop_file` between mutation batches,
each containing up to five callbacks. The runner explicitly keeps
`mutate_depth=5` and allows five times the per-input timeout plus five seconds
for that final batch and its statistics: 105 seconds normally, 305 for surface
knots. The process group is bounded by the startup allowance plus the requested
mutation duration plus that final-batch grace. Individual input limits stay at
20/60 seconds, and ignored stop requests are killed.
An early exit, startup overrun, missing final statistics or absent mutations
still fails the campaign. Reported peak RSS or slowest input exceeding its limit
also fails, even if the runtime exits successfully: a Linux seed took 64 seconds
against a 60-second alarm without a nonzero exit. Separating the phases fixes a reproduced Linux run
that completed its mutation budget but was killed before its last input and
final statistics finished; the failed evidence remains retained.

Surface-knot coefficient checks use exact continuity induction across the full
common raw-support partition. Equality on the preceding span and the known
knot multiplicities prove the lower coefficients; explicit integer comparisons
prove the rest. Periodic curves require all first-span coefficients because
they have no exterior zero polynomial. This remains a complete identity proof,
with adversarial comparisons against the exhaustive coefficient checker; it
does not substitute sampled points or approximate comparisons.

## Local reproduction

On a cargo-fuzz-supported Unix host with a C++ compiler:

```sh
rustup toolchain install nightly-2026-09-22 --profile minimal --component rustfmt
cargo install cargo-fuzz --version 0.13.1 --locked
python3 rust/tools/run_fuzz.py --seconds 300
# Reproduce a particular random campaign start (corpus state also matters):
python3 rust/tools/run_fuzz.py --target intersections --seconds 600 --seed 314159
```

The nightly toolchain is only for fuzzing; normal builds retain Rust 1.85 as
their minimum. `--toolchain` permits deliberate local toolchain experiments.
`target/fuzz-reports/summary.json` records the compiler, fuzzer, lock hash,
revision/dirty state, command, startup/mutation timing, budget completion,
replay/mutation execution counts, coverage edges and artifact
names. The raw log records the random seed. Cargo-fuzz 0.13.1 has no `--locked`
run option: the runner fetches with `--locked`, builds offline, and verifies
that the lock remains unchanged.

To replay and minimize a discovered failure (replace the example path):

```sh
cargo +nightly-2026-09-22 fuzz run intersections rust/fuzz/artifacts/intersections/crash-HASH --fuzz-dir rust/fuzz
cargo +nightly-2026-09-22 fuzz tmin intersections rust/fuzz/artifacts/intersections/crash-HASH --fuzz-dir rust/fuzz -- -max_total_time=120
```

For the surface-knot target, set `ASAN_OPTIONS=quarantine_size_mb=64` and include
`--sanitizer address --features asan-allocator` before the final `--` when
reproducing the campaign's allocator behavior.

Keep the original artifact. Investigate whether the defect is in the kernel,
the oracle, or its input contract. Add the minimized bytes under
`rust/fuzz/regressions/<target>/*.bin`, with a short explanation and a readable
ordinary Cargo regression. The runner automatically seeds those bytes on every
campaign. Fix the underlying defect before changing expected values or bounds.
Minimization is an explicit triage step, not a claim that CI automatically
understands or repairs failures. No lifetime reliability guarantee follows
from any finite campaign.

Slow-unit diagnostics also receive triage even when a campaign passes. The
checked-in [regression notes](fuzz/regressions/README.md) link saved inputs to
their originating run and ordinary exact fixtures. A release replay example
can time complete spline-intersection oracle checks separately from sanitizer
instrumentation; these timings are diagnostic, not a production latency gate.

If corpus replay starts exhausting the outer startup allowance, compact it with
`cargo +nightly-2026-09-22 fuzz cmin <target> --fuzz-dir rust/fuzz` and retain
the uncompressed artifact until the compacted corpus is validated. A stalled
campaign must be repaired, not counted as a successful mutation run.


Surface knot editing starts with up to 784 controls, preserving every complete
transverse coefficient identity after each edit. Its independent checker shares
basis construction and equation factorization across transverse fields rather
than repeating them per row. Every field and residual equation is retained.
Pairs of identical complete before/after control columns share one projection
only after comparing every integer in both columns. Any changed field creates
its own pair; tests corrupt every control component to verify rejection.
The fifteen-target daily workflow includes this campaign with the same startup,
mutation and memory limits; its larger per-input limit is documented above.
