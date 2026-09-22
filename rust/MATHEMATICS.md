# Mathematical foundation

Mathematical contracts define correctness; fuzzing and regression tests challenge
their implementation. OCCT remains a source and behavioral reference, with
independent references used to resolve disagreements. An application's saved
documents are not required to establish kernel correctness.

## Three separate questions

1. **Predicate:** which side of a line or surface is this represented point on?
   Its sign controls topology and needs a reliable decision, including zero.
2. **Construction:** where is the intersection, projection or offset? Computing
   new coordinates can round, be ill-conditioned, or have multiple solutions.
3. **Modeling tolerance:** are two entities close enough for a specific operation?
   This is a dimensional geometric policy, not a patch for arithmetic errors.

This distinction is standard in [CGAL's kernel design](https://doc.cgal.org/latest/Kernel_23/index.html#Kernel_23PredicatesandConstructions).
An exact sign for existing floating-point coordinates does not recover precision
lost when those coordinates were created. In particular, arbitrary rounded
translations/rotations do not necessarily preserve the exact orientation of
nearly collinear represented points.

## Implemented contract: 2D orientation

`predicates::orient2d(a, b, c)` returns clockwise, collinear or counterclockwise
according to the exact sign of

```text
D = (bx - ax)(cy - ay) - (by - ay)(cx - ax).
```

The domain is all finite IEEE-754 binary64 inputs, including signed zero and
subnormals. NaN and infinities return `Error::NonFinite`. There is no epsilon in
this predicate and no heap allocation. The ordinary Rust target arithmetic
assumption is IEEE binary64 round-to-nearest without fast-math reassociation;
nonstandard floating-point modes are not supported.

The usual path uses the determinant filter from
[Shewchuk's robust predicates](https://www.cs.cmu.edu/~quake/robust.html), specifically
the `orient2d` / `ccwerrboundA` calculation in his public-domain `predicates.c`.
The numerical sign is accepted only when its magnitude exceeds
`(3 + 16u)u * (abs(left) + abs(right))`, where `u = 2^-53`.

We conservatively restrict this filter to zero or coordinates with absolute
value in `[2^-400, 2^400]`. Nonzero coordinate differences are then at least
`2^-452` and at most `2^401`; nonzero products and their error bound remain
normal, and no product or sum overflows. A zero, uncertain result, or any input
outside that interval uses the exact fallback. This is deliberately narrower
than the predicate's public input domain.

### Exact fallback and its size bound

The fallback is an independently implemented integer calculation, not a port
of Shewchuk's floating-point expansion routines:

- Decode each finite `f64` into `sign * mantissa * 2^exponent`, with at most
  53 mantissa bits and exponent in `[-1074, 971]`.
- Expand the determinant into six signed coordinate products:
  `ax*by + bx*cy + cx*ay - ay*bx - by*cx - cy*ax`. This avoids rounding coordinate
  differences, including differences that would overflow `f64`.
- Each mantissa product is below `2^106` and fits `u128`. Express every term
  in units of `2^-2148`, requiring a nonnegative shift from 0 through 4090.
- Each shifted product is below `2^4196`. Even six products sum to less than
  `2^4199`, so 66 `u64` limbs (4224 bits) suffice for either signed subtotal.
- Accumulate positive and negative terms separately with integer carries, then
  compare from the highest limb. No subtractive floating-point cancellation,
  underflow or overflow enters this path. The input size and work are bounded.

The module and its tests document this reasoning. It is not a formal verification
of all compiled machine code or a correctness proof for the rest of the kernel.

### Where it is used

Polygon segment-crossing decisions use exact orientation signs. Ray classification
compares orientation and edge direction instead of constructing a rounded ray/edge
intersection. Euclidean proximity tests still apply the caller's linear tolerance
separately. Area/moment integration and point-to-segment distances remain ordinary
floating-point calculations with their existing validation and limits.

OCCT's [CSLib_Class2d](../src/FoundationClasses/TKMath/CSLib/CSLib_Class2d.cxx)
was consulted for the crossing and boundary-classification behavior. Its
normalization, tolerance treatment and acceleration grid are not ported here.
Rust's exact predicate is a deliberate numerical implementation choice; we do
not claim identical internal decisions to every OCCT version on all degeneracies.

## Implemented contract: spatial predicates

`orient3d(a,b,c,d)` returns the sign of
`((b-a) × (c-a)) · (d-a)`. Positive means above the oriented `abc` plane:
`orient3d(origin, X, Y, Z)` is positive. This is the **opposite sign** to
Shewchuk's `orient3d` convention. A degenerate defining plane returns coplanar.

The filter follows Shewchuk's `orient3d` / `o3derrboundA` calculation, with bound
`(7+56u)u * permanent` and the final sign negated. It is restricted to zero or
coordinate magnitudes in `[2^-200, 2^200]`. Nonzero coordinate differences are
at least `2^-252` and at most `2^201`; the three-factor arithmetic and its error
bound remain normal and finite. Uncertain results fall back to integers.

`in_sphere(a,b,c,d,q)` returns inside, boundary or outside the unique sphere
through `a,b,c,d`, independently of defining-point order. Coplanar defining
points return `Error::Degenerate`. It evaluates an exact lifted determinant;
it never constructs a rounded center or radius. There is currently no floating
filter on this predicate. Both APIs accept every finite binary64 coordinate,
reject non-finite values and apply no modeling tolerance.

### Integer representation and bounds

`exact.rs` represents each finite input as the signed integer `I(x)=x*2^1074`.
`I(x)` has fewer than 2,099 bits; a difference has at most 2,099 bits. Exact
vector differences, cross products and determinants then operate on integers.
The orientation determinant fits below `2^6300`; the degree-five sphere
determinant fits below the conservative bound `2^10510`. Intermediate sizes
are bounded by fixed input dimension, binary64 exponent range and polynomial
degree, with no user-controlled precision or convergence iteration.

These operations use `num-bigint` (locked to 0.4.8, default features disabled)
instead of introducing a second hand-written multi-precision implementation.
It and its integer-trait dependencies are Rust libraries; no native geometry
SDK is linked. The original allocation-free 2D predicate remains unchanged.

## Implemented contract: certified linear intersections

`Plane3::through_points` and `Triangle3::new` retain three finite, exactly
noncollinear defining points. Their exact integer normal defines the plane;
normalizing a rounded vector cannot move it. The supported operations are:

| Operation | Possible geometric results |
| --- | --- |
| `line_plane(p,q,plane)` | Disjoint parallel line, contained line, unique point |
| `segment_plane(p,q,plane)` | Disjoint, contained segment, unique point |
| `segment_triangle(p,q,triangle)` | Disjoint, unique point, coplanar overlap segment |

Segments and triangles are closed, including endpoints/edges. Equal segment
endpoints are a point at parameter zero; equal line-defining points are an
error. Triangles include coplanar edge overlap and vertex tangencies. Overlap
endpoints follow the input segment direction; triangle winding is irrelevant.

For plane determinant values `Dp` and `Dq`, the exact crossing parameter is
`t=Dp/(Dp-Dq)`. Zero and sign tests decide parallelism, containment and segment
membership before rounding. For triangle membership, drop a coordinate whose
exact normal component is nonzero and evaluate the three exact edge signs.
Coplanar segments are clipped against those half-planes with rational parameter
comparisons. No approximate intersection is fed back into a predicate.

For each coordinate, construct the integer numerator
`I(p_i)*den + (I(q_i)-I(p_i))*num`, with a positive denominator and scale
`2^-1074`. Coordinate numerators fit below a conservative `2^8410`; all rational
comparisons have bounded integer size. At most 63 exact comparisons over the
ordered positive finite binary64 bit patterns find the enclosing interval;
negative values use sign symmetry. No floating quotient is used.

`IntersectionPoint::bounds()` contains the exact mathematical point. Each
coordinate interval and `parameter()` is either a single exactly representable
value or two adjacent finite binary64 values. `position()` is a finite rounded
representative inside those bounds. It is **not** a new exact point guaranteed
to satisfy both input equations. Consumers must propagate the enclosure, refine
their representation, or explicitly apply a modeling tolerance. Repeatedly
discarding bounds and reclassifying rounded points loses this guarantee.

An infinite-line coordinate or affine parameter outside finite binary64 range
returns `Error::Unrepresentable`; it is never replaced by infinity or a guessed
parallel result. Segment constructions stay within their finite endpoints and
`t ∈ [0,1]`. Thin angles do not cause false parallel classifications. Enclosures
certify the represented inputs, not an uncertain physical measurement.

OCCT's `IntAna_IntConicQuad::Perform(gp_Lin,gp_Pln)` supplies the line/plane
reference and affine-substitution behavior. Its angular-tolerance parallelism
differs deliberately from these exact primitives. The live native comparison
checks 72 well-conditioned cases (intersection, parallelism, containment,
oblique planes, scales and endpoint parameters). Full exponent range,
near-parallelism and coplanar triangle clipping use independent rational
oracles; we do not claim parity with OCCT's tolerance decisions there.

## Implemented contract: quadratic algebraic roots

`polynomial::quadratic_roots(a,b,c)` solves `a*t²+b*t+c=0` for all finite
binary64 coefficients. Exact coefficient zeros determine the degree. The zero
polynomial returns `All`; a nonzero constant returns `None`; a linear polynomial
has one simple root. For a quadratic, the exact discriminant determines zero,
one double, or two distinct real roots. A double root appears once with
multiplicity two. There is no near-zero threshold or root-merging tolerance.

After making `a` positive, the roots are represented as
`(-b ± sqrt(b²-4ac))/(2a)`. This formula is **not evaluated in floating point**.
`algebraic.rs` retains expressions `(n+k*sqrt(D))/d`, with nonnegative integer
`D` and positive integer `d`. To compare against an exact rational value, clear
the positive denominator and examine `A+B*sqrt(D)`. Equal-sign terms have that
sign; opposite-sign terms are ordered by comparing `A²` and `B²*D` exactly.
Zero is decided exactly as well. Removing a common power of two from polynomial
coefficients reduces work without changing their roots.

`QuadraticRoot` retains this algebraic identity and multiplicity. `compare(x)`
orders it against any finite binary64 value. `bounds()` returns the smallest
finite binary64 enclosure, or `Unrepresentable` if the root exceeds that range.
An unrepresentable root still exists as an exact object and can be compared.
Distinct roots remain distinct even if their enclosures overlap. These closed
binary64 enclosures are not necessarily disjoint isolating intervals; general
algebraic root isolation is a separate capability (see the distinction in
[CGAL's algebraic kernel](https://doc.cgal.org/latest/Algebraic_kernel_d/classCGAL_1_1Algebraic__kernel__d__1.html)).
No cubic/quartic solver or general polynomial API is claimed.

## Implemented contract: curved primitive intersections

`Sphere3`, `Cylinder3` and `Circle3` accept finite centers, strictly positive
finite radii, and (where relevant) finite nonzero axis/normal vectors. Vectors
are interpreted exactly without normalization. A cylinder is its infinite
lateral surface, without end caps. A circle lies in the plane through its
center perpendicular to its supplied normal. These are query primitives;
they do not add spherical solids, trimmed cylinder faces, or arc profiles.

All three support infinite-line and closed-segment intersections. With
`w=p-center`, `d=q-p`, radius `r`, and cylinder axis `v`, substitute into:

```text
sphere:    |w+t*d|² - r² = 0
cylinder:  |(w+t*d) × v|² - r²*|v|² = 0
circle:    sphere equation AND (w+t*d) · normal = 0.
```

The coefficients are constructed with exact integers, including differences
that would overflow or lose precision in binary64. A coplanar line/circle
query uses the sphere quadratic; a line crossing the circle plane has one
rational candidate whose sphere equation is checked exactly before construction.
Segment membership compares exact algebraic parameters with zero and one before
any output bounds are calculated. Thus an out-of-segment unrepresentable root
cannot discard an otherwise valid hit.

Results distinguish empty, contained cylinder generator, one hit and two hits.
Each hit identifies a simple crossing, a true tangency, or a zero-length segment
on the primitive. A zero-length line is an error; an equal-endpoint segment is
a point at parameter zero. A segment entirely inside a sphere/cylinder with no
surface contact returns empty: this is surface intersection, not volume overlap.

Each point coordinate uses the exact affine expression `p_i+t*(q_i-p_i)` and
is bounded directly with integer comparisons, without substituting a rounded
parameter. This matters at large offsets: two parameter intervals may coincide
while the exact hit coordinates remain distinct and fully resolvable. Both
hits are retained, ordered by their exact parameter. Coordinates/parameters
outside finite binary64 range fail the entire requested construction explicitly.
Output representatives must not replace their uncertainty bounds in later
topological decisions.

All input dimensions and polynomial degrees are fixed. Cylinder coefficients
are degree four in lattice integers with at most 2,099-bit coordinate
differences; a conservative coefficient bound is 8,410 bits and discriminants
need fewer than 16,830 bits. After affine construction and binary64 comparisons,
the largest squared-comparison terms fit within a conservative 23,200-bit
bound. Each enclosure uses at most 63 bisection iterations. This bounds integer
work structurally; it is not a measured production latency or allocation limit.

OCCT reference paths are `IntAna_IntConicQuad`, `IntAna_Quadric`,
`gp_Sphere::Coefficients`, `gp_Cylinder::Coefficients`,
`IntAna2d_AnaIntersection::Perform(line,circle)` and
`math_DirectPolynomialRoots::Solve(A,B,C)`. OCCT's coefficient thresholds,
discriminant uncertainty band and radius epsilon can intentionally classify
near-degenerate cases differently. The live comparison covers 174 cases in the
common well-conditioned domain, including exact tangencies and native root
multiplicities. Independent mathematical oracles define the extreme and
near-tangent behavior; native results are never used to weaken those contracts.

## Independent and adversarial evidence

- `fixtures/orient2d.tsv`: 2,417 input triples with expected signs calculated by
  Python `fractions.Fraction.from_float`, using an exact rational determinant.
  Input bits are preserved. Cases include cancellation, neighboring floats,
  exponent boundaries, duplicates, signed zero, subnormals and maximum finite
  magnitudes. Cargo tests check all six point permutations: 14,502 comparisons.
- Another 10,000 seeded integer triples use a separate `i128` oracle and check
  exactly representable power-of-two scaling/translation and reflection.
- Unit tests exercise the certified fast path, required fallback paths, maximum
  accumulator carries and rejection of non-finite coordinates.
- `fixtures/predicates3d.tsv`: 824 orientation and 824 sphere cases from Python
  `Fraction` matrix elimination, including exact and perturbed coplanarity,
  sphere-boundary queries and extreme exponents. Orientation checks all 24
  permutations; sphere tests check four defining-point orders.
- `fixtures/intersections.tsv`: 963 line/plane, segment/plane and
  segment/triangle cases, with minimal coordinate and parameter enclosures
  calculated by a separate rational barycentric oracle and Python's rational
  conversion. Includes degenerate, unrepresentable, point, overlap and empty
  results. Tests also reverse segment endpoints and all triangle vertex orders.
- 256 generated prisms over 25 scales, with concave/convex radial polygons and
  optional holes, check volume/first-moment conservation under splitting,
  introduced cut area, winding and offset reversal, rigid-motion covariance,
  nonnegative inertia quadratic forms, topology, sampled bounds and boundary
  classification. Each has 32 further classification probes before/after motion.
  Their seed and case index make failures reproducible.
- `fixtures/quadratic.tsv`: 448 root cases from exact polynomial evaluation and
  vertex ordering in Python `Fraction`, independent of production radical
  comparisons. Covers degree reduction, repeated roots, near-double roots,
  cancellation, full exponent range and unrepresentable roots.
- `fixtures/curved.tsv`: 1,092 line/segment queries against circles, spheres and
  cylinders. The reference uses rational axial projection for cylinders and
  exact polynomial-sign comparisons for root/coordinate bounds. Separate tests
  cover adjacent-float tangencies, coincident parameter enclosures, endpoint
  clipping before overflow, direction rescaling and endpoint reversal.

These fixtures are bounded deterministic tests. Separately, [coverage-guided
fuzzing](FUZZING.md) mutates predicates, linear/curved intersections, polynomial
roots and standalone modeling
sequences under sanitizers, with retained corpora and daily campaigns.
Finite test evidence supplements the arithmetic argument; it does not prove a
general kernel correct. The OCCT comparison corpus remains a separate oracle.

```sh
cargo test --workspace --locked
cargo test --workspace --locked --release
python3 rust/tools/generate_predicate_fixtures.py --check
python3 rust/tools/generate_spatial_fixtures.py --check
python3 rust/tools/generate_curved_fixtures.py --check
```

Native CI executes both debug and optimized Rust tests. The rational fixture
generator is checked independently in the Linux validation job. Its default
write mode is for deliberate, reviewed corpus additions, never for making a
disagreement disappear.

## Next mathematical work

1. Add the distance/incircle comparisons required by subsequent algorithms.
   Preserve explicit units, domains and arithmetic bounds; do not reuse an
   arbitrary global epsilon.
2. Specify error and conditioning contracts for projections, curve/surface
   evaluation, root isolation and intersections. Use interval/error bounds or
   additional precision when ordinary arithmetic cannot establish the answer.
3. Strengthen topology invariants and tolerance propagation through each
   operation; sampled checks alone cannot certify full curve/surface agreement.
4. Extend the existing coverage-guided campaigns to each new capability,
   minimize failures, and expand resource/performance and independent-oracle
   checks as geometry and operation sequences become more complex.

General curve/surface intersections, higher-degree root isolation, Booleans,
generic topology history, STEP and meshing remain unimplemented. Certified
linear and quadratic primitives do not establish these capabilities or make
the rest of the kernel exact.
