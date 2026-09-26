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
This specialized radical API remains degree two; the general polynomial API
below handles degrees through 25 without closed-form cubic/quartic formulas.

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

## Implemented contract: rational spline curves

`curve::BSplineCurve3::new` stores a nonperiodic curve of degree 1..=25 with at most
4096 finite 3D poles. Its distinct knots must be finite and strictly increasing;
interior multiplicities are 1..=degree and end multiplicities 1..=degree+1.
The sum of multiplicities is `number_of_poles + degree + 1`. Weights default to
one; supplied weights must be finite and strictly positive, including positive
subnormals. Equal weights reduce geometrically to a polynomial curve, without
discarding unequal weights based on an epsilon.

With zero-based expanded knots `U` and `N` poles, the closed evaluation domain
is `[U[degree], U[N]]`, which must have positive width. This includes unclamped
knot vectors; the first/last distinct knots need not bound the domain.
`BezierCurve3` uses degree `N-1` and clamped knots 0,1. Nonperiodic
extrapolation and spline B-rep edges are not yet implemented. Exact Bézier
extraction/editing and rational curve knot editing are described below. These geometric
primitives do not change the current prism solid model.

The [de Boor recurrence](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/de-Boor.html)
operates on homogeneous poles `(w*x,w*y,w*z,w)`. Rust converts the represented
inputs to exact reduced rationals before products, differences or divisions.
For each local interpolation `H=(1-a)*L+a*R`, where `a` is affine in the original
parameter, the evaluator carries these derivatives:

```text
H'  = (1-a)*L'  + a*R'  + a'*(R-L)
H'' = (1-a)*L'' + a*R'' + 2*a'*(R'-L').
```

The positive-width active span is contained in every knot interval used by
the recurrence, so each divisor is strictly positive even at repeated knots.
Positive weights and the nonnegative partition of unity make the evaluated
homogeneous weight `W` strictly positive throughout the domain. For each
coordinate `X`, divide only after interpolation:

```text
P   = X/W
P'  = (X' - W'*P)/W
P'' = (X'' - W''*P - 2*W'*P')/W.
```

These exact identities give derivatives with respect to the original knot
parameter, not normalized arc length. They agree with the rational-derivative
contract studied in OCCT's `PLib::RationalDerivative`. The
[basis derivative identity](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-derv.html)
is used by the independent oracle, which evaluates basis functions instead of
interpolating poles and uses closed quotient formulas for the first two derivatives.

`DerivativeOrder` requests position, first or second derivative. Every requested
component gets its smallest finite binary64 enclosure, using exact rational
comparisons with at most 63 bisection iterations over ordered magnitude bits.
Unrepresentable requested derivatives return `Error::Unrepresentable`; a
position-only query can still succeed. Positive-weight positions lie in the
finite pole convex hull. Rounded output representatives remain approximate
and must not be reclassified as exact spline points.

At an interior knot, `KnotSide::Automatic` checks exact left/right jets when
the knot multiplicity permits a discontinuity in the requested derivatives.
Different derivatives return `DiscontinuousDerivative`. Exact agreement is
accepted even when the knot multiplicity alone would only guarantee lower
continuity. `Left` and `Right` explicitly choose one-sided values; selecting a
side outside the closed domain returns `OutOfDomain`. At a domain endpoint,
automatic evaluation takes the only interior limit. Position-only evaluation
is continuous under the permitted multiplicities.

The degree bound gives at most 325 local interpolation updates per one-sided
evaluation, involving at most 26 active poles; a discontinuity check can evaluate
both sides. Arithmetic uses `num-rational` 0.4.2 with `num-bigint`, independent
of ordinary floating-point intermediate ranges. Rational operand sizes depend
on the knot data and degree. These input/iteration bounds do not establish a
production latency or allocation budget; performance gates remain open.

## Periodic curves and rational tensor-product surfaces

`BSplineCurve3::new_periodic` and `KnotVector::new_periodic` require matching
end multiplicities `m` in `1..=degree`. Pole count is the sum of multiplicities
minus `m`; both periodic and nonperiodic directions require more poles than
their degree. This is a conservative supported-domain restriction, including
for periodic inputs that OCCT's pole-count formula might otherwise admit.

The fundamental domain is `[first knot, last knot]`, with both endpoints
representing the seam. Extend each end by `degree+1-m` knots, shifted by the
exact period, and cycle the poles in OCCT's order. Extended knots may be outside
finite binary64 range: they remain exact rationals internally. All finite
parameters are reduced as `u - floor((u-start)/period)*period` in exact rational
arithmetic. No floating division, remainder, overflow or tolerance snapping
chooses the span. The reduced parameter need not itself be representable by a
single `f64`. At the seam, Left selects the preceding period and Right selects
the following period. Automatic derivatives require exact agreement between
those jets; a closed position alone does not establish derivative continuity.

`BSplineSurface3` takes two independently validated `KnotVector`s and a U-major
control grid (`index = u*v_count+v`). Degrees are independently 1..=25, with at
most 4096 total finite poles and optional strictly positive finite weights.
`BezierSurface3` represents the clamped `[0,1]²` special case. Neither implies a
trimmed face, a solid, a projection/intersection operation or a regular surface:
singular and collapsed patches are valid evaluation inputs.

Evaluation first applies the differentiated homogeneous de Boor recurrence in
V, then U, retaining `(00,10,01,20,02,11)` partials. The recursive multivariate
quotient identity derived from `X=W*P` is:

```text
P_(a,b) = [X_(a,b) - sum binomial(a,i)*binomial(b,j)*W_(i,j)*P_(a-i,b-j)] / W
```

The sum covers `(i,j) != (0,0)`, `i<=a`, `j<=b`. Lower total orders are
calculated first. In particular, `P_uv=(X_uv-W_uv*P-W_u*P_v-W_v*P_u)/W`.
Positive tensor basis weights keep `W>0`. First order returns Du and Dv;
second order additionally returns Duu, Dvv and Duv. Every requested component
has its minimal finite enclosure. Queries fail atomically if a requested
component is unrepresentable. Derivatives use the original parameter units.

Each parameter direction has its own side selection. At intersecting knots,
automatic evaluation compares every relevant quadrant, including both periodic
seams. A mixed-derivative discontinuity cannot be hidden by agreeing position
and first derivatives. High multiplicity does not force rejection when the
actual requested jets agree exactly.

Each surface jet uses at most 26 V recurrences plus one U recurrence, each with
at most 325 homogeneous interpolation updates and six jet entries per component.
At most four side combinations are evaluated. Degree, pole and iteration limits
bound combinatorial work, not a production latency/allocation budget for exact
rational operands. Operational performance and cancellation remain open gates.

## Higher-degree roots and spline/plane intersections

`polynomial::Polynomial` accepts up to 26 finite binary64 coefficients in
ascending power order. Coefficients are exact represented numbers; no epsilon
reduces the degree. The zero polynomial has every point of the queried domain
as a root. Nonzero results are distinct real roots in exact increasing order,
with multiplicity retained. Queries cover the real line or a finite closed
interval, including endpoints. Root identity survives outside binary64;
only requesting a finite enclosure can fail as unrepresentable.

The solver clears denominators, removes positive integer content, takes
square-free gcd layers and isolates roots with a Sturm sequence. Every
pseudo-division elimination scales by a positive magnitude, preserving sign
variations even with negative leading coefficients. Gcds and Sturm chains use
a magnitude subresultant sequence (Collins; Cohen, *A Course in Computational
Algebraic Number Theory*, Algorithm 3.3.1). Each pseudo-remainder is scaled by
exactly `|lead|^(δ+1)` and divided exactly by `|g|·h^δ`. The terms differ from
classical subresultants only by sign, so these divisions are exact and every
chain element is a positive multiple of the Euclidean remainder. No integer
content is computed per remainder. A test checks every element's primitive
part and sign against the previous primitive Euclidean chain. Open-interval root counts
exclude endpoints, which are emitted once. An accepted isolating interval also
has nonroot endpoints. Iterative traversal avoids recursion depth from clustered
roots. Repeated gcd layers establish multiplicity.

`AlgebraicRoot::sign_at` evaluates another represented polynomial at the exact
root. Rational witnesses and sign-definite rational intervals are fast paths;
otherwise the Sturm–Tarski query uses the signed remainder sequence of P and
P′Q. This decides exact zero as well as positive/negative signs. The reference
is Li, Passmore and Paulson,
[Deciding Univariate Polynomial Problems Using Untrusted Certificates in Isabelle/HOL](https://wenda302.github.io/assets/pdf/rcf_jar.pdf),
§§5.1–5.3. This implementation is not formally verified. Retaining a defining
polynomial and isolating interval follows the algebraic-number model described
by the [CGAL Algebraic Kernel](https://doc.cgal.org/latest/Algebraic_kernel_d/index.html).

`AlgebraicRoot::compare_root` orders exact roots across different defining
polynomials, independently of multiplicity or binary64 views. Let beta have
square-free defining polynomial Q and a nonrational isolator (a,b). First
compare alpha exactly with a and b. If alpha is inside, Q(alpha)=0 proves
equality. Otherwise, Q changes sign exactly once on this interval, so the sign
of Q(alpha) equals Q(a) precisely when alpha < beta. Rational beta uses the
existing rational comparison. This gives a terminating equality decision;
repeated approximate refinement is not needed to separate true ties.

The independent comparison fixtures use irreducible factor identities and
canonical real-root indices to recognize equality, then VAS interval refinement
to separate distinct values. The 130 equation pairs include shared factors,
repeated roots, degree-25 equations, coefficients at binary64's exponent limits,
roots outside binary64, and distinct roots whose binary64 views coincide.
The `roots` fuzz target additionally compares each known-factor root against an
independently constructed minimal equation and checks cross-equation ordering.
General rational-function images and global curved-distance minimizers require
additional machinery and are not implied by root ordering alone.

`intersection::spline_plane` substitutes each rational span into a three-point
plane using exact homogeneous coefficients on local parameter [0,1]. Positive
weights ensure a nonzero denominator. All numerator roots are isolated; zero
polynomials become maximal closed overlap intervals. Overlap endpoints are not
emitted again as isolated points. Shared-knot hits merge only by exact knot
identity; proximity and equal floating enclosures never merge distinct roots.

Each point retains the algebraic parameter and homogeneous coordinates. Exact
sign queries against X−xW yield minimal coordinate bounds; original parameter
units are restored exactly. Actual plane sides decide crossing versus same-side
tangency; queried domain endpoints are boundary contacts. Nonsmooth knots retain
separate left/right contact orders. Periodic queries cover the closed fundamental
domain with distinct start/end parameter events. This convention and observed
OCCT differences are [explicitly reviewed](NATIVE_SPLINE_PLANE_DIVERGENCES.md).

`RootIsolationOptions::max_subdivisions` defaults to 65,536, shared across all
spans of a spline query. Exhaustion returns `Error::ComputationLimit` without
publishing partial results. This limits subdivision only; exact gcd/sign work,
coefficient bit growth and allocations do not yet have production latency/memory
budgets or cancellation. Degree/input limits bound combinatorial work, not
production response time.

`spline_plane_in` queries an explicit finite closed interval `[first,last]` with
positive length. Nonperiodic bounds must stay inside the curve domain. Periodic
queries can cross seams, start outside the fundamental period, and span multiple
turns. Parameters are retained in the caller's original units. Reversed bounds,
zero-length intervals, nonfinite bounds and nonperiodic extrapolation are typed
errors; curve sense reversal is a separate operation and is not implied.

Shifted periodic knots are rational numbers and need not be representable by
binary64. Span translation, range clipping, endpoint equality and ordering
therefore remain exact. Roots are isolated inside each clipped interval; they
are not accepted/rejected by rounded output bounds. Query endpoints have only
their interior contact side. A zero span contributes its clipped positive-length
overlap; when a query only touches its endpoint from a nonzero neighboring span,
that endpoint becomes an isolated boundary event.

`SplinePlaneOverlap` now retains exact rational endpoints and minimal enclosures.
`parameters()` supplies rounded representatives; `parameter_bounds()` and
`compare_parameter()` preserve reliable decisions. Separate overlaps or isolated
hits can have identical binary64 enclosures and remain distinct. This owning
type is `Clone` rather than `Copy`.

`SplinePlaneOptions` adds a default limit of 4,096 visited spans, independently
of the root-subdivision limit. Before enumerating a periodic interval, the exact
number of translated spans is counted: for base span `(a,b)` and period P, its
count is `ceil((last-a)/P) - floor((first-b)/P) - 1`. Exceeding the budget returns
`ComputationLimit` before span enumeration or any root solving. Even queries
from `-f64::MAX` to `f64::MAX` on tiny periods cannot start an unbounded loop.
This traversal bound still does not constitute a hard CPU or allocation limit.

## Spline/sphere and spline/cylinder intersections

`spline_sphere` and `spline_cylinder` query the complete closed fundamental
domain. Their `_in` variants accept the same explicit closed ranges as plane
queries, including multiple periodic turns. The `_with_options` variants use
`SplineSurfaceOptions` for span traversal and root-subdivision limits. A cylinder
is infinite and has no caps; these queries do not construct a spherical solid
or test overlap with the interior of a volume.

For homogeneous coordinates `(X,W)`, center/origin `c`, radius `r`, and a
nonzero cylinder axis `v`, let `D=X-c*W`. The exact scalar polynomials are:

```text
sphere:   D·D - r² W²
cylinder: (v·v)(D·D - r² W²) - (v·D)².
```

Positive spline weights make `W` strictly positive. The cylinder's `v·v` is
strictly positive, so clearing these denominators introduces no roots, and
the polynomial sign agrees with outside/inside surface sidedness. Unlike a
rounded unit-axis construction, this accepts any finite nonzero represented
axis, including subnormal components. Degree-25 curves produce polynomials
through degree 50. The internal exact isolator supports these equations;
the public `Polynomial` coefficient-input limit remains degree 25.

`SplineSurfaceIntersection` shares exact point identities, minimal coordinate
and parameter enclosures, one-sided contact orders and maximal overlap
intervals across planes, spheres and cylinders. The prior `SplinePlane*`
names remain aliases. Order 50 contacts are retained; a constant curve on a
surface occupies its entire parameter interval. Overlap endpoints are excluded
from isolated points, and distinct periodic parameter events are not merged
because they revisit the same position. Outside-domain, nonfinite, traversal
and subdivision failures remain atomic.

Before repeated coordinate sign queries for sphere/cylinder equations, each
nonrational root interval is tightened by at most 128 exact rational bisections.
Plane queries skip this unconditional refinement and use the adaptive sign
filter described below.
Its square-free defining
polynomial has exactly one simple root in that interval, so opposite endpoint
signs certify the retained half. This only improves the interval sign filter;
an undecided sign still uses Sturm–Tarski, including exact zero decisions.
The refinement has no tolerance or output-rounding effect. Dense degree-25
quadric fixtures exposed the need for this filter; those cases remain in the
ordinary debug/release tests. Span/subdivision limits still do not promise a
hard CPU or allocation bound.

Interval-Horner sign filtering carries integer numerators over a shared positive
denominator. For endpoint denominator LCM `D`, write `[a,b]=[A/D,B/D]`.
Starting with the leading coefficient over denominator 1, each Horner step
forms the minimum and maximum of `L*A,L*B,H*A,H*B`, then adds `c*D^k`,
where the previous bounds had denominator `D^(k-1)`. These are exactly the
same rational interval bounds as reduced-fraction arithmetic. Their signs
need no division because `D^k>0`. This avoids repeated large GCDs while retaining
all four products for intervals of either sign or crossing zero. Zero-containing
bounds still fall back to Sturm–Tarski. Two saved fuzz inputs exposed the cost:
a subnormal trim endpoint on a degree-seven plane query, and the sphere
intersection of `(2*t^25,0,0)`. Both remain in the corpus and exact fixtures;
this optimization changes neither tolerance nor the mathematical answer.

If the initial interval filter cannot decide a sign, algebraic sign queries now
try up to 64 exact sign bisections on a local copy of the isolator and repeat the
filter. The same single-simple-root invariant justifies this step. Unresolved
signs still use the complete Sturm–Tarski calculation on the refined interval;
no root identity, multiplicity, bound or exact-zero decision is approximated.
These bounded filter refinements are separate from root isolation's subdivision
budget. A saved degree-25 fuzz polynomial with a full-exponent query at its
root `-1` exercises the former expensive fallback and remains in the independent
root fixtures and fuzz corpus.

## Exact proximity of linear sets

`proximity::closest_points` supports every ordered pair of a point, infinite
line, closed segment, infinite three-point plane, or closed triangle. Finite
binary64 defining points are interpreted exactly. Zero-length segments are
singletons; zero-length lines and collinear planes/triangles are rejected.
No modeling tolerance, rounded normal, or near-parallel threshold is used.

The result retains one rational point on each operand, original affine
parameters and exact squared distance. Line/segment parameters use `a+t(b-a)`;
plane/triangle parameters use `a+u(b-a)+v(c-a)`. Segment parameters are in
`[0,1]`; triangle parameters satisfy `u>=0`, `v>=0`, `u+v<=1`. The API does
not claim uniqueness or enumerate a continuum of minima. Exact comparisons
and the rational result remain usable even when a point, parameter, squared
distance or unsquared distance has no finite binary64 enclosure. Each requested
conversion fails separately with `Unrepresentable`. Rounded representatives
are not exact incidence witnesses. Distances have input length units; squared
distances have squared units. There is no B-rep, curved-set or clearance-volume
query in this API.

The solver enumerates the nonempty faces of each bounded simplex (one point,
three segment faces, seven triangle faces), or the whole unbounded affine set.
For each pair, minimize `||delta + D*x||²`, where columns of `D` are the first
face directions and the negatives of the second face directions. Solve the
normal equations `DᵀD*x = -Dᵀdelta` exactly, with free variables set to zero.
Input coordinates have power-of-two denominators. Clearing their common
denominator leaves integer normal equations. Fraction-free
[Bareiss elimination](https://www.ams.org/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf)
retains integer coefficients; every division is exact. The last nonzero pivot
is a nonsingular-subsystem determinant, supplying a common denominator for
back substitution by Cramer's rule. Feasibility, witness construction and
candidate comparison keep integer numerators over positive shared denominators.
Only the selected result is converted to reduced rational values. This avoids
repeated large fraction reductions without changing rank, minima or tie order.
This system is consistent: `ker(DᵀD)=ker(D)`
and the right-hand side is orthogonal to that kernel. Keep solutions with
feasible bounded-face barycentric coordinates, then select the smallest exact
squared distance. A fixed tie order makes repeated queries deterministic.

Completeness includes singular systems. Take a minimum on its smallest active
faces. It is stationary on their affine hulls. If a null direction changes
bounded barycentric coordinates, move within the solution set until one reaches
a boundary, lowering the face dimension without changing distance. Repeat.
On the resulting enumerated faces, every null direction leaves bounded
coordinates unchanged, so the solver's particular solution is feasible.
When both faces are unbounded there are no feasibility inequalities. Collapsed
segments also have their singleton faces. Thus at least one enumerated
candidate realizes the global minimum. There are at most 49 systems of at most
four variables; arithmetic sizes depend on input bits. This finite bound is
not a measured production latency guarantee.

The independent witness checker uses first-order convex optimality, not another
face search or elimination. This is the supporting-hyperplane criterion in
Boyd and Vandenberghe's [Convex Optimization, §4.2.3](https://www.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).
For returned `p` and `q`, let `d=p-q`. Check membership exactly, then check
`d·(v-p)>=0` at every vertex of the first bounded set and `d·(v-q)<=0` at
every vertex of the second. For each unbounded set, check orthogonality to
every direction. These conditions extend to all feasible points by convexity
or affine linearity. For any other feasible `x,y`, put `h=(x-p)-(y-q)`.
Then `||x-y||²-||p-q||² = 2*d·h + ||h||² >= 0`, certifying global
minimality. Tests also reconstruct both points from the returned parameters.

Python `Fraction` fixtures independently use closed projections, cross-product
line distances, plane intersection, edge halfspaces and boundary candidates.
They check the exact squared distance and minimal enclosures. The square-root
reference seeds its enclosure using integer square root on a subnormal lattice;
production searches binary64 bit order with exact squared comparisons. Native
OCCT observations are a third check with [explicitly reviewed differences](NATIVE_PROXIMITY_DIVERGENCES.md).

## Complete intersections of linear convex sets

The [public contract](LINEAR_INTERSECTIONS.md) extends incidence to the entire
closed intersection for all 25 point/line/segment/plane/triangle pairings.
`intersection::LinearPrimitive3` reexports the proximity input type. The output
`LinearIntersection` distinguishes empty, point, segment, filled convex polygon,
infinite line and infinite plane. `ExactPoint3` retains rational coordinates;
requesting minimal binary64 bounds or a representative is separately fallible.
A result may therefore be usable for exact incidence despite an unrepresentable
coordinate. No topology identifiers, face orientation or history are inferred.

Each input becomes exact affine equalities plus closed linear inequalities.
Points fix three coordinates. Lines fix two independent equations; segments
also bound a coordinate with nonzero direction. Planes have one nonzero normal
equation; triangles add the three inward edge halfspaces. Normals are rational
cross products of the represented inputs, never rounded unit vectors.

Exact elimination of the combined equations gives an inconsistent system or
an affine space `x = o + B t` of dimension zero, one or two. Substitution reduces
the remaining inequalities to that space:

1. In dimension zero, test all inequalities at the unique point.
2. In dimension one, intersect every exact closed parameter bound. Equal bounds
   are one point; reversed bounds are empty. If any operand is bounded, both
   ends are bounded. Otherwise there are no inequalities and the entire line
   remains. A ray cannot arise from the supported operands.
3. In dimension two with no inequalities, the result is the entire plane.
   Otherwise a triangle operand makes the feasible polyhedron bounded. Every
   extreme point lies at two independent active boundary lines; enumerate
   all pairs (at most fifteen), solve each exactly and retain it only if all
   inequalities hold. A nonempty bounded polyhedron is the convex hull of its
   extreme points, including when it collapses to a point or segment. Thus an
   empty candidate list proves emptiness, and taking their hull gives the
   complete set rather than selected contact samples.

The finite hull removes exact duplicate and collinear vertices. A polygon's
cycle begins at its lexicographic minimum and selects the smaller direction;
it is not an oriented B-rep face. Segments use sorted endpoints. Lines set their
first nonzero direction coordinate to one and that origin coordinate to zero;
planes normalize their first nonzero normal coefficient to one. These unique
representations establish exact operand and defining-point order invariance.

At most six equalities, six inequalities and three affine unknowns are involved.
This bounds the combinatorial work, not a wall-clock production latency.
The supported binary64 range bounds input bit lengths. All eliminations,
feasibility decisions and hull orientations use exact rational arithmetic.

The independent reference constructs line intersections using cross products,
intersects triangle boundary edges with the other operand, and collects
contained endpoints. Triangle/triangle extreme points are input vertices
inside the other triangle or boundary intersections; for noncoplanar triangles,
the intersection segment endpoints occur on an input boundary. A separate
3D gift-wrapping hull establishes the complete reference result. Rust's oracle
uses Gram barycentric triangle membership; Python uses oriented halfspaces and
its own integer arithmetic. Neither reproduces production affine elimination
or its halfspace-pair enumeration. Tests compare full exact results, all
individual input permutations and operand swaps, and independently verify
minimal finite bounds. Fuzzing also checks agreement with the zero-distance
predicate on scaled-integer modes. Native observations remain supplementary;
[known nonresults and extra geometry](NATIVE_LINEAR_INTERSECTION_DIVERGENCES.md)
never weaken the mathematical checks.

## Exact Bézier extraction and editing

The [editing contract](BEZIER_EDITING.md) retains positive homogeneous rational
controls `H_i=(w_i*x_i,w_i*y_i,w_i*z_i,w_i)` and an increasing exact domain `[a,b]`.
Its local Bernstein parameter is `t=(u-a)/(b-a)`. No edit converts controls back
through binary64. Each extracted arc keeps its own one-sided endpoint jets.

For a nonempty spline span, the degree-p blossom evaluated at p-i copies of the
clipped lower parameter and i copies of the upper parameter is Bernstein control
i on that interval. Generalized de Boor stages evaluate these arguments in exact
rational arithmetic. The denominator at stage r, index j, is
`U[span-p+j+p-r+1]-U[span-p+j] > 0`; both arguments lie in the active span,
so each interpolation fraction belongs to `[0,1]`. Thus every output weight
stays positive. Shared lower-argument prefixes reduce duplicate blossom work.
Periodic spans use the original extended knot vector with exact translated
parameter bounds. Preflight counting rejects excessive turns before enumeration.

Subdivision uses the two boundary diagonals of the homogeneous de Casteljau
triangle. Its final entry is shared exactly by both results. For `t=n/d`, controls
are first lifted to integers over a common positive denominator D. Each stage
uses `(d-n)*A+n*B` and denominator `D*d^r`; only returned controls are reduced.
Evaluation uses the same integer recurrence. Elevation similarly keeps a common
denominator, multiplying it by p+1 at each increment. These are exact arithmetic
rearrangements, not rounding or projective rescaling of the published controls.
They avoid repeated fraction reductions exposed by degree-25 fuzz regressions.
Trimming composes
subdivisions while retaining original parameter units. Reversing control order
substitutes `a+b-u`. Elevation from p to p+1 uses
`H'_i = i/(p+1)*H_(i-1) + (1-i/(p+1))*H_i`, with unchanged endpoint controls.
These Bernstein identities preserve the entire homogeneous polynomial, hence
its rational geometry because the weight polynomial is positive.

Evaluation differentiates homogeneous Bernstein controls, including division by
the original domain length, then applies the rational quotient rule. Exact
points and derivatives remain available when a finite floating enclosure cannot
be represented. Position-only queries avoid unnecessary derivative conversion.
Input rationals are normalized and zero denominators rejected; interval checks
and degree limits precede edits. Degree <=25 and the caller's arc-count budget
bound combinatorial work. Arbitrarily large rational cut parameters and long
edit sequences do not yet have hard time/allocation guarantees.

Independent Python and Rust oracles expand Cox basis functions into power
polynomials, apply direct binomial affine substitutions, and compare complete
coefficients/controls. They do not repeat the production blossom or subdivision
recurrences. Rust's checker clears common denominators before linear polynomial
transforms and normalizes only their results. Exact jets, minimal enclosures,
commutation and source endpoint identities add further checks. The native
API observations are supplementary; [reviewed degree-25 differences](NATIVE_BEZIER_EDITING_DIVERGENCES.md)
never bypass these mathematical checks. General spline knot refinement/removal
is described next; attaching these curves to topology remains separate work.

## Exact rational B-spline knot editing

`ExactKnotVector` retains the same multiplicity/pole-order family over normalized
rationals. `ExactBSplineCurve3` stores positive homogeneous controls and supports
exact jets, rational-interval extraction, batch refinement and exact removal.
See [KNOT_EDITING.md](KNOT_EDITING.md) for input and failure contracts.

For one inserted knot `u`, let `k` be its last index in the old expanded vector
(or the preceding knot index for a new value), `s` its old multiplicity and `p`
the degree. For `i=k-p+1..k-s`, the new homogeneous controls satisfy

```text
alpha_i = (u-U_i)/(U_(i+p)-U_i)
Q_i = (1-alpha_i)*P_(i-1) + alpha_i*P_i.
```

Earlier controls are copied and later controls shifted one index. The closed
valid domain makes each interpolation coefficient lie in `[0,1]`; exact
insertion preserves positive weights. This includes controls inactive inside
an unclamped fundamental domain. They must not be discarded simply because
point samples cannot observe them.

Removal constructs the coarser knot vector and solves these equations from the
unaffected left anchor. For a strictly interior finite removal, every divisor
`alpha_i` is strictly positive. The recovered right anchor must equal its
unaffected neighbor in all four homogeneous components. Exact equality proves
inverse insertion; disagreement returns no curve. Final nonpositive weights
also return no curve. This is a homogeneous representation contract, stronger
than arbitrary equivalence of Euclidean rational functions or tolerance-based
simplification.

Periodic editing unrolls the cyclic knot/control functions into five periods.
Insertions update the relevant translated copies, including both outer domain
ends for a seam; removals act on four strictly interior translated copies.
Because supported curves have more poles than degree, the desired central
extended knot vector lies within the edited support. Its exact knot sequence
selects the canonical control block. Removing the last seam occurrence advances
the origin to the next distinct knot and retains the period. Parameter-to-point
correspondence is preserved, even where native behavior differs.

The independent checker expands Cox basis functions on every common span,
including the entire finite raw support for unclamped curves and one complete
period for periodic curves. A coefficient linear system recovers the complete
new controls or proves inconsistency, independently of inverse insertion.
The Rust checker clears denominators and uses fraction-free integer equations,
asserting exact divisions. Once the controls are uniquely determined, it checks
every remaining coefficient equation by exact substitution;
Python uses `Fraction` elimination, with a separate Greville reconstruction on
selected cases. No sampled-point comparison determines removal acceptance.
The 723 fixtures, 667 native observations and thirteenth sustained fuzz target
cover this scope. General spline elevation and surface knot editing remain
separate work. Arbitrary rational bit lengths still have no wall-clock bound.

## Exact tensor Bézier extraction and editing

`ExactBezierSurface3` retains a U-major homogeneous grid and two increasing
rational domains. Extraction raises each clipped boundary's multiplicity to
the degree (already clamped boundaries retain degree+1) by exact local knot
insertion in U and V. The active block is the Bernstein representation of that
rectangle. A reusable map on unit controls preserves the entire homogeneous
tensor polynomial; integer dot products apply it to grid rows after clearing
common denominators. This computes the same map as spline blossoming with fewer
interpolation stages. An unclamped local endpoint uses the last equal knot
index, not a span convention that assumes clamped end multiplicities.
Traversal checks both axis counts and their Cartesian product
before constructing geometry; the default limits are 4096 patches and
1,048,576 output controls. They bound combinatorial output, not rational bit
growth or wall-clock time.

Row-wise homogeneous de Casteljau, reversal and degree elevation provide exact
U/V subdivision, rectangular restriction and degree elevation. Transposition
exchanges controls, degrees and parameter domains. Constant-U/V substitution
returns an `ExactBezierCurve3` in the other original parameter's units, including
exact boundary curves. Every weight stays positive by convex combination.
Operations retain raw homogeneous controls without projective rescaling.

Homogeneous differentiation in both axes and the multivariate quotient rule
retain exact `(00,10,01,20,02,11)` jets. Binary64 enclosure is separately requested
and may fail without losing these values. Each patch owns its one-sided boundary
limits; neighboring patches need a separate continuity decision. Rational cut
parameters are normalized and checked before work. Arbitrary trim loops,
nonrational algebraic cuts, sewing and B-rep integration are not implied.

Independent Python and Rust Cox tensor power coefficients, affine binomial
substitutions and closed quotient formulas check complete controls/coefficients
and jets. Shared boundary curves, axis commutation, sub-ULP rational cuts,
derivative overflow and preflight exhaustion have additional direct tests.
See [the contract](SURFACE_EDITING.md) and the
[native observations](NATIVE_SURFACE_EDITING_DIVERGENCES.md).

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
- `fixtures/linear_sets.tsv`: 684 complete rational intersections, all 433 native
  inputs, full-exponent/subnormal and overflow cases, empty/point/segment/line/plane
  results, coplanar polygons with three through six vertices, and typed invalid
  definitions. Independent boundary algorithms and all defining-point orders
  are checked; the `linear_sets` target adds retained mutation campaigns.
- `fixtures/proximity.tsv`: 554 exact cases across all 25 ordered pairings,
  including all 370 native inputs, subnormal/full-exponent coordinates, singular
  minima, invalid geometry and separately unrepresentable outputs. Every valid
  result and its operand-swapped query pass the global optimality certificate;
  repeat queries must return the same witness. The dedicated fuzz target adds
  arbitrary binary64 mutations, exact coordinate permutations and vertex reversal.
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
- `fixtures/splines.tsv`: 1,059 independent exact basis-function/quotient cases
  for positions and first/second derivatives. Includes degree 25, rational and
  polynomial curves, unclamped domains, repeated knots, derivative-side
  decisions, extreme weights/coordinates and subnormal knot spans. Separate
  tests cover power-of-two weight invariance, coordinate permutation, original
  parameter units, malformed input and position-only queries despite derivative
  overflow, periodic seams and exact wrapping of extreme parameters. The live native corpus has 456 curve observations.
  One OCCT 7.6.3 second-derivative result exceeds the native comparison budget;
  the exact oracle confirms Rust and the version/input/value-pinned difference
  is reported separately, never counted as a match. See `VALIDATION.md`.

- `fixtures/surfaces.tsv`: 791 exact tensor-basis/closed-quotient cases including
  all U/V periodicity combinations, mixed derivatives, intersecting knots,
  degree 25 in both directions, unclamped domains and extreme exponents.
  Separate invariants transpose the U/V grid and rescale common weights.
  The native corpus has 210 surface observations. High-degree native evaluation
  discrepancies are independently verified and [reviewed separately](NATIVE_SPLINE_DIVERGENCES.md).

- `fixtures/bezier-editing.tsv`: 636 complete exact extraction/editing results,
  including all 546 native inputs. Independent homogeneous power coefficients
  determine every expected control. Additional cases cover full exponents,
  subnormal and adjacent knots, periodic multi-turn ranges and shifted sub-ULP
  arcs. Separate tests use rational cuts 2^-2048 apart, exact edit commutation,
  endpoint derivatives, malformed parameters and preflight traversal limits.

- `fixtures/real-roots.tsv`: 95 exact cases including degree 25, negative leading
  coefficients, multiplicities, clipping, algebraic sign queries, unrepresentable
  roots and distinct clustered roots with identical binary64 enclosures. SymPy
  irreducible factors and Vincent–Akritas–Strzebonski continued fractions provide
  independent isolation; rational interval evaluation decides reduced query signs.
- `fixtures/spline-plane.tsv`: 284 complete results from independent exact
  Cox–de Boor basis polynomials and that continued-fraction oracle. Includes all
  214 native inputs, random rational spans, periodic/unclamped domains, unequal
  contact orders, overlap chains, extreme weights and subnormal spans. All
  preexisting 1,850 spline jet fixture rows remain unchanged. The original 185
  intersection rows are also unchanged; 98 appended cases cover clipped contacts,
  multi-period queries, nonrepresentable shifted knots, subnormal periods,
  degree-25 clipping and closely clustered roots beside a trim boundary.
- `fixtures/spline-quadric.tsv`: 118 complete sphere/cylinder results, including
  all 92 native inputs, dense rational degree-25 curves, order-50 tangencies,
  contained rational arcs and generators, periodic/trimmed/unclamped ranges,
  subnormal radii/axes, overflowing squared magnitudes and coincident parameter
  enclosures. The independent basis/continued-fraction oracle uses the cylinder
  cross-product equation, while production uses dot-product projection.
- `fixtures/knot-editing.tsv`: 723 complete representations and operation flags,
  including 667 native inputs, exact rejection, rational weights, extreme
  domains, inactive unclamped controls and periodic seam origin changes. Every
  output is independently recovered from complete Cox coefficient equations.

These fixtures are bounded deterministic tests. Separately, [coverage-guided
fuzzing](FUZZING.md) mutates predicates, linear/curved intersections, polynomial
roots, spline evaluation/intersections and standalone modeling
sequences under sanitizers, with retained corpora and daily campaigns.
Finite test evidence supplements the arithmetic argument; it does not prove a
general kernel correct. The OCCT comparison corpus remains a separate oracle.

```sh
cargo test --workspace --locked
cargo test --workspace --locked --release
python3 rust/tools/generate_predicate_fixtures.py --check
python3 rust/tools/generate_spatial_fixtures.py --check
python3 rust/tools/generate_curved_fixtures.py --check
python3 rust/tools/generate_spline_fixtures.py --check
python3 rust/tools/generate_surface_fixtures.py --check
```

Native CI executes both debug and optimized Rust tests. The rational fixture
generator is checked independently in the Linux validation job. Its default
write mode is for deliberate, reviewed corpus additions, never for making a
disagreement disappear.

## Intersections after exact curve editing

Plane, sphere and infinite-cylinder queries also accept `ExactBSplineCurve3`
directly. Exact knot refinement and removal therefore compose with certified
intersection without rounding their rational knots, homogeneous controls or
trim interval. The shared engine forms each homogeneous span polynomial in a
local rational parameter, clears its positive weight denominator and isolates
every real root of the resulting implicit equation. Identically zero spans
produce maximal closed overlaps. Exact shared-knot identity merges duplicate
contacts while preserving left/right contact orders; roots with the same
floating enclosure remain separate.

`ExactSplineSurfacePoint` retains that algebraic identity and its original
rational parameter transformation. Comparisons against rational parameters and
coordinates require no binary64 result. Coordinate comparisons first use the
exact positive-weight control hull, then the sign of `H_c - x H_w` at the root
when the hull does not decide the answer. This remains valid when individual
controls exceed binary64 range but cancel to a finite contact.

All four homogeneous polynomials share one positive integer scale, cleared once
per retained span. A query for rational `x=n/d` therefore signs `d H_c - n H_w`
using integer coefficients, without repeated rational normalization. Before a
nontrivial algebraic sign query, a polynomial whose degree reaches that of the
square-free defining polynomial is reduced by positive pseudo-division. At a
root of the defining polynomial, this remainder has the original query's sign,
including exact zero. Interval filters and the complete Sturm-Tarski fallback
then operate on that lower-degree remainder. Neither optimization changes the
mathematical decision or substitutes sampling for a complete check.

The de Boor coefficient interpolation itself also keeps one positive
denominator per vector of polynomials. At each stage it clears the two input
denominators and the two affine mixing coefficients once, performs the same
linear-polynomial interpolation with integers, then cancels a common divisor
of the denominator and all coefficients. Final conversion to normalized
rationals recovers exactly the original homogeneous polynomials. Independent
Cox coefficient fixtures and all four curve representations check this change.

Bounded root refinement also tries rational candidates obtained from common
continued-fraction prefixes of the isolating interval. Subtracting a shared
integer part and reciprocating its positive fractional interval reverses the
bounds; undoing those maps produces a candidate in the original interval.
The kernel independently verifies both interval membership and an exact zero
of the defining polynomial before storing a rational root. A failed candidate
leaves the complete algebraic path unchanged. Candidates are tried initially
and every 16 bisections, with at most 32 eager bisections for planes and 128 for
quadrics. This avoids repeating an algebraic sign solve for contacts already
provably rational; it is not a tolerance or a rationality assumption.

Minimal floating enclosures are optional, fallible conversions. A parameter
outside finite binary64 range can still have a finite position, and a finite
parameter can still have an unrepresentable coordinate. An enclosure failure
does not erase any exact contact or overlap; converting the complete result
succeeds for every requested component or returns an error. Rational overlap
endpoints remain directly available. Span and root-subdivision limits still
bound counts, not arbitrary-precision operand size or wall-clock time.

The 274 additional Python `Fraction` fixtures independently expand known-factor
equations and check complete contacts, orders, coordinates, overlaps and minimal
bounds. The existing 402 independently certified intersection cases also run
through exact conversion, rational refinement and exact removal. Native
comparisons reuse their 306 original geometries across those representations;
this does not multiply the distinct native case count. The fourteenth fuzz
target adds mutated exact edits, complete coefficient identities, close roots,
periodic seams and extreme rational domains. See
[the exact intersection contract](EXACT_SPLINE_INTERSECTIONS.md).

## Next mathematical work

1. Add incircle comparisons and extend certified linear-set distances to the
   curved geometry required by subsequent algorithms.
   Preserve explicit units, domains and arithmetic bounds; do not reuse an
   arbitrary global epsilon.
2. Complete general spline elevation acceptance, add arbitrary face trimming and projections, and extend root isolation
   to general intersections. Use interval/error bounds or
   additional precision when ordinary arithmetic cannot establish the answer.
3. Strengthen topology invariants and tolerance propagation through each
   operation; sampled checks alone cannot certify full curve/surface agreement.
4. Extend the existing coverage-guided campaigns to each new capability,
   minimize failures, and expand resource/performance and independent-oracle
   checks as geometry and operation sequences become more complex.

General curve/surface intersections beyond spline/plane/sphere/cylinder, Booleans,
generic topology history, STEP and meshing remain unimplemented. Certified
linear/quadratic primitives and spline evaluation do not establish these capabilities or make
the rest of the kernel exact.

## Exact B-spline surface knot edits

`ExactBSplineSurface3` retains both rational knot vectors and all positive
homogeneous tensor controls. Combined refinement preflights both axes and the
Cartesian control limit. Exact removal must preserve every homogeneous
transverse row over full raw support; an inexact inverse returns `None` and
cannot partially update a surface. U/V periodic origins and their aliases are
independent. Retained isocurves can be passed to exact analytic intersectors.

The differentiated de Boor map is computed once per axis and applied through
integer dot products, with exact rational quotient rules through total order
two. Explicit knot quadrants and automatic continuity retain the existing
surface evaluation contract. See [surface knot editing](SURFACE_KNOT_EDITING.md)
for rational domains, limits, independent equation checks and exclusions.

## General exact degree elevation

Degree elevation preserves all four homogeneous spline functions on the
original active domain. Periodic multiplicities rise by the degree increment,
with the same cyclic origin and period. For a nonperiodic axis, the original
first/last active distinct-knot indices determine how many exterior flat knots
are removed; partial exterior multiplicities are retained. Inactive controls
are determined by a fully clamped zero-padded working extension and its exact
crop. Equality outside the original active interval is not implied.

Production uses Prautzsch rank averaging and shared exact refinement maps.
The independent oracles instead reconstruct every working control from every
Cox power-coefficient equation, require full rank, and verify all residuals.
They do not infer equality from samples. The Rust equation solver uses primitive
integer elimination for degree changes, skipping zero pivot columns while
retaining all constraints and checking integer divisions. The target pole count
and the full tensor product are checked before production control arithmetic.
See [contracts, native discrepancies and current gates](DEGREE_ELEVATION.md).

## Complete point-to-spline minimum sets

Closed point-to-rational-B-spline queries enumerate all isolated minimum
parameters and maximal constant-distance intervals. They retain periodic aliases
and singular stationary points. Positive homogeneous weights give squared
distance `N/W²`; its degree-at-most-73 stationary numerator `N'W-2NW'`, all knots and both
trim endpoints form a complete candidate set. The degree-at-most-25 common
polynomial of the three homogeneous spatial deltas first identifies every
zero-distance parameter. A real zero proves the local minimum without solving
the larger stationary equation; all-zero deltas give a whole minimum interval.
After the global distance reaches zero, only complete zero sets can contribute
further minimizers. Spans with no real zeros fall back to the stationary
equation until then. Shared knot hits are deduplicated by exact parameter
identity. Exact interval filters and
quotient-ring image polynomials decide global ties without approximate equality.
The image construction carries relative rational scales and full integer row
dependence witnesses. Independent Cox/VAS/resultant fixtures verify complete
minima, tight views and every coefficient of the image construction.
See [the full contract, proof structure and limits](SPLINE_PROXIMITY.md).

## Certified B-rep validation

A curve-on-surface deviation `D(t) = C(t_edge) - S(P(t))` for line/arc edges
against plane line/arc pcurves and cylinder line pcurves is a harmonic sum
`A0 + A1 t + Σ_ω (C_ω cos ωt + S_ω sin ωt)` with exact rational frequencies.
Equal frequencies are combined, then
`sup_[0,1] |D| <= sqrt(|A0|² + |A0+A1|²) + Σ_ω sqrt(λmax(Gram(C_ω, S_ω)))`.
Two terms at different frequencies `ω1 < ω2` are `Re(z1 e^{iω1 t} + z2
e^{iω2 t})` with `z = C - iS`; since `|e^{i(ω2-ω1)t} - 1| <= (ω2-ω1) t`, on
`[0, 1]` their sum is at most the ellipse bound of `z1 + z2` plus
`|z2| (ω2 - ω1)`. Terms adjacent by frequency use this pair bound whenever it
is smaller, so a circle against a pcurve whose period differs by a few ulps
(OCCT prints `2π` as `6.28318530717959`) certifies as it should.
Loop winding is the sign of the closed form of `∮ u dv - v du` over pcurves
closed by chords. Shell orientation is the sign of Green's volume integral
`Σ ∫ G dv` with `∂G/∂u = S·(S_u × S_v)`: `G = u(o·(x×y))` on planes and
`G = r[(o·(y×n)) sin u + (o·(x×n)) cos u + r·det(x,y,n)·u]` on cylinders.
Cavity containment uses exact ray/plane and ray/cylinder solves. A ray counts
only when every hit is certifiably farther than twice the tolerance from its
face's boundary. Every decision runs first in exactness-preserving binary64
intervals (TwoSum/FMA error signs) and then, if undecided, in rational
intervals; anything still undecided is reported as uncertified. See
[the validator's contract and evidence](BREP_VALIDATION.md).

Stored enclosures (M5) are the upper ends of those same enclosures: for a
vertex, `sqrt` of the squared distance to each curve end (and to a vertex
loop's surface); for a fin, the harmonic bound above; for a face, `sqrt` of
each squared UV gap with `du` scaled by the radius. Each is stored as the
next binary64 value above the interval's upper end, and at least `2^-80`. A
later check `d² ≤ b²` (or `sup |D| ≤ b`) is then strict in the tier that
produced it: `b² − hi(d²) ≥ 2 b ulp(b) > 0` exceeds the width of the
binary64 enclosure of `b²`, and `(2^-80)² = 2^-160` lies far above the
rational tier's `2^-192` grid. Tolerance decisions outside the validator
(T4) take binary64 inputs as exact rationals. Distances compare squared,
`|p − q|² ≤ (Σ t_i)²` for a threshold `Σ t_i ≥ 0` that is itself an exact
sum (for example `r + tol` or `max(a, b) − min(a, b) − tol`). A
point–segment distance uses the exact foot parameter `t = w·d`, with
`|w|²`, `|p − b|²` or `(w × d)² / |d|²`. These decide first in binary64
intervals, then exactly. The material-area screen `A ≤ tol · P / 2`
encloses `π` and the edge-length square roots in rational intervals, and
counts an undecidable case as degenerate.

## Complete spline/line and spline/segment preimages

Closed rational B-spline ranges against infinite lines or closed segments
return every isolated curve parameter and every maximal closed parameter
interval. With positive span weight `W`, `Δ=(X,Y,Z)-AW` and `D=B-A`, a
parameter lies on the line exactly when `Δ×D=0`. The primitive gcd of those
three components has degree at most the curve degree. A collapsed segment uses
`Δ=0` directly. Segment clipping is `0 <= Δ·D <= (D·D)W`, without division.
Contained spans are clipped by a complete sign decomposition. Every boundary
root is isolated, and each open cell's sign is the first nonzero derivative
sign at its left boundary, so no approximate midpoint decides topology.
Results from different spans are ordered by exact comparison of positive
affine root images, with no image resultant. Periodic aliases and backtracking
preimages remain distinct. Independent sample-point clipping fixtures and a
constructed-answer fuzz target verify complete sets. See
[the contract, native bridge and acceptance evidence](SPLINE_LINEAR_INTERSECTIONS.md).

This work also changed the shared root arithmetic. Polynomial content uses a
binary integer gcd. Rational-root recognition adds the rational root theorem
candidate `ceil(lead·lower)/lead`, which is unique once the isolator is
narrower than `1/|lead|`, and still requires an exact zero check. Both change
cost only, not results.

## Exact support of split and merged pieces

The history checker (`history::check`) decides whether a split child or a
merged parent lies on the support of the whole with exact rational arithmetic
on the binary64 inputs, within the body's resolution `t`:

* a point `p` is within `t` of the line through `a` and `b` when
  `|(b - a) x (p - a)|^2 <= t^2 |b - a|^2`;
* a point `q` is within `t` of the plane through `o` with normal `n` when
  `((q - o) . n)^2 <= t^2 |n|^2`, and of the axis through `o` with direction
  `n` when `|(q - o) x n|^2 <= t^2 |n|^2`;
* a plane piece (origin `o'`, normal `n'`) must face the whole's normal
  (`n' . n > 0`), and every end `p` of its line boundary edges, projected onto
  the piece's own plane as `q = p - ((p - o') . n' / n' . n') n'`, must lie
  within `t` of the whole plane. The signed distance to the whole plane is
  affine on the piece's plane, so its maximum over the face is attained on the
  projected boundary; a circular boundary edge must have a normal exactly
  parallel to both planes (`m x n' = m x n = 0`) and its projected centre
  within `t`;
* a cylinder piece must have an axis exactly parallel to the whole's, with
  `s = t - |r' - r| >= 0` and a point of its axis within `s` of the whole's
  axis (`|(p - o) x n|^2 <= s^2 |n|^2`): parallel axes at distance `d` and
  radii differing by `e` keep the surfaces within `d + e <= t` of each other
  everywhere. The builder's pieces share the whole's axis direction bit for
  bit.

Nothing else counts as the same support, so a split along a slightly rotated
plane or a merge of non-coplanar walls is reported, never absorbed.
