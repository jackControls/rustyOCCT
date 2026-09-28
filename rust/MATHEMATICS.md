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

General curve/surface intersections beyond spline/plane/sphere/cylinder and the
lines and circles of S7c.1, Booleans,
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

A curve-on-surface deviation `D(t) = C(t_edge) - S(P(t))` for line, arc and
ellipse edges against plane line, arc and axis-aligned ellipse pcurves,
cylinder line pcurves and cylinder sinusoids `v = a0 + a1 cos u + a2 sin u`
(whose image `O + a0 n + cos u (R x + a1 n) + sin u (R y + a2 n)` rotates at
the pcurve's own frequency, S8a.2) is a harmonic sum
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

## Cones and general mass properties (S3)

A cone face is `S(u, v) = O + ρ(v)(cos u x + sin u y) + v cos a n` with
`ρ(v) = R + v sin a`, as OCCT's `Geom_ConicalSurface`: `v` is arc length
along a ruling and the apex is at `v_a = -R / sin a`. Its pcurves live on the
universal cover in `(u, v)`, like a cylinder's. A UV gap scales `du` by
`|ρ(v)|` at the larger end, so every pcurve meets the apex at any `u`. A
point's distance to the (double) cone is the smaller of `|r cos a - R cos a -
z sin a|` and `|r cos a + R cos a + z sin a|` in the meridian half-plane,
enclosed in the same two tiers. A line pcurve's deviation is harmonic along a
ruling (`du = 0`, affine in `t`) and a parallel (`dv = 0`, one frequency);
anything else is not, and certifies only as a failure, like an arc pcurve on
a cylinder.

A plane's section of a cylinder (S8a.2) integrates in closed form. Along an
axis-aligned ellipse pcurve `u = cx + a cos t`, `v = cy + b sin t` the
twice-area term is `b (cx - ou)(sin t1 - sin t0) - a (cy - ov)(cos t1 -
cos t0) + a b (t1 - t0)`, the plane flux `a (cy - v0)(cos t0 - cos t1) + a b
((t1 - t0)/2 - (sin 2t1 - sin 2t0)/4)`, and a moment polynomial in `(u, v)`
a polynomial in `(cos t, sin t)` integrated by the Fourier expansions of the
arc's. Along a sinusoid `v = a0 + a1 cos u + a2 sin u` on a cylinder,
`-∫ v du = -(a0 Δu + a1 Δsin u - a2 Δcos u)`, `∫ u dv = Δ[a1 (u cos u - sin
u) + a2 (u sin u + cos u)]`, the twice-area term `Δ[a1 w cos u + a2 w sin u
- 2 a1 sin u + 2 a2 cos u] - (a0 - ov) Δu` with `w = u - ou`, the flux
`-∫ v (A sin u + B cos u + C) du` from the integrals of `1, sin, cos, sin
cos, cos², sin²`, and each moment monomial `v^k cos^p u sin^q u` a
polynomial in `(cos u, sin u)` integrated exactly. Containment parity along
an ellipse is its unit circle's after scaling each axis about the centre,
and a point is clear of an ellipse by `m` when its scaled distance from the
unit circle exceeds `m / b` (the scaling shrinks distances by at most `b`).

The apex is a pole: a vertex loop of the wall face. On a cone face whose edge
loops wind once in total (`|Σ w| = 1`), the first vertex loop closes the band
at the apex, so the windings need not balance. Its vertex must lie on the
surface and within tolerance of the apex (`pole_off_apex`). The periodic
area `-∮ v du` adds the pole's term `2π · W · v_a`, with `W` the winding the
pole closes and `2π` from the certified `π`. Orientation flux on a cone is
`-∮ F du` with `∂F/∂v = ρ(v) h(u)` and `h(u) = cR + c(O·x cos u + O·y sin
u) - s O·n` (with `c = cos a`, `s = sin a` and `O` the origin relative to
the reference). `F` is taken from the apex on a face with a pole, where it
vanishes (`F = ρ(v)² h(u) / (2s)`), and from `v = 0` elsewhere (`F = (R v + s
v² / 2) h(u)`). A nearly cylindrical cone's apex is far away (at `s = 2.4e-4`
it is `4096 R` along the ruling), and starting there would leave two large
terms of a band's loops to cancel. Along a line pcurve `F` is a quadratic in
the line's parameter times `cos` and `sin` of `u0 + w`, integrated in closed
form by the moments `∫ w^k cos(u0 + w) dw` and `∫ w^k sin(u0 + w) dw`, `k ≤
2`. Cavity containment on a cone face is not yet decided
(`uncertified_containment`).

General mass properties (REVIEW_NOTES.md U2) apply the divergence theorem to
every solid region's boundary:

* the volume is `(1/3) ∫ p·N`, the first moments `(1/2) ∫ p_i² N_i`, the
  second moments `(1/3) ∫ p_i³ N_i`, and the mixed ones `(1/2) ∫ p_i² p_j N_i`
  (for `xy`, `yz` and `zx`), all with `N = S_u × S_v` and `p` relative to a
  reference point of the body, so that a far body does not cancel;
* the area is `∫ |N|` and a face's centre `∫ p |N| / ∫ |N|`. `|N|` is 1 on a
  plane, `r` on a cylinder and `|ρ(v)|` on a cone, a polynomial in `v` of
  fixed sign on a face that does not cross the apex;
* on a plane every integrand is a polynomial in `(u, v)`; on a cylinder or
  cone a polynomial in `v` times `cos^a u sin^b u`. Green's theorem turns each
  face integral into `-∮ F du` over the face's loops, closed by their chords,
  with `F` the integrand's antiderivative in `v`;
* along a line pcurve this is a polynomial in the line's parameter times
  `cos` and `sin` of integer multiples of `u`. The exact Fourier expansion of
  `cos^a sin^b` (dyadic coefficients with numerators below `2^(a+b)`, held
  exactly in binary64 and memoized) reduces it to the moments `∫ w^j
  cos(f(u0 + w))`, which obey an exact recursion; all fourteen integrands of
  a line share its trigonometric values and moments. Along a plane arc the
  integrand is a trigonometric polynomial, integrated through the same
  expansion.

Every step runs in the validator's two tiers, binary64 intervals first and
then rational intervals, so the results are enclosures with certified error
bounds, not quadrature estimates. `Topology::mass_enclosure` returns them and
`Solid::mass_properties` reports a cone's midpoints. On all 546 identity
prisms the engine agrees with the prism closed forms to `2e-15` relative.
The 22 cones of `primitive-cases.txt` agree with the independent reference
(mpmath disc integrals) to `1e-12` of their scale, the bound for the rounded
frame, angle and slant the kernel stores (`tests/cones.rs`).

### Spheres

A sphere face is `S(u, v) = O + R (cos v e(u) + sin v n)` with
`e(u) = cos u x + sin u y`; `S_u × S_v = R² cos v q` with `q` the outward
radial direction, so `|N| = R² cos v` and the u-scale of a UV gap is
`|R cos v|`, zero at the poles `v = ±π/2`. A pole closing a band is on the
band's material side: the material lies left of the traversal, so a band
winding `+u` on a forward face closes at the north pole. Along a meridian
the deviation is one rotating term of frequency `dv` in the plane of `e(u0)`
and `n`; along a parallel, one of frequency `du` in the plane of `x` and
`y`.

The mass integrands (and the orientation flux `S·N`) are sums of
`cos^a u sin^b u cos^c v sin^d v`. `F(u, v)`, the antiderivative in `v`
from the pole on a face with one and from the south pole otherwise, is
never expanded: along a parallel `v = v0`, `F(u, v0)` is a trigonometric
polynomial in `u` whose coefficients are the exact integrals of
`cos^c sin^d` over `[lower, v0]`, integrated in `u` by the same Fourier
expansion; along a meridian `-∮ F du` is zero; along any other segment (a
chord closing a gap), its value lies in `-du` times `F` over the segment's
bounding box, since the mean value of `F` on the segment does. A whole
sphere has no loops: on the universal cover its region is bounded by the
two pole lines, of which only the north one contributes, `∫ F(u, π/2) du`
over one turn, with the face's sense. The independent reference integrates
its own closed form of the flux, `G(u, v) = R²[A(u)((v + π/2)/2 + sin 2v/4)
+ (O·n)(sin² v - 1)/2 + R (sin v + 1)]`, by quadrature.

### Tori

A torus face is `S(u, v) = O + (R + r cos v) e(u) + r sin v n`, a ring
torus (`R > r`); `S_u × S_v = r (R + r cos v) q` with `q = cos v e(u) +
sin v n` the tube's outward normal, so the UV-gap scales are `R + r cos v`
in `u` and `r` in `v`. Both parameters are periodic: a loop winding
`[wu, wv]` closes shifted by `(2πwu, 2πwv)`, and the windings in `v`
balance. Along a meridian (`du = 0`) the deviation is one rotating term of
frequency `dv` in the plane of `e(u0)` and `n` about the tube's centre;
along a parallel, one of frequency `du`.

The integrands are sums of `cos^a u sin^b u cos^c v sin^d v`, as on the
sphere. Loops winding in `u`, or none, take `-∮ F du` with `F` the
antiderivative in `v` from 0: exact along parallels, zero along meridians,
enclosed along chords. Loops winding in `v` (a wedge's wall) take `∮ G dv`
with `G` the antiderivative in `u` from 0, the same routine with the
parameters exchanged: exact along meridians, zero along parallels. The
whole torus has no loops; on the cover its region is the period square, and
by Green's theorem its integral is `∫ F(u, 2π) du` over one turn, with the
face's sense. The orientation flux is `r (R + r cos v)(A(u) cos v + (O·n)
sin v + R cos v + r)` with `A(u) = O·x cos u + O·y sin u`; over the whole
torus `A` integrates to nothing and the flux is `6π²Rr²`, three times the
volume `2π²Rr²`. The independent reference integrates closed forms of both
antiderivatives by quadrature (`cell_reference.py`), and
`primitive_reference.py` the partial revolution's moments in closed form
(Pappus's theorems for the tube's disc or the meridian region).

The 22 tori of `torus-cases.txt` agree with the reference to `1e-12` of
their scale (`tests/tori.rs`).

## Spline uses, areas and fluxes (S4b-d)

**Deviation by exact composition.** Where every factor of a use is
rational, `D(t) = C(s(t)) - S(P(t))` is a rational function on every common
piece: between the edge's knots (in the use's fraction, mirrored for a
reversed use), the pcurve's knots and a line pcurve's crossings of a spline
surface's knot lines. Each factor is extracted exactly as a Bézier piece
(homogeneous Bernstein coordinates on `[0, 1]`). On a plane
`S(u, v) = O + u X + v Y` with the stored axes, and `Ŝ = (O P_w + X P_u +
Y P_v, P_w)`; on a spline surface's patch of degrees `(du, dv)` over
`[u0, u1] × [v0, v1]`, with `U = (P_u - u0 P_w)/(u1 - u0)` and `V` likewise,
`Ŝ = Σ Q_ij C(du, i) C(dv, j) U^i (P_w - U)^(du-i) V^j (P_w - V)^(dv-j)`,
which is `S_h(P) P_w^(du+dv)`, a polynomial. Then `N = C_xyz Ŝ_w - Ŝ_xyz C_w`
and `W = C_w Ŝ_w` are polynomials, `D = N/W`, and on `[0, 1]`
`|N| <= max_i |N_i|` (a convex combination of the Bernstein coefficients)
and `W >= min_j W_j`. When every `W_j > 0` and `|N_i|^2 <= tol^2 W_j^2` for
all `i`, `j`, the piece is within tolerance; otherwise it is halved (de
Casteljau) up to ten times. Cancellation is exact because the pieces are
exact; the products run in binary64 intervals first and exactly when those
cannot decide. A sample `|N(t)|^2 > tol^2 W(t)^2` certifies a failure. The
measured enclosure of such a use is `max |N_i| / min W_j` over its pieces.

**Signed area.** On a piece `u = U/W`, `v = V/W`,
`u dv - v du = (U V' - V U')/W^2 dt`: the weights' derivatives cancel. With
equal weights the integral is exact (the mean of the Bernstein
coefficients); otherwise, on each of 64 subpieces with `1/W^2` in `[a, b]`
(from the coefficients of `W`), `∫ M/W^2` lies in `a ∫M ± (b - a) max |M_i|`.
The same enclosure gives `-∮ v du = -∫ V (U' W - U W')/W^3`, a plane's flux
with spline pcurves.

**Flux of a nonrational spline surface.** On a patch, in local coordinates,
`f = X·(X_ū × X_v̄)` is a tensor Bernstein polynomial (the common weight
cancels to `1/w^3`, a constant absorbed by using the poles divided by it).
Its antiderivative from `v̄ = 0` is `H_k = Σ_(j<k) f_j / (n + 1)` along each
row; with the full columns below the patch,
`G = Σ_below H(ū, 1) + H(ū, v̄)` is the antiderivative in `v` of the global
integrand from the domain's start, in local units of `u`, and by Green's
theorem the face's flux is `-∮ G dū` along its loops, each piece in one
patch, a polynomial (or, for a rational pcurve, a quotient by a power of
`P_w`) integrated as above. The reference integrates the same flux by
nested Gauss-Legendre quadrature of the exact surface derivatives.

**Taylor enclosures.** Every other spline use is bounded piece by piece,
in the piece's own parameter `τ` in `[0, 1]`:
`|D(τ)| <= |D(½)| + |D'(½)|/2 + sup|D''|/8`, with
`D' = C' - S_u P_u' - S_v P_v'` and
`D'' = C'' - (S_uu P_u'^2 + 2 S_uv P_u' P_v' + S_vv P_v'^2 + S_u P_u'' + S_v P_v'')`.
The terms at `½` are tight enclosures; `D''` is enclosed over the piece:
by the certified trigonometry for lines, arcs and the analytic surfaces,
and for splines by the hulls of exact Bézier derivative nets through the
quotient rule (a spline surface over the pcurve's box, the hull over the
patches it meets: the surface is C1 across a knot line (R4) and piecewise
C2, so the remainder bound holds). A piece is halved until the bound meets
the tolerance (fourteen times in binary64 intervals, six exactly); a
midpoint beyond it certifies a failure. For a valid use `D(½)` and `D'(½)`
are rounding-small and the remainder shrinks as `h^3`. The reference bounds
the same form on bisected pieces with `mpmath` interval de Boor
evaluation.

**Containment parity.** The `+u` ray from a point crosses a spline pcurve
an odd number of times exactly when a certain partition says so: a piece
whose control points lie certainly right of the point contributes whether
its ends lie on different sides of `v = p_v` (half-open, as for segments),
whatever its shape; a piece certainly left, above or below contributes
nothing; any other piece is halved exactly, up to twelve times, and a
pcurve through the point stays undecided. The reference integrates the
winding angle instead.

**Green integrals along spline pcurves.** On a cylinder, cone, sphere or
torus the flux and mass integrands have exact antiderivatives `F(u, v)` in
`v` (trigonometric in `u`, polynomial or trigonometric in `v`), but along a
rational pcurve `F(P(τ))` is not integrable in closed form. Each Bézier
piece of the pcurve is halved six times; on a piece, `F` is enclosed over
the box of its control points (their convex hull contains the piece) and
`P_u' = (U' W - U W')/W^2` over the hull of its Bernstein coefficients
divided by the square of the weights' lower bound, so `-∫ F du` lies in
their product: width `O(h^2)` per piece, `O(h)` in total. A gap-closing
chord is enclosed the same way.

**Flux of any nonperiodic spline surface.** When the exact route does not
apply (a rational surface, a pcurve piece across a patch boundary), the
antiderivative `G(u, v) = ∫_(v_0)^v f(u, s) ds` of `f = S·(S_u × S_v)` is
enclosed over a box `U × V` from 32 equal strips of the `v` domain: the
strips wholly below `V` contribute `Δs f(U × strip)`, the strip that `V`
reaches contributes `[0, v_hi - s_k] f(U × [s_k, v_hi])`, and a box below
the domain's start (a pcurve on the boundary patch's extension) the
negated part down to it. `f` is enclosed from the surface jets of the
Taylor enclosures (the quotient rule over exact Bézier nets). Along a line
pcurve, 32 pieces each give `-Δu G(box)`; along a spline pcurve the Green
integral above. Widths are `O(Δs + h)`: enough for the sign of a shell's
flux, which is all the orientation check needs.

**Scale, not position.** An enclosure of a rational quotient over a box,
`X(box)/w(box)`, has a width proportional to `|X/w|`, not to the body's
size, and the same holds for `u dv - v du` far from the UV origin. So every
shell's flux is `∮ (S - r)·N` with `r` a vertex of the shell (the
divergence of `S - r` is 3 whatever `r`), the strips lift the patches
translated to `r` exactly (`X - r w` on the homogeneous poles), a loop's
signed area is taken about a point of the loop, and a plane's flux and
mass integrals run in UV about a point of the loop or face (on a plane
`∮ du = 0` around a closed loop, and the mass integrands change variables
exactly). A body a few units from the origin but `10^-3` across is then
decided as at the origin.

**Mass properties with spline geometry.** The fourteen integrands of the
general mass properties (volume, first and second moments, area and face
centre, U2) follow the flux's two routes. On a nonrational patch the
position `p = X - r` (the poles minus the reference, a Bernstein partition
of unity) and `N = X_ū × X_v̄` are tensor polynomials, so the ten volume
and moment integrands (`p·N/3`, `p_i^2 N_i/2`, `p_i^3 N_i/3`,
`p_i^2 p_j N_i/2`) are too, and each is integrated exactly by the column
antiderivatives above. `|N|` is not polynomial: the area and face-centre
terms, and every term on a rational surface, use the strips, with the
integrands enclosed from the jets (`|N|` by the interval square root).
Along a spline pcurve on a cylinder or cone each term's antiderivative
`F` (a polynomial in `v` times `cos^a u sin^b u`) is enclosed over the
pieces' boxes by the Green integral. The reference integrates the same
properties by Green's theorem with nested Gauss–Legendre quadrature of the
exact surface jets, broken at knots. Since F8 every one of these
first-order routes is only the fallback of the certified quadrature
below.

**UV gaps on spline surfaces.** A spline surface has no length scale, so
the gap between consecutive pcurves' ends is measured in 3D: both ends
exact (a line's end, a spline's clamped end), each surface point evaluated
exactly by de Casteljau on the patch whose domain contains the end clamped
into the surface's domain. An end at most `1e-9` of the domain's width
outside it lies on that patch's polynomial extension, as OCCT evaluates a
file's fifteen-digit pcurve ends; farther out the gap is uncertified.

## Certified Gauss–Legendre quadrature of spline mass (F8)

The strips and Green integrals above are first order. Mass properties
(not validation, which needs signs only) now integrate every non-polynomial
term by a Gauss–Legendre rule whose remainder is bounded rigorously; the
first-order routes stay as the fallback.

**The rule and its remainder.** For `f ∈ C^(2n)[a, b]` and the `n`-point
rule with nodes `x_i` and weights `w_i` on `[-1, 1]`,
`∫_a^b f = (b - a)/2 Σ w_i f(m + (b - a) x_i / 2) + E`, `m` the midpoint,
with `E = (b - a)^(2n+1) (n!)^4 / ((2n + 1) ((2n)!)^3) f^(2n)(ξ)` for some
`ξ` in `[a, b]`. In the Taylor coefficient `f_(2n) = f^(2n)/(2n)!`,
`E = K_n (b - a)^(2n+1) f_(2n)(ξ)` with `K_n = (n!)^4 / ((2n + 1)
((2n)!)^2)` (`K_1 = 1/24`, the midpoint rule; `K_2 = 1/180`). Hence
`∫ f ∈ (b - a)/2 Σ W_i f(X_i) + K_n (b - a)^(2n+1) [f_(2n)](X)` for
enclosures `X_i ∋ x_i`, `W_i ∋ w_i` and `[f_(2n)](X)` of the coefficient
over the whole piece. As `(n!)^2/(2n)! ≈ √(πn)/4^n`, `K_n (2h)^(2n+1) ≈ π h
(h/2)^(2n)` on a piece of half width `h`: for an integrand analytic within
`ρ` of the piece, whose coefficients are about `M/ρ^(2n)`, the remainder is
about `π h M (h/(2ρ))^(2n)`. The kernel uses `n = 8`.

On a box `[a, b] × [c, d]` the tensor rule `Q = Q_x Q_y` satisfies
`∫∫ f - Q f = (I_x - Q_x)[∫ f dy] + Q_x[(I_y - Q_y) f]`. The first term is
the one-dimensional remainder of `F(x) = ∫_c^d f(x, y) dy`, whose
coefficient `F_(2n)(ξ) = ∫_c^d f_(2n,0)(ξ, y) dy` lies in `(d - c)
[f_(2n,0)](X × Y)`; the second is `Σ_i w_i` times the remainders along `y`
at the nodes `x_i`, each in `K_n (d - c)^(2n+1) [f_(0,2n)](X × Y)`, and the
weights are positive and sum to `b - a`. So
`∫∫ f ∈ Q f + K_n (b - a)^(2n+1) (d - c) [f_(2n,0)](X × Y) + K_n (b - a)
(d - c)^(2n+1) [f_(0,2n)](X × Y)`: only the pure coefficients along each
axis are needed, each from a univariate series with the other variable
held as an interval.

**Nodes and weights.** The nodes are the roots of the Legendre polynomial
`P_n` (`(k + 1) P_(k+1) = (2k + 1) x P_k - k P_(k-1)`). Each binary64
approximation `x̃` (Newton's iteration) is widened to `[x̃ - δ, x̃ + δ]`
until the exact rational values of `P_n` at the two ends have opposite
signs; the `n` brackets are disjoint, so each holds exactly one of the `n`
roots. Each is then bisected exactly, by the sign of `P_n` at its
midpoint, down to `2^-100`, and the weight `w = 2 / ((1 - x^2) P_n'(x)^2)`
(with `(1 - x^2) P_n' = n (P_(n-1) - x P_n)`) is enclosed over it in
rational interval arithmetic. Lifted into binary64 intervals the nodes and
weights are then adjacent values; brackets of `2^-52` had left node sums
`1e-13` wide. The rational enclosures are computed once per process.

**Coefficients by interval Taylor arithmetic.** A truncated series
`a = Σ_(k<L) a_k ε^k` stands for `a(x + ε)`. With `a_0`, `b_0` enclosures
and `k >= 1`: `(ab)_k = Σ_(j<=k) a_j b_(k-j)`; `q = a/b`: `q_k = (a_k -
Σ_(1<=j<=k) b_j q_(k-j)) / b_0`; `r = √a`: `r_0 = √a_0`, `r_k = (a_k -
Σ_(1<=j<k) r_j r_(k-j)) / (2 r_0)`; `c + i s = exp(i a)`: `c_k = -(1/k)
Σ_(1<=j<=k) j a_j s_(k-j)`, `s_k = (1/k) Σ_(1<=j<=k) j a_j c_(k-j)`. A
Bernstein polynomial of a series argument is evaluated by de Casteljau on
series. For every point `x` of `X` these recurrences produce the exact
coefficients at `x` from the exact coefficients of the inputs; interval
operations are inclusion isotone, so evaluating them on the inputs'
enclosures at `X + ε` encloses the coefficient at every `x ∈ X`, in
particular at `ξ`. The same evaluation certifies smoothness: a division
runs only when `b_0` excludes zero and a square root only when `a_0` is
positive, over all of `X`, so the integrand is analytic on the piece. The
integrands are written once over both kinds of number: the node values
evaluate on plain enclosures, the remainders on series.

**What is integrated.** (a) Along a rational spline pcurve piece on a
nonrational patch (the patch-exact route of the moments), `∫_0^1 M/W^k`,
`M` and `W` Bernstein. (b) Along a rational spline pcurve on a plane, and
any spline pcurve on a cylinder, cone, sphere or torus, `-∫_0^1 F(u(τ),
v(τ)) u'(τ) dτ`, `u = U/W`, `v = V/W`, `u' = (U'W - UW')/W^2`, `F` the
closed-form antiderivative in `v`: a polynomial in `(u, v)` (plane; the
pcurve is translated to the face's reference point exactly, on its
homogeneous poles, so a far face keeps its scale), polynomials in `v` times
`cos^a u sin^b u` (cylinder, cone), or `cos^a u sin^b u` times the Fourier
form
`α_0 (v - v_l) + Σ_f (α_f (sin f v - sin f v_l) + β_f (cos f v_l - cos f
v)) / f` of `∫_(v_l)^v cos^c sin^d` (sphere, torus). (c) On a spline
surface, the terms that are not tensor polynomials: every term on a
rational patch, the four `|N|` terms on a nonrational one. With `G(u, v) =
∫_(v_a)^v f(u, s) ds` from the domain's start, a boundary piece `(u(τ),
v(τ))` in the patch `[u0, u1] × [v0, v1]` contributes `-∫ G du =
-Σ_(Q below) ∫_(u(0))^(u(1)) ∫_(Q) f ds du - ∫_0^1 u'(τ) ∫_(v0)^(v(τ))
f(u(τ), s) ds dτ` (the first term by the change of variables `u = u(τ)`,
an exact differential whatever `u`'s monotonicity). In the patch's local
coordinates `ū = (u - u0)/Δu`, `v̄ = (v - v0)/Δv` every mass integrand is
`g(S)·N` or `g(S) |N|`, homogeneous of degree one in `N = S_u × S_v =
N̄ / (Δu Δv)` with `N̄ = S_ū × S_v̄`, so `f du dv = f̄ dū dv̄` exactly, `f̄`
the same integrand with `N̄`. Substituting `s̄ = σ v̄(τ)`, the piece's own
term is `J = -∫_0^1 ∫_0^1 ū'(τ) v̄(τ) f̄(ū(τ), σ v̄(τ)) dσ dτ`, and each
patch `Q` below contributes the same `J` of the line from `ū_Q(u(0))` to
`ū_Q(u(1))` at `v̄ = 1`. The swept region lies in the patch, where the
homogeneous surface `h(ū, v̄) / w(ū, v̄)` is analytic while `w > 0` (checked
over each box), so the tensor rule above applies with the series along
`τ` (`σ` an interval) and along `σ` (`τ` an interval). A piece must lie in
one patch: lines are split exactly where they cross knot lines, a spline
pcurve piece is located by its control points, a chord (a gap between
pcurve ends) by certain ends; anything else falls back to the strips.

On a plane the first form would also serve, `M/W^k` with `M` the
Bernstein expansion of `F(U/W, V/W) W^d (U'W - UW')`, but its degree is
the integrand's times the pcurve's and its node values in binary64 were up
to `6e-12` wide relative to themselves on the rounded corner's caps, where
`F(u, v) u'` from the pcurve's own coordinates stays at rounding.

**Adaptivity.** The first rule on a whole piece gives `s_k = Σ W |f_k|`
at its nodes, the piece's absolute integral of each integrand, and `r_k`,
the rounding width of its node sum. A piece or box of relative size `a`
(its length or area as a fraction of the piece's) is accepted when each
remainder is at most `a max(2^-44 s_k, 4 r_k)`; the remainders are formed
first and the node sum only for an accepted piece. Otherwise it is halved
across the direction with the larger remainder against that budget, up to
twelve halvings per direction, where it is accepted as it is. A series
undefined over a whole piece (an enclosure too wide to exclude a zero, as
`|N|^2` over a whole rational patch) sends the piece to its halves; after
twelve halvings in all it is taken for a singularity (a degenerate patch
edge where `|N|` vanishes), and after 2,048 pieces or boxes in one integral
the work is too; either way the rule gives up and the first-order route
runs. The
acceptance rule only decides the work: every accepted enclosure is sound
by the remainder above, and the result's width is reported, not assumed.
On the fixtures the widths come from rounding (`2^-40` and `2^-44` gave
the same enclosures), `7.6e-15` of the volume, `5.6e-15` of the area,
`2.3e-14` of `V^(1/3)` for the centroid and `1.6e-13` of `V^(5/3)` for the
inertia at worst, where the first-order routes gave `2e-3`–`5e-3` of the
nonrational walls' areas and `5e-2` to `0.6` on the rational corner.

## Profiles with circular arcs (S5)

**Sectors.** An arc of centre `c` from `a` to `b` contains the direction `w`
exactly when, with `u, v` its bounding directions turning counter-clockwise,
`u × w >= 0` and `w × v >= 0` for a minor sweep (`u × v > 0`), not both
negative for a major one, and `u × w >= 0` for a half turn: rational signs.

**Distances.** A point is within `t` of an arc when it is within `t` of an
end or its direction lies in the sector and `r - t <= |p - c| <= r + t`,
decided by the exact sum-of-terms predicates (`r` is an input). Between a
segment and an arc the minimum lies at an end of either, at the pair through
the foot of the centre on the segment (`||f - c| - r| <= t`, squared), or at
a crossing; between two arcs at an end, at one of the four pairs along the
line of centres (`|d ± r2 ∓ r1| <= t`, the distance between centres against
sums of radii, exact) or at a crossing. Crossings are irrational: they are
located in rational intervals with square roots, and an undecided location
counts as touching. Adjacent pieces through the shared point `v` meet again
at `v`'s reflection (along the segment, or across the line of centres),
rational, which must not lie on both pieces farther than `2t` from their
shared points; tangents exactly opposite are a doubling back.

**Moments.** With the Green forms `A = ½∮(x dy - y dx)`, `∫x = ∮x²/2 dy`,
`∫y = -∮y²/2 dx`, `∫x² = ∮x³/3 dy`, `∫y² = -∮y³/3 dx`, `∫xy = ∮x²y/2 dy` on
every segment alike (the polygon's symmetric per-edge forms agree only over
a closed loop), an arc `x = cx + r cos θ`, `y = cy + r sin θ` gives
polynomials in `cos θ, sin θ` of degree at most 4, integrated by the
reduction formulas. The area screen encloses `½(Σ p × q + Σ r²(φ - sin φ))`
in rational intervals, `φ` the certified sweep and `sin φ = u × v / |u||v|`.

**Ray parity.** An arc splits at the circle's top and bottom (exact points)
into parts monotone in `y`; a part between heights on either side of the
point's crosses once, right of the point when `±sqrt(r² - (y - cy)²) >
x - cx`, the sign by the half it runs on, decided by squaring.

## Measures of sheets and wires (S6)

A sheet's area and centre sum its faces' `|N|` and `p |N|` integrals from
the certified mass module, each face's area taken positive and its moment
with the same sign, and divide in the certified tiers. A wire's length and
first moment about the body's reference point `r` are closed forms per edge:
a segment `a b` has length `|b - a|` and moment `|b - a| ((a + b) / 2 - r)`;
an arc of radius `R` about `o` in the frame `(x, y)` from angle `a0` through
the sweep `s` has length `R |s|` and moment

    R sign(s) ( (o - r) s + R ((sin a1 - sin a0) x + (cos a0 - cos a1) y) ),
    a1 = a0 + s,

with `cos` and `sin` enclosed by the certified series; a whole circle has
length `2 pi R`, `pi` enclosed by its certified value, and its centre is
`o`. An acorn measures zero at its vertex. Every result is an enclosure; the
native comparison requires it to contain `BRepGProp`'s value up to `1e-9`
of the row's magnitude (OCCT's own integration error on spline faces).

## Analytic surface intersections (S7a)

Surfaces are exact point sets of their stored data, with axes normalised
exactly: `|p - o|^2 - ((p - o) . a)^2 / |a|^2 = r^2` for a cylinder, a
rational quadric; a cone's `rho cos a = +-(r cos a + h sin a)` (both
nappes) with `h` along the unit axis, transcendental through `cos a` and
`sin a`. Degeneracies are exact rational predicates: `n1 x n2 = 0` for
parallel planes, `a . n = 0` and `a x n = 0` for a cylinder's axis parallel
or normal to a plane, `|off|^2 / |n|^2` against `r^2` for tangency,
coplanar axes by a triple product, a cone's apex through a plane only when
the apex is rational (radius zero at the origin) or the plane contains the
axis.

* **Plane/cylinder.** The axis meets the plane at `c + a lambda`,
  `lambda = (o - c) . n / (a . n)`; the section is an ellipse with semi-minor
  `r` along `a x n` and semi-major `r |a| |n| / |a . n|` along
  `n x (a x n)`.
* **Plane/cone.** With the apex `V = c - a r / (|a| tan a)`, the plane's unit
  normal `n`, `cos b = n . a / |a|`, `sin b = |e1| / |a|` for the axis's
  projection `e1` on the plane, and `D = (V - o) . n`, a point
  `V - D n + s u1 + t u2` (`u1 = e1 / |e1|`, `u2 = n x u1`) lies on the cone
  when `|u|^2 cos^2 a = (u . a)^2 / |a|^2`, that is
  `A (s - s0)^2 + cos^2(a) t^2 + K = 0` with `A = cos^2 b - sin^2 a`,
  `s0 = -D cos b sin b / A` and `K = -D^2 cos^2 a sin^2 a / A` (the
  constant simplifies because `(cos^2 a - cos^2 b)(cos^2 b - sin^2 a) -
  cos^2 b sin^2 b = -cos^2 a sin^2 a`). `A > 0` is an ellipse with semi-axes
  `|D cos a sin a| / A` along `u1` and `|D sin a| / sqrt(A)` along `u2` (the
  first is the larger since `A <= cos^2 a`); `A < 0` a hyperbola with the
  same expressions as semi-transverse and semi-conjugate axes; `A = 0` is
  impossible. Through the apex (`D = 0`), `A > 0` leaves the apex and `A < 0`
  two generatrices `cos p u1 +- sin p u2`, `cos p = cos a / sin b`.
* **Sphere/sphere and parallel cylinders.** The radical plane (or line):
  the centre `c1 + d k`, `k = (|d|^2 + r1^2 - r2^2) / (2 |d|^2)`, radius
  squared `r1^2 - k^2 |d|^2`, all rational; tangency when it is zero.
* **Crossing cylinders of equal radius.** The two ellipses lie in the planes
  through the crossing point spanned by the common normal and a bisector
  `u1 +- u2`, with normal `u1 -+ u2`: semi-minor `r` along the common
  normal, semi-major `r` over the cosine between an axis and that normal.
* **Coaxial pairs.** With `h` along the common axis, a cylinder is
  `rho = r`, a cone `rho = c + t h` (`t = tan a`), a sphere
  `rho^2 = R^2 - (h - h0)^2`; circles lie where `rho1 = +-rho2`, solved as
  linear or quadratic equations in `h`. Equal slopes decide exactly (the
  same half-angle); constants agree exactly only when the rational and the
  `tan a` parts agree separately (the same surface); a rational apex on the
  other surface is a point.

## Procedural intersection curves (S7b)

On the ruled cylinder (axis `a` unit, radius `r`), in the frame whose `x` is
the unit common normal towards the other surface, a ruling
`P(u) + v a`, `P(u) = o + r (cos u x + sin u y)`, meets the other quadric
where `A v^2 + 2 B v + C = 0` and the curve is `v = (-B +- sqrt(D)) / A`,
`D = B^2 - A C`.

* **Two cylinders.** With the other axis `a2` (unit), `M = 1 - a2 a2^T`,
  `A = |M a|^2 = |a x a2|^2`, `B = M(P - o2, a)`, `C = |M (P - o2)|^2 - r2^2`.
  Lagrange's identity gives `A |M w|^2 - (M w . M a)^2 = |M a x M w|^2 =
  (w . (a2 x a))^2` for `w = P - o2`, and in this frame
  `w . (a2 x a) / |a2 x a| = r cos u - d` with `d` the axes' distance, so
  `D = A (r2^2 - (r cos u - d)^2)`. The curve exists where
  `d - r2 < r cos u < d + r2`; with `r <= r2` the upper bound always holds,
  so a loop is `cos u > c = (d - r2) / r`, two rings `c < -1`, a
  figure-eight `c = -1` and a tangent point `c = 1`.
* **A cylinder and a sphere.** `A = 1`, `D = R^2 - dist(c, ruling)^2 =
  R^2 - e^2 - r^2 + 2 e r cos u` with `e` the centre's distance from the axis
  (`x` points towards it): a loop is `cos u > c = (e^2 + r^2 - R^2) / (2 e r)`.
* **Enclosures.** `t = arccos c = atan2(sqrt(1 - c^2), c)` in rational
  intervals (`certified::atan2`); points `P(u) + v a` with `cos u`, `sin u`
  and the square root enclosed, binary64 intervals first.
* **A sphere and a cone (S7b.2).** On the cone's rulings through the apex
  `V`, `d(u) = cos h a + sin h (cos u x + sin u y)` is a unit vector, so
  `A = 1`, `C = |V - c|^2 - R^2` and, with `x` towards the centre,
  `B(u) = (V - c) . d(u) = b0 + b1 cos u`, `b0 = cos h (V - c) . a`,
  `b1 = sin h (V - c) . x`. `D = B^2 - C > 0` where `B > sqrt(C)` or
  `B < -sqrt(C)`: `cos u > k` or `cos u < k` with
  `k = (+-sqrt(C) - b0) / b1` (the inequality reversing with the sign of
  `b1`), an arc `[-arccos k, arccos k]` or `[arccos k, 2 pi - arccos k]`, or
  nothing, or the whole turn.
* **A cylinder and a cone (S7b.2).** The cone is `f(p) = cos^2 h |p - V|^2 -
  ((p - V) . a2)^2` about its apex (both nappes), so with the cylinder's
  ruling `A = cos^2 h - (a . a2)^2` is constant and `D = B^2 - A C` a
  trigonometric polynomial of degree two in `u`; `D' = 2 B B' - A C'` with
  `B' = Q(P', a)`, `C' = 2 Q(P', P - V)`. Over a piece `U` with middle `m`,
  `D(U) is in D(m) + D'(U) (U - m)` (the mean-value form), whose width is
  `|D'| |U|`: beside a nearly tangent loop, where `D` dips quadratically,
  the natural interval extension needs pieces as small as the square of
  their distance from the loop, the mean-value form only as small as the
  distance.
* **A torus and a plane or a sphere (S7b.3a).** On the meridian circle
  `p(t) = C + r (cos t e + sin t a)`, `C = o + R e`, a plane `n . (p - q)`
  is `f0 + alpha cos t + beta sin t` with `f0 = n . (C - q)`,
  `alpha = r n . e`, `beta = r n . a`; a sphere `|p - s|^2 - rho^2` has
  `f0 = |C - s|^2 + r^2 - rho^2`, `alpha = 2 r (C - s) . e`,
  `beta = 2 r (C - s) . a`. With `D = alpha^2 + beta^2 - f0^2` the unit
  vector `(cos t, sin t) = (-alpha f0 -+ beta sqrt(D), -beta f0 +- alpha
  sqrt(D)) / (alpha^2 + beta^2)` needs no trigonometry. Writing
  `c = m cos phi` (`m` the length of the plane normal's, or of `s - o`'s,
  component normal to the axis, `x` along it): for a plane
  `D = r^2 c^2 + r^2 (n . a)^2 - (d0 + R c)^2`, `d0 = n . (o - q)`; for a
  sphere `D = 4 r^2 (R - c)^2 + 4 r^2 w_a^2 - (K - 2 R c)^2` with
  `w = s - o`, `w_a = w . a`, `K = R^2 + |w|^2 + r^2 - rho^2`. Both have
  the leading coefficient `(r^2 - R^2)` (times 4) `< 0`; the plane's
  discriminant is `4 r^2 (d0^2 + (R^2 - r^2)(n . a)^2)` and the sphere's
  `16 r^2 ((K - 2 R^2)^2 + 4 (R^2 - r^2) w_a^2)`, zero only in the special
  cases. `D(phi) = P(m cos phi)` has double roots in `phi` only at
  `phi = 0, pi` (where `P(+-m) = 0`) or at the vertex, so every node and
  tangency is one of the exact classes.
* **Coaxial pairs.** In the half-plane `(rho, z)` the meridian circle
  `(rho - R)^2 + z^2 = r^2` meets a sphere's circle, another torus's
  meridian circle (the radical line; tangent when the centres' distance is
  `r1 + r2` or `|r1 - r2|`, rational squares), a cylinder's line
  `rho = r_c`, and a cone's lines `rho = s (r_c + g (z - z_a) tan h)`
  (`s = +-1` the nappe, `g = +-1` the axis's direction): the quadratic
  `(k1^2 + 1) z^2 + 2 k1 (k0 - R) z + (k0 - R)^2 - r^2 = 0` with
  `k0 = s (r_c - g z_a tan h)`, `k1 = s g tan h`.

## Traced torus curves (S7b.3b)

* **The field.** On the parameter torus `G(phi, t) = f(p)` with
  `f(p) = (p - q)^T Q (p - q) - k`, `Q = c2 - u u^T`; with `w = p - q`,
  `G_phi = 2 Q w . p_phi`, `G_t = 2 Q w . p_t` and
  `G_ab = 2 p_a^T Q p_b + 2 Q w . p_ab`, where `p_phi = rho e'`,
  `p_t = r (cos t a - sin t e)`, `p_pp = -rho e`, `p_pt = -r sin t e'`,
  `p_tt = -r (cos t e + sin t a)`, `rho = R + r cos t`.
* **Folds.** Over a box `B` with middle `m`, `G(B)` and `G_t(B)` are in
  their mean-value forms `G(m) + grad G(B) (B - m)`; a box where either
  excludes zero holds no fold. The Krawczyk operator
  `K = m - Y F(m) + (1 - Y J(B)) (B - m)`, `F = (G, G_t)`,
  `J = [[G_phi, G_t], [G_pt, G_tt]]`, `Y` the inverse of `J`'s midpoint,
  maps `B` into its interior only when `B` holds exactly one zero of `F`,
  regular; `K` then encloses it, and iterating narrows it. For a tangency
  the same operator of `grad G` with the Hessian certifies a unique critical
  point in its box.
* **A fold's box.** `G_phi` of one sign on the box makes each piece of the
  zero set a graph over `t`; with no zero on the top and bottom edges its
  pieces end on the sides, and two certified roots on one side and none on
  the other leave exactly one arc, turning at the fold. The two branches at
  `phi` from the fold are about `sqrt(2 |phi - phi_f| |G_phi / G_tt|)`
  apart, which sizes the box's height. A tangency's box: with a unique
  critical point inside (a saddle or an extremum) no closed piece can lie in
  it, so two and two side roots give two branches crossing at the tangency,
  none an isolated point.
* **Tracks.** A window `W` over a step `[phi0, phi1]` holds exactly one root
  for every `phi` of the step when `G_t` keeps a sign on the step times `W`
  and `G` has opposite certain signs on `W`'s ends for every `phi`, the
  latter by the mean-value form in `phi` (the natural extension overestimates
  it by the branch's own motion). The window is predicted from the slope
  `dt/dphi = -G_phi / G_t` with a margin of one and a half times the
  predicted motion.
* **Tangency.** In a rational basis `(U, V)` of the spine's plane the three
  conics are `E0 = |u U + v V|^2 - R^2`, `E1 = M (d + u U + v V) .
  (a x (u U + v V))` (`d = o - o2`, `M` the projector normal to the axis)
  and `E2 = |M (d + u U + v V)|^2 - k`, each quadratic in `v`. With
  `E0 = a2 v^2 + a1 v + a0` and `E1 = b2 v^2 + b1 v + b0`, the resultant is
  `(a2 b0 - a0 b2)^2 - (a2 b1 - a1 b2)(a1 b0 - a0 b1)`, of degree four in
  `u`, and the common root `v = (a2 b0 - a0 b2) / (a1 b2 - a2 b1)`; `E2`'s
  numerator at it vanishes on the algebraic root exactly when the critical
  point is a tangency. `V` is sheared (`V + s U`, `s` a small rational)
  until no critical point makes the denominator vanish (two critical points
  with one `u`).
* **Two tori (S7b.3b.2).** For the other torus's quartic
  `f = (|w|^2 + K)^2 - 4 R2^2 |M w|^2` (`K = R2^2 - r2^2`, `M` the projector
  normal to its axis) the jets use `grad f = 4 (|w|^2 + K) w - 8 R2^2 M w`
  and the Hessian form `4 (|w|^2 + K) v . z + 8 (w . v)(w . z) -
  8 R2^2 v . M z`. For the tangencies, with `s1` on the first spine,
  `w = s1 - o2`, `rho = |M w|` and `Q = k - |w|^2 - R2^2`, the distance to the
  second spine critical along it is `|w|^2 + R2^2 + 2 sigma R2 rho`; setting
  it to `k` gives `2 sigma R2 rho = Q`, the criticality along the first spine
  `w . T1 + sigma R2 (M w . T1) / rho = 0` becomes
  `E1 = (w . T1) Q + 2 R2^2 (M w . T1)` (a cubic: `w . T1` is linear, as
  `(u U + v V) . (a x (u U + v V)) = 0`), and the distance
  `E2 = Q^2 - 4 R2^2 rho^2`. Modulo `E0 = a2 v^2 + a1 v + a0` (`a2` a nonzero
  constant) both are linear in `v`, `c0 + c1 v` and `d0 + d1 v`; the common
  roots of `E0` and `E1` are the real roots of `a2 c0^2 - a1 c0 c1 + a0 c1^2`
  with `v = -c0 / c1`, and `E2` vanishes there exactly when `d0 c1 - d1 c0`
  does. The second spine's point is `o2 - sigma R2 M w / rho`, `sigma` the
  sign of `Q`.
* **A cone's rulings (S7b.4).** With `d(u) = cos h a + sin h (cos u x +
  sin u y)`, `d_u = sin h (cos u y - sin u x)` and `d_uu = -sin h (cos u x +
  sin u y)`: `A = Q(d, d)`, `A_u = 2 Q(d, d_u)`, `A_uu = 2 Q(d_u, d_u) +
  2 Q(d, d_uu)`, `B = Q(V - q, d)` with `B_u`, `B_uu` along, and `C` constant.
  Substituting `v = tan(t / 2)` into `A v^2 + 2 B v + C` and multiplying by
  `cos^2(t / 2)` gives `G = ((A + C) + (C - A) cos t) / 2 + B sin t`, of
  degree one in `t`; the factors are `alpha sin psi + beta cos psi` with
  `(alpha, beta) = (A, 2 B)` (the apex on the other surface, `C = 0`) or
  `(2 B, C)` (twins, `A = 0`). On a factor's chart the apex's rulings are
  the roots of `B`: `b1 cos u + b2 sin u = -b0`, two when
  `b1^2 + b2^2 > b0^2`, at `u = atan2(b2, b1) +- arccos(-b0 / sqrt(b1^2 +
  b2^2))`.
* **Points.** The interval Newton operator `N(T) = m - G(m) / G_t(T)`
  contains every root of `G(phi, .)` in `T`; intersected with `T` it
  converges quadratically.

## Lines and circles against surfaces (S7c.1)

Along a line `p0 + s d` with rational `p0`, `d`, a plane's, cylinder's,
sphere's and torus's implicit functions (`n . w`, `|w|^2 - (w . a)^2 / |a|^2
- r^2`, `|w|^2 - r^2`, `(|w|^2 + R^2 - r^2)^2 - 4 R^2 (|w|^2 - (w . a)^2 /
|a|^2)` with `w = p - o`) are rational polynomials in `s`. A real root's
multiplicity (from the repeated gcds with the derivative) is the order of
contact: two or more is a tangency.

A circle's points are `o + u U + v V` with `E0 = |U|^2 u^2 + 2 (U . V) u v +
|V|^2 v^2 - r^2 = 0`. Reducing the surface's function `F` modulo `E0` (as a
polynomial in `v`, `E0` monic up to the constant `|V|^2`) leaves
`c0(u) + c1(u) v`, equal to `F` on the circle. The resultant in `v` of `E0`
and `c0 + c1 v` is `a2 c0^2 - a1 c0 c1 + a0 c1^2`; its leading coefficient
in `v` being constant, the order of a root `u0` is the sum of the
intersection multiplicities of `E0 = 0` and `F = 0` over the points above
`u0`, and `F + g E0` has the intersection multiplicities of `F`. Where
`c1(u0) != 0` there is one point above `u0`, `v = -c0 / c1`, real, so the
root's order is that point's multiplicity; a shear `V + k U` (rational `k`)
moves the roots where `c1` vanishes. `c0 = c1 = 0` identically is
`F = 0` on the whole circle; `E0` being irreducible, a zero resultant
implies it.

The same holds for any plane conic `E0` whose square terms do not both
vanish in the sheared basis: an ellipse or a hyperbola of stored axes `x`,
`y` (S7c.2), `minor^2 xi^2 +- major^2 eta^2 = major^2 minor^2` with `xi = u
+ k v`, `eta = v`; a shear making the hyperbola's `v^2` term vanish (`k =
major / minor`) is skipped.

**Two plane conics (S7d.1).** In one plane the second conic's equation,
written in the first's plane by exact frame coordinates (the rows of
`[x y (x × y)]^-1`, rational), is a quadratic function there, so the
resultant above applies with it as the surface's function; `c0 = c1 = 0` is
the same conic. In planes crossing along `x0 + t d` (`d = n1 × n2`, `x0` the
rational solution of `n1 . x0 = n1 . o1`, `n2 . x0 = n2 . o2`, `d . x0 = 0`),
a common point lies on that line, where each conic's equation is a
polynomial `q_i(t)` of degree at most two: the points are the real roots of
`gcd(q1, q2)`. Both curves are tangent to the line there exactly when the
root is double in both, and then they share their tangent line; otherwise
their tangents differ (a conic's tangent at a simple root of `q_i` is not the
line, and a curve in another plane cannot share any other tangent).

**A spline against a plane conic (S7d.2).** With `D = X - o W` on a span,
`L = n . D` vanishes where the spline meets the conic's plane and `Q = W^2
E(X / W)` (the conic's quadratic equation made homogeneous: `|D|^2 - r^2 W^2`,
or `minor^2 xi(D)^2 +- major^2 eta(D)^2 - major^2 minor^2 W^2` with the
rational frame rows) where it meets the conic's quadric cylinder, so the
common points are the common roots, the roots of `gcd(L, Q)`. With `W > 0`,
`d/ds (n . (X / W - o)) = L' / W` and likewise for `Q` at a root: the spline's
tangent lies in the plane exactly when the root is multiple in `L` and along
the conic exactly when it is multiple in `Q`, both exactly when it is at
least double in the gcd.

**A spline against a cone.** With `h = (D . a) / (|a| W)` and `rho^2 =
(|D|^2 - (D . a)^2 / |a|^2) / W^2` for `D = X - o W`, the cone `rho^2 = (r +
h tan a)^2` is `W^2 F = Q0 + tau Q1 + tau^2 Q2 = 0` with rational `Q0 =
|D|^2 - (D . a)^2 / |a|^2 - r^2 W^2`, `Q1 = -2 r (D . a) W`, `Q2 = -(D . a)^2`
and `tau = tan a / |a|`. For rational `a != 0`, `e^{2 i a}` is transcendental
(Lindemann-Weierstrass), hence `tan a` and `tau` (`|a|` algebraic). A
polynomial identity `sum tau^i Q_i = 0` in the curve parameter makes each
coefficient a polynomial in `tau` with rational coefficients vanishing at
`tau`, so all of them vanish: `F = 0` on a span exactly when `Q0 = Q1 = Q2 =
0`. At a common root `s0` of the `Q`s the coefficient of `(s - s0)^m` in `F`
is `sum tau^i c_{i,m}`, zero exactly when every `c_{i,m}` is: the order of
`s0` in `F` is its least order in the `Q`s, its order in their gcd.

**Certified logarithms.** `ln x = k ln 2 + 2 atanh((m - 1) / (m + 1))` for
`x = 2^k m`, `m` in `[2/3, 4/3]` (the series argument at most `1/5`, its
remainder below `25/24 |y|^(2N+1) / (2N + 1)`), `ln 2 = 2 atanh(1/3)`;
`asinh x = ln(x + sqrt(x^2 + 1))`, odd. Both are monotone, so an interval's
image is the hull of its ends'.

## C1 of spline cells (R4)

A cell is C1 in its own parameterisation (U3 of `REVIEW_NOTES.md`). For a
B-spline of degree `p`, a knot of multiplicity `m` leaves the curve
`C^{p-m}` there, and one copy of the knot is removable with zero residual
exactly when the curve is `C^{p-m+1}` (Tiller's removal condition). A
binary64 spline has interior multiplicities `m <= p`, so it is C1 by
construction wherever `m < p`, and at `m = p` it is C1 exactly when one
removal succeeds. The kernel runs that removal
(`ExactBSplineCurve3::removable`, the exact homogeneous inverse insertion
with the resulting weights' signs not examined) on the homogeneous poles
`(w x, w y, w z, w)`: it decides C1 of the homogeneous curve, which implies
C1 of the rational curve `N / w`; a rational joint that is C1 only because
its weights compensate is reported. A periodic seam is an interior knot. A
removal the basis cannot hold (a periodic basis left with no more poles than
its degree) first inserts a simple knot in the middle of the longest span,
which changes neither the curve nor its continuity at the tested knot. A
surface is C1 across a knot line when every control row across it is,
since the surface is a combination of those rows with linearly independent
basis functions in the other parameter.

The independent reference (`spline_cell_reference.py`) never removes a
knot: it compares the homogeneous curve's position and first derivative
from both sides of each knot, by de Boor's recurrence on the two spans in
`Fraction`s. The seam's left side is the domain's end.

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

## Tessellation bounds (T-a)

A mesh triangle has nodes `X_i` and parameter points `p_i` on its face's
surface `S` (lifted into one period; a pole's `u` is free). Its reported
bound certifies a map from the flat triangle onto `S` that moves no point
farther (`TESSELLATION.md`); three maps are certified and the smallest bound
is kept.

**Linear map.** A point with weights `λ` goes to `S(p)`, `p = Σ λ_i p_i`.
With `d_i = p_i - p` (so `Σ λ_i d_i = 0`), Taylor's formula with integral
remainder gives

    Σ λ_i S(p_i) - S(p) = Σ λ_i ∫_0^1 (1 - s) D²S(p + s d_i)[d_i, d_i] ds.

The segments `p + s d_i` stay in the triangle. If `a >= |S_uu|`,
`b >= |S_uv|`, `c >= |S_vv|` over its parameter box, `|D²S[w, w]| <=
a w_u² + 2 b |w_u w_v| + c w_v²`. `Σ λ_i d_{u,i}²` is the variance of a
distribution on an interval of length `U` (the triangle's `u` extent), at
most `U² / 4`, and by Cauchy–Schwarz `Σ λ_i |d_{u,i} d_{v,i}| <= U V / 4`.
With `∫ (1 - s) ds = 1/2`:

    |Σ λ_i S(p_i) - S(p)| <= (a U² + 2 b U V + c V²) / 8,

and the flat point `Σ λ_i X_i` is within the largest gap `|X_i - S(p_i)|`
more. On a triangle with one parameter edge (`V = 0`) this is the familiar
`M h² / 8` of a chord.

**Coefficients.** With `E(u) = cos u x + sin u y` from the stored frame and
`σ` a bound on the spectral norm of `[x y n]` (Gershgorin on its Gram
matrix, about one), any combination `α x + β y + γ n` is at most
`σ |(α, β, γ)|`, so

| Surface | `a` | `b` | `c` | `|N_u|` | `|N_v|` |
| --- | --- | --- | --- | --- | --- |
| plane | 0 | 0 | 0 | 0 | 0 |
| cylinder `r` | `σ r` | 0 | 0 | 1 | 0 |
| cone `R`, `α` | `σ max |R + v sin α|` | `σ |sin α|` | 0 | `|cos α|` | 0 |
| sphere `R` | `σ R sup |cos v|` | `σ R sup |sin v|` | `σ R` | `sup |cos v|` | 1 |
| torus `R`, `r` | `σ (R + r sup |cos v|)` | `σ r sup |sin v|` | `σ r` | `sup |cos v|` | 1 |

over the box's `v` range; `|cos|` and `|sin|` are bounded from the box's
middle by their Lipschitz constant 1. The normal turns by at most
`|N_u| U + |N_v| V` over the triangle (the angle between unit normals is at
most the length of their path), with a `1e-12` relative allowance for the
frame's departure from orthonormality.

**Tangential correction.** On each analytic surface `S_uv` is parallel to
`S_u` (`S_u = g(v) E'(u)`, `S_uv = g'(v) E'(u)`), and the unit normal
`N ∝ S_u x S_v` is orthogonal to it, whatever the frame. So
`N(p) . S_uv(ξ) = (N(p) - N(ξ)) . S_uv(ξ)`, at most `t b` for the
triangle's normal turn `t` (at most 2): the normal part of the linear map's
deviation `e = Σ λ_i S(p_i) - S(p)` is at most
`(a U² + 2 t b U V + c V²) / 8`. Its tangential part `e_t = S_u δu + S_v δv`
has `|δu| <= |e| / |S_u|` and `|δv| <= |e| / |S_v|` (orthogonal
derivatives; the frames' defect of a few units in the last place is
covered by a `1e-9` relative allowance), and mapping the point to
`S(p + δp)` instead leaves `|e_n| + |S(p + δp) - S(p) - J δp|`, the second
term at most `(a' δu² + 2 b' |δu δv| + c' δv²) / 2` with the coefficients
over the box widened by `δv`. `|S_u|` is bounded below over the box (`r`;
`|R + v sin α|` at the nearer end when it keeps its sign; `R cos v` at an
end, concave on the band; `R + r (cos v_m - h)`); where it may vanish (a
pole in the box) the correction is not used. On cones and tori the mixed
term dominates elongated triangles, and this removes it to first order.

**Fan map at a pole.** A triangle with vertex `k` at a pole and the others
at `(u_m, v_m)`, `(u_n, v_n)`: the point with weights `λ` goes to
`S(u_μ, Σ λ_i v_i)`, `u_μ` interpolating `u_m`, `u_n` in the ratio
`λ_m : λ_n`. It is continuous (every `u` maps to the pole) and linear on the
opposite edge, so it agrees with a neighbour's linear map there. Writing
the flat point as `(1 - s) X_k + s Q_μ`, `Q_μ` on the opposite chord, its
distance from the image is at most `s |Q_μ - S(u_μ, v_μ)|`, the opposite
edge's interpolation, plus the chord of the meridian `u = u_μ` from the pole
to `v_μ`, at most `c V² / 8`. On a cone the meridian is a ruling and `c = 0`:
an apex fan meets its base segment's bound.

**Edges.** An arc of radius `r` and sweep `φ` has `|C''| <= σ r φ²` in its
fraction, so a segment of fraction step `Δt` deviates by at most
`σ r (φ Δt)² / 8` from its chord, plus its end nodes' distances from the
exact curve points; a line segment by those only. Its tangent turns by
`|φ| Δt`.

**Counts.** An edge takes the least uniform count `n` with
`σ r (φ / n)² / 8 <= 0.9 δ` and `|φ| / n <= θ`, and for each curved face
using it `(a Δu² + 2 b Δu Δv + c Δv²) / (8 n²) <= 0.45 δ` and
`(n_u Δu + n_v Δv) / n <= 0.45 θ` over its pcurve's extents: a triangle on a
boundary segment has at least the segment's extents, so the segment must
leave room for the triangle's own bound.

**Evaluation.** Node positions are the midpoints of outward-rounded binary64
enclosures of `S(p)` or `C(t)` (`Fast`, with its certified `cos_sin` and
`2π`), their gaps the enclosures' reach; a vertex's gap is its distance from
the enclosed exact point; every bound is summed and multiplied in the same
tier and rounded up.

## Conic sections and certified projections (S8d.2)

**A cone's plane section.** In the cone's frame the plane is `F = a u + b v
+ c w + d`, `|m|^2 = a^2 + b^2 + c^2`; the cone has half-angle `α` (`tan α =
(r1 - r0) / H`), apex `V` (virtual on a frustum) and axis `k` pointing from
`V` into the solid. With `N` the unit normal, `h = F(V) / |m|`, `e1` the
axis's direction projected on the plane (`κ1 = |k - (k.N) N|`, `κn = k.N`),
`e2 = N x e1` and the origin at `V`'s foot `O = V - h N`, a point `O + x e1 + y
e2` lies on the cone when `(x κ1 - h κn)^2 = (h^2 + x^2 + y^2) cos^2 α`:

    A x^2 + C y^2 + D x + F0 = 0,  A = cos^2 α - κ1^2,  C = cos^2 α,
    D = 2 h κ1 κn,  F0 = h^2 (cos^2 α - κn^2),

symmetric about `e1`. The kind is exact: `sign(A) = sign(H^2 |m|^2 - (a^2 +
b^2)(H^2 + (r1 - r0)^2))` on the stored binary64 data. `A > 0` is an
ellipse centred at `x = -D / 2A` with semi-axes `sqrt(G / A)` along `e1` and
`sqrt(G / C)` across (`G = D^2 / 4A - F0`, so the first is the major one);
`A < 0` a hyperbola with `G < 0`, transverse semi-axis `sqrt(G / A)` along
`e1`, conjugate `sqrt(-G / C)`, its branch the one where `x κ1 - h κn > 0`;
`A = 0` the parabola `x = -F0 / D - (C / D) y^2`, focal distance `|D| / 4C`.
A plane through the apex (`F(V) = 0` exactly, or within the resolution of a
frustum's virtual apex) cuts the rulings through its rim crossings. A
sphere's section is the circle of radius `sqrt(R^2 - d^2 / |m|^2)` about the
centre's foot.

**Crossings and arcs.** The plane meets an end circle of radius `r` at height
`w` where `k + r |(a, b)| cos(φ - atan2(b, a)) = 0`, `k = c w + d`: it misses,
touches or crosses as `k^2` is above, equal to or below `r^2 (a^2 + b^2)`
(exact). The crossings' parameters on the conic, sorted (cyclically on an
ellipse or circle), cut it into arcs; an arc belongs to the solid when its
midpoint's height lies between the ends. The wall's side on an arc's left
(about the outward normal) is the side the plane's normal points to from
`n_s x T` at its midpoint.

**Taylor jets.** A jet is the truncated Taylor series `sum c_k s^k` of a
function of the fraction about a base, each coefficient an enclosure in the
certified tiers. Sums, products, quotients (`d_k = (n_k - sum_{j<k} d_j
e_{k-j}) / e_0`), square roots, `exp` (and `cosh`, `sinh`), `cos` and `sin`
(by their coupled recurrences) and `atan2(y, x)` (the integral of `(x y' - y
x') / (x^2 + y^2)` from a base angle) follow the usual recurrences. With the
base an interval `I`, `c_k` encloses `f^(k)(t) / k!` for every `t` in `I`.
`Real::exp` halves its argument to `|x| <= 1/2`, sums 14 terms with the
remainder bound `2 |x|^15 / 15!` and squares back.

**Projections.** A `Projection` pcurve evaluates its edge's point jet (a
conic's closed form in `cos`, `sin`, `cosh`, `sinh` of the parameter) in
the surface's frame and inverts the surface map: `u = atan2(y, x)` (turned
by the recorded lift, `ref + atan2` of the vector rotated back by `ref`, so
the branch cut stays opposite), `v = z`, `z / cos α`, `atan2(z, ρ)` or
`atan2(z, ρ - R)` as the surface needs.

**Certified quadrature.** `∫_0^1 g` over a piece `[lo, hi]` with midpoint
`m` is the Taylor polynomial of order `n` about `m` integrated exactly over
`[lo - m, hi - m]` plus a remainder at most `sup_I |c_{n+1}| (|lo - m|^{n+2}
+ |hi - m|^{n+2}) / (n + 2)`, `c_{n+1}` the interval jet's next coefficient
over the piece. Pieces are bisected (at most 40 times) until every remainder
is within its share of `1e-12` (absolute for areas, fluxes and other sign
decisions; relative to each integrand's largest value at eight points, at
least one, for mass moments); several integrands share one evaluation of
the projection per piece. Interval jets over a piece overestimate their
high coefficients geometrically (about fourfold per order through the
quotient in `atan2`), so order 14 is used: the remainder's power of the
piece's width outruns the overestimation. The pieces gather where the
section's angle about the axis turns fast (near an apex or pole),
logarithmically in its distance.

**A torus's plane section (S8d.3).** With `alpha = a cos u + b sin u`, the
plane on the torus is `F = C + W cos(v - psi)`, `C = R alpha + d`, `(W cos
psi, W sin psi) = r (alpha, c)`; the section exists where `D = W^2 - C^2 =
-(R^2 - r^2) alpha^2 - 2 R d alpha + r^2 c^2 - d^2 >= 0`, between the roots
`alpha_{1,2} = (-R d -+ r sqrt(c^2 (R^2 - r^2) + d^2)) / (R^2 - r^2)`
(always real). `alpha` ranges over `[-rho, rho]`: both roots beyond it give
two loops over every `u` (about the axis), one inside a single `u`-interval
whose ends meet the branches `v = psi +- acos(-C / W)` in one contractible
loop, both inside two intervals, their loops about the tube when `C` has
opposite signs at the two roots (the branches meet at `psi` at one end and
`psi + pi` at the other) and contractible otherwise. Each comparison is the
sign of `x + y sqrt(e) - z sqrt(f)` on rationals, decided by squaring twice.
Over `v` the section is `u = atan2(b, a) +- acos(q / rho)`, `q = -(c r sin v
+ d) / (R + r cos v)`, valid for every `v` exactly when the loops wind about
the tube. A cap's branches over `u` are joined round each turning point
(where `D = 0`, a vertical tangent) by a graph over `v` from the points
where `|dv / du| = 1` nearest it (`dv/du = -F_u / F_v`); between them the
slope never vanishes, so the graph over `v` has no turning point there.

**Holes in wound faces.** An unwound loop in a face wound in `u` lies inside
when the signed crossings of the `+v` ray from its first point with the
other loops' line pcurves and chords, over every `u` alias (`+1` where a
piece runs in `-u`, `-1` in `+u`), plus one for a north pole, equal the
face's sign (`+1` forward): the region lies left of its boundary, so
directly below a piece running in `-u`.
## Tessellation bounds for splines (T-b)

T-a's deviation argument needs only that the map is `C1` with bounded
second derivatives on pieces: Taylor's formula with integral remainder
holds for `C1` functions whose derivative is absolutely continuous, so a
segment or triangle may straddle a knot where the cell is `C1` (checked
exactly, R4) and its second derivative jumps. What changes is where
`a >= |S_uu|`, `b >= |S_uv|`, `c >= |S_vv|` (and `|C''|`) come from.

**Cells.** A spline's range is cut into its exact Bézier pieces, and each
piece into cells by de Casteljau halving in the `Fast` tier: every control
of a cell is an interval containing the exact control of the exact
subdivision, since each halving is a convex combination evaluated with
outward rounding. A cell covers the parameter interval `[t_0, t_0 + L]` (a
surface cell `[u_0, u_0 + L_u] × [v_0, v_0 + L_v]`) and is a rational
Bézier curve of degree `p` in `τ = (t - t_0)/L`: `C = A/w`,
`A = Σ B^p_i(τ) w_i P_i`, `w = Σ B^p_i(τ) w_i`, `w_i > 0`.

**Translation.** For any point `c` (the cell's first control),
`A_c = A - w c = Σ B^p_i w_i Q_i` with `Q_i = P_i - c`, and
`A_c = w (C - c)`. `C - c` is the convex combination
`Σ (w_i B^p_i / w) Q_i` (positive weights), so `|C - c| <= R = max |Q_i|`,
and `w >= ω = min w_i`.

**Differences.** The derivative in `t` of a Bézier polynomial is
`p/L Σ B^{p-1}_i Δ_i`, `Δ_i` the first differences of its controls, and
its second `p (p-1)/L² Σ B^{p-2}_i Δ²_i`; each is a convex combination of
the differences and a norm is convex, so `|A_c'| <= p/L max |Δ(wQ)_i|`,
`|A_c''| <= p(p-1)/L² max |Δ²(wQ)_i|`, `|w'| <= p/L max |Δ w_i|` and
`|w''| <= p(p-1)/L² max |Δ² w_i|`. A tensor cell has the same bounds per
direction, and `∂_u ∂_v` of `Σ B^p_i B^q_j X_ij` is
`pq/(L_u L_v) Σ B^{p-1}_i B^{q-1}_j Δ_u Δ_v X_ij`.

**Quotients.** Differentiating `A_c = w (C - c)`:

    A_c'  = w' (C - c) + w C',
    A_c'' = w'' (C - c) + 2 w' C' + w C'',

so `|C'| <= D1 = (|A_c'| + |w'| R) / ω` and
`|C''| <= D2 = (|A_c''| + 2 |w'| D1 + |w''| R) / ω`. For a surface, from
`A_c = w (S - c)`:

    S_u  = (A_u - w_u (S - c)) / w,
    S_uu = (A_uu - 2 w_u S_u - w_uu (S - c)) / w,
    S_uv = (A_uv - w_u S_v - w_v S_u - w_uv (S - c)) / w,

and `S_vv` like `S_uu`: `D_u = (|A_u| + |w_u| R)/ω`,
`a = (|A_uu| + 2 |w_u| D_u + |w_uu| R)/ω`,
`b = (|A_uv| + |w_u| D_v + |w_v| D_u + |w_uv| R)/ω`,
`c = (|A_vv| + 2 |w_v| D_v + |w_vv| R)/ω`. With equal weights the `w`
differences vanish and the bounds are the control net's own differences.

**Segments.** A segment between parameters `u_k` and `u_{k+1}`, `h` apart,
maps `λ` to `C((1 - λ) u_k + λ u_{k+1})`; T-a's argument with `V = 0`
bounds its deviation from the chord of the exact points by `D2 h² / 8`,
`D2` the largest over the cells it meets, and the nodes' gaps add. The
unit tangent `T = C'/|C'|` has `|T'| <= |C''| / |C'|`, so it turns by at
most `D2 h / s` for a lower bound `s` of `|C'|` on the segment: `|C'|` at
the middle (by de Casteljau on the cell in the tier,
`C' = (A_c' - w' (C - c))/w`) less `D2 h / 2`.

**Normals.** With `m = S_u × S_v`, `N = m/|m|` has
`N_u = (m_u - N (N · m_u)) / |m|`, of length at most `|m_u| / |m|`;
`m_u = S_uu × S_v + S_u × S_uv`, so `|m_u| <= M_u = a D_v + D_u b` and
`|m_v| <= M_v = b D_v + D_u c` over the cells. On a triangle's parameter
box of extents `U`, `V` about its centre `o`,
`|m| >= μ = |m(o)| - (M_u U + M_v V)/2`, and when `μ > 0` the normal turns
by at most `(M_u U + M_v V)/μ` along any path in the triangle, as in T-a.
The analytic surfaces' tangential correction and fan map rely on their
structure (`S_uv` parallel to `S_u`, poles) and are not used.

**Numerators (implementation).** The first measurements (the fixtures'
rational quarter-circle wall at 34 times BRepMesh's nodes, 811 against 24,
and their sharp rational corner at 7,900 times, 354,559 against 45) showed
the quotient rule's triangle inequality and `|m_u|` losing to
cancellations the exact expressions keep. Three sharper
certificates are kept beside the ones above, the smallest bound winning;
each works on Bernstein polynomials of the cell (products, differences and
derivatives of enclosed coefficients in the tier, a polynomial bounded by
its largest coefficient):

* *Rational cells* (degrees at most 6 for curves, 3 for surfaces). `C' =
  N_1 / w²`, `C'' = N_2 / w³` with `N_1 = A' w - A w'`,
  `N_2 = (A'' w - A w'') w - 2 w' N_1`, so `|C'| <= max |N_1| / ω²` and
  `|C''| <= max |N_2| / ω³`; for a surface `S_u = N_u / w²`,
  `S_uu = ((A_uu w - A w_uu) w - 2 w_u N_u) / w³`,
  `S_uv = ((A_uv w + A_u w_v - A_v w_u - A w_uv) w - 2 w_v N_u) / w³` and
  `S_vv` likewise. On a ruled rational wall `A_uv w - A_v w_u` vanishes
  identically, and so does its enclosure, up to rounding.
* *Normal cones* (surface cells of degrees at most 6, rational ones at most
  3). `M = A_u × A_v` (nonrational) or `M = N_u × N_v` (rational) is a
  positive multiple of `S_u × S_v`, `w²` or `w⁴` times it in the cell's own
  parameters. With its value `o` at the cell's middle as axis, if every
  coefficient `M_k` has `M_k · o > 0`, every normal of the cell lies in the
  cone of half-angle `α = atan max_k |M_k × o| / (M_k · o)` (a positive
  combination of vectors in a convex cone stays in it), and
  `|M| >= min_k M_k · o / |o|`. Over a triangle's box, with `d` the normal
  at its centre, each normal is within `β_c + α_c` of `d` (`β_c` the angle
  between `d` and a cell's axis), so any two within twice the largest; and
  `|S_u × S_v| >= |M| / w_max^e`.
* *Rates.* Pointwise `|N_u| = |M × M_u| / |M|²` for any positive multiple
  `M` (the tangential part of `M_u` cancels in the cross product, and so
  does the derivative of the multiple), so on a cell
  `|N_u| <= max |M × M_u| / (min |M|)²` from the coefficients of the
  polynomial `M × M_u` (surface degrees at most 4, rational at most 2); a
  triangle's box takes the largest over its cells. For an edge the same
  holds with `N` a multiple of `C'` (`A'`, or `N_1`): the tangent turns
  over a cell's part `[τ_a, τ_b]` by at most
  `max |N × N'| (τ_b - τ_a) / μ²`, `μ` the enclosed `|N|` at the part's
  middle less `max |N'| (τ_b - τ_a)`, summed over the cells a segment
  meets.

A cell's multiple `M` depends on its parameter lengths (`L_u L_v` in its
own parameters) and on the weights; each cell's ratio is formed with its
own `M` before the largest is taken, so no two cells' multiples are
compared.

**Cusp corners.** At a corner `c` of the domain whose boundary rows leave
it in one direction, `S_u × S_v` vanishes: there `S_u` is a nonzero
multiple of the corner control's difference to its neighbour along `u`
(rational or not, `± p w_1 / (w_0 L_u)`), likewise `S_v`, so the test is
exact on the patch's controls. The unit normal's limit then depends on the
direction of approach (`Motor-c.brep` 378's four blends: 0.45 rad between
their two boundary rows), every bound above divides by a least `|M|` that
is zero on a box at `c`, and a boundary segment or triangle at the corner
never certifies. On the corner's patch and on its boxes `[c, c ± 2^-k L]`,
`k <= 3` (tighter third derivatives, until rounding takes over), with `y`
a box's own parameters from `c` and `M` the positive multiple above
(`A_s × A_t`, or `N_s × N_t`, one polynomial on the box, exactly zero at
`c`), Taylor's formula gives

    M(c + y) = G y + H[y, y] / 2 + E,  |E| <= T(y) / 6,

`G = [M_s M_t](c)`, `H` the second derivatives at `c` (enclosed, de
Casteljau at the corner), `T(y) = t_sss |y_s|³ + 3 t_sst y_s² |y_t| +
3 t_stt |y_s| y_t² + t_ttt |y_t|³` from the largest coefficients of the
third derivatives. `h(y) = G y + H[y, y] / 2` is quadratic, so over a
triangle `y_1 y_2 y_3` it lies in the convex hull of its quadratic Bézier
controls, the blossoms `h[y_i, y_j] = G (y_i + y_j) / 2 + H[y_i, y_j] / 2`.
A triangle at the corner is its rays `c + r d`, `0 < r <= 1`, `d` on the
opposite edge: `M / r = (1 - r) G d + r h(d) + E / r`, `|E / r| <=
r² T(d) / 6`, so the controls are the `G d_i` and the blossoms of the
opposite edge. Every normal is then within `asin ε` of a vector of the
controls' convex cone, `ε = T_max / (6 μ)` with `T_max` at the vertices'
largest `|y_s|`, `|y_t|` and `μ` the controls' least component along their
mean direction (a lower bound of the hull's norms). When the controls are
pairwise within 90°, the cone's directions are pairwise within the
widest angle `β` between two controls (for a fixed direction the
directions within `β <= 90°` of it form a convex cone, so the widest pair
is attained at controls), and the normal turns by at most
`β + 2 asin ε <= β + π ε` (`asin` convex on `[0, 1]`), the corner itself,
where no normal exists, excluded. A triangle's (or boundary segment's)
turn is the smallest of this and the bounds above, the corner turn only
for triangles inside one of the boxes. Along a row from the corner `β` is
zero and `ε` shrinks with the segment; a fan's triangles turn by their
sector's `β`.

**Counts.** An edge on a spline face needs, for each boundary segment with
parameter extents `Δu`, `Δv` between its pcurve's nodes,
`(a Δu² + 2 b Δu Δv + c Δv²)/8 <= 0.45 δ` and
`(M_u Δu + M_v Δv)/μ <= 0.45 θ` over the cells of the chord's box: T-a's
thin-triangle condition with local coefficients. A spline pcurve on an
analytic surface enters T-a's condition with the extents `H max |u'|` and
`H max |v'|` (`H` its parameter length, `D1` per coordinate over its
cells) and its `v` range from its control hull (a rational curve with
positive weights lies in the convex hull of its poles).
