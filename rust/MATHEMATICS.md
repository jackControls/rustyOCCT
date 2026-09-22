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
- 256 generated prisms over 25 scales, with concave/convex radial polygons and
  optional holes, check volume/first-moment conservation under splitting,
  introduced cut area, winding and offset reversal, rigid-motion covariance,
  nonnegative inertia quadratic forms, topology, sampled bounds and boundary
  classification. Each has 32 further classification probes before/after motion.
  Their seed and case index make failures reproducible.

These are bounded deterministic tests, not sustained coverage-guided fuzzing.
Finite test evidence supplements the arithmetic argument; it does not prove a
general kernel correct. The OCCT comparison corpus remains a separate oracle.

```sh
cargo test --workspace --locked
cargo test --workspace --locked --release
python3 rust/tools/generate_predicate_fixtures.py --check
```

Native CI executes both debug and optimized Rust tests. The rational fixture
generator is checked independently in the Linux validation job. Its default
write mode is for deliberate, reviewed corpus additions, never for making a
disagreement disappear.

## Next mathematical work

1. Extend exact predicates to 3D orientation and the distance/incircle decisions
   needed by the next algorithms. Specify units, domains and arithmetic bounds
   separately; do not reuse an arbitrary global epsilon.
2. Specify error and conditioning contracts for projections, curve/surface
   evaluation, root isolation and intersections. Use interval/error bounds or
   additional precision when ordinary arithmetic cannot establish the answer.
3. Strengthen topology invariants and tolerance propagation through each
   operation; sampled checks alone cannot certify full curve/surface agreement.
4. Add structured coverage-guided fuzzing, failure minimization, resource-limit
   tests and independent numerical references for these new contracts.

General certified constructions, exact 3D predicates, Booleans, generic topology
history, STEP and meshing remain unimplemented. Input restrictions and typed
errors remain necessary even with an exact orientation predicate.
