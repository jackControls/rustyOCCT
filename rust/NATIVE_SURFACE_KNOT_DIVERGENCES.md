# Native surface knot editing comparisons

The bridge captures all 1,748 complete OCCT results before Rust execution.
The original pre-implementation capture records OCCT 7.9.3 and source-reference
revision `3d097a0328e71b826377d4814ab05ec3c3d23871` separately. It includes every
supported degree in each axis, complete unclamped raw support, independently
periodic U/V knots, seams, redundant directions and adapted original surface
GTest insertion/removal inputs.

Every Rust observation must match the independent exact coefficient equations,
including operation flags, both complete knot vectors, domains and the full
homogeneous control grid. Native observations contain Euclidean poles and
weights; comparisons use absolute `1e-10` plus relative `2e-12` per field.
Native removal is explicitly requested with tolerance `1e-9`. No tolerance is
applied to Rust's exact homogeneous removal contract.

OCCT 7.9.3 produces 1,602 direct matches. The 146 reviewed differences all
follow insertion/removal round trips for which the independent equations
restore the original complete grid exactly:

- 140 cases: OCCT rejects at least one inverse removal after floating-point
  insertion. Exact Rust controls recover the original representation.
- Six cases: OCCT accepts all operations but its rounded controls differ by up
  to approximately `4.18e-10`, exceeding the bridge budget while remaining below
  the explicitly requested native removal tolerance.

These are differences, not native matches. The review file pins the exact input
hash, OCCT version, every native control bit and knot, exact-result hash and
difference fields. A review never permits Rust to disagree with independent
mathematics. New runtime observations or changed values require a fresh review.

Linux OCCT 7.6.3 independently verifies the same 1,748 inputs, with 1,601
direct matches and 147 reviewed differences: 142 rejected inverse removals and
five accepted round trips with control errors outside the bridge budget. Every
reviewed exact result is again the original complete homogeneous grid. The
largest accepted native control error is approximately `4.18e-10`. Each runtime
has its own complete observations and review pins; their results are not merged
into a single match count.

Run `python3 rust/tools/compare_surface_knots.py --occt-root /path/to/occt`.
`--strict-native` ignores reviews and fails every native discrepancy.
Fixtures use full controls with one shared integer denominator, preserving
all rational values without substituting hashes or sampled point checks.
