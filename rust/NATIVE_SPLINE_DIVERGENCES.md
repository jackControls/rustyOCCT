# Reviewed native spline numerical differences

The native comparison budget remains `1e-10 + 2e-12*abs(native)`. This is our
test budget, not a claim that OCCT promises it for every degree and knot vector.
The independent exact tests and fuzz oracles have no numerical exception.

The new corpus was captured from OCCT before implementing periodic Rust curves
and tensor-product surfaces. On macOS arm64 with installed OCCT 7.9.3:

| Family | Observations | Matching cases | Reviewed divergent cases |
| --- | ---: | ---: | ---: |
| Curves | 456 | 427 | 29 |
| Surfaces | 210 | 169 | 41 |

All these new disagreements involve degree 25. Every Rust jet agrees with an
independent exact rational reference: Cox-de Boor basis recursion, its derivative
identity, tensor basis products, and closed quotient formulas. This reference
does not interpolate control points as the production implementation does.
Periodic knots are independently indexed through an infinite cyclic sequence.

A particularly direct check uses a clamped surface corner. Its exact position
is the final control point, regardless of the positive weights:

| Case | Exact corner / native LocalD0 | Position returned by native LocalD2 |
| --- | --- | --- |
| `d25_1_p00_1` | `(24, 0, -4)` | `(23.999992373135385, 0, -3.999996208109506)` |
| `d25_25_p00_1` | `(24, 25, 5)` | `(24.000001021255564, 25.000000518163041, 5.000020032463568)` |

Thus even two native evaluation orders can return different positions at a
corner. Exact basis evaluation also verifies each derivative, including Duv;
the review does not infer correct derivatives merely from the corner position.
The source path is recorded in `SOURCE_MAP.md`: curve D2 evaluates with
`BSplCLib::Bohm` then `PLib::RationalDerivative`; surface D2 uses the bivariate
evaluation and `BSplSLib::RationalDerivative`. The evidence establishes native
floating-evaluation discrepancies, not a uniquely localized internal defect.

`fixtures/occt-spline-jet-divergences.json` pins each reviewed case's family,
native version, input SHA-256, all native output bits, and a digest of the exact
rational jet (canonical comma-separated `numerator/denominator` entries).
`tools/review_spline_jets.py` recomputes that exact jet from the original input
and accepts Rust only if every component is one of its two minimal enclosure
endpoints. Changed native bits, inputs, versions, exact answers or incorrect
Rust results cannot use the review. Guards explicitly test these failures.

Reports count reviewed cases separately from matches; no high-degree cases
were removed and no comparison threshold was widened. Run either comparison
with `--strict-native` to fail on reviewed differences too. Full inputs, native
and Rust observations, runtime versions and reports are retained in CI artifacts
under `target/spline-oracle` and `target/surface-oracle`. The installed native
runtime is distinct from source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`.

The original nonperiodic OCCT 7.6.3 case remains independently pinned in
`occt-spline-divergences.json` and described in `VALIDATION.md`.
