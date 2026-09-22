# Spline/plane native observations and exact contracts

The 132-case `GeomAPI_IntCS` corpus was captured on OCCT 7.9.3 before the Rust
intersector was implemented. It covers rational and polynomial Bézier curves,
clamped/unclamped/periodic B-splines, repeated contacts, contained spans, three
scales and two plane orientations. No native point is deleted or clustered by
the bridge. This is API differential testing, not additional unchanged DRAW
test coverage.

Every Rust result is independently checked: complete point/overlap counts,
contact type, one-sided orders, and minimal parameter/coordinate enclosures.
Python `Fraction` basis polynomials and SymPy 1.14.0 irreducible factorization
and continued-fraction isolation supply that answer. Production uses homogeneous
pole interpolation and integer Sturm/Sturm–Tarski sequences. Neither native
approximate roots nor Rust observations determine the oracle.

Local OCCT 7.9.3 has **97 matches and 35 reviewed differences**. A match uses
`1e-10 + 2e-12*abs(reference)`, the existing native comparison budget. This is
not an asserted OCCT guarantee. The registry pins version, input SHA-256, every
observed native value bit, and independent complete certificate SHA-256. A
review is accepted only after rechecking the complete exact Rust answer.
Changed or unreviewed observations fail. `--strict-native` fails on reviewed
differences too. Linux distribution observations need their own review.

## Reviewed families

| Family | Exact result / native observation |
| --- | --- |
| `factor_1_oblique_8.0` | One double root at 1/2. Native returns approximately 0.49999999745192014 and 0.500000000001005. Both observations are retained. |
| `factor_4_*` | Two double roots at 1/4 and 3/4. Native reports extra events near tangencies. `axis_1.0` returns four points, including three near 3/4. |
| `factor_7_*` (five reviewed cases) | One order-seven crossing at 1/2. Native parameter or position exceeds the budget; `axis_1.0` returns parameter 0.4999999966809783. |
| `rational_25_*` (five reviewed cases) | Exact basis substitution/root isolation confirms Rust. A native parameter or position exceeds the budget. In `axis_1.0`, the plane equation forces X=14.5 exactly; native gives 14.499999999638021. |
| `contained_*` | The numerator is identically zero on [0,1]. Rust returns that closed overlap; native reports neither points nor segments. |
| `partial_overlap_*` | The numerator vanishes identically on [1,2]. Rust returns an overlap with its endpoints included. Native reports endpoints 1 and 2 as isolated points and no segment. |
| `periodic_2_*` | Rust queries the **closed** domain [0,3], retaining seam parameter events 0 and 3. Native reports the start seam once. This is a declared parameter-domain convention difference, not a claimed OCCT defect. |

Factor-family control ordinates come from exact factored polynomials converted
to Bernstein form and scaled to integers before encoding. Their roots and
multiplicities survive binary64 input. The independent oracle nevertheless
recomputes every span from the actual encoded coordinates and weights.

## Reproduction

Install the modeling-algorithms SDK and the pinned test-only oracle:

```sh
python3 -m venv target/math-oracle-venv
target/math-oracle-venv/bin/python -m pip install --require-hashes -r rust/tools/math-oracle-requirements.txt
target/math-oracle-venv/bin/python rust/tools/generate_real_root_fixtures.py --check
target/math-oracle-venv/bin/python rust/tools/generate_spline_plane_fixtures.py --check
target/math-oracle-venv/bin/python rust/tools/compare_spline_plane.py
```

`target/spline-plane-oracle` retains inputs, native version, complete observations
and a report separating matches, reviewed differences and failures.
`--native-only` captures OCCT without building Rust. The production library
does not invoke C++, Python, SymPy or native OCCT.
