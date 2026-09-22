# Native proximity observations

The 370 inputs were captured from installed OCCT 7.9.3 before the Rust proximity
implementation. The pinned source reference remains
`3d097a0328e71b826377d4814ab05ec3c3d23871`; the runtime is not built from it.
The corpus covers all 25 ordered point/line/segment/plane/triangle pairings,
three scales and three coordinate frames, plus geometric degeneracies and the
original BUC60870 edge/vertex input. Lines and planes are unbounded; segments
and triangles include their boundaries. A collapsed segment is represented by
a native vertex, the same singleton geometric set.

Local OCCT 7.9.3 observations:

| API | Distance matches | Nonresults | Numerical differences | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| `BRepExtrema_DistShapeShape` | 316 | 54 | 0 | 0 |
| `gp_Pnt`, `gp_Lin`, `gp_Pln` distances | 120 | 0 | 6 | 244 |

Linux OCCT 7.6.3 was captured in [CI run 35724489035](https://github.com/jackControls/rustyOCCT/actions/runs/35724489035).
It has the same 316 B-rep matches, 54 nonresults and 412 witness pairs, and all
126 applicable affine distances match. Its angle diagnostics select the
parallel branch for the six oblique cases above. The 370 complete Rust rational
results are byte-identical on Linux and macOS. The 54 Linux nonresults were
reviewed after capture, with separate version pins; the initial CI failure was
the required rejection of unreviewed observations, not a waived Rust failure.

These are two overlapping API observations of 370 distinct inputs, not 740
independent cases. All 412 returned B-rep witness pairs are retained and checked
for operand ownership and distance. Rust returns one deterministic minimum
pair, not the native enumeration of all reported extrema. Every Rust result
passes independent analytic Fraction formulas, minimal distance enclosures,
and an exact membership/supporting-plane certificate.

## Unbounded B-rep nonresults

`IsDone()` is false for six families at all nine scale/frame combinations:
the general plane/plane pair (`FF`), parallel lines, coincident lines, parallel
planes, coincident planes, and a parallel line/plane pair. The B-rep algorithm
collects extrema over vertex/edge/face strata; this API does not produce a
witness for those inputs. Its nonresults remain explicit. The affine API
supplies an additional distance observation, but does not turn a B-rep
nonresult into a match. Rust independently certifies all 54 minima.

## Oblique affine parallelism

In frame 2, at scales 0.125, 1 and 32, the parallel-line and parallel-plane
families have exact squared distance `189 * scale²`. For example, at scale 1,
the line direction is `(4,4,1)`, and the separation `(3,-6,12)` is perpendicular
to it. The plane normal is proportional to `(1,-2,4)`, with separation three
times that vector. Both distances are `sqrt(189)`, approximately
`13.74772708486752`.

The native line API instead returns `10.606601717798213`, and the plane API
returns zero. Both use `gp_Dir::IsParallel(..., gp::Resolution())` to select
the parallel branch. The 7.9.3 diagnostic records identical represented unit
directions on both operands, yet `gp_Dir::Angle` returns approximately
`8.5011e-18` for the lines and `4.4686e-18` for the planes. Both exceed
`gp::Resolution() = 2.2250738585072014e-308`, so parallelism is false. The line
routine follows its cross-direction branch; the plane routine returns zero
from its nonparallel branch. These observations explain this runtime's branch
choice; no claim is made about every OCCT version/compiler combination.

Rust retains the original defining points and decides rank exactly. The
independent cross-product distance formulas and global optimality certificates
agree on `189 * scale²`. The native budget remains
`1e-10 + 2e-12*abs(native)`, without relaxation.

## Review enforcement

[`occt-proximity-divergences.json`](fixtures/occt-proximity-divergences.json)
pins the runtime version, API, case, exact input hash, native status/value bits
and independently computed squared distance for each reviewed observation.
The bridge certifies Rust before consulting reviews. Changed observations,
versions, inputs or incorrect Rust results fail. Deliberate corruptions test
those guards; `--strict-native` fails on reviewed nonresults/differences too.
Runtime-specific reviews do not authorize observations from another version.

```sh
python3 rust/tools/compare_proximity.py --occt-root /path/to/occt
python3 rust/tools/compare_proximity.py --reuse-capture --strict-native
```

Inputs, both native result streams, all Rust rational witnesses, runtime
version, native angular diagnostics and the comparison report are saved under
`target/proximity-oracle/` and uploaded by CI. The BUC60870 geometry is reused
in native and Rust queries; this does not claim execution of the original
GoogleTest or add a passing unchanged DRAW test.
