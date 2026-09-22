# Retained fuzz regressions

The runner seeds every `*.bin` under the matching target directory. Keep original
artifact bytes and names so campaign evidence remains traceable.

## Proximity: full-exponent triangle pairs

`proximity/slow-unit-1052e64729cba6e7060eca6bb910717b1678903a.bin` was saved during
corpus replay in [campaign 35724489063](https://github.com/jackControls/rustyOCCT/actions/runs/35724489063)
at `acfc37e7ef5044995b7cee997dac5a25baedc7a3`. The Linux sanitizer reported 14
seconds for the complete input and oracle, with no wrong answer or failure.
Its bytes are exactly the `random_TT_5` seed from `fixtures/proximity.tsv`:
two triangles with independently varying coordinate exponents, from subnormal
values through approximately `3.5e272`. That existing fixture independently
checks the exact minimum distance and global supporting-plane certificate.

Profiling showed repeated fraction reduction in candidate solves and point
evaluation. The solver now clears a common dyadic coordinate denominator,
uses integer-preserving Bareiss elimination and back substitution, and compares
homogeneous integer candidates before reducing only the chosen result. Native
corpus witnesses and all existing exact fixture expectations are unchanged.
Local five-run release replay medians, including both kernel and oracle, changed
from approximately 0.5445 seconds to 0.0836 seconds (6.5 times faster). These
are diagnostic timings on one host, not a general production latency guarantee.
The 20-second fuzz input limit is unchanged.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_proximity -- rust/fuzz/regressions/proximity/*.bin
```

## Spline intersections: expensive exact sign filtering

Both inputs came from [campaign 35715526133](https://github.com/jackControls/rustyOCCT/actions/runs/35715526133)
at `41b7b7d4f2fc153b0abb89855fd592f4fdac3ca3`. They were slow-unit
diagnostics, with no wrong answer, sanitizer failure, or timeout. The Linux
instrumented log reported 11 and 20 seconds respectively. These include the
mathematical oracle; they are not production latency measurements.

- `spline_intersections/slow-unit-782c568a8128da67116fb2dd30aacb47cd105d91.bin`
  encodes `(2*t^25,0,0)` on `[0,1]` against the unit sphere. There is exactly one
  crossing, at `t=2^(-1/25)`, position `(1,0,0)`. Ordinary Cargo coverage is
  `sphere_power_25_0` in `fixtures/spline-quadric.tsv`.
- `spline_intersections/slow-unit-b7d2cd2aaf90f34747c5d6472152631484e387ab.bin`
  encodes a rational degree-seven Bézier curve on `[-1,3]`, an oblique plane,
  and query `[2^-1074,2]`. Its local plane numerator is proportional to
  `t*(t-1/2)*(t-3/4)^4*(t-5/4)`, where `t=(u+1)/4`. The query contains a
  simple crossing at `u=1` and an order-four boundary contact at `u=2`.
  Ordinary Cargo coverage is `fuzz_subnormal_trim_degree7` in
  `fixtures/spline-plane.tsv`, independently generated with Python rational
  basis functions and continued-fraction roots.

Local release profiling put almost all replay time in the kernel, especially
fraction reduction in interval-Horner sign checks. Carrying integer numerators
over a shared positive denominator preserves the exact interval filter and its
Sturm–Tarski fallback. No geometric tolerance or expected result changed.

Replay the complete oracle without instrumentation or mutation from the repo root:

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_spline_intersections -- \
  rust/fuzz/regressions/spline_intersections/*.bin
```

Reported times include input construction, kernel work and all oracle checks.
They are diagnostic; correctness tests do not impose a machine-dependent
wall-clock threshold. The sanitizer campaigns retain their per-input timeout.

## Bézier editing: dense degree-25 periodic operations

`bezier_editing/periodic-degree25-split.bin` and
`bezier_editing/periodic-degree25-combined.bin` retain initial sanitizer slow
units `53ae219acb66d68bc0bbede17bec25d0b628f1cc` and
`a8115407babd486fe26e65a5ad9e43473fee7fb7`. They are original seeds, not minimized
crashes. Both passed the complete mathematical oracle, but initially took
11 and 16 seconds under AddressSanitizer. They exercise rational periodic
degree-25 extraction across two periods, splits, and combined
trim/elevate/reverse/split checks. The checker now clears common coefficient
denominators before binomial power transforms, rather than repeatedly reducing
fractions in every multiply/add. Geometry, assertions and timeout limits are
unchanged. Runtime includes independent oracle work, not just kernel work.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_bezier_editing -- rust/fuzz/regressions/bezier_editing/*.bin
```

## Polynomial roots: a wide isolator and a large query

`roots/slow-unit-0d9082d327a81ddd2e03cf8963cb57f12c9bf198.bin` came from
[campaign 35718215174](https://github.com/jackControls/rustyOCCT/actions/runs/35718215174)
at `61f2ab507cf6aa5ef8b0834834d25df20104e365`. Its 16-second instrumented
execution passed the oracle. The degree-25 defining polynomial is

```text
-(t+2)^11*(t+1)*t*(t-1)^2*(t-2)^2*(t²-2)*(t²-3)*(t²-5)*(t²+1).
```

It queries a degree-14 polynomial whose exact binary64 coefficients appear in
`fuzz_wide_isolator_query` in `fixtures/real-roots.tsv`. Profiling identified the
kernel's sign query at `t=-1`, with a wide initial isolator, as the expensive
step. This behavior also reproduced before the shared-denominator optimization.
Up to 64 exact bisections now give the interval filter another chance to decide
the sign; a still-undecided sign takes the complete Sturm–Tarski path. The
independent fixture verifies all eleven distinct roots, their multiplicities,
minimal bounds and query signs.

Replay the entire oracle with:

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_roots -- rust/fuzz/regressions/roots/*.bin
```
