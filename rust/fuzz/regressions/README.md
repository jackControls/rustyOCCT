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

The combined-edit seed then exceeded the unchanged 20-second limit on Linux
in [campaign 35740265339](https://github.com/jackControls/rustyOCCT/actions/runs/35740265339)
at `69c79163a5ef21ecf17a0f9e448d7eeb6849a08c` (the alarm fired at 27 seconds).
`periodic-degree25-trim.bin` and `unclamped-degree25-trim.bin` retain the other
reported Linux slow units; `mutated-degree25.bin` preserves the new slow unit
`7da76f940ea56d81abfc066fa5fba8ff44517a03` from the local ten-minute campaign.
Stage timing isolated substantial kernel subdivision/evaluation cost. Those
de Casteljau recurrences and elevation now use shared-denominator integer rows,
normalizing only returned values. The same combined input's local optimized
kernel-plus-oracle replay fell from 0.912 to 0.332 seconds. These measurements
are machine-specific diagnostics; no checks, input families or limits were
removed. The replay example prints the separate stage timings.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_bezier_editing -- rust/fuzz/regressions/bezier_editing/*.bin
```

## Tensor patches: degree-25 unclamped composed edits

`surface_editing/unclamped-degree25-combined.bin` is the complete 2,816-byte
seed saved as `slow-unit-a76a8bdec26f2a5742d5fb75e5865b37aaa9b452` by the first
local surface-editing campaign. It took 13 instrumented seconds and passed all
coefficient, jet, enclosure and commutation assertions. This is a retained
performance case, not a minimized crash. Both directions have degree 25 and
unclamped knots; the sequence restricts, reverses, exchanges axes and splits.

Initial profiling showed repeated spline blossoming on every row dominated
extraction. Reusing exact blossom maps reduced extraction from roughly 2.88 to
0.31 seconds locally. Sharing one de Casteljau triangle for position/D1/D2
further reduced jet work. The complete uninstrumented kernel-plus-oracle replay
fell from about 4.02 to 1.23 seconds; machine-specific timings are diagnostic,
not a latency guarantee. The 20-second input limit and all assertions remain.

The first Linux campaign at `253055ba` completed its full mutation budget, but
reported the same combined seed at **25 seconds** during replay and
`unclamped-degree25-extraction.bin` (`b91792ae56462056aaea07e89b28c271d87083b9`)
at 16 seconds. A successful run did not establish adequate time margin:
libFuzzer's alarm and slow-unit timing are not a hard real-time deadline.
The local ten-minute campaign also retained the 1,577-byte mutation
`mutated-unclamped-degree25.bin` (`151a86ed4e5e91dffeca5886c9bb2463933c0e98`)
at 10 instrumented seconds. Both retained inputs pass their complete oracles.

Exact local knot insertion now builds the identical extraction map with fewer
stages than repeated blossom queries. The independent polynomial checker also
reuses denominator-cleared Cox matrices and skips identity substitutions.
These changes retain every coefficient and assertion. Local full replays of
the combined, extraction and mutated cases take approximately 0.75, 0.37 and
0.80 seconds; the combined case's extraction stage is about 0.073 seconds.

The second Linux run at `7b32a108` did time out on the original combined input
after a 21-second alarm, before corpus replay or mutation completed. It is an
explicit failed campaign, not a successful fuzz run. It also retained the
clamped (`f9f28c9651fe1344f65fef3d35b6c9887e1fb641`, 17 seconds) and doubly
periodic (`a05dfa063c3e189f31e39b03f04374de608fdfd1`, 13 seconds) combined seeds.
All five full inputs are retained, not described as minimized crashes. The
independent surface fixture generator decodes their geometry and recomputes
every output control, including the previously missing degree-25 unclamped
ordinary regression. The source test names begin with `fuzz_`.

The third Linux campaign at `73af204c` passed, including its complete mutation
budget, but the retained mutation `151a86ed...` still took 20 instrumented
seconds. The complete local ten-minute campaign passed 1,285 mutation executions
without a new slow artifact. A green job alone did not resolve the timing risk.
Detailed profiling separated kernel work, complete-coefficient conversion,
independent jets and enclosure checks. The checker now carries one positive
integer denominator across tensor polynomial transforms and reuses exact
binomial/monomial matrices; it reduces only returned scalar values. Complete
coefficient equality still checks every component by exact cross multiplication.
No assertion, input family, timeout or memory limit is removed. On this local
host the retained mutation's uninstrumented replay fell from about 0.80 to 0.49
seconds. The exhaustive second-oracle pass checks all 744 fixtures, and is also
run by Linux release CI. Timings remain diagnostics rather than guarantees.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_surface_editing -- rust/fuzz/regressions/surface_editing/*.bin
```

## Exact knot editing: degree-25 removal equations

`knot_editing/constant-periodic-d25-removal.bin` retains the complete 16-byte
input from local artifact `timeout-511950b9670c399e097346e15a935fa51c700b3b`.
It was a seed-replay timeout, not a geometry disagreement or a minimized crash.
The ordinary `constant_periodic_d25_remove_existing` fixture and dedicated
`retained_degree25_removal_exercises_complete_integer_equation_solver` test
cover the same full coefficient system.

The initial sanitizer run exceeded the unchanged 20-second per-input limit.
Release profiling measured 6.61 seconds, including 4.38 seconds solving the
independent removal equations and 2.14 seconds checking full polynomial identity;
kernel removal itself took approximately 0.0013 seconds. The checker now keeps
Cox polynomials over common integer denominators and eliminates primitive
integer equations. Every coefficient and rejection assertion remains. On the
same local host, complete release replay fell to 0.57 seconds. The next campaign
completed its full seed replay and 60-second mutation budget. These diagnostic
timings do not establish a production latency guarantee.

The first Linux knot-editing campaign at `08c97766` passed its complete mutation
budget but reported 15- and 22-second seeds. The full inputs are retained as
`constant-periodic-d25-origin-removal.bin` (`a21ad3d7...`) and
`constant-periodic-d25-seam-roundtrip.bin` (`cd58aee8...`), from
[CI run 35822554768](https://github.com/jackControls/rustyOCCT/actions/runs/35822554768).
They also map to complete ordinary fixtures and the dedicated regression test.
These are unminimized performance cases; the green job did not remove their
latency risk.

The checker now uses fraction-free elimination with asserted exact divisions.
After obtaining full rank, it checks every remaining coefficient equation by
exact substitution instead of eliminating dependent equations again. Negative
weights and failed-removal consistency checks remain unchanged. Complete local
release replays improved from about 1.00/0.62 seconds to 0.24/0.31 seconds for
the origin/seam cases; every one of the 723 fixtures also passes the revised
second oracle. Linux sanitizer timing remains a separate validation gate.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_knot_editing -- rust/fuzz/regressions/knot_editing/*.bin
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
