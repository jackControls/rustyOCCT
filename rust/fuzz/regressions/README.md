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

At `1fcee25d`, the Linux Bézier campaign exhausted its 600-second startup
allowance replaying 381 retained seeds in
[run 35835064613](https://github.com/jackControls/rustyOCCT/actions/runs/35835064613).
It failed before initialization or mutation; the absence of a crash did not
make this a successful campaign. Profiling showed the independent Cox
polynomial checker dominated the retained periodic inputs. It now carries one
positive denominator per basis polynomial and skips identity substitutions.
The complete combined replay fell from about 0.30 to 0.16 seconds locally.
Every coefficient, derivative, enclosure and operation check remains. Linux
release CI now also verifies this second oracle against all 636 independently
generated Python fixtures.

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

`surface_editing/slow-unit-3e01ced2e0056a1159d07d2a69d7405e6dddf050.bin`
retains the complete 802-byte input reported at 21 seconds by
[run 35835064613](https://github.com/jackControls/rustyOCCT/actions/runs/35835064613).
That job passed, but this input exceeded the nominal 20-second target; the
libFuzzer alarm is not a hard real-time bound. It is a degree-two doubly
periodic surface with full-exponent binary64 controls and positive weights,
queried over three periods in both directions. This is an unminimized slow
input, not a mathematical disagreement. The ordinary
`fuzz_full_exponent_periodic` fixture independently recomputes all resulting
controls, and its complete second coefficient/jet oracle always runs.

Closed-form quotient derivatives in the independent tensor checker now clear
one common denominator and use integer arithmetic until the final scalar.
Their formulas remain separate from the kernel's recursive quotient rule.
Complete local replay fell from about 0.60 to 0.34 seconds. All 745 fixtures
are checked by both independent oracles; these diagnostic times are not
production bounds, and no assertion or campaign limit was relaxed.

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
second oracle. The next Linux campaign at `7b0fc6c6` completed corpus replay
and its full mutation budget, with a maximum reported input time of 13 seconds.
The three original removal cases no longer emitted slow-unit diagnostics.

That campaign retained a different 432-byte input as
`periodic-d25-refinement.bin` (`a363a30e7ba20a77d76980d6f4cddd2c7c71343a`),
from [CI run 35824211538](https://github.com/jackControls/rustyOCCT/actions/runs/35824211538).
It refines a degree-25 periodic rational curve with 28 poles and seam
multiplicity two, inserting `85/257` and `171/257` to multiplicity 25.
All full-support coefficient, extraction, jet, enclosure and periodic wrapping
checks passed. Local complete release replay took 0.433 seconds, including
0.213 seconds of edits and 0.196 seconds of extraction/jet/bound checks;
the sanitizer replay took 4.359 seconds. The Linux 13-second observation is
below the unchanged 20-second limit, with no failure or mathematical
disagreement. The full original input is preserved and the ordinary
`retained_periodic_degree25_refinement_preserves_complete_polynomial` test
replays its operations against the independent complete polynomial oracle.
These timings remain host-specific diagnostics, not production guarantees.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_knot_editing -- rust/fuzz/regressions/knot_editing/*.bin
```

`knot_editing/slow-unit-66a9ef4c53624190feb8541eef54ef38fde65a8b.bin`
retains the full 12-byte mutation reported at 14 seconds by the same
[run 35835064613](https://github.com/jackControls/rustyOCCT/actions/runs/35835064613).
It removes the origin of a degree-25 constant periodic curve with the smallest
positive binary64 knot spacing. Complete local release replay passed in about
0.52 seconds. The ordinary
`retained_subnormal_periodic_degree25_removal_preserves_complete_polynomial`
test verifies the independent removal equations, full-support polynomial
identity, exact derivatives and finite enclosures. The input was not minimized,
and fixed-input replay is separate from mutation coverage.

## Exact edited-curve intersections: high-degree weighted coordinates

The first local `exact_spline_intersections` campaign, in the working tree based
on `d2350d6e`, completed its 60-second mutation budget with 93 mutations and no
mathematical or sanitizer failures. It retained two unminimized 112-byte inputs:

- `exact_spline_intersections/degree25-negative-huge-domain.bin`
  (`5f281664175a36f607cb6443d18e1c23354d5004`), reported at 15 seconds.
- `exact_spline_intersections/degree25-positive-huge-domain.bin`
  (`0afe9b337cc4748968cb489a52eeeed5cd49aa57`), reported at 17 seconds.

Both represent degree-25 rational curves with independently known factored
cylinder contacts, arbitrary positive Bernstein weights and the cylinder axis
permuted to Y. They refine at fractions `85/257` and `171/257`, then remove the
first cut exactly. Their parameter domains are respectively
`[-2^2048, -2^2048 + 7/3]` and `[2^1024, 2^1025]`. Exact contacts remain usable
although finite parameter enclosures are impossible.

Complete release replay took about 1.80/1.77 seconds, mostly in coordinate sign
queries and minimal enclosures. The independent full polynomial identity check
took only about 0.01 seconds. The kernel now clears one shared denominator for
all homogeneous coordinates per span and reduces polynomial sign queries by
positive pseudo-remainders modulo the root's square-free defining polynomial.
These changes preserve signs, exact zeros and the complete algebraic fallback.
Complete release replay fell to about 0.88/0.83 seconds; individual sanitizer
replay took 7.17/7.20 seconds. These are host-specific diagnostics, not a kernel
latency guarantee, and fixed-input replay is not a fuzz campaign.

Both full inputs are added to every campaign's corpus. The ordinary
`retained_degree25_weighted_contacts_and_parameter_ranges` test independently
reconstructs their factors and weights, proves the complete edited function
unchanged, and verifies every contact, order, coordinate and minimal enclosure.
No degree, operand-size, assertion, timeout or RSS limit was relaxed.

The first Linux campaign at `c197a8cc` exhausted its 660-second outer deadline
while replaying the 330 seeds, so it failed without claiming any mutation
coverage. It also saved the 21-second `degree25-unit-domain.bin`
(`69f8157e36d0c5b2d9b386e589daa359010b04d6`) from
[CI run 35830708600](https://github.com/jackControls/rustyOCCT/actions/runs/35830708600).
This is the same weighted cylinder recipe on `[0,1]`. The local 600-second
campaign at that revision completed 1,160 mutations, but saved a separate
10-second plane input, `mutated-rational-contacts.bin`
(`a20f16ed32d2aa6ac05bb12e663766a033f9dafe`), whose mutated factors and weights
are retained unchanged. The ordinary test covers both additional full recipes.

Root refinement now recognizes rational candidates only after proving their
membership in the isolating interval and an exact zero of the defining
polynomial. Failed candidates retain the complete algebraic fallback. Plane
intersections perform up to 32 initial bisections with these checks; quadrics
retain their 128-step maximum. Complete release replay of the Linux input
improved from 0.83 to 0.43 seconds, and the mutated plane from 0.76 to 0.25
seconds. Local sanitizer replays took 4.35 and 2.50 seconds respectively; these
fixed-input runs are separate from mutation campaigns. Their assertions and
the campaign's corpus, time and memory limits remain intact. These measurements
still do not imply production latency bounds.

At `fd8fca65`, Linux finished startup in 593.19 seconds and completed 60.09
seconds of mutation, but the combined 660-second deadline killed the process
before its last input and final statistics finished. This remained a failed
campaign, with no invented mutation count, in
[run 35832787056](https://github.com/jackControls/rustyOCCT/actions/runs/35832787056).
Its 13-second degree-25 sphere seed on `[0,1]` is retained as
`degree25-sphere-unit-domain.bin` (`c8011cfcc8d4f68a3f2e52dc62433fe566de6813`)
and included in the same ordinary complete-contact test. The local campaign
at that revision completed 1,202 mutations over 600 seconds without a new slow
input. These local results do not replace the failed Linux result.

Further profiling separated the remaining cost: about 0.26 seconds in exact
de Boor span polynomials and 0.10 seconds in the implicit equation per retained
cylinder/sphere input. Shared positive denominators now avoid repeated rational
reductions during de Boor coefficient interpolation. The five retained full
release replays take about 0.07–0.22 seconds with the same complete oracles.

The runner also now enforces startup, mutation and shutdown independently:
600 seconds for build/replay, the full requested mutation duration, then a
25-second shutdown grace for the unchanged 20-second input timeout and final
reporting. Tests reproduce the observed 593.19-second boundary, reject late or
missing initialization, allow a final input to finish, and kill an ignored stop
request. This changes shutdown accounting; it does not extend the startup,
per-input, memory or geometry limits, or accept incomplete reports.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_exact_spline_intersections -- \
  rust/fuzz/regressions/exact_spline_intersections/*.bin
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


## Surface knots: degree-25 periodic origin removal

`surface_knots/periodic-degree25-origin-removal.bin` is the 12-byte structural
reproducer for an initial surface-knot fuzz startup timeout. Removing the four
trailing zero bytes from the saved 16-byte input leaves the decoded geometry
unchanged. It constructs a 28 by 28, degree-25 periodic surface, redundant in U
with nonconstant weighted V geometry, then deletes the U origin and checks the
entire result, patch coefficients, isocurves and wrapped exact partials.

The first complete release replay took about 7.2 seconds. Shared Cox matrices
and elimination reduced checker repetition; reusable differentiated de Boor
maps reduced kernel jet work. The same full release replay then took about
0.77 seconds on the development machine. These are diagnostic measurements,
not portable performance guarantees. Ordinary Cargo tests preserve the same
complete geometry and independent checks. No homogeneous coefficient or residual
check was dropped.

```sh
cargo run --manifest-path rust/fuzz/Cargo.toml --locked --release \
  --example replay_surface_knots -- rust/fuzz/regressions/surface_knots/*.bin
```


`surface_knots/unclamped-degree25-both-refinement.bin` retains the full 28 by
28 unclamped tensor that exceeded the initial 20-second instrumented budget.
Its 3,152 bytes omit only unused trailing data from the original 3,216-byte
input. The additional `unclamped-degree25-full-multiplicities.bin` raises both
new knot multiplicities to 25. Shared sparse Boehm maps and integer affine
blends remove repeated kernel reductions; common Cox matrices and exact integer
content cancellation preserve every independent raw-support equation.
Both complete cases are ordinary Cargo regressions. The full-multiplicity case
still takes about 23 seconds with instrumentation, so this new target has an
explicit 60-second per-input verification budget. Existing targets retain their
20-second budget, and all targets keep the 2 GiB process memory limit.

`degree23-input-at-campaign-memory-limit.bin` preserves the 432 used bytes of
the input active when the first long campaign crossed 2 GiB after 419 mutations.
That input passes on its own; it is not a standalone OOM reproducer. The failed
campaign's manifest, full log and corpus remain in the local verification
artifacts. A new sustained campaign must pass the unchanged memory gate before
this surface-knot milestone is considered verified.

A subsequent diagnostic replay of all 197 saved inputs reached 2,205 MiB RSS,
while sanitizer accounting after completed inputs remained near 25 MB of live
allocations. The runtime's usual allocator purge runs only during mutation,
after seed replay. The test-only surface harness now also purges freed allocator
memory between completed inputs, at most once per second. Both the failed
2 GiB campaigns and the separate 4 GiB diagnostic replay are retained; the
diagnostic is not a successful gating campaign. The fresh campaign still uses
the 2 GiB memory gate and full mutation budget.

Allocator cleanup alone did not fix macOS replay: it still exceeded the cap
with approximately 38 MB live and 247 MB quarantined, plus unused allocator
chunks. Ordinary release replay of the identical 197 inputs twice passed with
64,700,416 bytes peak RSS. The target now explicitly budgets its sanitizer
quarantine at 64 MiB; the shorter freed-memory detection window and reproduction
settings are documented in [FUZZING.md](../../FUZZING.md). The resulting local
campaign completed 600.01 seconds of mutation and 337 mutation executions after
replaying all 197 saved inputs, with a 1,018 MiB peak RSS and a slowest input of
25 seconds. The failed earlier runs remain evidence.

Linux run `35870767498` completed mutation but reported the retained full-
multiplicity input at 64 seconds against its 60-second alarm. The runner now
checks final per-input and RSS statistics even after a zero exit. Complete
before/after control-column pairs share polynomial projections only when every
integer in both columns is identical; deliberate corruption of every field is
still rejected. This reduced the local full-input replay from approximately
3.0 to 2.6 seconds without removing coefficient equations.

The same run exhausted the original fixed 600-second build/replay allowance in
five older targets as retained corpora grew. Startup now scales with the input
count, separately capped at 3,600 seconds, as specified in FUZZING.md. Every
saved input remains in replay, and requested mutation time and per-input limits
are unchanged. These runner changes require fresh complete CI campaigns.
