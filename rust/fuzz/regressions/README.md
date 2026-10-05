# Retained fuzz regressions

The runner seeds every `*.bin` under the matching target directory. Keep original
artifact bytes and names so campaign evidence remains traceable.

## Tessellation: the first spline campaign's slow units (T-b)

`tessellation/slow-unit-5f6821b6d0e321e663b20155e303bb9db0e1514c.bin` (a
spline prism) and `tessellation/slow-unit-ce04689a160ed2848b526283b48e2b2c51c7c209.bin`
(a spline sheet) were saved by the first local 600-second campaign after T-b
added spline bodies to the target, as slow-unit diagnostics (the slowest
input 19 seconds on a loaded host) with no wrong answer, sanitizer failure or
timeout. Local AddressSanitizer replay then took about 9 and 4 seconds per
execution: the prism's validation, its exact certified mass enclosure, two
tessellations and a projection oracle with the kernel's exact evaluation.
The target now checks a spline prism's volume against a quadrature of its
profile, re-tessellates only meshes of at most 500 triangles, projects with
its own binary64 Cox-de Boor evaluation, builds prisms on axis-aligned
frames, and draws spline deflections from a 64th of the size and angles
from 0.35 rad; the two inputs (their bytes now decode to other bodies) replay
in about 3.3 and 2.1 seconds under the sanitizer on a host at load 10,
within the `20 s / 2.6` triage rule. The 20-second input limit is unchanged.

## B-rep validation: a far, small spline prism

`brep_validation/slow-unit-47bdac9ab6a50070d76f4f45029d0419a26705eb.bin` was
saved during corpus replay in [campaign 36302492715](https://github.com/jackControls/rustyOCCT/actions/runs/36302492715)
at `ac92ec89`, a slow-unit diagnostic (14 seconds on Linux) with no wrong
answer, sanitizer failure or timeout. It is mutation 32 on the cubic bulge
(S4), moved by `p -> 2^-10 p + (-7.5625, 4.5, 0.0625)`, then a pcurve shifted
within the resolution: four validations (the determinism check, the
constructor, the mutation's) and the certified mass properties. Local
AddressSanitizer replay took about 4.6 seconds per execution, most of it in
the exact deviation composition of the moved spline uses (wide exact
rationals from the scale and the offset), within the `20 s / 2.6` triage
rule; reusing the report's constructed topology for the mass check brought it
to about 4.1 seconds. The 20-second input limit is unchanged.

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

Linux run `35880483859` retained a new 3,152-byte control-grid mutation as
`surface_knots/unclamped-degree25-mutated-grid.bin`
(`6390006ed8fa76be9aabd74d9dd9c8201b9ededb`). The runtime timed it out after 68
seconds against the unchanged 60-second limit. It passes all mathematical
checks in ordinary replay and is now included in the degree-25 Cargo regression.
The checker proves full polynomial equality by continuity induction along the
entire common raw-support partition: at a knot of maximum multiplicity m,
equality on the previous span forces the first p+1-m coefficients of the
degree-p difference to vanish. It explicitly checks the remaining coefficients
of every field. Nonperiodic sums start at zero outside their raw support;
periodic sums require a complete first-cell comparison. An adversarial test
compares this certificate with every explicit power coefficient, including
corrupted inactive controls, weights and constant shifts at periodic seams.
The complete saved-input replay decreased from 2.60 to 1.41 seconds locally;
instrumented replay took 26.35 seconds. These are diagnostic measurements,
and fresh Linux campaigns must still satisfy the 60-second input gate.

The same run killed the curve-knot campaign during shutdown without final
statistics. Inspection of the pinned libFuzzer `FuzzerLoop.cpp` showed that
`stop_file` is checked outside `MutateAndTestOne`, which can execute five
callbacks. The runner now explicitly pins that existing mutation depth and
reserves the full last batch: 105 seconds for 20-second inputs, 305 for
60-second inputs. This preserves input limits and mutation depth; it does not
turn the failed run into a pass. The original logs and corpus are retained.

Before these checker changes, clean revision `b413b83e` completed a second
local campaign with 600.10 seconds of mutation, 267 mutation executions after
296 retained inputs, a 1,021 MiB peak RSS and a 25-second slowest input. Linux's
resource failures remain separate from that successful local evidence.

At `f24b486a`, all 15 Linux campaigns passed, including the retained surface
timeout input (41-second slowest input, 440 MiB peak RSS for that campaign).
The larger local surface corpus had grown to 385 inputs and exhausted the
890-second startup estimate before mutation began. Its failed manifest and
log remain retained. The runner now calculates startup from the actual
20/60-second per-input allowance plus build time, with a one-hour absolute
cap, instead of assuming two seconds per saved input. It still requires a
complete separate mutation budget and positive mutation executions. No input
is discarded and no individual time or memory limit is raised.

## Degree-elevation rejection contract

The first degree-elevation seed replay stopped on
`crash-0e91d858c2385d16d407340416d35a9770047460`, before mutation started.
The test harness incorrectly expected `InvalidCurve` for degree reduction;
the public knot-axis contract returns `InvalidSpline` for invalid degree and
`LimitExceeded` for an oversized result. The seven-byte
`degree_elevation/degree-reduction-error.bin` retains the minimized trigger.
The checker now distinguishes both typed errors and still requires rejection.
The failed campaign is preserved as a failure, not counted as fuzzing time.

A separate dense degree-24-to-25 unclamped tensor replay with coordinates
scaled by 2^4096 took 100.985 seconds under ASan. Single-file libFuzzer replay
does not enforce its normal per-input alarm; the enclosing 95-second probe
was recorded as failed, and process inspection confirmed no leftover replay.
`degree_elevation/extreme-unclamped-both-24-25.bin` retains all 2,932 relevant
bytes (the unused input suffix was removed). Primitive elimination reduced
integer growth; reconstructing every source unit column before applying the
map to all transverse fields removed coordinate magnitude from factorization.
All unit-column equations and all resulting controls remain checked exactly.
The fixed instrumented replay and full campaigns must pass their original
60-second/2-GiB gates before acceptance.

Linux runs [35897946980](https://github.com/jackControls/rustyOCCT/actions/runs/35897946980)
and [35899308823](https://github.com/jackControls/rustyOCCT/actions/runs/35899308823)
both timed out on the identical
`timeout-1a8f85cf535d9243909bccfd017a763e7fa28bec` during seed replay. It raises
a 27-by-27 uniform periodic surface from degree 24 to 25 in both directions,
then checks the independent full grid, transposition and exact jets.
`degree_elevation/uniform-periodic-both-24-25.bin` retains the 2,932 consumed
bytes of the 3,216-byte input; only its unused suffix was removed.

Profiling identified repeated production control-map construction. The two
identical axes now share their exact map, and periodic degree elevation uses
the three periods needed by the complete canonical support instead of five.
The proof and independent unit-column boundary tests are described in
[`DEGREE_ELEVATION.md`](../../DEGREE_ELEVATION.md). Every resulting control and
all composition/derivative assertions remain enabled. On the local machine,
ordinary replay fell from 4.253 to 2.074 seconds; the ASan replay took 36.086
seconds. Single-file replay is timed separately because libFuzzer does not
enforce its usual alarm in that mode. Published Linux campaigns must still
pass the unchanged 60-second/2-GiB input gates.

## Degree elevation: a degree-19 periodic axis raised to 25

`degree_elevation/timeout-559a48007e6dd6414eb28692d5545d662d920ec1.bin` is
the original 21-byte artifact of the scheduled
[run 36868817257](https://github.com/jackControls/rustyOCCT/actions/runs/36868817257)
(87 s under AddressSanitizer against the 60 s limit). It is a surface: U a
uniform periodic degree-19 axis of 23 simple knots scaled by `2^4096`, V the
clamped degree-2 axis `0, 1, 4`, small integer controls, raised to degrees
25 and 2, then staged (an identity, then the same elevation again). The
scale is immaterial (the same input at scale 1 takes as long): nearly all of
the time was the production control map, Prautzsch's rank curves over three
periods, six degree steps of about 21 ranks each inserting about 63 knots,
every blend a `BigRational` product and sum reduced by `num_integer`'s
binary gcd. Every coefficient fits 64 bits there. The map is now computed in
machine words first, the knots moved by one positive affine map onto
integers (the insertion ratios are invariant under it), each sum and product
reduced as `BigRational`'s, and any overflow computes it again with
`BigRational` through `crate::rational`'s faster gcd; knot refinement's
blends use that too. The map, and so every result, is the same exact map:
the target's independent reconstruction, composition and jets pass
unchanged, and the kernel's tests compare the word and rational maps.

| Input | Release replay | AddressSanitizer (loaded Mac) |
| --- | ---: | ---: |
| `timeout-559a…` | 9.6 s → 1.1 s (map 4.7 s → 0.36 s) | 262 s → 14 s (733 → 28 Gcycles) |
| `slow-unit-25ac…` (55 s on Linux) | 4.0 s → 1.7 s | 62 s → 39 s (168 → 103 Gcycles) |
| corpus `41cb4b82…` | 5.6 s → 3.1 s | 100 s → 35 s (190 → 97 Gcycles) |
| corpus `8db50da3…` | | 48 s → 34 s (131 → 93 Gcycles) |

The run's other slow units now spend most of their time in the independent
oracle's Cox systems, which this change does not touch. With debug
assertions the target's 278 corpus inputs and 4 regressions replay without a
failure, the slowest in 3.4 s. The 60-second limit is unchanged.

## Spline/linear: sub-float known-factor root pairs

`spline_linear/slow-replay-9ff24f8f748d219d6d9db4edc3ee40340a02bad0.bin` and
`spline_linear/slow-replay-42dacee5663eef4c37bfd32d84bf2b9c29617986.bin` came
from a deterministic development replay of 3,000 random inputs, before any
libFuzzer campaign. Both are known-factor curves elevated to degree 24 or 25
with a sub-float root pair near `1/3`. Their complete release checks took 65.5
and 29.9 seconds, with no wrong answer. Profiling found schoolbook integer gcds
during polynomial content removal and late recognition of a large-denominator
rational root. After the binary gcd and rational-root-theorem candidate, they
take 0.79 and 0.60 seconds in release. The local AddressSanitizer target takes
6.7 and 5.0 seconds. Limits and assertions are unchanged.

`spline_linear/timeout-b85e667578fb887a586777f84b229f4bca17e25e.bin` is the
original artifact from the first clean local 600-second campaign at `da7be8be`.
It is a degree-25 known-factor line query. After 143.65 seconds of replay and
462 mutation executions, the sanitizer callback exceeded 20 seconds; the
release check then took 4.48 seconds. Every exact zero proof at one rational
parameter built a gcd with a 13,378-bit query polynomial. That root's degree-10
defining polynomial has a 615-bit leading coefficient. The fixed 64-step
refinement stopped at isolator width about `2^-604`, just short of the `2^-615`
needed by the rational-root-theorem candidate. A capped extra bisection to that
bound (at most 64 steps) recognizes the root before any gcd. The release check
now takes 1.25 seconds and the local sanitizer target 11.7 seconds.

## B-rep interop: a location index beyond the table

`brep_io/crash-9557e3ae7d308500f8f03642dc84fa23cc1e7438.bin` was found by the
first local smoke campaign of `brep_io`, before any commit of the target. A
mutated prism text referenced location 1 while the location table was empty.
The reader range-checked subshape locations but not those of edge
representations and faces, so the converter indexed the table out of bounds. The reader now range-checks every curve, pcurve, surface and
location reference and returns `BrepError::Reference`.

## B-rep interop: the harness's integer nudge overflowed

`brep_io/crash-faf6616697c9ef66c6d010852377229eb8a0c82c.bin` is the original
75-byte artifact of the scheduled
[run 36868817257](https://github.com/jackControls/rustyOCCT/actions/runs/36868817257):
`attempt to add with overflow` at `src/brep_io.rs:100`, in the target's own
`mutate`, before any text reached the reader. The mutation that turns a
number into an integer casts it with `as i64`, which saturates a huge or
infinite number at the end of `i64`'s range, and then added a nudge in
`-2..=2` under overflow checks. No kernel value is involved. The nudge now
saturates too (`saturating_add`), the same integer token for every other
number and the same bytes consumed, so every check on the mutated text is
unchanged; the input replays in 0.02 s with debug assertions.

## Surface editing: a completing input past the 20-second limit

`surface_editing/timeout-5cf6cbdd8250eeb22f43e0502d33c014522c1a2b.bin` and
`surface_editing/slow-unit-f240e0e8b3821b28ef9edb4177bac843e1213cdd.bin` come
from the push run at `3ba8c1f2`
([run 36236430909](https://github.com/jackControls/rustyOCCT/actions/runs/36236430909)).
Neither is a checked-in seed; both came from the evolving CI corpus. On the
Linux runner under AddressSanitizer the first exceeded the 20-second
per-input limit and the second took 14 seconds. Both complete without a
failure or disagreement. Locally on the development Mac:

| Input | Without sanitizer | AddressSanitizer (three runs) |
| --- | ---: | ---: |
| `timeout-5cf6…` | 4.06 s | 8.02, 8.20, 8.20 s |
| `slow-unit-f240…` | | 4.70 s |

8.2 s exceeds 20 s ÷ 2.6 ≈ 7.7 s, the local budget that corresponds to the
Linux limit, so this is a budget, not a runner fluke. `surface_editing` now
has the 60-second per-input budget of `surface_knots` and
`degree_elevation`. No assertion changed.


## Analytic intersections: a nearly tangent cylinder and sphere

`analytic_intersections/timeout-6b3c325edc7880d84c46ffe7154a5bfcdbde2276.bin`
was found by a local 120-second campaign before S7b.1 was committed: a
cylinder of radius 4.125 about the axis `(3, 3, 3)` and a sphere of the same
radius whose centre is offset from the axis by `8.25` along the frame's `x`
(rounded): nearly tangent from outside, a loop about `1e-8` wide. The first
implementation isolated the loop's ends by subdividing the discriminant in
binary64 intervals and then in rational ones, with the natural interval
extension: beside a nearly tangent loop that needs pieces as small as the
square of their distance from it, and the input did not finish in 600 s. The
discriminant has a closed form in the canonical frame, so the ends are now
`+-arccos c` by one certified arctangent; the input replays in 0.09 s
without a sanitizer.

## Analytic intersections: a traced curve's point at its track's end

`analytic_intersections/crash-a90f9c5a16e0f0758da7f17d68ce3a9ac2822040.bin`
was found by a local 300-second campaign before S7b.3b.1 was committed: a
cone and a torus (a kind byte of 224 or more) sharing an origin. The check
evaluates every track at `lo + (hi - lo) j / 4`; for `j = 4` the sum rounded
one unit beyond the track's end, and `TracedCurve::point_at` refused the
meridian as outside the track. A meridian within rounding of the range's end
is now clamped to it (and a branch between two boxes' sides at one meridian
keeps its root as a window).

## Analytic intersections: rational point evaluation on a traced curve

`analytic_intersections/timeout-5b31d67b3c1f9d25363f62841d8a4fe690bfd195.bin`
was found by the same campaign: a torus of major radius 0.5 and minor 0.25
and a cylinder of radius 0.25 through its centre on an oblique axis, with two
exact tangencies. Near the folds binary64 intervals resolve a track's root
only to about their rounding over `|G_t|`, so points fell back to rational
intervals, narrowed by bisection: dozens of rational cosines per point, over
20 s per input under AddressSanitizer. The rational tier now narrows the
binary64 bracket by the interval Newton operator (a few evaluations); the
input replays in 1.5 s without a sanitizer.

## Analytic intersections: steep branches of two thin tori

`analytic_intersections/timeout-248e3b316e70fe35093e805d7ae31e09c7a5ac09.bin`
was found by a local 600-second campaign before S7b.3b.2 was pushed: two
tori of major radii 4.375 and 3.125 with thin tubes (0.254 and 0.483) on
perpendicular axes. Their loops are narrow in the carrier's meridian angle
and tall round its tube, so beside the folds the branches climb steeply and
the tracks' windows, predicted from the slope alone, kept missing them; steps
fell below `1e-6` and were retried in rational intervals, 7.5 s without a
sanitizer. The windows are now predicted to second order (the branch's
curvature from implicit differentiation) and rational steps are tried only
below `1e-9`; the input replays in 0.8 s.

## Analytic intersections: a fold box's edge crossed beside a twin fold

`analytic_intersections/timeout-8dce95c9ef5b34d79c816425397a39d342103b2f.bin`
timed out in Linux CI's scheduled run (42 s locally without a sanitizer): a
torus of major radius 0.5 and minor 0.25 about `(3, -1, -1)` and a cylinder
of radius 0.25 through its centre about `(-1, -1, -1)`, 10 degrees off the
torus's equatorial plane (the family of `timeout-5b31…` above). The curve
crosses itself at two exact tangencies and has eight folds in four close
pairs, `3e-4` apart in the meridian angle: S-bends of a branch. The first
fold box on which `G_phi` keeps a sign is tall enough in `t` to reach the
twin fold's branch, which crosses its top or bottom. Binary64 intervals
descended to the edge subdivision's floor at the crossing, then rational
ones repeated the descent with a rational cosine at every fresh midpoint,
about 2 s per box. A certain sign change of `G` between exact points of an
edge now proves the crossing at once, and the rational tier is skipped (no
tier can free that edge); the next, smaller box succeeds in binary64 as
before, every certified result unchanged. The input replays in 2.3 s
(18.5 s before, measured together); the three slow units the run reported
beside it (1.4 to 1.7 s of CPU each, a torus and a cylinder, two cones, two
tori) are the traced graph's ordinary work and are unchanged.

## Analytic intersections: a plane an ulp off a torus's axis on Linux

`analytic_intersections/crash-6aaf4957eb5df03ee6bcfb941a50682840e8c058.bin`
crashed only in Linux CI's per-push run at `428349e8` ("a parameter inside
a component is on the curve"): a plane with normal `(-3, -1, 2)` through the
centre of a torus (major radius 1.375, minor 0.15) about `(-1, -1, -2)`,
both at `(6.875, 6.9375, 6.875)`. The normals are perpendicular, and on
macOS the stored ones stay exactly so: the plane contains the axis, two
meridian circles. glibc's `hypot` rounds correctly, and there each stored
normal is an ulp off macOS's and their dot product `2.3e-17`, so the plane
meets the torus in S7b.3a's two loops, each within `2.5e-18` of the meridian
angle `+-pi/2`, whose nearest binary64 value is `6.1e-17` away: their ends'
enclosures crossed, and every parameter the check took between them was off
the loop. A procedural curve's loop now has disjoint end enclosures with
both branches defined between them, or the intersection is
`ComputationLimit`; `toroidal.rs`'s tests build glibc's frames bit for bit
and turn either normal by an ulp or two on any host. On macOS the input
still gives the two circles.

## Analytic intersections: twin cones' points at infinity beside `pi/2`

`analytic_intersections/timeout-949542bf24afc86a8dae30e2bf314dc51483511c.bin`
timed out in Linux CI's scheduled run 36716623883 at `428349e8` (over the
target's 20 s under AddressSanitizer; 11 to 16 s locally without one):
two congruent cones (reference radius 0.125, half-angle `0.0708`) on
parallel axes `(1, 1, 1)`, the first at `(4.0625, 4.0625, 4.0625)` and the
second `0.75` along its frame's `x` (the harness's tangent mode). They meet
in the conic of their radical plane (S7b.4's twin factor), whose branches
run to infinity where `B` vanishes, on the rulings normal to the offset:
the rounded offset puts those crossings about three ulps inside `pi/2`'s
binary64 upper bound, a subdivision point of `roots_along`. Binary64
intervals cannot sign `G` there and descended to their floor; rational
intervals signed it at once, but then halved the rootless piece beyond it
(`G` monotone, of one sign at both ends, nearly zero at one) down to their
floor as well, four rational cosines per piece, and bisected each root with
a rational cosine per step: about 6 s per call, and the check makes two. A
piece on which `G_u` keeps a sign and whose ends' certain signs agree is
now excluded at once, and the rational tier takes what binary64 intervals
settle (a piece's exclusion, a point's sign, a root's bracket down to their
rounding) before its own evaluations, ending each root with interval Newton
steps and bisection. Every root is the same adjacent pair of binary64
values as before: 896 `roots_along` calls over the corpus, these
regressions and 600 random pairs of cones (two thirds of them on parallel
axes with equal half-angles) agree bit for bit. The input replays in
0.8 s (10.9 s before, measured together).

## Curve/surface: a circle nearly tangent to a cone

`curve_surface/timeout-fdb4ae721b33aab65287ff8213216b9633bd2975.bin` and
`curve_surface/slow-unit-5141896143c65fa4feb23b2d35ce6aa19345b013.bin` were
found by the first local campaigns of the S7c.1 target: a circle beside a
cone at the sum of its radius and the cone's (the harness's tangent mode),
touching the cone's section circle up to the rounding of the stored frames.
Binary64 intervals cannot settle the boxes around the near-tangency, and the
rational tier then subdivided the whole turn with rational cosines: over
20 s under AddressSanitizer. The rational tier now revisits only the boxes
binary64 leaves undecided, merged, under a budget of 16 boxes, with the
cone's apex and axis computed once; the inputs replay in 0.19 s and 0.24 s
without a sanitizer.

## Split: a round hole off the plane

`split/crash-39ab84da8a9b4d13d549164b53edee0890ef9fb2.bin` was found by the
first local campaign of the S8a target: a square with a round hole split by
a plane parallel to its axis that misses the hole. The hole's whole circle
enters the section as one piece from a placeholder vertex, which had no
stored position, and tracing its cycle's orientation looked it up. The
placeholder is now a point of the circle; the input replays in 0.008 s.

## Split: a plane passing a hole's vertex within the resolution

`split/crash-769f29ee776d0370d29778bce537ec268da15e07.bin` came from the
second campaign: a square with a square hole split through an outer vertex
by a plane that passes a hole vertex within the resolution without touching
it exactly. The hole of one piece then lies within the resolution of that
piece's boundary, which a profile rejects (`InvalidHole`); the split now
reports it `Degenerate`, as it does other sub-resolution pieces.

## Split: an oblique plane grazing a hole's circle

`split/crash-6d802bdd1a8deb1da573d4cdf1b35e71fcb7168e.bin` came from the
first S8a.2 campaign: in the tilted frame a plane meant to touch the round
hole's circle at the top cap crosses it at two points close together (the
frame's stored axes are not exactly orthogonal). The section's hole cycle
then starts within the resolution of its piece's outer boundary, which no
piece owns; the split now reports it `Degenerate`, as it does a piece's
hole within the resolution of its boundary.

## Split: a winding loop closing within rounding of a turn

`split/timeout-387b23f8d333d1b9c65ea4672d2af5b6c964cce5.bin` came from the
second S8a.2 campaign: an oblique plane cuts a square's round hole across
its top cap, so the hole's wall in the lower piece has an upper loop of an
arc and a sinusoid ending within rounding of its start a turn on. The
binary64 mass tier could not sign that closing chord's `du` and fell back to
rational intervals (6.3 s in release, a timeout under AddressSanitizer);
the chord is now enclosed over its box in either tier (0.3 s).

## Split: a cap's rim crossed within rounding of its tangent

`split/crash-69190d85572d409cd66f328d6e3a1fbeb4878c66.bin` came from the
first S8d.2 campaign: a cap (radius 3.5, from latitude 0.25 to its pole) cut
by a plane through its rim's point containing the rim's tangent there. The
point is rounded, so the plane crosses the stored rim at two points 1.4e-7
apart (above the resolution) around a segment whose sagitta is 7e-16. The
tiny rim arc's side and the section arcs' in-band status were decided by
midpoints within rounding of the plane and the rim, and both rim arcs fell
on one side (an open loop). Sides are now decided on the arc farthest from
the plane or the ends and propagated (a crossed rim's arcs lie on opposite
sides; arcs alternate at crossings, not at touches), a crossing whose segment
is thinner than the resolution is `Degenerate`, and a projection's
quadrature never falls back to rational jets (a fix of the same run: the
sliver had sent the mass to them for 204 s). The input replays in 0.01 s.

## Split: an oblique plane along a U's notch within rounding

`split/crash-98c8688e871a05f212a6cce95a17f93650cbe0d2.bin` came from the
second S8d.2 campaign (a prism input): in the tilted frame a U (sides 6.75
and 4.375) cut through the notch's inner corner by a plane whose trace on
the bottom cap runs along the notch's floor. The corner was rounded into the
frame, so the trace passes both ends of the floor within 1e-15 without
meeting either, and the section's cycles could not close. A stored profile
vertex off the trace by less than the resolution (exact on the stored data)
is now `Degenerate`, as a crossing that close to its segment's ends was. The
input replays in 0.002 s.

## Split: a cap cut near its pole (slow unit)

`split/slow-unit-4238fd63c6db929689705944f50d5785839ac17b.bin` was saved by
the clean 600-second S8d.2 campaign at `a4e1c9df`, a slow-unit diagnostic
(10 s under AddressSanitizer, 0.92 s in release) with no wrong answer: a cap
(latitude 0.25 to its pole) cut by the target's plane mode 7, whose circle
passes near the pole, where the projection quadrature refines while
the section's angle about the axis turns fast. It is within the split
target's 60-second input limit (S8d.2).

## Split: a torus section's angle on atan2's branch cut

`split/crash-412e2b4737d69bc31fd08b7e0b52e0a139a45e98.bin` was replayed from
the corpus by the first S8d.3 campaign: a thin torus (major 3, minor 0.25)
cut by a plane whose normal has no `y` part and a negative `x` part, so the
section's loops are graphs over `v` about `u = atan2(0, negative) = pi`, on
the binary64 tier's branch cut, where it declines. Every point of the rings
was then unknown and validation could certify no gap. The angle is now that
of the normal's trace turned back by its binary64 value, as the jets' angles
are. The input replays in well under a second.

## Split: tori cut in spiric sections (slow units)

`split/slow-unit-3965cc7684d5f1398a09a86ce4dfb856336444af.bin` and
`split/slow-unit-3c49096c4d0a4a3e60cdddef2eb662fe70cb62ff.bin` were saved by
the clean S8d.3 campaign at `9c21cd2a`, slow-unit diagnostics with no wrong
answer: whole tori cut obliquely (the target's plane modes 5 and 8), split
and moved rigidly, about 4 s each in release then, most of it certified
projection quadrature on the pieces' torus faces. Moving a piece's mass
properties with it instead of enclosing them again (`0df34152`) halved
them. Within the split target's 60-second input limit.

## Split: spline prisms cut obliquely (slow units)

`split/slow-unit-67ae28f45bf8c12f623afa46bd15d5a35c6d09d2.bin` and
`split/slow-unit-9ff692515f26ff227bc0ec1be47e543351fc3f83.bin` were saved by
the first S8b campaign at `66c32112` (no crash): a square with a lens hole of
two cubics given clockwise, and a quadratic bulge, in the tilted frame, cut
through a vertex at an angle (the target's spline prisms, plane mode 7), split
and moved rigidly. The first took 66 s under AddressSanitizer, beyond the
60-second input limit (7.4 s in release): exact rational halvings in the
separation screen, a cycle's orientation area from 256 exact evaluations,
and an oblique piece's certified mass enclosed again after a rigid motion.
The screen now halves outward binary64 boxes, the orientation area is a
binary64 polynomial integral, and a moved oblique piece takes its source's
mass, moved (`55bb0516`): 1.6 s and 0.7 s in release.
`split/slow-unit-e63ac50da1489df1eb606f9b800f18074062b053.bin` was saved by
the clean campaign at `5c50de15`: the lens hole again, cut at a dyadic
point at an angle (plane mode 6), 18 s under AddressSanitizer and 1.7 s in
release, within the limit.

## STEP: a torus band's ring starting a rounding below its seam

`step/crash-6b0ff1a4f7646a3f2fce7e959be0fc3033928e74.bin` was found by the
first clean campaign of the STEP-a target: the quarter-torus elbow with one
meridian's reference direction scaled and flipped, so the circle's
parameter starts on the tube's inner side. The import is valid (a torus
band between two rings wound in `v`), but its `.brep` text did not read
back: the writer reduced the `+v` ring's start, `-1.2e-16`, into `[0, 2 pi)`
for the latitude seam and left the rings' pcurves where they were, a period
below it. The writer now shifts a face's line pcurves into `[v0, v0 + 2 pi]`
as it does in `u`; the input replays in 0.06 s.

## STEP: a line pcurve on a spline surface that is not its edge's length

`step/crash-19bd352aea32baea4b0ecbdf14c4a2985ce70d18.bin` was found by the
first smoke run of the STEP-b target (60 s, local): the B-spline prism,
whose side surface's `v` runs over `[0, 1]` while its vertical edges are 4
long. The import is valid (the file's line pcurves over the fins'
fractions), but its `.brep` text did not read back: the writer wrote every
line pcurve at unit speed against its edge's parameter, right on planes and
along rulings and circles where the pcurve's length is the edge's, wrong
here (`(0, 4)` for `(0, 1)`). A line pcurve whose length is not its edge's
span now runs over the span; every other one is written as before. The
input replays in 0.05 s, and `tests/step.rs` round-trips the four spline
bodies.

## Boolean: a tool filling the object's hole over part of its height

`boolean/crash-769e2af11386b1251080b11d3dc937f416957d63.bin` was found by the
first campaign of the S9a target at `56f26e85`: a stadium inside a square
with a round hole, filling the hole over part of the square's height, fused
in the tilted frame. The heights differ, so the result is a stack (S9a.2);
the kernel returned the square unchanged because the fused profile's
boundary came from the square alone (the hole's pieces, inside the stadium,
were dropped), which it took for the stadium lying inside the square. The
volumes' identity caught it. Containment in 2D is now an empty exact cut,
and identical profiles are two empty cuts; the input is `OutOfDomain` and
replays in 0.03 s.

## Split: a zone's rim one ulp off its circle

`split/crash-674ebdd30b135f3807e5587ebe7e453e412265cb.bin` was found by a
600-second campaign of the split target at `56f26e85`: a zone of a radius
4.25 sphere between latitudes -1.375 and 0.75, split by an oblique plane
crossing its upper rim. The rim's arcs on the pieces lay on a circle one ulp
above the zone's rim, so the history check (debug build) found the split
arcs off their parent's support. The split's end heights took the
latitude's sine next to its cosine, which the compiler fused into one
`sincos` whose sine differs from `sin` in the last place on macOS; the
zone's builder took `sin` alone. Every revolved builder now takes a ring's
height and radius from out-of-line `scaled_sin` and `scaled_cos`; the input
replays in 0.5 s.

## Boolean: spline stacks' slow unit, a stack touching itself, a chained stack

`boolean/timeout-d0a3de29450a8c2ba8c2ce0f33b80aac23e5e11d.bin` was saved by
the 600-second campaign of the boolean target at `cfc641c9` as a timeout (26
seconds under AddressSanitizer against the 20-second limit): lens-hole
spline profiles as stacks, their three results measured (the certified mass
enclosure of spline walls), then each piece moved rigidly by rebuilding its
stack in the moved frame. A result's rigid motion now moves its stored
geometry, Bernstein degree elevation divides once by the degree in the
certified tier (it bracketed `i / n` from rationals), and the input replays
in 1.8 s in release and 27 s under the sanitizer: the target's per-input
limit is now 60 s, as the split target's.

`boolean/replay-7baa6762d444c7fb869dfa8d1d0f6f44fbaaf8f44ecae5253a71d003a84400ee.bin`
and `boolean/replay-01782c8ee258a903615e59c13d8449ca8de1d520.bin` are corpus
inputs that failed when S9b.2's chained stage was replayed over the corpus
with debug assertions (original corpus names). The first (at `cfc641c9`
too, before the chained stage): a square with a lens hole whose corner lies
on the edge line of a block standing on it; the stack touches itself at
that vertex and validation reported a non-manifold vertex (`InvalidTopology`)
where S9a's rules refuse (`Degenerate`); it replays refused. The second: a
tilted stack given to a Boolean again, whose stored model's cap triangles
spanned a run of collinear edges without their middle vertices (stored
vertices, not collinear exactly), leaving the model open; stored faces are
now triangulated on their own vertices (ears clipped), and
`tests/polyhedral_booleans.rs` holds the stack. Chaining larger results
took up to 165 s an input under the sanitizer (the stored models' exact
fragments), so the target chains first results of at most 12 faces (this
input's 21-face stack is no longer chained).

## Boolean: identical prisms with a spline hole given either way round

`boolean/crash-0072749f9f70085aacd7c92185839596132c5c59.bin` was found by the
600-second campaign of the boolean target at `e95fdfbc`: two identical
tilted squares with a lens hole of two cubics, the hole given
counter-clockwise in one and clockwise in the other, fused. The fused
profile traces the hole the other way round from one input, so the
result's spline walls are that input's reversed in u with the other sense
and its spline edges that input's curves reversed; the history checker
compared a wall's representation (its sense first) and a curve's
parameters, and found the merged support different. It now takes a spline
wall and its reversal in u with the other sense as one oriented surface,
a plane and its opposite normal with the other sense as one oriented
plane, and a spline piece on a whole's reversal at the piece's parameters
(`tests/booleans.rs`, `identical_prisms_with_a_spline_hole_fuse_into_one`);
the input replays in 0.24 s.


## Boolean: the campaign's slow unit at S9c.1

`boolean/slow-unit-7f11d0a5e88411231c5b1054a237eb1f79c8cccb.bin` was saved by
the 600-second campaign of the boolean target at `80ba3d9c` (S9c.1 merged),
which was otherwise clean: 363 mutation executions after a 2,607 s replay of
1,081 inputs, 35,013 edges, 486 MB peak. The input (16 bytes) takes 40 s
under AddressSanitizer and 2.8 s in a release build with debug assertions,
within the target's 60 s limit; it passes every check and is kept so the
replay covers it.

## Boolean: a sphere less a patch clear of its poles

`boolean/replay-1608293dab85fd14619a26a3b966e3fff7eb24e5.bin` is a corpus
input that S9d.1's sphere tool (the spline byte's top bits) turned into a
tilted U whose corner and bottom cap pass through a whole sphere's centre.
The fused sphere face was the sphere less a four-sided patch below its
stored pole, one contractible loop running as a hole, which the validator
refused (`loop_winding`: only a torus's face could be its surface less its
loops). A sphere's face may now be too, its mass the whole sphere less the
loops' integrals (`tests/sphere_booleans.rs`,
`a_sphere_less_a_patch_clear_of_its_poles`, checks the volumes).

`boolean/replay-709962d223ef329daa4670824ca32c3bce44136d.bin`, from the same
replay, fused a tilted prism with a sphere whose sections lie in its faces:
the history checker compared a split circle with its whole by their frames
(equal only when the circle's normal is the planes' exactly), and a
section's circle rounded from its rational basis, or a rim split at the
sphere's seam, differs from its whole's frame within rounding. Circles now
compare by their centres, normals and radii within tolerance, and a circle
in a plane by its centre's and axes' ends' distances, as an ellipse.

## Boolean: a frustum on a tilted prism's top within rounding

`boolean/replay-5b0afacf4f57b4067337eb44e6c433fc4be8ca0b.bin` is a corpus
input that S9d.3a's cone tool (the spline byte in `160..192`) turned into a
frustum standing on a tilted diamond prism's top: its frame's origin, the
prism's frame's point at the top's height rounded to binary64, lies 2.7e-16
off the top's exact plane, so the frustum's base disc is parallel to the
top and apart from it by far less than the resolution. The top cut the wall
in a ring beside its rim, and the sliver between them could not be nested
by the loops' binary64 images: the result was open (`InvalidTopology`). Two
parallel faces apart by no more than the resolution are now one plane
within it, `Degenerate` (`tests/cone_booleans.rs`,
`a_frustum_on_a_top_within_rounding_is_degenerate`).

## Boolean: a cone's apex at its base carrying a meeting

`boolean/replay-26c72abf2f24213912aa6c5455e5c8dd9e01c467.bin` is a corpus
input that S9d.3a's cone tool made a cone with its apex at its base (the
spline byte's radii `(0, r)`) against a prism with arcs, which S9d.3b.1
evaluates: its meeting with the prism's cylinder is a ring over the cone's
angle, `Curve3::Meet` with the cone as its carrier and radius zero at the
frame's origin. The validator took the radius for a cylinder's and refused
the curve (`degenerate_curve`); a cone's radius there may be zero (three
more corpus inputs alike). `tests/cones_booleans.rs`,
`a_cone_with_its_apex_at_its_base_carries_a_ring`, checks the volumes.

## Boolean: a frustum's caps an ulp off parallel on Linux

`boolean/replay-6d1fd0617f7790963cdb479676b89690a3868d96.bin` is a corpus
input that failed only on Linux CI (the scheduled and pushed runs at
`fc695afd`): a frustum over a tilted diamond prism's heights, its base frame
built from the prism's normal and so normalized again. `Vec3::normalized`
divides by the platform's `hypot`, which rounds the tilted normal's length to
1 on macOS and not on every Linux or Windows runner, so there the frustum's
caps are an ulp off parallel to the prism's bottom and top: not coincident,
and not refused as a sliver either, and the result's loops wound wrongly
(`InvalidTopology`, `loop_winding`). Two planar faces within the resolution
of each other across their boxes' overlap, parallel or not, are now one
plane within it (`Degenerate`); on macOS the caps stay exactly parallel and
coincident and the Boolean evaluates. `tests/cone_booleans.rs`,
`a_frustum_turned_by_an_ulp_on_a_top_is_degenerate`, turns the frame's normal
by an ulp on any host.

## Boolean and analytic intersections: Linux CI timeouts before the speed-ups

`boolean/timeout-f3138f9fc2095c3090ac8ee9b4382b117b163ac2.bin` (a frustum
tool against a prism, 60 s under AddressSanitizer on Linux CI at `b33fd65b`)
and `analytic_intersections/timeout-b4449b175f865da524262dc5bc9ed2287e10db1b.bin`
(over its 20 s limit there) were found by the per-push fuzz run before the
curved Booleans' exact arithmetic was sped up (Lehmer gcd, isolators kept per
field) and before the torus curve graph's fold boxes stopped at a proved
crossing. They replay in 0.9 s and 0.6 s after both (release, no
AddressSanitizer); kept so the per-push runs keep timing them.

## Boolean: a three-quarter wedge's holes, a cavity in a torus part

`boolean/replay-80c01a1abdc7d0548c96f72a9cf9ec4c2d575313.bin` and
`boolean/replay-0028ade227457ca9ef8fad7b8892eaa955687d3f.bin` are corpus
inputs whose spline byte and flags were set to S9d.4b.1's torus parts (the
byte in `148..160`) for a local replay with debug assertions of 6,579 such
variants before the part tool's first campaign. The first: a wedge of three
quarters of a turn coaxial with a tilted square frame with a square hole,
fused. The fused wall is a patch over more than half a turn with no loop
winding, and its holes' projection pcurves, lifted from the surface's
principal angles, lay a turn away in `u` from its outer loop
(`inner_loop_outside`); a torus face none of whose loops winds now shifts
each hole by whole turns into its outer loop on the cover
(`tests/torus_segment_booleans.rs`,
`a_three_quarter_wedge_fused_holds_its_holes`). The second (33 inputs
alike): a wedge in a turned frame fused with a rectangle, the result cut by
the turned box lying inside it: a cavity in a solid bounded by a torus face
with loops, whose containment the validator's rays cannot decide
(`uncertified_containment`, as for spheres and whole tori with loops). A ray
missing the whole sphere or torus now counts no hit, and a Boolean whose
result reports only an undecided containment is `ComputationLimit`
(`a_cavity_in_a_segment_is_undecided`). Both replay in under 3 s.

## Boolean: a torus band's certified integrals

`boolean/slow-unit-17e131e351b53a29f4966c6dd07b19837efa0626.bin` is a
corpus input (16 bytes): a torus band (the tube from `v = 0.5` to `2.25`
with its cap, the spline byte 158) in the tilted frame against a prism,
fused, cut and in common. Most of its time was the validator's and the
mass's certified integrals along the band's projections (REVIEW_NOTES.md,
the certified-integrals track): 5.3 s in a release build with debug
assertions and 63 s under AddressSanitizer on a host at load 10 to 20
before that track, 1.6 s and 17 s after. It passes every check and is kept so the
replay times it; a twin of it (`8997...`, another height) is not.

## Boolean: a sphere's pole at a cone's apex height, a loop on a cone's far nappe

`boolean/replay-775bd1c0c918547b46fc02bd6fa93e8a5b139c75.bin` and
`boolean/replay-e7003ca0ed37a9487a672376ba20ed2be1106cf9.bin` are corpus
inputs rewritten into S9d.3c's decode (the object a cone against a turned
zone tool, and a sphere against a turned cone tool) for a local replay with
debug assertions of 2,842 variants before S9d.3c's first campaign. The
first (25 variants alike): a zone's pole at the frustum's virtual apex's
height, tested against the frustum's `Curve3::Meet` pieces with the sphere,
reached the carrier's place of a point at its apex height (a radius of
zero: `a meeting's point off the apex`, a panic); such a point is on none of
a cone carrier's pieces (the apex itself is refused before). The second:
the sphere met the double cone in two loops, one on each nappe, and the
height graph's turning points on the far nappe (the carrier's radius
negative there) were placed on the wrong side of the angle `phi`, four
events at one place (`two turning points within rounding`); the side is
now `g / rho`'s sign, and the branch the sphere's gradient along the
ruling. Both replay in under 2 s.

## Boolean: a void under a band's hole, a loop through a sphere's pole

`boolean/replay-a8c53a2351252e92a262cbd7ca85b060ac023caf.bin` and
`boolean/replay-0ca818deb641a61baf12a95a7802f376d308a30e.bin` are corpus
inputs rewritten into S9d.4c's decode (a cap above latitude 1/4 against a
torus band `0.5..2.25`, and a whole sphere against a quarter wedge stood on
its side) for a local replay with debug assertions of 3,000 variants before
S9d.4c's first campaign. The first (seven variants alike, and a box or a
pipe over the band's hole the same way): the cap's plane lies just above
the band's upper end disc and covers its hole, so their fuse holds a thin
void bounded by the disc, the band's inner wall and the cap's plane, a
shell of both inputs' faces that the assembly made a solid of its own
(`shell_orientation`); a shell of both inputs' faces whose certified flux,
built alone, is turned inward is now a cavity of the result, whose
containment the validator's rays leave undecided (`ComputationLimit`, as
S9d.4b.1's). The second: the wedge's end half-plane, on the rounded
direction of its quarter turn, passes within rounding of the sphere's
poles, and the sphere face's loop through the pole had its winding from its
pcurves' changes alone, half a turn at the pole left out (`uv_gap`); a
loop's winding is now its last pcurve's end against its first's start.
Both replay in under 2 s.

## Boolean: a cavity beside a lens hole's spline walls

`boolean/crash-1962418372712ecacc6f188064195b85e80eb14b.bin` was found by
the scheduled run 36868817257 (at `0dbd7c44`), mutated from
`replay-26c72abf…`: a square with a lens hole of two cubics in the tilted
frame, less a holed square inside it in 2D and over the middle of its
height. The cut leaves a closed void around the lens, a valid solid, but the
validator's rays do not decide spline faces, so the void's containment is
`uncertified_containment`. The polyhedral and curved results map a result
failing only on that (or on `uncertified_loop_winding`) to
`ComputationLimit` (S9d.4b.1, S9d.4b.2a); S9a.2's stacks did not and
returned `InvalidTopology`. Both now take one rule (`undecided` in
`solid/boolean.rs`); `tests/booleans.rs`,
`a_cavity_beside_spline_walls_is_undecided`, checks it in the tilted and the
axis-aligned frame, the fuse being the object and the common the tool. It
replays in 0.15 s.

## Split: a spiric section filling AddressSanitizer's stack depot

`split/oom-dfe03a750d28f053865e6e878831548356994524.bin` is the input the
scheduled run 37008675181 (at `0dbd7c44`) was running when the process
crossed the 2 GiB RSS gate: a whole torus (major 1.625, minor 0.875) in the
tilted frame cut by an oblique plane, a spiric section. It is not an
unbounded allocation: alone it runs in 0.28 s and 8 MB with debug
assertions, and 49 variants with the torus's axes or the plane's normal
turned by one or two ulps (the Linux `hypot` lesson) all split alike in
0.2 s within 21 MB; at the gate the run held 50 MB live and 59 MB
quarantined, and its RSS had climbed steadily from 75 MB through 2,770
inputs. The rest was AddressSanitizer's stack depot, which keeps every
distinct allocation and free stack for the process's life: this input
alone records 985,416 distinct 30-frame stacks (248 MB; 493 MB peak RSS
alone under the sanitizer), its certified projections' rational interval
arithmetic reached through deep and varied call paths, and 1,400 corpus
inputs replayed locally left 3.9 million (941 MB of depot, 1,350 MB RSS).
`split` now keeps five-frame stacks as `analytic_intersections`,
`curve_surface`, `curve_curve` and `boolean` do (`SHORT_STACK_TARGETS` in
`tools/run_fuzz.py`): 14,849 stacks and 255 MB for this input, 88,036 stacks
(10 MB) and 434 MB RSS for the 1,400. The 2 GiB gate, the quarantine and
the input limits are unchanged.

## Boolean: spline walls in any position (S9f.1)

Three inputs replayed before S9f.1's kernel was committed, none from a
campaign. `boolean/replay-ec02b1bf118dc78762cd2737908f5eeb5c9e1f83.bin`
and `boolean/replay-d8d4cf6a46ae12c670d47796801753bbb0316362.bin` are
spline variants of the corpus (the decode's frame byte set to a turned,
leaning, tilted or side tool, the spline byte to a spline object or tool,
the byte at index 14 to R4's knot for every other object; 2,852 replayed
with debug assertions): a bulge prism apart from a leaning box at the
height of its joints, whose fuse, cut and common came out open
(`InvalidTopology`), the `+u` ray's crossings with the spline counted by
the opposite predicate to the chords' at a joint at the point's height; and
a tilted box's hole wall in the plane of a lens hole's corner, which the
kernel took as a plane crossing the prism along the joint's edge (a
non-manifold cut) where it touches it there (`Degenerate` now, the
arriving and leaving segments on one side of the plane).
`boolean/replay-197f00ac2a7f93dcd90d3a794ededf7d90f12497.bin` is a corpus
input whose index 14 byte now decodes R4's knot: in one frame S9a.2 cut the
spline across its knot of multiplicity two and the rounded poles next to it
left the result's profile off C1 (`InvalidCurve`); pieces are now cut from
the curve with the knot removed once (`spline::piece`), as 58 other corpus
inputs found. `tests/spline_any_booleans.rs` checks all three
configurations. They replay in 0.2 s, under 0.1 s and 0.3 s with debug
assertions (the third in 4 s under AddressSanitizer).

## Boolean: a wall within rounding of tangency to the chained stage's sphere

`boolean/crash-d238291d9edbe60570ae31479609845872d763c0.bin` is from the
local campaign at `6582f379`, mutated from a corpus input. It decodes two
stadiums in the shared tilted frame (the object `2.75` by `1.0` over
heights `0..1.75`, the tool `1.25` by `2.25` stacked on its top), and its
chained byte gives the cut's one solid, the object itself, to the
`GIVEN_BALL` sphere of radius 1.25 about the turned box's frame. The
sphere's centre lies `1.25` from the stadium's flat wall exactly, but the
frames' rounding put the wall 4e-16 inside the sphere: the plane's section
was a circle of radius 1e-8, and the chained cut and common failed with
`InvalidTopology(degenerate_curve)` (the plain stadium against the sphere
too). A plane crossing a sphere within the resolution of tangency (a cap
no higher than the resolution) is now `Degenerate`, as S9d.2c's circles
within the resolution of tangency are; a plane missing it by less stays a
miss (`tests/sphere_booleans.rs`,
`a_wall_crossing_within_the_resolution_of_tangency_is_degenerate`). It
replays in 0.3 s with debug assertions, and 348 single-byte mutations of
it replay without a failure.

## Boolean: two cylinders' axes within rounding of parallel

`boolean/crash-3e989fac80f9cef1f1a6abe1d874ec76855724c3.bin` is a replayed
variant of the boolean target with its shipped switches, which panicked in
the validator at `7199e06a` and `93e6fcd0` (`a conic or section
evaluates`). It decodes a square `8.5` by `8.5` with a round hole of radius
`1.375` over heights `0..2.25` in the tilted frame and a torus band
(`-2.5..-0.25`, radii `2.0625` and `1.03125`) in the same frame; the fuse's
first solid, the holed prism itself, goes to the `GIVEN_ROUND` cylinder of
radius 1.25, whose frame's normal is the tilted frame's normalized again,
an ulp off it. The two cylinders' models crossed within rounding of
parallel, so S9c.2b.1's quartic gave a meeting of an ulp's sweep of the
carrier's angle, stored on axes exactly parallel where it has no point.
Two cylinders whose axes are within `10^-12` of parallel without being
parallel are now `Degenerate` unless certainly apart within their faces'
bounds, where they evaluate as before (`meet::cyl_pair`;
`tests/turned_booleans.rs`,
`cylinders_within_rounding_of_parallel_are_degenerate_unless_apart`, the
holed prism, a rod and the given result, and pairs apart). The torus band
is incidental: the holed prism alone against the cylinder panicked too.
It replays in 0.6 s with debug assertions, and 363 single-byte mutations
of it (every value of the chained byte among them) replay without a
failure.

## Boolean: a spline wall's and a cylinder's axes within rounding of parallel

`boolean/crash-a3fb9da3e3d83ef1cc3b272a30e6716eed8145b2.bin` was found by
an audit of the other pair kinds after the input above, decoding variants
of the boolean target with its shipped switches; at `b903f3fa` it failed
(`unexpected error invalid topology: degenerate_curve`), as did 15 more of
the 2,304 searched (the object's spline profile, its sizes, its height and
the chained operation varied; two of them with `PrecisionLoss`), none
since. It decodes a square
`5` by `5` with the split target's lens hole (two cubics reaching `x =
+-1.25`, `y = +-0.9375`) over heights `0..0.5` in the tilted frame, and a
rectangle tool above it (`1.5..2.5`), so the cut's first result is the
prism itself; the chained byte gives it to the `GIVEN_ROUND` cylinder of
radius 1.25, whose frame's normal is the tilted frame's normalized again,
an ulp off it. The cylinder crosses the lens's spline walls on an axis
within rounding of theirs: along the walls' rulings its function's `w^2`
coefficient is within rounding of zero, and the meeting within the faces a
sliver of the run no binary64 edge holds. A spline wall and a cylinder
whose axes are within `10^-12` of parallel without being parallel are now
`Degenerate` unless certainly apart within their faces' bounds, as two such
cylinders are (`spline_crossing::meeting_with`;
`tests/spline_crossing_booleans.rs`,
`a_spline_wall_and_a_cylinder_within_rounding_of_parallel_are_degenerate_unless_apart`).
It replays in 0.1 s with debug assertions, and 605 single-byte mutations of
it (every value of the chained byte among them) replay without a failure.

## Boolean: a meridian loop through one pole of the chained stage's sphere

`boolean/replay-21928984a54f016cf7e3e987dcdd3f07df37c15c.bin` is a corpus
input (its original name kept) that failed the debug-assertion replay at
`86b1d834`, the merge of S9e.4b.3a's kernel with `GIVEN_MET` switched on
(`unexpected error invalid topology: uv_gap`); it passed at either parent,
`GIVEN_MET` off at the one with S9e.4b.3a. It decodes a rectangle `4.25`
by `0.5` over heights `0..2.25` in the tilted frame and a torus band
(`0.5..2.25`, radii `1.3125` and `0.65625`) in the tilted frame offset by
an amount that rounds; the chained byte gives the cut's first result to the
`GIVEN_BALL` sphere about the middle of its first torus section
(`GIVEN_MET`), which lies on the bar's wall. The wall holds the sphere's
axis: its section is a meridian through the north pole, inside the bar,
and the sphere face's loop runs up it, through the pole and down the other
side, then closes along the bar's bottom circle. S9e.4b.3a's rule for a
loop through a pole set its winding to none without lifting its pcurves,
so its closing fin was a turn off its first unless the loop happened to
start at the pole; the cut failed `uv_gap`. A box with a wall through a
sphere's axis fails alike in 76 of 80 frames and operations at `86b1d834`
(`tests/sphere_booleans.rs`,
`a_meridian_loop_through_one_pole_closes_on_its_first_fin`). The loop's
pcurves from that pole on are now lifted by the turn, and the rule takes
only the pole the face would close at (`assemble.rs`); the chained
operations evaluate as before S9e.4b.3a, their volumes the same
(`tests/given_met_booleans.rs`,
`a_ball_about_a_section_on_a_wall_through_its_axis`). It replays in 3.1 s
with debug assertions, and 540 single-byte mutations of it (every value of
the chained byte among them) replay without a failure.

## Boolean: a cone section given on the other cone's disc

`boolean/replay-6fab9d419900cda3aad1257fe2ec8e8709205f20.bin` is a corpus
input (its original name kept) that failed the debug-assertion replay with
`CONE_PAIRS` switched on (`unexpected error invalid topology:
vertex_off_curve`); off, it decodes a prism against the frustum. On, it
decodes a frustum of radii `2.625` and `1.3125` over heights `0..0.75` in
the axis-aligned frame and a frustum of radii `1.3125` and `0.65625` in the
frame turned to the normal `(0, 3, 4)`; the chained byte gives their fuse's
first result to the `GIVEN_BALL` sphere about the middle of its first
meeting (`GIVEN_MET`). The fuse has an edge where one frustum's end disc
cuts the other's wall, a cone section on the wall's cone. Given to another
Boolean, that edge's faces are the wall and the disc, each on its own
input's model (`view`), and both models are cones: the section's curve was
rounded on the first model with a cone among the edge's faces' views,
which for the disc is the other frustum, whose frame turned the section's
plane into another one, a hyperbola off the section's vertices. 28 of the
459 cone-pair variants of the corpus (every third input moved onto a cone
object against the cone tool) failed alike, every one in the chained
stage, its partner the box, the cylinder or the sphere, and one more was
refused as a plane within the resolution of the wrong cone's apex. The
section's curve is now rounded on the cone it lies on (its frame and radii,
`ConeSec::lies_on`); `tests/given_met_booleans.rs`,
`a_box_across_a_cone_section_on_another_cones_disc`, gives the fuse of a
frustum and a leaning cone across its wall to the turned box (it failed
`vertex_off_curve` before). The input replays in 1.6 s with debug
assertions.
