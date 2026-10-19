# Work budgets and cancellation (S10)

Every operation can run under a budget: a number of work units, an
optional wall-clock deadline and an optional cancellation token. An
operation that exceeds its budget returns `Error::BudgetExhausted` with a
partial diagnosis, and another thread can cancel it. This is step S10 of
`REVIEW_2026-10-18.md`. The evidence that it is needed before any
interactive use is the torus near-hang: a rod crossing a torus's inner
equator by `1e-9` ran for over half an hour (`torus_hole_rod_booleans.rs`).
The meeting's own verification budget (`VERIFY_BOXES`) fixed that case,
but no caller could bound or stop an operation.

## Refined decisions

**API.**
- `rusty_occt::Budget` holds three things: work units (`Budget::work(n)`;
  `Budget::unlimited()` counts work but never refuses for it), an optional
  deadline (`with_deadline(Instant)`, or `with_wall_time(Duration)` from
  now) and an optional `Cancel` token (`with_cancel`).
- `Cancel` wraps an `Arc<AtomicBool>`. `Cancel::cancel()` may be called
  from any thread.
- The operations that already take a `Context` take the budget there, in
  the new field `Context::budget` or through `Context::with_budget(&budget)`.
  These are Booleans, plane and height splits, stacks, imports
  (`Solid::imported_in`), constructions and transforms.
- Operations without a `Context` run under `Budget::run(|| ...)`. These are
  validation, mass enclosures, tessellation, the intersection and
  proximity queries, and `.brep`/STEP import.
- `Budget::run` installs the budget, on the calling thread, for every
  kernel operation the closure calls. It is generic over any error type
  `From<rusty_occt::Error>`, so the DRAW worker's own error type works.
- `Context::new` leaves `budget` `None`. Every existing signature and call
  is unchanged.

**The refusal: a new variant, not `ComputationLimit`.** The review asks for
"`ComputationLimit` with its partial diagnosis". The kernel already uses
`ComputationLimit(&'static str)` for its own fixed limits (`VERIFY_BOXES`,
`folds`' boxes, tessellation's node count, an undecided comparison), and
those differ from a budget refusal in two ways:
- A fixed limit is a final refusal of the input; a budget refusal is the
  caller's choice, and a larger budget may succeed.
- Internal code catches some `ComputationLimit`s to try an alternative. The
  seam attempts continue on the seam's `ComputationLimit`, and the chain
  and form attempts of imports continue on refusals. A budget refusal must
  never be taken for one of those.

The budget's refusal is therefore `Error::BudgetExhausted(Exhaustion)`.
`Exhaustion` has three fields:
- `cause`: `Cause::Work`, `WallTime` or `Cancelled`;
- `stage`: the work-unit point that refused, a `&'static str` from the
  table below;
- `work`: the units counted, the refused unit included.

Its `Display` reads "operation budget exhausted: work units at a torus
meeting's verification box after 1001 work units". The kernel's fixed
limits stay `ComputationLimit`, and their meaning is unchanged. This is a
lasting choice of the public error type, recorded for the user. Making the
refusal `ComputationLimit` instead would mean changing that variant's
payload at about 190 construction sites and narrowing the internal
catches.

**What a work unit is.**
- A unit is one iteration of a loop the kernel already has, at the points
  where work can grow with the input. In the Booleans:
  - an arrangement's face pairs, pierces, sections, segments and face loops;
  - the near-miss rules' pairs;
  - a torus meeting's verification boxes (generalising `VERIFY_BOXES`),
    turning-point boxes, tracing steps and piece plans;
  - seam, chain and form attempts;
  - assembly's faces, edges, vertices, cavity trials and solids;
  - the polyhedral engine's clips, fragments and groups;
  - the profile arrangement's pairs.
- Elsewhere:
  - a split's section segments and loops;
  - a certified integral's pieces and boxes;
  - a spline deviation's pieces;
  - a ray's faces;
  - tessellation's rounds, edges, faces and refined triangles;
  - an intersection's root pieces, Newton steps and tracing steps.
- The table under "Work-unit points" lists each point with its stage name.
- A unit is counted, and the budget checked, at the head of the iteration,
  before the work it stands for. A unit that refuses is counted but not
  done.
- Work units are deterministic: the same operation on the same inputs
  counts the same units on every host, in every run, whatever the process
  computed before. Two rules make this hold:
  - No unit is counted inside a computation the kernel keeps across
    operations or fills lazily. Those are the field arithmetic's memos
    (`num.rs`), kept roots, resultants and eliminants, generators, the
    projection and wall jets, kept trigonometric enclosures and an edge's
    lazily drawn samples. A kept value's later use would otherwise count
    nothing.
  - Five kept computations do contain units. They are the curved engine's
    arrangements (`curved/mod.rs`'s `ARRANGED`, `pieces.rs`'s `ARRANGED`
    and `INNER`, `splits.rs`'s `ARRANGED`) and the validator's kept
    integrals (`mass.rs`'s `KEPT`). Each records the units its computation
    counted, run by run of one point (`budget::recorded`), and charges them
    again, in order, on every later use (`Recorded::replay`). A use
    therefore refuses where the computation would have, at the same unit
    and the same stage.
- A computation during which the budget refused is never kept. Its value
  may be that refusal, or a consequence of it, such as an undecided
  integral.
- Work-unit refusals are part of the contract. With a budget of `n` units,
  an operation that counts `w` units unlimited behaves in one of two ways:
  - when `w <= n`, it succeeds with the identical result and counts `w`;
  - when `w > n`, it is refused for work at unit `n + 1`, with
    `Exhaustion::work == n + 1` and `Budget::spent() == n + 1`.
- The check is exactly as granular as the count: every unit is checked as
  it is counted, so no operation exceeds its budget by more than the one
  refused unit. The time a unit stands for is bounded by the placement of
  the points. That stretch is measured as the longest wall time between
  two checks of an operation, over the corpora (below).

**Results.** No partial or invalid body is ever returned:
- Once a budget has refused, every later unit under it refuses at once
  without counting.
- `Budget::run` (which every `Context` operation goes through) replaces
  whatever the operation returns with `BudgetExhausted` once the budget
  has refused, whether that is a result, a different error or a body built
  on an undecided verdict.
- Where a loop's function cannot return the refusal (validation's
  `Option`, `Verdict` and `Vec<Issue>` returns), the unit makes it return
  "undecided" (`None`, `Verdict::Unknown`). That makes the validation
  fail, never pass, and the operation's end turns it into the refusal.

**Wall time and cancellation.**
- The deadline and the token are checked at the same points as the work.
  A deadline costs one `Instant::now()` a unit, and only when one is set.
- A refusal for wall time (`Cause::WallTime`) or for cancellation
  (`Cause::Cancelled`) is not deterministic: it depends on the host, its
  load and the other thread. It says so, and `Exhaustion::deterministic()`
  is `false`.
- How promptly a cancellation lands is bounded by the longest stretch
  between two checks. That is measured, not promised.

**Threads.**
- `Budget` and `Cancel` are `Send` and `Sync`. A `Budget`'s clones share
  one count and one refusal, so an application can keep a clone on another
  thread to read `spent()` or `exhausted()`.
- A budget, once refused, refuses every later operation run under it.
- A budget applies on the thread that runs the operation. The kernel runs
  each operation on one thread.
- Several threads may share one budget. The total stays bounded, but each
  operation's count is then deterministic only when one operation runs at
  a time.
- `Context` stays `Copy`. It holds a reference to the budget.

**Nesting.** An operation inside another runs under the outer budget when
its own context has none. It runs under its own budget when it has one: a
Boolean's inner Booleans and imports with `Context::new` (no budget) count
toward the outer budget.

**Defaults.** Without a budget, nothing is counted and nothing is refused.
A unit then costs one thread-local read, or a thread-local slow path while
a kept computation is being recorded. Every existing result, fixture and
history is unchanged. The cost of counting must stay under about 1% of the
boolean corpus's CPU (measured below).

**Not counted.** Parsing `.brep` and Part 21 text is linear in the file and
counts no units. An import's constructions and validations count theirs.
Constructions without loops (frames, primitives' topology) count none.

## Work-unit points

Each point counts one unit per iteration, under the stage name an
`Exhaustion` reports (`grep -rn 'budget::unit' rust/kernel/src` lists them).

| Where | Stages |
| --- | --- |
| Curved Booleans, arrangement (`curved/graph.rs`) | a face pair's meeting, coincidence, near miss and section; an edge pair on one surface; an edge's pierce of a face; a crossing's point; an edge's split; a section's segment; a face's pieces; a face loop's step; a hole's refinement; a stored point's split |
| Curved Booleans, near misses (`curved/near.rs`) | a face pair's near miss, grid row, Newton seed and critical point; an edge's near miss of a face or a vertex; a vertex's near miss of a face; vertex and edge pairs' near misses |
| Torus meetings (`curved/torus_curved.rs`) | a torus meeting's verification box (each box `VERIFY_BOXES` also draws), turning-point box, tracing step, piece plan and component |
| Attempts | a seam attempt (`curved/mod.rs`, `given.rs`, `pieces.rs`); a chain attempt and a form attempt (`solid/imported.rs`) |
| Assembly (`curved/assemble.rs`) | a result face's loops, a result vertex's neighbours, a cavity's trial, a result solid, edge and face |
| Polyhedral Booleans (`boolean/polyhedra.rs`) | the near-miss pairs; a face's clip; a fragment's side; a fragment edge's points; an edge chain; a face group; a result solid |
| Profiles and stacks (`profile/boolean.rs`, `boolean/stack.rs`) | a profile vertex pair, segment pair, segment's vertices and pieces, piece pair, piece's side, loop and loop pair; a stack level; a face group; a result solid |
| Splits (`solid/split.rs`, `split/oblique.rs`) | a split piece; a section's boundary, segment and loop; a section loop's step; an oblique piece |
| Imports (`solid/imported.rs`, `imported/turned.rs`) | a face pair's prism; a face grouping round; a turned part pair |
| Validation (`topology/validate*`, `jet.rs`) | an integral's piece and box; a quotient integral's and a Green integral's piece; a Green integral's fin; a face's integrals and flux; a spline deviation's piece; a ray's face; a ray crossing's piece |
| Tessellation (`tessellation.rs`, `tessellation/spline.rs`) | a tessellation round; a spline cell; an edge's discretisation; a spline edge's doubling; a face's mesh; a mesh grid; a mesh triangle |
| Intersections (`intersection/*`) | a root's piece; a curve root's piece; an edge's piece; a turning-point box; a tracing step; a Newton step |

## Evidence

**Budget fixtures.**
- The slowest distinct inputs of three corpora at `50afeef4` are kept under
  `rust/fixtures/budgets/`:
  - `boolean`: eight inputs, 3.1 to 4.5 s each in a release replay;
  - `split`: eight inputs, 1.3 to 1.4 s each;
  - `brep_validation`: four inputs, 1.0 to 1.2 s each.
- `rust/fixtures/budget-work.tsv` lists them with those seconds and pins the
  work units of every metered operation of each:
  - `boolean`: fuse, cut and common, 1,669 to 4,769 units;
  - `split`: the split, 454 to 2,293 units;
  - `brep_validation`: each validation, 10,884 to 17,027 units.
- `rust/fuzz/examples/replay_budgets.rs --check` runs them with their
  targets' oracles and S10's invariant. It passes on macOS and on a Linux
  host (glibc 2.39, rustc 1.99.1) with the same counts. The fuzz
  workflow's `boolean` job runs it on CI's Ubuntu runner.
- `tests/budgets.rs` pins two kernel counts, which every kernel job checks
  on its host:
  - two boxes' fuse: 477 units;
  - the tilted torus and a rod crossing its inner equator by `1e-6`:
    51,200 units.

**The fuzz invariant** ([FUZZING.md](FUZZING.md), "Budgets (S10)").
- The `boolean`, `split`, `brep_validation` and `tessellation` targets each
  run one long operation an input twice:
  1. under a budget drawn from the input's hash, log-uniform below `2^24`;
  2. unlimited.
- The two runs must agree in one of two ways:
  - work within the budget: the same result and the same work;
  - otherwise: refused for work with exactly one unit past the budget.
- Replayed over the local corpora and regressions on the Linux host, every
  check holds:

| Target | Checked | Within the budget | of them, no work | Refused | Budget equal to the work |
| --- | --- | --- | --- | --- | --- |
| `boolean` | 1,482 | 950 | 199 | 532 | 7 |
| `split` | 3,378 | 2,653 | 756 | 725 | 38 |
| `brep_validation` | 3,170 | 2,599 | 327 | 571 | 41 |
| `tessellation` | 1,289 | 855 | 0 | 434 | 1 |

- The same replays with debug assertions on macOS, one process an input,
  have no failure. They cover:
  - `boolean`: 1,482 inputs;
  - `split`: 3,571 inputs;
  - `brep_validation`: 3,217 inputs;
  - `tessellation`: 1,315 inputs.

**Cancellation from another thread.**
`another_thread_cancels_a_long_boolean_promptly` (`tests/budgets.rs`) fuses
the tilted torus with the `1e-6` rod. That Boolean takes about a second in
release. A second thread cancels it after 50 ms, and the operation must
return `Cause::Cancelled` within 500 ms of the cancellation. On the
development Mac it returned in 11 ms in release and in 21 ms at
`opt-level = 2` with debug assertions. Afterwards the same fuse, unlimited,
evaluates on the same thread: the cancelled arrangement was not kept.

**The DRAW adapter's `cpulimit`.** On the Rust backend the worker installs
a wall budget with that deadline for every later command (`draw_bridge.tcl`,
`draw_worker.rs`; [UPSTREAM_TESTS.md](UPSTREAM_TESTS.md)'s "Deliberate
limits"). `test_draw_bridge.py` checks two things:
- commands under it evaluate;
- `cpulimit 0` stays unsupported.

The registered cases take their groups' limits of 60 s or more, so the 30 s
contract and the ledger are unchanged.

## Measurements

**The longest stretch between two checks** (`Budget::run`'s start, every
unit and its end). It bounds how late a deadline or a cancellation lands.
It was measured with an uncommitted instrumented build, in release, on the
Linux host, over every corpus input and regression:

| Target | Median | 90th percentile | 99th percentile | Longest |
| --- | --- | --- | --- | --- |
| `boolean` | 2.4 ms | 16 ms | 66 ms | 114 ms (a stack's cells after its profiles' loops) |
| `split` | 0.6 ms | 3.6 ms | 54 ms | 58 ms (a spline integral's setup) |
| `brep_validation` | 0.1 ms | 0.9 ms | 3.3 ms | 11 ms |
| `tessellation` | 0.4 ms | 19 ms | 67 ms | 167 ms (after a grid face's mesh) |

Every stretch is under a fifth of a second. Before the validation's
integrals had units, the boolean corpus's longest stretch on the development
Mac was 1.1 s, from a result's last integral box to the operation's end.

**Cost.**
- Without a budget a unit is a thread-local read. With an unlimited budget
  installed it is a thread-local slow path and an atomic add.
- The boolean corpus counts 4.0 million units in its fuse, cut and common
  (1,482 inputs, 14,784 at most for one input), against about 750 s of CPU.
- The cost was measured on the development Mac (load about 7.5, other
  work), with release replays of the boolean corpus, back to back and
  twice. The three variants are:
  - `50afeef4`;
  - S10 with the old harness (no budget anywhere);
  - S10 with every metered operation run once under an unlimited budget
    (every unit counted).

| Run | `50afeef4` | S10, no budget | S10, counting |
| --- | --- | --- | --- |
| first | 768.5 s | 754.4 s | 754.3 s |
| second | 752.6 s | 751.6 s | 753.5 s |

- The cost is within the runs' noise. The second round, the quieter one,
  shows -0.1% without a budget and +0.1% counting, well under the 1%
  allowed.
- The release suite is also within noise:
  - wall: 726.9 s becomes 717.4 s;
  - CPU: 2,299.5 s becomes 2,280.0 s;
  - tests: 857 become 861, with `tests/budgets.rs`.

## Open items

- **For the user:** the refusal is a new variant, `Error::BudgetExhausted`,
  rather than `ComputationLimit` with a payload ("Refined decisions"). Both
  are mechanical to change while nothing outside the repository uses them.
- Memory budgets are not counted. The 2 GiB fuzz limit and the kernel's
  fixed size limits remain the only bounds on memory.
- Parsing `.brep` and Part 21 text counts no units.
- The longest stretches left are a stack's cells after its profiles'
  loops, a spline integral's setup and a tessellation's checks after its
  last triangle, each under 0.2 s on the Linux host. Points can be added
  there if an application needs faster cancellation. Each addition changes
  the pinned counts and must be recorded again (`replay_budgets --write`,
  `tests/budgets.rs`).
