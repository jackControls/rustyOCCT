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
