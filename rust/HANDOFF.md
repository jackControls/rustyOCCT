# Handoff: where to pick up

**R15, the green-gate rule** (`REVIEW_2026-10-18.md`): a red workflow is
the next task, always. No "done", "verified" or "accepted" bullet may be
written while either workflow is red at that revision. This header states
the CI state of the head commit. A push that turns a workflow red is
followed only by commits that address it or revert it.

**CI state of the head** (as of the branch point of X4 and X5):

- CI at `36f83779`: green. "Rust kernel" run 39190992103 (Ubuntu,
  Windows, macOS and Rust 1.85), "Rust geometry fuzzing" run 39190991866,
  "Rust Boolean captures on Linux" run 39197100462.
- CI at `50afeef4` (the head of `rust-kernel`: X6's CI changes, merged at
  `5e861738`): pending. Fill in each workflow's run and result here when
  confirmed; a later head adds its own line.

## Pick up

1. Read [STATUS.md](STATUS.md): what exists and at which revision, the
   open items in order, the decisions answered, how to run the gates.
2. Confirm the head's CI above (`gh run list -R jackControls/rustyOCCT
   -L 6`); a red workflow comes before anything else (R15).
3. X6's open gate, a week of scheduled runs within budget: the review's
   Status and [FUZZING.md](FUZZING.md)'s sharded minimisation and CI wall
   times.
4. Then the next step, S10 (cancellation and work budgets), and S11 to
   S15 in U11's order: the review's "Next development steps".
5. How a step is done, verified and committed:
   [CONTRIBUTING-KERNEL.md](CONTRIBUTING-KERNEL.md).
6. Contracts and evidence per capability: [BOOLEAN.md](BOOLEAN.md) (its
   comparisons, Linux records and open items), [SPLIT.md](SPLIT.md),
   [TESSELLATION.md](TESSELLATION.md), [VALIDATION.md](VALIDATION.md),
   [UPSTREAM_TESTS.md](UPSTREAM_TESTS.md), [FUZZING.md](FUZZING.md),
   [IDENTITY_AND_HISTORY.md](IDENTITY_AND_HISTORY.md),
   [TOPOLOGY_MODEL.md](TOPOLOGY_MODEL.md) and
   [MATHEMATICS.md](MATHEMATICS.md).
7. The history of decisions and evidence up to the review of 18 October:
   [REVIEW_NOTES.md](REVIEW_NOTES.md), frozen (no new bullets).

## This page's former items

Code comments, tools and the Linux workflow cite items of this page as it
stood up to `50afeef4` (`git show 50afeef4:rust/HANDOFF.md`). Where their
content is now:

- item 2, the Boolean captures' Linux records and the steps to take one:
  BOOLEAN.md's "Linux records (CI)"; the comparisons' table: its
  "Comparisons";
- item 3, S9e.4b's sub-steps and the version-3 `.brep` header:
  BOOLEAN.md's sections and "Open items and refusals", REVIEW_NOTES.md's
  bullets per sub-step and "The version-3 `.brep` header, done";
- item 5, the near misses, (d) the shallow lens among them: BOOLEAN.md's
  near-miss paragraphs, REVIEW_NOTES.md's "The shallow lens, done";
- item 6, the torus meeting's verification budget and the tilted torus's
  validation: REVIEW_NOTES.md's bullets of those names (S9);
- "How a step is done", "Verification", "Working rules and lessons":
  CONTRIBUTING-KERNEL.md; "User decisions answered": STATUS.md.
