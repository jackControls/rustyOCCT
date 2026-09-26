# Review notes and next steps

Reviewer's findings, decisions and the recommended order of work, written for
the implementing agent to read and execute. It complements
`IDENTITY_AND_HISTORY.md` and `TOPOLOGY_MODEL.md`; it does not replace their
contracts or acceptance records. When a step below is accepted, record the
evidence in the guide that owns it and mark the step here as done.

**Reviewed state:** head `3ba8c1f2` on 2026-09-26. M0–M5 and T1–T2 are
implemented, pushed and accepted; both workflows are green at `0913b5e2`; the
full local release suite passes at the head.

## Findings

* **F1. One untriaged fuzz timeout.** The push run at `3ba8c1f2` (run
  36236430909) failed its `surface_editing` job with artifacts
  `timeout-5cf6cbdd8250eeb22f43e0502d33c014522c1a2b` and
  `slow-unit-f240e0e8b3821b28ef9edb4177bac843e1213cdd`. The input is not a
  checked-in seed; it came from the evolving CI corpus. It completes: 4.0
  seconds locally on the development Mac without a sanitizer. Under
  AddressSanitizer on the Linux runner, which the notes measure at about 2.6
  times slower, it crossed the 20-second per-input limit. The later green run
  at that commit was the scheduled campaign, not a fix, and nothing in the
  repository records a triage. `FUZZING.md` requires one.
* **F2. T2 captured native observations after implementation.** Every other
  milestone captured first. The T2 notes record the break candidly and
  explain what it does and does not affect; the independent reader was
  written by the same author after debugging the Rust reader. Acceptable
  because disclosed; it must stay unique.
* **F3. Retained upstream evidence is measured and currently zero.** The
  ledger labels 35,766 assertions: 22,920 model-independent, 0
  mapped-and-verified, 12,846 lost, of which 8,953 are lost only because no
  original case yet exercises an existing mapping natively. Only 177 of the
  17,879 cases do nothing but restore a data file and run checks (160 in
  `bugs`, 9 in `heal`; validity in 129, properties in 111), and most of those
  files carry splines, cones or spheres. Open Cascade publishes a public
  testing dataset of more than 2,500 shapes on its developer portal, separate
  from the sources, and its test manual states that many cases use data that
  is confidential and available only at Open Cascade.
* **F4. Cross-platform allowances are accumulating.** Two captures now
  tolerate last-bit libm differences between macOS and Linux at `2^-46` of the
  case's size. Each is documented in its own notes; there is no single list.
* **F5. Open items carried in the guides, all still open:** the D11 continuity
  check before any spline edge or face; the H8 rule for composites across
  algorithm levels; `on_modify` exercised only by hand-written histories;
  poles and degenerate edges not representable, so cones, spheres, tori and
  splines cannot be imported; the `.brep` writer records the resolution, not
  each enclosure; no `restore` command in the DRAW adapter.
* **F6. Resolved since the M3 review:** fused parents are ordered axially
  (`1429a803`), so ids no longer depend on call order.

## Decisions

* **R1. F1 is a budget, not a defect.** Re-measure the input locally *with*
  AddressSanitizer. If it exceeds 20 s ÷ 2.6, move `surface_editing` onto the
  60-second per-input budget in `run_fuzz.py` (`TARGET_INPUT_SECONDS`, beside
  `surface_knots` and `degree_elevation`); if it does not, the runner was
  slow and the budget stays. Either way, keep the artifact as a slow-unit
  regression under `rust/fuzz/regressions` with the measured times, and add
  to `FUZZING.md` the rule that a timeout whose input completes is triaged by
  this measurement. Never re-run to green.
* **R2. Test data is fetched by a separate, pinned step, never by the bridge,
  never by CI, and never committed.** `rust/tools/fetch_occt_test_data.py` downloads the
  public dataset into an ignored `target/occt-test-data`, verifies a recorded
  SHA-256, and prints the file inventory. The runner distinguishes
  `not-fetched` from `private-data` (a file absent from the public dataset);
  both are non-passes. Local use only; CI never downloads it (U1).
* **R3. `restore` is a thin call into the T2 reader and converter.** It
  reports unsupported constructs by name as the interop already does; a case
  whose file is unsupported is `Unsupported`, not a failure. Every restored
  and checked case that passes natively moves its assertions to
  mapped-and-verified through the existing ledger rules.
* **R4. Continuity is C1 in the cell's own parameterisation, not G1.** It is
  exactly decidable and a G1-only curve loses nothing by a vertex at the
  knot. A B-spline is C1 at an interior knot of multiplicity below its
  degree by construction; at multiplicity at or above the degree it is C1
  exactly when the kernel's exact knot removal reduces that multiplicity by
  one with zero residual. The same test on the control net along a knot line
  decides each surface direction, and it applies to pcurves. Deliverables:
  `edge_not_c1`, `pcurve_not_c1`, `face_not_c1` in `Topology::check`, always
  certified; an independent `Fraction` check by exact removal; two fuzz
  mutations, a broken knot that must be reported and an exact repeat that
  must stay valid. It is the first deliverable of S4, since nothing in a
  topology can fail it before then.
* **R5. Composites record every step's level.** Replace the single `level`
  of a composite by `steps: Vec<(OperationId, OperationKind, AlgorithmLevel)>`,
  one entry for a single operation, concatenated by `History::then`; a
  composite's level *is* that list. Composites are never replayed, so nothing
  else changes and derivations stay level-free. The checker requires a
  non-empty ordered list. Test it now by composing a hand-written history at
  a synthetic level 2 with one at level 1; no callable second level is
  needed. Closes the H8 open item.
* **R6. Surfaces of revolution come before splines, in the order cone,
  sphere, torus.** The cone's apex is the simplest pole and becomes a vertex
  loop; the sphere is the first zero-loop face; the torus is the first face
  wound in `v`, closing that T1 gap. Each ships with: the `Surface` variant
  with certified evaluation; its pcurve convention on the universal cover;
  validator invariants and `cell_reference.py` rules for pole loops and `v`
  windings; a primitive builder with derivations, roles and a history; the
  count synthesizer's one degenerate edge per pole loop; the converter and
  writer mapping to OCCT's conical, spherical and toroidal surfaces with
  degenerate edges; fixtures, fuzz extensions of `brep_validation`,
  `identity`, `history` and `brep_io`; native bridges on `BRepPrimAPI`
  history, `BRepCheck` and `BRepGProp`; derived DRAW cases for `pcone`,
  `psphere` and `ptorus` counts; docs. General mass properties arrive with
  the cone (U2).
* **R7. Native observations before implementation is mandatory.** The T2
  exception stays the only one. Any future exception needs the user's
  approval before the first line of kernel code.
* **R8. List every platform allowance in one place.** `VALIDATION.md` gains a
  table of the captures that tolerate cross-platform rounding, each with its
  bound and reason, updated whenever one is added.

* **R9. Viewer commands get their own recorded status, on both backends.**
  A local trial with the dataset present ran `bugs/modalg_1/buc60684` on
  native DRAW: `restore`, `checkshape`, `prism` and `checkprops` all ran, and
  the case was then latched `unsupported` because its final `checkview`
  calls `smallview`. 7,428 upstream cases touch the viewer, almost always as
  a final image dump after the geometric assertions. Add a status
  `viewer-skipped` meaning "every geometric assertion was evaluated; the
  image commands were recorded, not run". It is not a pass, it never counts
  a stripped file, and the ledger counts the evaluated assertions exactly as
  it does for a pass. Files are never edited.
* **R10. Data lookup descends one level, as upstream's `locate_data_file`
  does.** The bridge's `locateData` searches only the directories given, so
  the dataset's `brep`, `geom`, `iges`, `msv`, `others`, `step` and `xbf`
  subdirectories had to be passed one by one. Match upstream: search each
  `--data-dir` and its immediate subdirectories, so `--data-dir` takes the
  dataset root.

## Next steps, in order

### S1 — small closures (no new geometry)

R1 triage, R5 composite steps, R8 allowance table. Gate: kernel and fuzz CI
green at one revision; `history_contracts.rs` covers the two-level
composition; the regression and the measurements are committed.

### S2 — `restore` and the public test data

Local trial on 2026-09-26 with the dataset in `target/occt-test-data`:
`buc60684` resolves its file once the subdirectories are given (R10); native
DRAW runs every geometric assertion and is latched by `smallview` (R9); the
Rust adapter reports `unsupported DRAW signature: restore ...`, which R3
removes. So S2's order is R10, R9, R3, the two runner statuses, the fetch
script, then the ledger.

R2 fetch script, R3 `restore`, the two new runner statuses, and the ledger
re-recorded. Report the count of restore-only cases that became
mapped-and-verified and the count reported `Unsupported` by construct name,
so S3 can be prioritised by what the corpus actually contains. Gate: the
original-test bridge green on the Rust adapter and native DRAW, ledger
recorded, no automatic download anywhere in the bridge.

### S3 — cones, spheres and tori in the cell model

R6, one surface per accepted revision, each with the full shipping pattern
and its own acceptance record in `TOPOLOGY_MODEL.md`. General mass
properties land with the cone: surface integrals with certified quadrature
and an error bound that feeds the enclosure contract (U2). Expected upstream
effect: `pcylinder`, `pcone`, `psphere`, `ptorus` and `box` are the most
common constructions in self-contained cases, so count and property
assertions start moving to mapped-and-verified.

### S4 — the continuity check, then spline edges and faces

R4 first, as its own accepted record, then rational spline edges and faces
in topology with the exact modules that already exist, pcurves on the
universal cover, the converter and writer mapping, and the interop's
unsupported list shrinking accordingly.

### After S4

`PORTING.md` step three: general face trimming, curve/surface intersection
and projection, then Booleans. Not scoped here.

## User decisions (answered 2026-09-26)

* **U1. Test data on CI: yes, if the licence allows; local disk is scarce.**
  Located and fetched on 2026-09-26. The dataset is an official OCCT GitHub
  release asset, not a portal download: `opencascade-dataset-7.9.0.tar.xz`
  under the `V7_9_0_beta1` and `V7_9_0_beta2` releases of
  `Open-Cascade-SAS/OCCT` (earlier versions under `V7_8_0` and `V7_7_0`).
  OCCT's own public CI actions download exactly this file for its test runs,
  and Gentoo's `sci-libs/opencascade` package downloads it for `USE=test`
  under the package licence `Open-CASCADE-LGPL-2.1-Exception-1.0`.
  * URL: `https://github.com/Open-Cascade-SAS/OCCT/releases/download/V7_9_0_beta2/opencascade-dataset-7.9.0.tar.xz`
  * size 98,739,184 bytes; SHA-256
    `a92ed91c3271c299287c1c404bb9454d463251094bfc848acc44aee60e6a026c`;
    SHA-512 equal to Gentoo's Manifest entry
  * extracted: 3,388 files, 346 MB, directories `brep geom iges msv others
    step xbf`; 1,494 `.brep`, 875 `.rle`, 322 `.stp`, 201 `.draw`, 58 `.igs`
  * no `LICENSE` or `README` inside the archive, and no licence sentence in
    the release notes, the announcement or the test manual. The basis for use
    is that it is an asset of an LGPL-2.1-with-exception release, published
    and consumed by upstream's own CI. That is a reasonable basis for
    fetching it to run tests; it is not a basis for committing or
    redistributing any file from it. The user decides whether that basis is
    enough (see the reply of 2026-09-26).
  * coverage against `tests/`: 2,924 of 6,963 literal referenced file names
    are present (42%); 5,600 of 11,556 data-dependent cases have every file
    (48%), 514 some, 5,442 none. By group, fully covered: boolean 1,291 of
    2,078, bugs 1,112 of 3,283, heal 668 of 1,186, offset 593 of 1,786,
    sewing 508 of 700, thrusection 248 of 248, blend 263 of 292, lowalgos
    104 of 215. The remainder is the confidential set the manual describes.
  * older public portal archives (`shapes_7.5.0.tgz`, `opencascade-dataset-7.x`
    under `dev.opencascade.org`) are gone: every such URL now returns the
    portal's HTML page.
  * local copy: `target/occt-test-data/opencascade-dataset-7.9.0` (ignored);
    22 GiB were free locally, so it was not necessary to use `jackgpu`, which
    was on the subnet but not answering.
  **Decided 2026-09-26: local only, no CI.** The fetch script pins this URL,
  size and SHA-256, extracts into `target/occt-test-data` and never commits
  a file; it is run by a developer, never by a workflow. CI continues to
  report every data-dependent case as `not-fetched`, and the ledger records
  the local run separately with the dataset's SHA-256 so a local
  mapped-and-verified count is never mistaken for a CI one. A missing file
  is `private-data` when its name is absent from the archive inventory and
  `not-fetched` when the archive is absent. Revisit CI use only when the
  user says so.
* **U2. Mass properties: certified error bounds.** General mass properties
  arrive with the cone as surface integrals with certified quadrature and a
  bound that feeds the enclosure contract; the prism closed forms stay.
* **U3. Continuity: C1 in the cell's parameterisation**, after the
  explanation below; G1-only joints are handled by exact reparameterisation
  where the speed ratio is rational and by a vertex otherwise.
* **U4. Cadence: continuous.** The agent does not pause for permission
  between milestones; it stops only for a decision listed in this file or a
  new one it records here first.
* **U5. Fuzz budget: approved** for `surface_editing` if R1's sanitizer
  measurement calls for it.

### Why C1 rather than G1

Two curves meeting at a knot are G1 when their unit tangents agree: the shape
has no corner. They are C1 when the derivative vectors agree, direction and
magnitude, in the shared parameter: the parameter speed does not jump. G1 is
the property the geometry needs; C1 is a property of the parameterisation.
Every C1 joint is G1; a G1 joint whose speeds differ by a factor is not C1.

The kernel requires C1 for three reasons. It is exactly decidable: for a
rational B-spline, C1 at an interior knot is precisely "the exact knot
removal already in the kernel reduces the knot's multiplicity by one with
zero residual", a rational computation with no tolerance, while G1 asks
whether two derivative vectors are positive multiples of each other, which
still needs the same machinery plus a ratio. It costs nothing geometrically:
a G1 joint with a rational speed ratio becomes C1 by scaling the knot
interval on one side exactly, which the exact knot editing supports, and an
irrational ratio, which imports cannot produce from binary64 data, would get
a vertex, adding an entity but no shape. And it keeps every parameter-based
algorithm honest: Newton steps, projections and arc-length integrals that the
intersection and proximity modules already use assume derivative continuity
within a span sequence, so a C1 cell is the contract they were written for.
The cost is that imported curves joined with different speeds show one more
vertex than in the source system when the ratio is not rational; the
converter records that as `Imported` provenance, never as an approximation.

## Status

* S1 — implemented at `7e463cb2`; gate pending CI. R1: the CI input takes
  8.0–8.2 s locally under AddressSanitizer (4.1 s without), above
  20 s ÷ 2.6, so `surface_editing` has the 60-second budget; both inputs are
  regressions with their times in `fuzz/regressions/README.md`, and
  `FUZZING.md` states the triage rule. R5: `History::steps`, the
  `steps_invalid` check and the two-level test in `history_contracts.rs`.
  R8: the table in `VALIDATION.md`.
* S2 — implemented; gate pending CI.
  * R10 deviation: upstream's `locate_data_file` searches every
    subdirectory level breadth-first (skipping dot-directories), not one
    level. The bridge now matches upstream exactly, including the case's own
    `data` folder and the case variables of `_run_test` (`casename`,
    `imagedir` and the rest, whose absence failed two cases).
  * R9 and R3 as decided. `restore` also reports a saved curve or surface
    (another DRAW object type) as `unsupported`. The reader accepts a root
    location glued to stray bytes and the worker decodes leniently, as OCCT
    reads nothing after the root: two upstream files need it.
  * Survey of the 192 restore-only cases with the dataset: native evaluates
    33, Rust 2, both registered (`bug27264_1` pass, `buc60769`
    viewer-skipped); `buc60684` is native viewer-skipped, Rust unsupported
    (`prism` without `Copy`). The ledger has its first mapped-and-verified
    assertion (`bug27264_1`'s `checknbshapes`). 103 cases need private
    data.
  * For S3's priority, what the 31 natively evaluated cases Rust cannot yet
    run need (a case may need several): free faces 22, B-spline curves on
    surfaces 18, B-spline curves 16, B-spline surfaces 13, rectangular
    trimmed surfaces 9, Bézier surfaces 3, tori 3, cones 2, spheres 2,
    extrusion surfaces 2; two need Booleans. Free faces and B-spline
    surfaces dominate the data-dependent corpus. Cones, spheres and tori
    remain first under R6 because `pcone`, `psphere` and `ptorus` dominate
    the self-contained cases; the data argues for free faces (sheet bodies)
    early in S4.
  * U1 is not fully resolved: CI does not fetch the dataset until the user
    confirms the licence basis recorded under U1. Locally the fetch script
    verifies the pinned archive; data cases report `not_fetched` on CI,
    which the contract accepts.
* S3 — pending
* S4 — pending
