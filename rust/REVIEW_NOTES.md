# Review notes and next steps

Reviewer's findings, decisions and the recommended order of work, written for
the implementing agent to read and execute. It complements
`IDENTITY_AND_HISTORY.md` and `TOPOLOGY_MODEL.md`; it does not replace their
contracts or acceptance records. When a step below is accepted, record the
evidence in the guide that owns it and mark the step here as done.

**Reviewed state:** head `3ba8c1f2` on 2026-09-26. S1–S4 are done
(2026-09-27, see Status); what remains open for the user is listed there. M0–M5 and T1–T2 are
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
* **R11. The cone's history capture follows its builder (recorded
  2026-09-26 by the implementing agent; accepted by the user on
  2026-09-27).** The cone capture at `dde086c1` observed `BRepPrimAPI_MakeCone`'s
  structure, `BRepCheck` and `BRepGProp` before any kernel cone code, but no
  history: `MakeCone` reports none and DRAW `pcone` saves none. The history
  the kernel reports for a cone is that of revolving its meridian, which
  OCCT gives through `BRepPrimAPI_MakeRevol`, and that was not captured
  before the kernel's cone builder and its derivations were written
  (uncommitted at the time). The `MakeRevol` capture is taken before any of
  that is committed and before the Python enumeration and the role
  correspondences exist; its metadata says the builder existed. The sphere
  and the torus capture history with everything else, before any kernel
  code. If the user does not accept this, the cone's history bridge is
  labelled like T2's (captured after implementation) and nothing else
  changes.
* **R12. Retained fuzz corpora grow without bound (recorded 2026-09-26 by
  the implementing agent; the user chose weekly minimisation on
  2026-09-27, see Status).** The
  push campaign at `0c94aa53` failed `surface_knots` because its corpus
  replay no longer fit the 3,600-second startup cap: 2,766 s for 564 inputs
  at `1a76d29e`, over 3,600 s for 577 at `0c94aa53`, with no input over its
  60-second limit (slowest 32 s). Nothing in `surface_knots` changed; the CI
  cache retains every input libFuzzer adds. Interim, as R1 did for
  `surface_editing`: `surface_knots` gets a 7,200-second startup cap
  (`TARGET_MAX_STARTUP_SECONDS`) and the fuzz job 195 minutes. The lasting
  fix is to minimize retained corpora (`-merge=1` into a fresh directory,
  for example weekly) so replay stays bounded; it costs one more replay per
  minimization and drops inputs that add no coverage. The same push also
  showed the cone checks slowing `identity`, `history` and `brep_io`
  replay (about 2 s per input under Linux AddressSanitizer). That was the
  cone's mass engine, not a budget: it is now 20 to 40 times faster
  (shared trigonometric moments, memoized exact dyadic Fourier
  coefficients).

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

Decisions for R4, recorded before its code (2026-09-26):

* **Representation first.** R4 adds `Curve3::BSpline` (a `BSplineCurve3`
  over its whole domain, the edge fraction mapped affinely onto it),
  `Curve2::BSpline` (a new `BSplineCurve2`, the planar counterpart, for
  pcurves) and `Surface::BSpline` (a `BSplineSurface3`, `u` or `v` periodic
  as its basis is). Nothing builds or imports them yet: the reader still
  reports B-spline records unsupported until the rest of S4.
* **The criterion is homogeneous.** A binary64 B-spline has interior
  multiplicities at most its degree, so a knot of lower multiplicity is C1
  by construction and only a knot of multiplicity equal to the degree
  (degree 1 included) is tested: one exact homogeneous removal
  (`ExactBSplineCurve3::remove_knot`, the surface's `remove_u_knot` /
  `remove_v_knot` for a whole knot line) with zero residual. That is C1 of
  the homogeneous curve, which implies C1 of the rational one; a rational
  joint whose weights alone jump is reported, as R4 states. A periodic
  spline's seam knot is interior and tested the same way; a spline ring
  edge must be periodic (`ring_edge_open` otherwise). Every knot of the
  curve or surface is tested, trimmed away or not: the cell's
  parameterisation is the whole curve or surface.
* **Until the rest of S4, spline geometry is never called valid by the
  other geometric checks.** Each check that cannot yet certify spline
  geometry reports its `uncertified_*` kind (vertex on curve, pcurve on
  edge, UV gap, loop winding, containment, shell orientation), and mass
  properties are `None`. The continuity issues are exact, so R4's fixtures
  and fuzz mutations require them exactly and leave the rest of the report
  unchanged by the mutation.
* **The independent check does not remove knots.** The reference compares
  the homogeneous curve's left and right first derivatives at each tested
  knot in `Fraction`s (de Boor on each side), the equivalent condition.

Decisions for the rest of S4, recorded before its code (2026-09-26):

* **Order.** S4a native capture, before any kernel code that certifies
  spline geometry: the native construction protocol and
  `occt_brep_check_oracle.cpp` learn B-spline curves, 2D curves and
  surfaces; spline cases with OCCT rows (planar faces bounded by spline
  edges, ruled spline walls, spline pcurves on planes and cylinders, and
  their mutations) record BRepCheck verdicts, tolerances, per-use
  deviations and `BRepGProp` properties. S4b certified deviation, UV gaps
  and vertex gaps for spline uses. S4c loop winding and containment with
  spline pcurves. S4d face flux and mass properties with spline pcurves and
  surfaces. S4e interop (B-spline records in the reader, converter and
  writer; periodic spline surfaces' windings with their domain's period).
  Each lands with its reference, fixtures, fuzz coverage and docs; S4 is
  accepted once as a whole, like S3.
* **Deviation (S4b).** Where the composition is rational (a plane or spline
  surface, a line or spline pcurve, a line or spline edge), the kernel forms
  the homogeneous difference `N = C_h w_SP - (S∘P)_h w_C` exactly on every
  common Bézier piece (split at the knots of both curves and at the
  pcurve's crossings of the surface's knot lines) and bounds `|D| <= max
  |N_i| / min (w_C w_SP)_j` from Bernstein coefficients: exact cancellation,
  no subdivision. Otherwise (an arc edge or pcurve, or a periodic analytic
  surface under a spline pcurve) it bounds `D` by a second-order Taylor
  enclosure, `|D(t_m)| + |D'(t_m)| h/2 + sup|D''| h²/8`, over pieces
  refined until it is within tolerance, with the hulls of exact Bézier
  pieces for the spline factors and the certified trigonometry for the
  analytic ones. A failure is certified by an exactly evaluated sample, as
  today. The reference is independent: exact `Fraction` composition in the
  rational case, `mpmath` interval arithmetic otherwise.
* **UV gaps on spline surfaces** are measured in 3D, `|S(a) - S(b)|` from
  exact evaluations, instead of by a length scale.
* **Spline edges and pcurves carry a range (decided 2026-09-27, before
  S4e).** In the `data/occ` corpus 999 spline edges use their curve's whole
  domain, but 200 use a sub-range and 251 more are trimmed curves; trimming
  a spline exactly yields poles binary64 cannot hold. So `Curve3::BSpline`
  and `Curve2::BSpline` hold a `SplineSpan`, the curve and a closed range
  inside its domain (the whole domain when the kernel builds it), the
  fraction mapped affinely onto the range. R4 then tests the knots strictly
  inside the range: the cell's parameterisation is the range. A surface is
  still tested over its whole domain.
* **Periodic spline surfaces stay unsupported (decided 2026-09-27, during
  S4e).** Their windings with the domain's period were planned for S4e,
  but in the `data/occ` corpus one solid of 77 has a periodic spline surface
  (`PeriodicBSplineSurface`) and one a seam on a nonperiodic closed spline
  surface (`SeamOnBSplineSurface`); every other spline surface (419 of 420
  records) is nonperiodic and unseamed. The reader reports both by name, the
  validator does not certify a periodic spline face, and windings on spline
  surfaces wait for a case that needs them.
* **Mass properties with spline geometry (S4d, decided 2026-09-27, before
  its code).** The S4a capture has no `BRepGProp` row, so native properties
  come first: `compare_brep.py --family spline --capture-properties`
  records OCCT's volume and area (each with its error estimate at a
  requested relative precision of `1e-12`), centre of mass and inertia
  about it for every spline model (`fixtures/occt-spline-properties`),
  before any kernel code integrates a spline surface or a spline pcurve on
  a cylinder or cone; every later run must reproduce them. Then:
  * along a spline pcurve on a cylinder or cone, each term's `-∫ F du` by
    the Green integral already used on spheres and tori;
  * on a nonrational spline surface whose pcurve pieces each lie in one
    patch, the ten volume and moment integrands are tensor Bernstein
    polynomials of the patch (`p` and `N = X_ū × X_v̄` from the poles),
    integrated exactly by the flux's column antiderivatives; the four
    terms with `|N|` (area and face centre) are enclosed by the flux's
    strips;
  * on any other nonperiodic spline surface, all fourteen terms by the
    strips.
  The Green integral and the strips are first order, so those enclosures
  are wide (their widths are recorded); a higher-order certified quadrature
  is later work, not S4. The bridge requires every kernel enclosure of a
  spline model valid on both sides to contain OCCT's value up to OCCT's own
  error estimate and `1e-9` relative (amended after the capture to `1e-8`:
  OCCT's area of `spline_bulge` is `1.7e-9` relative from its closed form
  while estimating `2e-16`, `occt-spline-properties/NOTES.md`). The reference integrates the same
  properties of every valid spline fixture case by nested Gauss–Legendre
  quadrature (`cell_reference.py`), and `brep_validation.rs` requires the
  kernel's enclosures to contain them.

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

## Review of S1–S4 (2026-09-27, head `a7301d68`)

All four steps are done and recorded, with captures before code except the
cone's `MakeRevol` history (R11, accepted by the user). Both workflows are
green through `cde6d254`; the head's runs are queued; the full local release
suite passes at the head. The kernel is 29.5k lines with 14.7k lines of
tests, 20.4k lines of Python references and bridges, 162 fixture files and
24 fuzz targets. Cones, spheres and tori exist with poles, zero-loop faces
and `v` windings; spline cells exist with the C1 check; general mass
properties are certified; 60 of the 77 `data/occ` solids import (58
certified); the ledger has its first 2 mapped-and-verified assertions.

Findings that shape the next moves:

* **F7. CI time is now dominated by corpus replay.** The fuzz workflow takes
  74–82 minutes per push; `surface_knots` alone replays 614 retained inputs
  for 4,286 s of its 4,422 s and mutates for 60 s. The kernel workflow takes
  64 minutes because of the Windows job. Weekly minimisation (R12) helps
  but cannot change the ratio: replay of exact tensor inputs is the cost.
* **F8. Spline mass enclosures are first order.** Areas of nonrational
  spline walls are enclosed within 0.2–0.6%, the rational corner's volume
  within 5%. Recorded honestly and decided as later work; not yet
  production-grade for measurement. (Addressed 2026-09-28 by the certified
  quadrature of the parallel track: within `1e-12` of each property's
  scale on every spline fixture.)
* **F9. Profiles still accept only polygons and circles.** Arcs exist in
  the topology and the validator's fixtures but not in `Boundary`, so the
  application's first job, Extrude of a sketch with arcs and fillets, cannot
  run.
* **F10. The data corpus is blocked mostly by sheet bodies.** The S2 survey
  of restore-only cases native evaluates but Rust cannot: free faces 22,
  spline curves on surfaces 18, spline curves 16, spline surfaces 13,
  trimmed surfaces 9. S4 removed the spline blockers; free faces, that is
  bodies without a solid region, are the largest remaining one, and the
  model already allows them (D1, D9).

## Next moves (2026-09-27)

Ordered by dependency and value. Each ships the full pattern and records its
acceptance in the guide that owns it.

### S5 — arcs in profiles and prisms

`Boundary` gains circular-arc segments beside lines, tangent or not, for
outlines and holes; prism walls on arcs are partial cylinder faces with two
vertical edges (no wrap), caps get arc edges. Ids and history keep the
segment roles; the interop maps arcs both ways; `profile`-built cases in
the DRAW adapter; the `MakePrism` bridge extended; the neutral generator's
existing arc prisms become builder fixtures. Value: the Extrude job takes
real sketches. Gate: polygon and circle ids byte-identical again.

Decisions for S5, recorded before its code (2026-09-27):

* **Input.** `Boundary::path(points, segments, tolerance)`: segment `i` runs
  from point `i` to point `i + 1` (cyclically) and is a `Segment::Line` or a
  `Segment::Arc { center, radius, ccw }` turning counter-clockwise or
  clockwise about the profile normal. Both ends must lie within tolerance of
  the arc's circle, the sweep is the turn from the first end's direction to
  the second's, in `(0, 2π)`, and a path has at least two segments. A path
  of lines only is exactly `Boundary::polygon` (same stored points, same
  ids). The radius is an input, so most arc decisions keep the exact
  sum-of-terms predicates of `decide.rs`.
* **Orientation and labels.** The stored path runs counter-clockwise by its
  signed area, arcs' bulges included; a clockwise input is reversed as a
  polygon is (points `1..` reversed, each segment traversed backwards, arcs'
  directions flipped) and its labels map the same way.
* **Validity.** Every segment's chord and every radius exceed the tolerance.
  Non-adjacent segments stay farther apart than the tolerance, decided
  exactly where the distance is a sum of input terms and in rational
  intervals otherwise, an undecided case counting as touching (a validity
  screen, as `area_is_degenerate`). Adjacent segments may not double back
  at their shared point, may not meet again, and neither's far end may lie
  within tolerance of the other (the polygon rule).
* **Prisms.** An arc segment gives circular-arc bottom and top edges, the
  vertical lines at its ends and a partial cylinder wall (one unwound loop
  of line pcurves on the cylinder's cover, sensed so its normal leaves the
  material), with the same roles and ids a line segment gets; caps take
  circular-arc pcurves. Profile moments add each arc's circular segment in
  closed form.
* **Evidence first.** The identity reference learns arc segments, the
  native `MakePrism` probe learns arc wires, and arc prisms (fillets,
  notches, arcs over half a turn, two-arc lenses, holes with arcs, reversed
  input, labels, tilted frames, transforms) are captured natively
  (`fixtures/occt-arc-prism-preimplementation`) before any kernel arc code;
  `compare_history.py --family arc` compares them afterwards.

### S6 — sheet and wire bodies

Bodies without a solid region: a planar face from a profile, a face on any
supported surface with loops, wires. `Topology::check` classifies them
(D9); the history treats a face body like any other; the `.brep` converter
imports free faces and shells, which unblocks the largest group of data
cases; `mkface`, `mkplane` and `mkedge` in the adapter. Value: ledger, and
the tools that later operations need (a face as a splitting tool).

Decisions for S6, recorded before its code (2026-09-27):

* **Classes (D9).** `Topology::class()` computes `Solid` (every bounded
  region solid, every face between a solid and a void, no wire edges or
  acorn vertices), `Sheet` (faces only, no solid region: each face between
  void regions, possibly the same one), `Wire` (wire edges only), `Acorn`
  (acorn vertices only) or `General`. Nothing stores it.
* **Sheets.** An open sheet's faces have both sides in one shell of the
  infinite void: the radial alternation holds with one fin per boundary
  edge, so an open boundary is not a free edge, and there is no
  orientation flux to decide. A closed shell without a solid bounds a
  bounded void region, oriented and contained as a solid's shell is. Every
  face-level check (loops, pcurves, gaps, windings, continuity, enclosures)
  is unchanged. A shell is still one connected component; a body of
  several components is several bodies (as import makes one body per
  solid).
* **Wires and acorns.** A wire body is one shell of the infinite void
  listing connected wire edges (no fins), each with its curve's checks; an
  acorn body lists one vertex. Their vertex gaps are certified as edges'.
* **Builders.** `Body::face_from_profile` (a planar face with the profile's
  loops) and `Body::wire_from_boundary` (the boundary's edges and
  vertices), with new operation kinds `MakeFace` and `MakeWire` and roles
  `Face`, `Edge` and `Vertex` from the profile's elements, in the identity
  encoding and its independent reference; rigid motions keep their ids.
  Mass properties of a sheet are its area and centre, of a wire its length
  and centre.
* **Interop.** The reader turns free faces and shells into sheet bodies,
  free wires and edges into wire bodies and free vertices into acorns; the
  writer writes them back as OCCT faces, shells, wires and vertices;
  synthesized counts follow OCCT's: a sheet of one face counts as a free
  face (no shell), of several as one shell; a wire of one edge as a free
  edge (no wire), of several as one wire (so an imported one-face shell or
  one-edge wire counts differently, a reviewed difference if it occurs).
* **Evidence first.** The independent validator (`cell_reference.py`)
  learns sheets, wires and acorns, the generator adds such cases and their
  mutations, and native `BRepCheck`, counts and `BRepGProp` properties of
  them, and the corpus's free shapes read natively, are captured before any
  kernel code accepts them.
* **DRAW.** `plane`, `cylinder`, `line` and `circle` objects, `mkface` and
  `mkedge` on them, `mkplane` of profile faces with arcs, and
  `checkshape`, `nbshapes`, `sprops` and `lprops` of sheets and wires.

### S7 — intersections of the analytic family (PORTING step 3)

Curve/curve, curve/surface and surface/surface for planes, cylinders,
cones, spheres and tori, certified, with the existing exact modules
extended: plane/quadric intersections are exact conics; quadric/quadric
intersections are space curves that are not splines. Decision U7 below
chooses their representation. Every result carries enclosures; every
degeneracy (tangency, coincidence, containment) is a declared case with a
fixture; native `IntTools`/`GeomAPI` bridges and the `lowalgos` upstream
group. This is the mathematical core Booleans stand on.

Decisions for S7, recorded before its code (2026-09-27):

* **Order.** S7a: surface pairs whose intersection is a point, lines or a
  conic (plane/plane, plane with a cylinder, cone or sphere, sphere/sphere,
  and the quadric pairs `IntAna_QuadQuadGeo` solves in closed form:
  coaxial, parallel-axis and common-apex configurations). S7b: the other
  quadric pairs and every torus pair, as procedural curves (D13). S7c:
  curve/surface for lines, circles and the S7a conics against every
  analytic surface. S7d: curve/curve. The `lowalgos` upstream group closes
  the step.
* **Exact geometry, exact degeneracy.** A surface is the exact point set its
  stored binary64 data define: a plane through the stored origin with the
  stored normal, a cylinder, cone, sphere or torus about the line through
  the stored origin along the stored normal, normalised exactly (the
  projector `n n^T / |n|^2` is rational), with the stored radii and angle.
  (The validator's distances use the stored normal as if unit; the two
  differ by rounding, far inside any tolerance, but only the normalised
  sets are the intended quadrics exactly.) A cone's `cos a` and `sin a` are
  transcendental for every nonzero binary64 angle, so a plane is never
  exactly parallel to a generatrix: an exact parabola never occurs, and a
  plane/cone intersection is an ellipse, a hyperbola, lines through a
  rational apex, a point or a circle. A configuration is degenerate (parallel, coaxial, tangent, coincident,
  through an apex) only when it is exactly: every such decision is an exact
  predicate on the stored data, as in `decide.rs`. A near-degenerate input
  gets the generic result (a very eccentric ellipse, two very close lines)
  with its enclosures. Merging within a tolerance belongs to the face
  algorithms of S8 and S9, which have bounded faces and a resolution.
  `IntAna_QuadQuadGeo` snaps with angular and linear tolerances; where it
  does, the difference is a reviewed native difference.
* **Results.** Surfaces are the topology's `Surface` values. A result is
  empty, coincident, a finite set of points, or curves: lines, circles,
  ellipses, parabolas and hyperbolas (S7a), procedural curves (S7b). Each
  curve's canonical parameters (a point and unit directions, radii or
  semi-axes, a focal distance) are returned in binary64 with certified
  intervals containing the exact values; points likewise. Canonical forms
  fix every sign and ordering (for example an ellipse's major axis first,
  its direction's first nonzero coordinate positive), so an independent
  reference can check the intervals.
* **Evidence first.** An independent reference (`analytic_intersection_
  reference.py`, mpmath at 60 digits from the same binary64 inputs)
  classifies every fixture case and computes its canonical parameters; a
  native `IntAna_QuadQuadGeo` and `GeomAPI_IntSS` capture of the same cases
  is taken before any kernel intersection code. Fixtures declare every
  degeneracy class with at least one exact and one near case.

Decisions for S7b, recorded before its code (2026-09-27):

* **Ruled parameterisation.** Every quadric pair S7a leaves contains a
  cylinder or a cone. The curve is parameterised on it (the thinner
  cylinder, else the cylinder, else the first cone) by the angle `u` of a
  ruling in an exactly orthonormal frame about its axis: substituting the
  ruling into the other quadric gives `A v^2 + 2 B(u) v + C(u) = 0`, so the
  curve is `v = (-B +- sqrt(D)) / A` where `D(u) = B^2 - A C >= 0`, two
  branches joined where `D` vanishes. This is D13's certified
  parameterisation: a point is evaluated in rational intervals with an
  enclosure; the pcurve on the ruled surface is `(u, v)` exactly and on the
  other surface its projection.
* **Topology by exact predicates.** For two cylinders with crossing axes
  the classes are decided exactly from the squared axis distance `d^2`
  against `(r1 + r2)^2` and `(r1 - r2)^2`: empty, a tangent point, one loop,
  a figure-eight (internal tangency, a node) or two rings around the thinner
  cylinder (equal radii at `d = 0` are S7a's ellipses); a cylinder and a
  sphere likewise from the distance of the centre to the axis against the
  radii. The roots of `D` (a loop's ends) are then simple: each is certified
  by a sign change of `D` and a derivative bounded away from zero, and their
  number must be the class's, or the result is `ComputationLimit`.
  *Amended in implementation:* in the canonical frame `D` has a closed form
  (Lagrange's identity for two cylinders, the centre's distance for a
  sphere), so a loop's ends are `+-arccos c` by one certified arctangent. The
  subdivision first implemented took over 600 s on a nearly tangent pair the
  fuzz target found (`fuzz/regressions/README.md`); S7b.2's cones, which have
  no such form, keep subdivision with a mean-value enclosure and a budget.
* **Order.** S7b.1: cylinder/cylinder and cylinder/sphere. S7b.2: pairs with
  a cone (transcendental coefficients: the classes come from certified root
  counts, exact only where the apex is rational). S7b.3: tori, by their
  meridians (a plane or a sphere meets a meridian circle in a quadratic;
  cylinders, cones and tori in a quartic). Each sub-step keeps the evidence
  order: an independent 60-digit reference and a native `GeomAPI_IntSS`
  capture (sample points on both surfaces, curve counts) before its kernel
  code.

Decisions for S7b.2, recorded before its code (2026-09-27):

* **Sphere and cone** (the centre off the axis): parameterised on the cone's
  rulings through its apex `V`, `V + v d(u)` with `d(u) = cos a a + sin a
  (cos u x + sin u y)` unit and `x` towards the centre, so `A = 1`,
  `C = |V - c|^2 - R^2` constant and `B(u) = b0 + b1 cos u` a sinusoid. The
  apex inside the sphere (`C < 0`): two rings, one on each nappe. Outside
  (`C > 0`): loops where `|B| > sqrt(C)`, one around `u = 0` and one around
  `u = pi`, each end `arccos` in closed form, or two rings when `B` keeps a
  sign beyond `sqrt(C)`. The apex on the sphere (`C = 0`, possible only for a
  rational apex) makes one root the apex itself and the other switch
  branches where `B` vanishes; it stays `NotConic` until a component for it
  exists.
* **Cylinder and cone** (axes not coaxial): parameterised on the cylinder,
  where `A = cos^2 a - (a1 . a2)^2` is a nonzero constant (zero would make the
  cylinder's axis a generatrix direction, impossible for a binary64 angle).
  `D(u)` is a trigonometric polynomial of degree two without a convenient
  closed form: its roots are isolated by certified subdivision with the
  mean-value enclosure and a work budget (S7b.1's lesson), and the
  components come from the certified count and signs: loops between roots
  where `D > 0`, two rings where `D > 0` throughout. A tangency (a double
  root) cannot be certified and is `ComputationLimit`; it never occurs
  exactly with a transcendental half-angle.
* **Two cones** (not coaxial): parameterised on the first, `A(u)` vanishes
  where the rulings are parallel and a branch goes to infinity; they stay
  `NotConic` until unbounded components exist (a later part of S7b).
* **Evidence first**, as for S7b.1: the reference (same parameterisation and
  frames, 80-digit roots by dense sampling) and a `GeomInt_IntSS` capture of
  these pairs before their kernel code.

Decisions for S7b.3, recorded before its code (2026-09-27):

* **Split.** S7b.3a: a torus with a plane or a sphere (every position), and
  every coaxial torus pair (a plane normal to the axis, a sphere centred on
  it, a coaxial cylinder, cone or torus). S7b.3b: a torus with a cylinder, a
  cone or a torus off its axis, where a meridian meets the other surface in
  a quartic. Unbounded components (two cones, the apex on a sphere) stay a
  later part of S7b.
* **Meridians.** A torus (axis `a` normalised exactly, major `R`, minor
  `r`) is `C(phi) + r (cos t e(phi) + sin t a)`, `C = o + R e`,
  `e(phi) = cos phi x + sin phi y` in an exactly orthonormal frame, `x`
  along the component normal to the axis of the plane's stored normal or of
  the direction from the torus's origin to the sphere's centre. A plane or a
  sphere restricted to a meridian circle is `f0 + alpha cos t + beta sin t`,
  so each `phi` gives the points of two branches, `t = atan2(beta, alpha) +-
  arccos(-f0 / sqrt(alpha^2 + beta^2))`, evaluated algebraically (the unit
  vector `(cos t, sin t)` is `(-alpha f0 -+ beta sqrt(D), -beta f0 +- alpha
  sqrt(D)) / (alpha^2 + beta^2)`), joined where `D = alpha^2 + beta^2 - f0^2`
  vanishes. `D` is a quadratic `P(c)` in `c = m cos phi` with rational
  coefficients and a negative leading one (`r < R`), `m` the length of that
  normal component (the square root of a rational).
* **Classes by exact predicates.** The signs of `P(m)`, `P(-m)` (numbers
  `u + v sqrt(q)`, decided exactly), of `P`'s discriminant and the vertex's
  position against `+-m` give every class: empty, one or two tangent points,
  one loop around `phi = 0` or `pi`, two loops (mirror images), two loops
  touching at `phi = 0` or `pi` (a node: two components sharing an end), two
  rings, a figure-eight with its node at `0` or `pi`, and Villarceau's two
  circles (`P(m) = P(-m) = 0`: two loops `[0, pi]`, `[pi, 2 pi]` sharing both
  ends). Loop ends are `arccos` of a root over `m` in closed form. A
  component's arcs between its ends and nodes are what S8 and S9 need; which
  arc continues which smoothly through a node is not recorded.
* **Exact special cases.** A plane containing the axis: two meridian
  circles. A plane normal to the axis or a sphere centred on it: circles
  about the axis (two, one tangent, none). A sphere containing a meridian
  circle (the centre in the equatorial plane on the tangent line of a
  meridian circle's centre, `rho^2 = r^2 + |w|^2 - R^2`) stays `NotConic`.
  Coaxial pairs meet in circles about the axis where their meridians meet in
  the half-plane: a cylinder (a vertical line), a cone (two lines through the
  apex, certified; never exactly tangent), a torus (two circles, exact
  classes; the same torus is `Same`).
* **Evidence first**, as before: the reference extended with the meridian
  parameterisation (80-digit `D` by evaluating the other surface on the
  meridian circle, roots by dense sampling, classes by its own exact
  predicates) and a `GeomInt_IntSS` capture of the S7b.3a cases before any
  kernel torus intersection code.

Decisions for S7b.3b, recorded before its code (2026-09-27):

* **Meridians again.** A torus and a cylinder, a cone or another torus off
  its axis are parameterised on the torus's meridians (of two tori the first
  by stored data), `(phi, t)` on the flat parameter torus, in the frame whose
  `x` is the other axis's component normal to the torus's axis (the
  direction to the other's origin when the axes are parallel). The curve is
  the zero set of `G(phi, t) = f(p(phi, t))`, `f` the other surface's
  implicit function (a quadric: `G` of degree two in `cos t`, `sin t`, up to
  four points per meridian; a torus: four, up to eight). The
  parameterisation is a diffeomorphism, so the curve's topology is that of
  `G = 0` on the flat torus; no closed form separates the branches.
* **The curve as a graph.** Vertices are the tangencies (singular points);
  edges are arcs between them, or closed smooth components without one.
  Each arc is a chain of tracks (graphs `t_i(phi)` over a `phi`-range)
  joined at folds, where a component turns in `phi` (`G = G_t = 0`). A
  closed smooth component keeps its winding numbers on the torus. S8 and S9
  consume the arcs as edges.
* **Certified construction.** Folds are found by subdivision of the
  parameter torus (binary64 intervals, a box they cannot settle again in
  rational intervals; the mean-value form excludes a box where `G` or `G_t`
  certainly does not vanish) and certified by the Krawczyk operator of
  `(G, G_t)` (a unique regular solution); a work budget bounds it (S7b.1's
  lesson). Each fold, and each tangency, gets a box on whose top and bottom
  edges `G` does not vanish and whose side edges carry the certified simple
  roots of its local picture (two on one side and none on the other for a
  fold, two and two for a crossing, none for an isolated point), with a
  unique critical point of `G` inside for a tangency. Between them the
  branches are followed by chains of certified windows (a sign change of `G`
  across the window for every `phi` in its step, `G_t` of one sign on it),
  which fixes how branches continue (through `t = pi` too) and gives every
  point evaluation its window. A box or window that cannot be certified
  within the budget is `ComputationLimit`.
* **Tangency exactly.** A singular point of the curve is a tangency of the
  surfaces. A torus is the pipe of radius `r` about its spine circle and a
  cylinder the pipe of radius `r_c` about its axis, so they touch exactly
  where the distance between the spine and the axis has a critical value
  `r + r_c` or `|r - r_c|`: three conics in a rational basis of the spine's
  plane, decided by exact resultants and real root isolation
  (`polynomial::real`), the basis sheared until distinct critical points
  have distinct coordinates. The contact is an isolated point or two
  crossing branches by the certified sign of the Hessian of `G` there; a
  higher contact, and an axis meeting the spine with `r = r_c`, are
  `ComputationLimit`. A cone is never exactly tangent to a torus for a
  binary64 half-angle (its `cos` and `sin` are transcendental). Two tori
  touch where their spines have a critical distance `r1 +- r2` (S7b.3b.2).
* **Order.** S7b.3b.1: a torus and a cylinder or a cone. S7b.3b.2: two tori.
  Each keeps the evidence order: an independent reference (the same
  meridians; every meridian's roots from 80-digit polynomial roots in
  `z = e^{it}`; the critical meridians as the real roots of the resultant of
  `z^2 G` and `z^2 G_t`, a trigonometric polynomial found from exact samples;
  components by continuity between meridians placed about them; tangencies
  by a Groebner basis in sympy) and a `GeomInt_IntSS` capture before the
  kernel code.

Decisions for S7b.3b.2, recorded before its code (2026-09-27):

* **The same graph.** Of two tori the carrier is the first by stored data;
  the other's implicit quartic `(|w|^2 + R2^2 - r2^2)^2 - 4 R2^2 |M2 w|^2`
  restricted to a meridian circle is still of degree two in `cos t`,
  `sin t` (`|w|^2` is affine along a circle), so folds, boxes and tracks are
  S7b.3b.1's with the quartic's gradient and Hessian in the jets.
* **Tangency exactly.** Two tori touch where their spines have a critical
  pair at the distance `r1 + r2` or `|r1 - r2|`. For a point `s1` of the
  first spine (`u`, `v` in a rational basis of its plane, the circle `E0`)
  the distance to the second spine is critical along it at
  `D = |w|^2 + R2^2 + 2 sigma R2 rho` (`w = s1 - o2`, `rho = |M2 w|`,
  `sigma = +-1`); `D = k` gives `2 sigma R2 rho = Q`, `Q = k - |w|^2 - R2^2`,
  and the criticality along the first spine becomes the cubic
  `E1 = (w . T1) Q + 2 R2^2 (M2 w . T1)` and the distance the quartic
  `E2 = Q^2 - 4 R2^2 |M2 w|^2`, `sigma` the sign of `Q`. Both reduce modulo
  `E0` (its `v^2` coefficient is a nonzero constant) to `c1 v + c0` and
  `d1 v + d0`; the critical points are the real roots of
  `a2 c0^2 - a1 c0 c1 + a0 c1^2` with `v = -c0 / c1`, and a tangency where
  `d0 c1 - d1 c0` vanishes on the algebraic root. `Q = 0` there (the point
  on the other's axis) and spines meeting with `r1 = r2` are
  `ComputationLimit`, as are tori touching along a curve.
* **Evidence first**: the reference with the torus's quartic and its
  tangencies by a Groebner basis of both spines' points (four coordinates,
  five equations), and a `GeomInt_IntSS` capture of 14 pairs before the
  kernel code.

Decisions for S7b.4, recorded before its code (2026-09-28):

* **Scope.** Two cones not coaxial, and a cone whose rational apex lies
  exactly on a sphere or a cylinder: the pairs S7b.2 left `NotConic` or could
  not certify. Two cones with one rational apex meet in their common
  generatrices (the certified roots of `A(u)`, S7a's lines), or the apex
  alone. Their curves may be unbounded
  (a branch escapes where a ruling of the first cone is parallel to a
  generatrix of the second) or pass through the apex.
* **The projective ruling.** On the first cone's rulings through its apex
  (of two cones the first by stored data), `V + v d(u)`, the other quadric is
  `A(u) v^2 + 2 B(u) v + C`; with `v = tan(t / 2)`,
  `G(u, t) = ((A + C) + (C - A) cos t) / 2 + B sin t` is smooth on the torus
  `(u, t)` mod `2 pi` and its zero set is the curve with its points at
  infinity (`t = pi`). S7b.3b's traced graph applies unchanged (its field
  becomes a chart: the torus's meridians or the cone's rulings); a component
  crossing `t = pi` is unbounded, and its crossings of infinity are counted
  (it may cross and return, winding zero).
* **Factors.** The apex exactly on the other surface (`C = 0`): `G = sin(t/2)
  (A sin(t/2) + 2 B cos(t/2))`, the apex for every ruling and the curve
  `A sin psi + 2 B cos psi = 0` (`psi = t / 2`, period `pi`), one point per
  ruling, through the apex where `B` vanishes (a crossing there, or the apex
  isolated when `B` keeps a sign; a double root of `B` is
  `ComputationLimit`). Parallel cones of equal half-angles (`A = 0`
  identically): `G = cos(t/2) (C cos(t/2) + 2 B sin(t/2))`, the circle at
  infinity common and the conic `2 B sin psi + C cos psi = 0` in their radical
  plane. Both are traced on the `(u, psi)` torus with the same machinery.
* **Tangency.** A cone's `cos` and `sin` are transcendental for every
  nonzero binary64 half-angle, so a tangency away from the apex needs a
  symmetric configuration. *Amended in implementation:* such configurations
  occur exactly: two congruent cones whose axes cross at a point equidistant
  from their apexes (S7a's `kk_crossing`) meet in two conics crossing where
  the cones touch, whatever the half-angle. Those nodes are not decided
  exactly: the certification cannot settle them and the result is
  `ComputationLimit`, which S7a's comparison accepts for that `not_conic`
  row.
* **Evidence first**: the reference (the same projective parameterisation
  with `A`, `B`, `C` found by evaluating the other surface on the ruling at
  `v = -1, 0, 1`, the scan, critical meridians and components of the torus
  reference, the factors in closed form) and a `GeomInt_IntSS` capture of 13
  pairs before the kernel code.

Decisions for S7c, recorded before its code (2026-09-28):

* **Scope and order.** S7c.1: the topology's analytic edges (lines, circles
  and arcs of a `Frame3`) against every analytic surface (plane, cylinder,
  cone, sphere, torus): what S8's split needs beside the existing
  spline/plane intersections. S7c.2: S7a's ellipses and hyperbolas (the
  edges S8 and S9 create) against every analytic surface, and spline edges
  against cones and tori. Procedural and traced curves against surfaces
  come with S9, where they become edges.
* **Exact curves.** A line is the exact set through its two rational points,
  parameterised `p0 + s (p1 - p0)`; a circle the exact point set of its
  stored frame (the plane through the origin normal to the stored normal,
  the stored radius), parameterised by the angle from the stored x axis as
  `Curve3::Circle`. A segment or an arc keeps the points within its range;
  S7c.1 intersects the whole curve.
* **Methods.** A line against a plane, cylinder, sphere or torus: the
  surface's function along it is a polynomial of degree one, two or four
  with rational coefficients, its real roots isolated exactly with their
  multiplicities (`polynomial::real`). A circle against them: in a rational
  basis `(U, V)` of its plane the surface's function reduced modulo the
  circle's conic `E0` is `c0(u) + c1(u) v`; the points are the real roots of
  the resultant `a2 c0^2 - a1 c0 c1 + a0 c1^2` with `v = -c0 / c1` (the basis
  sheared as in `tangency.rs`), their multiplicities the contacts, exactly;
  the angle is a certified `atan2` in the stored frame. A cone's `cos` and
  `sin` are transcendental: a line's quadratic and a circle's trigonometric
  polynomial of degree two in intervals, a root certified by a sign change
  and a derivative of one sign (a tangency there is `ComputationLimit`).
* **Degeneracies, exactly.** The function's being zero identically along the
  curve (`c0 = c1 = 0`, or the zero polynomial) is `Contained`: a line in a
  plane or along a cylinder, a circle in a plane, on a sphere, on a coaxial
  cylinder, a torus's meridian or parallel. A double root is a tangency.
* **Results.** `curve_surface(curve, surface)` gives `Empty`, `Contained`, or
  points, each with its curve parameter and point enclosed and its contact
  (crossing or tangent).
* **Evidence first**: an independent reference (sympy: the exact polynomial
  along a line, `real_roots` with multiplicities; the circle's two
  equations by a Groebner basis, a tangency where their gradients are
  parallel; 80-digit roots for cones) and a native `GeomAPI_IntCS` capture
  of 36 cases before the kernel code.

Decisions for S7c.2, recorded before its code (2026-09-28):

* **Conic curves.** S7a's ellipses and hyperbolas become edges in S8 and
  S9 with stored binary64 data, as OCCT's `Geom_Ellipse` and
  `Geom_Hyperbola`: `intersection::Conic` is an ellipse or a hyperbola of a
  `Frame3` with semi-axes `major` and `minor`, the exact point sets
  `o + major cos t x + minor sin t y` and `o + major cosh t x + minor sinh t
  y` (the stored axes `x`, `y`; one branch of the hyperbola), `t` the
  parameter of `ElCLib` (an ellipse's in `(-pi, pi]`, a hyperbola's real).
  S8 gives them their `Curve3` variants with the same data.
  `conic_surface(conic, surface)` gives `Empty`, `Contained` or points, as
  `curve_surface`.
* **Conic methods.** In the stored axes, sheared, a point `o + xi x + eta y`
  is on the curve exactly when `minor^2 xi^2 +- major^2 eta^2 = major^2
  minor^2` (and `xi > 0` on the hyperbola's branch, a certified sign): the
  circle's resultant method of S7c.1 applies unchanged, with exact
  multiplicities and containment against planes, cylinders, spheres and
  tori. The parameter is a certified `atan2(eta / minor, xi / major)` or
  `asinh(eta / minor)` (a certified logarithm joins `certified.rs`). Against
  a cone an ellipse is the circle's trigonometric polynomial of degree two;
  a hyperbola's, times `4 z^2` with `z = e^t`, a quartic in `z > 0` with
  interval coefficients, its roots isolated by certified subdivision
  between Cauchy bounds of it and its reciprocal, `t = ln z`; a tangency
  there is `ComputationLimit`.
* **Splines against tori.** `spline_torus`: the torus's quartic in
  homogeneous form, `(|D|^2 + (R^2 - r^2) W^2)^2 - 4 R^2 W^2 (|D|^2 -
  (D . a)^2 / |a|^2)` with `D = X - o W`, of degree `4 p` on each span,
  through the shared exact isolator of `spline_sphere` and `spline_cylinder`:
  every root with its contact orders, and every span on the torus as an
  overlap.
* **Splines against cones.** A cone is the set `rho = |r + h tan a|` about
  its stored axis, so on a span `W^2 F = Q0 + tau Q1 + tau^2 Q2` with exact
  polynomials `Q0 = |D|^2 - (D . a)^2 / |a|^2 - r^2 W^2`, `Q1 = -2 r (D . a)
  W`, `Q2 = -(D . a)^2` and one irrational `tau = tan a / |a|`. `tau` is
  transcendental (`a` a nonzero rational, Lindemann-Weierstrass), so a span
  lies on the cone exactly when `Q0 = Q1 = Q2 = 0`, and a root common to all
  three (their exact gcd: a rational apex) is an exact point, a tangency
  (the cone's apex is singular). The other roots are isolated by certified
  subdivision with `tau` in intervals, each by a sign change and a
  derivative of one sign; one that cannot be certified is
  `ComputationLimit`. `spline_cone` gives points (parameter, point,
  contact) and overlaps (whole spans, exact knots).
* **Evidence first**: the reference extended (sympy: the conics' two
  equations by Groebner bases, a spline's exact span polynomials against a
  torus with `real_roots` and multiplicities; 80 digits for cones), fixtures
  of ellipses, hyperbolas and splines against every surface with each class,
  and a native `GeomAPI_IntCS` capture on `Geom_Ellipse`, `Geom_Hyperbola`
  and `Geom_BSplineCurve` before the kernel code.

Decisions for S7d, recorded before its code (2026-09-28):

* **Scope and order.** S7d.1: pairs of the analytic curves (a line segment's
  whole line, a circle or an arc's whole circle, S7c.2's ellipses and
  hyperbolas' branches): `curve_curve(a, b)` over
  `intersection::AnalyticCurve` (`Edge(Curve3)` or `Conic(Conic)`) gives
  `Empty`, `Coincident` (the same point sets), or points, each with both
  curves' parameters and the point enclosed and its contact (tangent where
  the curves share their tangent line there, else crossing), sorted by the
  first curve's parameter. S7d.2: spline edges against these curves and
  each other, decided with its own decisions after S7d.1.
* **Planes, exactly.** A circle lies in the plane through its origin normal
  to its stored normal; an ellipse or a hyperbola in the plane through its
  origin spanned by its stored axes (normal `x × y`, rational). A line is
  `p0 + s (p1 - p0)`.
* **Methods.** Two lines: exact linear algebra (coincident, parallel, skew,
  or one rational point). A line and a plane conic: a line crossing the
  plane meets it in one rational point, on the conic exactly when the
  conic's equation vanishes there (a crossing); a line in the plane
  substitutes into the conic's equation, a quadratic with rational
  coefficients (a double root a tangency); a line parallel to the plane
  off it misses. Two plane conics in one plane: S7c's resultant, the second
  conic's equation written in the first's plane by exact frame
  coordinates (a multiple root a tangency, `c0 = c1 = 0` coincidence, a
  hyperbola's branch by a certified sign). In two crossing planes: both
  conics' equations along the planes' common line are quadratics with
  rational coefficients; the points are the real roots of their gcd, a
  tangency where the root is double in both. Parallel planes: empty.
* **Evidence first**: an independent reference (sympy: each curve's
  rational parameterisation substituted into the other curve's plane and
  conic equations, the common real roots of the two polynomials by their
  gcd; a tangency where the tangent directions are parallel), fixtures of
  every pair and class (crossings, tangencies coplanar and across planes,
  coincidence, parallel and skew lines, the other branch, misses), and a
  native `IntTools_EdgeEdge` capture on edges of those curves before the
  kernel code.

Decisions for S7d.2, recorded before its code (2026-09-28):

* **Scope.** Rational B-spline edges against circles, ellipses and
  hyperbolas' branches (`spline_curve(spline, curve)`; against lines
  `spline_line` already serves). Two spline edges meet in 3D only at
  isolated, non-generic points whose exact decision needs resultants of
  degree up to `2 p q` in the span parameters (1,250 at degree 25): they
  come with S9, where spline edges meet through the faces that share them,
  as S7c left procedural and traced curves to S9.
* **Method.** On each span, with the homogeneous polynomials `(X, W)`, the
  conic's plane gives `L = n . (X - o W)` (degree `p`) and its equation
  `Q = W^2 E(X / W)` (degree `2 p`), both exact. A span lies on the conic
  exactly when `L = Q = 0` (an overlap between exact knots); otherwise the
  points are the real roots of `gcd(L, Q)` in the span, with the exact
  isolator, a hyperbola's branch by a certified sign. The spline's tangent
  is the conic's there exactly when the root is at least double in both
  (in the gcd): `L'` vanishing puts the tangent in the plane, `Q'` along the
  conic. Knots shared by two spans give one point.
* **Results.** Points with the spline's parameter, the conic's, the point
  and the contact, and overlaps as spline parameter ranges.
* **Evidence first**: the reference extended (exact span polynomials into
  the conic's implicit equations in `(X, Y, Z)`, sympy's gcd and
  `real_roots`, tangency by parallel tangent directions), fixtures of each
  class (crossing through the plane, in-plane crossings and tangencies,
  a span on the circle as an exact rational quarter circle, a partial
  overlap, the hyperbola's other branch, misses), and a native
  `IntTools_EdgeEdge` capture before the kernel code.

### S8 — general planar split and face trimming (SplitBody job)

Split any supported solid by an arbitrary plane: face/plane intersection
curves, loop splitting on the universal cover, region classification, the
`Split`, `Generated` and `Deleted` relations M3 introduced but on real
geometry. Native `BRepAlgoAPI_Splitter` bridge; upstream `bsplit` cases.
The first general topology-changing algorithm, and the rehearsal for S9.

Decisions for S8, recorded before its code (2026-09-28):

* **Order.** S8a: prisms of line and arc profiles (polygons, circles,
  stadiums, holes): planar faces and cylindrical walls, so the new edges are
  segments and the walls' sections (lines, circles or ellipse arcs). S8b:
  spline prisms (a plane crossing a linear extrusion of a spline meets it in
  an affine image of the profile, an exact spline). S8c: the cylinder, cone
  and sphere primitives (conic arcs, a cone's lines through its apex). S8d:
  tori (D13's procedural edges). S8e: the upstream `bsplit` group and the
  sheets and wires. Each sub-step keeps the evidence order.
* **Operation.** `Solid::split_by_plane(operation, plane)` (`plane` a
  `Frame3`: origin and normal) returns the pieces as solids, those against
  the normal first, then along it, each in a deterministic order (by the
  lowest id of its input faces), with the operation's history. A plane
  missing the solid returns it unchanged. A piece that is itself a prism of
  the input's frame (a plane normal to the axis: M3's height split; a plane
  parallel to it: the prism of the profile's section by the plane's line,
  its ids renamed from the input by provenance) is built as one, keeping the
  prism's exact queries (S8a.1); any other piece is a general body built
  through `TopologyParts` and validated before it is returned (S8a.2 on:
  planes oblique to a prism's axis, and the other families).
* **Exact decisions, binary64 geometry.** Which side of the plane every
  vertex lies on, whether an edge or a face crosses, touches or lies in the
  plane, and where boundaries meet the section are exact predicates on the
  stored data (S7's intersections). New vertices, edges and pcurves store
  binary64 data rounded from those exact or enclosed values; as everywhere
  in the topology (TOPOLOGY_MODEL.md) each fin carries its certified
  enclosure and the body's validation certifies every new entity within the
  resolution. D13's exact representation applies to procedural edges
  (S8d), whose pcurves are their exact projections.
* **New edge and pcurve kinds.** `Curve3` gains ellipse and hyperbola arcs
  (S7c.2's `Conic` data with a start angle or parameter and a sweep) and,
  in S8d, procedural edges. `Curve2` gains, on a planar face, ellipse and
  hyperbola arcs in the plane's frame, and on a cylinder the plane section's
  graph `v = a0 + a1 cos u + a2 sin u` over `u = u0 + sweep f` (on a cone
  the quotient of two such, S8c). Validation, `.brep` I/O (OCCT's
  `Geom_Ellipse`, `Geom_Hyperbola`, and B-spline pcurves where OCCT has no
  analytic kind), tessellation and mass properties take each kind as it
  lands.
* **Degeneracies, exactly.** A face lying in the plane is not cut: it and
  its solid stay on the side the rest of the solid is on (a cap in the
  plane: the solid is unchanged). A plane tangent to a wall along a ruling,
  or through a vertex or along an edge without crossing the solid there,
  cuts nothing at that contact. Pieces meet the plane only in their cut
  faces: one planar face per region of the solid's section, with its holes.
* **History.** A face, edge or region crossing the plane is `Split` into one
  child per piece it reaches (ordinal by piece); entities on one side are
  `Unchanged` in their piece; each cut face is `Generated` from every face
  its boundary runs along, each cut edge from the face it cuts and each cut
  vertex from the edge it cuts, as M3's height split does. The height split
  of M3 stays, and on its inputs the general split gives the same pieces.
* **Evidence first**: an independent reference (mpmath: each piece's volume
  and centre by quadrature of the clipped height over the profile, its area
  as its caps', walls' and cut faces' by the same slicing; the pieces'
  face, edge and vertex counts come from the native capture and the
  kernel's validation, not from a second arrangement), fixtures of every class
  (planes crossing caps and walls at angles, parallel and normal to the
  axis, through a vertex, containing an edge, tangent to an arc wall,
  missing, in a cap), and a native `BRepAlgoAPI_Splitter` capture (each
  piece's volume, area, centre and counts) before the kernel code.

Decisions for S8c, recorded before its code (2026-09-28):

* **Order changed: S8c before S8b.** S8b splits spline prisms, and no
  builder makes one yet: a profile holds lines and arcs (S5). It needs
  spline profile segments first (their validation as simple closed
  curves by certified spline proximity, exact area moments, the extruded
  wall as an exact degree-(p, 1) spline surface), whose tessellation and
  precise mass are the running parallel tracks T-b and F8. S8b follows
  them; S8c needs nothing new but its own evidence.
* **Scope.** The cylinder primitive is a prism of a circle, split by S8a.
  S8c.1: a cone, frustum or sphere zone by a plane normal to its axis, the
  pieces the same primitive between exact heights rounded to binary64
  (M3's height split for revolved solids: the cut face a disc or, through
  the apex, the apex itself); a whole sphere by any plane, its pieces two
  caps whose frame's axis is the plane's normal (a sphere has every axis).
  S8c.2: a plane containing the axis, the pieces half-solids bounded by
  meridians (a cone's rulings, a sphere's great half-circles, both line
  pcurves) built as general bodies. Other planes cut a cone in conics and
  a zone in circles whose pcurves are transcendental graphs over the
  angle (OCCT approximates them by splines): they wait for D13's
  procedural edges (S8d) and are `OutOfDomain` until then.
* **History.** As S8a: pieces renamed by provenance, the lateral face
  `Split` into faces of the same surface. The history checker's support
  test learns cones, spheres and tori (a piece lies on the whole's surface
  when the axes are parallel and the apex, centre or centre circle within
  the tolerance, whatever their stored frames' x axes).
* **Evidence first**: `split_reference.py` extended to cones and spheres
  (each side's volume, area and centre by quadrature of the clipped
  revolved profile), fixtures of every class, and `occt_split_oracle.cpp`
  extended with `BRepPrimAPI_MakeCone`/`MakeSphere` and a capture before
  the kernel code.

Decisions for S8e, recorded before its code (2026-09-28):

* **Scope.** S6's bodies by a plane: `Body::split_by_plane(operation,
  plane)` returns the pieces with their sides (below first) and the
  history, as `Solid::split_by_plane` does. A body is a planar sheet
  (`Body::face_from_profile`) or a closed wire (`Body::wire_from_boundary`)
  on its frame's plane, of lines, arcs, circles and (S8b) splines; the
  plane's trace on the body's plane is the line `a u + b v + d = 0` in the
  frame's coordinates, its coefficients exact rationals of the stored data.
* **Parallel planes.** A plane parallel to the body's (`a = b = 0`) cuts
  nothing: the body is returned itself, every entity `Unchanged`, on the
  side it lies on, or `Below` when it lies in the plane (a face in the
  plane stays whole, S8's rule).
* **Sheets.** The face's pieces are the profile's section by the trace
  (S8a.1's and S8b.3's `Section`, unchanged: exact sides, crossings,
  tangencies and chords); each piece is a planar sheet of its profile on
  the body's frame, renamed by provenance as a prism's caps are: a whole
  edge or vertex keeps its id, parts are `Split` children with canonical
  ordinals, a vertex on the plane is split into one copy per piece, each
  cut edge (a chord) is `Generated` from the face and each cut vertex from
  the edge it cuts; the face is `Split` into one child per piece. A trace
  missing or touching the face returns the body.
* **Wires.** A closed wire's crossings with the trace (the same exact
  section of its boundary, without chords) cut it into runs; each maximal
  run of consecutive pieces on one side is one open wire (a new
  construction: a path of segments with two free ends, its edges and
  vertices on the frame's plane, validated as a wire body), a run's order
  the boundary's stored order from the run's first piece. An edge lying
  along the trace joins the run it continues in stored order; a wire
  touching the trace without crossing it returns the body. History: whole
  edges keep their ids, parts are `Split`, a crossing's vertex is
  `Generated` from the edge it cuts and a stored vertex on the plane
  between two runs is split into one copy per run; the wire's region (the
  void) is unchanged.
* **Upstream.** The `bsplit` group (`boolean/splitter/A5`, `B5` and the
  cases of `bugs/*` that call `bsplit`) is registered in
  `upstream-draw.json` as capability sentinels (their shapes need the
  general builder of S9; native DRAW must pass them where the dataset is
  present). A derived case `split_plane` drives the adapter's
  `bclearobjects`, `baddobjects`, `baddtools`, `bfillds` and `bsplit` with a
  plane face (`plane` then `mkface`) as the tool on a prism, a sheet and a
  wire, checking `checkshape`, `checknbshapes` and `checkprops` of the
  result on both backends.
* **Evidence first**: the reference extended to sheets (each side's area
  and centroid from the profile's section) and wires (each side's length
  and centroid, exactly for lines and arcs, by quadrature for splines),
  fixtures of every class (crossing, through a vertex, along an edge,
  tangent to an arc or a hole, missing, parallel, containing the body, in
  both frames; lines, arcs, circles, holes and splines), and a native
  `BRepAlgoAPI_Splitter` capture of sheets and wires split by a plane face
  (each piece's area or length, centre and counts) before the kernel code.

Decisions for S8b, recorded before its code (2026-09-28):

* **Unblocked.** T-b tessellates spline edges, pcurves and faces and F8
  encloses spline mass properties within 1e-12 of their scale; S8b now
  follows S8d, before S8e.
* **S8b.1, spline profile segments.** `Segment::Spline` holds a planar
  nonrational B-spline (`BSplineCurve2`, degree 1 to 7) whose first and last
  poles are its path points exactly and which is C1 inside (R4's exact test);
  a rational spline is `OutOfDomain` (its moments are not polynomial).
  `Segment` stops being `Copy`. Validation keeps its rules: a boundary is a
  simple closed curve, each pair of pieces farther apart than the resolution
  except adjacent ones at their shared point, decided by a certified
  separation screen: exact Bernstein subdivision of each spline piece with
  control-hull boxes against lines, arcs and spline pieces, in the binary64
  tier then the rational one, an undecided pair treated as touching
  (`SelfIntersection`, as `decide::area_is_degenerate` does); adjacent
  pieces leave their shared point in directions decided exactly and are
  separated beyond a ball about it. Area moments are exact Bernstein
  integrals; point location uses the crossing parity and the certified
  point-to-spline distance.
* **S8b.1 amendments, from its implementation (2026-09-28).** (a) A
  spline segment is also screened against itself: its Bézier arcs are cut
  (exactly, at halves) into pieces whose control polygons turn through
  less than a half-turn, each therefore monotone along a direction and
  simple, and the pieces are screened pairwise, consecutive ones exempt
  only in disjoint direction sectors from their shared point. A cusp's
  halves double back there (`SelfIntersection`); a turn that does not
  narrow within twelve halvings is `OutOfDomain`. (b) The screen
  subdivides in outward binary64 intervals: each control point's box
  holds the exact one, so a part lies in their hull, and a part keeps its
  exact ends where they are its original's, so a sector from a shared
  point is the cone of the other boxes' corners. (Exact rational halvings
  at depth 48 made a lens hole's oblique split take 66 s under the
  sanitizer in S8b's first campaign; with them, a sampled orientation area
  and an oblique piece's mass recomputed under rigid motion it took 7.4 s
  in release, now 2.2 s.) (c) Point location takes the crossing parity of the `+x`
  ray with exact `y`-monotone pieces, and a point within the resolution of
  a spline, or undecided, is on the boundary, as for arcs.
* **S8b.2, spline prisms.** A spline segment's wall is the exact degree-(p,
  1) B-spline surface over its knots and the heights (its poles the
  profile's lifted by the frame, rounded); its cap edges are `Curve3`
  splines, its cap pcurves the profile spline in the cap's frame, its wall
  pcurves lines in the surface's parameters. Mass (F8), tessellation (T-b),
  validation and `.brep` records (S4) apply unchanged; roles and ids are
  the arcs' (identity encodes elements, not geometry).
* **S8b.3, splits.** Normal to the axis: the height split. Parallel: the
  profile's section by the plane's line, crossings by exact spline/line
  preimages, each piece's spline the exact restriction at the crossing's
  rounded parameter with its poles rounded (within the resolution of the
  original). Oblique: S8a.2's footprint and crease, the crease an affine
  image of the footprint's spline (exact, rounded) with a spline pcurve on
  the wall (identity in `u`, affine in `v`).
* **S8b.3 amendments, from its implementation (2026-09-28).** (a) Roots
  come with exact multiplicities: a simple root crosses, an even one
  touches (the side kept; a touch inside the solid pinches as an arc's
  does), an odd one of three or more or a tangency at a knot is
  `Degenerate`, and a spline along the trace `OutOfDomain`; sides are exact
  signs between distinct roots. (b) A crease part's ends take their
  vertices' heights exactly (the edge's end poles and the pcurve's end `v`),
  so no pcurve leaves the wall's `v` range by rounding. (c) The identity in
  `u` rounds (Greville abscissae are averages of knots), so a pcurve's
  hull crosses the wall's knot lines by a rounding step: the exact
  deviation bound (S4b) now takes the patch holding a piece but for slivers
  at most `2^-20` of its size, adding a certified bound of each
  neighbouring patch's departure from its polynomial over the sliver's box
  (both nets re-expressed exactly there; the difference lies in the hull of
  theirs); an undecided exact bound falls back to the Taylor enclosure. (d)
  Green's exact path (S4d) takes a patch whose box a pcurve piece's curve
  provably keeps to when its control polygon leaves it (a crease nearly
  touching a cap): each coordinate against each bound exactly nonnegative
  on `[0, 1]`. Without it the strips' enclosure was wide (a fuzz seed:
  0.13 relative). (e) The independent history check compares spline supports:
  a spline edge's exact Bézier arcs with its parent's over its range, a
  wall's rows with its parent's extrusion across its direction, a planar
  piece's spline boundary by its poles, each within the tolerance. (f) The
  `.brep` writer writes an edge whose span runs against its curve along it,
  used the other way (S8b.2's prisms of profiles given clockwise were
  unwritable).
* **Evidence first**: the independent reference with spline profiles
  (mpmath Bernstein areas and moments; prism volumes and split sides by
  slicing), fixtures of every class (splines crossing, touching and
  tangent to the plane's line, holes, mixed with lines and arcs), and a
  native capture (faces bounded by `Geom2d_BSplineCurve` edges, prisms,
  `BRepAlgoAPI_Splitter`) before the kernel code.

Decisions for S8d, recorded before its code (2026-09-28):

* **S8d.1, exact.** A torus (whole, band or wedge) by a plane normal to
  its axis cuts the tube in parallels, circles whose pcurves are lines of
  constant `v`; by a plane containing the axis in the tube's meridian
  circles, lines of constant `u`. The pieces are general bodies on the
  input's own torus surface with planar cut faces (annuli, discs); a whole
  torus by a plane containing its axis gives two half-turn wedges, the S3
  construction. (A torus v-segment is the revolved region between the
  tube's arc and the axis, as `BRepPrimAPI_MakeTorus` makes it, so a torus
  cut normal to its axis does not give one.)
* **S8d.2, D13's engine, with cones and zones by any plane.** A pcurve may
  be a `Projection`: the exact inverse of its face's surface map applied to
  the fin's own edge, lifted continuously on the cover from a stored start.
  Its deviation is zero by definition when it projects the fin's own edge
  onto the face's own stored surface; vertex gaps and degeneracy come from
  certified evaluation. Integrals along it (twice-areas, periodic areas,
  orientation fluxes, mass moments) come from certified adaptive
  quadrature: interval Taylor enclosures of the integrand on subintervals
  from enclosures of the curve's derivatives (the surface inverse
  differentiated implicitly), subdivided until the width asked for; ray
  crossings and boundary clearance from certified monotone pieces. A
  cone's and a sphere's plane sections are conics and stay explicit (D13):
  `EllipseArc`, new `HyperbolaArc` and `ParabolaArc`, and `CircularArc`,
  with `Projection` pcurves on the cone or sphere.
* **S8d.3, spiric sections.** A torus's plane section is `Curve3::Section`:
  a graph over the torus's angle between its turning points, `v = atan2(B,
  A) +- acos(C / |(A, B)|)` with `A, B, C` affine in `cos u, sin u`,
  evaluated with certified enclosures; its pcurves are `Projection`s. Tori
  by any plane.
* **S8d.3 refined, before its code (2026-09-28).** A plane cuts a whole
  torus in two loops winding once about the axis (the pieces tube bands,
  their walls wound in `u`), two loops winding once about the tube (the
  pieces C-shaped, their walls wound in `v`), or one contractible loop (a
  cap whose wall is a disc on the torus, and the rest, whose wall is the
  torus less that disc: unwound loops only, covering both periods), or it
  touches or misses the tube. The transitions (Villarceau circles, the
  lemniscate where a plane parallel to the axis touches the inner equator,
  tangency at a point) are `Degenerate`. `Curve3::Section` is the section
  as a graph over `u` (`v = atan2(B, A) +- acos(-C / |(A, B)|)`, `A, B, C`
  affine in `cos u`, `sin u`) or over `v` (`u` likewise), whichever keeps
  the slope within one; a loop is cut where its slope crosses one, a vertex
  there, so no edge reaches a turning point and each stays analytic: its
  jets, D13's quadrature and the tessellation's interval second
  derivatives apply unchanged. The validator learns the torus less discs
  (a face on a torus whose edge loops are all unwound with the inner
  sign), its mass the whole torus less the holes' Green integrals.
  Evidence first: `split-spiric-cases.txt` (bands, caps, C-shapes, a tilted
  frame) with the reference and a native capture before
  `solid/split/spiric.rs` exists.
* **Tessellation and interop.** Segment bounds come from second-derivative
  enclosures (closed forms for conics, interval evaluation otherwise). The
  `.brep` writer approximates a `Projection` or a `Section` by a B-spline
  whose certified bound becomes the edge's tolerance (D13's interchange
  approximation); until it does, those are `Unwritable`.
* **Evidence first**: `split_reference.py` extended to tori (each slice an
  annulus cut by the plane's line) and the conic and circle cases of S8c's
  fixtures, a native `BRepPrimAPI_MakeTorus` capture before the kernel
  code; the certified quadrature checked against mpmath integrals.

### S9 — Booleans for the analytic family (Combine job)

Fuse, cut and common: intersect faces (S7), split (S8), classify by
regions, assemble shells and regions, report complete histories with the
split and merge relations. Native `BOPAlgo` bridge; the `boolean` group has
1,291 fully covered data cases and thousands of self-contained ones, the
ledger's largest lever. Tangent and coincident faces are declared cases from
the first fixture, never deferred.

Decisions for S9, recorded before its code (2026-09-28):

* **Operations.** `Solid::fuse(operation, other)`, `Solid::cut(operation,
  tool)` and `Solid::common(operation, other)` return the result's solids
  (each a maximal connected solid region, ordered by the lowest input face
  id it keeps, none for an empty result) and the operation's history
  (`Fuse`, `Cut`, `Common`). A result that is itself a prism of one of the
  inputs' frames (one profile over one height range) is built as one,
  keeping the prism's exact queries; any other is a general body built
  through `TopologyParts` and validated before it is returned. A rigid
  motion of a general body moves its stored geometry (frames and points,
  rounded), measures its enclosures again and keeps every id.
* **Order.** S9a: two prisms whose frames have bitwise-equal axes and whose
  origins differ by a vector with binary64 coordinates in them (the
  coordinates decided exactly), of lines, arcs and circles: the Boolean is
  a stack of height slabs (between the four caps' heights), each slab's
  profile the 2D Boolean of the profiles present in it, consecutive slabs
  with equal profiles merged (M3's fuse). S9a.1: results that are one
  prism (every common; a cut whose tool spans the object's heights; a fuse
  of equal height ranges, or whose slabs all have one profile). S9a.2: the
  other stacks, general bodies whose caps between slabs are the regions
  where consecutive profiles differ (implemented). S9b: prisms of line profiles in any relative position
  (polyhedra: exact plane arrangements). S9c: arc walls in any position
  (cylinders against planes and each other, S7 and D13's curves). S9d:
  cones, spheres and tori. Splines join each sub-step when their pairwise
  intersections are exact (S9a.2: spline profiles, spline/line from S8b.3,
  spline/conic from S7d.2, spline/spline by resultants).
* **2D Booleans (S9a).** The two profiles' boundaries are arranged
  exactly: every crossing of two segments is decided on the stored data
  (lines by rational crossings, a line and an arc by a quadratic surd, two
  circles by their radical line), rounded to binary64 as a new vertex; a
  crossing within the resolution of a vertex, or two crossings within it
  of each other, is `Degenerate` as in S8. Every piece of either boundary
  is classified against the other profile at a point strictly inside it
  (exact side predicates, or the certified `Profile::classify`): inside,
  outside, or on its boundary (a shared piece, with the same or the
  opposite direction). Fuse keeps each profile's pieces outside the other
  and the shared pieces of the same direction, cut keeps the object's
  pieces outside the tool, the tool's inside the object reversed and the
  shared pieces of opposite directions, common keeps each profile's pieces
  inside the other and the shared pieces of the same direction; the kept
  pieces are traced into cycles (S8's rule: the next piece at a vertex the
  first clockwise from the incoming one), counter-clockwise ones outer
  boundaries and clockwise ones holes, each validated as a profile.
* **Degeneracies, exactly and from the first fixture.** Coincident
  boundary pieces (collinear overlapping segments, arcs of one circle),
  tangent arcs, a vertex of one on the other's boundary, identical
  profiles, one inside the other, disjoint and touching profiles, and
  caps at equal heights are fixtures, never deferred. Pieces thinner than
  the resolution are `Degenerate`; a result touching itself at a point or
  along an edge (two solids sharing an edge) is `Degenerate` until the
  kernel holds non-manifold bodies.
* **History.** An input face, edge or vertex kept whole keeps its id
  (`Unchanged`, or `Modified` when its geometry or bounding ids changed);
  one kept in parts is `Split`; faces and edges of both inputs lying on one
  another (coplanar caps, coincident walls) are `Merged` into one result
  entity; an input entity not kept is `Deleted`; new edges and vertices
  where the inputs' faces meet are `Generated` from the faces they lie on.
  The region is `Merged` (fuse), `Split` or `Modified` as the solids are.
  The independent history check applies unchanged.
* **Evidence first** per sub-step: an independent reference (mpmath: the
  result's volume, area and centre by slicing both solids at every height
  and clipping their sections exactly, in Fractions for lines and by
  Green's theorem for arcs), fixtures of every class above in both frames,
  and a native `BRepAlgoAPI_Fuse`/`Cut`/`Common` capture (the result's
  validity, solids, counts, volume, area and centre) before the kernel
  code; then the kernel, a probe and `compare_boolean.py` with
  fingerprinted reviews, a `boolean` fuzz target, and the upstream
  `boolean` group's self-contained cases through the DRAW adapter
  (`bfuse`, `bcut`, `bcommon`, `bop`/`bopfuse`...).
* **S9a.2's spline profiles, decisions recorded before their code
  (2026-09-28).** Either profile may hold S8b's nonrational B-spline
  segments (degree at most 7, open, simple). Meetings are exact on each
  Bézier span: a line by the roots of its equation on the span (S8b.3's),
  a circle by the roots of `(x - c_x)^2 + (y - c_y)^2 - r^2` on it (degree
  `2p`, S7d.2's substitution), another span by resultants: each span's
  implicit equation (its Bézout matrix's determinant, exact) on the
  other's parametrization gives a polynomial of degree `p q` in that
  parameter, both ways, and each root in one span is paired with the one
  root in the other whose certified points' boxes meet (two candidates, or
  none, is `ComputationLimit`). A root of even multiplicity where the two
  curves stay on one side of each other is a tangency and cuts nothing
  (S9a.1's amendment (b)); a meeting within the resolution of a stored
  vertex is that vertex, and one within it of another meeting is
  `Degenerate`, as S9a.1's. Two spline segments share pieces only when they
  are one curve (equal degree, knots and poles, in either direction); any
  other overlap (an identically vanishing resultant, or a line or circle
  containing a spline's span) is `OutOfDomain`. A spline piece between two
  cuts is its segment's restriction (S8b.3's knot insertion), classified at
  the fraction `0.4453125` of its parameter range by the certified
  `Profile::classify`, and a traced result joins consecutive pieces of one
  segment into its restriction between their outer ends (a boundary keeps
  no vertex where it does not turn). In stacks a spline piece's wall is
  S8b's degree-`(p, 1)` surface, joined with the walls of the same piece
  across slab heights and with those of the adjacent pieces of the same
  segment across their shared end. History follows S9a's rules, supports
  by S8b's spline arms (`arcs_within`). Evidence first: the reference
  extended to spline boundaries (each atom's area and moments by Green's
  theorem over the classified pieces, the spline meetings by 40-digit root
  finding independent of the resultants), spline fixtures of every class
  above (a spline crossing a line, an arc and another spline, a tangency, a
  shared spline, a spline hole, stacks with spline walls, both frames) and
  a native capture before the kernel code; then the kernel, the probe and
  `compare_boolean.py`, the `boolean` target's spline profiles and a
  campaign.
* **S9b, decisions recorded before its code (2026-09-28).** Polyhedra in
  any relative position: two solids each a prism of a line profile
  (polygons with polygon holes) in any frame, or a Boolean's stack or a
  plane's piece of such a prism (every face planar, every edge a line).
  * **Exact models.** Each input is decided on its construction's exact
    model, never on its rounded vertices: a point of a prism is `o + u x +
    v y + w n` in rationals from the frame's stored binary64 origin and
    axes and the profile's and heights' binary64 values, a face's plane the
    exact plane through its model points (a cap's normal `n`, a wall's `(q
    - p) x n` for its segment `pq`), a plane piece's cut plane its stored
    exact coefficients; so every model vertex lies exactly on its faces'
    planes. An input without such a model (an imported body, a sheet) is
    `OutOfDomain` until S9c's general faces.
  * **Arrangement.** Every face of one input is intersected with every
    face of the other exactly: two non-parallel planes meet in a rational
    line, clipped to each face's polygon by exact parameters along it; a
    face's pieces are the exact 2D arrangement of its polygon with those
    segments in its plane (projected on its normal's largest coordinate
    plane, rational throughout), coplanar faces arranged with each other's
    polygons. Each piece is classified against the other solid at a
    rational interior point (exact ray parity along a direction chosen off
    every edge, or on the other's coplanar face with the same or opposite
    orientation) and kept as S9a keeps pieces; kept pieces sharing a plane
    and orientation join into maximal faces, collinear edges into one, as
    OCCT's unified result.
  * **Rounding.** Result vertices are the exact model's points rounded to
    binary64, faces' planes their exact planes rounded (origin and normal),
    each solid validated with measured enclosures; a rounding the
    validation refuses is `PrecisionLoss`. Degeneracies are S9a's: a
    result touching itself along an edge or at a vertex, a piece thinner
    than the resolution (measured on the exact model), `Degenerate`.
  * **Results.** General bodies holding their exact model (the inputs'
    constructions, the operation and the solid's index), so a rigid motion
    rebuilds them in the moved frames as S8's pieces and S9a.2's stacks do,
    and a result is an input again; histories by S9a's rules (a face
    continues the input faces whose plane and orientation it shares and
    whose polygon it overlaps, an edge or vertex those it lies on, new
    edges and vertices generated from the faces meeting there).
  * **Evidence first.** A reference in Fractions independent of the
    arrangement: each profile cut into convex pieces (ear clipping with
    exact orientation), each prism a union of convex prisms, and every
    pair's intersection by exact half-space clipping, giving the common's
    volume and moments exactly; fuse and cut by inclusion and exclusion;
    areas by clipping each face polygon against the other solid's convex
    pieces (inside, outside, or on a coplanar face with either
    orientation); solid counts by union-find over convex cells sharing
    positive area. Fixtures of every class (rotated boxes, a tilted tool,
    coplanar faces of either orientation, an edge on a face, edges crossing,
    a vertex on a face, identical, nested, touching and apart solids, a
    stack and a plane piece as inputs, both inputs rotated), a native
    `BRepAlgoAPI_Fuse`/`Cut`/`Common` capture before the kernel code, then
    the kernel, the probe and `compare_boolean.py`, the `boolean` target's
    rotated frames, the DRAW survey's 647 cases of frames with different
    axes and a campaign.
* **S9c, decisions recorded before its code (2026-09-28).** Arc walls in
  any position: prisms whose profiles hold arcs and circles as well as
  lines (and S9b.1's polyhedra), in any relative position.
  * **Sub-steps.** S9c.1: every pair of faces meets in lines, circles or
    ellipses or not at all (S7a's items): plane/plane, plane/cylinder, and
    two cylinders parallel (lines), coaxial (circles) or of equal radii with
    crossing axes (two ellipses). S9c.2: two cylinders meeting in S7b.1's
    procedural curves (D13). A pair of a later sub-step is `OutOfDomain`.
  * **Faces in their parameters.** Each input face is split in its own
    parameters (a plane's frame, a cylinder's cover `(u, v)`) by its
    boundary and every curve where the other input's faces meet its
    surface (S7a's items, restricted to the other face by exact or certified
    tests): on a plane lines, circles and ellipse arcs, on a cylinder lines
    (generatrices), circles (`v` constant) and a plane's section (S8a.2's
    sinusoid `v = a0 + a1 cos u + a2 sin u`). The crossings of two such
    curves in one face are S7d's certified curve/curve points (a conic pair,
    a sinusoid and a line or another sinusoid by the exact equation in
    `cos u`, `sin u`), each rounded once to binary64 on both faces' curves,
    as S9a's arrangement rounds its crossings. Each piece is classified
    against the other input at a certified interior point (exact for
    planes, the certified `classify` for cylinders), kept by the set
    function as S9b.1 keeps fragments, and pieces sharing a surface and
    orientation are joined into maximal faces, their edges joined where they
    run straight on (one line, circle or ellipse) between the same two
    faces.
  * **Geometry.** An edge's 3D curve is S7a's canonical item rounded once
    (a line segment, a circular arc, an ellipse arc; an ellipse's frame
    sharing its plane's `x` axis when it lies in an input plane, as S8a.2's
    cut faces), its pcurves the exact forms above; vertices rounded; every
    solid validated with measured enclosures, `PrecisionLoss` where the
    rounding fails it. Tangencies (a plane tangent to a cylinder along a
    generatrix, two cylinders tangent) and results touching themselves are
    `Degenerate`, as in S8 and S9a.
  * **History.** S9a's and S9b.1's rules: a face continues the input faces
    its pieces lie on facing their way, an edge or vertex those it lies on,
    new ones are generated from the faces meeting there.
  * **Evidence first.** An independent reference: each result's volume,
    area and centre from slicing the two solids in parallel planes (each
    slice's region bounded by segments and ellipse arcs, its area and
    moments by Green's theorem in closed form at 40 digits, the slices'
    breakpoints where the section's structure changes found as roots,
    integrated by Gauss-Legendre between them) and areas from the surfaces'
    own parameterisations over each result face's region in its parameters,
    checked against closed forms (a cylinder cut by planes, Steinmetz solids
    of equal radii, coaxial and parallel cylinders), inclusion and exclusion,
    and sampling; fixtures of every class (a tilted cylinder through a box,
    a box corner in a cylinder, perpendicular and oblique equal cylinders,
    parallel and coaxial cylinders, a cylinder through a hole, planes along
    a generatrix and tangent cylinders declared `degenerate`, both frames);
    a native `BRepAlgoAPI` capture before the kernel code; then the kernel,
    the probe and the comparison, the `boolean` target's arcs in turned
    frames, the DRAW survey and a campaign.
* **S9c.2, decisions recorded before its code (2026-09-28).** Two
  cylinders meeting in S7b.1's procedural curves (D13).
  * **Sub-steps.** S9c.2a: cylinders in exact frames (S9c.1's circular
    measure: the world's axes permuted or reversed), so crossing axes are
    perpendicular, of any radii and offset, whose section meets no cap's
    circle: every vertex is then a quadratic surd (a generatrix meets a
    cylinder in a quadratic), only the section curves between vertices are
    quartics. Where the section crosses a cap's circle the vertex is a
    nested surd (the cap's plane meets the other cylinder in generatrices at
    a surd offset, the circle meets those in `sqrt(q + q' sqrt k)`): S9c.2b's,
    `OutOfDomain` in S9c.2a. S9c.2b: cylinders in
    turned frames (affine models not circular in a common measure): a
    circle meeting the other cylinder, parallel cylinders' generatrices at
    an ellipse's and a circle's crossings, and the quartic curves; its
    points are algebraic numbers of degree up to four held as isolating
    intervals of their exact polynomials, refined on demand, two not
    separated at `1e-30` of the case's size `Degenerate` (S9a's crossings
    within the resolution). A pair whose curve comes within the resolution
    of a node (nearly equal radii, nearly crossing axes: the stored turned
    frames of `steinmetz_oblique` and `steinmetz_tilted` make the two
    ellipses of S9c.1 a quartic pair of loops about `1e-8` apart) is
    `Degenerate`. S9c.2b's decisions are refined by its evidence before its
    code. A pair of a later sub-step stays `OutOfDomain`.
  * **Classes, exactly.** In the frame `(a, b, e = a x b)` of the two axes
    `a` and `b`, a point of both walls has `eta = p.e` with `beta^2 =
    rA^2 - (eta - eA)^2` and `alpha^2 = rB^2 - (eta - eB)^2` (`beta`, `alpha`
    its coordinates across `A`'s and `B`'s axes along `b` and `a`), so the
    curve's class comes from the two intervals `[eA - rA, eA + rA]` and
    `[eB - rB, eB + rB]`: one strictly inside the other, two rings about
    the thinner cylinder (each a graph over its angle, whole turns); the
    two overlapping in part, one loop; an end of one on the other's end or
    interior end (a tangency: outside, or inside with a node, the
    figure-eight) `Degenerate`; apart, none; equal intervals with crossing
    axes are S9c.1's ellipses.
  * **Curve.** `Curve3::Meet`: the meeting of two cylinders as a graph over
    the first's angle (S7b's ruled parameterisation): the ruling at `u =
    start + sweep f` meets the second cylinder at `v = (-B + s sqrt(B^2 -
    A C)) / A`, `s` the branch. An edge's range keeps `B^2 - A C > 0`
    strictly (no turning point) so it is analytic: its jets, D13's
    quadrature and the tessellation's interval second derivatives apply as
    to S8d.3's sections; its pcurves are `Projection`s (on its own carrier
    evaluated by its angle); the `.brep` writer refuses it (`Unwritable`)
    until D13's interchange approximation, as sections. A turning point of
    one cylinder's graph (its ruling tangent to the curve) is regular on the
    other's, so a loop is four graphs: over `B`'s angle about the two
    points where `alpha = 0`, over `A`'s about the two where `beta = 0`,
    switched at rational points of the thinner cylinder's angle between
    them (`beta` then a surd of one field). Rings without a vertex are split
    at a rational point, as S9c.1's circles, tried again elsewhere when a
    meeting falls there. Switch and split vertices stay in the result
    (OCCT's intersection splines are split elsewhere: reviewed counts).
  * **Everything else is S9c.1's:** vertices where an edge meets a face,
    positions on a graph by its carrier's angle (pairs of surds), pieces
    traced in each face's parameters and kept by exact membership, joins,
    rounding, validation, histories; tangencies and results touching
    themselves `Degenerate`.
  * **Evidence first.** Fixtures of every class in exact frames (a thin
    pipe through a thick one with crossing and offset axes: two rings;
    a partial bite: one loop; equal radii offset: one loop; a pipe ending
    inside the other; a hole's wall crossed by a pipe; internal and external
    tangency declared `degenerate`) and S9c.2b's in turned frames (oblique
    and skew unequal cylinders, a turned pipe through a box's round hole);
    the reference checked against closed forms (the common of crossing
    perpendicular cylinders, `4 int sqrt((rA^2 - (eta - eA)^2) (rB^2 -
    (eta - eB)^2)) d eta`, and Legendre's `8 rA / 3 ((rA^2 + rB^2) E(k) -
    (rA^2 - rB^2) K(k))`, `k = rB / rA`, for crossing axes) and Monte Carlo;
    a native `BRepAlgoAPI` capture before `solid/boolean/curved/
    procedural.rs` exists; then the kernel, the comparison, the DRAW survey
    and a campaign.
  * **Evidence (2026-09-28).** `generate_procedural_boolean_fixtures.py
    --check`: 28 fixtures (22 S9c.2a in exact frames, 6 S9c.2b in turned
    ones; 22 solid, 6 degenerate with reasons: internal and external
    tangency, every operation) from S9c.1's reference unchanged (its slicing
    takes any radii and offset); the common of perpendicular cylinders by a
    quadrature in `eta` (volume, moments, both walls' and caps' areas, so all
    three operations) and Legendre's form (2.2e-41 against the quadrature;
    `16 r^3 / 3` at equal radii) agree with the reference within 7.9e-41,
    turned frames' forms (over `sin phi`, lenses) within 5.0e-16, S9c.1's
    per-pair checks within 1.3e-40 and 3.2 standard errors, no near
    coincidence; the curved fixtures regenerate byte for byte. The capture
    `occt-boolean-procedural-preimplementation` (`compare_procedural_boolean.py`,
    the kernel `unsupported` on all 28): every result valid with the
    reference's solids, the degenerate ones included; 5 matches without a
    quartic edge, 23 reviewed: BRepGProp's default integration on faces
    bounded by OCCT's approximated intersection curves misses by up to
    7.8e-5 in exact frames and 3.0e-4 in turned ones, the same results
    measured adaptively within 2.7e-9; six solids' counts change when
    unified. S9c.2a's kernel next.
  * **S9c.2a implemented** (`solid/boolean/curved/procedural.rs`,
    `Curve3::Meet`): perpendicular cylinders in exact frames as the
    decisions describe; all 22 exact-frame fixtures as declared (16 inside
    the reference, the 6 tangencies refused), the 6 turned-frame ones
    `OutOfDomain`, every history checked, results deterministic and moved
    rigidly; `compare_procedural_boolean.py` 5 matches and 23 reviewed (the
    11 kernel results whose counts differ added to their measure reviews:
    the loops' switch vertices, OCCT's own splits of its rings), S9a's,
    S9b's and S9c.1's comparisons unchanged. `Curve3::Meet` was written
    while the evidence agent took the capture (the edge type only; the
    Boolean module came after the capture). Amendments, from the
    implementation: (a) a cap's circle meeting the other cylinder is a
    nested surd `x + y sqrt(E)` (`x`, `y`, `E` in one quadratic field),
    decided exactly against the edge's arc and the face and refused only
    where it lies on both; (b) the validator certifies a band's holes
    bounded by sinusoids and projections (the signed `+v` ray at every `u`
    alias; projections by certified pieces), which also fixes the DRAW
    survey's `ZC5` fuse and cut of S9c.1 (area 19.7221 as native); (c) a
    hole no piece holds (a sliver whose binary64 image turned the wrong
    way) is `Degenerate`, fixing the survey's `T7` and `Y2` (six cases,
    now refused as a piece thinner than the resolution); (d) the `boolean`
    fuzz target stands the tool on its side (an exact frame) by a byte's
    top bit, so its arcs reach S9c.2a. The campaign at `9932a232` was
    clean: 289 mutation executions after a 2,720 s replay of 1,205 inputs,
    35,566 edges, 468 MB peak, the slowest input 16 s under
    AddressSanitizer.
  * **DRAW survey of S9c.2a (2026-09-28, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 822 (804
    before), none failing (8 before), all registered (`bopfuse_simple/ZE3`,
    a sentinel until now, and 17 new). Of S9c.1's 48 cylinders S9c.2 takes,
    16 evaluate (`ZE3` to `ZE6` of the four `bop*_simple` grids: radius 0.5
    through radius 1 at right angles, two quartic rings; areas as native
    DRAW's to its six digits), 16 are a tangency (`ZE7` to `ZF1`, the thin
    axis moved 0.5 off: walls touching at a point) and 16 are S9c.2b's
    turned frames (`ZD9`, `ZE2`, `ZF2`, `ZF3`). S9c.1's failures:
    `bopfuse_simple` and `bopcut_simple` `ZC5` evaluate (band holes of
    sinusoid and projection pcurves certified), the six `T7` and `Y2` are
    refused as a piece thinner than the resolution (a hole no piece holds).
    Sentinels added: `bopfuse_simple/ZD9` (S9c.2b), `ZE7` (tangency), the
    six `T7`/`Y2`. The ledger does not change (`checkprops -s`,
    `checkshape`).
  * **S9c.2b refined, before its code (2026-09-29).** Two sub-steps.
    S9c.2b.1: cylinders in any frames whose axes cross (affine models,
    elliptic in the world), with every vertex still a quadratic surd (an
    input's generatrix against the other cylinder). Each cylinder's
    discriminant `D_K` over its angle is a quartic in `t = tan(theta / 2)`
    with rational coefficients; its real roots (the turning points of `K`'s
    graph) are isolated exactly (Sturm sequences, rational intervals).
    `D_K > 0` all round: two rings over `K`; otherwise each interval of
    `D_A >= 0` holds one loop, whose turning points of either kind are
    ordered along it by certified enclosures, and a rational point of `A`'s
    angle is placed between each adjacent pair of different kinds, the
    loop's pieces graphs over `B`'s angle about `A`'s turning points and
    over `A`'s about `B`'s. Every piece is then verified exactly: no root of
    its carrier's discriminant within its range (Sturm counts at the
    switch points' angles, surds of one field). A double root, or a local
    extremum of `D_K` whose value lies within `(A_K res / 2)^2` of zero (two
    branches within the resolution: a near node, the stored turned frames'
    Steinmetz pairs, or a loop thinner than the resolution) is `Degenerate`.
    A cap's circle against the other cylinder is a quartic in its `t`: no
    root within the arc's range (exact counts), or every root's point
    certainly outside the face (isolating intervals refined to `1e-30`,
    interval membership), else `OutOfDomain` (S9c.2b.2). S9c.2b.2: those
    points and parallel cylinders in turned frames (generatrices at an
    ellipse's and a circle's crossings) as certified algebraic numbers,
    refined by its evidence before its code. S9c.1's `steinmetz_oblique`
    and `steinmetz_tilted` fixtures become declared `degenerate` (near
    nodes: the reference integrates the stored models' pair of loops about
    `1e-8` apart, the kernel refuses them within the resolution).
  * **S9c.2b.1 evidence (2026-09-28).** `generate_turned_boolean_fixtures.py
    --check`: 15 fixtures of two cylinders in turned frames with crossing
    axes, every vertex a quadratic surd (13 solid: a bite and skew unequal
    cylinders as one loop, a pipe through as the object, a pipe ending
    inside and a tilted pipe across a box's round hole as rings; 2
    `degenerate` near nodes of equal radii) from S9c.1's reference and
    S9c.2's closed forms unchanged: the common over `sin phi`, the
    perpendicular pair with its cap and areas, Cavalieri and Legendre's
    form within 4.2e-16 (ideal frames), S9c.1's per-pair checks within
    6.9e-41 and 2.4 standard errors; the classes from the intervals as
    declared; no cap's circle crosses the other cylinder's model (checked
    numerically along each circle). The stored models of the declared nodes
    have their intervals' ends 2.4e-17 and 8.3e-17 apart (two rings about
    1e-8 apart); where a stored `x` or `y` lies along the common
    perpendicular exactly (`steinmetz_oblique`, `steinmetz_tilted`) they are
    equal exactly, a double tangency meeting in two conics, not the loops
    1e-8 apart stated above: `Degenerate` either way. The capture
    `occt-boolean-turned-preimplementation` (`compare_turned_boolean.py`,
    the kernel `unsupported` on all 15): every result valid with the
    reference's solids; the nodes match within 1.5e-10 (OCCT's section two
    ellipses), 13 reviewed: BRepGProp's default integration misses by up to
    3.3e-4, the same results measured adaptively within 9.1e-10; four
    solids' counts change when unified. S9c.2b.1's kernel next.
  * **S9c.2b.1 implemented** (`solid/boolean/curved/turned.rs`): cylinders
    with crossing axes in turned frames as the refined decisions describe;
    all 15 turned fixtures as declared (13 inside the reference, the 2 near
    nodes refused), the procedural set's `oblique` and `skew` inside the
    reference, S9c.1's `steinmetz_oblique` and `steinmetz_tilted` refused
    and now declared `degenerate` (the evidence's finding: their extents
    are equal exactly, a double tangency), every history checked, results
    deterministic and moved rigidly; `compare_turned_boolean.py` 2 matches
    and 13 reviewed (the 5 loops' counts with their switch vertices added),
    the procedural, curved and earlier comparisons without failure.
    Amendments, from the implementation: (a) the chart for a loop's first
    cylinder has its antipode where `D < 0` (an axis point, else a point
    between two roots), so every loop's `t` range is finite; (b) turning
    points of the second kind are placed on the loops from isolators
    narrowed by 160 bisections (their midpoints' binary64 views order the
    events; the exact verification of every piece is what certifies the
    result); (c) a cap's circle is refused only where a root of its quartic
    lies on the edge's arc, not yet where it lies off the other face
    (S9c.2b.2's interval membership). The campaign at `c28c633c` was
    clean: 360 mutation executions after a 3,137 s replay of 1,286 inputs,
    36,227 edges, 448 MB peak, the slowest input 21 s.
  * **DRAW survey of S9c.2b.1 (2026-09-28, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 822, as
    before, none failing, all registered; no status changes. Of S9c.2a's 16
    cylinders in turned frames none evaluates: `ZD9` and `ZE2` of the four
    `bop*_simple` grids (equal, axes meeting at right angles) are
    `Degenerate` by decision (a near node), `ZF2` and `ZF3` (parallel) stay
    refused for S9c.2b.2 under S9c.2's general message. `bopfuse_simple/ZD9`
    re-purposed as the near node's sentinel; `bopfuse_simple/ZF2` added for
    S9c.2b.2. The ledger does not change. The bridge self-test's S9c.2 gap
    (a turned cylinder beside another, in fact apart) is decided since
    S9c.2b.1, as native DRAW decides it; it is now checked by volume, and
    the gap is ZF2's parallel pair.
  * **S9c.2b.2, decisions before its code (2026-09-29).** Every vertex
    S9c.2b.1 refuses lies on a cap's circle at an angle whose half-angle
    tangent `alpha` is a real algebraic number of degree at most four (a
    root of the circle's quartic against the other cylinder; against a
    parallel one, of its section's conic): its coordinates, and its angle
    on any cylinder of either input, lie in `Q(alpha)`. The surds'
    coefficients generalize from the rationals to such a field: a number is
    `a + b sqrt(d)` with `a`, `b` in `Q` or in one `Q(alpha)` (polynomials
    in `alpha`, `alpha` a square-free defining polynomial and a rational
    isolating interval: `polynomial/real.rs`) and `d` rational, its sign
    exact (Sturm-Tarski at `alpha`, and S9c.1's tower rule for the surd).
    S9c.2a's nested surds are the case `alpha` = the circle's `t`. Two
    numbers of different fields (two such vertices compared along one
    curve) are ordered by enclosures refined to `1e-40` of the case's size;
    unseparated there they count as equal, which the arrangement refuses
    (`Degenerate`, two meetings within rounding). Parallel cylinders in
    turned frames meet in generatrices through their cross-sections'
    crossings (algebraic `t` on the first's circle), their points ordered
    by height. Evidence first: fixtures where a section crosses a cap's
    circle (exact and turned frames: a pipe through the rim of a thicker
    cylinder, a tilted pipe entering through a cap's rim, a box's round
    hole crossed at its rim) besides the existing parallel ones
    (`parallel_hole`, S9c.1's `parallel_cylinders`), and a native capture
    before `solid/boolean/curved/algebraic.rs` exists.
  * **S9c.2b.2 evidence (2026-09-28).** `generate_capped_boolean_fixtures.py
    --check`: 18 fixtures of two cylinders whose quartic section crosses a
    cap's circle within both faces, all solid: 9 in exact frames (`XY`
    against `SIDE`: a pipe across a cap's rim, two rings; a bite at the
    rim, one loop; a pipe ending partway through the wall; a round hole
    crossed at its top rim) and 9 oblique in turned ones (`TILTX`, `LEAN`
    and `TILT` pipes across the rim, entering through it into a blind hole,
    ending in the wall, crossing a hole's rim), from S9c.1's reference and
    S9c.2's `Perpendicular` unchanged: the exact pairs' closed forms (the
    common clipped by both inputs' caps, areas included; the hole by a
    circular segment less the hole's cylinder) within 6.4e-41, S9c.1's
    per-pair checks within 5.9e-41 and 2.9 standard errors; the classes
    as declared; every pair has a cap's circle crossing the other cylinder
    on its face (checked numerically along each circle, 20 crossings at
    least 2.0 inside the other face, slopes at least 0.25). Perpendicular
    pairs in stored turned frames were dropped: a cap plane parallel to the
    other's generatrices within rounding is a near coincidence of the
    slicing (breakpoints 1e-17 apart). The capture
    `occt-boolean-capped-preimplementation` (`compare_capped_boolean.py`,
    the kernel `unsupported` on all 18: `(S9c.2b)` in exact frames,
    `(S9c.2b.2)` in turned ones): every result valid with the reference's
    solids; 2 match, 16 reviewed: BRepGProp's default integration misses by
    up to 1.9e-5, the same results measured adaptively within 1.0e-9 in
    exact frames and 1.9e-8 in turned ones (a small result's approximated
    section, unchanged at `Eps = 1e-12`); eight solids' counts change when
    unified. S9c.2b.2's kernel next.
  * **S9c.2b.2 implemented** (`solid/boolean/curved/algebraic.rs`, surds
    over `Q(alpha)` in `curved/num.rs`): every cap circle's crossing is an
    algebraic vertex, parallel cylinders in turned frames meet in
    generatrices; all 18 capped fixtures inside the reference, the
    procedural set's `parallel_hole` and S9c.1's `parallel_cylinders`
    inside theirs, no cylinder fixture left `OutOfDomain`; every history
    checked, results deterministic and moved rigidly;
    `compare_capped_boolean.py` 18 reviewed (ten results' counts: the same
    faces as OCCT's unified result, edges and vertices where each splits
    its section curves), every other Boolean comparison without failure.
    The number type was generalized while the evidence agent took the
    capture (the fields only; `algebraic.rs` came after the capture).
    Amendments, from the implementation: (a) S9c.2a's nested surds and
    their refusals are gone: every cap circle's crossing, in exact or
    turned frames, is a root of the circle's quartic; (b) the antipode of
    the circle's chart base, where the quartic drops a degree, is a
    rational vertex; (c) the bridge self-test's three cylinder gaps are
    evaluated now (their values are the fixture comparisons'). The campaign
    at `7e1281ad` was clean: 314 mutation executions after a 2,925 s
    replay of 1,379 inputs, 37,038 edges, 510 MB peak, the slowest input
    30 s.
  * **DRAW survey of S9c.2b.2 (2026-09-28, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 830 (822
    before), none failing, all registered: `ZF2` and `ZF3` of the four
    `bop*_simple` grids (equal parallel cylinders, one moved and turned by
    -120 or 120 and 60 degrees), their areas native DRAW's to its printed
    digits and the closed forms' within 6e-15. No other status or reason
    changes; S9c.2's general message no longer occurs. The
    `bopfuse_simple/ZF2` sentinel evaluates and is registered so; with its
    reason gone, none replaces it; the other sentinels are refused as
    before. 749 cases stay refused. The ledger does not change.

* **S9d, decisions recorded before its code (2026-09-29).** Cones,
  spheres and tori against the prisms S9a to S9c take, and against each
  other.
  * **Sub-steps.** S9d.1: a sphere (S3's `sphere_*`: whole, a cap or a
    zone between two latitudes) against polyhedral prisms (every face a
    plane) in any position: every pair of faces meets in a line (planes) or
    a circle (a plane and the sphere), every vertex a quadratic surd (a
    line against the sphere; the sphere's circle against a plane meets it
    on a line through the circle's plane). S9d.2: a sphere against prisms
    with arcs (a cylinder and a sphere: S7b.1's procedural curve, a graph
    over the cylinder's angle, `Curve3::Meet`'s kind) and two spheres (a
    circle). S9d.3: cones (conic sections, S8d.2; procedural against
    cylinders, spheres and cones, S7b.2). S9d.4: tori (spiric sections,
    S8d.3; procedural otherwise). A pair of a later sub-step is
    `OutOfDomain`.
  * **Models.** Each input on its exact model (S9c.1's rule): a sphere the
    set `|p - o|^2 <= r^2` of its stored centre and radius, a zone's end
    planes at its stored heights along the frame's stored normal; a
    membership test exact (a quadratic and linear signs, pushes as S9c.1's
    at first order: along the sphere's tangent plane a push keeps to it).
  * **Faces in their parameters.** The sphere's face is traced in S3's
    `(u, v)` (longitude, latitude) as S9c.1 traces a cylinder's in `(u,
    w)`: split at a rational seam meridian (a half-plane through the axis
    at a rational point of the equator, each input's seam its own, tried
    again at another when a meeting falls on it); a pole is a vertex only
    where a whole sphere's piece holds it, a loop round a pole wound in
    `u`. Plane sections are circles (centre the centre's projection,
    radius `sqrt(r^2 - d^2)` a surd): a vertex where two meet is a
    quadratic surd, placed on each by its exact `(cos, sin)` about its own
    centre.
  * **Geometry.** A section is `Curve3::Circle` or `CircularArc` on the
    plane's frame (radius rounded once), its pcurves a circle on the plane
    and a `Projection` on the sphere; tangencies (a plane tangent to the
    sphere, a vertex of a prism on the sphere, an edge tangent to it) and
    results touching themselves `Degenerate`.
  * **History.** S9a's and S9b.1's rules.
  * **Evidence first.** A reference independent of the arrangement:
    slices by parallel planes (each slice a polygon's convex pieces
    against a disc, its area and moments by Green's theorem over segments
    and arcs in closed form, breakpoints where a slice passes a vertex or
    touches the sphere as exact roots, Gauss-Legendre between), areas from
    each face's own parameters, checked against closed forms (spherical
    caps and zones, a sphere and a half-space, octants, a box inside, a
    sphere inside), inclusion and exclusion and Monte Carlo; fixtures of
    every class in exact and turned frames (a sphere through a box face,
    a box corner in a sphere, a sphere through a box edge, a zone against
    a box, a box through a zone's cap, a sphere inside, touching declared
    `degenerate`); a native `BRepAlgoAPI` capture before
    `solid/boolean/curved/sphere.rs` exists; then the kernel, the
    comparison, the DRAW survey and a campaign.
  * **S9d.1 evidence (2026-09-28).** `sphere_boolean_reference.py` slices
    a sphere, cap or zone (the end planes `w = R sin(latitude)` in the
    frame's affine coordinates, heights as stored) and a prism ear-clipped
    into convex pieces by planes normal to the zone's axis: discs cut by
    the zone's lines against the pieces' sections, the operations apart by
    classified boundaries (Green's theorem over segments and arcs), the
    common again by clipping, breakpoints exact quadratics' roots; the
    sphere's face by Archimedes, planar faces and end discs in closed form
    in their planes, solids by convexity (exact distances to polyhedra, runs
    of faces' boundaries outside the sphere).
    `generate_sphere_boolean_fixtures.py --check`: 30 fixtures (25 solid, 1
    empty, 4 degenerate; 5 in turned frames), every class the decisions
    list: a box's face, corner (exact, `TILT`, `LEAN` against a `TILTX`
    sphere) and edge through a sphere, a post through a zone's flat end, a
    box through a zone's end and band, a hemisphere crossed through disc
    and dome, discs coplanar with a box's face of either orientation,
    cavities both ways, a turned bar and a slab cutting in two, an L prism
    severed, an octant; a tangent face, vertices on the sphere and a
    tangent edge declared `degenerate`. Caps and zones only in exact frames
    at the origin, where the affine end plane and the kernel's cap plane
    (through `o + h n`, normal to `n`) are one plane exactly; the reference
    takes the affine reading, as S9c.1's prisms' caps. Closed forms of every pair
    within 9.3e-40 in exact frames and 2.2e-16 in turned ones, inclusion
    and exclusion 9.2e-41, the area identity 3.3e-40, every face's classes
    1.8e-40, a second direction 1.6e-40, Monte Carlo 2.7 standard errors,
    no near coincidence. The protocol takes a sphere on either side (its
    identity rows; the native oracle a `sphere` row by
    `BRepPrimAPI_MakeSphere`). The capture
    `occt-boolean-sphere-preimplementation` (`compare_sphere_boolean.py`,
    the kernel `unsupported` on all 30, `OutOfDomain("a Boolean of a solid
    with curved faces or edges in any position (S9c)")`): every result
    valid with the reference's solids, all 30 match within 3.3e-9 (no
    review), three solids' counts change when unified. S9d.1's kernel next.
  * **S9d.1 implemented** (`solid/boolean/curved/sphere.rs`): spheres, caps
    and zones against polyhedral prisms in S9c's arrangement, as the
    decisions describe; all 30 fixtures as declared, every history checked,
    results deterministic and moved rigidly; `compare_sphere_boolean.py` 30
    matches, counts included, no review. Amendments, from the
    implementation: (a) a whole sphere's split axis is a generic rational
    direction moving with the seam (a cap's or zone's is its frame's): with
    the frame's axis a box's face through the centre put a pole on the box
    at every seam (the octant); (b) hemispheres are traced in their
    stereographic projection (the orthographic one compressed the rim below
    the sampled boundary's chord error, and a hole near it fell outside);
    (c) a section circle meeting no seam inside both faces is a ring, given
    one vertex at a rational point; (d) a result's sphere face whose loops
    wind once holds its pole as a vertex loop, as S3's caps do, and a
    meridian's or parallel's pcurve is a line (a loop through a pole needs
    no `u` there); (e) the validator decides a ray against a whole sphere
    (a box's cavity in a sphere), rays against sphere pieces staying
    undecided; (f) a Boolean result's bounds hold its spheres' boxes (edges
    do not bound a sphere's face; the classifier's prefilter refused a
    fused sphere's south pole); (g) the `boolean` fuzz target makes its tool
    a sphere, a cap or a zone by the spline byte's two top bits, and
    accepts S9d's refusals; (h) a sphere's face whose first loop runs as a
    hole (none winding) is the sphere less its loops, as a torus's (S8d.3),
    its integrals the whole sphere's less the loops': the first replay with
    sphere tools found a tilted U's corner at a sphere's centre fused into
    such a face (`fuzz/regressions/README.md`); (i) the history checker
    compares circles by centres, normals and radii within tolerance and a
    circle in a plane as an ellipse (the same replay: a rim split at the
    seam, a section's circle rounded from its basis, differ from their
    wholes' frames within rounding). The campaigns at S9d.1 and at
    S9d.2 (`0f2cec1c`) replayed their corpora (1,418 inputs) under
    AddressSanitizer without a failure but passed their startup hour before
    mutating (exit 124); the corpus replayed with debug assertions at each
    step instead, and S9d.3a's campaign (below) is the clean one.
  * **DRAW survey of S9d.1 (2026-09-29, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 831 (830
    before), none failing, all registered. Of the 101 Booleans of a sphere
    refused before as solids other than prisms none evaluates: 33 are a
    unit box whose corner is at a unit sphere's centre, its far corners on
    the sphere (`Degenerate` by decision: 21 reported as a meeting at every
    seam tried, the arrangement taking a vertex on a sphere's face for a
    seam conflict and retrying; 12, the box quarter-turned, as a tangency
    of its wall), 60 a sphere against a cylinder with the cylinder first
    (S9d.2's `OutOfDomain`), and 8 the same with the sphere first and its
    equator the cylinder's rim, refused before the pair (a meeting at
    every seam tried). The one sphere against a box whose vertices stay
    off it, `bopcommon_simple/ZP9`, was refused by the adapter (a
    `psphere` on a DRAW `plane`, recorded until now as a `pcylinder` on a
    plane); the adapter places a sphere on a plane's frame now, and the
    case evaluates, its area native DRAW's to its printed digits and the
    closed form's within 9e-16 relative. The `bcommon_simple/A1` sentinel
    is refused as a vertex on the sphere now; three sentinels are added
    (`bopfuse_simple/ZH5`, S9d.2; `bfuse_simple/A4`, a tangent wall;
    `boptuc_simple/ZH5`, the sphere first). 748 cases stay refused. The
    ledger does not change.

  * **S9d.2 refined, before its code (2026-09-29).** A sphere against a
    prism with arcs, and two spheres. Two spheres meet in a circle (their
    radical plane's section; tangent spheres `Degenerate`); a cap's circle
    meets a sphere where it meets the sphere's circle in its plane (two
    circles in a plane: a quadratic surd). A cylinder and a sphere whose
    centre lies on the cylinder's axis meet in circles of latitude
    (parallels of the sphere, rational heights); otherwise in S7b.1's quartic:
    on the cylinder's `(u, w)`, `cos(u - u0) = q(w)` with `q` quadratic in
    the height, so the curve is a graph over `u` away from its turning
    points (the ruling tangent to the curve, `D(u) = 0`) and over `w` away
    from its others (`q(w) = +-1`), each piece analytic: `Curve3::Meet`
    with the sphere as its other quadric over `u`, and a new graph over the
    height, `u = u0 +- acos(q(w))`, switched at rational heights between the
    two kinds of turning point, verified exactly as S9c.2b.1's pieces (Sturm
    counts of each carrier's discriminant). A sphere tangent to a cylinder
    (a double root) or within the resolution of it is `Degenerate`.
    Evidence first: fixtures of every class (two spheres crossing, nested and
    apart; a sphere through a cylinder with its centre on the axis and off
    it, rings and a loop; a stadium or a round hole against a sphere; a
    pipe through a sphere; tangencies declared `degenerate`) and a native
    capture before `solid/boolean/curved/spheres.rs` exists.
  * **S9d.2 evidence (2026-09-29).** `spheres_boolean_reference.py` slices
    a sphere or cap and a prism of lines, arcs and circles (or two
    spheres) in a rational affine chart of the slices (the prism's `(u,
    v)`): a sphere's section an ellipse (a circle in an exact frame), the
    prism's its profile's preimage, every boundary a segment or an arc `C
    + A (cos t, sin t)`; the operations apart by classified boundaries
    (two arcs meet at the unit roots of a quartic), Green's theorem over
    them exactly, the sphere's face by Archimedes; breakpoints exact roots
    (two sections tangent where `det(lambda Q1(s) + Q2(s))` has a double
    root in `lambda`, or concentric sections coincide); flat walls, caps
    and discs in closed form in their planes, cylindrical walls by their
    angle; solids by components followed through the slices.
    `generate_spheres_boolean_fixtures.py --check`: 33 fixtures (28 solid,
    1 empty, 4 degenerate; prisms in `TILT` twice, whole spheres in
    `TILTX` and `LEAN`), every class the refined decisions list: two
    spheres crossing, nested, apart; a pipe through a sphere (the ring, the
    ends), in `SIDE`, a box with a coaxial hole (five solids), a hemisphere
    and a pipe; a rod through the sphere (two rings), a bite (one loop), a
    cylinder ending inside it; a stadium; tangent spheres and cylinders
    declared `degenerate`. Closed forms (the lens as two caps, the ring `pi
    h^3 / 6`, a coaxial cap, the holed box, an off-axis cylinder by the
    lens of two discs along its axis) within 9.2e-41 in exact frames and
    4.2e-17 in turned ones, inclusion and exclusion 9.2e-41, the area
    identity 2.3e-40, every face's classes 1.8e-40, a second direction
    (the prism cut obliquely) 1.8e-40, Monte Carlo 3.1 standard errors, no
    near coincidence, cap circles at least 0.48 from tangency. No protocol
    or oracle change. The capture `occt-boolean-spheres-preimplementation`
    (`compare_spheres_boolean.py`; the kernel `unsupported` on all 33:
    `OutOfDomain("a sphere against a cylinder or a sphere (S9d.2)")`
    against a prism, S9b.2's `(S9c)` refusal for two spheres): every result
    valid with the reference's solids, 18 match, 15 reviewed (BRepGProp's
    default integration on approximated quartics, up to 2.4e-6; adaptively
    within 3.4e-9), seven solids' counts change when unified. S9d.2's
    kernel next.
  * **S9d.2a implemented** (`solid/boolean/curved/spheres.rs`): two
    spheres, and spheres against prisms with arcs where the cylinder's
    rulings meet the sphere all round (rings); 24 of the 33 fixtures as the
    reference, the 9 loops `OutOfDomain`, every history checked, results
    deterministic and moved rigidly; `compare_spheres_boolean.py` 12
    matches and 21 reviewed (ten kernel results' counts: OCCT's unified
    faces, edges and vertices split elsewhere), every other comparison
    unchanged. Amendments, from the implementation: (a) S9d.2 is split:
    S9d.2a takes rings (a coaxial sphere's parallels among them, as
    `Curve3::Meet` rings) and two spheres, S9d.2b the loops (graphs over the
    height, `u = u0 +- acos(g(w) / rho)`, exact where the cylinder's frame
    is, rational heights giving surd points) and circles of unequal axes;
    (b) `Curve3::Meet` and the arrangement's pieces take a sphere as the
    other quadric; (c) a whole sphere's split axes are two rows of a
    rational rotation, so its great circle's basis has equal axes and its
    points against a cylinder are surds over `Q(alpha)`; (d) only a real
    repeated root of a discriminant is a tangency (a coaxial sphere's is
    constant, its chart's quartic `D (1 + t^2)^2`); (e) the validator tries
    a band's hole at two more points of its first fin, as other holes; (f) a
    meeting on a hemisphere is a seam's only on the split, so a prism's
    vertex on a sphere is `Degenerate` at once (the DRAW survey's box
    corners). The campaign: see S9d.1's (its replay passed the startup
    hour) and S9d.3a's.
  * **S9d.2b implemented** (`solid/boolean/curved/spheres.rs`,
    `Curve3::Rise`): a sphere and a cylinder in an exact frame meeting in
    loops, pieces over the height about the rulings' tangencies and over the
    angle about the circles', switched at rational angles, each verified
    exactly; 32 of the 33 fixtures as the reference (the bite in `TILT`
    `OutOfDomain`: a turned cylinder's circle is not round in its model, so
    its height graph has no closed form), every history checked, results
    deterministic and moved rigidly; `compare_spheres_boolean.py` 12 matches
    and 21 reviewed. `Curve3::Rise`'s `phi` and `rho` are its binary64
    `atan2` and `hypot` of `(alpha, beta)`, constants of the curve (an
    interval `atan2` on the branch cut failed its jets). The campaign: see
    S9d.3a's.
  * **DRAW survey of S9d.2 (2026-09-29, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 851 (831
    before), none failing or timing out, all registered. Of the 60 spheres
    against cylinders refused before as S9d.2's `OutOfDomain` (`ZH5` to
    `ZJ3` of the four `bop*_simple` grids: a `pcylinder` of radius 4 and
    height 8, a sphere of radius 4, 2 or 6 centred on its top cap) 20
    evaluate: radius 6, the cap inside the sphere and a coaxial ring on
    the wall, unturned or turned by quarter turns; areas native DRAW's to
    its printed digits and the closed forms' within 9.1e-16 relative. The
    other 40 are refused for reasons the decisions did not name: radius 4
    (24, and the 8 of `boptuc_simple` refused before as a meeting at every
    seam tried), the rim on the sphere's equator, a tangency along a
    circle whose discriminant vanishes identically, as a computation limit
    (`negative_chart` gives up before the tangency test: `Degenerate` by
    the decisions); radius 2 turned a quarter turn about x (16), the
    cap's plane through the sphere's poles, as `PrecisionLoss` (a box's
    wall through an upright sphere's centre alike: S9d.1's domain too; the
    same sphere unturned, or tilted 30 degrees, evaluates). The 21 unit
    boxes with a corner at a unit sphere's centre are refused as a vertex
    of one input on the other's face now (amendment (f)). Sentinels:
    `bopfuse_simple/ZH5` re-purposed (the computation limit),
    `bcommon_simple/A1` updated, `boptuc_simple/ZH5` dropped (its new
    reason `ZH5`'s), `bopfuse_simple/ZI4` added (`PrecisionLoss`). 728
    cases stay refused. The ledger does not change. After the survey: a
    discriminant vanishing identically is tested first (`sphere_cyl`), so
    the 32 of radius 4 are a tangency between the inputs, `Degenerate`, as
    the decisions say (`ZH5`'s purpose updated). S9d.1's follow-up, done
    after S9d.3a: a section through a stored sphere's pole, or within the
    resolution of it (a sphere turned by a binary64 quarter turn), gets a
    vertex there (exactly on the section, at its rational direction nearest
    the pole) that assembly keeps, so a meridian's pcurves are lines on
    either side and the loop turns half a turn at the pole (the validator's
    `u` gaps scale by `cos v`); a section through a pole off the axis's
    planes is `OutOfDomain` (its pcurves' `u` undefined at the pole).
    `ZI4` evaluates, its area native DRAW's (`tests/sphere_booleans.rs`,
    `sections_through_a_spheres_poles`); the survey after S9d.3a's
    registered 12 of the 16 (`ZI4`, `ZI6`, `ZI7` of the four grids, the
    sentinel `bopfuse_simple/ZI4` re-purposed) and found `ZI5` wrong,
    registered after its fix (the DRAW survey of S9d.3a, below).

  * **S9d.3 refined, before its code (2026-09-29).** Cones and frusta (S3's
    `cone_*`: radii `bottom` and `top` over heights `0..h` on the stored
    frame) in S9c's arrangement. Sub-steps: S9d.3a against polyhedral
    prisms; S9d.3b against prisms with arcs, spheres and other cones (S7b.2's
    procedural curves, graphs over a cone's angle about its ruling's
    tangencies, `Curve3::Meet` with a cone carrier); a pair of a later
    sub-step `OutOfDomain`. The model: the set between the end planes where
    the point's distance to the axis is at most the radius interpolated
    linearly in height, exact on the stored axes (affine where they are
    not orthonormal); its wall split into two halves by a plane through the
    axis at a rational point (each input its own seam, retried), its apex a
    vertex (virtual unless the cone's own). A plane meets the wall in a
    conic that is a graph over the cone's angle: the ruling at `(cos u, sin
    u)` meets the plane at `v(u)`, a rational function of `(cos, sin)`, so a
    rational angle gives a rational point, a line meets the cone in a
    quadratic surd, and a conic is split where the plane is parallel to a
    ruling (a parabola's or hyperbola's ends at infinity lie beyond the
    frustum's end planes, which bound the pieces). A plane through the apex
    (lines) or tangent to the cone along a ruling is `Degenerate`, as S8d.2
    refuses them. Edges: the sections as S8d.2's `EllipseArc`,
    `HyperbolaArc`, `ParabolaArc` or `CircularArc`, their pcurves
    `Projection`s on the cone. Evidence first: fixtures of every conic class
    (a frustum cut by a box's faces in ellipses, parabolas and hyperbolas, a
    box through a cone's apex region, a cone inside a box, planes through the
    apex and tangent along a ruling declared `degenerate`, exact and turned
    frames), an independent reference and a native capture before
    `solid/boolean/curved/cone.rs` exists.
  * **S9d.3a evidence (2026-09-29).** `cone_boolean_reference.py` takes a
    cone or frustum on the decisions' model in its chart (the stored axes
    as rationals, the radius linear in the height, a zero radius an apex)
    and slices it and a prism (S9d.1's convex pieces) by planes normal to
    the axis, each slice a disc of radius `r(s)`, or tilted from it within
    the cone's complement, each an ellipse made a disc by coordinates
    orthonormal for the quadric's form: S9d.1's classification of the
    boundaries (a line within 1e-30 of tangency taken as tangent: a face
    tangent along a ruling touches every slice's circle), the wall by its
    element `r(w) N(theta)` over the arcs inside, breakpoints exact
    quadratics' roots (a face's line tangent to the section the
    discriminant of the quadric along it), planar faces by lines of
    constant height across them, end discs in closed form, solids by
    convexity (`K - P` the components of the sets `int K n H_f` joined
    where `int K n H_f n H_g` has volume). `generate_cone_boolean_
    fixtures.py --check`: 30 fixtures (25 solid, 1 empty, 4 degenerate; 6
    in turned frames), every class the refined decisions list: a frustum
    cut in an ellipse, a parabola clipped by both end planes and a
    hyperbola; a box through a cone's apex region and a box's corner below
    it; a cone inside a box, a box inside a frustum; a half-space normal to
    the axis; a slab and a bar cutting in two; coplanar discs of either
    orientation; a frustum in `TILT` below a box's face, a cone in `LEAN`
    with its apex in a box (an oblique cone), a frustum in `TILTX` inside a
    box and across a slab; a wall through the apex (lines), a face tangent
    along a ruling and a box's vertex on the wall declared `degenerate`.
    Closed forms of every pair (an aligned box by rectangles inside the
    sections' discs and the walls' hyperbolic chords, a half-space by
    circular segments, a slab as two, a frustum) within 1.2e-40 in exact
    frames and 1.2e-16 in turned ones, inclusion and exclusion 9.2e-41, the
    area identity 1.8e-40, every face's classes 9.2e-41, a second direction
    (every section an ellipse) 1.8e-40, Monte Carlo 2.7 standard errors, no
    near coincidence. The protocol takes a cone on either side (its
    identity rows; the native oracle a `cone` row by
    `BRepPrimAPI_MakeCone`). The capture
    `occt-boolean-cone-preimplementation` (`compare_cone_boolean.py`, the
    kernel `unsupported` on all 30, `OutOfDomain("a Boolean of a solid with
    curved faces or edges in any position (S9c)")`): every result valid
    with the reference's solids, all 30 match within 1.4e-8 (BRepGProp's
    default integration; adaptively within 5.0e-9, every section an exact
    conic), no review, thirteen solids' counts change when unified.
    S9d.3a's kernel next.
  * **S9d.3a implemented** (`solid/boolean/curved/cone.rs`): cones and
    frusta against polyhedral prisms, all 30 fixtures as the reference (26
    solid results within its enclosures, the 4 degenerate refused), every
    history checked, results deterministic and moved rigidly;
    `compare_cone_boolean.py` 25 matches and 5 reviewed (OCCT's conical
    faces keep their seams: one edge and one vertex more per solid), every
    other comparison unchanged. Amendments, from the implementation: (a)
    the wall is not split: it is a graph over the plane of `(u, v)` (`k`
    never zero), one face traced in that projection, each rim one closed
    edge with a vertex at a rational seam point (a meeting there tried
    again at another), the apex a point inside the face; a result's wall
    whose loops wind once closes at the apex, a vertex loop continuing the
    input's; (b) a section's rounding is S8d.2's (`cone_conic`, factored
    out of its setup), a plane normal to the axis a circle about it, and
    its pcurves on the cone `Projection`s (a circle about the axis a line of
    constant `v`); (c) a line's double root against the cone beyond the
    solid's ends is no meeting (a line through a frustum's virtual apex);
    (d) the validator casts rays against cones (hits on the surface's
    nappe, the apex one crossing more up the ray when the face closes at
    it) and clears a closing chord of no length within rounding by its
    start (a ring's one fin: the frustum's cavity), the history checker
    takes hyperbolas and parabolas in a plane's pieces, a rigid motion moves
    them, and a body's bounds hold every vertex and those conics' extremes;
    (e) the DRAW bridge's cone case (a box's edge along the cone's axis
    through its apex) is a tangency now: a cone standing in the box
    evaluates instead, and a cone against a cylinder is S9d.3b's; (f) the
    fuzz target's tool is a cone or frustum when its spline byte lies in
    `160..192`; (g) two parallel faces apart by no more than the resolution
    are one plane within it, `Degenerate` (the corpus replay with debug
    assertions: a frustum on a tilted prism's top, its origin rounded 2.7e-16
    off the top's plane, left a sliver whose loops' binary64 images could not
    be nested, an open result; `fuzz/regressions/README.md`). The
    campaigns at `6f51406c` (S9d.3a merged), all clean to the end of their
    600 seconds of mutation: `boolean` on U6's sampled replay (418 of 1,419
    inputs, 1,426 s; the full replay passes the startup hour under
    AddressSanitizer, and minimising the corpus did too, exit 124, the
    corpus unchanged), 217 mutation executions, 40,143 edges, the slowest
    input 56 s, 565 MB; `brep_validation` (the validator's rays against
    cones and chords of no length), 2,508 mutation executions, 20,344
    edges; `split` 1,373, 30,723 edges; `tessellation` 1,980, 19,809
    edges. No new artifact.
  * **DRAW survey of S9d.3a (2026-09-29, `UPSTREAM_TESTS.md`).** Of the
    boolean group's 1,802 self-contained cases Rust evaluates 887 (851
    before), none failing or timing out, all registered. Of the 138 cases
    of a cone refused before as solids other than prisms, 72 are a box of
    side 4 and a `pcone` frustum (`ZF5` to `ZH4` of the four `bop*_simple`
    grids) and 36 of them evaluate: a frustum on the box's central
    vertical standing on it, inside it or through its bottom face (circles
    on its faces), and a frustum of the box's height whose axis is 2
    outside a wall (hyperbolas on the walls across it; two turned 30
    degrees about the axis); areas and volumes native DRAW's to its printed
    digits and closed forms' within 8.1e-16 relative. The other 36 are
    `Degenerate` for what their geometry holds (a wall through the
    frustum's axis with the box's corners on a rim tangent there to two
    edges, or a base circle on the top face tangent to two of its edges),
    reported as the first the arrangement meets: a tangency between the
    inputs 24, a vertex of one input on the other's face 6, a plane
    through a cone's apex 4, an edge of one input meeting an edge of the
    other 2. The 66 cones against cylinders are S9d.3b's `OutOfDomain`.
    The 32 cylinders' rims on a sphere's equator are a tangency between
    the inputs, as S9d.2's correction says. The kernel's solids other than
    prisms are now 16 tori against cylinders (`ZL2` to `ZL5`); the
    adapter's 61 are unchanged. Sentinels added: `bopfuse_simple/ZJ4`
    (S9d.3b's `OutOfDomain`), `ZG8`, `ZG2`, `ZG4` and `boptuc_simple/ZG2`
    (the four `Degenerate` reasons), `bopfuse_simple/ZL2` (a torus); none
    changes. 692 cases stay refused. The ledger does not change. After
    the near-coplanar refusal and S9d.1's pole follow-up, a survey again:
    only the 16 `PrecisionLoss` spheres change (none is refused by either
    change); 12 evaluate, areas and volumes native DRAW's and the closed
    forms' within 1.1e-15 relative, and are registered; `ZI5` (the
    sphere's stored axis `-y`, reference direction `-z`) is wrong: its
    fuse, cut and common keep the areas but their volumes are off by `32
    pi / 9` (the common 27.9253 against `16 pi / 3`), and its `btuc`
    fails validation (`uncertified_shell_orientation`) where native DRAW
    evaluates; a box in place of the cylinder alike. Not registered
    then; Rust evaluated 902, 899 registered, one failing.
    The fix, after it: not the pole follow-up's, but S9d.1's own and older
    (a box's vertical edge along an upright sphere's axis gave a quarter
    ball 14.58 for `8 pi / 3`, at the commit before the follow-up too): a
    sphere face's closing chord at a pole between two meridians, its ends'
    `v` a rounding apart, was enclosed over its `u` hull (`sph_lines`), the
    volume's enclosure 26 wide and its midpoint (the reported volume) far
    off, the validator's orientation sign undecided (`ZI5`'s `btuc`). Such a
    chord now takes the parallel's exact value widened by the `v` spread
    (`MATHEMATICS.md`); `ZI5`'s four operations and the wedges give the
    closed forms (`tests/sphere_booleans.rs`,
    `wedges_through_a_spheres_poles`). The fixture tests held the
    reference inside each enclosure but not the enclosure narrow: every
    curved Boolean fixture test now requires volume and area enclosures
    within `1e-9` relative, all of them passing. After it (59455fd7) a
    survey again: only the 16 former `PrecisionLoss` cases differ from
    S9d.3a's, all evaluating; `ZI5`'s four agree with native DRAW's
    volumes and centres of gravity and with the closed forms (1.1e-15
    relative), and are registered: Rust evaluates 903, each registered,
    none failing, 676 refused. A volume audit, since the registrations
    check areas: all 903 registered cases of the group run again with
    `vprops` on both backends; Rust's volumes and centres agree with
    native DRAW's to its printed digits or 1e-6 relative except in eight,
    where native DRAW is off: seven of the crossed cylinders `ZE3` to
    `ZE6` (native volumes up to 7.8e-6 relative from the closed form,
    Rust's within 2.0e-15, native centres up to 4.0e-6 off the axes) and
    `bopcut_simple/A9` (native `vprops` of two slabs puts their centre at
    `(0.53125, 0.515625, 0.46875)`, its solids' own centres averaging to
    Rust's `(0.5, 0.5, 0.5)`).
  * **S9d.3b evidence (2026-09-29).** `cones_boolean_reference.py` takes a
    cone or frustum against a prism with arcs, a sphere, cap or zone, or
    another cone, each input its exact surfaces in world coordinates, and
    slices both in S9d.2's chart along a direction in which every quadric's
    section is an ellipse (a cone's axis where the other allows, else a
    prism's or a rational combination, checked exactly): a cone's section an
    ellipse cut by its end planes' lines, the operations apart as S9d.2's,
    breakpoints from every edge of one input against every surface of the
    other (a line's quadratic, a circle's resultant in its plane), planes
    tangent to sections, sections tangent and shrinking to points. Every
    face but the sphere's is swept by lines on it (a cone's rulings, a
    cylinder's and a flat wall's generatrices, parallel lines on a planar
    face), cut at every surface's roots and classified, breakpoints the
    roots of each surface's leading coefficient and discriminant along the
    lines and of every two surfaces' resultant; solids as S9d.2's, same-key
    components followed by their centroids, sections vanishing on both
    sides of a breakpoint not joined. `generate_cones_boolean_fixtures.py
    --check`: 40 fixtures (35 solid, 1 empty, 4 degenerate; 4 in turned
    frames), every class the refined decisions list: a coaxial pipe
    (circles), a cylinder through the wall off the axis, rods across a
    frustum and a cone's tip and around its apex, a box with a coaxial hole,
    a stadium; coaxial spheres (the apex region among them), a sphere off
    the axis (a loop), a dome; two cones tip in tip, frusta of opposite
    slopes, a frustum on another's disc, crossing cones, a cone in `LEAN`,
    parallel axes; a cylinder tangent along a ruling, an inscribed sphere,
    apexes touching and equal cones declared `degenerate`. Each expected row
    names its sub-step if S9d.3b is split as S9d.2 was: `S9d.3b.1` (36)
    where some input's rulings meet the other's quadrics transversally
    wherever their curve runs (circles, graphs over a cylinder's or a cone's
    angle), `S9d.3b.2` (4: spheres off the axis, the inscribed sphere) where
    both inputs' rulings are tangent to the other on their faces, by the
    rulings' exact discriminants. Closed forms of 22 of the 27 pairs (the
    lens of two discs along the axis for parallel cylinders, spheres and
    cones, coaxial or not; the strip of a rod across the axis; the holed box
    as S9d.3a's box less the core; equal cones) within 1.7e-40 in exact
    frames and 4.4e-17 in turned ones, inclusion and exclusion 1.2e-41, the
    area identity 1.1e-40, every face's classes 5.0e-41, the cone's wall two
    ways 3.6e-41, a second direction 4.6e-41 (19 pairs), Monte Carlo 2.9
    standard errors, no near coincidence, every edge and vertex at least
    0.0126 from tangency with or incidence on the other's surfaces. No
    protocol or oracle change. The capture
    `occt-boolean-cones-preimplementation` (`compare_cones_boolean.py`, the
    kernel `unsupported` on all 40, `OutOfDomain("a cone against a prism
    with arcs, a sphere or a cone (S9d.3b)")`): every result valid, the
    reference's solids in 39 (the inscribed sphere's cut one solid
    natively), 23 match (circles within 4.1e-14; two cones of equal
    half-angle with parallel axes meet in a plane, OCCT's hyperbolas within
    6.1e-10), 17 reviewed (BRepGProp's default integration on approximated
    quartics, up to 4.1e-6; adaptively within 6.0e-9; and that count),
    eleven solids' counts change when unified. S9d.3b's kernel next.
  * **S9d.3b refined, before its code (2026-09-29).** S9d.3b is split as
    the evidence names its cases. S9d.3b.1: a cone against a cylinder, a
    sphere or another cone where one input's rulings meet the other's
    surface twice all round: rings over that input's angle (S9c.2a's
    graphs), the carrier the first input (a cylinder before a cone) whose
    discriminant `D = B^2 - A C` is positive all round and whose ruling's
    leading coefficient `A` has no real root (a cylinder's `A` constant, a
    cone's a quadratic form in its angle, the ruling `b + k w` along `n + k
    (cos x + sin y)`, `D` then a quartic form); coaxial pairs are such
    rings too (their `D` constant). The other quadric takes a cone's term:
    `sum (g_i . p - e_i)^2 - (r + t (h . p - e_h))^2`, `t` the cone's slope
    (zero for a cylinder or a sphere), in the rings, a cap's or a rim's
    circle against it (`algebraic.rs`, a quartic's roots) and a sphere's
    circles against it (`spheres.rs`). Two cones whose quadrics differ by
    an affine function (parallel axes, equal slopes: `A` zero for every
    ruling) meet on that plane: S9d.3a's plane section of one, kept where
    the other holds it. `Curve3::Meet` takes a cone carrier and a cone as
    the other quadric (`half_angle`, `other_half_angle`, zero for a
    cylinder), its jets, bounds, projections and the validator's checks
    with them. A pair where no input's rulings meet the other all round
    (loops: the spheres off a cone's axis) is S9d.3b.2's, `OutOfDomain`; a
    repeated real root of `D` or of `A` where the curve runs is a tangency,
    `Degenerate`. The kernel in `solid/boolean/curved/cones.rs`.
  * **S9d.3b.1 implemented** (`solid/boolean/curved/cones.rs`): cones
    against cylinders, spheres and cones meeting in rings over a carrier or
    on a plane; of the 40 fixtures 32 are the reference's (its solids
    within the kernel's enclosures, one empty), 4 are refused as declared
    (`Degenerate`) and 4 left to S9d.3b.2 (`OutOfDomain`), every history
    checked, results deterministic and moved rigidly;
    `compare_cones_boolean.py` 21 matches and 19 reviewed (the capture's
    measures, and entity counts: OCCT's seams and curve splits), every
    other comparison unchanged. Amendments, from the implementation: (a)
    `turned.rs`'s forms take any degree (a cone carrier's `D` a quartic
    form, its chart polynomial an octic), and a form is decided to vanish,
    or not to, by its chart polynomial (`cos^2 + sin^2 = 1` unreduced in
    its terms); (b) a carrier whose `A` has real roots (its rulings along
    the other's asymptotic directions, branches to infinity: a cone in
    `LEAN` against a frustum, `cones_lean`) is S9d.3b.2's, as the decisions
    said, though the evidence counted it S9d.3b.1's; (c) `A` within rounding
    of zero (`rod_ruling`: a cylinder 1e-16 off a cone's ruling) and `D`
    vanishing identically (`ball_inscribed`: a sphere inscribed in a cone)
    are `Degenerate` (the latter declared so, taken from S9d.3b.2); (d) a
    carrier's apex, or a frustum's virtual apex, on the other quadric is
    `Degenerate`; (e) a ring on a cone's wall gets a vertex at a rational
    angle (no seam crosses it); (f) the validator takes `Curve3::Meet`'s
    radius zero at a cone's frame's origin (its apex there: the corpus
    replay's cone with its apex at its base, four inputs; a regression and
    `a_cone_with_its_apex_at_its_base_carries_a_ring`); (g) a circle on a
    cone's quadric beyond its ends (a pipe's cap on the frustum's cone
    extended) is no meeting, not a tangency; (h) the DRAW bridge evaluates
    a coaxial pipe cut from a frustum, and guards a torus's refusal (DRAW
    cannot give a cone and a sphere of their own ids). Campaign: the boolean
    campaign at `29278966` (600 s, a sampled replay) clean, 921 runs, the
    slowest input 14 s under AddressSanitizer, no new artifact;
    `brep_validation` there clean too, 4,498 runs, the slowest input 16 s
    of its 20.
  * **DRAW survey of S9d.3b.1 (2026-09-29, `UPSTREAM_TESTS.md`).** Of the
    66 cases refused as S9d.3b's `OutOfDomain` (a `pcylinder` of radius 4
    and height 8 and a cone or frustum, `ZJ4` to `ZL1` of the four
    `bop*_simple` grids; `bcut_simple/G9`, `H3`), 36 evaluate on Rust, none
    wrong: coaxial frusta on the cap, inside the cylinder or through its
    caps (circles), and frusta whose axis crosses the cylinder's at right
    angles (two quartic rings); volumes, areas and centres of gravity the
    closed forms' and `cones_boolean_reference.py`'s within 6.5e-16
    relative (centres 2.0e-15), native DRAW's to its printed digits on the
    coaxial ones and up to 7.8e-6 relative off on the crossing ones
    (BRepGProp on approximated quartics). 28 are registered
    (`bopfuse_simple/ZJ4`, the sentinel for S9d.3b's `OutOfDomain`,
    re-purposed). `ZK7` and `ZK8` (a frustum of radii 6 and 1 across the
    cylinder, its wide end through both caps) are right but take 63 to
    118 seconds alone on the debug worker, past the contract's 30 a case;
    five timed out at the survey's 120 seconds, eight cases at once; none
    is registered. The other 30 are refused, none as S9d.3b.2's loops: a
    plane through a cone's apex (`ZJ5`: the bottom cap's plane through a
    frustum's virtual apex) 4, a tangency between the inputs (`ZK1`: a
    frustum's base rim the cylinder's top rim) 4, a piece thinner than the
    resolution (`ZK2` to `ZK4`: that frustum turned about its axis) 12, a
    plane within rounding of a cylinder's direction (`ZK9`, `ZL1`: a
    frustum's end disc across the cylinder after a quarter turn about y)
    8, and a stack with a cone's wall given to another Boolean (`G9`,
    `H3`, S9c) 2. The rounding ones are the adapter's: it turns a prism by
    quarter turns exactly but a cone by the kernel's rotation (a quarter
    turn's cosine 6.1e-17, the frustum's axis `(1, 0, 6.1e-17)`); `ZK1`
    turned a whole turn is refused alike. Rust evaluates 934 (903 before),
    931 registered, none failing, five timing out; 640 refused (676).
    Sentinel added: `bopfuse_simple/ZK9` (a plane within rounding of a
    cylinder's direction); the others are refused as before. The ledger
    does not change. After the adapter's exact quarter turns for every
    solid and S9d.3b.2, a survey again: `ZK9` and `ZL1` (8) evaluate,
    `cones_boolean_reference.py`'s within 1.4e-15 relative, and are
    registered (the `ZK9` sentinel re-purposed; `bopfuse_simple/ZK1`, the
    coincident rims, the new one); `ZK2` to `ZK4` are a tangency between
    the inputs as `ZK1` is; no case reaches S9d.3b.2's code; `ZK7` and
    `ZK8` take 23 to 34 seconds on the debug worker (60 to 92 before the
    exact turn; 1.7 to 4.9 on a release build) and stay unregistered.
    Rust evaluates 946, 939 registered, one timing out; 632 refused. The
    volume audit of every registered case again: Rust's values unchanged,
    the turned spheres the closed forms' within 7.8e-16; a full contract
    run holds but for two timeouts under a load average of 16, which pass
    run again.
  * **S9d.3b.2 implemented** (`solid/boolean/curved/cones.rs`,
    `spheres.rs`): a cone and a sphere in loops (S9d.2b's graphs over the
    height, the carrier's circle of radius `b + k w`: `u = phi +- acos(g(w)
    / (rho(w) rho))`, `Curve3::Rise` with the cone's `half_angle`; the loop
    builder takes any ruled carrier in an exact frame), and a carrier whose
    rulings reach the other's asymptotic directions (its `A`'s simple real
    roots, algebraic directions over `Q(alpha)`): its turn split there into
    open pieces of both branches, the finite one's switch at `w = -C / 2B`,
    every root of the ruling's quadratic taken as `C / (-B - s sqrt(D))`
    where `(-B + s sqrt(D)) / A` cancels (exactly, in binary64 and in the
    jets). `cones_lean` evaluates as the reference; `ball_side`'s sphere
    touches the frustum's base plane inside its disc, a tangency
    (`Degenerate` by S9d.1's rule; the evidence checked edges and vertices
    against the other's surfaces, not faces against faces), declared so in
    the generator and the expected rows now (the reference's rows and the
    native capture unchanged), so the loop is checked on the same pair with
    the sphere's radius 0.875
    (`a_sphere_off_a_cones_axis_meets_it_in_a_loop`, against
    `cones_boolean_reference.py`'s rows for it); `ball_tilt` (a turned
    frame) stays `OutOfDomain`, as S9d.2b's loops in turned frames. A
    loop's section pieces reach the arrangement (the first build dropped
    them: the cone and sphere found apart). The DRAW survey: `13bb636a`
    (`UPSTREAM_TESTS.md`, ZK9 and ZL1 registered). Campaign: the boolean
    campaign at `29278966` (600 s, a sampled replay) clean, 921 runs, the
    slowest input 14 s under AddressSanitizer, no new artifact;
    `brep_validation` there clean too, 4,498 runs, the slowest input 16 s
    of its 20.
  * **S9d.4 refined, before its code (2026-09-29).** Tori (S3's
    `torus_*`) in S9c's arrangement. Sub-steps: S9d.4a, a whole torus
    (the full tube and turn) against polyhedral prisms; S9d.4b, against
    prisms with arcs, spheres, cones and tori (S7b.3's certified traced
    curves); a v-segment or a wedge, or a pair of a later sub-step,
    `OutOfDomain`. The model: `(|p|^2 + R^2 - r^2)^2 <= 4 R^2 (u^2 + v^2)`
    in the stored frame's coordinates as rationals (affine where its axes
    are not orthonormal), membership exact at first order as the others'.
    Its wall is traced in the torus's own `(u, v)` in four patches, cut at
    a meridian through a rational point of the equator (each input its
    own seam, retried) and at the outer and inner equators (`v = 0`, `v =
    pi`): the seams are circles with rational centres, axes and points,
    the patches injective in `(u, v)`. A plane meets the wall in a spiric
    section: on `alpha = a cos u + b sin u` the tube's circle at `u` meets
    it where `r (alpha cos v + c sin v) = -(R alpha + d)`, a graph over `u`
    whose point at a rational `u` lies in one quadratic field (the root of
    `r^2 (alpha^2 + c^2) - (R alpha + d)^2`, a quadratic form in `(cos u,
    sin u)`), or over `v` symmetrically; S8d.3's classes decide (loops about
    the axis over `u`, about the tube over `v`, a cap's loop in two graphs
    of each kind switched at rational parameters between its turning
    points), each piece verified exactly (no turning point inside it). A
    plane tangent to the torus (a double root) or within the resolution
    of it is `Degenerate`. Edges are S8d.3's `Curve3::Section`, their
    pcurves `Projection`s on the torus. Evidence first: an independent
    reference (slices normal to the axis are annuli, two circles against
    the prism's convex pieces; areas from each face's own parameters),
    fixtures of every class (a box through the hole, bands, a cap cut from
    the tube's outside and its inside, loops about the tube, a box inside
    the tube, exact and turned frames, tangencies declared `degenerate`),
    and a native capture before `solid/boolean/curved/torus.rs` exists.
  * **S9d.4a evidence (2026-09-29).** `torus_boolean_reference.py` takes a
    whole torus on the decisions' model in its chart and a prism (S9d.1's
    convex pieces) and slices them two ways: normal to the axis, each slice
    an annulus (S9d.3a's classification on both circles, the operations
    their combinations, the wall by the circles' angles with the element `r
    rho / q`), and by the meridian half-planes about the axis, each section
    the tube's disc against a convex polygon (the volume by the cylindrical
    element, the wall in its own `(theta, phi)`); breakpoints exact
    polynomials' roots (quartics through Yun's square-free factorization:
    edges and faces' lines meeting the torus, a face's line tangent to a
    slice's circle; in `theta` the decisions' spiric quadratic in `alpha`),
    planar faces by their own parameters as S9d.3a's, solids followed
    through the meridians (the common's convex sections, `D - C` by the
    faces' outer half-planes, `C - D` by runs of the section's boundary
    outside the disc). `generate_torus_boolean_fixtures.py --check`: 35
    fixtures (30 solid, 1 empty, 4 degenerate; 6 in turned frames), every
    class the refined decisions list: a bar through the hole across the
    whole torus (loops about the tube) and a strip along the equator cut in
    three; a slab normal to the axis both ways and a half-space above the
    equator (loops about the axis); a cap cut from the tube's outside, a box
    through the hole cutting four, caps cut from the tube's inside by a
    box's vertical edges and a wedge's edge (contractible loops over two
    faces); a half-space through the axis (loops about the tube); a box
    inside the tube, the torus inside a box, a box in the hole; a torus in
    `TILT` above a plane of a Villarceau plane's inclination (one
    contractible loop), in `TILTX` across a slab between its saddle levels
    (two loops about the tube in each plane), a bar in `LEAN`; a face
    tangent along the top circle, a wall tangent to the inner equator, a
    face on a Villarceau plane and a box's vertex on the torus declared
    `degenerate`. Closed forms of 19 of the 20 pairs (an aligned box by
    rectangles inside both circles, half-spaces by circular segments of
    both) within 3.3e-40 in exact frames and 1.5e-16 in turned ones (the
    degenerate pairs 2.3e-29), inclusion and exclusion 1.1e-40, the area
    identity and every face's classes 3.0e-35 (1.3e-39 but for the vertex on
    the torus), the meridians as a second direction 1.6e-40, the wall two
    ways 2.1e-37, Monte Carlo 2.8 standard errors, no near coincidence
    outside the degenerate pairs, every other pair's vertices at least 0.18,
    faces' planes 0.031 and edges' crossings (a sine) 0.54 from tangency.
    The protocol takes a torus on either side (`encode_boolean_case`; the
    kernel's test support already read its rows; the native oracle a `torus`
    row by `BRepPrimAPI_MakeTorus(gp_Ax2, R, r)`; every older capture
    reproduces). The capture `occt-boolean-torus-preimplementation`
    (`compare_torus_boolean.py`, the kernel `unsupported` on all 35,
    `OutOfDomain("a Boolean of a solid with curved faces or edges in any
    position (S9c)")`): every result valid with the reference's solids, 16
    match (circles within 4.7e-16), 19 reviewed (BRepGProp's default
    integration on the spiric sections, which OCCT builds as B-splines,
    every section but circles: up to 6.3e-6; adaptively within 1.8e-9 but
    for a small cap, 1.2e-8, and the Villarceau plane, 2.0e-8), nine solids'
    counts change when unified. S9d.4a's kernel next.
  * **S9d.4a implemented** (`solid/boolean/curved/torus.rs`): a whole
    torus against polyhedral prisms, all 35 fixtures as the reference (30
    solid within the kernel's enclosures, one empty, the 4 degenerate
    refused), every history checked, results deterministic and moved
    rigidly; `compare_torus_boolean.py` 11 matches and 24 reviewed (the
    capture's measures, and entity counts: OCCT's seams and B-spline splits,
    the kernel's caps in four exact pieces), every other comparison
    unchanged. Amendments, from the implementation: (a) the tube's seams are
    the parallels at a rational angle `v0` and `v0 + pi` (moving with the
    retried seam), not the equators: a prism's face on the equatorial plane
    would lie along a seam at every retry; (b) a point's angles are placed
    exactly as unit directions, `u` along `(l_u, l_v) / rho` and `v` along
    `(rho - R, l_w) / r` with `rho = (|l|^2 + R^2 - r^2) / 2R`, rational in
    the point's field (a line's pierces are a quartic's roots, over
    `Q(alpha)`); (c) rings over `u` or `v` are decided from each
    discriminant's sign before any root is isolated (a plane normal to the
    axis has `D_v` a square's negative, one through it `D_u`, their double
    roots no tangency); (d) a torus face's loops lie on one sheet of the
    cover by their material's side (a reference loop winding in `u`, the
    others in the turn above it on a forward face, below it on a reversed
    one; in `u` likewise for loops about the tube); (e) the validator
    counts a ray's crossings of a whole torus (no loops) exactly, the
    positive roots of its quartic (a cavity in a torus's tube); (f) a
    rigid motion moves `Curve3::Section` edges, and the history checker
    takes them in a plane's pieces by points along them; (g) the fuzz
    target's tool is a whole torus when its spline byte lies in
    `144..160`; (h) the DRAW bridge evaluates a torus inside a box and
    guards a torus against a cylinder (S9d.4b). Campaign: the boolean
    campaign at `29278966` (600 s, a sampled replay) clean, 921 runs, the
    slowest input 14 s under AddressSanitizer, no new artifact;
    `brep_validation` there clean too, 4,498 runs, the slowest input 16 s
    of its 20.
    Performance, every result unchanged (the fuzz replay's 1,432
    boolean inputs hashed alike before and after): the corpus's slow torus
    inputs replay (debug assertions, no ASan) 1.7 to 3.3 times as fast
    (`17e131e3` 8.1 s to 3.3 s, `4274e084` 6.8 s to 2.1 s, `1b405929` 11 s
    to 4-5 s, `ad0dbe59` 8.1 s to 4.6 s) and CI's slow sphere input
    `fab20f09` 5.0 s to 2.3 s, from a Lehmer gcd for the curved Booleans'
    and the rational enclosures' arithmetic (`rational.rs`: `num_integer`'s
    Stein gcd is quadratic even against a small operand), bisection of
    isolators over one denominator, each algebraic field's narrowed
    isolators kept, a plane's spiric section found once for the four
    patches, surds' known fields not tested for squares again, binary64
    interval products from their two extreme corners where the factors'
    signs fix them (bit for bit the four corners'), and the mass integrands of a
    torus or sphere face sharing their `cos^a u sin^b u` jets; the
    certified mass integrals of the results' torus faces now take most of
    what remains.
  * **DRAW survey of S9d.4a (2026-09-29, `UPSTREAM_TESTS.md`).** The
    Boolean group holds no torus against a box, so no case evaluates
    newly: its tori are 16 `pcylinder`s with a coaxial torus (`ZL2` to
    `ZL5` of the four `bop*_simple` grids), refused as S9d.4b's
    `OutOfDomain` (a torus against a curved face) instead of S9c's, and
    three copies of one torus in `bopfuse_simple/ZP6`, the adapter's
    shared ids. No case reaches a segment or a wedge, none fails, and
    none times out (`bopcommon_simple/ZK8` evaluating within the survey's
    120 seconds this time, unregistered as the other `ZK7` and `ZK8`).
    Rust evaluates 947, 939 registered; 632 refused. The sentinel
    `bopfuse_simple/ZL2` now guards S9d.4b's `OutOfDomain`. The ledger
    does not change.
  * **CI on `rust-kernel` (2026-09-29).** The "Rust kernel" workflow had
    failed since S7 (last green `a608a3c2`), unnoticed; each cause:
    (a) `a_frustum_on_a_top_within_rounding_is_degenerate` failed on Linux,
    Windows and the 1.85 job: a frame normalized again divides by the
    platform's `hypot`, which rounds the tilted normal's length to 1 on
    macOS but not on every Linux or Windows runner, so the frustum's base
    was an ulp off parallel to the top and escaped 6f51406c's refusal of
    parallel planes within the resolution. Planar faces within the
    resolution of each other across their boxes' overlap, parallel or not,
    are now `Degenerate` (`a_frustum_turned_by_an_ulp_on_a_top_is_degenerate`
    turns the normal by an ulp on any host). The fuzz workflow's Linux-only
    crash `6d1fd061` (a frustum's caps an ulp off the prism's bottom and
    top, `loop_winding`) is the same and is kept as a regression. The
    `properties_baseline` assertion in those logs is the pre-T1
    regeneration step's expected one, not a failure. (b) The B-rep job
    stopped at captures without a Linux record: the boolean, split
    spiric, spline and sheet, and spline tessellation captures now have
    theirs (CI run 36556520655, each reviewed against the macOS capture in
    its `record.json`); the spline tessellation's Linux rows differ in
    nodes on 8 of 24 rows, so their reviewed differences need the next
    run's Linux fingerprints. (c) STEP-b's six reviews carried only macOS
    fingerprints: four Linux rows differ in their last bits, and have
    their own reviews now. (d) The upstream job's 40 minutes no longer
    covered the fixture generators' checks (31 minutes): they run in their
    own job in four groups, and the eight S9c-S9d Boolean generators,
    never checked in CI before, are the last two (the cones' reference
    alone took 31 minutes on a loaded host). (e) Fuzzing: the timeouts
    are slow units, a sphere (`fab20f09`, 13.6 s without AddressSanitizer)
    and a torus (`17e131e3`, 9.8 s) against prisms, and S7b's torus curve
    graph on an analytic input (`8dce95c9`, 42 s): the Booleans' number
    arithmetic 1.7 to 3.3 times faster (S9d.4a's note above; what remains
    is mostly the results' certified validation), and the torus curve
    graph's fold boxes stopping at a proved crossing (2.3 s, S7b.3b.1). (f) `compare_brep_io.py` has failed on every platform since
    F8's quadrature (first at `2f60b19f`): seven of the hammer's free
    spline faces now have enclosures narrower than OCCT's default
    BRepGProp error. OCCT's adaptive integration (Eps 1e-12) lies inside
    five of them and within the tool's 1e-9 of a sixth; on face 225 both
    of OCCT's integrations agree 2e-9 relative outside the enclosure: its
    pcurves leave a 2.5e-10 gap in the parameters, which the kernel closes
    by a chord and BRepGProp through the face's `UMin`. An independent
    40-digit Green integration (`free_face_measure.py`) lies inside the
    kernel's enclosure on all seven; they are reviewed differences now,
    each with its independent measure (`occt-brep-io-divergences.json`,
    OCCT's adaptive probe run only where a default measure falls outside).
    (g) With the generators moved out, the upstream job reached its runs
    for the first time since S8e and failed one: `bugs/heal/bug33171_1`
    (S8e's sentinel) fails on OCCT 7.6.3, Ubuntu's package and the Linux
    CI's native DRAW, whose unified shape is invalid (the bug the test was
    written for); a case may now state a native release's own outcome
    (`expected_occt_by_version`, a failure there recorded as that release's
    `known_failure`, the pinned release's expectation unchanged). The job's
    40 minutes did not cover the upstream runs (28 minutes) with the bridge
    tests and the native comparisons: 60 now. (h) With its records in
    place the split comparison ran on Linux for the first time since
    S8d.3: 26 of its reviews carried only macOS fingerprints (Linux ones
    added, the same differences), and `cap_oblique`'s review was stale on
    both hosts: since S9c.1 (`1c2f8afd`, found by bisection) the kernel's
    synthesized counts include a wound loop's seam split, OCCT's 3 faces,
    6 edges and 4 vertices, so only its measures still differ. The spline
    tessellation's Linux rows needed one review (`spline_dome/coarse`, the
    same understated deflection as on macOS). (i) The unoptimized test
    step (debug assertions, overflow checks) ran for hours with S9d's
    exact curved Booleans (`torus_segment_booleans` alone 459 s on this
    host): it keeps both checks optimized now (`CARGO_PROFILE_DEV_OPT_LEVEL`
    2, 40 s), and the test jobs have time limits. At `4d5d3d5b` every job
    and comparison passed, the B-rep job's in 60 minutes 3 seconds against
    its 60 (90 now); the fuzz workflow was green there too. Both
    workflows green at `b96730a7`, the first since `a608a3c2`. (j) The fuzz
    workflow's Linux-only `analytic_intersections` crash `6aaf4957` (at
    `428349e8`) is (a)'s `hypot` again, a plane through a torus's centre
    an ulp off containing its axis whose two loops lie within `2.5e-18` of
    the meridian angle `+-pi/2`, their ends' enclosures crossing with no
    binary64 parameter between them, so a procedural curve's loop now has
    disjoint end enclosures or the intersection is `ComputationLimit` (a
    limit of binary64 enclosures, not a degeneracy, which S7 decides only
    exactly; `toroidal.rs`'s tests take glibc's frames bit for bit).
    (k) The scheduled fuzz run 36716623883 (at `428349e8`) timed out on
    `analytic_intersections` `949542bf`, two congruent cones on parallel
    axes whose crossings of infinity (S7b.4's twin factor) lie three ulps
    from the subdivision point `pi/2`, where binary64 cannot sign `G` and
    the rational tier halved the rootless piece beside it to the floor:
    `roots_along` now excludes a piece on which `G_u` keeps a sign and
    the ends' certain signs agree, and its rational tier takes binary64's
    exclusions, signs and root brackets first and ends each root with
    interval Newton steps and bisection, every root the same adjacent
    binary64 pair as before (10.9 s to 0.8 s, `fuzz/regressions/README.md`).
    (l) The scheduled fuzz run 36868817257 (at `0dbd7c44`) found boolean
    `19624183`, a cut leaving a void beside a lens hole's spline walls that
    the validator's rays cannot decide (`uncertified_containment`), refused
    as `InvalidTopology` by S9a.2's stacks where the polyhedral and curved
    results say `ComputationLimit` (S9d.4b.1, S9d.4b.2a): every result now
    takes one rule (`undecided` in `solid/boolean.rs`, `tests/booleans.rs`).
    (m) The scheduled fuzz run 37008675181 (at `0dbd7c44`) stopped split at
    the 2 GiB RSS gate on `dfe03a75`, a torus's spiric section that alone
    takes 8 MB and 0.28 s (ulp-turned frames alike) but records 985,416
    distinct allocation stacks, 248 MB of AddressSanitizer's never-freed
    stack depot, with 50 MB live at the gate: split now keeps five-frame
    stacks as boolean does (`SHORT_STACK_TARGETS`; 1,400 corpus inputs
    1,350 MB RSS before, 434 MB after; no limit changed).
    (n) The scheduled fuzz run 36868817257's `brep_io` crash `faf66166`
    was the harness's own overflow, a number cast to `i64` (saturating at
    its end) and then nudged by `-2..=2`, no kernel value: the nudge
    saturates too, every check on the mutated text unchanged. (o) The same
    run's `degree_elevation` timeout `559a4800` (87 s under
    AddressSanitizer) raises a degree-19 periodic axis to 25, nearly all of
    it the production control map's `BigRational` blends: the map is
    computed in machine words on the knots moved affinely onto integers
    (the insertion ratios unchanged), again with `BigRational` on any
    overflow, the same exact map (9.6 s to 1.1 s in release, 262 s to 14 s
    under the sanitizer on a loaded Mac, `fuzz/regressions/README.md`).
  * **S9d.4b refined, before its code (2026-09-29).** Two sub-steps.
    S9d.4b.1: tori other than whole ones (S3's v-segments between two
    latitudes and wedges of a partial turn) against polyhedral prisms: the
    model S9d.4a's with the segment's end faces (at a latitude `v1`, the
    cone, cylinder or plane swept by the tube's point there; a wedge's
    meridian discs), its patches cut at the ends instead of a seam where an
    end lies. S9d.4b.2: a whole torus against prisms with arcs, spheres,
    cones and tori: the meeting of a torus with a quadric of revolution
    about another axis is a graph over the torus's `u` (the tube's circle
    at `u` against the quadric: a quartic in the half-angle tangent of
    `v`), or over the other surface's angle where it is ruled (its ruling
    against the torus: a quartic in the height), a point at a rational
    parameter algebraic of degree four (`Q(alpha)`, one generator per
    point, compared across fields by enclosures as S9c.2b.2's); pieces cut
    at their turning points' neighbourhoods by switches between the two
    kinds of graph, each verified exactly (the discriminant of the quartic
    in its parameter, Sturm counts); coaxial pairs meet in circles. The
    curve is D13's procedural cell, `Curve3::Toric` (the torus, the other
    surface, the parameter and the branch), its jets from the implicit
    function theorem on the two surfaces; a tangency or a double root is
    `Degenerate`. Evidence first per sub-step: an independent reference
    (S9d.4a's slices with the other solid's sections), fixtures of every
    class, a native capture before the kernel code.
  * **S9d.4b.1 evidence (2026-09-29).** One correction to the decisions
    first: S3 (`Topology::torus`, as `BRepPrimAPI_MakeTorus(gp_Ax2, R, r,
    low, high, angle)`) builds a v-segment as the region between the tube's
    arc and the axis, revolved, its end faces planar discs normal to the
    axis at the heights `r sin(latitude)` (`scaled_sin`, stored), each from
    the axis to the arc's end: not a cone, cylinder or plane swept by the
    tube's point, and not a band of the tube. So the model is S9d.4a's torus
    with the end planes `w = z` at the stored heights, the arc's ends where
    they meet the tube on the side of the latitude's cosine; in the meridian
    half-plane the section `M` holds a point with an odd number of the arc's
    points at its height to its right (the disc `[0, R +- q]` between the
    end heights, the tube's cap `[R - q, R + q]` beyond), its boundary
    counter-clockwise exactly when S3's wall is forward. The outer and inner
    halves' end planes `w = +-r` are tangent to the torus along their rings
    (the wall meets its discs tangentially along the rings, and a plane
    crossing such a ring meets the disc in a line and the wall in a spiric
    curve tangent to it at the ring); OCCT builds
    the inner half inside out (S3's capture), so the oracle reverses a part
    of negative volume. A wedge is the whole tube in the sector from the
    half-plane of `x` to that of the chart direction `(cos angle, sin
    angle)` rounded (the kernel's end disc): the half turn's end 1.2e-16 off
    the plane of `x`, the quarter's 6.1e-17 off that of `y`, so a prism face
    on those ideal planes is within rounding of an end (the fixtures keep
    coplanar faces to the exact start half-plane and the segments' end
    planes). `torus_segment_boolean_reference.py` slices both ways with `M`
    (normal slices: signed discs, a wedge's sector clipped, the common by
    clipping; meridian half-planes: `M`'s arc carrying all of Green's
    integrals, a wedge's section empty outside its turn), measures every
    face in its own parameters (the end discs exactly and a second way,
    `same`/`opp` for prism faces on their planes) and counts solids by
    sweeping the meridian half-planes (bands of height per section, samples
    between breakpoints joined one to one or refined, across a breakpoint
    through its own section, through the axis). 29 fixtures
    (`generate_torus_segment_boolean_fixtures.py --check`: 25 solid, 1
    empty, 3 degenerate; the outer half, the inner half, an upper band
    `0.5..2.25`, a half and a quarter turn; boxes and slabs across the tube,
    across an end disc, through the axis and the hole, a box standing on an
    end disc and one against a wedge's start disc, coplanar and opposite;
    `TILT` and `TILT2`): closed forms of 17 of 19 pairs within 8.2e-41 in
    exact frames (1.3e-16 turned, 4.3e-16 where a rounded end cuts the box),
    every operation two ways 6.4e-41, the wall two ways 7.0e-37, the end
    discs two ways 2.9e-41, the area identity (`+ 2 opp`) and face classes
    2.8e-40, Monte Carlo 3.3 standard errors, least margins 0.063 (faces),
    0.072 (vertices), sines 0.28 (edges) and 0.17 (rims). The native row
    `torus ... R r LOW HIGH ANGLE` for a segment or wedge
    (`BRepPrimAPI_MakeTorus` with latitudes and turn; a whole torus keeps
    its row, every older capture reproducing).
    `occt-boolean-torus-segment-preimplementation`
    (`compare_torus_segment_boolean.py`, keyed on
    `solid/boolean/curved/torus_segment.rs`, the kernel `unsupported` on all
    29 with `OutOfDomain("a Boolean of a torus segment or wedge (S9d.4b)")`):
    every result valid with the reference's solids, 16 match (lines and
    circles within 2.6e-16), 13 reviewed (BRepGProp's default integration on
    B-spline sections, up to 1.7e-5; adaptively within 5.7e-9 but for two
    results whose approximated sections bound a region 3.7e-8 and 3.2e-8
    off), nine solids' counts change when unified. S9d.4b.1's kernel next.
  * **S9d.4b.1 implemented** (`solid/boolean/curved/torus_segment.rs`):
    torus v-segments and wedges against polyhedral prisms, all 29 fixtures
    as the reference (25 solid within the kernel's enclosures, one empty,
    the 3 degenerate refused), every history checked, results deterministic
    and moved rigidly; `compare_torus_segment_boolean.py` 15 matches and 14
    reviewed (the capture's 13 measures, seven with entity counts: OCCT's
    B-spline sections split at its own points, the kernel's sections in
    exact pieces over `u` and `v`; and `ub_hole_cut`, a native match, by its
    counts alone), every other comparison unchanged. Amendments, from the
    implementation: (a) a segment's membership by the bands between its
    critical heights (its ends and `+-r`): in each the tube's outer and
    inner points at a height lie on the arc or not all along (decided once,
    exactly), so a point's side is exact in `t^2 - R^2`, the torus's quartic
    and its height (the tube's disc, `t < R + q`, `t < R - q` or nothing),
    pushed at first order as the whole torus's; a wedge's by the torus and
    its half-planes (their common within a half turn, their union beyond);
    (b) the wall in two patches, cut at the ends instead of a seam where an
    end lies (a segment's at its meridian seam, a wedge's at its parallel
    seam); a segment's rims its end planes' rings over `u` (S9d.4a's
    `TorusSec`; the halves' tangent `w = +-r` circles of radius `R`), a
    wedge's its end half-planes' rings over `v` (surds of `cos^2 + sin^2` of
    its rounded end direction), model edges now as conics and circles are,
    meeting a plane where the two planes' line meets the torus on the rim's
    branch (a repeated root on the other ring or circle no contact), rounded
    as the input's circles; (c) an inside-out segment's wall (the inner
    half) has its material outside the tube: its normals and its
    parameters' orientation turn; (d) the decisions' tangency is
    `Degenerate` on the part's wall only: a plane tangent to the torus, or
    within the resolution of it, off the wall is no contact, the part's
    section then its two branches over the wall's range (`v` for a segment,
    `u` for a wedge) where that graph's discriminant is positive all over a
    range holding the wall's exactly (Sturm counts, no extremum within the
    resolution of zero), or nothing where it is negative (`oh_side_common`:
    the box's wall `x = 1` touches the inner equator, off the outer half;
    `qw_box_*`: the wall `y = -1` touches it behind the quarter turn's
    start; the evidence's margins count only contacts on the wall); (e)
    `Solid::classify` decides segments and wedges
    (`decide::torus_part_location`: membership exact, the distances from the
    wall within its range and from the end discs in rational intervals),
    which every result of one needs (a result's classification asks its
    inputs'); (f) a torus face none of whose loops winds holds its holes
    inside its outer loop on the cover (a three-quarter wedge's fused wall);
    (g) a ray missing a whole sphere or torus crosses a face on it with
    loops no times, and a Boolean whose result the validator refuses only
    for an undecided containment (a cavity in a solid bounded by a sphere or
    torus face with loops, whose rays are not decided, as before) is
    `ComputationLimit`, not `InvalidTopology`; (h) the fuzz target's tool is
    a v-segment or a wedge when its spline byte lies in `148..160` (its bits
    2 and 3: the outer or inner half, a wedge of a quarter, a half, three
    quarters or 2 radians of a turn, or a band `0.5..2.25` or `-2.5..-0.25`,
    which by the flags' bits 3 and 4; `144..148` stay whole tori): 68 of the
    corpus's 1,419 boolean inputs now decode to parts; (i) the DRAW bridge
    evaluates a half turn inside a box. The corpus and the regressions
    replay with debug assertions without a failure, as do 6,579 variants of
    corpus inputs made parts (their spline byte and flags set), which found
    (f) and (g) (two kept as regressions, `fuzz/regressions/README.md`).
    Campaign: the boolean campaign at `85104dc3` (600 s, a sampled replay)
    clean, 1,040 runs, the slowest input 17 s under AddressSanitizer.
    DRAW survey: S9d.4b's, below (no case reaches a segment or a wedge).
  * **S9d.4b.2 evidence (2026-09-29).** Corrections to the decisions first:
    (a) their degree four holds for a quadric and for a torus about a
    parallel axis (its function along a tube's circle is of degree two in
    the tube's angle), but a torus about another axis meets a tube's circle
    in a trigonometric polynomial of degree four, a polynomial of degree
    eight in the half-angle tangent: two such tori's points at a rational
    parameter are algebraic of degree eight (`tori_ring`, `tori_cross`); (b)
    "a double root is `Degenerate`" means a singular point of the meeting
    (the surfaces tangent), not a graph's double root at a turning point:
    every pair but the coaxial ones and those apart has turning points (up
    to 12 in one face's family of curves; the stadium's events 34 over its
    faces), which the graphs' switches take; (c) two tori of parallel
    axes whose top, bottom or equator circles lie at one height and cross
    are tangent there (the first placement of `tori_side`, found by its
    quadrature's failure to converge): the kernel's tangency test must see
    it; (d) a prism whose flat wall meets its arc tangentially (a stadium)
    gives the torus's curves a surface switch on that edge, where the wall's
    and the cylinder's roots coincide to second order.
    `torus_curved_boolean_reference.py` takes each input's exact model (the
    torus S9d.4a's, the cone S9d.3a's, the sphere S9d.1's, a convex prism of
    segments and arcs S9d.2's) and sweeps every face of both by two families
    of circles and lines (the tori's tube circles and parallels, the
    sphere's meridians and parallels, a cylinder's and a cone's rulings and
    sections, a flat wall's generatrices and rows, a planar face's chords in
    two directions), the other input's surfaces along each curve a
    trigonometric polynomial of their degree (from `2d + 1` samples, their
    real roots the unit-circle roots of `z^d f`) or a line's polynomial,
    pieces classified inside and measured by the divergence theorem (area,
    volume and first moments in closed form but the area element's square
    root, 24-point Gauss-Legendre), the outer parameter integrated between
    events (where the curves' structure, their classes with the bounding
    surfaces, changes: a scan of 360 curves and bisection to 1e-35, 1e-20
    where two surfaces tangent along an edge make it rounding's to call); a
    tangency's double root split by rounding dropped; every quadrature
    node's structure its interval's. Solids by a sweep of the torus's normal
    slices (S9d.4a's annuli), rays from the axis joined where they overlap.
    44 fixtures (`generate_torus_curved_boolean_fixtures.py --check`: 10
    fuses, 16 cuts, 18 commons; 39 solid, 2 empty, 3 degenerate; the torus
    of radii 5/2 and 1): coaxial cylinders through the hole and around the
    tube, coaxial spheres about the centre and on the axis, a coaxial cone
    and a coaxial torus (circles); a rod across the torus (the common two
    solids, the rod less the torus three), a vertical cylinder and a frustum
    through the tube (graphs over their angle), a stadium across the tube, a
    sphere on the tube's top and one swallowing the tube, a cone in `LEAN`
    with its apex in the tube, tori linked apart, ringing the tube, of a
    parallel axis (the common two solids) and linked crossing the tube; a pin
    and a ball in the hole apart; the pipe and a sphere with the torus in
    `TILT`; `degenerate` a cylinder tangent along the outer equator, a
    sphere touching at a point and a coaxial torus touching along a circle.
    Checks: closed forms of 23 of the 27 pairs (one quadrature along the
    axis of the annulus against lenses of signed discs about a parallel
    axis, a rod's strip by rectangle-in-disc antiderivatives, sums apart)
    within 7.5e-39 in exact frames, 4.6e-17 turned; every operation two ways
    (the two families) 6.6e-35, both inputs from their faces 1.4e-39,
    inclusion and exclusion 1.1e-40, the area identity 1.4e-39, Monte Carlo
    2.8 standard errors, the degenerate pairs 3.8e-21 (closed forms
    3.1e-24); no near coincidence outside them, least margins 0.42 (the sine
    between surfaces where they meet), 0.066 (surfaces apart), 0.62 and
    0.031 (edges crossing the torus, and from tangency), 0.13 (the apex). No
    protocol row is new: the kernel's test support and the native oracle
    already build every input on either side (every older capture
    reproducing). `occt-boolean-torus-curved-preimplementation`
    (`compare_torus_curved_boolean.py`, keyed on
    `solid/boolean/curved/torus_curved.rs`, the kernel `unsupported` on all
    44 with `OutOfDomain("a torus against a curved face (S9d.4b)")`): every
    result valid with the reference's solids, 24 match (circles within
    1.6e-13), 20 reviewed (BRepGProp's default integration on B-spline
    sections, up to 9.8e-6; adaptively within 4.8e-9, unchanged at 1e-12),
    thirteen solids' counts change when unified. The kernel will need the
    torus's classification against cylinders, spheres, cones and tori, the
    graphs over `u` and over a ruled surface's angle with their switches,
    fields of degree eight for two tori, the three tangencies refused (along
    circles and at a point), linked solids in a fuse, and pieces across a
    stadium's tangent edges. S9d.4b.2's kernel next.
  * **S9d.4b.2 split, before its code (2026-09-29).** Two sub-steps by the
    evidence's corrections: S9d.4b.2a, a whole torus against prisms with
    arcs, spheres and cones (points algebraic of degree four, coaxial pairs
    in circles; `boolean-torus-curved` cases but the tori pairs), and
    S9d.4b.2b, two tori (degree eight about different axes, `tori_*`), its
    cases `LATER` in 2a's tests. In both, only a tangency of the surfaces
    (a singular point of the meeting) is `Degenerate`; a graph's turning
    point is a switch between graphs, not a refusal.
  * **S9d.4b.2a implemented** (`solid/boolean/curved/torus_curved.rs`): a
    whole torus against prisms with arcs, spheres and cones, all 32
    fixtures but the tori as the reference (29 solid within the kernel's
    enclosures, one empty, `pipe_equator` and `ball_touch` refused), the 12
    `tori_*` `LATER` (S9d.4b.2b's `OutOfDomain("a torus against a torus
    (S9d.4b.2b)")`), every history checked, results deterministic and moved
    rigidly (`tests/torus_curved_booleans.rs`, 33 s in the dev profile at
    opt-level 2); `compare_torus_curved_boolean.py` 19 matches and 25
    reviewed (the capture's 15 measure reviews of 2a's cases now with the
    kernel's entity counts, five native matches reviewed by their counts
    alone: OCCT's B-spline sections split at its own points and faces split
    at the torus's seams, the kernel's exact pieces), every other
    comparison unchanged. Amendments, from the implementation: (a) the
    graphs are over the torus's own `u` and `v` for every quadric, not over
    a ruled surface's angle: a sphere is not ruled, a curve is tangent to at
    most one of the meridian and the parallel at a regular point, and the
    quadric's function `G(u, v)` is of degree two in each angle, so a point
    at a rational parameter is algebraic of degree four either way; (b) the
    decisions' branch is a window of the other angle (rational directions)
    holding exactly one root over the piece's range, verified by Sturm
    counts on its ends, one isolated root inside and certified boxes clear
    of `G` or of its derivative in the other angle (no double root), as
    S7b.3b's tracks; (c) the components are traced numerically from the
    exact roots on a line of `u` between each two critical values (the
    discriminant in `v`, degree 24 in a chart of `u`) and every such root
    must lie on a verified piece, so none is missed; the symmetric fixtures
    put turning points of two branches at one `u` (the bore's top and
    bottom loops), which the discriminant's double roots do not refuse:
    each critical line is instead certified regular (boxes clear of `G`,
    `G_v` or `G_u`), a point where none clears a tangency within `1e-10`
    radians; a tube's circle on the quadric is `OutOfDomain`; (d) the
    switches lie where the meeting's slope in the angles is one, not midway
    between turning points of different kinds: a graph's series converges
    only as far as its own turning points, and switches near them made the
    validator's jets wide and a sphere's variant take 411 s; (e) coaxial
    pairs (`G` independent of `u`) are rings over `u`, stored as circles
    about the axis; (f) `Curve3::Toric` holds the torus, the other quadric,
    the parameter, its range and the window; its jets come from the
    implicit function theorem term by term (each coefficient linear in the
    new one, the functionals and `cos`, `sin` of the other angle continued
    by their recurrences), the base's root by interval Newton in its
    mean-value form over the base, its point from the root's change of sign
    in the window; rigid motion, tessellation rates, the history check and
    the writer's refusal as `Curve3::Meet`'s; (g) a cap's, a rim's or a
    whole sphere's great circle meets the torus at the roots of a polynomial
    of degree eight, a seam's tangency with a torus retried at another
    seam; (h) `Qd::to_f64` of a number whose 96-bisection enclosure is wide
    (a field's large coefficients) is its value at the isolator's middle,
    exactly, at more bisections until two agree (`ball_tilt_cut`'s places
    were off by a coarse enclosure, a vertex's angle on the wrong side of
    another); (i) a projection pcurve takes 64 or 256 anchors where 16 leave
    its lift unpinned (a section passing within 0.003 of a sphere's pole),
    and a Boolean whose result the validator refuses only for undecided loop
    winding (a torus face wound both ways, a torus knot) is
    `ComputationLimit`, both found by a fuzz replay (`tests/
    torus_curved_booleans.rs`, `a_section_near_a_spheres_pole_projects`,
    `a_torus_face_wound_both_ways_is_undecided`; not kept as fuzz
    regressions: 16 to 17 s each with debug assertions); (j) the fuzz
    target's object is a sphere or a cone against a whole torus by the
    flags' bits 5 and 6 (1 and 2): none of the corpus's 21 whole-torus
    inputs decodes differently, those against prisms with arcs now
    evaluating; (k) the DRAW bridge evaluates a coaxial pipe through a
    tube, guards a half turn against a cylinder, and the upstream sentinel
    `bopfuse_simple/ZL2` (a coaxial `pcylinder` through the tube)
    evaluates, registered so, none replacing it. The corpus and the
    regressions (1,435 inputs) replay with debug assertions without a
    failure, the slowest 5.5 s, as do 4,110 variants of corpus inputs made a
    whole torus against their prisms, spheres and cones (their spline byte
    and flags set), which found (d) and (i); the slowest variant, a ball
    about a torus in a turned frame, takes 26 s there (its sphere faces'
    certified integrals along the meetings). Native DRAW confirms the coaxial count reviews' reading:
    `pipe_hole_fuse`'s unified solid has 2 caps of one edge and 4 faces of
    three, the cylinder's two bands and the torus's band as two faces split
    on its outer equator (the kernel's 5 faces, 7 edges and 4 vertices with
    that band whole). Campaign: the boolean campaign at `f4ea7a2a` (600 s,
    a sampled replay) clean, 787 runs; its slowest input 41 s of the 60
    under AddressSanitizer (a frustum against a tilted prism, 3.2 s
    without: the exact arithmetic's allocations cost about twelvefold
    there). DRAW survey: S9d.4b's, below.
  * **S9d.4b.2b implemented** (`solid/boolean/curved/torus_curved.rs`): two
    whole tori, all 12 `tori_*` fixtures as the reference (10 solid within
    the kernel's enclosures, `tori_link_common` empty, `tori_kiss` refused),
    every history checked, results deterministic and moved rigidly
    (`tests/torus_curved_booleans.rs`, 95 s in the dev profile at opt-level
    2 on a loaded host, the parallel fixtures' certified integrals the most
    of it); `compare_torus_curved_boolean.py` 15 matches and 29 reviewed,
    the kernel supported on all 44 (the eight tori whose unified counts
    differ reviewed by them: coaxial and ringing pairs whose bands cross
    OCCT's seams, by native DRAW's faces, the others' B-spline sections
    split at OCCT's points), every other comparison unchanged and every
    other fixture's kernel rows bit-identical. Corrections to the evidence,
    from the implementation: (a) a tube's circle meets another torus in a
    trigonometric polynomial of degree two, not four: along a round circle
    `|l|^2` is affine in its angle's cosine and sine (a circle meets a torus
    in at most four points), so where both stored frames are exactly
    orthonormal `G` is of degree two in each angle as a quadric's (every
    fixture) and S9d.4b.2a's discriminant applies; degree four in each
    angle, points of degree eight, arise only where a stored frame is not
    exactly orthonormal (axes rounded by a turn), whose discriminant (degree
    112) is out of reach: there the turning points in `u` are enclosed by a
    certified subdivision (boxes clear of `G` or `G_v`, the rest under
    `1e-6` clear of `G_u`, one clear of none at `1e-10` a tangency), lines
    between them seed the traces, and every turning point's box must lie in
    a verified piece over `v`
    (`turned_tori_meet_in_fields_of_degree_eight`); (c) only extreme circles
    at one height are tangent: parallel tori whose top circles, or a top and
    a bottom, lie at one height and cross touch there (`Degenerate`), but
    with their equators at one height (`tori_side`'s first placement) the
    surfaces cross at 43 degrees and the meeting only turns in `u` there: a
    valid common of two solids, 7.74729, as native DRAW's (OCCT 7.9.3: the
    fuse 69.359, the cut 41.6007; `parallel_tori_at_one_height`).
    Amendments: (i) `Curve3::Toric` gains `other_minor` (the other torus's
    minor radius, zero for a quadric), its jets `S^2 - 4 R2^2 P` from the
    other torus's local coordinates' series; (ii) products in `Q(alpha)`
    reduce by `x^j mod p` kept over one denominator (an integer product and
    one reduction per coefficient), a field's sign is tried by a binary64
    enclosure over its generator's isolator first, and a meeting keeps its
    tangents by point: a turned pair's Boolean from 33 s to 15 s, every
    result unchanged; (iii) the fuzz target's object is a whole torus
    against a whole torus tool by the flags' and the heights byte's top bits
    with the flags' bits 5 and 6 clear (no corpus input of 1,435 decodes to
    one: with the flags' bit alone one did, and took 182 s under ASan), the
    tool then in its offset frame and one operation checked by its volume's
    bounds, its history and its motion; even so two tori take 46, 58 and 111
    s under ASan for variants at the median, the third quartile and the
    ninth decile of those that meet (a host at load 12), over the target's
    60 s: the campaign will meet them (the validator's and the mass's
    certified integrals along the meetings on both tori, which S9d.4b.2a's
    rows fix bit for bit); (iv) the DRAW adapter builds every `ptorus` with
    one set of ids, and a Boolean of two (transforms keep ids) is its `a
    solid other than a prism sharing ids`: two tori do not evaluate there
    yet. The corpus and the regressions (1,435 inputs) replay with debug
    assertions without a failure, the slowest 8.0 s, as do 2,626 variants of
    corpus inputs made two tori (their spline byte, flags and heights byte
    set; the median 2.2 s, the slowest 124 s on a loaded host) and, before
    (iii), 324 made two tori with all three operations, 150 of them chosen
    in turned frames (points of degree eight; the slowest 409 s). The fuzz
    target's torus pairs (`TORUS_PAIRS`) are off: under AddressSanitizer
    they took 46 s at the median and 111 s at the ninth decile of
    intersecting variants, over the target's 60, most in the validator's
    and mass's certified integrals along the meetings; since that parallel
    track (below) about 3 s, 22 s and 51 s at the slowest of 44 on this
    host, still too close to the limit at the Linux runners' 2.6 times
    (the ninth decile about 57 s), so they stay off until the degree-eight
    arrangement arithmetic is faster too; the kernel's tests and the 2,950
    replayed variants cover them meanwhile. DRAW survey: S9d.4b's, below.
    Campaign at `2efa1ec7` (600 s, a sampled replay): 910 runs to the end
    of both budgets, but one slow unit over the limit, 63 s under
    AddressSanitizer (`fd9294c7`, a whole torus in a turned frame against
    a prism with arcs, 9.8 s without; 30 s since the certified integrals'
    track). Campaign: the boolean campaign at `afa36c7f` (600 s, a sampled
    replay) clean, 1,173 runs, the slowest input 52 s under
    AddressSanitizer (S9a.2's spline stack `d0a3de29`, the host loaded). The torus pairs off, as above.
  * **DRAW survey of S9d.4b (2026-09-30, `UPSTREAM_TESTS.md`).** With
    S9d.4b.1, S9d.4b.2a and S9d.4b.2b, and the adapter's curved primitives
    fixed first (the parallel track below): each `pcone`, `psphere` and
    `ptorus` is built under an operation of its own, and a copy of one
    sharing ids with the other argument is built again from its
    constructor's numbers in its frame, as a prism is extruded again. 36
    cases evaluate newly, none wrong: `ZL2` to `ZL5` of the four
    `bop*_simple` grids (16, S9d.4b.2a's: a `pcylinder` of radius 4 and a
    coaxial torus of radii 4 and 1, the wall cutting the tube in two
    circles) and `ZM1`, `ZM2`, `ZM4` to `ZM6` (20, the adapter's: two
    coaxial `pcone`s, the narrow one on, inside or through the wide one);
    volumes, areas and centres of gravity closed forms' within 2.6e-15
    relative, native DRAW's to its printed digits; the 35 not yet in the
    manifest registered, `bopfuse_simple/ZL2` among them (the tori 8.9 s
    at most on the debug worker, the frusta 0.3). The other 41
    refused for shared ids reach the kernel and are `Degenerate`: a cone's
    apex on the other input's surface 21 (frusta on one cone or of one
    virtual apex, `ZL6` to `ZL9`, `ZM3`, `boptuc_simple/ZN2`), a tangency
    between the inputs 19 (coincident rims, `ZM7` to `ZN1`; bases on one
    plane, `ZN2`), a torus tangent to the other input's surface 1
    (`bopfuse_simple/ZP6`'s three copies of a torus about perpendicular
    axes, their tubes touching). No case reaches S9d.4b.1's segments and
    wedges; none fails or times out. Rust evaluates 983 (947), 975
    registered (940); 596 refused (632). The sentinel
    `bopcommon_simple/ZL6` (two cones sharing ids) is re-purposed as the
    frusta's common apex; no sentinel is added (`ZP6`'s native area is a
    known failure on Linux). The volume audit of the last survey's 947
    cases is unchanged bit for bit. The `gdml_public` tori stay refused by
    both hosts (Rust's first refusal `compound result`, then a `ptorus` on
    a DRAW `plane`). A full contract run holds; the ledger does not
    change.

  * **S9's remaining scope (2026-09-30), recorded before its code.** An
    inventory of the kernel's refusals tagged with S9's sub-steps against
    these notes: every named sub-step has its implemented bullet, but these
    remain, each evidence first, in this order. S9d.2c: a sphere's circles
    of unequal axes against a cylinder (S9d.2a's amendment (a) gave them to
    S9d.2b, which built its loops only) and sphere-cylinder loops in turned
    frames (the bite in `TILT`). S9d.3c: cone-cylinder and cone-cone loops,
    and a cone's loops in turned frames (`ball_tilt`), the rest of
    S9d.3b.2's "every pair where no input's rulings meet the other all
    round". S9d.4c: a sphere's cap or zone circle (a surd radius) against a
    torus, and torus v-segments and wedges against prisms with arcs,
    spheres, cones and tori (S9d's scope, in neither S9d.4b sub-step).
    S9e: curved inputs in any position, the general faces S9b's decisions
    deferred to S9c (a Boolean's result with curved faces, a stack with arc
    walls, an imported curved body as an input: the curved analogue of
    S9b.2). S9f: splines joining S9b to S9d where their pairwise
    intersections are exact (the Order decision), assessed first. By design
    and staying refused: a tangency between the inputs, a cavity beside
    several solids, spline segments along one curve of different forms, an
    arc ending off its circle, a section through a sphere's pole off its
    meridians, a torus's tube circle on the other surface, a result
    touching itself (until non-manifold bodies). Two refusals are dead code
    (`meet.rs`'s cylinders meeting in other curves, needing an exactly
    orthonormal oblique frame; `cones.rs`'s torus guard, which `graph.rs`
    never reaches). Acceptance (U6) needs the full boolean replay on a
    schedule run, which exceeds the startup hour under AddressSanitizer in
    one process: it is sharded across schedule jobs now (the parallel
    track below). The S9a.2 and S9b.1 "Pending: the campaign"
    lines were closed by the later clean boolean campaigns (`326ad26c`,
    `29278966`, `85104dc3`, `f4ea7a2a`).
  * **S9d.2c refined, before its code (2026-09-30).** Why each is refused.
    (1) A cap's or zone's own circles (its rims, its split's great circle)
    lie on bases of unequal lengths wherever its stored axis is not of unit
    length exactly (every turned frame: a rim's `y = n * x` has `|y|^2 =
    |n|^2 |x|^2`), so `c + s (cos x + sin y)` with one surd `s` is no longer
    the circle and `circ_quadric` refuses; a whole sphere's circles keep
    equal axes (two rows of a rational rotation, S9d.2a (c)), and so do an
    exact frame's. (2) On a turned cylinder's stored axes the model's circle
    at height `w`, `o + r (cos x + sin y) + w n`, is an ellipse in the world,
    and the sphere's function on it, `F(u, w) = |d + w n + r (cos x + sin
    y)|^2 - R^2` (`d = o - c`), is of degree two in `(cos u, sin u)`, not
    linear as in an exact frame: `u(w)` has no closed form (a quartic in the
    half-angle tangent at each height), which `RiseCrv::at`, its branch and
    `loops`' exact check (the quartic `dw` in `w`) assume. Decisions. (1)
    Such a circle is `c + dx x + dy y` with `xx dx^2 + yy dy^2 = r2`
    (`Circ`'s places are already `(dx, dy)`, and its order and angles take
    unequal axes: plane sections always had them): in coordinates `(s, t)`
    turned from `(dx, dy)` by a rational rotation (the identity, then
    `(3, 4, 5)`, `(5, 12, 13)`, `(8, 15, 17)`), both it and the other
    quadric's function are quadratics in `t`; their resultant is a quartic
    in `s` whose simple real roots give each meeting's `t` as the common
    root `(a2 b0 - a0 b2) / (a1 b2 - a2 b1)`, a point of `Q(alpha)`. A
    repeated root or a vanishing denominator (two meetings of one `s`, or a
    tangency) tries the next rotation; failing all four is a tangency
    (`Degenerate`), a resultant vanishing identically the circle on the
    surface (as `circle_quadric`). A turned frame's rounding never leaves an
    exact tangency, so a crossing within the resolution of one (a rim
    tangent to a cylinder but for the axes' rounding) is `Degenerate` too:
    the meetings with the quadric offset by `+-2 rho res` (`rho` its radius
    at the circle's centre) are counted, and a count that differs from the
    quadric's own means an extremum of its function along the circle within
    that band. Cones take the same function (their other quadric's radius
    term), so a turned cap's circle against a cone is met the same way. (2)
    In a turned frame a loop's pieces over the angle stay `Curve3::Meet`'s
    (the ruling's quadratic is exact in any affine frame), and a piece over
    the height is the one root of `F(., w)` in its branch's half-turn (the
    side `plus` of the rational direction `(alpha, beta) = 2 (x . d, y .
    d)` in the cylinder's local coordinates, S9d.2b's branch test), for `w`
    over its range: its point at a rational height the quartic's root there
    (`Q(alpha)`, degree four). Each such piece is verified exactly: (a) no
    root of `F` on the half-turn's boundary directions over the range, `E(w)
    = rho2 P(w)^2 - l(w)^2` (`F = P + l / sqrt(rho2)` there, a quartic in
    `w`; in an exact frame `-rho2 dw`); (b) no double root of `F(., w)`
    over the range, the discriminant of its quartic in `t` (degree twelve in
    `w`); (c) one root in the half-turn at a rational height inside; then
    the root is analytic in `w` (the implicit function theorem). The height
    graph's turning points (the ellipse tangent to the sphere) are that
    discriminant's real roots: in an exact frame it is `-4 disc(Q_w)
    |Q_w(i)|^4` (`F (1 + t^2)^2 = (1 + t^2) Q_w(t)`, `Q_w(i) = 2 r (alpha +
    i beta)` constant and nonzero), so a turned frame's perturbation keeps
    the pair near `+-i` apart and every real root a real tangency. A point's
    height is its local `w` (the inverse frame's third row, the axis in an
    exact frame), its branch and the cylinder's membership by its local
    `(u, v)`: every exact-frame result unchanged. The topology's
    `Curve3::Rise` keeps its closed form: on the stored axes its binary64
    reading lies on the model's cylinder and within rounding of the sphere
    (as `Meet`'s other frame's axes are read), an edge's range clear of
    turning points by its switches, so no curve type, validator, mass or
    exchange change. A cone's loops in a turned frame stay refused as
    S9d.3c's (`OutOfDomain`). Evidence first: an extension of S9d.2's
    reference taking a cap's end planes as the kernel reads them (through
    `o + h n` normal to the stored axis, the same plane in an exact frame)
    and slicing along a turned cap's axis (the prism cut obliquely), a
    closed form where one exists and a two-way check (a cap and its
    complement against the whole sphere sliced along the prism's axis),
    fixtures of both classes with a declared degenerate (a turned rim
    tangent to a cylinder within rounding), and a native capture before
    `solid/boolean/curved/spheres_turned.rs` exists.
  * **S9d.2c evidence (2026-09-30).** `spheres_boolean_reference.py` now
    reads a cap's end planes as the kernel does (`AxisSphere`: through `o +
    h n` normal to the stored axis; the same tuples in an exact frame, so
    S9d.2's and S9d.3b's references write their fixtures unchanged) and,
    where they are not normal to the prism's axis, slices along the cap's
    axis with the prism cut obliquely: the prism's caps classified in their
    own `(u, v)` against the ball's ellipse and zone half-planes (a flat
    wall's `ball_region`), a cylindrical wall's zone bound varying along
    its generatrices (events where it meets the sphere, the rim crossing
    the wall, a quartic, or the prism's ends, a quadratic), and components
    of one key (two crescents of a circle less an oblique ellipse)
    numbered by their centres along a fixed direction.
    `generate_spheres_turned_boolean_fixtures.py --check`: 18 fixtures (17
    solid, 1 degenerate): hemispheres in `TILT`, `LEAN` and `TILTX` against
    a coaxial pipe (the rim and the split crossing its wall between its
    rings), the bite (the equator across its loop, all three operations)
    and the rod (across its rings); loops in `LEAN` (the bite, all three),
    `TILT` (its fuse and common; S9d.2's fixture is the cut), `TILTX` (a
    long loop round a thick cylinder) and `R125` (turned about the axis); a
    hemisphere in `TILT` against the bite in `LEAN`; and a hemisphere's rim
    tangent to a cylinder within rounding declared `degenerate`. Caps are
    hemispheres only (a height of exactly 0, whatever the platform's
    `sin`). Closed forms (a whole sphere against a turned cylinder by
    S9d.2's lens along the axis; the hemisphere against the coaxial pipe by
    circular segments along the pipe's axis, its faces by central symmetry
    and its disc inside the pipe a circle and an ellipse) within 8.1e-17, a
    cap and its complement against the whole sphere sliced along the
    prism's axis within 1.7e-41, inclusion and exclusion 1.7e-41, the area
    identity 1.4e-40, every face's classes 9.2e-41, a second direction
    (whole spheres) 1.4e-40, Monte Carlo 2.2 standard errors, no near
    coincidence but the declared pair's (its rim's crossings 4.4e-22 of
    the radius squared from tangency), every other cap circle at least 0.16
    from tangency; Python 3.9 and 3.12 write the same files. The capture
    `occt-boolean-spheres-turned-preimplementation`
    (`compare_spheres_turned_boolean.py`; the kernel `unsupported` on all
    18: 8 turned circles of unequal axes, 10 turned loops): every result
    valid with the reference's solids, 1 match, 17 reviewed (BRepGProp's
    default integration on approximated quartics, up to 1.4e-5; adaptively
    converged within 1.5e-8, `dome_lean_bite_common` 4.0e-8 of a small
    result, 1.0e-8 absolute against its area times its edges' tolerance),
    five solids' counts change when unified. S9d.2c's kernel next.
  * **S9d.2c implemented** (`solid/boolean/curved/spheres_turned.rs`,
    `spheres.rs`): a turned cap's circles of unequal axes against a
    cylinder (`circ_ellipse`: the resultant in rotated coordinates, points
    of degree four, a crossing within the resolution of tangency
    `Degenerate`) and a sphere's loops with a turned cylinder (`Height`:
    `F`'s terms, its discriminant's real roots the height graph's turning
    points, the boundary quartic `E`, each piece over the height verified
    and its points at rational heights the quartic's roots on its branch),
    as the refined decisions describe; `RiseCrv` reads heights, branches and
    the cylinder's membership in the carrier's local coordinates (every
    exact-frame result unchanged). All 18 fixtures as the reference (17
    within the kernel's enclosures, `rim_tangent` refused as a crossing
    within the resolution of tangency), every history checked, results
    deterministic and moved rigidly; `Curve3::Rise`'s closed form on the
    turned frames' stored axes lies within `1e-12` of both surfaces
    (`tests/spheres_turned_booleans.rs`, 39 s in the dev profile at
    `opt-level` 2 as CI runs it; the module's own tests check the
    resultant against S9d.2a's surds on a round circle and the
    discriminant's real roots against `dw`'s in an exact frame).
    `compare_spheres_turned_boolean.py` 0 matches and 18 reviewed (the
    measures, and fifteen results' edges and
    vertices where the kernel splits its loops at its switches and OCCT its
    intersection curves at their own points; the faces OCCT's unified
    ones); `compare_spheres_boolean.py` 12 matches and 21 reviewed as
    before, `bite_tilt_cut` now within the reference (its review adds its
    counts); every other comparison unchanged (`ball_tilt_common` still
    `OutOfDomain`, now `(S9d.3c)`). Amendments and corrections: (a) a turned
    cap's circle against a cone stays refused as S9d.3c's (the decisions
    said cones would be met the same way, but no evidence covers them); (b)
    the evidence's `dome_tiltx_rings` does not cut its rings: the rod's
    rings lie on either side of `TILTX`'s equator and the split's great
    circle crosses the rod (the generator's description corrected; the case
    is still of the class); (c) the loop builder's guard for a sphere
    centred on the carrier's axis is a `ComputationLimit` (unreachable: such
    a pair meets in rings), not S9d.2b's refusal. The `boolean` fuzz target
    decodes both configurations already (a turned cap tool against the
    object's arcs, a sphere tool in the tilted frame shared with the
    object's arcs): no decode change. The corpus (1,421 inputs) and the 14
    regressions replay with debug assertions without a failure (the slowest
    7.9 s), and 2,842 variants rewritten into the two configurations (the
    object a stadium or a square with a round hole, the tool a cap or zone
    in a turned frame or a sphere, cap or zone in the tilted frame shared
    with the object) likewise (the slowest 22.8 s): 2,475 of them evaluate
    at least one operation (7,417 of their operations), the others refused
    as documented (a turned cap's plane within the resolution of the
    object's cap, tangencies, a vertex on the other's face), none as
    S9d.2b's `OutOfDomain`. Campaign: the boolean campaign at `afa36c7f` (600 s, a sampled
    replay) clean, 1,173 runs, the slowest input 52 s under
    AddressSanitizer (S9a.2's spline stack `d0a3de29`, the host loaded). DRAW survey: that of
    S9d.2c, S9d.3c and S9d.4c, below (no upstream case's status changes).
  * **S9d.3c refined, before its code (2026-09-30).** Why each is refused.
    (1) `cones.rs`'s `cone_pair` takes a ruled carrier (a cylinder before a
    cone) whose ruling's leading coefficient `A` has no real root and whose
    discriminant `D = B^2 - A C` is positive all round (two rings over its
    angle) or negative all round (apart), a sphere as the other quadric
    (S9d.2b's loops), or a carrier whose `A` has simple real roots and whose
    `D` is nowhere negative (S9d.3b.2's open pieces, split where a branch
    runs off to infinity); every other pair of two ruled faces, a cone
    against a cylinder or a cone whose `D` changes sign over both carriers'
    angles, is `OutOfDomain("a cone meeting a curved face in a loop
    (S9d.3b.2)")`: its curve turns over both angles, so no graph over one
    angle covers it. S9c.2b.1's `turned::crossing` covers such loops for two
    cylinders (graphs over both angles, switched between turning points of
    different kinds), but reads `A` as a cylinder's constant and the other
    quadric as a cylinder. (2) `plane_of`'s refusal: a carrier whose `A`
    vanishes identically is taken for two cones whose quadrics differ by an
    affine function (their asymptotic cones one cone, their quadratic parts
    proportional); a cylinder exactly along a cone's ruling makes the
    cylinder's `A` zero as well (and the cone's `A`, nonnegative, a double
    root), its parts not proportional, and falls through as a loop. The
    curve is then a graph over the cylinder's angle with one finite root per
    ruling, running to infinity where `B` vanishes. (3) `spheres::loops`
    refuses a turned cone (`ball_tilt`, since S9d.2c with an S9d.3c
    message): S9d.2c's height graph `F(u, w)` and its checks were written
    for a cylinder's constant radius. (4) `meet.rs` refuses a turned cap's
    circle (of unequal axes) against a cone for want of evidence (S9d.2c's
    amendment (a)). Decisions. (1) Two ruled faces whose `D` changes sign
    over both carriers meet in loops met as S9c.2b.1 meets two turned
    cylinders', generalised to ruled carriers (any affine frames, exact or
    turned: the ruling's quadratic is exact on the stored axes): carrier 0
    the cylinder (its `A` constant), else operand 0's cone; each interval of
    `D_0 > 0` (in a chart based at a negative point) one component, its
    plus branch from one end to the other and its minus branch back; the
    other carrier's turning points (`D_1`'s real roots, where its ruling
    touches the first quadric, `w = -B_1 / A_1`) placed along it by binary64
    views (the chart's `t` and the branch); switches at rational angles of
    carrier 0 between adjacent turning points of different kinds (points
    with one surd, `sqrt(D_0)`); the pieces between switches graphs over
    carrier 0's angle about carrier 1's turning points and over carrier 1's
    angle about carrier 0's, each `Curve3::Meet` with that carrier (a
    cone's `half_angle`, the other's `other_half_angle`): no new curve type,
    validator, mass or exchange change. Each piece is verified exactly: no
    root of its carrier's `D` strictly inside its range (Sturm counts at its
    ends, surds for carrier 1's ranges), both switches and an exact interior
    point (on carrier 0's branch at a rational angle) on its branch; the
    switches are ordered exactly along the component (rational angles on
    known branches), so the verified pieces are its arcs between
    consecutive switches and cover it once. Degrees: `D` a quadratic form
    over a cylinder's angle (a quartic in the chart's `t`), a quartic form
    over a cone's (an octic); the curve a quartic. A carrier with `D`
    positive all round keeps S9d.3b.1's rings. (2) Where the rulings reach
    the other's asymptotic directions (`A`'s simple real roots: two cones
    whose direction cones cross; a cone and a cylinder never, but along a
    ruling) and `D` changes sign over both, a component runs through
    infinity: the directions where `A_0` vanishes are cuts on the branch
    that runs off there (`-sign B_0`; `D_0 = B_0^2 > 0` at them), the
    finite branch's point there (`w = -C_0 / 2 B_0`, in `Q(alpha)` of degree
    four) a switch, as S9d.3b.2's open pieces, and carrier 1's `A_1` roots
    alike; a piece may end at a cut (its range ending at the algebraic
    direction, its branch the one running off there), verified as in (1)
    with no root of `A` inside either. The pieces are then the arcs between
    consecutive switches and cuts, the first point at infinity along each
    the one its verified piece reaches. (3) A cylinder exactly along a
    cone's ruling (the cylinder's `A` zero, the cone's a double root), and
    any repeated real root of `A`, is `Degenerate`: the limit of S9d.3b.1's
    "a ruling within rounding of the other's direction", whose branch at
    infinity a turned frame's rounding puts on either side. (4) A turned
    cone's loops with a sphere: S9d.2c's height graph with the cone's radius
    `rho(w) = r + k w`: `F(u, w) = |d + w n + rho(w) (cos x + sin y)|^2 -
    R^2` has coefficients of degree two in `w` (the circle's terms
    `rho(w)^2 x . x` and `2 rho(w) x . (d + w n)`), the half-turns' boundary
    quartic `E(w) = rho2 P(w)^2 - l(w)^2` (`P` and `l` quadratics now), the
    discriminant in `w` of degree at most twelve (its real roots the height
    graph's turning points), and a piece's point at a rational height the one
    root on its branch's half-turn (`Q(alpha)`, degree four) where `rho(w) >
    0`; `Curve3::Rise` with the cone's `half_angle` keeps its closed form as
    the binary64 reading. (5) A turned cap's circle against a cone: S9d.2c's
    resultant with the cone's radius term (as its decisions said), the band
    of its tangency test `2 rho res` with `rho` a bound of the radius term
    over the circle (for a cylinder or a sphere its constant, as before).
    (6) `Degenerate`: a repeated real root of `D` (a tangency), an extremum
    of `D` within the resolution of zero (S9c.2b.1's `near_node`, `A`'s
    sampled least size its scale), two turning points, or a turning point
    and a cut, within rounding along a component, a cone's apex on the other
    quadric (S9d.3b.1), and a turned cap's circle within the resolution of
    tangency to a cone (S9d.2c's rule). The kernel in
    `solid/boolean/curved/cones_loops.rs` (the ruled loops) with
    `spheres_turned.rs`'s height graph taking a slope. Evidence first: the
    cones' reference (`cones_boolean_reference.py`, S9d.2c's reading of a
    turned cap's end planes) on fixtures of cone-cylinder and cone-cone
    loops in exact and turned frames (a rod grazing a frustum's wall, a rod
    across a cone's tip region, a rod in `LEAN`, a thin cone across a
    frustum's wall, a wide cone whose direction cone crosses the frustum's),
    turned cones' loops against spheres, and turned caps against cones, with
    a cylinder along a cone's ruling exactly and a turned rim tangent to a
    cone within rounding declared `degenerate`; closed forms where they
    exist (a rod across a cone's axis by S9d.3b's strip, a whole sphere
    against a turned cone by its lens along the cone's ideal axis, a
    hemisphere against a coaxial cone by circular segments along the axis),
    two-way checks (an input split in two along its axis against the whole,
    a cap and its complement against the whole sphere, a second slicing
    direction), and a native capture before `cones_loops.rs` exists.
  * **S9d.3c evidence (2026-09-30).** `cones_boolean_reference.py` now
    takes a turned cap (S9d.2c's `AxisSphere` reading of its end planes, its
    circle's edge and its disc spanned by vectors exactly in its plane; the
    frame's axes in an exact frame, so S9d.3b's rows are unchanged, its
    `--check` passing). `generate_cones_loops_boolean_fixtures.py --check`:
    31 fixtures (29 solid, 2 degenerate; 9 in exact frames): cone-cylinder
    loops (`rod_graze`, all three, and `rod_bitten`: a rod across a
    frustum's axis grazing its wall; `tip_graze` beside a cone's apex;
    `rod_lean`; `tilt_rod`, the frustum in `TILT`), cone-cone loops
    (`cones_graze`, all three; `cones_turned`, `TILT` against `LEAN`;
    `cones_asymptotic`, a cone in `LEAN` whose direction cone crosses the
    frustum's, the curve through infinity), turned cones' loops with spheres
    (`ball_tilt`, all three, S9d.3b's pair; `ball_lean`; `ball_r125`),
    turned hemispheres against cones (`dome_tilt_cone`, coaxial;
    `dome_lean_frustum`; `dome_cone_turned`, both turned), and declared
    `degenerate` a cylinder exactly along a cone's ruling (`rod_ruling_exact`,
    the cone's slope `fl(0.6) / fl(0.8)` of `LEAN`'s stored axis, bit
    checked in the frames' file) and a turned rim tangent to a cone within
    rounding (`rim_tangent_cone`). Closed forms (a rod across the axis by
    S9d.3b's strip, a whole sphere against a turned cone by the lens, a
    hemisphere against a coaxial cone by circular segments along the axis)
    within 2.2e-40 in exact frames and 1.3e-16 in turned ones; two-way
    checks: an input cut in two along its axis against the whole within
    5.3e-42 (six pairs), a hemisphere and its complement against the whole
    sphere within 1.4e-41; inclusion and exclusion 1.2e-41, the area
    identity and every face's classes 6.9e-41, the cone's wall two ways
    2.7e-41, a second direction 1.7e-41 (twelve pairs), Monte Carlo 2.5
    standard errors, no near coincidence but the declared pairs', every
    other edge and vertex at least 0.0156 from tangency with or incidence on
    the other's surfaces, every plane at least 0.05 of the radius from
    tangency to the other's spheres and parallel cylinders (a check added
    when the kernel found two first fixtures' planes tangent to a sphere,
    `Degenerate` by S9d.1's rule: `ball_r125`'s sphere on the frustum's base
    inside its disc, `dome_lean_frustum`'s top plane on the sphere off both
    faces; both moved clear and captured again before the kernel's commit);
    Python 3.9 and 3.12 write the same files; the
    generator's check a CI group of its own (`cones-loops`, 25 CPU minutes
    locally). The capture `occt-boolean-cones-loops-preimplementation`
    (`compare_cones_loops_boolean.py`; the kernel `unsupported` on all 31:
    17 as S9d.3b.2's loops, 8 as a turned cone's loop, 6 as a turned cap's
    circle against a cone): every result valid with the reference's solids,
    7 match (the turned hemispheres), 24 reviewed (BRepGProp's default
    integration on approximated loops, up to 5.8e-6, 4.5e-4 on the declared
    ruling; adaptively converged within 1.7e-8), five solids' counts change
    when unified. S9d.3c's kernel next.
  * **S9d.3c implemented** (`solid/boolean/curved/cones_loops.rs`,
    `cones.rs`, `spheres.rs`, `spheres_turned.rs`, `meet.rs`): cone-cylinder
    and cone-cone loops over both carriers' angles in any frames (a
    component through infinity cut at the rulings' asymptotic directions, its
    finite branch switched at `-C / 2B`), a cylinder exactly along a cone's
    ruling `Degenerate`, a turned cone's loops with a sphere by S9d.2c's
    height graph with the cone's radius, a turned cap's circle against a
    cone by S9d.2c's resultant, as the refined decisions describe. All 31
    fixtures as the reference (29 within the kernel's enclosures, each at
    most `1e-9` wide; `rod_ruling_exact` refused as "a cylinder along a
    cone's ruling", `rim_tangent_cone` as S9d.2c's crossing within the
    resolution of tangency), every history checked, results deterministic
    and moved rigidly; the loops' `Curve3::Meet` pieces run over both
    carriers and lie within `1e-9` of both surfaces, the turned cones'
    `Curve3::Rise` pieces within `1e-12`, the turned rims' crossings on the
    rim's plane (`tests/cones_loops_booleans.rs`, 24 s in the dev profile at
    `opt-level` 2; the module's tests check Sturm's count with a root at a
    range's end and the boundaries' order, `spheres_turned.rs`'s a cone's
    turning points against `dw`'s with the apex's left out).
    `compare_cones_loops_boolean.py` 5 matches and 26 reviewed (the
    capture's 24 measures, and the kernel's edges and vertices where it
    splits its loops at its switches: `dome_cone_turned`'s two now reviewed
    for their counts alone); `compare_cones_boolean.py` 21 matches and 19
    reviewed as before, `ball_tilt_common` now within the reference (its
    review adds its counts); every other comparison unchanged.
    Amendments and corrections: (a) the evidence had two pairs with a plane
    tangent to a sphere (`ball_r125`'s sphere on the frustum's base inside
    its disc, `dome_lean_frustum`'s frustum's top plane on the sphere off
    both faces), refused by S9d.1's rule wherever the touch lies: the
    generator now checks every plane against the other's spheres and
    parallel cylinders, both pairs moved clear, and the capture was taken
    again before the kernel's commit (its own commits, `c3eeaa31` and
    `555a5974`); (b) an asymptotic direction's binary64 place along a
    component is its generator's enclosure (a first build took the
    isolator's middle, far off before refinement: `cones_asymptotic`'s cut
    misplaced, a piece verified across an `A` root and refused); (c) a point
    at a cone carrier's apex height is on none of its `Curve3::Meet` pieces
    (`MeetCrv::on`: the replays reached `place` with a zone's pole there, a
    panic, 25 variants; a regression); (d) the height graph's turning points
    on a cone's far nappe lie on the other side of `phi` (the side by `g /
    rho`, the branch by the sphere's gradient along the ruling rather than
    the axis): a sphere meeting both nappes had four events at one place
    (`two turning points within rounding`, or no event of the other kind),
    S9d.3b.2's exact-frame loops alike (a regression); (e) S9d.3b's
    `ball_tilt_common` evaluates (`tests/cones_booleans.rs`'s `LATER` empty),
    and S9d.2c's test of the S9d.3c refusals is gone; (f) the near-node test
    of a carrier whose `A` has real roots takes `A`'s largest sampled size.
    The `boolean` fuzz target decodes the object a sphere against a cone or
    sphere tool and a cone against a sphere, cap or zone tool (the flags'
    bits 5 and 6, as against a whole torus), a cone against a cone tool off
    (`CONE_PAIRS`: two cones' certified integrals, the result's mass and
    validation, take up to 80 s with debug assertions alone, 14 of 711
    variants above 20 s, past the target's 60 s under the sanitizer). The
    corpus (1,421 inputs) and the 16 regressions replay with debug
    assertions without a failure (the slowest 13.1 s), and 2,842 variants
    rewritten into the new configurations (the object's arcs against a
    turned or sideways cone tool, a cone against a cone with `CONE_PAIRS`
    on, a sphere against a cone, a cone against a turned cap or zone)
    likewise (the slowest 79.5 s, a cone pair; a sphere and a cone 25.7 s):
    2,349 of them evaluate at least one operation (7,045 of their
    operations), the others refused as documented (a spline profile in any
    position, S9c; two parallel faces within the resolution of one plane, a
    plane through a cone's apex, tangencies, thin pieces, a cone's apex on
    the other's surface, a ruling within rounding of the other's direction,
    a near node), none as S9d.3b.2's or S9d.3c's `OutOfDomain`. Campaign: the
    boolean campaign at `428349e8` (600 s, a sampled replay) clean, 928
    runs, the slowest input 17 s under AddressSanitizer (cone pairs off,
    `CONE_PAIRS`). DRAW survey: that of S9d.2c, S9d.3c and S9d.4c, below
    (no upstream case's status changes).
  * **S9d.4c refined, before its code (2026-09-30).** Why each is refused.
    (1) `torus_curved::circ_torus` meets a sphere's own circle with a torus
    only on a basis of equal lengths whose scale to the radius is rational
    (`c + s (cos x + sin y)`, `s` rational: a whole sphere's great circle in
    any frame, a hemisphere's rim in an exact one), as a conic
    (`conic_torus`: the torus's quartic along it, an octic in the
    half-angle tangent). A cap's or zone's rims have the radius `sqrt(R^2 -
    h^2)`, a surd wherever `R^2 - h^2` is no square (every rim but a
    hemisphere's), and in a turned frame every circle of a cap or zone lies
    on a basis of unequal lengths (S9d.2c); both are refused with
    `OutOfDomain("a sphere's circle of a surd radius against a torus
    (S9d.4b)")`. With a surd scale the octic's coefficients lie in
    `Q(sqrt(sigma))`, and its norm, of degree sixteen, holds every crossing
    twice (the parameter's and its antipode's). (2) `curved::attempt`
    refuses a torus v-segment or wedge against anything but a polyhedral
    prism (`OutOfDomain("a torus segment or wedge against a curved face
    (S9d.4b)")`): S9d.4b.1 meets a part's rims (S9d.4a's `TorusSec` rings:
    a segment's end plane's ring over `u`, of radius `R +- sqrt(r^2 -
    z^2)`; a wedge's end half-plane's ring over `v`, about `R (c, s) /
    sqrt(c^2 + s^2)` on its rounded end direction) with planes only
    (`rim_plane`, the two planes' line against the torus), and
    `meet::edge_surface` holds no case of a rim against a cylinder, sphere,
    cone or torus (`unreachable`); S9d.4b.2 built the wall's meeting with a
    quadric or a torus for whole tori. Decisions. (1) A sphere's circle `c +
    dx x + dy y`, `xx dx^2 + 2 xy dx dy + yy dy^2 = r2` (`r2` rational, the
    basis's lengths any: a rim in an exact or a turned frame, a split's great
    circle), against a torus: the torus's function on the circle's plane is
    a quartic in `(dx, dy)` (its local coordinates affine in them); in
    coordinates `(s, t)` turned by a rational rotation (S9d.2c's, then
    further Pythagorean ones) the circle is a quadratic in `t` of constant
    leading coefficient (positive: the basis's Gram form turned), the
    quartic reduced modulo it is linear in `t`, `r1(s) t + r0(s)`, and their
    resultant `a2 r0^2 - a1 r0 r1 + a0 r1^2` is a polynomial of degree at
    most eight in `s`, whose simple real roots give one crossing each, `t =
    -r0 / r1`: a point of `Q(alpha)` of degree at most eight, as
    `conic_torus`'s. A repeated root, or `r1` vanishing at a root (two
    crossings of one `s`, or a tangency), tries the next rotation; failing
    all is a tangency (`Degenerate`); a resultant vanishing identically is
    the circle on the torus (`Along`, as `conic_torus`'s). A crossing within
    the resolution of a tangency is `Degenerate`, as S9d.2c's: the counts
    with the torus's function offset by `+-delta` must be its own, `delta = 8
    R r (R + r) res` (on the torus its gradient in local coordinates is `8 R
    r rho`, `rho <= R + r`). A round circle with a rational scale keeps
    `conic_torus` (every S9d.4b.2 row unchanged). The cap's discs meet the
    torus in S9d.4a's spiric sections (a plane tangent to the torus
    `Degenerate`), its wall in S9d.4b.2a's meeting with the whole sphere,
    kept where the cap's face holds it. (2) A part's wall meets a quadric or
    a torus in S9d.4b.2's meeting with the whole torus's surface, traced and
    verified over the whole torus, its pieces kept where the part's patches
    hold them (S9d.4b.1's face membership, `part_in_face`); its end discs
    meet quadrics in S9c's and S9d's plane sections and another torus in
    S9d.4a's spiric sections; what is new is a rim against a curved surface.
    Along a rim one angle is fixed at a quadratic surd: a segment's `v` at
    `(s sqrt(q2), z) / r` (`q2 = r^2 - z^2`, `s` the sign of the latitude's
    cosine), a wedge's end `u` at `(c, s) / sqrt(n2)` (`n2 = c^2 + s^2` of
    its rounded direction; the start's `u = (1, 0)` rational). The other
    surface's function on the torus's angles, S9d.4b.2's `G(u, v)`, is there
    a form in the free angle with coefficients in `Q(sqrt(d))`, `A + sqrt(d)
    B`; in a chart whose antipode is no crossing its polynomial `p = A +
    sqrt(d) B` has its roots among the norm's, `A^2 - d B^2`, rational and of
    degree eight (`G` of degree two in each angle: a quartic over
    `Q(sqrt(d))`; sixteen for two tori in frames not exactly orthonormal,
    S9d.4b.2b's degree four in each angle; two exact tori meet in degree
    four, as a quadric does). The rim's crossings are the norm's real roots
    at which `p` vanishes, decided exactly in `Q(alpha)(sqrt(d))` (the
    conjugate's are the other ring's, or the opposite half-plane's
    circle's), each point in that tower (the fixed angle's surd with the
    free angle's field); a root of `p` where `p'` vanishes too is a
    tangency, `Degenerate`. A rim of rational radius (the halves' circles of
    radius `R` at `w = +-r`, stored as conics) keeps S9d.4b.2's conic
    crossings; a `d` that is a square is rational. No curve type,
    validator, mass or exchange change: the rims keep `TorusSec`, the
    meetings `Curve3::Toric`, the cap's circles `Circ`. (3) `Degenerate`: a
    singular point of a torus's meeting with the other surface anywhere on
    the whole torus (S9d.4b.2's regularity is certified over the whole torus
    before its pieces are cut to the part; a tangency off a part's wall is
    refused, not taken as no contact as S9d.4b.1 (d) takes a plane's: a
    quadric or torus tangent off the wall is no part's construction, as the
    halves' end planes are), a rim or a sphere's circle tangent to the other
    surface, a sphere's circle within the resolution of tangency to a torus,
    a plane tangent to a torus (S9d.4a's; on a part's own wall S9d.4b.1's),
    and every earlier rule (a cone's apex on the other's surface; a tube's
    circle on the other surface by design). A half's end plane is tangent
    to its torus along its rim: a surface crossing the rim meets the disc
    and the wall in curves tangent there, on opposite sides of the crossing
    (one curve through a smooth edge), as S9d.4b.1's planes do. The kernel
    in `solid/boolean/curved/torus_parts.rs` (a sphere's circle's resultant
    and the rims' crossings), `mod.rs`'s refusal gone. Evidence first: an
    independent reference (`torus_parts_boolean_reference.py`: S9d.4b.2's,
    every face of both inputs swept by two families of curves, with
    S9d.4b.1's part model and a cap's or zone's, their end discs by chords
    and their walls over their ranges; solids by the torus's normal slices),
    fixtures of sphere caps and zones against tori (coaxial ones meeting in
    circles, caps and zones across the tube with rims of surd radius, turned
    frames) and of the outer and inner halves, a band and wedges against a
    cylinder, a sphere, a cone and a torus (coaxial ones, across the rims
    and end discs, turned frames), a cap's rim within the resolution of
    tangency to a torus and a part's rim tangent to a cylinder declared
    `degenerate`; closed forms where the other's sections are discs about
    the torus's axis (coaxial pairs; a wedge's by its turn and its end
    discs), two-way checks (a segment and its complement against the whole
    torus, a cap and its complement against the whole sphere, both families
    of every face), and a native capture before
    `solid/boolean/curved/torus_parts.rs` exists.
  * **S9d.4c evidence (2026-09-30).** `torus_parts_boolean_reference.py`
    takes S9d.4b.2's reference (every face of both inputs swept by two
    families of curves against the other's surfaces and membership, the
    divergence theorem on the classified pieces) with two more inputs:
    S9d.4b.1's part model (a segment's end planes at the stored heights,
    its wall's tube circles over the arc `[v_lo, v_hi]` and parallels over
    that range, oriented by its `sigma`, its end discs from the axis to the
    arc's ends; a wedge's half-planes on the rounded end direction, its wall
    over the turn and the tube's discs in its half-planes) and a cap or zone
    as S9d.2c's `AxisSphere` reads it (end planes through `o + h n` normal
    to the stored axis, its wall's meridians and parallels about the unit
    axis between the ends' latitudes, its discs by chords); solids by the
    part's (else the torus's) normal slices, a thin region sheared between
    two slices joined through a chain of intervals at levels between them
    (bisected: overlaps alone, as S9d.4b.2's sweep joins, counted the tip of
    a half's ring outside a coaxial dome apart). Heights are `r
    sin(latitude)` rounded once (`sin_rn`, asserted equal to the host's
    `sin` for every latitude used). `generate_torus_parts_boolean_fixtures.py
    --check`: 53 cases (25 pairs; 7 fuses, 23 cuts, 23 commons; 51 solid, 2
    degenerate; 45 in exact frames) on the torus of radii 5/2 and 1: a zone
    and a cap of a coaxial sphere (circles), a cap on the tube's top with a
    rim of surd radius across the tube (all three), a zone about `x` across
    the outer side, a hemisphere in `LEAN` (its rim on unequal axes) and a
    zone against the torus in `TILT`; the outer and inner halves against
    coaxial pipes, a coaxial cone, sphere and lower hemisphere, a sphere
    about a point of the outer half's top rim and a small torus ringing the
    tube's outer side; the band `0.5..2.25` against a sphere, a frustum and
    a ringing torus across its rims of surd radius, and in `TILT` against a
    coaxial pipe; a quarter wedge against a pipe (all three) and a small
    torus of a parallel axis across its end disc, a half wedge against a
    sphere across its start disc and rim, a three-quarter wedge against a
    frustum across its end disc, a quarter wedge in `LEAN` against a sphere
    across its start rim; `degenerate` a cap's rim within the resolution of
    tangency to the outer equator and a quarter wedge's start rim touching
    a sphere. Checks: closed forms of the 9 coaxial pairs (one quadrature
    along the axis of both inputs' radial intervals, the walls by their own
    elements, a wedge by its turn and its meridian discs) within 6.9e-40 in
    exact frames and 3.8e-17 turned; two-way checks, a half and the other
    half (weighted by their orientation) against the tool as the whole torus
    against it (S9d.4b.2's reference) and a cap or zone and the rest of its
    sphere against the torus as the whole sphere, within 4.3e-42 (five
    pairs; a band's rest of the tube bounds its region with a boundary
    crossing itself, whose parity set is no signed complement); every
    operation two ways 1.5e-31, both inputs from their faces 8.3e-35,
    inclusion and exclusion 1.8e-35, the area identity 8.3e-35, Monte Carlo
    2.5 standard errors; no near coincidence outside the declared pairs,
    every other pair's surfaces meeting at a sine of at least 0.38, edges
    (the rims too) crossing the other's surfaces at 0.29, and every plane at
    least 0.0059 of the case's size from tangency to the other's spheres and
    tori (a check the kernel's rules need: a plane tangent to a sphere or a
    torus is refused wherever it touches). Python 3.9 and 3.12 write the
    same files; the generator's check a CI group of its own (`torus-parts`,
    50 CPU minutes locally). Fixture corrections before the capture, from
    the reference's and a first kernel's runs: a cylinder tangent to the
    outer half's rim also touched the outer equator (a surface tangency),
    a small torus about `x` across the quarter wedge's end disc lay within
    rounding of its equatorial plane (the end's rounded direction: a near
    node), a frustum across the three-quarter wedge's end disc had its axis
    in the end plane (a plane through its apex), a zone's sphere passed
    through the torus's top circle with a parallel in its tangent plane
    (the reference's quadrature did not converge), a small torus crossing
    the outer half's rims, where its end planes are tangent to the torus,
    left the reference's quadrature nodes off their intervals' structure,
    and events of two crossings symmetric about a wedge's rim's centre fell
    within rounding of each other: all moved clear. The capture
    `occt-boolean-torus-parts-preimplementation`
    (`compare_torus_parts_boolean.py`, keyed on
    `solid/boolean/curved/torus_parts.rs`; the kernel `unsupported` on all
    53, 16 as a sphere's circle of a surd radius against a torus and 37 as
    a torus segment or wedge against a curved face): every solid count the
    reference's, 26 match (the coaxial pairs within 2.4e-14), 27 reviewed:
    25 BRepGProp's default integration on B-spline sections (up to 7.4e-6;
    adaptively converged within 2.3e-8), and `qw_lean_ball`'s two results
    wrong natively (the cut 4.7e-5 off in area, the common 1.6e-4 in volume
    and invalid under `BRepCheck_Analyzer`); fourteen solids' counts change
    when unified. S9d.4c's kernel next.
  * **S9d.4c implemented** (`solid/boolean/curved/torus_parts.rs`,
    `torus_curved.rs`, `meet.rs`, `mod.rs`, `assemble.rs`, `num.rs`): a
    sphere's circle of a surd radius or on unequal axes against a torus by
    the resultant of its quadratic and the torus's quartic in rationally
    rotated coordinates (eight rotations, the band test with `8 R r (R + r)
    res`), and torus v-segments and wedges against prisms with arcs,
    spheres, cones and tori, their rims against a curved surface by the
    norm of the other surface's function at the rim's fixed surd angle (the
    points in `Q(alpha)(sqrt(d))`), `mod.rs`'s refusal gone, as the refined
    decisions describe. All 53 fixtures as the reference (51 within the
    kernel's enclosures, each at most `1e-9` wide; `rim_touch` refused as "a
    sphere's circle within the resolution of tangency to a torus",
    `qw_rim_touch` as "a torus part's rim tangent to the other input's
    surface"), every history checked, results deterministic and moved
    rigidly, the coaxial pairs in circles and the others in `Curve3::Toric`
    pieces, a band, a quarter wedge and a cap turned with their tools
    keeping the reference's volumes (`tests/torus_parts_booleans.rs`, 48 s
    in the dev profile at `opt-level` 2). `compare_torus_parts_boolean.py`
    16 matches and 37 reviewed (the capture's 27, 23 of them with the
    kernel's entity counts, and ten native matches by their counts alone:
    coaxial bands crossing OCCT's seams, sections split at OCCT's points
    against the kernel's exact pieces); every other Boolean comparison
    unchanged in its counts (boolean 45/0, its splines 33/13, polyhedral
    43/2, curved 42/2, procedural 4/24, turned 2/13, capped 0/18, sphere
    30/0, spheres 12/21, cone 25/5, cones 21/19, torus 11/24, torus segment
    15/14, torus curved 15/29, spheres turned 0/18, cones' loops 5/26).
    Amendments and corrections: (a) a rim's fixed angle is read from the
    rim's own point at the free angle's `(1, 0)` (`TorusSec::at`: a
    segment's `v` direction over `r`, a wedge's `u` direction over the
    point's distance from the axis), whichever branch the rim is; (b) a
    fuse's void between both inputs' faces (a band's end disc and inner
    wall under another input's face over its hole) was assembled as a solid
    of its own and failed validation (`shell_orientation`; S9d.4b.1's box
    over a band's hole the same, found by the replays here): a shell of
    both inputs' faces whose certified flux, the shell built alone, is
    turned inward is now a cavity of the result, whose containment the
    validator's rays leave undecided (`ComputationLimit`, as S9d.4b.1 (g));
    (c) a loop through a sphere's pole (a quarter wedge's end half-plane,
    on the rounded direction of its turn, within rounding of a sphere's
    axis) wound by the sum of its pcurves' changes, a tie of half a turn at
    the pole rounded away (`uv_gap`): a loop's winding is now its last
    pcurve's end against its first's start; refusing a section within the
    resolution of a pole but off it was tried first and dropped, as it
    refused S9d.1's `wedges_through_a_spheres_poles` (the DRAW grids'
    `ZI5`), which the winding by ends keeps; (d) surds over algebraic
    fields (`Qd::sign`, `tower_sign`, `mixed_dot_sign`) try a binary64
    enclosure before their exact products (a turned zone against a turned
    quarter wedge about 16 to 12 s; every other comparison's rows unchanged); (e)
    the DRAW bridge's guard of a half turn against a cylinder (S9d.4b)
    evaluates now, moved to its common's volume, `torus_parts_boolean_
    reference.py`'s 0.53575919723822447; (f) the `boolean` fuzz target
    decodes a sphere or cone object against a v-segment or wedge tool and a
    cap or zone object (the heights byte's top bit) against a whole torus or
    a part, in the axis-aligned frames; in the tilted or a turned frame only
    with `TURNED_PARTS`, off: of 3,000 replayed variants the 1,201 caps and
    592 parts in turned frames took 2.1 s and 0.6 s at the median, 7.5 s
    and 4.6 s at the ninth decile and 26 s and 50 s at the slowest (debug
    assertions, no sanitizer), against 2.3 s at the slowest of 457 in the
    axis-aligned frames; parts against prisms with arcs, refused before,
    evaluate in every frame (4.5 s at the slowest of 750). The corpus (1,423
    inputs) and the 17 regressions replay with debug assertions without a
    failure (the slowest 2.9 s), as do the 3,000 variants with the committed
    decode (the slowest 4.6 s; 24 to 36 s under AddressSanitizer for the
    five slowest, a host running the test suite beside it) and with
    `TURNED_PARTS` on (the slowest 50 s), which found (b) and (c) (two kept
    as regressions, `fuzz/regressions/README.md`). No DRAW upstream case's
    status changes before the survey (none reaches a part or a cap against
    a torus). DRAW survey: that of S9d.2c, S9d.3c and S9d.4c, below (no
    upstream case's status changes). Campaign: the boolean
    campaign at `07ea3b3e` (600 s, a sampled replay) clean, 975 runs, the
    slowest input 23 s under AddressSanitizer (turned parts off,
    `TURNED_PARTS`).
  * **DRAW survey of S9d.2c, S9d.3c and S9d.4c (2026-09-30,
    `UPSTREAM_TESTS.md`).** The 1,802 cases of the Boolean group run again
    on both backends after S9d.2c, S9d.3c, S9d.4c and the certified
    integrals' speed-up (a parallel track below): no case evaluates newly,
    none changes its refusal, none fails or times out, and the sentinels
    are refused as before. The group's `psphere`s and `ptorus`es are whole
    solids (no cap, zone or torus part) and its cones meet the `pcylinder`
    in circles or rings, not loops, so none of the three sub-steps shows.
    The speed-up does: `ZK7` and `ZK8` of the four grids (8, S9d.3b.1's
    frustum of radii 6 and 1 across the cylinder, its wide end through
    both caps), right since S9d.3b.1's survey but 23 to 34 seconds alone on
    the debug worker at a load of 5 and up to 120 loaded, too near or past
    the contract's 30, take 8.7 to 12.9 now (the
    worker before the speed-up 18.9 against 8.9, side by side at a load of
    1.7); their measures are `cones_boolean_reference.py`'s within 1.3e-15
    relative (native DRAW's within 7.8e-6), and they are registered.
    Rust evaluates 983, 983 registered (975), 596 refused, none failing
    and none timing out (the other counts as before). The volume audit of
    the 983 cases: native DRAW's values unchanged, Rust's bit for bit on
    854 and within rounding on 129 (at most 2.2e-15 relative in volume,
    7.2e-16 in area, 1.4e-15 in the centres): 118 by the speed-up's
    binary64 series products (the worker built at its parent reproduces
    the last audit bit for bit, at it these values) and 11 centres of
    `ZI4` to `ZI7` (the section through a sphere's poles) by up to
    2.6e-16 with S9d.4c's kernel (the worker at its predecessor gives the
    last audit's); no bug. Against the closed forms and references the
    turned spheres are within 1.4e-15 (7.8e-16 before), the crossed
    cylinders' areas 1.0e-15 (6.7e-16), the crossing frusta 2.8e-15, the
    torus and coaxial frusta 2.2e-15; Rust and native DRAW disagree on the
    same 35 as before, native off in each. `bopfuse_simple/ZP6` stays a
    torus tangent to the other input's surface; the `gdml_public` tori stay
    refused by both hosts (Rust at `compound result`, then a `ptorus` on a
    DRAW `plane`; native at `add` or `wire`). A full contract run holds
    (the slowest Boolean case 9.4 seconds); the ledger does not change.

  * **S9e refined, before its code (2026-09-30).** Why each is refused.
    (1) `polyhedra::build` gives a pair to the curved arrangement only when
    both inputs are constructions with exact models (`curved::applies`:
    prisms, spheres, cones, tori); a Boolean's result
    (`Construction::Polyhedron`, whose inputs it keeps), a stack
    (`Construction::Stack`, its two profiles and heights in the object's
    frame) or a plane's piece goes to S9b.2's stored model, which refuses
    any face but a plane and any edge but a line (`OutOfDomain("a Boolean of
    a solid with curved faces or edges in any position (S9c)")`); a
    polyhedral one against a prism with arcs is refused by S9b.1's
    `prism_model` (`"a Boolean of a prism with arcs and a solid other than
    a prism (S9c)"`), and `model::Prism::new` refuses every other
    construction (`"a Boolean of a solid other than a prism with arcs in
    any position (S9c)"`). (2) What such a body stores is not an exact
    model: its faces keep the inputs' stored surfaces (a wall's cylinder on
    a frame at its circle's rounded centre, not the prism's affine cylinder
    on its stored axes), its vertices and edge curves are rounded once from
    the arrangement's exact points and curves, so no stored vertex lies on
    its faces' stored surfaces exactly, and what the inputs' models decide
    exactly (a profile's line tangent to its arc, a cap shared by both
    inputs, a section through a vertex) does not survive the rounding; a
    model re-derived from the stored surfaces (an edge as its faces'
    meeting, a vertex as their common point) meets a smooth join or a
    shared surface as a tangency within rounding. `Curve3` stores rounded
    data too: a `Meet`, `Rise` or `Toric` edge its two stored surfaces and
    angles, a conic its rounded frame, none of them the arrangement's exact
    curve. Decisions. (1) *The given model.* A Boolean's result given to
    another Boolean is decided on its construction's exact model: the
    arrangement of its inputs' exact models that built it, run again (its
    inputs as kept by the result, the seams tried in the same order, the
    thread's arrangements cache), whose kept pieces are the body. Every face
    of the model is an input face holding kept pieces, on that input's
    exact surface in that input's frame (the arrangement's faces are then
    no longer on one model's frame: each face reaches its surface's frame
    and data through its input's model, a view); its region is where the
    input's exact `in_face` holds and the other input's sides at the point,
    pushed off the face both ways, are those the first operation keeps
    (`On` where a push stays on the other's boundary: the first
    arrangement's sections); its orientation is the kept pieces' (a cut's
    tool faces reversed). Membership is the first operation's set function
    over both inputs' exact membership with the same pushes. Every edge of
    the model is an arrangement edge on a kept piece's boundary (an input
    edge's exact line or arc, or a section of two input faces: S9c.1's
    lines, circles and ellipse arcs), every vertex an arrangement vertex
    (rational, or a quadratic surd), placed exactly: the stored vertices
    are not used. (2) *Its names.* The re-run's assembly gives the result's
    slots in the stored order (one deterministic assembly), so each model
    face carries the stored ids of the result faces holding its pieces
    (several where a slot through an input face left it in parts), each
    arrangement edge the id of the result edge it lies in (none inside a
    result face: a full circle's seam), each vertex the stored id where it
    is a result vertex (none where the result's edge runs on through it).
    The re-run is checked against the stored topology: equal slot counts
    and every re-run vertex rounding to the stored one bit for bit, or
    within the resolution after a rigid motion (a moved result moves its
    inputs, whose models are exact in the moved frames), else
    `ComputationLimit("a given result rebuilt differently")`; a result that
    was one of several solids of its Boolean is S9e.2's (its region is not
    its construction's set function alone). (3) *The second operation.* The
    given model enters `graph::arrange_shared` as either input unchanged in
    kind: its edges are pierced by the other input's faces and the other's
    edges by its faces (a line or a conic against a plane or a cylinder:
    S9c.1's exact meetings, so a new vertex on a given section is the
    meeting of three surfaces without a new computation), its faces meet
    the other's in S9c's sections restricted to both regions, its pieces
    are classified by the other's membership and the other's pieces by the
    given model's; tangencies, a vertex of one on the other's face, edges
    meeting and results touching themselves are S9c.1's `Degenerate`,
    including where the other input meets a given model's section at a
    point where the first result has no face (the region test's `On`). (4)
    *Names of the second result.* Faces, edges and vertices continue the
    given result's entities by S9a's rules over the given model's ids, so
    the second operation's history is over the first result's ids and
    chains onto the first's history: a given face kept whole keeps its id,
    in parts is `Split`, one on a face of the other input `Merged`; a new
    edge or vertex is `Generated` from the given result's faces and edges it
    lies on. A result piece on a model face carrying several ids continues
    the result face holding it: the second arrangement's pieces of that
    face joined across the other input's sections alone (never across a
    given edge) form one result face's part, one of whose pieces holds a
    given edge, whose side names it. A result edge along a given edge keeps
    that edge's stored circle or ellipse frame (its new ends' angles on it);
    other edges are rounded as S9c.1's. (5) *Sub-steps.* **S9e.1**: a
    Boolean's result from the curved arrangement whose inputs are prisms of
    lines, arcs and circles in any frames (S9c.1's pairs: every pair of
    faces meeting in lines, circles or ellipses), the only solid of its
    Boolean, given to another Boolean as object or tool with a prism of
    lines, arcs and circles in any position or with another such result,
    every pair of their faces again S9c.1's; a given edge at an irrational
    point (a vertical through a surd crossing) against a cylinder is
    `ComputationLimit`, as a section's line is in S9c.1. **S9e.2**: a stack
    with arc walls (S9a.2) and a polyhedral result (S9b) given with arcs,
    and a result of several solids: their stored topology comes from
    another assembly (the stack's slabs, S9b's fragments) or holds only
    part of the construction's region, so the re-run's arrangement is
    matched to the stored topology geometrically (each stored vertex the
    arrangement vertex within the resolution of it, one to one; edges by
    their ends and curves; faces by their edges) and a point's solid is
    decided by exact parity against the matched faces' regions; DRAW's
    `bcut_simple/L3` to `L6` (a pocketed disc, a stack, cut by a turned
    cylinder) are S9e.2's. **S9e.3**: given results with spheres, cones
    and tori (their inputs' ball, funnel and ring through the views) and
    with procedural edges (`Meet`, `Rise`, cone and torus sections,
    `Toric`) in the given model: such an edge against the other input's
    faces meets three surfaces off S9c.1's conics (resultants, certified
    isolation along the procedural curve), and deeper chains (a given
    result whose input is itself a result); DRAW's `bcut_simple/G9` and
    `H3` (a cone fused on a cylinder, then cut) are S9e.3's. **S9e.4**:
    imported bodies, a `Solid` from an imported topology of planes,
    cylinders, spheres, cones and tori (a new constructor), decided on its
    stored surfaces as exact affine models (the stored frames' binary64
    axes in rationals, as a prism's): an edge the meeting of its two faces'
    stored surfaces (the branch through its stored curve's midpoint, unique
    within the resolution), a vertex the point common to its faces' stored
    surfaces within the resolution of the stored one (unique, else
    `Degenerate`); stored faces tangent along an edge or on one surface
    across it (smooth joins, faces split at a seam) `Degenerate` until such
    edges keep exact data, and turned frames taken as their affine models
    as S9c's prisms are. Staying refused: spline faces and edges (S9f), and
    every rule by design listed in S9's remaining scope. (6) *Evidence
    first, S9e.1.* The case protocol chains: a `then OP ID` row (`swapped`
    after it: the first result is the tool) after the tool's rows and a
    third solid's rows, the first result's one solid the second's argument.
    An independent reference (`chained_boolean_reference.py`): S9c.1's
    slicing reference generalized to three prisms and a set function of
    three (slices in planes holding every axis, so the three axes take at
    most two directions; each slice's convex pieces of the first operation,
    then of the second against the third's parallelograms; every face of
    each input swept against both others, its pieces classified against
    each and kept where the function differs across the face, a face shared
    by two inputs counted once), nothing from the kernel; checks against
    S9c.1's pair reference for `X = A op1 B` (`V(X ∪ C) + V(X ∩ C) = V(X) +
    V(C)`, `V(X ∖ C) = V(X) - V(X ∩ C)`, `area(X ∪ C) + area(X ∩ C) =
    area(X) + area(C)` where no faces coincide), closed forms (boxes and
    coaxial cylinders) and Monte Carlo; fixtures of every class
    (`generate_chained_boolean_fixtures.py`: a box less a tilted hole, then
    less a vertical hole crossing it; a pin fused through a hole, then cut;
    a result as the tool; a third prism on a result's cap, coplanar; a
    plane through a result's ellipse edge; a slot leaving an input face in
    two result faces; exact and turned frames; declared `degenerate`: a
    third cylinder tangent to a result's cylinder, a third prism's edge
    through a result's vertex); a native capture (OCCT's both operations,
    the first result's solid the second's argument) before
    `solid/boolean/curved/given.rs` exists; then the kernel (the given model
    in `given.rs`, views in the arrangement), the probe and the comparison
    keyed on that file, histories chained, the `boolean` target's first
    curved results given to a second operation, the DRAW survey and a
    campaign.
  * **S9e.1 evidence (2026-09-30).** The case protocol chains: a `then OP
    ID` row (`swapped`) and a third prism's rows after the tool's
    (`identity_reference.encode_chained_case` and `native_chained_case`;
    `occt_boolean_oracle.cpp` runs the second Boolean on the first result's
    one solid, `tests/support/boolean_protocol.rs` likewise, the probe's
    rows the second's). `chained_boolean_reference.py` is S9c.1's reference
    generalized to three prisms and the chained set function: slices in
    planes holding every axis (at most two directions among the three), the
    first operation's convex pieces clipped from the first two sections'
    parallelograms and the second's against the third's (a convex piece's
    half-planes its edges), breakpoints over the three sections' lines
    (parallel lines of two sections coinciding, three lines of at least two
    sections concurrent, filtered on the closures of the sections they
    belong to); every face of each prism swept against both others at once
    (S9c.1's sweep against each, their crossings merged, breakpoints where
    a crossing of one meets a crossing of the other), a piece kept where the
    chained function differs across the face, one on faces of several
    prisms counted by the first of them. `generate_chained_boolean_
    fixtures.py --check`: 30 cases (10 chains; 10 fuses, 10 cuts, 10
    commons; 23 solid, 1 empty, 6 degenerate) in `XY`, `SIDE` and `TILT`: a
    box less a tilted hole halved by a box whose wall holds the hole's axis
    (the wall crossing the hole's ellipse edges), that result as the tool
    of a slab crossing the hole (the cut two solids), a box fused with a
    tilted pin stepped by a box whose wall holds the pin's axis and as the
    tool of a slab above it, a box on the holed top face (coplanar faces,
    edges crossing on them; the common empty), a box less a groove of
    radius 5/4 across its top (its top face in two result faces) drilled
    through one of them, a column of the hole's frame whose cap cuts the
    hole in a circle, a quarter cylinder bored coaxially; `degenerate` a
    third cylinder tangent to the groove along a generatrix and a third
    box whose top edge lies on the groove's wall (the line `y = 2.25, z =
    3`, `0.75^2 + 1^2 = 1.25^2`). Checks: closed forms (the groove drilled,
    the quarter bored) within 1.3e-40; S9c.1's pair reference for the first
    result `X` in `V(X u C) + V(X n C) = V(X) + V(C)` and `V(X - C) = V(X)
    - V(X n C)` (swapped `V(C - X)`), moments too, within 1.3e-41, and
    `area(X u C) + area(X n C) = area(X) + area(C)` in the 9 chains
    without a face of the third on another's within 7.0e-41; every face's
    classes their closed-form area within 2.6e-41; Monte Carlo 2.7
    standard errors; no near coincidence (the box on the holed top moved to
    `z` 7.5: its top plane's events at the tilted hole's end cap fell 3.9e-16
    from the hole's ellipse crossing its edge, the rounded `TILT` axes
    separating an ideal coincidence). `test_chained_boolean_reference.py`
    checks three boxes (one in `SIDE`) by their grid cells exactly for
    every chain, with shared planes of both orientations, a third prism
    far from the result, and the common's symmetry. Python 3.9 and 3.12
    write the same files; the generator's check a CI group of its own
    (`chained`, 26 s locally). The capture
    `occt-boolean-chained-preimplementation` (`compare_chained_boolean.py`,
    keyed on `solid/boolean/curved/given.rs`; the kernel `unsupported` on
    all 30, its first Boolean evaluating and the second refused as "a
    Boolean of a solid with curved faces or edges in any position (S9c)"):
    every result valid with the reference's solid count, all 30 matching
    (volumes within 5.7e-10, areas 8.0e-10, centres 5.9e-10 of the size;
    exact frames 2.9e-16), no review; seven solids' counts change when
    unified. S9e.1's kernel next.
  * **S9e.1 implemented** (`solid/boolean/curved/given.rs`, `model.rs`'s
    `Prism::view` and `given`, `graph.rs`, `assemble.rs`'s provenance
    `Made` and names, `mod.rs`): a Boolean's result of prisms with lines,
    arcs and circles, the only solid of its Boolean, given to another
    Boolean as object or tool with such a prism or another such result, on
    its given model as the refined decisions describe. All 30 fixtures as
    the reference (24 within the kernel's enclosures, each at most `1e-9`
    wide; the 6 declared degenerate refused: the tangent third cylinder as
    a tangency between the inputs, the third box's edge as an edge of one
    input on the other's face), every second history complete over the
    first result's ids (each of them resolved), results deterministic and
    moved rigidly, the first results translated with their third prisms
    keeping the reference's volumes, the grooved top's two result faces
    each named by its own part, the stored frames the reference's bit for
    bit (`tests/chained_booleans.rs`, 0.7 s at `opt-level` 2).
    `compare_chained_boolean.py` 24 matches and 6 reviewed for their
    entity counts (a diagnostic build of the oracle printing the unified
    solids' vertices, not committed, finds every native-only vertex on the
    hole's or pin's stored cylinder at its seam angle, or at the third
    box's imprint on the holed top, OCCT keeping those after unifying, and
    no kernel vertex missing natively); every other Boolean comparison
    unchanged in its counts (boolean 45/0, its splines 33/13, polyhedral
    43/2, curved 42/2, procedural 4/24, turned 2/13, capped 0/18, sphere
    30/0, spheres 12/21, cone 25/5, cones 21/19, torus 11/24, torus segment
    15/14, torus curved 15/29, spheres turned 0/18, cones' loops 5/26,
    torus parts 16/37). The refusals name their sub-steps now: a stack, a
    plane's piece or a polyhedral result with curved faces or against arcs
    (S9e.2), a result of several solids (S9e.2), a result of spheres,
    cones, tori or with procedural edges, deeper chains and a result
    against a sphere, cone or torus (S9e.3). Amendments and corrections,
    from the implementation: (a) parallel circular cylinders' relation
    (`CylPair::Parallel`) holds the second cylinder's circle in the
    first's frame, and an arc of the second input's cylinder was met
    against the first's with it as its own: its crossings were missed and
    the result left open (`InvalidTopology("an open Boolean of arcs in any
    position")` for a cylinder's cap circle across a parallel one's wall in
    exact frames, either input first; latent since S9c.1, no fixture
    reaching it): the relation is taken from the arc's side now
    (`tests/curved_booleans.rs`, the lens's closed form); (b) the
    independent history check takes an ellipse split into arcs as one curve
    (each arc's centre and axes' ends on the whole's), as it takes a
    circle: a given result's section edge cut by the second Boolean was
    reported `split_support_differs`; (c) a model face holding several
    result faces names each second-arrangement piece through its pieces
    joined across the edges they share and a model edge on the group's
    boundary (no point-in-face test needed); (d) a rigid motion by a
    rotation rounds the frames apart, and a box's wall parallel to the
    tilted hole's axis then lies within rounding of its direction (S9c.1's
    refusal, the first Boolean's too), so the moved inputs are checked
    under a translation by binary64 steps. The `boolean` fuzz target's
    chained stage gives first results with curved faces to the turned box
    too (`GIVEN_CURVED`, on): replaying the corpus (1,425 inputs) and the 19
    regressions with debug assertions, no failure, the slowest 3.6 s; of
    390 curved first results the stage reached, 26 evaluate on their given
    models, the rest refused as documented (199 results of spheres, cones
    or tori, S9e.3; 104 stacks, S9e.2; 34 a plane within rounding of a
    cylinder's direction; others S9c's). DRAW: of the boolean group's 25
    cases giving a result to another Boolean that the scan finds outside
    S9b.2's registered pockets, 17 need the public dataset, `bcut_simple/
    L3` to `L6` (stacks) are S9e.2's, `G9` and `H3` (a cone fused on a
    cylinder) S9e.3's, `bopcut_simple/ZQ1` and `bopfuse_simple/ZP6` refused
    as before: no status changes, the ledger does not (`L3`'s sentinel
    purpose now names S9e.2). Campaign: the boolean
    campaign at `51c08edf` (600 s, a sampled replay) clean, 1,101 runs, the
    slowest input 22 s under AddressSanitizer. DRAW survey: that of S9e.1
    and S9e.2, below (no status changes beyond S9e.2's rollex).
  * **S9e.2 refined, before its code (2026-09-30).** Why each is refused
    today: a stack given with an arc (or with curved walls of its own) goes
    to `polyhedra.rs`, whose `stored_model` refuses its cylinders and whose
    `prism_model` refuses the partner's arcs ("a Boolean of a prism with
    arcs and a solid other than a prism (S9e.2)"); an S9b.1 result
    (`Construction::Polyhedron` of two line prisms, built by
    `polyhedra::build`) likewise; one solid of several of an S9e.1-class
    result reaches `given::model`, which refuses it ("a Boolean's result of
    several solids given to another Boolean (S9e.2)"). Decisions. (1)
    *Which inputs.* Given to another Boolean as object or tool with S9e.1's
    partners (a prism of lines, arcs and circles in any position, or another
    given result): (a) a stack (`Construction::Stack`, S9a.2) with a
    cylindrical wall, or of planar faces only when its partner has an arc or
    a cylindrical face; (b) an S9b.1 result (both its inputs line prisms)
    under the same condition; (c) one solid of a result of several solids
    of S9e.1's class, of (a) or of (b). A stack or S9b.1 result with a
    polyhedral partner stays S9b.2's stored model (exact on planar faces,
    as now). A result whose inputs are not both prisms (a result of a
    stack, a stored model or a result: deeper chains) stays S9e.3's; a
    plane's piece (`Clipped`, `Half`) is a general body's, decided on its
    stored surfaces with S9e.4's imported bodies (its refusal names S9e.4).
    (2) *The construction's model is the curved arrangement's.* A stack is
    `A op B` of its two profiles' prisms (both on its frame: the object's,
    the tool's profile translated into it, over its heights `ha` and `hb`),
    an S9b.1 result `A op B` of its two kept prisms: the point set is the
    inputs' exact set function whichever arrangement built it. The given
    model re-runs S9c.1's arrangement of the two prisms' exact models
    (same-frame prisms included: walls parallel, caps coplanar, both seams
    apart on a shared circle) with S9e.1's seams and cache, its kept pieces
    the body, its faces, edges, vertices, regions and membership exactly
    S9e.1's; stack.rs's slabs and polyhedra.rs's fragments are not used
    (their exact data is 2D or triangulated, not S9c.1's faces, edges and
    pieces). (3) *A geometric match names it.* A stack's or S9b.1 result's
    stored slots come from another assembly, so the re-run's assembly (its
    solids, faces joined across shared edges, edges joined where they run
    straight on, closed curves as rings: S9c.1's rules, which S9a.2's and
    S9b.1's unified results follow too) is matched to the stored topology:
    each stored vertex the one re-run vertex within the resolution of it,
    one to one and onto; each stored edge the re-run edge between the
    matched ends whose points at a quarter, a half and three quarters of its
    range lie within the resolution of the stored edge's at those fractions
    (either direction; a ring's points within the resolution of the stored
    ring's curve); each stored face the re-run face bounded by the matched
    edges (the same set) on a surface of the same kind; the given solid the
    one re-run solid so matched (a stack's `index` numbers its own
    assembly's solids, not the re-run's). A stored entity unmatched or
    matched twice, or a re-run entity unmatched, is
    `ComputationLimit("a given result rebuilt differently")`: the
    assemblies differ (a seam, a join), which a fixture shows first. An
    S9e.1-class result keeps its slot check (one assembly; the solid its
    `index`). The model's faces, edges and vertices carry the matched stored
    ids as S9e.1's do, and its edges the matched stored curves. (4) *One
    solid of several.* Such a solid holds only part of its construction's
    region, and deciding a point's solid by parity along a ray would leave
    the quadratic surds (a ray's crossing of a cylinder) or need a ray
    parallel to every axis. The given model is instead the whole
    construction (every re-run solid's faces, edges and vertices; the
    construction's regions and membership), and the second arrangement's
    pieces are sorted by solid once it is built, by adjacency, exactly:
    each piece of a given face by a given edge on its group (the pieces of
    that model face joined across the edges they share; every group is
    bounded by given edges of one solid), each piece of the other input
    with a side inside the construction by a section or given edge on it
    (the given face's pieces along that edge), else by its neighbours
    across the other input's own edges inside the construction. Pieces of
    the other solids' given faces are dropped; a piece of the other input
    inside another solid is outside the given one (both sides). A shell of
    the other input inside the construction meeting none of its faces (the
    partner strictly inside one solid of several) has no adjacency:
    `ComputationLimit("a solid inside one of a given result's several
    solids")`. The arrangement holding every solid, a degeneracy between the
    partner and another solid of the given result (a tangency, an edge on a
    face) is S9c.1's `Degenerate` too, though the Boolean with the given
    solid alone is not degenerate: a conservative refusal, declared by a
    fixture. (5) *Routing.* `curved::applies` takes a stack or an S9b.1
    result when it has a cylindrical face or its partner an arc or a
    cylindrical face (a given result), with S9e.1's partners only (a
    sphere, cone or torus partner is S9e.3's). (6) *Evidence first.* The
    case protocol picks one solid of several: `then OP ID [swapped] [solid
    X Y Z]` takes the first result's solid holding the point strictly inside
    (`Solid::classify`; natively `BRepClass3d_SolidClassifier`), exactly one
    of them; without `solid` the first result must be one solid, as in
    S9e.1. The independent reference (`chained_boolean_reference.py`,
    generalized to a fourth prism): a selector `D`, a box sharing an axis
    direction with the others, the chain `op2(op1(A, B) and D, C)`
    (swapped `op2(C, op1(A, B) and D)`), with checks that `D` separates the
    first result (no area of `D`'s faces bounds `op1(A, B) and D`, the pick
    point inside it, the pair identities with `V(X_k)` from the three-prism
    chain `op1(A, B) and D`). Fixtures (`generate_given_boolean_fixtures.py`)
    of every class: DRAW's `bcut_simple/L3` to `L6` in the rollex's
    geometry (a disc less a pocket cut across its rim over the top part of
    its height, a stack; then a cylinder whose bottom lies on the pocket's
    floor, turned as `L3` and `L4` and in the stack's frame as `L5` and
    `L6`); a box with a boss (a stack) drilled across; S9b.1 results (a box
    and a turned box; a box cut by a tilted box, bored through the slanted
    face); a box severed by a cylinder (two solids of S9e.1's class), one of
    them drilled, as object and as tool; a box severed by a tilted slab (two
    S9b.1 solids), one of them against a cylinder; declared `degenerate`: a
    third cylinder tangent to a stack's arc wall along a generatrix, a third
    prism tangent to the solid of a severed box not picked. A native
    capture before `solid/boolean/curved/matched.rs` (the match and the
    sorting by solid) exists, the comparison `compare_given_boolean.py`
    keyed on it, `test_given_boolean_reference.py`, the generator's check a
    CI group (`given`); then the kernel.
  * **S9e.2 evidence (2026-09-30).** The case protocol picks one solid of
    several: `then OP ID [swapped] solid X Y Z`
    (`identity_reference.encode_chained_case` and `native_chained_case`
    with a pick point; `occt_boolean_oracle.cpp` takes the solid
    `BRepClass3d_SolidClassifier` finds the point inside, exactly one,
    `tests/support/boolean_protocol.rs` the one `Solid::classify` does).
    `chained_boolean_reference.py` takes a fourth prism, the selector `D`:
    each slice's first pieces clipped to `D`'s parallelograms before the
    second operation, every face swept against the three others (their
    crossings merged pairwise), a piece on several prisms' faces counted by
    the first in the order object, tool, third, selector, and
    `separation`, the area of `D`'s faces with the first result on either
    side. `generate_given_boolean_fixtures.py --check`: 36 cases (12
    chains; 12 fuses, 12 cuts, 12 commons; 30 solid, 6 degenerate; 18
    stacks, 6 S9b.1 results, 12 picked solids) in `XY`, `SIDE`, `DOWN`,
    `TILT` and `R125`: DRAW's rollex (a disc of radius 60 less a pocket of
    radius 40 across its rim over the top 6 of its height 20, then a
    cylinder of radius 30 whose bottom lies on the pocket's floor, crossing
    the pocket's wall and the rim) in `DOWN` (`L3`, `L4`) and in the stack's
    frame (`L5`, `L6`); a box with a boss bored coaxially (closed forms),
    sliced by a tilted slab (its common two solids) and as a tilted box's
    tool; a box fused with a box turned about its axis, bored through both;
    a box less a tilted box (a slanted top) bored through the slanted face;
    a box severed by a cylinder across it (two solids of S9e.1's class), the
    lower drilled and as a box's tool; a box severed by a tilted slab (two
    S9b.1 solids), the lower bored through its slanted face; `degenerate` a
    third cylinder internally tangent to the rollex's rim along a
    generatrix and a third box whose wall touches the severed box's other
    solid along a generatrix (the conservative refusal). Checks: closed
    forms within 3.3e-41; the pair identities for the given solid `X`
    (S9c.1's pair reference, or for a picked solid the three-prism chain
    `op1(A, B) and D`, one solid while the pair counts two) within 7.7e-42
    and the area identity within 8.2e-41; the selectors' separation exactly
    zero and each pick point inside its solid; face classes 2.2e-41; Monte
    Carlo 3.1 standard errors; no near coincidence outside the declared
    chains (the lower piece's bore first had its top `z = 5` and the slab's
    far plane meeting its wall where the selector's plane met the box's
    top, two triples concurrent at one slice: the bore raised to 5.5). The
    rollex's cut has area 30152.95448 (DRAW's `checkprops -s 30153`).
    `test_given_boolean_reference.py` checks four boxes (a box severed by a
    slab of the `SIDE` frame, the selector holding one part) by their grid
    cells exactly as object and tool, the separation of a selector across
    the first result (its face's area there), and a selector around the
    whole first result giving S9e.1's chain; S9e.1's fixtures are written
    unchanged by the generalized reference. The generator's check is a CI
    group of its own (`given`, 50 s locally on six workers). The capture
    `occt-boolean-given-preimplementation` (`compare_given_boolean.py`,
    keyed on `solid/boolean/curved/matched.rs`; the kernel `unsupported` on
    all 36, its first Boolean evaluating and the second refused as S9e.2's):
    every result valid with the reference's solid count, all 36 matching
    (volumes within 5.5e-9, areas 5.7e-9, centres 3.9e-10 of the size, a
    cylinder's piece under the tilted planes; exact frames 3.3e-15), no
    review; the rollex's results lose faces, edges or vertices when
    unified. S9e.2's kernel next.
  * **S9e.2 implemented** (`solid/boolean/curved/matched.rs`, `given.rs`'s
    construction and whole-construction model, `mod.rs`'s routing): a
    stack with a cylindrical wall (or of planes against arcs), an S9b.1
    result of line prisms given with arcs and one solid of a result of
    several, given to another Boolean as object or tool with a prism of
    lines, arcs and circles or another given result, on their
    construction's curved arrangement as the refined decisions describe.
    All 36 fixtures as the reference (30 within the kernel's enclosures,
    each at most `1e-9` wide; the 6 declared degenerate refused, each "a
    tangency between the inputs"), every second history complete over the
    given solid's ids (the other solid's ids neither in the result nor in
    its history), results deterministic and moved rigidly, given solids
    translated with their third prisms keeping the reference's volumes, a
    partner strictly inside one solid of several `ComputationLimit`, the
    stored frames the reference's bit for bit (`tests/given_booleans.rs`,
    1.6 s at `opt-level` 2). `compare_given_boolean.py` 36 matches, OCCT's
    unified entity counts on every evaluated case, no review; every other
    Boolean comparison unchanged in its counts (boolean 45/0, its splines
    33/13, polyhedral 43/2, curved 42/2, procedural 4/24, turned 2/13,
    capped 0/18, sphere 30/0, spheres 12/21, cone 25/5, cones 21/19, torus
    11/24, torus segment 15/14, torus curved 15/29, spheres turned 0/18,
    cones' loops 5/26, torus parts 16/37, chained 24/6). The stored
    slots of every stack and S9b.1 fixture match the re-run's assembly one
    to one (S9a.2's and S9b.1's unified results keep the curved
    assembly's vertices, edges and faces). Amendments and corrections,
    from the implementation: (a) pieces of both inputs on one surface
    joined into one result face, where the other input's piece is reversed
    beside a kept one (the rollex's cut: the cylinder standing on the
    pocket's floor, its bottom disc reversed over the material beside the
    pocket, joined with the floor), ran their loops about their own faces'
    opposite normals, and the face's loops did not close
    (`InvalidTopology("a joined face's loops do not close")`, latent since
    S9c.1, no pair of prisms having reached it): each piece's loops now
    run about the joined face's normal and every vertex's way on is chosen
    about it (`assemble.rs`); (b) a stack moved rigidly was given the
    height range `[0, 0]` (S9b.2's `Stack::moved`), so it classified every
    point off that height outside, and a result with a moved stack among
    its inputs classified its own vertices off its boundary (the moved
    results' check here; latent since S9b.2): the moved stack keeps its
    range; (c) a matched stored circle or ellipse whose frame turns against
    the conic's parameter is read the other way (`Given::flip`), since
    another assembly's frame need not follow the conic's (no fixture's
    does: S9a.2's and S9b.1's frames follow it as the curved assembly's
    do, so the reading is unexercised); (d) the
    remaining refusals of a plane's piece name S9e.4 (`polyhedra.rs`'s
    `prism_model` and `stored_model`, `model::Prism::new`), a given result
    with a sphere, cone or torus S9e.3 as before. The `boolean` fuzz
    target's chained stage now gives stacks with arc walls and a result's
    first solid of several too, and by the chained byte's upper bits its
    partner is a cylinder of radius 1.25 on the turned box's frame
    (`GIVEN_ROUND`, on): replaying the corpus (1,425 inputs) and the 19
    regressions with debug assertions, no failure, the slowest 4.4 s on a
    loaded machine; of the corpus's chained stages with a curved first
    result or the cylinder partner, 45 evaluate (28 curved first results
    against the box, 13 of them a result's first solid of several; 17
    against the cylinder), the others refused as documented (a result of
    spheres, cones or tori or a deeper chain, S9e.3; spline stacks, S9f;
    S9c.1's degeneracies). DRAW (`UPSTREAM_TESTS.md`): of the 25 cases
    giving a result to another Boolean that S9e.1's scan found unregistered
    (`L3` the sentinel), `bcut_simple/L3` to `L6` (the rollex) evaluate on
    both backends, Rust's `checkprops -s 30153` reading 30152.9544786437837
    (the reference's within 7.6e-17 relative) and the volume audit's
    volume within 3.2e-16 of the reference's: registered (987 registered,
    `L3` no longer a sentinel); the others as before (`G9` and `H3` S9e.3's,
    the dataset's restored arguments and constructs, private data, `ZQ1`,
    `ZP6`), none failing or timing out; the ledger does not change.
    DRAW survey: that of S9e.1 and S9e.2, below. Campaign: the boolean
    campaign at `6c221525` (600 s, a sampled replay) clean, 1,230 runs,
    the slowest input 30 s under AddressSanitizer (S9a.2's spline stack
    `d0a3de29`).
  * **DRAW survey of S9e.1 and S9e.2 (2026-09-30, `UPSTREAM_TESTS.md`).**
    The 1,802 cases of the Boolean group run again on both backends after
    S9e.1 and S9e.2 (the public dataset, 120 seconds a case): native
    DRAW's statuses unchanged, Rust's changed only where S9e.2's run of
    the chained cases found: `bcut_simple/L3` to `L6` (the rollex)
    evaluate, as registered there, and `G9` and `H3` (a frustum fused onto
    a cylinder, then cut) are refused as S9e.3's, not as S9c's solid with
    curved faces in any position. No other case evaluates newly or changes
    its refusal, none fails or times out, the sentinels are refused as
    before (`L3` no longer one). Rust evaluates 987 (983 before), 987
    registered, 592 refused (596), 223 unsupported on both. The volume
    audit of the 987: on the 983 of the last audit Rust's values bit for
    bit and native DRAW's to the digit, so S9e.1's and S9e.2's kernels and
    amendments change no registered case and no worker at an earlier
    commit was needed; the rollex's S9e.2's audit's bit for bit, the
    references' (`rollex_turned_cut`, `rollex_flat_cut`) within 3.2e-16
    relative in volume, 7.5e-17 in area and 1.5e-16 in the centre, native
    DRAW's within its printed digits. Rust and native DRAW disagree on the
    same 35 as before, native off in each. `bopfuse_simple/ZP6` and the
    `gdml_public` tori refused as before. A full contract run holds (the
    slowest Boolean case 9.8 seconds, the rollex 1.6); the ledger does not
    change; no kernel change.
  * **S9f refined, before its code (2026-09-30).** Why splines stop today:
    a spline prism in another frame, or whose offset or heights round,
    reaches `polyhedra::stored_model` ("a solid with curved faces or edges
    in any position (S9c)") or the curved engine's model refusal ("a spline
    profile in a Boolean of prisms in any position (S9c)", both to be
    relabelled S9f); no spline-surface intersection exists, every spline
    capability being a curve against a surface or 2D and rounded (S8b.3,
    S9a.2). Sub-steps, each evidence first, after S9e: S9f.1, spline prisms
    against polyhedral prisms in any position (a profile of lines, arcs and
    splines against line profiles, same-axis pairs whose offset rounds
    included): spline walls enter the curved engine (`Seg::Spline`,
    `Surf::SplineWall`, `Crv::Spline`); a plane meets a wall `o + S_x(t) x
    + S_y(t) y + w n` in `a(t) + b w = 0`, a crease `w = -a(t)/b` (an
    affine image of the profile, an exact B-spline, S8b.3's `plane_image`)
    or, parallel to the axis, generatrices at the roots of a degree-`p`
    polynomial; every vertex in `Q(alpha)` of degree at most `p`; exact
    membership at rational points or on a wall at a known parameter
    (S8b.3's and S9a.2's polynomials reused, their roundings not). S9f.2a,
    exactly parallel axes: spline walls against line, arc and spline walls
    in generatrices at exact 2D crossings (degrees `p`, `2p`, `pq`, the
    last up to 49). S9f.2b, cylinders with crossing axes: `A w^2 + 2B(t) w
    + C(t)` of degrees 0, `p`, `2p` along a ruling, branches as
    `MeetCrv` over the spline parameter, turning points and vertices of
    degree `2p`, loops by verified window graphs (as `torus_curved`), a new
    D13 curve kind. S9f.3, spheres and cones the same way. Refused, with
    these reasons: tori (a quartic in `w`, its turning points' discriminant
    of degree `12p`, traced curves only), spline walls against spline
    walls with crossing axes (degree `q` in `w`, no closed-form section, a
    discriminant of degree up to `q(q-1)p` and the implicit curve's
    extraneous branches; a wall of degree at most 2 on every span is a
    quadric and joins S9f.2b), spline results and imported spline bodies
    as inputs (S9e's general faces), rational and periodic profile splines
    (refused already). Risks: R4 under rounding (an interior knot of
    multiplicity `p` in a turned frame, its poles rounded: a fixture first,
    then exact knot removal before lifting or `PrecisionLoss`), U10 not
    reached (profiles and creases stay nonrational), a plane within rounding
    of a wall's axis (`Degenerate`, S9c.1's (c)), tangencies along a
    generatrix or at a knot (`Degenerate`), fields of degree 7 to 49 in
    `K` where S9a.2 rounded. Evidence to reuse: the 46 `boolean-spline-*`
    profiles in new frames, `boolean_reference.SplinePair`'s slicing and
    40-digit meetings, the spline capture's oracle and its reviewed
    BRepGProp errors, S8b's oblique-split fixtures for creases.
    * **S9f.1 evidence (2026-09-30), before its code (`BOOLEAN.md`).**
      `curved_boolean_reference.py` takes spline walls (S9c.1's assertion
      dropped): each spline segment is its exact Bezier spans, a line meets
      a span at the real roots of a degree-`p` polynomial (monotone runs,
      Newton steps at 40 digits, no resultant), an event between a spline
      crossing and a line's crossing moving linearly is an exact degree-`p`
      root, a tangency one of degree `p - 1`, knots are profile vertices;
      a spline span's wall is swept by its generatrices, and every
      operation's volume and moments are computed a second way by the
      divergence theorem over the classified face pieces.
      `generate_spline_any_boolean_fixtures.py --check` writes 38 cases
      (12 fuses, 13 cuts, 13 commons; 32 solid, 1 empty, 5 degenerate):
      S9a.2's bulge, dome, blob, wave, capsule and lens hole and a new
      `kink` in new frames, oblique creases, generatrices of planes
      parallel to the axis, tool edges piercing spline walls, coplanar caps
      of either orientation, a same-axis pair whose offset rounds
      (`bulge_offset`, both in `TILT`), turned (`R125`, `TURN`), leaning and
      tilted tools, the blob as a leaning tool; `degenerate` a plane tangent
      along the dome's apex generatrix, one tangent at the kink's knot, and
      a `TILTX` cap plane 8.9e-17 from a `TILT` blob's axis. Checks: the
      divergence theorem, inclusion and exclusion and the area identity
      within 1.2e-40; S8b's split reference on half-space boxes (creases
      and generatrices) 1.4e-41 exact, 1.9e-18 turned; S9a.2's `SplinePair`
      on `bulge_offset` (its offset taken exactly) and on its own 17
      one-spline fixtures 2.2e-17 (the 10 in `XY` to every digit); margins
      outside the declared pairs at least 0.019. The capture
      `occt-boolean-spline-any-preimplementation`
      (`compare_spline_any_boolean.py`, keyed on
      `solid/boolean/curved/spline_walls.rs`; the probe `unsupported` on
      all 38): every result valid with the reference's solids, 23 matches,
      15 reviewed (BRepGProp's default integration up to 6.7e-7, the wave
      1.0e-3; Green's theorem over OCCT's own faces and pcurves, a
      diagnostic build, within 1.3e-8 of the reference on all 38).
      Corrections to the decisions: (a) R4 bites before any Boolean: a
      profile spline C1 exactly at an interior knot of multiplicity `p`
      whose poles lose C1 when lifted into a turned frame (the knot at (5,
      4) along (-2, 1), poles (7, 3), (5, 4), (3, 5), in `TILT`: the knot's
      pole 2^-53 off its neighbours' midpoint) is refused by S8b's
      extrusion itself (`Solid::extrude_with`:
      `InvalidTopology("edge_not_c1")`, a panic of the protocol's
      construction), so the exact knot removal before lifting belongs to
      S8b.2's prisms, and a Boolean fixture of it waits for that; the
      fixture `kink_lean` keeps R4's knot of multiplicity two in `TILT` with
      poles that stay C1 lifted (its apex, horizontal), so the curved
      engine meets the knot as a C0 knot of the wall's surface
      (`test_spline_any_boolean_reference.py` checks both liftings); (b) the
      kernel refuses from two places, both to be relabelled S9f: a spline
      prism with lines only reaches `polyhedra::stored_model` (32 cases),
      one whose profile also holds an arc (the capsule) the curved engine's
      model (6); S9f.1 must route both to the curved engine; (c) "a plane
      within rounding of a wall's axis" arises from stored frames alone
      (`TILT` against `TILTX`: `n_A . (x_B * y_B)` is `-8.9e-17`, exactly zero
      for `TILT` against `SIDE`), so the kernel's test is S9c.1's (c) on the
      exact models, not on the frames' names; (d) S8b's oblique-split
      fixtures serve as closed forms (half-space boxes) rather than as
      Boolean fixtures. S9f.1's kernel next.

  * **S9e.3 refined, before its code (2026-10-02).** Why each is refused
    today. (1) `given::construction` takes a Boolean's result only when
    both its inputs are prisms: a result of a sphere, cone or torus (or of
    their caps, zones and parts) and a result one of whose inputs is itself
    a result, a stack or a plane's piece (a deeper chain) are refused as
    "a Boolean's result of solids other than prisms given to another
    Boolean (S9e.3)" (the fuzz target's ~180 chained stages refused as
    S9e.3's, DRAW's `bcut_simple/G9` and `H3`). (2) `given::model` keeps
    only lines and conics as the model's edges ("a Boolean's result with
    procedural edges given to another Boolean (S9e.3)") and reverses only a
    plane's or a cylinder's side for a cut's tool face. (3) `curved::build`
    refuses a given result whose partner is a sphere, cone or torus ("a
    Boolean's result given to another Boolean with a sphere, cone or torus
    (S9e.3)"). (4) Behind the refusals: `Given::view` goes down one level
    (a leaf is a primitive's model); the given model knows an edge's own
    surface only for a cylinder's conic (`own`), while a sphere's circle
    meets a sphere through its own sphere (`edge_surface`'s `own_ball`, the
    model's `ball`, none on a given model); `same_quadric` compares two
    faces' surfaces on their model's frame, but a given model's faces lie
    on their leaves' frames; `edge_surface` has no meeting for a section of
    two faces as an edge (`Meet`, `Rise`, `Toric`, a plane's section of a
    cone, `unreachable!`: only input edges were pierced); `assemble.rs`
    rounds a section's procedural curve through the arrangement's models by
    the curve's operand, which for a given edge are the first Boolean's.
    Decisions. (1) *Leaves of every kind.* A given result's construction is
    re-run whichever arrangement decided it: its inputs' exact models are
    `model_of`'s (prisms of lines, arcs and circles; spheres, caps and
    zones; cones and frusta; whole tori, v-segments and wedges; and given
    results themselves, below), its arrangement S9c's or S9d's
    (`arrange_shared` with every S9d pair it holds), named as S9e.1's (the
    re-run's slots the stored ones, checked) or matched as S9e.2's. Its
    faces are its inputs' faces holding kept pieces, on their exact
    surfaces through `view`, which now goes down to the primitive model
    (a leaf that is itself given passes its face to its own leaf); a face's
    outward normal is its leaf's composed with each level's orientation
    (`behind`), so a sphere's, cone's or torus's face kept as a cut's tool
    is reversed there and only there: no rule reads a side from a given
    face's own `Surf` (regions, normals and membership go through the
    views and the set functions), which now copies the leaf's surface
    unchanged for every kind. Membership and regions recurse: a level's
    set function over its two leaves' memberships, each leaf's own
    (`Prism::member`, `in_face`). (2) *Edges of every kind.* Every
    arrangement edge on a kept piece's boundary is a model edge with its
    first arrangement's exact curve and places: lines, conics (`Ang`), a
    sphere's circles (`Circle`), a plane's sections of a cone (`Cone`) and
    of a torus and the parts' rims (`Torus`), and the meetings of two
    curved faces (`Meet`, `Rise`, `Toric`; S9d.3b's cone pairs' `Meet`).
    Each model edge records its carrier: the model face on whose surface
    its place is measured (a cylinder for a conic, as now; a sphere for a
    sphere's circle, whose own sphere `edge_surface` needs; the cone for a
    cone section; the torus for a torus section or rim; the curve's own
    carrier for a meeting), through which the second arrangement reaches
    the leaf's data. (3) *Deeper chains.* A given result's inputs may be
    given results (of S9e.1's, S9e.2's or this class), each its own given
    model re-run and named as above, so the model is the construction tree
    evaluated level by level, each level's arrangement over its two
    leaves' models; the cache of arrangements keeps every level's re-run
    while the next is built. The tree is limited to three Booleans (the
    given result's and two below it: a result given twice more):
    deeper, `ComputationLimit("a given result's construction deeper than
    three Booleans")`, each level's re-run growing the next's fields and
    time (a kernel test shows the limit).
    (4) *A given result against a sphere, cone or torus.* The partner's
    faces meet the given model's through the views in the S9d pairs (the
    leaf surfaces' relations: a sphere against planes, cylinders, spheres
    and cones, a cone against those, a torus against every surface), the
    partner's edges meet the given faces as S9d's input edges meet such
    surfaces, and the second arrangement is S9d's; `build`'s refusal is
    removed. (5) *Meetings of given edges* (the second arrangement's new
    vertices on them). A given edge meets the partner's faces by S9c.1's
    and S9d's edge meetings through its carrier: a line or conic against
    every surface, a sphere's circle against planes, cylinders, spheres
    (its own sphere the carrier's ball), cones and tori, a torus's section
    or rim against planes (`rim_plane`: the line of the two planes against
    the torus) and against quadrics and tori (`rim_far`). One meeting is
    new and of the same degree as S9d.3a's: a plane's section of a cone
    against another plane meets it on the line of the two planes against
    the cone (`line_cone`, a quadratic), the roots kept on the section's
    branch and placed by the cone's angle (a section through the cone's
    apex is already refused). A given edge that is the meeting of two
    curved surfaces (`Meet`, `Rise`, `Toric`, a cone pair's meeting), or a
    plane's section of a cone against a curved face, met by a partner's
    face is three surfaces two of them curved, a quartic and more along
    the curve (resultants and certified isolation along the procedural
    curve, the refined decisions of S9e): **S9e.3b**,
    `OutOfDomain("a given result's meeting of two curved faces met by
    another face (S9e.3b)")`; such an edge the partner's faces do not reach
    stays in the model with its first arrangement's curve. (6) *Names and
    stored curves.* Names as S9e.1's; a result edge along a given edge
    keeps that edge's stored curve: a circle or ellipse its frame (its new
    ends' angles on it, S9e.1's `given_arc`), a given edge of any other
    curve kept whole its stored curve (reversed where the result runs it
    the other way), a piece of one cut by the partner rounded from the
    exact curve on its carrier's leaf (the curve's own model data, not the
    second arrangement's operands). (7) *Degenerate.* S9d's rules in the
    second arrangement (tangencies, a vertex of one on the other's face,
    edges meeting, a plane tangent to a sphere or a torus wherever it
    touches, a plane through a cone's apex, a torus's tube circle on the
    other's surface, a result touching itself), the given model's faces
    and edges as an input's. DRAW's `G9` and `H3` are such a case: the
    `pcylinder` of radius 1 about `(5, 0)` touches the frustum's top circle
    of radius 6 at `(6, 0, 4)` from inside (a rim tangent to the third's
    wall, where both seams lie: DRAW's comment), so they stay refused, now
    as "a tangency between the inputs" (S9c's rule), and OCCT's area
    727.481 is not compared. (8) *What stays refused, and why.* S9e.3b's
    meetings (5); deeper trees (3); configurations the `boolean` fuzz
    target switches off for time (`TORUS_PAIRS`, two whole tori;
    `CONE_PAIRS`, two cones; `TURNED_PARTS`, parts, caps and zones against
    tori in turned frames) are not given in its chained stage either (no
    first result is built), though the kernel takes them in its tests;
    the chained stage keeps its limit of 12 faces on the first result; a
    plane's piece (S9e.4); splines (S9f); and every rule by design of S9's
    remaining scope. (9) *Fuzzing.* The chained stage's first results of
    spheres, cones and tori reach the given model with the turned box and
    cylinder (`GIVEN_CURVED`, `GIVEN_ROUND`); by the chained byte's next
    bit the partner may be a sphere of radius 5/4 about the turned box's
    centre (`GIVEN_BALL`: a given result against a sphere), the per-input
    time under AddressSanitizer kept well under the target's 60 s (else
    the switch off, as the others). (10) *Evidence first, S9e.3a.* The
    case protocol chains further: a case may hold several `then` rows,
    each with the next solid's rows, the previous result's one solid (or
    picked solid) given to the next Boolean (the Python encoder, the
    native oracle and the Rust protocol). An independent reference
    (`chained_curved_boolean_reference.py`): S9d.4b.2's face sweeps
    (`torus_curved_boolean_reference.py`: every face of every input
    covered by two families of circles and lines, cut along each curve at
    the other inputs' surfaces' roots, measured by the divergence theorem)
    generalized to several inputs and a set function of their
    memberships, each piece classified by the other inputs' memberships at
    its point pushed off the face both ways (so faces of several inputs on
    one surface need no special case: a piece bounds the chain where the
    function differs across the face, and is counted by the first input
    whose boundary holds it), nothing from the kernel; checks: both
    families of every face, the inputs' closed forms, the pair identities
    for the first result `X` and the third (`V(X u C) + V(X n C) = V(X) +
    V(C)`, `V(X - C) = V(X) - V(X n C)`, moments, and the area identity
    where no face of the third lies on another's), closed forms of coaxial
    chains, S9e.1's slicing reference on chains of prisms (two engines),
    Monte Carlo; solids counted by a grid at two resolutions against the
    declared count. Fixtures (`generate_given_curved_boolean_fixtures.py`)
    of every class: results with a sphere (a domed box cut by a tilted box
    across the dome's circle, drilled through the circle), with a cone
    (DRAW's `G9` with the `pcylinder` moved clear of the top rim, the
    frustum's section by the cylinder's top crossed by a box; a box with a
    conical countersink), with a torus (a box grooved by a torus, cut by a
    tilted box across the groove's circles), with a procedural edge the
    partner does not reach (a cylinder fused on a sphere, `Rise`, then cut
    by a box clear of the meeting), deeper chains (a box less two crossing
    holes then cut by a tilted box; a domed box drilled then cut), given
    results against a sphere, a cone and a torus; declared `degenerate`:
    `G9` itself and a box tangent to a given result's sphere. A native
    capture before `solid/boolean/curved/chain.rs` (the given model's
    carriers, recursion and the cone section's meeting) exists, the
    comparison `compare_given_curved_boolean.py` keyed on it,
    `test_given_curved_boolean_reference.py`, the generator's check a CI
    group of its own (`given-curved`); then the kernel, its tests
    (`tests/given_curved_booleans.rs`), the fuzz target's stage, the DRAW
    survey of the chained cases and a campaign. **S9e.3b** after it, evidence
    first: the meetings of (5).
  * **S9e.3a evidence (2026-10-02).** The case protocol chains further
    Booleans: a case may hold several `then` rows, each followed by its
    solid's rows, the previous Boolean's one solid (or picked solid) the
    next one's argument (`identity_reference.encode_chained_case` and
    `native_chained_case` with `more`; `occt_boolean_oracle.cpp`'s stages;
    `tests/support/boolean_protocol.rs`'s `Case::more` and `given_by`).
    `chained_curved_boolean_reference.py` is S9d.4b.2's reference
    (`torus_curved_boolean_reference.py`: its inputs, curve families, roots
    along a curve and quadrature) generalized to several solids and a set
    function of their memberships: every face of every solid swept by its
    two families against the surfaces of all the others, each piece
    classified by the others' memberships at a golden point (pushed off the
    face by 1e-34 of the size both ways only for a solid with a surface on
    the face's, where its roots vanish: faces of several solids on one
    surface need no special case), a piece bounding the chain where the set
    function differs across the face, counted by the first solid whose
    boundary holds it; one sweep serves every expression of the solids (the
    given result, the chain, each solid). Two findings while building it:
    roots of another solid's surface at a curve's end within rounding (a
    face's edge on another solid's plane, as where boxes share edges) are
    dropped, else the structure flips on noise; and an event the scan
    misses (a feature narrower than its spacing) shows as a quadrature node
    of another structure, bisected from that node against the interval's
    middle (`refined`; none was needed in the fixtures). Solids by rays
    along `z` at two resolutions (S9d.4b.2's binary64 `FloatModel`, a
    prism's lower cap at its own height: S9d.4b.2's puts it at zero, where
    all its prisms start) against declared counts.
    `generate_given_curved_boolean_fixtures.py --check`: 48 cases (16
    chains; 16 fuses, 16 cuts, 16 commons; 42 solid, 6 degenerate; classes
    9 `sphere`, 15 `cone`, 6 `torus`, 3 `procedural`, 6 `deep`, 9
    `partner`) in `XY`, `SIDE`, `TILT` and `R125`: a box with a dome (a
    sphere of radius 3 about its top face's centre) with a `TILT` slab
    across the dome's circle and a cylinder through it (swapped); DRAW's
    `G9` body (a cylinder of radius 9 and height 3 fused with a frustum of
    radii 7 and 6 on its base, coplanar bottoms) with the `pcylinder` moved
    to `x = 4.5`, a box crossing the frustum's section by the cylinder's
    top and a coaxial bore; a box with a countersink (a frustum from radius
    1 to 3 at the top face) drilled coaxially; a box grooved by a torus of
    radii 3 and 1 (its top at `z = 1/2`: circles of surd radii) with a
    `TILT` slab across the groove and a cylinder through its outer circle;
    a sphere fused with a peg off its axis (a `Rise`) cut by a box clear of
    it; deeper chains: a box less a hole along `x` and a vertical hole, then
    a `TILT` slab, and the domed box drilled, then an `R125` box whose wall
    crosses the dome clear of the drill; a box less the hole along `x`
    against a sphere about its top face's centre, a frustum standing in it
    and a torus whose tube crosses the hole; `degenerate` `G9` itself and a
    box touching the dome's top. Checks: coaxial closed forms (rings about
    the axis by Simpson's rule between their kinks, horizontal faces and
    walls where the set changes) within 5.9e-41; the two families within
    3.3e-35; each solid's closed form 1.5e-40; the pair identities for the
    given result and the last solid within 4.1e-41 and the area identity
    1.5e-40 (15 chains); Monte Carlo 2.4 standard errors; every surface met
    at a sine of at least 0.24 on the scanned curves and every family's
    events at least 2.6e-5 apart outside the declared cases.
    `test_given_curved_boolean_reference.py` checks boxes by their grid
    cells exactly (frames of different axes, faces on one plane both ways,
    a fourth box), a hemisphere less a coaxial cylinder by closed forms and
    S9e.1's quarter bored against the slicing reference. Python 3.9 and
    3.12 write the same files; the generator's check a CI group of its own
    (`given-curved`, 23 CPU minutes, 8 on four workers). Fixture
    corrections before the capture, from the reference's runs and a first
    kernel's (outside the repository): the groove's slab at `w = 1` held a
    chord tangent to the torus's bottom circle along its parabolic
    direction (the reference's quadrature did not converge) and at `w =
    1.4` passed through a point of the torus's top circle (the kernel found
    a result thinner than the resolution), and at 1.5 to 2.1 its common
    kept a sliver of volume 0.03 whose ray counts differed between
    resolutions: moved to 1.8 to 2.4; the peg's wall held the sphere's axis
    (S9d.1's refusal of a section through a pole off the meridians): moved
    off it; the groove's drill left a crescent too thin for the rays:
    widened. The capture `occt-boolean-given-curved-preimplementation`
    (`compare_given_curved_boolean.py`, keyed on
    `solid/boolean/curved/chain.rs`; the kernel `unsupported` on all 48, its
    first Booleans evaluating and the one given a result refused as S9e.3's):
    every result valid with the reference's solid count, 33 matching
    (volumes within 1.6e-8, areas 1.0e-8, centres 2.1e-9 of the size; exact
    frames 1e-15), 15 reviewed (BRepGProp's default integration on faces
    bounded by approximated meetings of two curved faces, up to 6.4e-6;
    adaptively within 2.4e-9, unchanged at 1e-12); `G9`'s cut has DRAW's
    area 727.481; 18 results' counts change when unified. S9e.3a's kernel
    next.
  * **S9e.3a implemented** (`solid/boolean/curved/chain.rs`, `given.rs`'s
    leaves of every kind, `graph.rs`'s edge places and meetings, `meet.rs`,
    `cone.rs`, `torus.rs`, `assemble.rs`'s carriers and stored curves,
    `model.rs`, `mod.rs`'s routing, `history.rs`'s split supports): a
    given result's construction re-run whichever arrangement decided it
    (spheres, caps, zones, cones, tori and their parts, and given results,
    as leaves; views down to the primitive model, normals and the
    parameters' sign composed through each level's orientation), every
    first-arrangement edge a model edge with its exact curve and places
    (`Given::places`), deeper chains to three Booleans
    (`chain::MAX_DEPTH`), a given result against a sphere, cone or torus,
    as the refined decisions describe. All 48 fixtures as the reference
    (42 within the kernel's enclosures, each at most `1e-9` wide; the 6
    declared degenerate refused, each "a tangency
    between the inputs": `G9`'s rod touching the frustum's top circle and
    the box touching the dome), every stage's history complete over its
    inputs and chaining (each stage's given solid one of the previous
    stage's results by its ids), results deterministic and moved rigidly,
    given results translated with their last solids keeping the
    reference's volumes, a tree of four Booleans `ComputationLimit`, a
    peg's `Rise` met by a box S9e.3b's, the stored frames the reference's
    bit for bit (`tests/given_curved_booleans.rs`, 14 s at `opt-level`
    2). `compare_given_curved_boolean.py` 25 matches and 23 reviewed: the
    15 native measures reviewed before, and 22 results whose entity
    counts differ from OCCT's unified ones (the kernel's meetings of two
    curved faces in exact pieces switched at rational points and its
    periodic faces whole, where OCCT splits its approximated meetings at
    its own points and its faces at their seams; the holes' seam vertex
    OCCT keeps); every other Boolean
    comparison unchanged in its counts (boolean 45/0, its splines 33/13,
    polyhedral 43/2, curved 42/2, procedural 4/24, turned 2/13, capped
    0/18, sphere 30/0, spheres 12/21, cone 25/5, cones 21/19, torus 11/24,
    torus segment 15/14, torus curved 15/29, spheres turned 0/18, cones'
    loops 5/26, torus parts 16/37, chained 24/6, given 36/0; S9f.1's 38
    unsupported). Amendments to the decisions, from the implementation:
    (a) a plane's section of a cone normal to its axis is given as the
    circle it is (`Crv::Conic` on the cone's frame, its angle the cone's:
    the same places), so every surface meets it as a conic (`G9`'s body:
    the frustum's circle at the cylinder's top, met by a rod and a box);
    (b) a torus's plane section meets a curved face by S9d.4c's `rim_far`
    only where it is a circle of the torus at a fixed angle (a plane normal
    to or holding the axis, as a part's rims are); another plane's section
    met by a curved face is S9e.3b's: the fuzz replay found `rim_far`
    taking a general section (a box's face across a torus, then a sphere)
    and the second result open (`InvalidTopology`, three corpus inputs
    with the sphere partner); (c) a given edge of S9e.3b's curves is not
    met where one of its faces' surfaces is apart from the other face's
    (a plane parallel to a cylinder's axis beyond its radius, a plane
    beyond a sphere, a curved pair's relation `Apart`), so a meeting the
    partner does not reach stays given (`dome_deep`'s `Rise`); (d) a result
    edge over the whole of a given edge (every model edge of its stored
    edge, none split) keeps the stored curve, oriented by its ends; a piece
    of a cone's or a torus's section is rounded from the exact curve on its
    primitive model (`chain::curve_model`, also for every section of the
    second arrangement, whose carrier's frame and data are its view's);
    a given sphere's circle keeps its stored frame as a given conic does;
    (e) the independent history check takes a torus's plane section, a
    hyperbola and a parabola split by another Boolean as one curve (the
    same torus or frame and plane or size, whatever its range), as S9e.1's
    took an ellipse. The `boolean` fuzz target's chained stage now gives
    first results of spheres, cones and tori, and by the chained byte's
    next bit its partner is a sphere of radius 1.25 about the turned box's
    centre (`GIVEN_BALL`, on): replaying the corpus (1,425 inputs) and the
    19 regressions with debug assertions, no failure, the slowest 5.0 s on
    a loaded machine (the chained stage's operations 4.3 s at the
    slowest); of the corpus's 1,238 chained operations 659 evaluate (409
    before), none is refused as S9e.3's (468 before), 156 as S9e.3b's (a
    meeting of two curved faces, or a torus's general section, met by the
    partner), the others as documented (spline stacks, S9c's arcs off
    their circles, S9c.1's degeneracies, a plane within rounding of a
    cylinder's direction). DRAW (`UPSTREAM_TESTS.md`): of the 20
    unregistered cases giving a Boolean's result to another Boolean,
    surveyed on both backends, `bcut_simple/G9` and `H3` are refused now
    as a tangency between the inputs (the decisions' refusal; before, as
    S9e.3's), the others as before (5 private data, 6 a restored shape
    given to a Boolean, 5 constructs the adapter does not read, `ZQ1`'s
    wire, `ZP6`'s torus tangent to the other input); none evaluates newly,
    fails or times out, so none is registered; the ledger does not change.
    The local campaign at `6582f379` found a tilted stadium's flat wall
    4e-16 inside the tangent plane of the `GIVEN_BALL` sphere by the turned
    frame's rounding, met in a circle of radius 1e-8 the validator refused
    (`degenerate_curve`), so a plane crossing a sphere within the
    resolution of tangency (a cap no higher than it) is now `Degenerate`
    as S9d.2c's circles are (`sphere::plane_section`, also for two spheres'
    radical plane;
    `fuzz/regressions/boolean/crash-d238291d9edbe60570ae31479609845872d763c0.bin`,
    `tests/sphere_booleans.rs`).
    Campaign: the boolean campaign at `b2765f20` (600 s, a sampled
    replay) clean, 989 runs, the slowest input 22 s under
    AddressSanitizer (the first, at `6582f379`, found `d238291d`, a wall
    crossing a sphere within rounding of tangency: S9d.1's sections now
    refuse it). DRAW survey: that of S9e.3a and S9f.1, below (no status
    changes; `G9` and `H3` refused as a tangency, as above).
  * **S9f.1 implemented** (`solid/boolean/curved/spline_walls.rs`,
    `model.rs`'s `Seg::Spline`, `Surf::Spline` and `Crv::Spline`, `meet.rs`,
    `graph.rs`'s places and joints, `assemble.rs`'s spline curves and wall
    pcurves, `mod.rs`'s routing, `polyhedra.rs`; R4 in `topology.rs` and
    `split/spline.rs`; the validator's slivers in `spline_flux.rs`,
    `quadrature.rs` and `spline_deviation.rs`): spline prisms (profiles of
    lines, arcs and splines) against prisms of lines in any position,
    same-axis pairs whose offset rounds included, as the refined decisions
    describe. A spline segment is its exact Bézier arcs run as the profile
    runs (its run parameter the curve's own, mirrored for a reversed span);
    a plane meets its wall in the crease `w = -a(tau)/b` (the rounded curve
    S8b.3's `plane_image` of the piece, its wall pcurve S8b.3's
    `wall_pcurve`, both ends at their vertices' heights) or, parallel to
    the axis, in generatrices over the exact roots of `a` (degree `p`); a
    line meets the wall at its trace's roots and a cap edge or crease meets
    a plane at the plane's function's roots along it. Every root is the
    arc's parameter, rational or a generator of `Q(alpha)` (S9c.2b.2's
    `Gen`), one generator per root of one polynomial (its primitive
    coefficients and index), so every vertex on a wall is in `Q(alpha)` of
    degree at most `p` at a known parameter; membership is exact (a point
    of a segment found again at its generator or by `gcd(S_x - u, S_y -
    v)`, its side from the run's tangent; another rational point by the
    `+u` ray's crossings at the roots of `S_y - v`). Both refusal sites are
    routed to it (`polyhedra::stored_model` keeps a spline prism against a
    plane's piece, relabelled S9f; `model.rs`'s refusal is gone) and what
    stays refused is labelled S9f: spline walls against spline walls or a
    prism with arcs (S9f.2), a sphere or cone (S9f.3), a torus, a given
    result or a plane's piece, and a Boolean's result with spline walls
    given to another Boolean. All 38 fixtures as the reference (33 within
    the kernel's enclosures of volume, area and centre, each at most `1e-9`
    wide; the 5 degenerate refused: the dome's apex generatrix, the kink's
    knot, `TILT`'s blob against `TILTX`'s caps 8.9e-17 off its axis), and
    `boolean-spline-any-r4-*` (6 cases written after the capture, no native
    rows): R4's rounded knot in `TILT` creased across the knot by an `XY`
    box whose `x` planes, exactly parallel to `TILT`'s axis, cut it in
    generatrices, and the blob in `TILT` against a `SIDE` box (its caps
    exactly parallel to the axis: generatrices, where `TILTX`'s are
    refused); every history complete, results deterministic and moved
    rigidly (`tests/spline_any_booleans.rs`, 9 tests, 15 s at `opt-level`
    2, 15 s in release). `compare_spline_any_boolean.py` 22 matches and 16
    reviewed: the 15 native measures reviewed before (their reasons now
    state the kernel's outcome) and `capsule_stand_cut`'s counts (OCCT's
    unified result keeps the capsule's top spline edge split in four where
    the standing box's footprint crosses it, three vertices found by a
    diagnostic build printing the unified solid's vertices; the kernel
    joins pieces of one curve); every other Boolean comparison unchanged,
    `compare_split.py` 72/56 and `compare_brep.py --family spline` 10/3.
    Amendments to the decisions and the evidence, from the
    implementation: (a) R4 is fixed in the lifting by the preferred
    option: every lift or placement of a profile spline (a prism's walls and
    cap edges, a crease, every piece of S8b.3's splits, S9a.2's profile
    Booleans and stacks, and S9f.1's results: `topology::lifted_spline`,
    `spline_wall`, `placed_spline`, `spline::piece`) first removes one copy
    of each interior knot of multiplicity `p` where the profile is exactly
    C1 (`topology::c1_reduced`: the profile's own poles less the knot's,
    the same curve in the same parameter, C1 by construction), and a lift
    still not C1 is `PrecisionLoss`; a general body's rigid motion removes
    it where the motion's rounding breaks C1 (`moved_parts`). R4's profile
    of the evidence now extrudes in `TILT` (`InvalidTopology("edge_not_c1")`
    before); the stored prisms of profiles with such knots change in every
    frame (their poles less the knot's), no other stored geometry; (b) a
    crease's pcurve is the identity in `u` by rounded Greville abscissae,
    so its Bézier arcs end a rounding step across the wall's knot lines
    (`kink_lean`, the wave) and its control polygon may leave the wall's
    `v` range by a cap (the wave's crease near the top): the deviation
    bound already took slivers (S8b.3's (c)); now the deviation also boxes a
    control polygon past the domain by the curve's exact nonnegativity
    against the bound, and the exact Green path and the certified
    quadrature take both, the slivers integrated on the patch's polynomial
    with their error added (at most `2 (deg + 1) δ` of `u` travel times
    `|G|` of both patches over the sliver's box): the kink's common's
    volume enclosure was 14% wide (the strips), the wave's area 1.7%;
    (c) a plane holding a spline joint's vertical edge (between two spline
    walls) across its face (`blob_lean`: the box's wall `y = 4` through the
    blob's C1 joint at `(8, 4)`, an edge on a face, refused by S9c.1's rule)
    takes the edge's parts inside it as its face's, and meets the walls
    there, when the arriving and leaving segments lie on the plane's two
    sides; touching (both on one side, or a tangent) is `Degenerate` (a
    replayed variant found a lens hole's corner in a tilted hole wall,
    non-manifold before); (d) the ray's crossings with a spline count by
    the chords' predicate (an end strictly above the point): counted the
    other way a joint at the point's height flipped the parity (a replayed
    variant, a box apart from a bulge at its joints' height, an open
    result); (e) an edge of one input in the plane of a face of the other
    but off that face stays refused (S9c.1's rule), so the rounded-knot
    fixture's box stands off the plane of the prism's base edge; (f) a
    vertex found twice keeps one point, and an edge's place there is
    recomputed from that point (its field); (g) a line or curve over a
    spline tangent to the other's face is `Degenerate` with its own label
    (the dome's tangent fixtures refuse there first, before their
    sections). The `boolean` fuzz target decodes spline prisms against its
    turned, leaning, tilted and side tools through S9f.1 (`SPLINE_WALLS`,
    on; S9f's refusals accepted), and the byte at index 14, at or above 128,
    makes the object's spline R4's knot (`knot_profile`; earlier inputs
    decode as before but for that byte). Replays with debug assertions: the
    corpus (1,426 inputs) and the 23 regressions, none failing, the slowest
    6.9 s; 2,852 spline variants of the corpus (the object or the tool a
    spline prism, the tool turned, leaning, tilted or on its side, R4's
    knot for every other object), none failing after (c), (d) and (a)'s
    pieces (59 corpus inputs decoding R4's knot had found S9a.2's pieces
    across it refused, `InvalidCurve`), median 0.29 s, the slowest 2.95 s
    (23 s under AddressSanitizer on this host, where the corpus's slowest input takes 29 s); three kept as regressions
    (`fuzz/regressions/README.md`). 
    Campaign: the boolean campaign at `b2765f20` (600 s, a sampled
    replay) clean, 989 runs, the slowest input 22 s under
    AddressSanitizer (the first, at `6582f379`, found `d238291d`, a wall
    crossing a sphere within rounding of tangency: S9d.1's sections now
    refuse it). DRAW survey: that of S9e.3a and S9f.1, below (no case of
    the group reaches S9f.1's kernel).
  * **DRAW survey of S9e.3a and S9f.1 (2026-10-03, `UPSTREAM_TESTS.md`).**
    The 1,802 cases of the Boolean group run again on both backends after
    S9e.3a, S9f.1 and `b2765f20` (a plane crossing a sphere within the
    resolution of tangency `Degenerate`) (the public dataset, 120 seconds
    a case): no status changes on either backend. Rust evaluates 987, all
    registered, 592 refused, 223 unsupported on both, as before. The one
    reason changed is S9e.3a's run's: `bcut_simple/G9` and `H3` refused as
    a tangency between the inputs (the rod's circle on the frustum's top
    circle), not as S9e.3's. No case evaluates newly, none fails or times
    out, the sentinels are refused as before. S9f.1's kernel reaches no
    case (the group's spline solids come from `nurbsconvert`, which the
    adapter does not run); `b2765f20` turns no case `Degenerate` (the
    group's spheres are cut through their centres, `ZI4` to `ZJ3`, or by
    `ZP9`'s wall 1.053 from the centre of a sphere of radius 7.5, the
    box's other planes missing it by 9.5 or more) and changes no refusal.
    The volume audit of the 987: Rust's values the last audit's bit for
    bit and native DRAW's to the digit, so no worker at an earlier commit
    was needed; the same 35 disagreements, native off in each (three
    cases first ran with an empty script, the audit script's two backend
    threads rewriting one copy; run again they agree, and the script now
    writes one per backend). `bopfuse_simple/ZP6` and the `gdml_public`
    tori refused as before. A full contract run holds (the slowest Boolean
    case 13.5 seconds on a loaded machine, a worker at the last survey's
    commit as slow at that load); the ledger does not change; no kernel change.
  * **S9e.3b refined, before its code (2026-10-03).** Why each is refused
    today. A given result's edge on a meeting of two curved faces
    (`Crv::Meet`: a cylinder's or a cone's rulings against another quadric;
    `Crv::Rise`: a cylinder or cone and a sphere over the height;
    `Crv::Toric`: a torus and a quadric or another torus) or on a plane's
    section of a cone or of a torus other than a circle of the torus
    (`Crv::Cone`, `Crv::Torus`) is a model edge since S9e.3a, but where a
    face of the other input reaches it `meet::edge_surface` refuses it
    ("a given result's meeting of two curved faces met by another face
    (S9e.3b)"): the meeting of the edge's two surfaces with the partner's is
    three surfaces, two or three of them curved, and no S9 pair computes it
    (S9e.3a's only such meetings are a cone's section against a plane, on
    the line of the two planes, and a torus's circle at a fixed angle,
    S9d.4c's `rim_far`). Where one of the edge's faces' surfaces is apart
    from the partner's face (`graph.rs`'s `unmet`: a curved pair `Apart`,
    `chain::planes_apart`) the edge is kept as given; and `assemble.rs`
    refuses a result edge over part of such a given edge of a meeting of two
    curved faces (it keeps only a whole one's stored curve). In the fuzz
    corpus 156 of 1,238 chained operations stop there; DRAW's chained cases
    do not reach it. Decisions. (1) *Classes.* (a) *A plane among the
    three*: a given `Meet`, `Rise` or `Toric` met by the partner's plane, and
    a given section of a cone or a torus met by the partner's cylinder,
    sphere, cone or torus; (b) *three quadrics*: a given `Meet` or `Rise`
    met by the partner's cylinder, sphere or cone. (2) *Representation.*
    Every new vertex is a point whose three coordinates lie in one
    `Q(alpha)` (S9c.2b.2's fields: `alpha` a real root of a square-free
    eliminant `E` with rational coefficients, isolated by rational
    intervals, one generator per root, numbers polynomials in `alpha`; a
    `Qd` without its surd): no tower, because both parameters of the point
    are found in `Q(alpha)`. The two parameters: in (a), affine coordinates
    `(s, t)` of the plane, `p = p0 + s e1 + t e2` (`e1`, `e2` rational
    vectors spanning it: cross products with its normal), the other two
    surfaces restricted to the plane polynomials `G(s, t)` of degree 2 (a
    quadric's conic) or 4 (a torus's spiric curve); in (b), a ruled one of
    the three (the given curve's carrier first) parameterized by the
    half-angle tangent `x` of its angle in a chart rotated by a rational
    base direction (S9c.2b.1's `Chart`) and the place `y` along the ruling,
    `p = (A(x) + y B(x)) / Q(x)` (a cylinder's circle point and axis over
    `1 + x^2`; a cone's apex and generatrix direction), the two other
    quadrics restricted quadratics in `y` with coefficients of degree at
    most 4 in `x`. `E` is the resultant in the second parameter: of degree 4
    for two conics, 8 for a conic and a spiric curve, 8 for three quadrics
    (Bezout's counts; the chart's spurious factor `1 + x^2` removed). At
    each real root the fibre's two polynomials in the second parameter have
    their gcd over `Q(alpha)` taken by Euclid's algorithm (leading
    coefficients' exact signs at `alpha`): of degree one, the parameter is
    `-g0 / g1` and the point a rational map of the two; of degree zero, no
    point there (a root of the leading coefficients); of two or more, two
    points share the fibre (a plane parallel to a cylinder's rulings, a
    partner parallel to the carrier) and the projection is retried: the two
    parameters' roles exchanged, then rational shears `s + c t` (`c` = 1,
    -1, 2, 1/2, ...) of the plane's coordinates or, in (b), another ruled
    surface of the three as carrier; none separating the points is
    `ComputationLimit("a given meeting's points not separated by a
    projection (S9e.3b)")`. The chart's antipode (a rational direction) is
    tested exactly and the chart's base moved where it holds a point. (3)
    *Exactness.* Every real root of `E` is isolated over the whole line, so
    no meeting is missed; each point found is verified exactly (the three
    surfaces' functions zero at it, signs in `Q(alpha)`), and the given
    edge's curve keeps it only where its own exact test holds (`on`: both
    its surfaces, its branch, its window and range), so a point of the
    other branch, of another piece or of the same surfaces' other component
    is dropped; its place on the curve (the carrier's angle, a height, a
    torus angle's direction) is computed from the point in its field, and
    the second arrangement takes it as any vertex of an edge met by a face
    (on the face, inside it; its sections through it by their own exact
    tests; orders on a curve between numbers of different fields by
    enclosures, as S9c.2b.2's). (4) *Degenerate.* The partner's surface
    tangent to the given curve at a meeting on the edge: the three gradients
    linearly dependent there (`det(grad F1, grad F2, grad F3) = 0`, exactly:
    three surfaces through a point with dependent normals, which a double
    root of `E` at a point is), `Degenerate("a given meeting of two faces
    tangent to a face of the other input (S9e.3b)")`; two meetings of one
    given edge with one surface within the resolution of each other (a
    crossing within the resolution of a tangency, as S9d.4c's rims)
    likewise; a meeting at the given edge's end, S9's rule (a vertex of one
    input on the other's face); the given curve on the partner's surface
    (`E` zero identically and the curve on it), `Along`, taken as faces on
    one surface where one of the edge's faces lies on it, else S9's rule (an
    edge of one input on a face of the other); and S9d's rules for every
    pair in the second arrangement, unchanged. (5) *What stays refused, and
    why.* By cost, a torus among three curved surfaces (a `Toric` met by a
    curved face, a `Meet` or `Rise` met by a torus) and two tori's quartics
    in one plane (a torus's general section met by a torus, two tori's
    `Toric` met by a plane): eliminants of degree 16 (32 for two tori and a
    quadric) in the chart's variable, every vertex in a field of that degree
    and every later predicate a product there (S9d.4b.2b's two tori, the
    lightest of these, are already switched off in the fuzz target for
    time): `OutOfDomain("a given result's meeting met by a third surface of
    degree sixteen or more: a torus among three curved surfaces, or two
    tori in one plane (S9e.3b)")`; a projection that separates no points
    (2), a limit; everything else as before. (6) *Stored curves.* A result
    edge over part of a given edge of a procedural curve keeps the given
    edge's stored curve with the piece's ends' parameters on it (a `Rise`'s
    heights, a `Meet`'s or a `Toric`'s angle on the stored carrier, as
    S9e.1's `given_arc` for conics), a cone's or a torus's section's piece
    rounded as S9e.3a's (`chain::curve_model`); names and histories as
    S9e.3a's (a piece Modified from its given edge). (7) *Fuzzing.* The
    chained stage's partners (the turned box, the cylinder, `GIVEN_BALL`'s
    sphere) reach these meetings; a switch `GIVEN_MET` (on), switched off if
    the per-input time under AddressSanitizer comes near the target's 60
    s. (8) *Evidence first.* Fixtures (`generate_given_met_boolean_fixtures.py`)
    from S9e.3a's chained reference unchanged (`chained_curved_boolean_reference.py`:
    it decides the chain by its solids' memberships, so a partner across a
    given meeting needs nothing new), with an independent check of the
    meetings themselves (`given_met_reference.py`): each pair of surfaces of
    the given result's solids traced along the families of one's faces as
    roots of the other, the partner's surfaces' signs along the traced
    curve bracketed and the triple points refined by Newton's method in 40
    digits, kept where the given result's boundary has an edge there (its
    set function differs among the four quadrants about the curve) and the
    partner's face holds it; every chain outside the declared cases must
    hold such a point with every triple point's sine (the normals'
    determinant over their lengths) at least 0.05, and each declared
    degenerate chain its triple point with dependent normals within 1e-30.
    Chains of every class: a sphere with a peg (`Rise`) met by a `TILT`
    slab, by a wall through the peg's axis (parallel to its rulings), by a
    pipe, by a ball and by a frustum; two crossing rods (`Meet` rings) met
    by a `TILT` slab, by a wall parallel to a rod and by a ball; a frustum
    pierced by a pipe (a cone's and a cylinder's `Meet`) met by a slab; a
    torus with a rod through its tube (`Toric`) met by a `TILT` slab and a
    wall; a frustum cut obliquely (a cone's section) met by a rod and a
    ball; a torus cut by a wall off its axis (a spiric section) met by a
    pipe and a ball; declared `degenerate`: a box's face tangent to the
    peg's `Rise` at its lowest point, and a pipe through that point whose
    normal there is dependent on the sphere's and the peg's (no two of the
    three surfaces tangent). A native capture before
    `solid/boolean/curved/triple.rs` (the three surfaces' meetings) exists,
    the comparison `compare_given_met_boolean.py` keyed on it,
    `test_given_met_boolean_reference.py`, the generator's check a CI group
    of its own (`given-met`); then the kernel, its tests
    (`tests/given_met_booleans.rs`), the fuzz target's switch, the DRAW
    survey of the chained cases and a campaign.
  * **S9e.3b evidence (2026-10-03).** `generate_given_met_boolean_fixtures.py
    --check`: 50 cases of 17 chains (45 solid, 5 declared `degenerate`;
    given edges of classes `rise` 20, `meet` 12, `toric` 6, `cone` 6,
    `spiric` 6; partners' surfaces plane 24, cylinder 9, sphere 14, cone 3)
    in `XY`, `SIDE`, `TILT` and `TILTX`: a sphere of radius 5 fused with a
    peg of radius 1 about `(0, 2)` (a `Rise` loop from `(0, 3, 4)` up to
    `(0, 1, sqrt 24)`) met by a `TILTX` slab, by a box whose wall `y = 2`
    holds the peg's axis (parallel to its rulings, swapped), by a pipe along
    `x`, a ball and a frustum along `x`; a rod of radius 2 along `x` fused
    with a rod of radius 1 along `z` off its axis by 1/2 (two `Meet` rings)
    met by a `TILT` slab, a wall through the thin rod's axis and a ball; a
    frustum of radii 3 and 1 fused with a pipe through its wall (a cone's
    and a cylinder's `Meet`) met by a wall; a torus of radii 3 and 1 fused
    with a rod through its tube (`Toric` loops) met by a `TILTX` slab and a
    wall through the rod's axis; the frustum less a `TILT` box (an elliptic
    section) met by a rod and a ball; the torus less a box beyond `x = 5/2`
    (a spiric section) met by a pipe and a ball; declared `degenerate`: a
    box whose bottom face `z = 4` touches the `Rise` at its lowest point and
    a ball of radius 13/16 about `(0, 15/4, 69/16)` through that point, its
    normal there dependent on the sphere's and the peg's though no two of
    the three surfaces touch (fuse and cut: its common's pieces touch at the
    point, which rays count unsteadily). The rows are S9e.3a's chained
    reference's: the two families within 1.1e-35, each solid's closed form
    1.6e-36, the pair identities for the given result and the partner within
    9.2e-41 and the area identity 1.8e-40 (15 chains), Monte Carlo 3.2
    standard errors, the solid counts the declared ones by rays at two
    resolutions, every surface met at a sine of at least 0.41 on the
    scanned curves and every family's events at least 3.7e-6 apart outside
    the declared cases. `given_met_reference.py` finds the meetings
    themselves independently: for two solids of the given result, each face
    of one swept by its families, the roots of the other's surfaces along
    each scanned curve traced across the face, each partner surface's sign
    changes between traced points refined by Newton's method on the curve's
    two parameters in 40 digits, kept where the given result's set holds one
    or three of the four quadrants about the meeting (an edge of it) and the
    partner's membership changes across its surface: 34 triple points, each
    found at least twice (from both solids' faces or both families), every
    chain outside the declared cases holding at least one of its class on
    the partner's surface of its kind, every sine (the unit normals'
    determinant) at least 0.38; each declared tangency at its point on all
    three surfaces within 1e-30 with its normals' determinant zero, on the
    given result's edge and the partner's face. Amendment to the decisions'
    "unchanged": the declared tangencies' sections touch on the faces
    through the point, where the chained reference's structure flipped on
    rounding noise (a crescent between two section curves classified at its
    golden point once shorter than 1e-37, an event found again and again)
    and its quadrature met endpoint singularities converging slowly; for
    them alone (`MERGE`, `QUAD`, off by default) roots of different surfaces
    within 1e-30 of a curve's range are one and the quadrature takes 1e-22
    of the size to the fourth, their checks at 1e-18 (their rows rounded to
    binary64; the kernel refuses them); S9e.3a's fixtures regenerate byte for
    byte. Fixture corrections before the capture, from the reference's runs:
    slabs grazing a torus's top or leaving thin caps on a rod were moved
    (solid counts unsteady between the rays' resolutions); a frustum
    partner standing in the peg (thin crescents between near-parallel
    walls) laid along `x`; a wall through `(3, 0)`, tangent to the torus's top
    circle (a family curve tangent to it, the quadrature not converging),
    moved to `x = 3.1` with the rod; a rod internally tangent to the
    frustum's base circle moved off it; a pipe tangent to the meeting (its
    section's fourth-order contact with the sphere's parallels) and a ball
    with the tangency at its pole replaced by the ball above; slabs cut to
    their solids' extent (Monte Carlo samples the inputs' common box).
    Python 3.9 and 3.12 write the same files; `test_given_met_boolean_reference.py`
    checks the triple points in closed form (`(+-1, 2, 2 sqrt 5)` on the
    wall through the peg's axis, `(0, 3/2, +-sqrt 7 / 2)` and `(0, -1/2,
    +-sqrt 15 / 2)` on the wall through the thin rod's), on all three
    surfaces with their sines for the elliptic section and the rod, the
    declared tangency, and a meeting of the two surfaces off the given
    result's edges; the generator's check a CI group of its own (`given-met`).
    The capture `occt-boolean-given-met-preimplementation`
    (`compare_given_met_boolean.py`, keyed on `solid/boolean/curved/triple.rs`;
    the kernel `unsupported` on all 50, its first Booleans evaluating and the
    one given a result refused as S9e.3b's): every result valid with the
    reference's solid count, 12 matching, 38 reviewed (BRepGProp's default
    integration on faces bounded by approximated sections, up to 4.8e-6;
    adaptively within 3.2e-8, unchanged at 1e-12, the smallest results the
    farthest); 12 results' counts change when unified. S9e.3b's kernel next.
  * **S9e.3b implemented** (`solid/boolean/curved/triple.rs`, `meet.rs`'s
    routing, `graph.rs`'s clip boxes, `assemble.rs`'s pieces of given
    meetings, `chain.rs`'s `piece`; the curves' `surfaces`; the validator's
    holes in torus bands wound in `v`, `polyhedra.rs`'s face boxes,
    `history.rs`'s split meetings, `num.rs`'s kept exact signs): a given
    result's edge on a meeting of two curved faces (`Meet`, `Rise`,
    `Toric`), or on a cone's or a torus's general plane section, met by a
    face of the other input, as the refined decisions describe. With a
    plane among the three, the other two restricted to the plane's rational
    affine coordinates (an integer basis of cross products with its normal,
    its point on the axis of the normal's largest component); with three
    quadrics, the two others along a ruled one's rulings (the curve's
    carrier, else the partner's) over a rational chart whose antipode holds
    no common point; the restrictions made primitive, their resultant by
    Sylvester determinants at integers and Newton's interpolation; at each
    real root of its square-free part the fibre's gcd over `Q(alpha)`, every
    point in one `Q(alpha)` of degree at most 8, verified on all three
    surfaces and kept by the curve's own test; dependent gradients there
    `Degenerate`, two meetings within the resolution likewise; a torus among
    three curved surfaces and two tori in a plane `OutOfDomain` by cost.
    All 50 fixtures as declared (42 within the kernel's enclosures, each at
    most `1e-9` wide; the 8 declared degenerate refused: `peg_touch` and
    `peg_kiss` as the partner tangent to the given meeting, `cone_cut_rod`
    by S9d.3a's rule, below), every stage's history complete and chaining,
    results deterministic and moved rigidly, given results translated with
    their last solids keeping the reference's volumes, and three more: the
    elliptic section met by the rod lowered clear of the apex, the crossing
    rods' meeting met by a wall along both rods' rulings (a fibre of two
    points, the projection sheared) and a torus's meeting with a ball
    refused by cost, each consistent by the kernel's pair identities
    (`tests/given_met_booleans.rs`, 9 tests). `compare_given_met_boolean.py`
    8 matches and 42 reviewed: the 38 native measures reviewed before, and
    the kernel's entity counts against OCCT's unified ones in 32, 4 of them
    alone (its exact
    meeting pieces and whole periodic faces, OCCT's split approximations
    and seams); every other Boolean comparison unchanged in its counts.
    Amendments to the decisions and the evidence, from the implementation:
    (a) `cone_cut_rod` is declared `degenerate` after the capture: the rod's
    top cap's plane `z = 6` holds the frustum's apex, S9d.3a's refusal of a
    plane through a cone's apex, which the kernel found (its inputs, rows
    and native observations unchanged; the generator's check and the
    capture's reproduce); (b) the validator decides a hole in a torus face
    wound in `v` (a plane's cut across the tube leaves a band) by the `+u`
    ray's signed crossings with the other loops, every `v` alias, the
    band's right side running in `+v` (the `+v` rule of a band wound in `u`
    with the parameters exchanged): before, such a hole was always
    `UncertifiedContainment`, a `ComputationLimit` (`torus_rod_tilt`'s fuse
    and cut, the rod's hole in the band); (c) a Boolean result's bounds take
    each torus face's whole box and each cylinder's or cone's circles
    between its edges' heights (sampled, padded), as a sphere's before: a
    chained result's vertices classify against the given result, whose
    bounds from its edges alone left them outside (`torus_cut_*`'s and
    `torus_rod_wall`'s rigid motions); (d) the independent history check
    takes a split meeting (`Meet`, `Rise`, `Toric`) as one curve: the same
    surfaces and branch whatever its range; (e) a result edge over part of
    a given meeting keeps the stored curve from its start's parameter in
    the stored frame turning (or rising) by the places' sweep, and over
    part of a torus's section at a fixed angle stored as a circle an arc of
    that circle (the replay's history check found such pieces rounded as
    spiric sections); (f) a vertex at a pole of a sphere face it bounds
    stays in the result (the replay: a given meridian joined across the
    pole, its pcurve undefined there); (g) a budget: an eliminant whose
    coefficient bound (each restriction's degree times the other's largest
    coefficient, in bits) passes 4,096 bits is `ComputationLimit("a given
    meeting's eliminant past its budget of coefficient bits (S9e.3b)")`
    (the replay: torus sections in turned frames against the turned
    cylinder, eliminants of 8,000 to 11,000 bits, 50 to 170 s an input with
    debug assertions; the fixtures' at most 1,300); (h) for time: the
    eliminant of three surfaces is kept for every piece of the meeting and
    every face on the partner's surface, a simple root of it needs no check
    of its fibre, and each root's point is enclosed in binary64 intervals
    (Euclid's algorithm on the fibre in intervals, the restrictions scaled
    by a power of two) so a point certainly outside the edge's and the
    face's boxes is never constructed; the exact signs a generator decides
    by Sturm-Tarski are kept by polynomial (a vertex's surfaces are tested
    on every section through it).
    The fuzz target's chained stage may centre its partner on the first
    result's first meeting of two curved faces or torus section by the
    chained byte's next bit (`GIVEN_MET`), off: with it on, the corpus's
    slowest chained operations reaching such an edge took 60 to 71 s an
    input under AddressSanitizer on the Mac (8 to 10 s with debug
    assertions: vertices in fields of degree eight through the second
    arrangement), past the target's 60 s; off, a partner whose bounds meet
    the sampled box of a `Meet`, `Rise`, `Toric` or torus section edge of
    the first result is not given it. Replays with debug assertions of the
    corpus and the regressions (1,454 inputs, one process each), none
    failing: as committed (`GIVEN_MET` off) the slowest input 4.8 s (as at
    the base, not S9e.3b's; 34 s under AddressSanitizer), 684 of 1,176
    chained operations evaluating, 16 still reaching S9e.3b's path past the
    filter and refused as its limits or cost (the slowest such input 21 s
    under AddressSanitizer); with it on, 767 of 1,296 evaluating (659 of
    1,238 at S9e.3a), the slowest input 10.2 s, S9e.3b's refusals 16
    projections separating no points, 14 eliminants past the budget, 4
    meetings within the resolution and 2 by cost (156 refused as S9e.3b's
    before); with every chained partner centred on a meeting where the
    first result has one (`FORCE_MET`, 168 operations, 98 of them
    evaluating), the slowest input 15.2 s, 753 of 1,296 evaluating in all.
    Every check of HANDOFF's verification passes: formatting, clippy, the
    1.85 check, the test suite, every comparison unchanged with 0 failures
    (`compare_split.py` and `compare_brep.py --family spline` too), the
    reference tests, the generator's check and the ledger. DRAW
    (`UPSTREAM_TESTS.md`): the 20 unregistered cases giving a Boolean's
    result to another Boolean, surveyed on both backends, as after S9e.3a
    (the kernel reaches none of them: 5 private data, 6 a restored shape
    given to a Boolean, 5 constructs the adapter does not read, `ZQ1`'s
    wire, `G9` and `H3` a tangency between the inputs, `ZP6` a torus tangent
    to the other input); none evaluates newly, so none is registered and no
    volume audit is due; the ledger does not change. DRAW survey: that of
    S9e.3b, S9f.2a and S9f.2b.1, below (no status or refusal changes; the
    volume audit's values bit for bit).
    Campaign: the boolean campaign at `b0b9adc6` (600 s, a sampled
    replay) clean, 995 runs, the slowest input 21 s under
    AddressSanitizer (`GIVEN_MET` off).
    Amendment (the parallel track "The degree-eight arrangement's
    arithmetic", branch `given-met-speed`): `GIVEN_MET` is on, the
    corpus's slowest chained operations reaching a given meeting about a
    third of their time before; the 8 "points not separated by a
    projection" of the replay with it on (two cylinders with parallel axes
    meeting in two rulings a sphere crosses at equal heights) are separated
    by a projection along a ruled surface's place sheared by its angle's
    chart, tried after the ones above.
  * **S9f.2a refined, before its code (2026-10-03).** Why it is refused
    today: `curved::spline_pairs` refuses a spline prism against a prism
    holding arcs or circles ("a spline prism against a prism with arcs in
    any position (S9f.2)") or splines ("spline walls against spline walls
    in any position (S9f.2)"); behind that refusal `meet::section` and
    `meet::edge_surface` have no spline wall against a curved face
    (`spline_curved`), `SplineSeg::locate` finds a point of `Q(alpha)` on
    a segment only at its own arcs' generator, and a spline profile's
    membership is undecided at an irrational point off its segments
    (`OnProfile::Undecided`: none arose against S9f.1's polyhedral
    partners). Decisions.
    (1) *Domain.* Two prisms whose stored normals are exactly parallel
    (`n_A x n_B = 0` in rationals: bitwise equal or opposite, which a rigid
    motion of both keeps), at least one holding a spline segment and the
    other an arc, a circle or a spline (a prism of lines only is S9f.1's in
    any position), in frames that differ: turned about the axis (`TURN`'s
    quarter turn and a half turn exactly, `R125`'s rounded rotation,
    `TILT2` against `TILT` about the tilted axis), or with equal axes and an
    origins' offset that is not binary64 in the object's frame. The
    relation to S9a.2: S9a.2 decides two profiles in one exact frame
    (bitwise-equal axes, a binary64 offset, the routing unchanged): one 2D
    arrangement of profiles, its spline crossings found by the same
    resultants (`implicit_on`, degree `p q`) but rounded to binary64
    parameters at the profile level and paired by certified boxes, and
    stacks of slabs. Here the tool's profile lies in the object's frame
    under a rational affine map that is neither the identity nor a
    binary64 translation (a circle becomes an ellipse with rational
    implicit coefficients unless the map is a similarity), so no profile
    arrangement can be rounded first: the walls are faces of the curved
    engine (S9f.1's `Surf::Spline` and S9c's cylinders), meeting in
    generatrices over exact 2D crossings, every vertex exact and rounded
    once at assembly. Crossing axes stay refused: against a cylinder
    S9f.2b's, against a spline wall by design.
    (2) *Meetings*, all walls parallel, so every wall against wall section
    is generatrices: (a) a spline wall (degree `p`) against a cylinder
    wall: the cylinder's frame coordinates `(u', v')` of the spline's
    profile point are rational affine functions of `(u, v)` (the axis
    carries no `u'` or `v'`), so on each Bézier arc `E(tau) = (u'(S(tau)) -
    c_u)^2 + (v'(S(tau)) - c_v)^2 - r^2`, degree `2 p`, exact; generatrices
    over its roots strictly inside the spline's run, every vertex there in
    `Q(alpha)` of degree at most `2 p`, `alpha` the arc's own parameter
    (S9f.1's generators); (b) a spline wall against a spline wall: on each
    pair of arcs (degrees `p` and `q`, the second's control points mapped
    exactly into the first's frame), `R(tau) = f(S(tau))` with `f(X, Y) =
    Res_sigma(B_x(sigma) - X, B_y(sigma) - Y)` the second arc's implicit
    equation (Sylvester, exact, of degree `q`; S9a.2's polynomial), degree at
    most `p q`; at each root `alpha` in `[0, 1]` the second arc's parameter
    is an element of `Q(alpha)`, `sigma = -b(alpha) / a(alpha)` from the
    first subresultant `a sigma + b` of `B_x(sigma) - X` and `B_y(sigma) - Y`
    at the point (`a`, `b` polynomials in `(X, Y)` from exact evaluations
    and interpolation, once per arc), the crossing kept when `sigma` lies
    in `[0, 1]`; `a = 0` there (a node or cusp of the second arc's curve)
    is `ComputationLimit`. One field per crossing whichever way it is
    found: always the object's (operand A's) spline parameter's; (c) a cap
    edge or crease over a spline segment against the other's cylinder or
    spline wall: the same polynomials along its segment (it projects onto
    the segment), its place the run parameter; an arc's or circle's cap
    edge (a conic) against a spline wall: (a)'s roots, the point on its
    own cap and its angle `((u' - c_u) / r, (v' - c_v) / r)`; a spline cap
    edge of the second input against the first's spline wall: (b)'s
    crossings in the first's field at the second's run parameter; (d) line
    walls against spline walls S9f.1's (generatrices of degree `p`), arcs
    against lines or arcs S9c's (parallel cylinders: circles, or S9c.2b.2's
    lines through algebraic points).
    (3) *Places and membership.* A point of another field on a spline
    segment (a crossing found on the other wall) is found on it exactly:
    its own generator's arc first (S9f.1), then on each arc whose control
    box holds it the exact sign of the arc's implicit equation at the
    point, and where that vanishes the subresultant's parameter, on the arc
    when in `[0, 1]` with the arc's point there the given one. An
    irrational point off every segment (a plane's generatrix on one spline
    wall, or a crossing with a cylinder, tested against the other spline
    profile) is decided at a rational point of a box about it certified
    free of the profile's boundary (its enclosure narrowed until every line
    misses the box by its corners' sides, every circle by distance, every
    spline arc by its control boxes under exact subdivision), by S9f.1's
    ray there; no such box within 2^-64 of the point's size leaves it
    undecided (refused where it matters, as S9f.1's).
    (4) *Cost.* The fields' degrees are at most `2 p <= 14` against a
    cylinder and `p q` against a spline (up to 49 for two septics); arcs
    whose degrees' product exceeds 16 are refused (`ComputationLimit("two
    spline walls of degrees whose product exceeds 16 (S9f.2a)")`), checked
    by the kernel's tests; S9a.2's and the boolean target's profiles are of
    degree at most three (products at most 9).
    (5) *Degenerate.* Walls tangent along a generatrix: a root of `E` or `R`
    of multiplicity above one inside both runs (a spline touching or
    osculating a circle or another spline), as S9f.1's tangent planes and
    S9a.2's tangent splines; a touch at a knot of either (a root at a knot
    of multiplicity above one on either arc beside it, or the second's
    parameter at its knot with the curves not crossing); a vertical edge of
    one on the other's curved wall (a crossing at a segment's end) is the
    engine's edge on a face, as before; a crossing within the resolution of
    another or of a vertex the engine's own. Coincident walls (`R`
    identically zero on a pair of arcs: one curve in both frames, as when a
    symmetric profile is turned onto itself) are refused as
    `OutOfDomain("spline walls of both inputs on one surface in different
    frames (refused, S9f)")` (S9a.2 shares one curve in one frame only); `E`
    cannot vanish on an arc (a polynomial arc is no circle or ellipse).
    Walls within rounding of tangency are not refused exactly: their
    results validate or are `PrecisionLoss`; fixtures keep a margin.
    (6) *Stays refused*: crossing axes (S9f.2b, by design for two spline
    walls), spheres and cones (S9f.3), tori, given results with spline
    walls, coincident spline walls, rational and periodic profile splines.
    (7) *Evidence first*: `curved_boolean_reference.py`'s parallel slicing
    and face sweeps extended to spline chords against arcs' and splines'
    chords, their events the 2D crossings (a spline against a circle mapped
    into its frame by the exact roots of its degree-`2 p` polynomial, two
    splines by S9a.2's reference's subdivision in fractions and 40-digit
    Newton steps, no resultant); `generate_spline_parallel_boolean_fixtures.py`
    with spline walls against a stadium's line and arc walls, against discs
    and a holed square, a capsule against a disc (its arc against the
    disc's), spline against spline (turned, quarter and half turns, rotated
    about the tilted axis, rounding offsets), cubics against cubics and a
    quartic against a cubic, and declared degenerate tangencies (a dome
    touching a disc, a dome touching a half-turned dome); checked two ways
    (slicing and the divergence theorem), by inclusion and exclusion, the
    area identity, a second slicing axis, S9a.2's `SplinePair` where the
    map is exact (quarter and half turns, the tool's profile turned
    exactly; rounding offsets taken exactly), margins (crossing sines,
    near misses, vertices); a native capture before
    `curved/spline_parallel.rs` exists (`compare_spline_parallel_boolean.py`
    keyed on it); then the kernel in that file, its tests, and the boolean
    target's spline variants against arcs and splines in the frame turned
    about the axis (`SPLINE_PARALLEL`) and a rounding offset in the tilted
    frame.
    * **S9f.2a evidence (2026-10-03), before its code (`BOOLEAN.md`).**
      `curved_boolean_reference.py` takes arcs, circles and splines against
      spline walls on exactly parallel axes: the 2D crossings of the
      profiles' curves projected along the axis (`parallel_crossings`, in
      the exact affine map between the frames: a span against a circle at
      the real roots of the exact degree-`2 p` polynomial, a tangency a root
      of even multiplicity; two spans by subdivision in fractions and
      40-digit Newton steps, on the distance's gradient where tangent, no
      resultant) are breakpoints of the slicing and of every face sweep
      (S9f.1's assertions against a surd or a second spline chord now
      limited to crossing axes). `generate_spline_parallel_boolean_fixtures.py
      --check` writes 43 cases (15 fuses, 14 cuts, 14 commons; 39 solid, 4
      degenerate): spline walls against discs (`R125`, a quarter turn and
      a half turn `FLIP`, a new frame stored bit for bit), a stadium's arcs
      and lines (`TILT2` about `TILT`'s axis), a holed square, a capsule's
      spline against a disc; spline against spline (bulge and dome, blob
      and lens, the wave against its half-turned mirror, two domes about
      the tilted axis, a quartic hump against a cubic lens, two capsules);
      rounding offsets in `TILT`; and `degenerate` a disc and a half-turned
      dome each tangent to the dome along its apex generatrix. Checks: the
      divergence theorem, inclusion and exclusion and the area identity
      within 1.3e-40, a second slicing axis 1.3e-40, S9a.2's `SplinePair`
      on the 33 cases whose map is an exact turn or an offset (the tool's
      profile turned exactly into the object's frame) 2.4e-39 in `XY` and
      2.2e-17 in `TILT`, margins at least 0.06 outside the declared pairs
      (crossing sines 0.47, near misses 0.37, vertices 0.13), the declared
      pairs' crossing sines below 5.7e-21; Python 3.9 and 3.12 write the
      same files; `test_spline_parallel_boolean_reference.py`; the
      generator's check a CI group of its own (`spline-parallel`). The
      capture `occt-boolean-spline-parallel-preimplementation`
      (`compare_spline_parallel_boolean.py`, keyed on
      `solid/boolean/curved/spline_parallel.rs`; the probe `unsupported` on
      all 43, the 20 against arcs and the 23 against splines refused by
      `spline_pairs`): every result valid with the reference's solids, 31
      matches, 12 reviewed (BRepGProp's default integration up to 1.9e-7,
      the waves 8.6e-4; adaptively within 2.3e-9 in volume but for the
      waves; Green's theorem over OCCT's own faces and pcurves, a
      diagnostic build, within 1.1e-8 of the reference on all 39 results);
      four solids' counts change when unified. No correction to the
      decisions from the evidence. S9f.2a's kernel next.
  * **S9f.2a implemented** (`solid/boolean/curved/spline_parallel.rs`;
    `mod.rs`'s `spline_pairs` and `parallel_axes`, `meet.rs`'s sections and
    edge meetings, `spline_walls.rs`'s `foreign_param` and the arcs'
    implicit forms, `model.rs`'s `rational_proxy`, the operand on each
    `SplineSeg`; the validator's steep lines in `mass.rs`): spline walls
    against arc, circle and spline walls on exactly parallel axes as the
    refined decisions describe. Each Bézier arc keeps its implicit equation
    and inversion in its own frame, made once (`Implicit`: the Sylvester
    resultant and the first subresultant's `a`, `b` interpolated exactly on
    an integer grid; a coordinate of degree one inverts directly); a spline
    against a cylinder at the roots of `E` on the spline's arcs
    (`cylinder_roots`), two splines at the roots of `f(m(S))` on the
    object's segment's arcs, the other's parameter `-b / a` in the root's
    field (`crossings`, swapped when the tool's segment is given first, so
    the object's field always holds the crossing). Both refusal sites gave
    way (`spline_pairs` now refuses crossing axes only, relabelled S9f.2b
    and S9f; `spline_curved` names S9f.2b and S9f.3). All 43 fixtures as
    the reference (39 within the kernel's enclosures of volume, area and
    centre, each at most `1e-9` wide, with OCCT's unified counts on every
    solid; the 4 degenerate refused: a disc and a half-turned dome tangent
    to the dome along its apex generatrix); every history complete,
    results deterministic and moved rigidly
    (`tests/spline_parallel_booleans.rs`, 8 tests, 15 s at `opt-level` 2 and
    16 s in release: a mirrored dome refused,
    crossing axes refused with their labels, a quintic against a quartic a
    limit and a quartic against a quartic consistent, a disc through the
    wave's knot as S9a.2's one-frame Boolean; and the module's 3: the
    implicit forms and inversion, a crossing found from both sides in one
    field, coincident arcs refused). `compare_spline_parallel_boolean.py`
    31 matches and 12 reviewed, unchanged; every other Boolean comparison
    unchanged, `compare_split.py` 72/56 and `compare_brep.py --family
    spline` 10/3. Amendments to the decisions, from the implementation:
    (a) a crossing at a knot of either curve is kept only where the legs
    beside the knot lie on the other curve's two sides (`crosses_at_knot`,
    the knot's legs the arcs' end derivatives), at knots of both refused,
    as a touch at a knot; (b) coincident spline walls of two prisms meet
    first at their arches' ends, a vertex of one on the other's face
    (`Degenerate`), before any resultant vanishes; the module's test
    reaches the coincident arcs; (c) the validator's certified mass
    integrals on a cylinder or cone took a steep line (its `du` at most
    `2^-20` of its `dv`: the gap between a vertical edge's pcurve and a
    section's at a stadium arc's start in `TILT2`, their angles an ulp apart
    at `-7.4e-16`) by the expansion in its slope's powers, which left
    `dome_stadium_tilt2_common`'s centre enclosure `10^26` wide; such a line
    is now enclosed as `-du F` over its box (`rev_lines`, `BREP_VALIDATION.md`).
    The `boolean` fuzz target's spline variants against the stadium and
    the round hole or a spline in the frame turned about the axis evaluate
    through S9f.2a (`SPLINE_PARALLEL`, on), and in the tilted frame the byte
    at index 15, at or above 128, offsets a spline variant's tool by an
    amount that rounds (the tilted frame's own normal and `x`: normalized
    again they turned by an ulp, crossing axes). Replays with debug
    assertions: the corpus (1,430 inputs) and the 24 regressions, none
    failing, the slowest 6.3 s; 1,908 S9f.2a variants of every third
    corpus input (the object, the tool or both a spline prism against the
    stadium, the round hole or a spline, turned about the axis; and in the
    tilted frame at an offset that rounds), none failing, 5,606 of their
    operations evaluating (the rest refused as documented: the chained
    stage's spline results, S9f; sub-resolution results and contacts),
    median 0.31 s, the slowest 2.2 s (27 s under AddressSanitizer on this
    host, where the corpus's slowest input takes 29 s).
    Campaign: the boolean campaign at `089c5fbe` (600 s, a sampled
    replay) clean, 860 runs, the slowest input 44 s under
    AddressSanitizer on a host at load 10 (S9e.3b's agent working).
    DRAW survey: that of S9e.3b, S9f.2a and S9f.2b.1, below (no case of
    the group reaches S9f.2a's kernel).
  * **S9f.2b refined, before its code (2026-10-03).** Why it is refused
    today: `curved::spline_pairs` refuses a spline prism against a prism
    with arcs or circles whose axis crosses its own ("a spline prism
    against a prism with arcs on crossing axes (S9f.2b)"); behind it
    `meet::section`'s `spline_cyl` and `edge_surface`'s `wallcrv_cyl` and
    `conic_wall` need the exact parallel map between the frames
    (`spline_parallel::map2`, `None` on crossing axes: "a spline wall
    against a cylinder on crossing axes (S9f.2b)"), and no curve of the
    engine or of the topology is a spline wall's meeting with a quadric.
    Decisions.
    (1) *Domain.* Two prisms whose stored normals are not exactly parallel
    (`n_A x n_B != 0` in rationals), one holding a spline segment (with
    lines and arcs), the other arcs or circles (with lines) and no spline:
    its cylinder walls meet the spline walls in curves; its line walls are
    S9f.1's planes, arcs against arcs S9c's (and S9c.2's procedural
    meetings), unchanged. Spline walls against spline walls on crossing
    axes stay refused by design ("S9f refined"); spheres and cones are
    S9f.3's.
    (2) *The meeting.* Along the wall's ruling at the run parameter `tau`,
    `X = o + S_x(tau) x + S_y(tau) y + w n` (the spline prism's exact model,
    S9f.1's arcs), the cylinder's function in its own frame's exact rows
    (S9c's `Other`: `sum_i (g_i . X - e_i)^2 - r^2`, an elliptic cylinder in
    the world where the stored axes are not exactly orthonormal) is `F =
    A w^2 + 2 B(tau) w + C(tau)`, `A = sum_i (g_i . n)^2` a positive
    constant on crossing axes, `B` of degree `p` and `C` of degree `2 p` on
    each Bézier arc, so the branches are `w = (-B +- sqrt(D)) / A`, `D = B^2
    - A C` of degree `2 p`. Its roots are the turning points (the wall's
    ruling tangent to the cylinder, a graph over `tau` turning back, the
    curve regular there as a graph over the cylinder's angle).
    (3) *Sub-steps.* S9f.2b.1 (this step): every piece of the meeting a
    graph over `tau` with a fixed branch, between the segment's ends and
    the turning points: each branch kept over each maximal range of the run
    where `D > 0` (a turning point outside either face ends a range there,
    the pieces beside it outside too); a turning point inside both faces
    (loops, a meeting turning back inside the faces) is S9f.2b.2's,
    `OutOfDomain("a spline wall's meeting with a cylinder turning back
    inside the faces (S9f.2b.2)")`. S9f.2b.2: those loops, each piece about
    a turning point a graph over the cylinder's angle with the spline's
    parameter the root of the trace's equation in a verified window (as
    S9d.4b.2's `Toric`), switched at rational points, and the tower fields
    of (6).
    (4) *Representation.* Exact, in the curved engine: `Crv::WallMeet`,
    the spline segment on its model's frame, the other cylinder (`Other`),
    the branch and its range of `tau` (two numbers, each the segment's end
    or a root of `D`), placed by `tau` (`Pos::T`); a point is on it when its
    profile point lies on the segment within the range (S9f.1's `locate`),
    it lies on the cylinder and `A w + B(tau)` has the branch's sign
    (exactly); a piece's midpoint at a rational `tau` between its ends is
    in `Q(sqrt(D(tau)))`; its tangent `S' + w' n`, `w' = -(2 B' w + C') / (2
    (A w + B))`. Rounded, in the topology: `Curve3::WallMeet`, a new D13
    kind: the wall face's own stored surface (the spline's `BSplineSurface3`
    of degree `p` in `u`, one in `v`, as S9f.1 stores a wall), the other
    cylinder's stored frame and radius, the branch's sign and the `u` range
    (`start`, `sweep`); at `u` the surface's ruling `L(u) + v M(u)` (`L`
    and `M` the surface's two pole rows' combinations, `M` along the axis
    within rounding) meets the stored cylinder where `a v^2 + 2 b v + c =
    0`, and `v = (-b + s sqrt(b^2 - a c)) / a`, or `c / (-b - s sqrt(...))`
    where that cancels less (`Curve3::Meet`'s form). On its own wall its
    pcurve is its own `(u, v)` (a `Projection` that reads its parameters:
    deviation zero by definition, no inverse of the spline surface); on the
    cylinder a `Projection` by the cylinder's inverse, as `Meet`'s. Its
    validity (a nonrational wall of two pole rows, a sign, a range inside
    the surface's `u` domain, the discriminant positive at sampled
    fractions), points, jets, quadrature, tessellation bounds, rigid motion
    (the surface's poles and the other frame move), the writer's refusal
    (`Unwritable`, as `Meet`) and the importer's, history's,
    `curve_curve`'s and `curve_surface`'s refusals are each module's
    `Meet` handling extended.
    (5) *Certified bounds split at the knots.* The wall is C^(p-1) across
    its knots (the second derivative jumps for `p = 2`), so a jet of the
    curve comes from one knot span's polynomial: the integrals along it
    (the validator's Green integrals and mass terms on the wall and on the
    cylinder) are split at the knots' fractions, each piece integrated on
    its span's polynomial with the rounding slivers at the splits bounded
    (their width times both spans' integrands over the sliver, as S9f.1's
    crease pieces); a jet over a base across a knot is the union of both
    spans' jets only to order `p - 1`, else `None`; the tessellation's
    bounds take their pieces split at the knots.
    (6) *Degrees and fields.* Vertices: the meeting on the spline prism's
    caps, `F(h, tau) = 0` at a cap's height, degree `2 p` (`wallcrv_cyl`:
    the cap edge's curve against the cylinder); on the cylinder prism's
    caps, its cap plane's crease `w = h0 + h1 S_x + h2 S_y` in `F`, degree
    `2 p` (`conic_wall`: the cap's circle against the wall, its angle from
    the point); the cylinder prism's vertical edges against the wall,
    degree `p` (S9f.1's `line_wall`); the spline prism's vertical edges
    (knots, ends) against the cylinder, quadratic surds over `Q`. Every
    vertex lies in `Q(alpha)` of degree at most `2 p` (`alpha` the arc's own
    parameter, S9f.1's generators) or in `Q(sqrt(d))`; turning points are
    roots of `D`, degree `2 p`. A cylinder's cap plane holding the wall's
    axis direction exactly (perpendicular axes) meets the wall in
    generatrices, where its circle's points are in `Q(alpha)(sqrt(delta))`,
    `delta` in `Q(alpha)`, a tower the engine's numbers do not hold: such a
    cap circle meeting a spline wall is S9f.2b.2's (`OutOfDomain`); a rod
    through the wall with both caps outside is S9f.2b.1's.
    (7) *Degenerate.* The cylinder tangent to the wall (a root of `D` of
    multiplicity above one inside the run: an isolated point or a node of
    the meeting); a turning point at an interior knot (the cylinder tangent
    to the wall's ruling there) inside both faces' closures; a vertex's
    polynomial with a multiple root (an edge tangent to the other's face,
    S9f.1's and S9f.2a's labels); a turning point on a face's boundary; a
    vertex of one input on the other's face and crossings within the
    resolution, the engine's own rules. A meeting within rounding of
    tangency is not refused exactly (it validates or is `PrecisionLoss`);
    fixtures keep a margin.
    (8) *Stays refused*: S9f.2b.2's loops and towers, spline walls against
    spline walls on crossing axes, spheres and cones (S9f.3), tori, given
    results with spline walls, rational and periodic profile splines. The
    fields' degrees are at most `2 p <= 14` (S9f.2a's bound); no new limit.
    (9) *Evidence first*: `curved_boolean_reference.py`'s crossing-axes
    slicing and face sweeps taking a spline chord against a cylinder's
    (surd) chord: the slicing's events where a spline chord, a cylinder's
    chord and a height line are concurrent are the roots of the exact
    degree-`2 p` polynomials of (6) in the span's parameter, a cylinder
    wall's generatrices swept against a spline wall meet its crossings
    there too, and a spline wall's generatrices against the cylinder are
    S9c.1's surds over a polynomial family; `generate_spline_crossing_boolean_fixtures.py`
    across crossing angles (perpendicular rods through and beside spline
    walls, leaning and tilted cylinders whose meetings run cap to cap, a
    cylinder's cap inside the spline prism, holes, the capsule's arc and
    spline), S9f.2b.2's loops (refused by this step's kernel) and declared
    degenerate pairs (a rod touching the dome's apex, a turning point at
    the wave's knot); checked two ways (slicing and the divergence
    theorem), by inclusion and exclusion, the area identity, the faces'
    classes against their closed forms, perpendicular commons by the
    product of the profile's chord length and the disc's height chord
    (independent of the slicing's polygons), margins (vertices' crossing
    sines, turning points from caps and edges, near tangencies); Python
    3.9 and 3.12 the same files; a native capture before
    `curved/spline_crossing.rs` exists (`compare_spline_crossing_boolean.py`
    keyed on it); then the kernel in that file, its tests, and the boolean
    target's spline variants against the turned, leaning and tilted arcs
    (`SPLINE_CROSSING`).
    * **S9f.2b evidence (2026-10-03), before its code (`BOOLEAN.md`).**
      `curved_boolean_reference.py` takes arcs and circles against spline
      walls on crossing axes: the slicing's planes hold both axes, and an
      event where a spline chord, a cylinder's chord and a height line are
      concurrent is the walls' meeting at a cap's height, the real roots
      of an exact polynomial of degree `2 p` in the span's parameter
      (`wall_cylinder_events`: the cylinder's equation along the ruling,
      `w` fixed at the spline prism's cap or linear in `S(tau)` at the
      cylinder prism's); a cylinder wall's generatrices swept against a
      spline wall meet its implicit crossings at the same points and where
      their projection touches a span; the caps' sweep lines on spline cap
      edges and cylinder cap circles at the same roots; the turning points
      (`wall_cylinder_turns`, the roots of `B^2 - A C`) need no event of
      their own (the trace's tangencies). `generate_spline_crossing_boolean_fixtures.py
      --check` writes 34 cases (11 fuses, 11 cuts, 12 commons; 30 solid, 4
      degenerate; 28 of S9f.2b.1, 6 of S9f.2b.2): tilted rods along the
      bulge's and capsule's walls (the capsule's crossing its double knot),
      a leaning rod along the wave's top across its knot, a perpendicular
      rod covering the dome, a steep cylinder holding most of the blob as
      the object (`STEEP`, the normal (0, 5, 12) / 13), a tilted rod ending
      inside the blob (its cap circle on the wall), a stadium's arc and
      edges against the bulge, a tilted ring's hole around the dome, every
      meeting a graph over the spline's parameter inside the faces, its
      turning points at least 0.56 outside a face; loops (S9f.2b.2: a
      perpendicular rod through the bulge's wall, a tilted rod through the
      lens's); and `degenerate` a rod touching the dome's apex (a point
      contact) and one whose meeting with `knot`'s wall turns back at its
      knot of multiplicity two. Checks: the divergence theorem, inclusion
      and exclusion and the area identity within 9.4e-41, every face's
      classes 4.5e-41, the four perpendicular pairs' commons as `int L(y)
      W(y) dy` (the profile's chord along `x` times the disc's height
      chord) 7.1e-44; margins outside the declared pairs at least 0.011
      (vertical edges piercing curved walls, vertices 0.037 from the
      other's faces, creases 0.073 from the caps, meetings crossing caps
      at 0.11, turning points 0.56 outside a face, the loops' 1.5 inside
      both); the same files under Python 3.9 and 3.12;
      `test_spline_crossing_boolean_reference.py` (the meeting at a cap
      and its turning points against closed forms, the tower field's
      points, the product of chords, the degenerate margins); the
      generator's check a CI group of its own (`spline-crossing`). The
      capture `occt-boolean-spline-crossing-preimplementation`
      (`compare_spline_crossing_boolean.py`, keyed on
      `solid/boolean/curved/spline_crossing.rs`; the probe `unsupported` on
      all 34, refused by `spline_pairs`): every result valid with the
      reference's solids, 13 matches, 21 reviewed (BRepGProp's default
      integration up to 1.5e-6, `knot_turn_cut` 1.3e-5, the wave 2.9e-4;
      a diagnostic build's adaptive BRepGProp and Green's theorem over
      OCCT's own faces and pcurves, the better within 4.0e-9 in volume,
      1.4e-9 in area and 8.6e-9 in the centre on all 21); five solids'
      counts change when unified. No correction to the decisions from the
      evidence. S9f.2b.1's kernel next.
  * **S9f.2b.1 implemented** (`solid/boolean/curved/spline_crossing.rs`;
    `mod.rs`'s `spline_pairs`, `meet.rs`'s sections and edge meetings,
    `spline_parallel.rs`'s cap edges and conics on crossing axes, `graph.rs`'s
    placements, `assemble.rs`'s curves and the own wall's pcurves; the
    topology's `Curve3::WallMeet` with `topology/validate/wall_meet.rs` and
    its handling in `projection.rs`, `quadrature.rs`, `spline_flux.rs`,
    `mass.rs`, `spline_taylor.rs`, the tessellation's bounds, the writer,
    the reader, history, `curve_curve`, `curve_surface` and `chain.rs`):
    spline walls against cylinder walls on crossing axes as the refined
    decisions describe. `Crv::WallMeet` holds the segment, the cylinder's
    exact rows, the branch and its range of `tau`; its sections come from
    `D`'s exact roots on the segment's arcs, each turning point classified
    by its point on the ruling (`w = -B / A`) against both faces; vertices
    from S9f.1's `roots_of` on degree-`2 p` polynomials (`wallcrv_cyl`,
    `conic_wall`). All 34 fixtures as the reference (the 24 results of
    S9f.2b.1 within the kernel's enclosures of volume, area and centre, each
    at most `1e-9` wide; the 4 degenerate refused; S9f.2b.2's 6
    `OutOfDomain` naming S9f.2b.2); every history complete, results
    deterministic and moved rigidly (`tests/spline_crossing_booleans.rs`, 7
    tests, 42 s in release and 41 s at `opt-level` 2 with debug
    assertions: also a cap circle on the dome's wall along its axis refused
    as S9f.2b.2's and the same rod through the whole prism evaluated,
    crossing spline walls refused, and two found by the fuzz variants
    below; the module's 4: the curve on both surfaces with jets enclosing
    its points and slopes, jets across a knot to the wall's continuity
    and `None` past it, the pieces' exact splits, the own wall's pcurve).
    `compare_spline_crossing_boolean.py` 13 matches and 21 reviewed (the
    capsule's three reviews now name its counts too: OCCT splits the
    meeting's ellipse arcs on the capsule's caps at the rod's seam
    generatrix, `x = 0.2`, at `(0.2, -0.65, 0)` and `(0.2, 0.85, 2)`, and in
    the fuse stops the rod's seam at the hole, which `occt_counts` does not
    synthesize for a hole that does not wind); every other comparison
    unchanged, `compare_split.py` 72/56, `compare_brep.py --family spline`
    10/3, `compare_brep_io.py` 6835/7. Amendments to the decisions, from
    the implementation: (a) a jet over a base across a knot is the union of
    both spans' jets up to order `k + 1`, `k` the wall's least continuity
    (`p - 1` for R4's walls), not `p - 1`: Taylor's remainder holds with
    `f^(k)` absolutely continuous, its next derivative between the spans'
    values; integrals are split at the knots' fractions exactly (rationals
    of the stored `start` and `sweep`), so no sliver is left to bound; (b)
    `a`, `b`, `c` and `d` are exact Bernstein polynomials in a span's `ū`,
    made once per wall and cylinder by the product rule, enclosed at a point
    within a few units in the last place of their coefficients: formed from
    enclosed factors (the foot and the direction each from two rows), `d`
    was a hundred times wider, and a meeting's end `6.7e-7` in `u` short of
    a turning point (a lens wall under a leaning holed slab, a fuzz variant)
    left its `v` `2.2e-10` wide past the wall's cap, its chord's ends in no
    patch (`uncertified_shell_orientation`), now `3.5e-12`; over a range
    (remainder boxes) the factors' polynomials are combined instead; (c)
    the quadrature halves a wall piece toward a turning point (at most 40
    times), without a try where `d`'s binary64 values over the piece differ
    by more than a factor of two (a failing sweep spends its whole budget
    first: twenty failing tries took minutes), so that end takes some
    twenty halvings, one sweep each; (d) a closing chord on a spline
    surface may pass the domain's edges (`chord_patch`: the boundary patch's
    polynomial is the only one there), still within `2^-40` of a patch's
    width at its interior sides; (e) a turning point outside the faces but
    within the resolution of both is `Degenerate` with the label of one on
    a face's boundary (R4's knot under a leaning stadium: its ruling at
    `u = 1/2` tangent to the arc's cylinder on the top cap but for the
    frame's rounding, a fuzz variant whose meeting's end left its loop
    winding undecided); (f) the binary64 tier's conversions of a span's
    exact numbers are made once (most of an evaluation otherwise). The
    `boolean` fuzz target decodes a spline prism against the stadium or
    the round hole leaning, tilted or on its side through S9f.2b.1
    (`SPLINE_CROSSING`, on; before, refused as S9f.2b's), spline against
    spline there still refused. Replays with debug assertions: the corpus
    (1,430 inputs) and the 24 regressions, none failing, the slowest
    11.4 s on a host at load 14; 1,431 S9f.2b variants of every third
    corpus input (the object or the tool a spline prism against the stadium
    or the round hole, leaning, tilted or on its side), none failing, 3,232
    of their operations evaluating (the rest refused as documented: the
    chained stage's spline results, S9f; S9f.2b.2's cap circles along the
    wall's axis and turning points inside the faces; tangencies and
    sub-resolution contacts, 8 turning back on a face's boundary), median
    0.34 s, the slowest 5.7 s (47 s under AddressSanitizer on this host at
    load 9 to 14, where the corpus's two slowest inputs took 47 and 52 s).
    Campaign: the boolean campaign at `6c77655a` (600 s, a sampled
    replay) clean, 902 runs, the slowest input 50 s under
    AddressSanitizer at load 6 to 8 (`1b405929`, an existing corpus
    input, near the 60 s limit). DRAW survey: that of S9e.3b, S9f.2a and
    S9f.2b.1, below (no case of the group reaches S9f.2b.1's kernel).
  * **DRAW survey of S9e.3b, S9f.2a and S9f.2b.1 (2026-10-03,
    `UPSTREAM_TESTS.md`).** The 1,802 cases of the Boolean group run again
    on both backends at `b0b9adc6` (the public dataset, 120 seconds a case,
    four at once): every status and every refusal's reason is the last
    survey's (`507b8054`) on both backends. Rust evaluates 987, all
    registered, 592 refused, 223 unsupported on both, as before. No case
    evaluates newly, none fails or times out, the sentinels are refused as
    before. The three kernels reach no case: the group's spline solids come
    from `nurbsconvert` (96 cases), which the adapter does not run, so no
    spline wall meets another prism's wall (S9f.2a, S9f.2b.1); of its
    results given to another Boolean only `bcut_simple/G9` and `H3` reach
    the kernel, refused as a tangency between the inputs, and no given
    result the adapter makes holds a meeting of two curved faces met by
    the partner (S9e.3b, as its run of the chained cases found). The
    volume audit of the 987: both backends' values the last audit's bit
    for bit (Rust) and to the digit (native DRAW), the coaxial tori `ZL2`
    to `ZL5`, the spheres sectioned through their poles `ZI4` to `ZI7` and
    every case with cylinder or cone faces among them, so S9e.3b's
    validator for holes in torus bands wound in v, its bounds taking torus,
    cylinder and cone faces' bulges and its pole vertices, and S9f.2a's
    mass integrals taking a steep line by its box, move no registered
    case's values and no worker at an earlier commit was needed; the same
    35 disagreements, native off in each. `bopfuse_simple/ZP6` and the
    `gdml_public` tori refused as before. A full contract run holds (the
    slowest Boolean case 12.8 seconds on a loaded machine, 13.5 in the
    last survey); the ledger does not change; no kernel change.
  * **S9e.4 refined, before its code (2026-10-03).** Why it is refused
    today. A body without a construction has no `Solid` at all: the
    `.brep` converter (`occt_brep::import`) and STEP's reader give a cell
    `Topology` with the body's resolution (the largest OCCT tolerance), and
    nothing builds a solid from it; DRAW's adapter keeps a restored solid as
    a `Shape::Body`, which `boolean_argument` refuses ("an argument other
    than a solid the adapter made"), and S9b.2's stored model takes only a
    Boolean's results and a plane's pieces. What such a body stores (OCCT's
    `.brep`, `TopTools_ShapeSet` and `GeomTools`): surfaces and curves with
    17 significant digits (OCCT's binary64 values exactly), vertices,
    tolerances and edge ranges with 15, so a stored vertex lies off its
    faces' stored surfaces by up to about `5e-16` of its size (a tilted
    box's corner `1.4` beside its edge's line ending at
    `1.4000000000000004`); the converter builds every frame again with
    `Frame3::new` (the stored axes the kernel's normalization of OCCT's,
    within an ulp, with the platform's `hypot`). The S9e text's plan (every
    edge the meeting of its faces' stored surfaces, every vertex their
    common point within the resolution, faces' regions and the solid's
    membership decided exactly on general faces) needs an exact point
    membership and face regions for faces of any shape on curved surfaces,
    which nothing in S9 has (each model's `member` and `in_face` are its
    construction's: a profile's, a sphere's latitudes, a cone's heights).
    Decisions. (1) *Sub-steps.* **S9e.4a** (this step): an imported solid
    whose stored topology is one of the kernel's constructions (a prism of
    lines, arcs and circles; a sphere, a cap or a zone; a cone or a
    frustum; a whole torus) is that construction, read off its stored
    surfaces, and decided on the construction's exact model; **S9e.4b**:
    every other imported body (a general body on its stored surfaces, the
    S9e text's plan, with an exact membership of its own), and a plane's
    piece (`Clipped`, `Half`) against curved faces, which S9e.2 deferred
    here (its refusals keep naming S9e.4). (2) *The constructor.*
    `Solid::imported_with(operation, topology, resolution)` (and `_in`):
    the topology's entities renamed under the operation as
    `Topology::from_parts` names them (kind `External`, role `External`,
    each slot its ordinal) but with the operation's id, so two imported
    bodies have distinct ids; the history every entity `Generated` with no
    parents, as a primitive's. It recognizes the construction (3), builds
    it with ids of its own (an operation derived from the import's), and
    matches it to the stored topology (4); a body recognized as none is
    `OutOfDomain("an imported solid other than a prism, a sphere, a cone or
    a torus (S9e.4b)")`, spline faces or edges `OutOfDomain(... (S9f))`. (3)
    *Recognition*, each number read from the stored data once: a *prism*:
    the first pair (in stored face order) of plane faces whose outward
    normals are opposite, every other face a plane parallel to their normal
    or a cylinder whose axis is (candidates within `1e-9` in the sine; the
    match decides), the first of the pair the bottom cap; the frame the
    bottom cap's stored plane frame (turned `(x, -y, -n)`, exactly, where
    its normal leaves the material: the prism's normal points into it); the
    profile the bottom cap's loops: each vertex's local coordinates in the
    frame's exact affine map, rounded once to binary64; a stored line a
    line; a stored arc its circle (its centre's local coordinates rounded
    once, the stored radius, counter-clockwise where its sweep, its frame's
    normal against the prism's and its fin's sense agree); a ring a circle;
    the outer boundary the loop of the greatest area, the others holes; the
    heights `0` (the bottom cap's plane is the frame's own, exactly) and the
    top cap's stored origin's local height rounded once; `Solid::extrude`
    of them. A *sphere, cap or zone*: one spherical face and at most two
    planar faces, each bounded by one ring of it, normal to its stored axis;
    the stored frame and radius exactly, each disc's latitude `asin(h / R)`
    of its ring's centre's local height `h` (a missing end a pole, `-+pi /
    2`). A *cone or frustum*: one conical face and one or two discs normal
    to its axis, the other end the apex vertex; the stored frame with the
    axes kept bit for bit and the origin moved along the axis to the lower
    end's height where that is not zero, each end's ring's stored radius
    (`0` at the apex) and the upper end's height above the lower. A *whole
    torus*: one toroidal face without loops; the stored frame and radii
    exactly, the whole tube and turn. A torus's v-segment or wedge, a
    cavity or several shells, any other face: S9e.4b's. (4)
    *Verification.* The construction's topology is matched to the stored
    one by S9e.2's geometric match (`matched.rs`): every vertex the one
    stored vertex within the resolution, one to one; every edge the stored
    edge between the matched ends through its points at a quarter, a half
    and three quarters (a ring by its distance from the stored ring); every
    face the stored face bounded by the matched edges on a surface of the
    same kind. With the converter's validation (every stored edge on its
    faces' stored surfaces within the resolution, `pcurve_off_edge`), every
    stored surface then lies within twice the resolution of the
    construction's along every edge, and the construction's surface is the
    stored one exactly where it was read (a sphere's, a torus's, a cone's
    frame, a cap's plane). A body unmatched is S9e.4b's (`OutOfDomain`).
    (5) *The exact model and the stored edges.* The body is decided on its
    construction's exact model: S9c's prism model on the stored cap's frame
    (its affine axes in rationals), S9d's sphere, cone and torus models on
    the stored frames. No stored edge is trusted as an exact curve: every
    edge is the construction's (a prism's lines and arcs between its
    rounded profile points, its verticals, a sphere's or a cone's rings at
    their heights), and the stored curves only name and verify it. So a
    profile's line tangent to its arc (a smooth join) and a periodic face
    split at a seam into faces (two arcs of one circle, two cylinder faces)
    are the construction's profile data, decided exactly as S9c decides
    them, not refused (in S9e.4b, where an edge is its faces' meeting, they
    stay refused). (6) *In a Boolean.* `polyhedra::build` takes an imported
    input as its construction (with the construction's ids), so every
    pair's engine and every rule are the construction's; the result's
    components' plans name the imported body's entities through the match
    (each construction entity the stored entity matched to it), so the
    history is over the imported body's ids, and the result keeps the
    imported solid as its input (classification, rigid motion). A result
    of an imported solid given to another Boolean (S9e.1 to S9e.3) re-runs
    its first arrangement with the same construction. The construction's
    ids must be apart from the partner's (`InvalidLabel` otherwise, never
    met but by a forged operation id). (7) *Refused.* A prism's arc whose
    ends, rounded into the cap's frame, are off its circle exactly (S9c's
    requirement: a frame whose rotation rounds the profile's points, a
    profile turned by 30 degrees) `OutOfDomain("an imported prism's arc
    whose ends round off its circle in its cap's frame (S9e.4b)")`; S9's
    rules unchanged (`Degenerate`: an imported cylinder tangent to the
    partner's plane, a face of the partner touching an imported sphere at a
    point). (8) *Its other queries.* Classification is the construction's
    (within the resolution); mass properties the stored topology's
    certified enclosure (the body's own measure); bounds the construction's
    widened by the stored edges'; a rigid motion moves the stored topology
    (its stored geometry, as S9b's results) and the construction with it
    (rebuilt in the moved frame), the match kept (both keep their ids). (9)
    *The reader.* OCCT 8.1's `BRepTools::Write` writes format version 3
    under the copyright line `(c) Open Cascade`, which the reader refuses
    (its header check expects `Matra-Datavision`; 27 files of the public
    dataset carry the newer line): the fixtures are written in version 1
    (`TopTools_FormatVersion_VERSION_1`, no triangulations); the header
    check is left to the import track, since widening it changes which
    dataset files restore (a survey of its own). (10) *Fuzzing.* The
    `boolean` target's object written by the kernel's writer and read back
    (`occt_brep::write`, `read`, `import`), imported and given the same
    Boolean, its volumes those of the object's own result within the
    resolution (`IMPORTED`, a switch). (11) *Evidence first.* Bodies
    written by OCCT (`occt_boolean_oracle.cpp`'s `write` blocks:
    `BRepPrimAPI_MakeBox`, `MakeCylinder`, `MakeSphere`, `MakeCone`,
    `MakeTorus` and `MakePrism` of a profile) under
    `rust/fixtures/imported/`, generated, never from `data/`; cases whose
    inputs may be a `brep PATH` row (`identity_reference`'s `Case.brep`);
    the reference the constructions OCCT was given, through S9e.3a's
    chained reference (`chained_curved_boolean_reference.py`, pairs and one
    chain), S9d.1's sphere reference for a zone, coaxial closed forms,
    Monte Carlo and the pair identities, with each file's stored vertices
    checked within `1e-12` of the case's size of the construction's
    (`generate_imported_boolean_fixtures.py --check`): boxes (one in the
    `TILT` frame), cylinders (one along `x`), a prism with tangent arcs, a
    cylinder of two half faces, a plate with a hole, a sphere, a zone, a
    frustum, a cone with its apex and a torus, against prisms of lines and
    arcs, a sphere and another imported body, as object and as tool, and
    one chain; declared `degenerate` an imported cylinder tangent to a
    box's wall, a box touching an imported sphere and a cylinder tangent to
    an imported box's wall; declared `unsupported` the prism with arcs
    turned by 30 degrees. A native capture before `solid/imported.rs`
    exists (`compare_imported_boolean.py` keyed on it, the kernel's probe
    reading every file and reporting `unsupported`),
    `test_imported_boolean_reference.py`, the generator's check a CI group
    (`imported`); then the kernel, its tests (`tests/imported_booleans.rs`:
    enclosures within `1e-9` of the reference, degenerate cases refused,
    histories over the imported bodies' ids, determinism, rigid motion),
    the DRAW adapter's restored solids as Boolean arguments, the fuzz
    switch, the DRAW survey and a campaign.
  * **S9e.4 evidence (2026-10-03).** The bodies are OCCT's own output:
    `occt_boolean_oracle.cpp` takes a block `write NAME PATH` (one solid's
    rows, `BRepTools::Write` in format version 1 without triangulations)
    and new rows `box` (`BRepPrimAPI_MakeBox(gp_Ax2, DX, DY, DZ)`),
    `cylinder` (`BRepPrimAPI_MakeCylinder(gp_Ax2, R, H)`) and `brep PATH`
    (`BRepTools::Read` of the file's one solid, under the directory
    `OCCT_BOOLEAN_FIXTURES` names); `compare_imported_boolean.py
    --write-bodies` writes the 13 bodies of `boolean-imported-bodies.txt`
    to `rust/fixtures/imported/` (generated, nothing from `data/`): a box
    (`MakeBox`), a box in the `TILT` frame, a cylinder (`MakeCylinder`), a
    cylinder along `x`, a prism of lines and an arc tangent to them, a
    cylinder of two half faces, a plate with a hole (`MakePrism`, whose top
    shares its bottom's records under a translation), a sphere, a
    hemisphere (a cap), a frustum, a cone with its apex, a torus and the
    tangent profile turned by 30 degrees. The case protocol takes an
    imported input as one row `brep PATH` (`identity_reference.Case.brep`,
    `encode_case` and `native_case`; `tests/support/boolean_protocol.rs`'s
    `input` reads and converts the file, and since no constructor takes it
    reports the case `OutOfDomain`). `generate_imported_boolean_fixtures.py
    --check`: 69 cases of 23 groups (57 solid, 9 declared `degenerate`, 3
    `unsupported`; 21 of class `prism`, 15 `cylinder`, 12 `sphere`, 9
    `cone`, 6 `torus`, 3 `both`, 3 `chain`), each imported body against
    prisms of lines and arcs (`TILT` slabs, coaxial and crossing rods,
    boxes), a sphere or another imported body, as object and as tool, and
    one chain (the imported box less a rod, then with a `TILT` slab);
    `degenerate` an imported cylinder whose generatrix a box's wall
    touches, a box touching the imported sphere's top point and a rod
    touching the imported box's wall; `unsupported` the turned profile
    against a box. The reference is the constructions OCCT was given:
    S9e.3a's chained reference for 21 groups, S9d.1's sphere reference
    (`generate_sphere_boolean_fixtures.evaluate`, its own checks) for the
    hemisphere against a `TILT` box, the coaxial sections
    (`generate_given_curved_boolean_fixtures.Coaxial`) for the plate with
    a hole (with a Monte-Carlo estimate of its own). Checks, relative to
    the case's size: closed forms (coaxial sections for five groups, the
    hemispheres for two) within 4.3e-41; the two families 8.6e-32; each
    solid's closed form 2.8e-40; the pair identities 2.3e-41 and the area
    identity 1.9e-40; the sphere reference's checks within 1.4e-40; Monte
    Carlo 2.4 standard errors; solid counts by rays at two resolutions;
    every meeting's sine at least 0.45 and events at least 5.3e-5 of their
    range apart outside the declared cases (an event found twice within
    1e-25 counted once: a tangent profile's joint is both a vertex and its
    arc's extreme along a chord). Every file read independently
    (`stored_records`: `brep_io_reference`'s records and locations, each
    face and vertex instance as placed): its faces' kinds the
    construction's and every vertex within 3.7e-32 of the size on the
    construction's surfaces (OCCT's 15 digits). Correction from the
    evidence: `torus_curved_boolean_reference.Prism.closed` added a
    polygon's moment terms (`(x_i + x_j) cr / 6`) to the arcs' Green's
    forms (`int x^2 dy / 2`): the two agree only summed over segments
    alone, so a profile whose arc's centre lies off the origin's lines
    (the tangent profile) had its moments off by `cr / 6` per segment; the
    segments now take the arcs' forms (its fixtures' profiles, stadiums
    about the origin, were not affected: a check, the fixtures unchanged;
    every generator importing it is checked again). Fixture corrections
    from the reference's runs: the solid counts of three groups (the slab
    leaves the tangent profile's arc whole; the rod takes the frustum's
    narrow end; a slab below the ball replaced by one across its middle),
    a sphere touching the cone's apex moved to hold it inside (the
    apex at its centre, the meeting at right angles), a tangent box's
    common declared empty. `test_imported_boolean_reference.py` checks two
    boxes by their grid cells exactly, the coaxial sections against the
    chained reference, the hemispheres against a quadrature of their
    sections, the files' reader on the hemisphere's and the tangent
    profile's files, and the case list and its protocol rows. The
    generator's check is a CI group of its own (`imported`); Python 3.9
    and 3.12 write the same files. The capture
    `occt-boolean-imported-preimplementation` (`compare_imported_boolean.py`,
    keyed on `solid/imported.rs`; the kernel's probe `unsupported` on all
    69, every file read and converted): every result valid, 61 matching
    (volumes within 5.1e-9, areas 3.4e-9, centres 4.5e-9 of the size), 8
    reviewed: the crossing cylinders and the torus with the rod through its
    tube, BRepGProp's default integration on faces bounded by approximated
    sections (up to 2.8e-5; adaptively within 1.7e-9, unchanged at
    1e-12), and the declared tangent fuses, which OCCT keeps as two solids
    touching along a line where the reference's rays count one; 14
    results' counts change when unified. S9e.4a's kernel next.
  * **S9e.4a implemented** (`solid/imported.rs`, `Construction::Imported`,
    `polyhedra.rs`'s `substituted`, `given.rs`'s construction,
    `model.rs`'s profile test, `Frame3::flipped`; the DRAW adapter): an
    imported solid decided on the construction its stored surfaces give,
    as the refined decisions describe. All 69 fixtures as declared (57
    within the kernel's enclosures, each at most `1e-9` wide; the 9
    declared degenerate refused as tangencies; the 3 turned profiles
    `OutOfDomain` as S9e.4b's, their arcs' ends rounding off their
    circles), every history complete over the imported bodies' stored ids
    and none naming a construction's own, results deterministic and moved
    rigidly, both inputs translated keeping the reference's volumes, every
    body its construction (its closed-form volume within `1e-9`, points
    classified, its stored vertices on its boundary, two imports' ids
    apart), STEP solids (a box, cylinders one tilted, a sphere, a
    hemisphere, a cone, a frustum, a torus, an L prism, a plate with a
    hole) imported and cut consistently, a cavity and a spline prism
    refused, the stored frames the reference's bit for bit
    (`tests/imported_booleans.rs`, 9.5 s at `opt-level` 2).
    `compare_imported_boolean.py` 54 matches and 15 reviewed (the native
    measures as captured; entity counts: the kernel's whole periodic faces
    and exact meeting pieces against OCCT's seams and split approximations,
    the imported seam-split cylinder's two faces against OCCT's unified
    one), every enclosure within the reference with the comparison's
    `1e-12` slack; every other comparison unchanged (boolean 45/0, its
    splines 33/13, polyhedral 43/2, curved 42/2, procedural 4/24, turned
    2/13, capped 0/18, sphere 30/0, spheres 12/21, cone 25/5, cones 21/19,
    torus 11/24, torus segment 15/14, torus curved 15/29, spheres turned
    0/18, cones' loops 5/26, torus parts 16/37, chained 24/6, given 36/0,
    given curved 25/23, spline any 22/16, spline parallel 31/12, spline
    crossing 13/21, given met 8/42; split 72/56, brep spline 10/3, brep_io
    6,835/7, step 23/6 on STEP-b's SDK). Amendments and corrections, from
    the implementation: (a) the history checker requires a `Generated`
    relation's parents, so an import's entity is generated from a label of
    its stored slot (its kind above its ordinal), not from none; (b) a
    construction the stored data cannot build (a profile whose rounded
    points touch, a degenerate height) is none of the kernel's, S9e.4b's
    `OutOfDomain` (a dataset prism of `bcut_complex/I6` raised the
    profile's self-intersection before); (c) the curved engine's profile
    test took a point on an arc's chord as on the arc's side of it: two
    arcs of one circle (a seam-split cylinder) share their chord run either
    way, a point on it counted in both circular segments, and the box's cap
    piece inside the circle was kept, the fuse left open ("an open Boolean
    of arcs in any position") wherever the curved engine met such a prism
    in frames of equal axes (latent since S9c.1, a kernel prism of two half
    circles against a box with an offset that rounds reproducing it): a
    point on a chord is displaced alike for every chord, `(eta, eps)` as
    the ray's half-open crossings take it (`tests/curved_booleans.rs`). The
    `boolean` fuzz target's object is written by the kernel's writer, read
    back, imported and given the chosen operation again (`IMPORTED`, on):
    replaying the corpus (1,430 inputs) and the 24 regressions with debug
    assertions, no failure, the slowest 20 s on a host at load 11 to 37
    (11.8 s instrumented); of the imported operations 645 evaluate, 237 are
    refused as documented (42 arcs off their circles once rounded,
    S9e.4b's; 195 degeneracies of the converter's frames normalized again:
    faces within the resolution of one plane, thin faces, planes within
    rounding of a cylinder's or a spline wall's direction, tangencies) and
    476 spline objects are not recognized (S9f). With `GIVEN_MET` on the
    replay is clean too (the slowest 31 s); with all four off switches on at
    once one corpus input (`6fab9d41`) fails `vertex_off_curve` before the
    imported stage, each switch alone clean and S9e.4's engine change
    reverted the same: an interaction of the off switches, not S9e.4's, open.
    The `brep_io` corpus (1,469 inputs) and its 2 regressions replay clean
    (the slowest 0.6 s). Pending: the DRAW survey, the campaign.
  * **DRAW survey of S9e.4a (2026-10-03, `UPSTREAM_TESTS.md`).** The
    cases restoring a shape and giving it to a Boolean, never surveyed
    before (the adapter refused restored arguments; the Boolean surveys ran
    the self-contained cases): 1,814 of every group, none registered, run on
    the Rust adapter (the public dataset, 120 seconds a case) and on native
    DRAW where the import reaches them. 561 load private data, 598 give a
    Boolean a restored shape other than one solid, 123 stop at constructs
    the reader does not represent, 276 at the validator's rejection of a
    restored shape, 8 time out at their first restores, 77 at other
    commands the adapter does not read; the import reaches 171: 84 bodies
    none of the kernel's constructions and 31 prisms whose arcs round off
    their circles (S9e.4b), 3 with spline faces (S9f), 37 refused by S9's
    rules (tangencies, faces within the resolution of one plane, thin
    faces), and 16 evaluate on both backends with every check (boxes,
    hexahedra, wedges, prisms with cylindrical walls and cylinders fused,
    cut and intersected, `bfuse_complex/N3` a result of imported solids
    fused again, `bugs/modalg_6/bug21427`). Their volume audit (`vprops`
    and `sprops` before each `checkprops`): Rust's values native DRAW's to
    its printed digits in all 16. Registered (1,089 cases; the contract
    holds on both backends within 30 seconds, Rust 0.2 to 5.6 s on a loaded
    host); the ledger records `F5`'s and `Q2`'s `checknbshapes` confirmed
    natively (mapped-and-verified 2 to 4, lost 12,844 to 12,842). The
    survey found the refusal of a profile its rounded points make touch
    (`bcut_complex/I6`, now S9e.4b's: amendment (b) above). Campaign: the
    boolean campaign at `7199e06a`, S9e.4a's and S9f.2b.2's together (600 s,
    a sampled replay, `IMPORTED` and `SPLINE_CROSSING` on) clean, 845 runs,
    the slowest input 47 s under AddressSanitizer at load 12 to 17
    (`17e131e3`, an existing corpus input, the torus against prisms).
  * **S9e.4b refined, before its code (2026-10-03).** Why each class is
    refused today, from the DRAW survey of S9e.4a (the 114 restore cases
    refused as S9e.4b's, the solids of their 86 files listed by the
    import's reader: faces, edges, shells) and the fuzz target's imported stage.
    (a) *Prisms whose arcs round off their circles* (31 cases, the fuzz
    replay's 42): `polyhedra::substituted` refuses them ("an imported
    prism's arc whose ends round off its circle in its cap's frame
    (S9e.4b)") because S9c's model (`model::Prism::new`) takes each arc
    between its profile's binary64 points and needs both on the circle
    `|p - c| = r` exactly, and an imported profile's points are the stored
    vertices' local coordinates rounded once in the cap's frame, its arcs'
    centres likewise and their radii stored: a vertex of a cap turned about
    or tilted from the world's axes (`dee_turn`, every `CTO9xx` part whose
    frame is not the world's), or at an angle whose cosine is irrational
    (OCCT's 15 digits), lies off the circle by the rounding of its
    coordinates (up to 7e-14 of the radius in the surveyed `CTO`, `cts`
    and `pro` files, 5e-11 in `bug27948_a.brep`, whose stored vertices are
    OCCT's approximations within its tolerance). Their joints: a line and
    an arc (tangent or crossing) or two arcs of one circle in 25 cases;
    two arcs of different circles meeting (fillet chains, lenses) in 6
    (`bcut_complex/E8`, `P4`, `bfuse_complex/D5`, `E1`,
    `bugs/modalg_2/bug4993_1`, `_2`). (b) *Bodies none of the kernel's
    constructions* (83 cases, `imported::recognize`'s "an imported solid
    other than a prism, a sphere, a cone or a torus"): sphere pieces cut by
    planes (`so1` to `so7`: a sphere with one or three plane faces; 38
    cases), a hollow sphere of two spherical shells
    (`case_8_solid_repaired.brep`; 23 cases, each a Boolean with a wire,
    refused as such anyway), polyhedra other than prisms (hexahedra with
    slanted faces, a pyramid's frustum, a twisted hexahedron:
    `buc60803a`, `b`, `OCC578_w1`, `w2`, `pro9481b`, `CTO900_pro12559a`;
    7), prisms with walls along two directions (`CTO900_fra50089-part`,
    `CTO900_pro9476-part`, `CTO904_cts20370-part`, `bug28773_2`; 4), and
    cones, cylinders, spheres and tori together or a cone cut by planes
    (`OCC485a`, `cts21128c`, `d`, `shading_132`; 11). (c) *A plane's
    piece* (`Clipped`, `Half`) against curved faces, which S9e.2 deferred
    here (`curved::build`'s and `polyhedra.rs`'s refusals naming S9e.4).
    Decisions. (1) *Sub-steps*, ordered by what the survey's cases and the
    fuzz target need and by the machinery each takes: **S9e.4b.1** (this
    step): the prisms of (a) whose joints are a line and an arc or two arcs
    of one circle (25 cases, the fuzz replay's 42), on S9e.4a's
    construction with each arc's ends taken onto its circle (2);
    **S9e.4b.2**: polyhedra other than prisms (7 cases), every face its
    stored plane (the frame's origin and normal in rationals, an exact
    plane), every vertex the common point of its faces' planes (three or
    more planes: a rational point, unique and within the resolution of the
    stored one, else `Degenerate`), every edge the meeting of its two
    faces' planes between its vertices, decided as S9b.2's stored model
    decides a polyhedral result (faces' regions by their exact polygons,
    membership by parity); **S9e.4b.3**: a plane's pieces of a sphere, a
    cylinder or a cone (one curved face and planes: `so1` to `so7`,
    `shading_132`; 39 cases) and S9e.2's deferred `Clipped` and `Half`
    against curved faces: every edge a plane's meeting with the curved
    surface (a circle, an ellipse, a cone's conic) or with another plane,
    every vertex their common point (a quadratic surd, unique within the
    resolution of the stored one), each face's region on its stored
    surface bounded by those exact curves, membership by the curved
    surface's side and the planes' half-spaces (parity along a ray where
    the body is not convex); **S9e.4b.4**: the S9e text's plan in full:
    profiles whose joints are two circles' common points (a profile
    segment's ends in a quadratic field, which S9c's `Seg` with rational
    ends does not hold; 6 cases), prisms with walls of two directions (4)
    and bodies of several curved surfaces (11), and a profile its rounded
    points make touch (`bcut_complex/I6`). Not S9e.4b's: the hollow
    sphere's 23 Booleans with wires (no solid argument), spline bodies
    (S9f). The reader's header check (OCCT 8.1's version-3 `.brep` under
    `(c) Open Cascade`, 27 dataset files) stays open on the import track:
    widening it changes which dataset files restore, a survey of its own,
    and none of the 114 cases reads such a file. (2) *S9e.4b.1's
    representation.* The construction is S9e.4a's (the profile's points
    the stored vertices' local coordinates rounded once in the bottom cap's
    stored frame, each arc's centre likewise, its radius stored), its
    topology and the match to the stored topology unchanged. Its exact
    model takes each arc's end onto the arc's circle: an end off the
    circle `|p - c| = r` (in rationals) becomes the circle's rational point
    at its half-angle tangent `s`, `c + r ((1 - s^2), 2 s) / (1 + s^2)`
    where the end's `dx = p_x - c_x` is not negative, `c + r (-(1 - s^2),
    2 s) / (1 + s^2)` where it is, with `s` the binary64 rounding (once,
    to nearest) of `dy / (r + |dx|)` computed exactly (`|s| <= 1` within
    rounding); an end on its circle is kept. A joint of a line and an arc
    takes the arc's point, the line's end with it (its wall's plane passes
    through the point); a joint of two arcs of one circle (equal centres
    and radii as rationals) that circle's point, the same for both; a
    joint of two lines keeps its rounded point. The model's ends are then
    rationals of up to about 220 bits (`s` squared over the centre's and
    radius's bits) where they were binary64, and move by at most the end's
    distance from its circle plus `r 2^-52` (the half-angle's rounding):
    within the resolution, since the profile's validation holds every arc's
    ends within the resolution of its circle, so the model's vertices stay
    within the resolution of the stored ones (`check_slots`). A joint
    tangent in the body OCCT was given is tangent in no binary64 data (two
    circles' or a circle's and a line's common tangent points are
    irrational in general): after the ends are taken onto the circle, the
    line meets the circle at the joint at an angle within rounding of
    tangency and again within rounding of the joint, outside both faces,
    and the arrangement decides the joint by its exact turn as any other
    (S9c's profile test takes a smooth joint and a turning one alike). The flag that
    asks for it is the imported profile's (`Profile`'s `rounded_arcs`, set
    by `imported::prism`); a kernel profile's arc ending off its circle
    stays refused by design (S9c: the caller placed the end, it is not a
    rounding). The step's code is its own module
    (`solid/boolean/curved/snapped.rs`), called by `model::Prism::new` for
    a flagged profile; `polyhedra::substituted`'s refusal goes. (3)
    *Degenerate and refused.* S9's rules unchanged: a partner's face within
    the resolution of the imported prism's plane wall (coplanar in the
    construction: two faces within the resolution of one plane), a
    partner's face tangent to its cylinder wall, a vertex on a face;
    `OutOfDomain("an imported prism's arcs of two circles meeting at a
    joint (S9e.4b.4)")` where an end is shared by arcs of different circles
    and lies off either. S9e.4a's three `dee_turn` cases, declared
    `unsupported` as S9e.4b's, are this step's and are declared solid with
    its kernel. (4) *Its other queries* are S9e.4a's (classification by
    the construction, mass by the stored topology, rigid motion
    rebuilding the construction, whose flagged profile is taken onto its
    circles again in the moved frame: the same local data, the same
    points). (5) *Fuzzing.* The `boolean` target's `IMPORTED` stage reaches
    such prisms already (the object written, read back and imported:
    its 42 refusals); no new stage. (6) *Evidence first.* Bodies OCCT
    writes (`occt_boolean_oracle.cpp`'s `write` blocks, `MakePrism` of
    profiles whose world coordinates are the turned frames' roundings)
    under `rust/fixtures/imported/`: a stadium and a seam-split circle in
    the `TILT` frame, a rounded rectangle (four fillets) in the `R125`
    frame, a profile whose arc crosses its lines at an angle in the
    `TURN30` frame, and two arcs of different circles meeting (a lens,
    declared `unsupported`, S9e.4b.4's); cases against boxes, a crossing
    rod and a sphere, as object and as tool, both inputs imported, a chain,
    a box sharing the stadium's flat wall's plane and a box tangent to its
    arc wall (declared `degenerate`); the reference the constructions OCCT
    was given through S9e.3a's chained reference with S9e.4a's checks, each
    file's stored vertices on the construction's surfaces and its arcs'
    ends, in its cap's stored frame as written, off their circles
    (`generate_imported_arcs_boolean_fixtures.py --check`, a CI group
    `imported-arcs`, `test_imported_arcs_boolean_reference.py`); a native
    capture before `snapped.rs` exists (`compare_imported_arcs_boolean.py`
    keyed on it, the kernel's probe refusing every case); then the kernel
    and its tests (`tests/imported_arc_booleans.rs`: enclosures within
    `1e-9` of the reference, degenerate cases refused, histories over the
    stored ids, determinism, rigid motion).
  * **S9e.4b.1 evidence (2026-10-03).** Five bodies OCCT wrote
    (`boolean-imported-arcs-bodies.txt`, `MakePrism` of profiles whose
    world coordinates are their turned frames' roundings, written by
    `compare_imported_arcs_boolean.py --write-bodies` to
    `rust/fixtures/imported/`): `slot`, a stadium (lines tangent to its
    arcs) in the `TILT` frame; `halves_turn`, a circle as two arcs (a
    joint of two arcs of one circle) in the `TURN30` frame; `rounded`, a
    rectangle with four fillets in the `R125` frame; `notch`, an arc
    meeting its lines at an angle in the `TURN30` frame; `lens`, two arcs
    of different circles meeting, in the `TURN30` frame.
    `generate_imported_arcs_boolean_fixtures.py --check`: 36 cases of 12
    groups (27 solid, 6 declared `degenerate`, 3 `unsupported`; 30 of
    class `prism`, 3 `both`, 3 `chain`): the slot against a box across its
    arc end, as object and as tool, a vertical rod through it and a ball
    across its arc wall (a sphere's meeting with its cylinder); the split
    circle against a box across its seam; the rounded rectangle against a
    `TILT` slab (two solids cut); the notch against a box across its
    arc's crossing joint; the notch and the slot both imported (two solids
    cut); the slot less the rod, then with the box; declared `degenerate`
    a `TILT` box on the plane of the slot's flat wall (coplanar in the
    construction) and a `TILT` box whose wall touches the slot's arc wall
    along a generatrix; declared `unsupported` the lens against a box
    (S9e.4b.4). The reference is the constructions OCCT was given through
    S9e.3a's chained reference with S9e.4a's checks
    (`generate_imported_boolean_fixtures.evaluate_chain`), relative to the
    case's size: the two families within 4.1e-32, each solid's closed
    form 1.5e-40, the pair identities 4.4e-38 and the area identity
    1.5e-40, Monte Carlo 3.2 standard errors, solid counts by rays at two
    resolutions, every meeting's sine at least 0.30 and events at least
    3.2e-3 of their range apart outside the declared cases, quadrature
    estimates 9.1e-33. Every file read independently (`stored_records`):
    its faces' kinds the construction's, every stored vertex within
    3.9e-17 of the size on the construction's surfaces, and its arcs'
    ends, their stored vertices' local coordinates in the construction's
    frame (the exact affine map of its stored binary64 axes, of which the
    cap's stored frame the kernel reads is a rounding) rounded once, off
    their circles exactly: 4 of the slot's 8, all of the others' (4, 16,
    4 and 4).
    `test_imported_arcs_boolean_reference.py` checks the constructions'
    arcs ending on their circles exactly, the slot's volume and area
    against its closed form (the volume times the stored axes'
    determinant), the slot less the box by a quadrature of its sections,
    the off-circle test, and the case list and its protocol rows. The
    generator's check is a CI group of its own (`imported-arcs`); Python
    3.9 and 3.12 write the same files. Corrections from the evidence: the
    plan's off-circle check is made in the construction's frame (the
    cap's stored frame is the converter's normalization again of OCCT's
    17 digits, which the reference does not emulate); the lens was moved
    from the `TILT` frame, where its joints' coordinates round back onto
    both circles (S9e.4a would take it), to `TURN30`; the split circle
    likewise from `TILT` to `TURN30` (in `TILT` its ends, off its circle
    in the construction's frame, round back onto it in the cap's stored
    frame, and S9e.4a's kernel took it: found by the capture's check that
    the kernel refuses every case before this step's code); and the notch
    from `(2, 0, -1/2)` to `(9/4, 1/8, -1/2)`: its vertical edge at `x = 0`
    touched the slot's arc wall along its generatrix there, an edge
    tangent to a face that the reference's events, the same height found
    twice within 1e-25 and counted once, did not flag (a trial run of the
    cases on a draft of the step's code, not committed, refused it as two
    result vertices within the resolution). The capture
    `occt-boolean-imported-arcs-preimplementation`
    (`compare_imported_arcs_boolean.py`, keyed on
    `solid/boolean/curved/snapped.rs`; the kernel's probe `unsupported` on
    all 36, S9e.4a refusing every body's arcs): every result valid, 30
    matching (volumes within 1.3e-8, areas 1.2e-8, centres 4.2e-9 of the
    size), 6 reviewed: the ball across the slot's arc wall and the notch
    and the slot crossing, BRepGProp's default integration on faces bounded
    by approximated sections (up to 1.4e-5; adaptively within 3.2e-9,
    unchanged at 1e-12), and the declared tangent fuse, which OCCT keeps as
    two solids touching along a line where the reference's rays count one;
    5 results' counts change when unified. S9e.4b.1's kernel next.
  * **S9e.4b.1 implemented** (`solid/boolean/curved/snapped.rs`,
    `Profile`'s `rounded_arcs`, set by `imported::prism`;
    `model::Prism::new` takes a flagged profile's points through it;
    `polyhedra::substituted`'s refusal and `imported::arcs_on_circles`
    gone), as the refined decisions describe. All 36 fixtures as declared
    (27 within the kernel's enclosures, each at most `1e-9` wide; the 6
    declared degenerate refused, the coplanar box as two faces within the
    resolution of one plane and the touching box as a tangency; the lens
    `OutOfDomain` as S9e.4b.4's), every history complete over the imported
    bodies' stored ids, results deterministic and moved rigidly, both
    inputs translated keeping the reference's volumes and turned too where
    no face of one input is exactly parallel to the other's cylinder (a
    turn rounds such a pair into "a plane within rounding of a cylinder's
    direction", S9's rule), every body its construction (closed-form
    volumes, points classified, stored vertices on the boundary), the
    stored frames the reference's bit for bit
    (`tests/imported_arc_booleans.rs`, 15 s at `opt-level` 2 with debug
    assertions on a host at load 29). `compare_imported_arcs_boolean.py`
    27 matches and 9 reviewed (the native measures and the tangent fuse as
    captured; entity counts: the ball's meeting with the slot's arc wall in
    exact pieces, the split circle's two faces kept), every enclosure
    within the reference with S9e.4a's `1e-12` slack; every other
    comparison unchanged (`compare_imported_boolean.py` 54/15 with its
    `dee_turn` cases refused as declared; the 28 others of the verification
    as `HANDOFF.md`'s table, `compare_step.py` 23/6 on STEP-b's SDK), the
    tools' unit tests (297) and the ledger unchanged. Amendments and
    corrections, from the implementation: (a) S9e.4a's three `dee_turn`
    cases, declared `unsupported` because S9e.4a refused the turned profile's arcs first,
    are degenerate under S9's rules (the profile's corner at the origin on
    the box's corner: "a vertex of one input on the other's face"), which
    the kernel now reaches: declared `degenerate` (CORNER in
    `generate_imported_boolean_fixtures.py`; the reference's rows
    unchanged, `compare_imported_boolean.py` still 54 matches and 15
    reviewed); (b) a turned rigid motion of a pair with a face exactly
    parallel to the other's cylinder axis (an axis-aligned box against the
    slot in the `TILT` frame) is S9's degenerate case, not this step's, so
    the motion test turns only pairs without one. The `boolean` fuzz
    target's `IMPORTED` stage needs no new code: replaying the corpus
    (1,430 inputs) and the 24 regressions with debug assertions, no
    failure, the slowest 24.6 s on a host at load 25 to 48; of the imported
    operations 663 evaluate and 220 are refused as documented (S9's
    degeneracies of the converter's frames normalized again, one torus
    meeting's budget), against 645 and 237 with S9e.4a's 42 refusals of
    arcs off their circles, none of which remains, and no joint of two
    circles is met. A trial of the DRAW survey's 31 restore cases with such
    prisms on the Rust adapter (not the survey: no volume audit, nothing
    registered): 7 evaluate with every check (`bcut_complex/H3`, `K8`,
    `bfuse_complex/C9`, `E9`, `I6`, `N1`, `N9`), 17 are refused by S9's
    rules (12 faces within the resolution of one plane, 3 edges meeting, 2
    tangencies), 1 by a cylinders' section within the resolution of a node
    and 6 as S9e.4b.4's (two circles at a joint). DRAW survey: that of
    S9f.2b.2 and S9e.4b.1, below.
    Campaign: the boolean campaign at `788f8861`, with S9e.4b.1, S9e.4b.2,
    S9f.3a, S9f.3b and the near-parallel guards (600 s, a sampled replay,
    `SPLINE_SPHERE` and `SPLINE_CONE` off) clean, 909 runs, the slowest
    input 45 s under AddressSanitizer at load 3 to 7 (`e36969f1`, an
    existing corpus input).
  * **DRAW survey of S9f.2b.2 and S9e.4b.1 (2026-10-03,
    `UPSTREAM_TESTS.md`).** At `93e6fcd0` (`s9c2-kernel` with S9e.4a,
    S9f.2b.2 and S9e.4b.1 merged; the public dataset, 120 seconds a case,
    four at once). The 1,802 self-contained cases of the Boolean group on
    both backends: every status and every refusal's reason the last
    survey's (`b0b9adc6`), 987 evaluating and registered, 592 refused, 223
    unsupported on both, none failing or timing out, the sentinels refused
    as before, `bopfuse_simple/ZP6` and the `gdml_public` tori too. The
    1,814 cases restoring a shape for a Boolean on the Rust adapter, and
    the 171 the import reaches on native DRAW too: against S9e.4a's survey
    only the 31 prisms whose arcs round off their circles move, as S9e.4b.1's
    trial found: 7 evaluate on both backends with every check
    (`bcut_complex/H3`, `K8`, `bfuse_complex/C9`, `E9`, `I6`, `N1`, `N9`:
    prisms of the `CTO9xx`, `cts` and `pro` series with one or two
    cylindrical walls fused with or cut by cylinders, boxes and prisms of
    planes), 18 are refused by S9's rules (12 faces within the resolution of
    one plane, 3 edges meeting, 2 tangencies, `bug29807_b1` a cylinders'
    section within the resolution of a node) and 6 as S9e.4b.4's (two
    circles at a joint: `bcut_complex/E8`, `P4`, `bfuse_complex/D5`, `E1`,
    `bug4993_1`, `_2`); `bcut_complex/I6`, failing in S9e.4a's run, is
    refused as S9e.4b's (that survey's amendment). Of the 171: 23 evaluate
    on both backends, 55 are refused by S9's rules, 84 are bodies none of
    the kernel's constructions, 6 S9e.4b.4's and 3 spline bodies (S9f); 24
    are unsupported natively too (`checksection`, `bopargcheck`). No other
    status or reason moves; the same 8 time out at their first restores.
    S9f.2b.2's loops and towers reach no case: the Boolean group's spline
    solids come from `nurbsconvert`, which the adapter does not run, and the
    3 restored spline bodies (`bcut_complex/L9`, `O1`, `bfuse_complex/N7`)
    are refused as imported spline faces (S9f) before any wall meets a
    cylinder. The volume audit (`vprops` and `sprops` before each
    `checkprops`): the 7 new cases' volume, area and centre native DRAW's to
    its printed digits; the 1,003 registered cases' values the last audits'
    (`b0b9adc6`'s 987, S9e.4a's 16) bit for bit on both backends, so
    S9e.4b.1's snapped arcs and S9f.2b.2's kernel move no registered value;
    the same 35 disagreements, native off in each. The 7 are registered
    (1,096 cases; the contract holds for them on both backends within 30
    seconds, Rust 1.8 to 6.5 s on a host at load 25 to 31); the ledger does
    not change (no `checknbshapes` among them). A full contract run of the
    manifest holds on both backends with the dataset (30 seconds a case),
    the slowest Boolean case 16.2 seconds (`boptuc_simple/ZK8`, 12.8 in
    the last survey, on a host at load 7 to 31), the restore cases 0.2 to
    5.0 s. No case fails, crashes or panics; no kernel change.
  * **S9e.4b.2 refined, before its code (2026-10-04).** Why it is refused
    today: `imported::recognize` takes a body of plane faces and line
    edges only as S9e.4a's prism (a first pair of opposite planes, every
    other face along their normal), so a pyramid, a pyramid's frustum, a
    wedge with slanted faces or a result of boxes in different frames is
    `OutOfDomain("an imported solid other than a prism, a sphere, a cone
    or a torus (S9e.4b)")`: the DRAW survey's 7 cases (`bugs/modalg_1/
    buc60803`, `bug102_1`, `bug102_2`: `buc60803a` and `b`, two frustums
    of a pyramid, the second standing on the first's top, fused;
    `bugs/modalg_2/bug578_1`, `_2`: `OCC578_w1` and `w2`, frustums turned
    by about a right angle about `x` on either side of `z = 0.5`, fused
    and cut from a box; `bopfuse_complex/K5`: `pro9481b`, a frustum on
    the box `pro9481a`'s top; `bfuse_complex/D9`: `CTO900_pro12559a`, a
    prism along `x` with some walls drafted by 2 degrees, fused with the
    prism `CTO900_pro12559b` on its cap). What the files store: every face
    a plane, every edge a line between two vertices (15 or 17 digits for
    surfaces, 15 for vertices); every vertex has three faces, and lies off
    their planes' common point by up to 3.0e-13 (`CTO900_pro12559a`);
    where two bodies were built to share a face (`buc60803a`'s top is
    `b`'s base, `pro9481a`'s top `pro9481b`'s base, `CTO900_pro12559a`'s
    cap holds `b`'s), the shared vertices are the same binary64 points in
    both files, but each file's planes meet near them at points 5.3e-15 to
    7.2e-15 apart (`buc60803`; a scratch measurement, the dataset not
    committed). Decisions. (1) *The exact polyhedron is its stored
    vertices*, not its planes' common points (S9e.4b refined's plan,
    amended): each vertex its stored binary64 point as a rational, each
    edge the segment between its two stored vertices, each face the
    polygon of its loops' stored vertices, exactly planar where they are
    coplanar (an axis-aligned face, a triangle, a face OCCT built from
    exact points), else the triangles of that polygon in its projection on
    the stored normal's largest coordinate plane (ears clipped, holes
    bridged; S9b.2's stored model of a Boolean's result, which decides
    such faces already: each triangle exactly planar, within the
    resolution of the stored plane, the face's fragments joined back by
    the face), and the solid's membership by an exact ray's parity over
    those triangles. Why not the planes: a vertex of four or more faces (a
    pyramid's apex) has no common point of its rounded planes in general,
    and where two bodies share a face their planes' points differ though
    their stored points agree, so a Boolean of the two would leave a
    sliver along every shared edge (refused by S9b's rule as a face
    thinner than the resolution) where the stored data meet exactly; the
    stored points keep every coincidence the files store, and lie within
    the resolution of their faces' stored planes (the converter's
    validation), so the model is the body within its resolution. Numbers:
    rationals of binary64 (S9b.2's), no new field or degree. (2)
    *Recognition.* A body whose faces are all planes and edges all lines,
    of one solid region and one shell, is S9e.4a's prism where its
    recognition, construction and match succeed (a box or a prism of
    lines, unchanged), else a polyhedron (this step), decided on its
    stored topology: no construction and no match, its entities its stored
    ones, so the Boolean's history is over the stored ids directly. (3)
    *Degenerate.* A face whose stored vertices fold it (a triangle of its
    clipping facing against the stored plane's outward normal, or none
    found of the polygon's area: a sliver face rounded over)
    `Degenerate("an imported face folded by its stored vertices")`; S9's
    rules unchanged in the Boolean (S9b's exact arrangement: a partner's
    face within the resolution of an imported face and not on it exactly,
    as a box in a turned frame flush with a face whose stored corners are
    roundings, leaves a face thinner than the resolution; solids touching
    at a point or along an edge; a face meeting itself). An ill-conditioned
    vertex (three planes nearly dependent) needs no rule: its stored point
    is taken, not computed. (4) *Its other queries.* Classification
    `Boundary` within the resolution of a face's triangles (in binary64),
    else by the exact parity; mass the stored topology's certified
    enclosure (S9e.4a's); bounds the stored edges'; a rigid motion moves
    the stored topology (as S9b's results), the model its moved vertices.
    A result with an imported polyhedron given to another polyhedral
    Boolean is S9b.2's stored model of that result. (5) *Refused.* An
    imported polyhedron against a solid with curved faces or edges (a
    prism with arcs or circles, a sphere, a cone, a torus, a spline prism
    or a result of them) and a result of one given to a Boolean of curved
    faces: `OutOfDomain("an imported polyhedron against curved faces
    (S9e.4b.4)")` (the curved engine decides on constructions; a general
    body's faces in its arrangement are the S9e text's plan, S9e.4b.4); a
    polyhedron with a cavity or several shells (none among the surveyed
    cases) `OutOfDomain("an imported polyhedron with a cavity or several
    shells (S9e.4b.4)")`. (6) *Fuzzing.* The `IMPORTED` stage imports the
    object, always a construction; the chained stage's first result (at
    most 12 plane faces, a polyhedron in any position) is written, read
    back, imported and cut by the chained box too, its volume the chained
    cut's. (7) *Evidence first.* `occt_boolean_oracle.cpp` takes rows
    `wedge` (`BRepPrimAPI_MakeWedge(gp_Ax2, DX, DY, DZ, XMIN, ZMIN, XMAX,
    ZMAX)`: a pyramid where its top is a point, a frustum) and
    `polyhedron` (faces of given points sewn into a solid), and a `write`
    block taking a Boolean of two solids (its result's one solid written).
    Bodies OCCT writes under `rust/fixtures/imported/`: a tetrahedron and
    an octahedron (points, rounded in a turned frame or exact), a pyramid
    (its apex of four faces), a frustum and a wedge with slanted faces in
    turned frames, results of boxes (a box less a tilted box, a box fused
    with a turned box over other heights), and the survey's shapes from
    exact points (two frustums sharing a face, a frustum on a box's top, a
    partly drafted prism with a prism on its cap, two frustums at right
    angles fused and cut from a box); cases against boxes and slabs,
    another imported body, as object and tool, and chains; declared
    `degenerate` a slab in the pyramid's frame on its base's plane;
    declared `unsupported` a box with a cavity. The reference
    (`imported_polyhedra_boolean_reference.py`, independent of the
    kernel): each body the construction OCCT was given as convex cells
    (hulls of exact points: a wedge's local corners on the frame's stored
    axes, a polyhedron's binary64 points, a Boolean's cells by S9b's
    polyhedral reference), its Booleans' volumes three ways (inclusion
    and exclusion, the result's cells, the divergence over its boundary),
    areas two ways (the inputs' faces classified, the result's cells'
    boundary), prismatoid closed forms, Monte Carlo against the inputs'
    half-spaces, solid counts by cells meeting in area, clearances
    outside the declared cases; `generate_imported_polyhedra_boolean_
    fixtures.py --check` (a CI group `imported-polyhedra`, each file's
    stored vertices on the construction's boundary and the body no prism),
    `test_imported_polyhedra_boolean_reference.py`; a native capture before
    `solid/boolean/polyhedra/imported.rs` exists
    (`compare_imported_polyhedra_boolean.py` keyed on it, the kernel's
    probe refusing every case); then the kernel and its tests
    (`tests/imported_polyhedra_booleans.rs`).
  * **S9e.4b.2 evidence (2026-10-04).** `occt_boolean_oracle.cpp` takes the
    rows `wedge` and `polyhedron` and a `write` block of two solids about a
    `boolean` row (its result's one solid written). Fifteen bodies OCCT
    wrote (`boolean-imported-polyhedra-bodies.txt`, written by
    `compare_imported_polyhedra_boolean.py --write-bodies` to
    `rust/fixtures/imported/`): `tetra`, a tetrahedron of points rounded in
    the `TURN30` frame; `octa`, an octahedron of exact points (every vertex
    of four faces); `pyramid` (`MakeWedge` in the `TILT` frame, its apex of
    four faces), `truncated` (a pyramid's frustum in the `TURN30` frame)
    and `wedge` (slanted on three sides, in the `R125` frame); `notched`, a
    box less a box in a skew frame, and `ell`, a box fused with a turned box
    over other heights (results OCCT computed); `hollow`, a box with a
    cavity; and the survey's shapes from exact points: `steps_low` and
    `steps_high` (`buc60803a`, `b`: the second's base the first's top),
    `pedestal` (`pro9481b`: a frustum on a kernel box's top, its corners
    decimals of 15 digits), `draft` (`CTO900_pro12559a`: a prism along `x`
    of a profile with a reflex corner, two walls drafted) with `ridge`
    (`CTO900_pro12559b`, `MakePrism` on part of its far cap, S9e.4a's
    prism), `vane_up` and `vane_down` (`OCC578_w1`, `w2`: frustums in frames
    turned by right angles about `x`, their bases one square). The
    reference (`imported_polyhedra_boolean_reference.py`): every solid as
    convex cells in exact Fractions (a hull of each construction's points:
    a wedge's corners on the frame's stored axes, a polyhedron's binary64
    points; a Boolean's cells by S9b's polyhedral reference), a Boolean's
    result as S9b's cells. `generate_imported_polyhedra_boolean_fixtures.py
    --check`: 48 cases of 16 groups (40 solid, 4 empty, 1 `degenerate`, 3
    `unsupported`; 33 of class `polyhedron`, 9 `both`, 6 `chain`), each body
    against boxes, slabs and rods in the `XY`, `TILT` and `TILTX` frames or
    another imported body, as object and tool, and two chains (the vanes
    fused, then a box with them, swapped; the notched box less a rod, then a
    `TILT` slab); declared `degenerate` the fuse of a `TILTX` slab on the
    pyramid's base's plane (its cut and common the construction's: two
    stored base corners on the slab's plane, two off it on the pyramid's
    side by 7e-16), `unsupported` the box with a cavity against a slab.
    Checks, relative to the case's size: volumes three ways (inclusion and
    exclusion, the result's cells, the divergence over its boundary) and
    the pair identities on the cells' volumes exactly equal, areas two ways
    within 1.8e-40 and the area identity 6.0e-40, every imported body's
    closed form (prismatoid, a tetrahedron's determinant, an octahedron's
    diagonals) exactly, Monte Carlo 3.3 standard errors, solid counts by
    cells meeting in area; outside the declared group every vertex and edge
    of one input at least 4.0e-3 of the size from the other's faces and
    edges where not on them (196 exact contacts, all between bodies whose
    stored points are their construction's exactly: the shared faces), and
    crossing faces at sines of at least 0.36. Every file read
    independently (`stored_records`): its faces planes, every stored vertex
    within 6.6e-16 of the size of the construction's boundary (the exact
    bodies' stored points their construction's bit for bit), and no
    construction a prism (but `ridge`'s). `test_imported_polyhedra_
    boolean_reference.py` checks hulls, a wedge's corners and volumes, cubes'
    Booleans and a union's boundary in closed form, the margins on contacts
    made on purpose, the prism test and the case list and its protocol
    rows; the generator's check is a CI group of its own
    (`imported-polyhedra`); Python 3.9 and 3.12 write the same files.
    Corrections from the evidence: the pyramid's frustum named `truncated`
    (S9e.4a's cone body is `frustum.brep`, found by the files' check); the
    skew box moved off `(8, 8, 2.5)`, where its edge passed exactly through
    the rods' edges (found once the margins measured edges too: a trial of
    the cases on a draft of the step's code, not committed, refused it as a
    face thinner than the resolution and a non-manifold vertex), the
    octahedron's slab and the tetrahedron moved off vertices within
    rounding of the other's faces; the flush slab's cut and common, which
    the same trial evaluated within the reference, declared as such. The
    capture `occt-boolean-imported-polyhedra-preimplementation`
    (`compare_imported_polyhedra_boolean.py`, keyed on
    `solid/boolean/polyhedra/imported.rs`; the kernel's probe `unsupported`
    on all 48, S9e.4a refusing every body but `ridge` and the pairs with
    it): every result valid, 47 matching (volumes within 8.9e-16, areas
    5.5e-16, centres 2.5e-16 of the size: plane faces, BRepGProp exact), 1
    reviewed: the declared flush fuse, which OCCT makes one solid sharing
    the pyramid's base within its tolerance where the reference's rounded
    frames keep two apart (its volume the reference's, its area less twice
    the base square); 2 results' counts change when unified (the drafted
    prism and the ridge). S9e.4b.2's kernel next.
  * **S9e.4b.2 implemented** (`solid/boolean/polyhedra/imported.rs`,
    `imported::Recognized`, `polyhedra.rs`'s stored model for an imported
    polyhedron, `given.rs`'s refusal), as the refined decisions describe:
    an imported body of planes and lines is S9e.4a's prism where its
    recognition, construction and match succeed, else a polyhedron decided
    on its stored vertices (S9b.2's stored model: each face's outward
    normal from the region behind it, its polygon of stored vertices cut
    into exactly planar triangles, each edge between its stored vertices,
    membership by parity), its history over its stored ids directly; it
    classifies `Boundary` within the resolution of a face's triangle and
    by parity off them, and moves with its stored topology. All 48
    fixtures as declared (40 within the kernel's enclosures, each at most
    `1e-9` wide, and 4 empty; the flush fuse refused, `Degenerate("a
    face using an edge both ways")`; the cavity `OutOfDomain` as
    S9e.4b.4's), every history complete over the imported bodies' stored
    ids, results deterministic and moved rigidly, both inputs translated
    and turned keeping the reference's volumes (the cases without exact
    contacts, which a turn's rounding breaks), every body its
    construction (its exact volume within `1e-9`, points classified, its
    stored vertices on its boundary), a pyramid against a cylinder and a
    result of it given with a cylinder `OutOfDomain` as S9e.4b.4's, the
    stored frames the reference's bit for bit
    (`tests/imported_polyhedra_booleans.rs`, 10 s at `opt-level` 2 with debug
    assertions on a host at load 12 to 21).
    `compare_imported_polyhedra_boolean.py` 47 matches and 1 reviewed (the
    flush fuse, as captured), every enclosure within the reference with
    S9e.4a's `1e-12` slack, every entity count OCCT's after unifying;
    `compare_imported_boolean.py` 54/15 and
    `compare_imported_arcs_boolean.py` 27/9 unchanged (and `compare_step.py`
    23/6 on STEP-b's SDK, polyhedral 43/2, given 36/0, chained 24/6,
    boolean 45/0), the tools' unit tests (308) and the ledger unchanged.
    Amendments, from
    the implementation: (a) an imported face whose stored vertices are
    coplanar exactly is decided on that plane, its fragments joined with
    any face's on it (a construction's or another import's), as a
    construction's are: S9b.2 joins a stored face's fragments by the face
    alone, which left the drafted prism's walls apart from the ridge's on
    their shared planes in the fuse (17 faces where OCCT's unified result
    has 12; a kernel result given again keeps S9b.2's rule); (b) an
    imported prism of lines against an imported polyhedron is decided on
    its stored vertices too, not S9e.4a's construction, whose corners,
    re-derived from its rounded local coordinates, can miss the vertices
    the two files share by an ulp (`bopfuse_complex/K5`, a frustum on an
    imported box's top, refused as a face thinner than the resolution
    before; found by the trial below; the fixtures' `draft_ridge`
    exercises the path, every result unchanged). The `boolean` fuzz
    target's `IMPORTED` stage also writes the chained stage's first result
    of plane faces, reads it back, imports it and cuts it by the turned box
    again, its volume the chained cut's within `1e-9`: replaying the corpus
    (1,430 inputs) and the 25 regressions with debug assertions, no
    failure, the slowest 23.5 s on a host at load 28 to 37; 189 such first
    results reach it, 188 cut within the chained cut's volume and 1
    refused (a cylinder partner, S9e.4b.4's). A trial of the DRAW survey's
    7 restore cases with such polyhedra (not the survey: no volume audit,
    nothing registered), on the Rust adapter and native DRAW: 4 evaluate on
    both with every check (`bugs/modalg_1/buc60803`, `bug102_1`,
    `bug102_2`, two frustums sharing a face, fused; `bopfuse_complex/K5`
    since amendment (b)); 3 are refused: `bugs/modalg_2/bug578_1` and `_2`,
    whose frustums `OCC578_w1` and `w2` are turned about 1.3e-6 and 4.0e-6
    from right angles, so their bases lie 6.6e-7 to 2.0e-6 apart, above the
    resolution: their fuse is two solids, which the adapter refuses as the
    next Boolean's argument (one solid); `bfuse_complex/D9`, a
    corner the two files share stored 1e-13 apart (`z = 56.5616376719611`
    and `56.561637671961`), a face thinner than the resolution (S9's rule).
    DRAW survey: that of S9e.4b.2, S9f.3a and S9f.3b, below (the trial's 4
    cases registered, the other 3 refused as it found). Campaign: the
    boolean campaign at `788f8861`,
    with S9e.4b.1, S9e.4b.2, S9f.3a, S9f.3b and the near-parallel guards
    (600 s, a sampled replay, `SPLINE_SPHERE` and `SPLINE_CONE` off) clean,
    909 runs, the slowest input 45 s under AddressSanitizer at load 3 to 7
    (`e36969f1`, an existing corpus input).
  * **S9e.4b.3 refined, before its code (2026-10-04).** Why each is refused
    today, from the files the DRAW cases restore (read by the converter;
    the dataset is not committed). (a) *A plane's pieces of a sphere*
    (`so1` to `so7`, 38 cases): `imported::recognize` takes a body of one
    sphere face and planes only as S9e.4a's sphere, cap or zone, every
    plane face a disc normal to the stored axis bounded by one ring. `so1`
    and `so4` (the hemisphere of radius 10 about the origin above `z = 0`,
    the cap above `z = 5`) have such a disc, but OCCT split its rim into two
    half circles at the stored sphere's seam and opposite it (two arcs and
    two vertices after the converter), the stored pole a vertex loop of the
    sphere face. `so2`, `so3` and `so5` are spherical wedges: the ball
    above `z = 0` (`z = 5`) between two planes through its axis at 45 and
    135 degrees (20 and 110), their meetings with the sphere arcs between
    three-face vertices, the stored pole the vertex where the two planes'
    line meets the sphere; `so6` and `so7` are `so2` and `so3` turned (about
    `x` by 20 degrees, about a skew axis), the stored sphere's frame turned
    with them, so its pole lies within 1e-15 of the planes' line's point on
    the sphere, `so7`'s vertices OCCT's approximations off their surfaces by
    up to 1.3e-9. Every DRAW case is a Boolean of two of these, all on the
    one sphere of radius 10 about the origin: both inputs' sphere faces on
    one surface, which S9d.2 refuses (`spheres::sphere_sphere`,
    `Degenerate("two spheres about one centre")`). (b) *A cone's*
    (`shading_132`, `bugs/modalg_1/buc60926`): the frustum of radii 1 and 2
    over height 2 less the quarter between two planes through its axis, so
    not convex, its planes through the frustum's virtual apex (S9d.3a's
    `Degenerate("a plane through a cone's apex")`); its case fuses it with
    a unit sphere through two of its rim's vertices (a vertex of one input
    on the other's face, S9's rule). (c) *S9e.2's `Clipped` and `Half`*
    against curved faces: the refusals naming S9e.4 in `model::Prism::new`,
    `given.rs` and `curved::build`. Decisions. (1) *Sub-steps.*
    **S9e.4b.3a** (this step): an imported body of one sphere, cylinder or
    cone face and plane faces that is its primitive common the half-spaces
    of its planes, against prisms, spheres, cones, tori and other such
    bodies whose curved faces do not share its surface, on the model of
    (2); **S9e.4b.3b**: S9e.2's `Clipped` and `Half` on the same model (a
    split's primitive and its plane in the solid's frame; a `split` row in
    the case protocol and the oracle's common with a half-space);
    **S9e.4b.3c**: what the DRAW cases need beyond it: two pieces of one
    sphere (both inputs' faces on one sphere, as S9c's coincident cylinder
    walls), a ring split by stored vertices (`so1`, `so4`: the match taking
    a construction's ring as the stored arcs), pieces not convex in their
    planes (`shading_132`: the primitive less a hull) and two plane faces on
    one plane. Without S9e.4b.3c no DRAW case evaluates: this step's trial
    expects all 39 refused, the so-pairs as two inputs on one sphere. (2)
    *The model: a first Boolean.* The body is the given model (S9e.1,
    S9e.3a) of `primitive common hull`. The primitive: a whole sphere on
    the stored frame and radius exactly (its axis free, S9d.1's whole
    sphere); a cylinder, the circle of the stored radius about the stored
    frame's origin extruded along its normal over the body's span past both
    ends by a quarter of it (the span the heights of its stored edges'
    bounds' corners, rounded once); or a cone, the stored frame moved along
    its axis to the lower height and its radii `radius + w tan a` at both
    heights rounded once, the apex where the range reaches it. The hull:
    the convex body of each plane face's stored plane `o + u x + v y` (the
    stored axes as rationals, its normal `x * y` exactly, as S9e.4a reads a
    cap), on the side of the face's region, bounded by a cube about the
    primitive (its faces never meeting it). It is a new leaf model
    (`curved/pieces.rs`'s `Hull`): its faces the planes' convex polygons of
    the three planes' points lying in every half-space (rationals), its
    edges the lines between them, its membership every half-space (pushed
    in turn as the other models'), a face's region the other half-spaces.
    The first arrangement (S9c.1's engine on the two models, the
    primitive's seams tried in turn) and its assembly give one solid,
    matched to the stored topology by S9e.2's geometric match; the given
    model (S9e.1's: the primitive's and the hull's faces holding kept
    pieces, the arrangement's edges on their boundaries, exact vertices,
    the first Boolean's set function as membership) carries the stored
    ids, so a Boolean's history is over them directly (no construction ids,
    no names). Numbers: the engine's own (three planes' points rationals;
    two planes' line meeting a sphere, cylinder or cone at a quadratic
    surd; a plane's section a sphere's `Circ`, a cylinder's conic over its
    angle, a cone's `ConeSec`): no new field or degree. (3)
    *Recognition.* After S9e.4a's construction fails, a body not planar
    (S9e.4b.2's) of one solid region and one shell, one sphere, cylinder or
    cone face and plane faces, no spline edge, is a plane piece; its model
    is built and matched on import (`pieces::check`). Unmatched (not convex
    in its planes, a stored vertex splitting a ring)
    `OutOfDomain("an imported plane piece other than its primitive common
    its planes' half-spaces (S9e.4b.3c)")`; two plane faces on one plane
    `OutOfDomain("an imported plane piece with two faces on one plane
    (S9e.4b.3c)")`. (4) *Degenerate and refused.* S9's rules unchanged in
    both arrangements (a plane through a cone's apex, real or virtual; a
    section through the stored sphere's pole off its meridians; faces within
    the resolution of one plane; tangencies; a vertex on a face). One
    amendment of S9d.1's pole rule: a section within the resolution of the
    stored sphere's pole takes as its pole vertex a vertex of both faces
    already within the resolution of the pole, where there is one (two
    planes' line through a turned frame's pole, `so6` and `so7`: otherwise
    a second vertex within rounding of the first, refused as a result
    thinner than the resolution). Both inputs' faces on one sphere where
    either is a piece: `OutOfDomain("faces of both inputs on one sphere
    (S9e.4b.3c)")` (S9d.2's `Degenerate` otherwise, unchanged). (5) *Its
    other queries.* Classification: the primitive's, then each plane's side,
    within the resolution; mass: the stored topology's certified enclosure
    (S9e.4a's); bounds: the primitive's between the body's axial ends (a
    whole sphere's box), widened by the stored edges'; a rigid motion moves
    the stored topology and reads the piece off it again. A result of a
    piece given to another Boolean is S9e.3a's (its first arrangement's
    leaf the piece's given model). (6) *Fuzzing.* The `boolean` target's
    chained stage's first result of one sphere, cylinder or cone face and
    planes is written, read back, imported and cut by the turned box again,
    its volume the chained cut's (S9e.4b.2's stage widened). (7) *Evidence
    first.* Bodies OCCT writes (`occt_boolean_oracle.cpp`'s `write` blocks:
    a primitive's rows, `boolean common` or `cut`, a box's rows) under
    `rust/fixtures/imported/`: balls' wedges between three planes through
    their centre (in the world's frame and a tilted one, the stored pole on
    the planes' line, the tilted one's within rounding of it), a wedge
    above a plane off the centre (`so5`'s), a ball cut by one tilted plane
    and a band between two; a cylinder cut by two oblique planes and a
    cylinder's wedge between two oblique planes meeting across it; a
    frustum cut obliquely (an ellipse) and a cone with its apex cut by a
    tilted plane; cases against boxes, slabs, rods, spheres, another piece
    and an S9e.4a body, as object and as tool, and a chain; declared
    `degenerate` a box flush with a piece's plane face and a box touching a
    piece's sphere at a point; declared `unsupported` two pieces of one
    sphere and a ball less a box's corner (S9e.4b.3c). The reference: the
    constructions OCCT was given through S9e.3a's chained reference (a
    piece its first Boolean: `(P common B) op C`, swapped, two pieces
    `(P common B) op (Q common D)`, chained), with S9e.4a's checks and each
    piece's closed form where it has one, each file's stored vertices on
    the construction's surfaces, its faces one curved and planes, and the
    body no S9e.4a construction (`generate_imported_pieces_boolean_
    fixtures.py --check`, a CI group `imported-pieces`,
    `test_imported_pieces_boolean_reference.py`); a native capture before
    `curved/pieces.rs` exists (`compare_imported_pieces_boolean.py` keyed on
    it, the kernel's probe refusing every case); then the kernel and its
    tests (`tests/imported_piece_booleans.rs`: enclosures within `1e-9` of
    the reference, degenerate cases refused, histories over the stored ids,
    determinism, rigid motion).
  * **S9e.4b.3a evidence (2026-10-04).** Eight bodies OCCT wrote
    (`boolean-imported-pieces-bodies.txt`, `write` blocks of a sphere's row,
    `boolean common` (`bitten`: `cut`) and a box's or prism's rows, written
    by `compare_imported_pieces_boolean.py --write-bodies` to
    `rust/fixtures/imported/`), every one a sphere's piece: `octant`, the
    ball of radius 5 about `(5, 5, 4)` common the corner of its frame's
    axes (three planes through the centre); `octant_tilt`, the same in
    another frame; `upper`, the wedge above a parallel's plane between two
    meridian planes (`so5`'s); `lune`, the ball between two meridian
    planes; `half`, the ball on one side of a meridian plane; `zone_wedge`,
    between two meridian planes and two parallels' planes; declared
    S9e.4b.3c's: `octant_low`, the ball of `octant` below a parallel's plane
    (S9e.4a's cap) and `bitten`, a ball less a box's corner. The frames are
    rational rotations (`SKEW` to `SKEW4`: normals `(8, 4, 1)`, `(4, 4, 7)`,
    `(10, 11, 2)`, `(4, 1, 8)` over 9 or 15), each plane a meridian plane or
    a parallel's plane of its sphere's stored frame.
    `generate_imported_pieces_boolean_fixtures.py --check`: 45 cases of 15
    groups (33 solid, 6 declared `degenerate`, 6 `unsupported`; 33 of class
    `sphere`, 9 `both`, 3 `chain`): each piece against boxes, a `TILT` slab,
    rods upright and along `y`, a ball and an upright cone across its sphere
    face, as object and as tool, two pieces of two spheres, a piece and
    S9e.4a's imported box, and a chain (the octant less a rod, then with a
    `TILT` slab); declared `degenerate` a box in the turned octant's frame on
    its base plane (within the resolution of one plane in the file) and a
    box touching the half ball's sphere inside its face; declared
    `unsupported` (S9e.4b.3c) the octant and `octant_low` (one sphere) and
    the bitten ball. The reference is the constructions OCCT was given
    through S9e.3a's chained reference, each piece its first Boolean (`(P
    common B) op C`, swapped, `(P common B) op (Q common D)`, chained), with
    S9e.4a's checks relative to the case's size: the two families within
    7.5e-37, each solid's closed form 2.0e-40, the pair identities on the
    last Boolean's arguments 4.4e-38 and the area identity 9.9e-41, each
    piece's closed form 1.3e-41 (the corners' solid angles by Van Oosterom
    and Strackee, the lune's dihedral angle, the half), Monte Carlo 3.5
    standard errors, quadrature estimates 1.7e-32, solid counts by rays at
    two resolutions and each piece one solid, every meeting's sine at least
    0.31 and events at least 9.2e-6 of their range apart outside the
    declared groups. Every file read independently (`stored_records`): its
    faces one sphere and planes, every stored vertex within 7.6e-16 of the
    size on the construction's surfaces, and every construction no S9e.4a
    construction (a plane not normal to its sphere's axis) but `octant_low`.
    `test_imported_pieces_boolean_reference.py` checks the closed forms
    (exact axes against an eighth, a quarter and a half of the ball; turned
    axes' solid angles against the faces' spherical excess), the octant on
    the chained reference against its closed form, the S9e.4a test and the
    case list and its protocol rows. The generator's check is a CI group of
    its own (`imported-pieces`); Python 3.9 and 3.12 write the same files.
    Corrections from the evidence, amending the refined decisions' plan (7):
    (a) only spheres' pieces can be written and read: OCCT writes a
    sphere's section by a plane other than a meridian plane or a parallel's
    plane of its stored frame with a B-spline pcurve the converter does not
    certify on the sphere (`UncertifiedPcurveOffEdge`: a ball cut by one
    tilted plane, a band between two), and a cylinder's or a cone's oblique
    section as an ellipse record the `.brep` reader does not read
    (`Ellipse`, left to the import track: four such bodies written and
    refused before any check), so the planned cylinder's and cone's pieces
    are the kernel's own split pieces' topologies imported in the kernel's
    tests instead (no reference); (b) the reference sweeps a sphere by
    meridians about the world's `z` and its parallels, and a plane holding
    the world's `z` through the centre lies along a whole meridian there:
    the first bodies, on the world's axes, kept their groups running for
    over an hour each, so every body is on a turned rational frame (none of
    whose axes is a world axis or normal to one, nor one of the kernel's
    whole sphere's own axes, whose split would run through the planes'
    line); (c) a frame normalized differently by macOS's `hypot` than by a
    correctly rounded one (a first `SKEW3` of normal `(2, 10, 11)`) cannot
    carry a kernel-built partner (`frames.tsv`'s bits), so the frames are
    chosen among those both round alike; (d) OCCT splits a sphere face its
    seam crosses into two faces (two faces on one sphere, S9e.4b.3c's), so
    the half ball's sphere has its frame's `x` reversed; (e) the zone's
    wedge was symmetric about its centre, its walls' generatrices leaving
    the sphere at both ends at one parameter (two events 1.9e-18 apart),
    now between `-3/2` and `5/2`; (f) partners moved where a box's face was
    tangent to the zone's sphere, two balls were internally tangent (the
    pieces' pair) and a rod's or a pair's result was a sliver whose rays'
    counts disagreed between resolutions (the turned octant moved and its
    frame changed; the pair the octant and the turned octant). The capture
    `occt-boolean-imported-pieces-preimplementation`
    (`compare_imported_pieces_boolean.py`, keyed on
    `solid/boolean/curved/pieces.rs`; the kernel's probe `unsupported` on
    all 45, S9e.4a refusing every body): every result valid, 35 matching
    (within 1.6e-8 of the reference, the worst the chain's cut, its faces
    bounded by OCCT's approximated sections), 10 reviewed: the rod's and the
    cone's meetings with a piece's sphere (`tilt_rod`, `lune_cone`,
    `half_rod`, 9 cases, volumes up to 5.4e-7 relative off by BRepGProp's
    default integration, within 3.2e-9 measured adaptively at 1e-10 and
    1e-12 by a diagnostic build) and the declared touching fuse, which OCCT
    keeps as 2 solids sharing the point where the reference's rays count 1
    (its volume and area the reference's as printed); 6 results' counts
    change when unified. S9e.4b.3a's kernel next.
  * **S9e.4b.3a implemented** (`solid/boolean/curved/pieces.rs`,
    `imported::Recognized::Piece`, `model.rs`'s hull faces and edges,
    `given::built`), as the refined decisions describe: an imported body
    of one sphere, cylinder or cone face and plane faces that is none of
    S9e.4a's constructions is a plane piece, its primitive (the whole
    sphere on its stored frame and radius, or a cylinder or a cone on the
    curved face's stored frame reaching a quarter of the body's span past
    its ends, a cone clamped at its apex) common the half-spaces of its
    plane faces' stored planes; the model is S9e.1's given model of that
    Boolean, its second input a hull leaf model (the planes' convex
    polygons bounded by a cube about the primitive, their corners three
    planes' rational points, membership every half-space), the arrangement
    tried at the seams in turn as a whole sphere's, matched to the stored
    topology on import (S9e.2's match), its history over the stored ids
    directly; it classifies by the model and moves with its stored
    topology, read off it again. All 45 fixtures as declared (33 within the
    kernel's enclosures, each at most `1e-9` wide; `tilt_flush` refused as
    `Degenerate("two faces within the resolution of one plane")`,
    `half_touch` as a tangency; the octant against `octant_low`
    `OutOfDomain("faces of both inputs on one sphere (S9e.4b.3c)")`, the
    bitten ball refused on import as S9e.4b.3c's), every history complete
    over the imported bodies' stored ids and none from a primitive or a
    hull, results deterministic and moved rigidly, both inputs translated
    and turned keeping the reference's volumes, every body a piece (its
    volume its closed form, a point in it inside and its mirror in the
    centre outside, its stored vertices on its boundary), the kernel's own
    split pieces of a cylinder and a frustum by oblique planes and a zone's
    halves imported as pieces, their Booleans with a box and a ball obeying
    the pair identities, the octant against a whole torus, the stored
    frames the reference's bit for bit (`tests/imported_piece_booleans.rs`,
    40 s at `opt-level` 2 with debug
    assertions on a host at load 27 to 41, 19 s in release). `compare_imported_pieces_boolean.py` 32 matches and 13
    reviewed (the 10 captured; with the kernel, entity counts in 6:
    `lune_ball_fuse`, the cone's three and the chain's fuse and cut, the
    same faces as OCCT's unified result, each splitting its sections at its
    own points and seams), every enclosure within the reference with
    S9e.4a's `1e-12` slack; every other comparison unchanged
    (`compare_imported_boolean.py` 54/15, `compare_imported_arcs_boolean.py`
    27/9, `compare_imported_polyhedra_boolean.py` 47/1, given 36/0, given
    curved 25/23, given met 8/42, chained 24/6, and the rest of HANDOFF's
    table and `compare_step.py` 23/6 on STEP-b's SDK), the suite and the
    tools' unit tests (318) passing and the ledger unchanged. Amendments, from the
    implementation: (a) a section within the resolution of a stored
    sphere's pole takes a vertex of both faces already within the
    resolution of the pole as its pole vertex (`graph.rs`): two planes'
    line through a turned frame's pole meets it at a rounded point, and the
    arrangement made a second vertex there (`so6`'s and `so7`'s shape
    refused as a result thinner than the resolution); (b) a
    sphere face's loop through a pole whose pcurves turn half a turn there
    (a meridian circle through both poles, the half ball's plane) winds
    none, the face closing on it without a pole vertex loop beside it
    (`assemble.rs`); (c) a given model's circle matched to a stored circle
    or arc whose frame turns against it is read the other way
    (`Given::flip`, S9e.2's rule for lines, now for circles too: the
    pieces' rims); (d) faces of both inputs on one sphere where either is a
    piece are refused before the arrangement, `OutOfDomain` as
    S9e.4b.3c's (two pieces of one sphere meet along their whole common
    sphere, which the arrangement does not decide); (e) a piece is never
    S9b.2's substitution of a stored polyhedron (`polyhedra.rs`), its model
    the curved engine's; (f) evidence correction (a): only spheres' pieces
    come from `.brep` files, the cylinder's and the cone's tested on the
    kernel's own split pieces. The `boolean` fuzz target's `IMPORTED`
    stage also imports the chained stage's first result of one sphere,
    cylinder or cone face and plane faces (a piece) and cuts it by the
    turned box again, its volume the chained cut's within `1e-9`: replaying
    the corpus (1,432 inputs) and the 26 regressions with debug
    assertions, no failure, the slowest 13.5 s on a host at load 24 to 41;
    90 such first results reach it, 17 imported and cut within the chained
    cut's volume, 63 not written by the kernel's writer (a projection or
    sinusoid pcurve of a section by the turned box), 6 read back with
    pcurves off their edges and 4 refused on import as S9e.4b.3c's (2 not
    their primitive common their planes, 2 with two faces on one plane).
    A trial of the DRAW survey's 39 S9e.4b.3 restore cases (not the
    survey: nothing registered), on the Rust adapter: none evaluates, each
    refused as declared, S9e.4b.3c's or S9's: 18 (every case restoring
    `so1` or `so4`, whose rims OCCT split at the seam) as a piece other
    than its primitive common its planes' half-spaces, 20
    (pairs of `so2`, `so3`, `so5`, `so6` and `so7`) as faces of both
    inputs on one sphere, and `buc60926` `Degenerate("a
    plane through a cone's apex")` (`shading_132`, a three-quarter frustum
    whose planes pass through its virtual apex); native DRAW passes all 39
    (`viewer_skipped`), and before S9e.4b.3a all 39 were refused as an
    imported solid other than a prism, a sphere, a cone or a torus.
    Pending: S9e.4b.3b (S9e.2's `Clipped` and `Half` on this model) and
    S9e.4b.3c (pieces of one sphere, split rims, pieces not convex in their
    planes) and the Linux record of the capture. DRAW survey: that of
    S9e.4b.3a, S9e.4b.3b and the fuzz fixes, below (the 39 refused as the
    trial found, none evaluating). Campaign:
    the boolean campaign at `dd51af05` (600 s, a sampled replay, `IMPORTED`,
    `SPLINE_SPHERE` and `SPLINE_CONE` on) clean, 838 runs, the slowest input
    60 s under AddressSanitizer at load 6 to 8 (`e36969f1`, an existing
    corpus input reaching none of S9e.4's or S9f's code, 91 G instructions
    with debug assertions; 45 and 52 s in the last two campaigns at lower
    load), at the target's limit.
  * **S9e.4b.3b refined, before its code (2026-10-04).** Why it is refused
    today: a split piece (S8a.2's `Clipped`, a prism of lines, arcs and
    circles split by a plane oblique to its axis; S8c.2's and S8d.2's
    `Half`, a cone, frustum, cap or zone split by a plane through its axis
    or across it; S8d.1's and S8d.3's `Half`, a whole torus's piece) given
    to a Boolean reaches `polyhedra::build`, since `curved::applies` takes
    none: `stored_model` refuses its curved faces ("a Boolean of a solid
    with curved faces or edges in any position (S9e.4)") and `prism_model`
    a partner's arcs ("a Boolean of a prism with arcs and a solid other
    than a prism (S9e.4)"); a result of one given to another Boolean is
    refused by `given.rs` ("a Boolean's result of a plane's piece given to
    another Boolean (S9e.4)"), one beside a given result by
    `curved::build`. Decisions. (1) *Which.* A split piece where either
    input has a face other than a plane, against a prism, a sphere, a cone,
    a torus, an imported piece (S9e.4b.3a), another split piece or a given
    result, as object or tool, and a result of one given to another Boolean
    (S9e.3a's chain, its leaf the piece's model). A `Clipped` of a prism of
    lines against a polyhedral partner stays S9b.2's stored model (exact on
    planar faces, as now). (2) *The model: S9e.4b.3a's, from the split.* The
    given model of the piece's primitive common the half-spaces of its
    planes, read off its construction (no recognition): a `Clipped`'s
    primitive the prism (its profile on the piece's frame between its
    heights) and its plane the split's `a u + b v + c w + d` in that frame,
    the piece's side by its stored sign; a `Half`'s primitive the cone or
    frustum on its frame, for a zone or cap the whole sphere on its frame
    with its ends' parallels' planes `w = start` and `w = end` (S9e.4b.3a's
    model of a sphere's piece: a plane through the axis then meets the
    sphere through its stored poles as an imported half's does, where the
    cap's own model would take that section through its pole vertex), a
    torus's the whole torus; its plane the one the split built it on (a
    plane within a quarter of the resolution of holding the axis taken
    through it, as S8c.2 builds its halves about the trace; a torus's
    within the resolution of normal to its axis at the band's rounded
    height, S8d.1's), the side by its index (the first below). Each plane
    enters the world exactly through the frame's stored axes, `F(X) = m .
    (X - o) + d` with `m` the coefficients through the inverse of the
    axes' matrix (where the split took them from a plane's frame, its
    stored normal exactly); the hull (`pieces.rs`'s `Hull`, given its planes
    exactly) bounds them by the cube about the primitive, its faces'
    parameters the piece's stored faces' frames. A profile not convex across
    the plane leaves several pieces on a side: the arrangement's several
    solids, the piece the one S9e.2's match finds (its whole construction
    given, `keep_solid`). No new field or degree: the planes' points are
    rationals, their sections the engine's. (3) *Match.* S9e.2's geometric
    match of the arrangement's assembly to the piece's stored topology,
    whose ids name the model, so a Boolean's history is over the piece's
    ids directly (`ComputationLimit("a given result rebuilt differently")`
    where they differ, which a fixture shows first). (4) *Degenerate and
    refused.* S9's rules in both arrangements, unchanged: a cone's or
    frustum's half by a plane through its axis passes through its apex,
    real or virtual, `Degenerate("a plane through a cone's apex")` (S9d.3a's
    rule, as the S9e.4b.3a trial's `buc60926`); faces within the resolution
    of one plane (a box on a piece's cut face); tangencies; a vertex on a
    face; a closed section's own vertex in the model (at its conic's angle
    zero, no stored vertex of the piece's) on the other input's face, a
    seam conflict no seam of the second arrangement moves, `Degenerate("a
    meeting at every seam tried")` as for S9e.4b.3a's pieces. Refused: a
    torus's piece by a plane oblique to its axis (S8d.3's spiric pieces,
    their sections the split's procedural curves, which no rule of S9e.2's
    match takes) `OutOfDomain("a torus's piece by a plane oblique to its
    axis against curved faces (refused, S9e.4b.3b)")`; a spline prism's
    piece `OutOfDomain("a spline prism's split piece against curved faces
    (S9f)")`; a split zone and a sphere, zone or piece of the same sphere
    `OutOfDomain("faces of both inputs on one sphere (S9e.4b.3c)")`
    (S9e.4b.3a's rule). (5) *Its other queries*: unchanged, S8's
    (classification, mass, bounds, rigid motion); the model is built again
    from a moved piece (its key the construction, frame and topology). (6)
    *Fuzzing.* No fuzz target gives a split piece to a Boolean: the
    `split` target's oblique pieces of prisms with arcs and its revolved
    pieces are each given to a Boolean with a box holding the solid (the
    piece's own volume) and a turned box through its centroid (the pair
    identities), a switch of their own (`PIECE_BOOLEANS`). (7) *Evidence
    first.* A row `split ox oy oz nx ny nz xx xy xz below|above` after a
    solid's rows (`identity_reference.Case.split`, both protocols): the
    kernel splits the solid by the plane through the frame's origin normal
    to its normal and takes the one piece on that side; the oracle takes
    the solid common `BRepPrimAPI_MakeHalfSpace` of the plane on that side.
    The reference is the construction through S9e.3a's chained reference, a
    piece `P common H` with `H` a box on the plane's frame on the kept side
    reaching past the solid (`generate_split_pieces_boolean_fixtures.py
    --check`, a CI group `split-pieces`,
    `test_split_pieces_boolean_reference.py`): a cylinder's piece below a
    leaning plane (closed form), a prism of lines and an arc and a box cut
    across both caps, a frustum's and a zone's piece across their walls, a
    zone's half through its axis, a torus's band, against boxes, balls,
    rods, a cone and another piece, as object and as tool, and a chain;
    declared `degenerate` a frustum's half through its axis and a box on a
    piece's cut plane, `unsupported` a zone's half against its own ball. A
    native capture before `curved/splits.rs` exists
    (`compare_split_pieces_boolean.py` keyed on it, the kernel's probe
    refusing every case); then the kernel and its tests
    (`tests/split_piece_booleans.rs`).
  * **S9e.4b.3b evidence (2026-10-04).** The protocols take a split piece:
    after a solid's rows a row `split ox oy oz nx ny nz xx xy xz
    below|above` (`identity_reference.Case.split` in `encode_case` and
    `native_case`; `tests/support/identity_protocol.rs`'s `CaseSpec.split`,
    `build` keeping the one piece `Solid::split_by_plane` leaves on that
    side). `generate_split_pieces_boolean_fixtures.py --check`: 48 cases of
    16 groups (39 solid, 6 declared `degenerate`, 3 `unsupported`; 21 of
    class `prism`, 21 `revolved`, 3 `both`, 3 `chain`) on 8 pieces: a
    cylinder of radius 3 below a leaning plane across its wall
    (`cyl_low`, its closed form), S9e.4a's `dee` below a tilted plane across
    both caps, a box above a leaning plane across both caps (plane faces
    only, against curved partners), a frustum below a tilted plane across its
    wall (an ellipse), a zone on `SKEW2` on one side of the plane through its
    axis (S8c.2's half), a zone on `SKEW` above a tilted plane (a circle
    across its wall), a whole torus above a plane normal to its axis (S8d.1's
    band) and a frustum's half by a plane through its axis; each against
    boxes, balls, rods along the world's `y` and upright, a cone along the
    zone's axis inside its band, as object and as tool, the cylinder's and
    the dee's pieces together, and a chain (the cylinder's piece less a rod,
    then with a level slab across it); declared `degenerate` the frustum's
    half through its axis against a box (its virtual apex) and a box on the
    cylinder piece's cut plane, `unsupported` (S9e.4b.3c) the zone's half
    against its own whole ball. Each split is checked to be the kernel's
    piece of its class (`split_class`, exactly on the stored axes: a
    prism's plane oblique to its axis, a revolved solid's within a quarter
    of the resolution of its axis or not, a torus's normal to its axis).
    The reference is the construction OCCT is given through S9e.3a's
    chained reference, each piece `P common H`, `H` a box on the plane's
    frame on the kept side reaching past the solid (each solid within four
    fifths of its reach), a zone's piece its whole sphere common `H` common
    the slab between its ends' parallels at the heights the kernel stores
    (`r sin(lat)`, the sine and the product each rounded once: the chained
    reference takes whole spheres): the two families within 6.2e-36 of the
    case's size, each solid's closed form 9.7e-36, the cylinder's piece's
    closed form 1.9e-41 (its volume the disc's area times the plane's height
    over the axis, its area the disc, the wall and the ellipse `pi r^2 |m| /
    |m_z|`, `m` its box's face's normal), the pair identities on the last
    Boolean's arguments 2.6e-38 and the area identity 3.8e-41, Monte Carlo
    2.8 standard errors (50,000 points a group), quadrature estimates
    1.9e-32, solid counts by rays at two resolutions and each piece one
    solid, every meeting's sine at least 0.083 and events at least 1.3e-5 of
    their range apart outside the declared groups.
    `test_split_pieces_boolean_reference.py` checks the cylinder's piece's
    closed form against a direct integration (a plane normal to the axis
    the cylinder's own), the chained reference on the piece against it, the
    splits' classes (a plane normal to a prism's axis, off a zone's centre
    or oblique to a torus's axis not that class), the half-space boxes and
    the case list and its protocol rows. The generator's check is a CI group
    of its own (`split-pieces`, 15 to 27 minutes on four workers locally); Python
    3.9 and 3.12 write the same files. Corrections from the evidence,
    amending the refined decisions' plan (7): (a) the chained reference
    takes whole spheres only, so a zone's piece is its sphere common two
    boxes (above); (b) half-space boxes reaching far past the solids made a
    group take up to 50 minutes (the chain's), so each box reaches just past
    its piece's solid and the chain's slab is level, and Monte Carlo takes
    50,000 points a group; (c) partners moved where a ball's pole touched a
    frustum's cap's plane, the chain's slab's plane met the cut plane on the
    cylinder's wall, a ray count's resolution split a thin neck (a box
    lowered, a cone along the zone's axis, the band's torus on the world's
    axes and its ball smaller), and the pair is the cylinder's and the dee's
    pieces (a zone's half beside them took the longest of all). The oracle
    takes the `split` row as the solid common `BRepPrimAPI_MakeHalfSpace` of
    the plane's face on the kept side (its reference point the origin a unit
    along or against the normal). The capture
    `occt-boolean-split-pieces-preimplementation`
    (`compare_split_pieces_boolean.py`, keyed on
    `solid/boolean/curved/splits.rs`, S9e.4a's slack; the kernel's probe
    `unsupported` on all 48, every split piece refused as S9e.4's): every
    result valid with the reference's solid count, 33 matching, 15 reviewed
    (BRepGProp's default integration on faces bounded by the approximated
    sections of the ball and the rods with the cylinder's piece, the ball
    with the frustum's, the rod with the zone's, the ball with the band
    and the chain's common: volumes up to 3.6e-6 relative off, within
    1.4e-8 measured adaptively at 1e-10 and at 1e-12 by a diagnostic
    build); 20 results' counts change when unified. S9e.4b.3b's kernel
    next.
  * **S9e.4b.3b implemented** (`solid/boolean/curved/splits.rs`,
    `curved::applies` and `model_of`, `pieces.rs`'s hull of exact planes,
    `given.rs`'s chains, `split::Clipped::parts`, `split::Half::parts` and
    `built_plane`), as the refined decisions describe: a split piece where
    either input has a curved face is the given model of its split's
    primitive common the half-spaces of its planes, read off the split (a
    prism's oblique piece its prism and plane; a cone's piece its cone, a
    zone's or cap's its whole sphere and its ends' parallels' planes, a
    torus's band its whole torus, each with the plane it was built on),
    each plane taken exactly into the world through the frame's axes, the
    hull's faces on the piece's stored faces' frames, matched to the piece's
    stored topology, the history over its ids. All 48 fixtures as declared
    (39 within the kernel's enclosures, each at most `1e-9` wide;
    `cone_axis` refused as `Degenerate("a plane through a cone's apex")`,
    `cyl_flush` as two faces within the resolution of one plane,
    `one_sphere` `OutOfDomain("faces of both inputs on one sphere
    (S9e.4b.3c)")`), every history complete over the pieces' ids and none
    from a primitive or a hull, results deterministic and moved rigidly,
    both inputs translated and turned keeping the reference's volumes (a
    plane exactly along a rod's axis, `cyl_rod`'s and `block_rod`'s, is
    within rounding of it once turned: S9's `Degenerate`, so those are not
    moved), the stored frames and the splits' planes the reference's bit for
    bit, every fixture's piece its vertices on its boundary and its common
    with a box holding it itself; beside the fixtures, a U profile's two
    pieces above one plane (each one solid of several of its construction,
    S9e.2's sorting) against a rod, a hemisphere's halves through its pole
    against a box (their commons the cap's), a cylinder's piece against
    S9e.4a's imported octant and a given result as object and tool, a
    torus's spiric piece refused (`tests/split_piece_booleans.rs`, 29 s at
    `opt-level` 2 with debug assertions, 31 s in release with debug
    assertions under the emulated correctly rounded `hypot` on a host at
    load 10 to 20). `compare_split_pieces_boolean.py` 26 matches and 22
    reviewed (the 15 captured; with the kernel, entity counts in 14, the same
    faces as OCCT's unified result, each splitting its sections at its own
    points and seams), every enclosure within the reference with S9e.4a's
    `1e-12` slack; every other comparison unchanged
    (`compare_imported_pieces_boolean.py` 32/13,
    `compare_imported_boolean.py` 54/15, `compare_imported_arcs_boolean.py`
    27/9, `compare_imported_polyhedra_boolean.py` 47/1, given 36/0, given
    curved 25/23, chained 24/6), the suite and the tools' unit tests (323)
    passing. Amendments, from the implementation and the split target's new
    stage (its replay found the last three, each an imported piece's too):
    (a) a zone's or cap's piece is its whole sphere with its ends' planes
    (S9e.4b.3a's model), not the cap: a plane through the cap's axis meets
    it through its stored pole, which the cap's own model took as a seam
    conflict at every seam; (b) a loop through a stored sphere's pole that
    winds none (S9e.4b.3a's rule) but does not start there has its fins from
    the pole on lifted by the turn it made, and the face's other loops
    lifted again to lie with it (`assemble.rs`'s `pole_lift`: a cap's half
    against a box across the seam, the last pcurve not closing on the
    first, `uv_gap`, or a hole on the wrong sheet, `inner_loop_outside`);
    (c) a stored circle matched as a ring runs the result's way, turned over
    where it runs the other (`given_arc`: a zone's rim as a parallel's plane's
    section of the whole sphere, that plane's normal out of the piece,
    `loop_winding`); (d) a cone's section stored as a ring runs the result's
    way about the cone's axis (`ring_about`: a frustum's oblique piece, its
    wall's loops winding the same way, "a cone's wall winding without an
    apex"); (e) a closed section's own vertex in a model (at its conic's
    angle zero) on the other input's face is a seam conflict no seam of the
    second arrangement moves: `Degenerate("a meeting at every seam
    tried")`, as for S9e.4b.3a's pieces (a box's face through the frustum's
    section's angle-zero point; no fixture). The `split` fuzz target gives
    each oblique piece of a prism with arcs and each piece of a cone, zone
    or cap by a plane not normal to its axis to Booleans
    (`PIECE_BOOLEANS`, on): its common with a box holding it is itself, and
    the first piece's fuse, cut and common with a box turned about its
    centroid obey the pair identities (a torus's pieces are not given, their
    arrangements taking seconds; the second piece not turned, for time: the
    corpus's slowest input took 62 s under AddressSanitizer with both; with
    the first only the five slowest others take 13 to 18 s and an existing
    slow regression 28 s, at load 10 to 19). Replaying the split corpus and
    its regressions (3,563 inputs; two failing inputs now regressions) with
    debug assertions, no failure: 890 inputs reach the stage with 1,718 pieces,
    1,404 commons with the holding box and the identities of 699 of the 868
    turned boxes evaluate, the rest refused as documented (a vertex of one
    input on the other's face, a spline prism's piece, S9f's, a plane
    through or within the resolution of a cone's apex, meetings within
    rounding, pieces or results thinner than the resolution, tangencies),
    the slowest input 3.8 s on a host at load 10 to 20; the boolean
    target's corpus and regressions (1,459 inputs) replayed with debug
    assertions, no failure, the slowest 15.6 s. Pending: the capture's Linux
    record. DRAW survey: that of S9e.4b.3a, S9e.4b.3b and the fuzz fixes,
    below (no case's status or reason moves with it). Campaigns: at `6d26886d` the `boolean`
    target's clean (886 runs, the slowest input 28 s under AddressSanitizer
    at load about 20), the `split` target's found a torus cap's round end
    lifted the wrong way (S8d.3's, latent since `d7049aac`) and, after its
    fix, a plane within rounding of a frustum's virtual apex reaching the
    arrangement through `PIECE_BOOLEANS` (both fixed, notes below); at
    `cc7ea7ff`, with those fixes and `CONE_PAIRS` on, both clean: `boolean`
    996 runs, the slowest 13 s, `split` 1,731 runs, the slowest 19 s.
  * **DRAW survey of S9e.4b.3a, S9e.4b.3b and the fuzz fixes (2026-10-05,
    `UPSTREAM_TESTS.md`).** At `d665df29` (`s9c2-kernel` with S9e.4b.3a,
    S9e.4b.3b, the loops' certified integrals, the validator's curve points
    guarded, the degree-eight arithmetic and the sheared projections, the
    slowest inputs' speedup, the prism walls' frames (`Frame3::at`), the
    loop through one pole and S9e.4b.3b's pole lift, the three cone-pair
    fixes, the narrower `Meet` and `Toric` integrals, a torus cap's round
    ends and a plane within rounding of a cone's apex; the public dataset,
    120 seconds a case, four at once). The 1,802 self-contained cases of the
    Boolean group on both backends: every status, every refusal's reason
    and every error the last survey's (`12b6c176`) field for field, 987
    evaluating and registered, 592 refused, 223 unsupported on both, none
    failing or timing out, the sentinels, `bopfuse_simple/ZP6` and the
    `gdml_public` tori refused as before: the new guards (a ring on a
    floor, a plane within rounding of a cone's apex) refuse no case that
    evaluated, and no case meets the "points not separated by a
    projection" limits, as before. The
    1,814 cases restoring a shape for a Boolean on the Rust adapter, and
    the 171 the import reaches on native DRAW too: no status moves, and 41
    refusals' reasons move, each from "an imported solid other than a
    prism, a sphere, a cone or a torus (S9e.4b)", native DRAW evaluating
    every one (`viewer_skipped`). S9e.4b.3a's 39 are refused as its trial
    found: 18 (every case restoring `so1` or `so4`) as a piece other than
    its primitive common its planes' half-spaces, 20 (pairs of `so2`,
    `so3`, `so5`, `so6`, `so7`) as faces of both inputs on one sphere, both
    S9e.4b.3c's, and `bugs/modalg_1/buc60926` as a plane through a cone's
    apex. Two CTO parts beyond the trial now reach the pieces' recognition:
    `bcut_complex/G4` (the part a box with a cylindrical boss, radius 52.93
    over 100, one cylinder face and plane faces: refused as a piece other
    than its primitive common its planes, S9e.4b.3c's) and
    `bcut_complex/I6` (the tool a block less a half cylinder, radius 35,
    whose wall is tangent to the block's own faces `y = 0` and `y = 70`:
    `Degenerate("a tangency between the inputs (S9c)")`, raised in the
    piece's own first arrangement, the cylinder primitive against its hull,
    not between the case's inputs: the tool alone against a box 1,000 away
    is refused alike; the reason names the wrong pair, the refusal stands).
    Of the 171: 27 evaluate on both backends, 58 are refused by S9's rules,
    39 as S9e.4b.3c's, 36 are bodies none of the kernel's constructions, 6
    S9e.4b.4's, 3 spline bodies (S9f) and 2 arguments of several solids;
    24 are unsupported natively too. The same 8 time out at their first
    restores. S9e.4b.3b's split pieces move no status or reason. The
    volume audit (`vprops` and `sprops` before each `checkprops`): native
    DRAW's values the last audit's bit for bit on all 1,014 registered
    cases; the Rust adapter's statuses the same, and 115 cases' values
    moved, 49 volumes, 34 areas and 250 centre coordinates, every one
    within rounding (volumes and areas at most 6.6e-16 relative, a centre
    at most 1.6e-15 of the solid's size): 113 cases at the loops' certified
    integrals' merge (`cfeab65d`: 30 volumes, 2 areas, 218 coordinates)
    and 60 at the off switches' merge, its narrower `Meet` integrals
    (`f87d150d`: 36 volumes, 33 areas, 154 coordinates; two cylinders or a
    cylinder and a cone, `bop*_simple/ZE3` to `ZE6`, `ZI8` to `ZJ3`, `ZK5`
    to `ZL1`), 58 in both (a Rust audit at `770bcdbc`, and the merges
    between on five cases, locate them). None crosses native DRAW's printed
    digits: 81 of the 115 agree to them as before, the other 34 are among
    the same 35 disagreements, each flagging the same fields, native off in
    each. No case evaluates newly, so none is registered (1,100 cases) and
    the ledger does not change. A full contract run of the manifest holds
    on both backends with the dataset (30 seconds a case), the slowest
    Boolean case 6.0 seconds (`bopfuse_simple/ZK8`, on a host at load 5 to
    10; `bopcommon_simple/ZK8` 16.9 in the last survey), the restore cases
    0.1 to 2.6 s, the rollex 0.4 to 0.7. No case fails, crashes or panics;
    no kernel change.
  * **S9e.4b.3c refined, before its code (2026-10-05).** Why the survey's
    38 restore cases of `so1` to `so7` are refused today, read off the
    files by the converter (the dataset is not committed): every one is a
    Boolean of two pieces of the sphere of radius 10 about the origin, its
    stored frame's origin and radius bit for bit alike in all seven files.
    The 20 pairs of `so2`, `so3`, `so5`, `so6` and `so7` reach
    `pieces::one_sphere` (`OutOfDomain("faces of both inputs on one sphere
    (S9e.4b.3c)")`, ahead of S9d.2's `Degenerate("two spheres about one
    centre")` in `sphere_sphere` and `circ_sphere`); the 18 with `so1` or
    `so4` are refused on import (`OutOfDomain("an imported plane piece
    other than its primitive common its planes' half-spaces (S9e.4b.3c)")`):
    OCCT split each one's rim into two half circles at the stored sphere's
    seam and opposite it (vertices at `x = +-10`, `+-8.66`), where the
    re-run's rim is one ring, so S9e.2's match by counts fails. What each
    pair needs beyond that: `so1` and `so4` (8 cases: `bcommon_complex/B2`,
    `bcut_complex/C5`, `C6`, `bfuse_complex/B5`, `bugs/modalg_2/bug413_1`,
    `_2`, `bugs/moddata_1/bug183_2`, `_3`) and `so4` and `so2` (5: `B5`,
    `D2`, `D3`, `B8`, `bug183_5`) meet in general position on the sphere:
    their circles cross or nest, no vertex of either lies on the other's
    circles or vertices, the wedge's pole vertex inside the cap's sphere
    face (at its stored pole); `so1` and `so2` (5), `so2` and `so3` (4) and
    `so5` and `so2` (4) are exact incidences on one sphere and one frame:
    the wedges' equator arcs on the hemisphere's rim circle, both wedges'
    base faces on `z = 0`, their corners at the centre and their pole
    vertices one point, `so5`'s and `so2`'s axis edges along one line; the
    12 with `so6` or `so7` (OCCT's turned copies of `so2` and `so3`) put
    each wedge's corner within `1e-15` of the other's at the centre (three
    rounded planes' rational point, not on it), so slivers thinner than the
    resolution. Decisions. (1) *Sub-steps*, by what unlocks the most cases
    and by the machinery each takes: **S9e.4b.3c.1** (this step): faces of
    both inputs on one sphere in general position, and an imported sphere
    piece's rim split by stored vertices (13 cases); **S9e.4b.3c.2**: exact
    incidences on one sphere, a vertex of both (a corner at the centre, a
    pole), a circle of both (an arc on the other's rim), a line of both
    (two wedges' axis edges) and plane faces on one plane with overlapping
    edges, the arrangement's vertices merged across the inputs and an edge
    of both one arrangement edge of both inputs' faces, as S9a's polygons
    on one plane (13 cases); **S9e.4b.3c.3**: a piece not convex in its
    planes (`shading_132`, the bitten ball) and two plane faces on one
    plane. The 12 with `so6` or `so7` stay refused by S9's rules (a result
    or piece thinner than the resolution, two faces within the resolution
    of one plane): a corner of each within rounding of the other's is not
    one vertex. (2) *Faces on one sphere: S9c.1's faces on one surface.*
    Two inputs' spheres are one when their exact models' centres and radii
    are equal rationals (the stored frames' origins and radii); their
    sphere faces are then faces on one surface as equal cylinders' and
    coplanar planes' are (`Arr::coinc`): each holds the other's edges that
    lie inside it, a piece of either is classified by the other's
    membership pushed off the sphere both ways, and the kept pieces of both
    facing one way join (the assembly unchanged). Each hemisphere keeps its
    stereographic chart (S9d.1), where both inputs' edges are traced; no
    new field or degree. What the arrangement lacks for spheres: where two
    circles of the inputs cross on the sphere. A rim crossing the other's
    rim is found as now (its pierce of the other's plane face, on that
    face's boundary), but a whole sphere's split great circle bounds no
    plane face, so every circle of one input on a face on the sphere is
    met with every such circle of the other: the first circle's points in
    the second's plane (`Circ::meet_plane`: its place a quadratic surd, the
    point in one quadratic field), a vertex of both where strictly inside
    both edges (`VKey::Circles`, named by both edges). Two pieces whose
    primitives are whole spheres would split them along one great circle
    (each piece's own first seam), so a piece's or a split zone's model
    under a partner on its sphere is built from the second arrangement's
    seam first, and the two great circles differ and move with its retried
    seams. A sphere, cap or zone (S9d.1) and an imported or split piece of
    its sphere are taken alike (any two inputs on one sphere). (3) *Rims
    split by stored vertices.* An imported sphere piece's stored vertex
    within the resolution of a section circle's edge of its own first
    arrangement, strictly inside it, and of none of its vertices splits
    that edge at the circle's point in the direction of the vertex's place
    rounded once (`Circ::at`, scaled onto the circle: a point of one
    quadratic field within the resolution of the stored one), a vertex the
    assembly keeps, as a pole (`VKey::Stored`); the re-run's rim is then the
    stored arcs and S9e.2's match takes it, the history over the stored
    ids. A stored vertex within the resolution of two edges is `Degenerate`.
    (4) *Degenerate and refused.* S9's rules unchanged off the sphere;
    circles of the inputs tangent on the sphere (`meet_plane`'s double
    root) and spheres within the resolution of one and not one
    (`Degenerate("two spheres within the resolution of one sphere")`) are
    `Degenerate`; concentric spheres of different radii stay S9d.2's
    `Degenerate("two spheres about one centre")`. An exact incidence on one
    sphere, a circle of both or a vertex of one at an end of the other's
    circle, is `OutOfDomain("a vertex or a circle of both inputs on one
    sphere (S9e.4b.3c.2)")`, found before the pierces so that a split great
    circle in it (retried at another seam) never hides it; S9e.4b.3c.3's
    refusals unchanged. Caps and zones on one frame meet at their stored
    poles (both splits hold the axis): S9e.4b.3c.2's. (5) *Its other
    queries* unchanged. S9e.4b.3a's `one_sphere` (the octant and
    `octant_low`) and S9e.4b.3b's (the zone's half and its own whole ball)
    become solid cases. (6) *Fuzzing.* The `split` target's
    `PIECE_BOOLEANS` gives a zone's or cap's piece also to a whole ball on
    its sphere in a turned frame: their common is the piece, their fuse the
    ball. (7) *Evidence first.* Bodies OCCT writes on one sphere in turned
    rational frames (the reference's families are the world's meridians):
    the hemisphere and the cap above a parallel's plane with their rims
    split as `so1`'s and `so4`'s (a `write` block's row `divide`,
    `ShapeUpgrade_ShapeDivideClosedEdges` on the written solid: each closed
    edge divided in two), and S9e.4b.3a's octant (`so2`'s wedge) in the
    cap's frame; cases of two pieces of one sphere (the hemisphere and the
    cap either way, the cap and the octant, a ball and the cap, the cap and
    a wedge above a parallel's plane in another frame, a chain), declared
    `degenerate` a ball within the resolution of the cap's sphere, declared
    `unsupported` (S9e.4b.3c.2) the hemisphere and the octant on one frame
    and two octants about one centre; S9e.3a's chained reference with each
    piece's closed form, its shared surfaces leaving the area identity
    (`generate_one_sphere_boolean_fixtures.py --check`, a CI group
    `one-sphere`, `test_one_sphere_boolean_reference.py`); a native capture
    keyed on the refusal it removes (`compare_one_sphere_boolean.py`); then
    the kernel and its tests (`tests/one_sphere_booleans.rs`).
  * **S9e.4b.3c.1 evidence (2026-10-05).** Five bodies OCCT wrote
    (`boolean-one-sphere-bodies.txt`, `write` blocks of a sphere's row,
    `boolean common` and a box's or prism's rows, written by
    `compare_one_sphere_boolean.py --write-bodies` to
    `rust/fixtures/imported/`), every one a piece of S9e.4b.3a's octant's
    sphere (radius 5 about `(5, 5, 4)`): `sphere_hemi`, the ball above its
    equator's plane on `SKEW4`, and `sphere_cap`, above the parallel's
    plane at `3/2` on `SKEW4`, each with its rim divided as `so1`'s and
    `so4`'s (a `write` block's new row `divide`:
    `ShapeUpgrade_ShapeDivideClosedEdges` with one split point on the
    written solid, the rim two arcs between the stored seam's point and the
    opposite one); `sphere_octant2` and `sphere_octant4`, the corners of
    `SKEW2`'s and `SKEW4`'s axes (`so2`'s wedges); `sphere_wedge`, above the
    parallel's plane at `1` between two meridian planes of `SKEW` (`so5`'s).
    `generate_one_sphere_boolean_fixtures.py --check`: 22 cases of 9 groups
    (16 solid, 1 empty, 2 declared `degenerate`, 3 `unsupported`; 14 of
    class `pieces`, 5 `sphere`, 3 `chain`): the hemisphere and the cap
    either way (`so1` and `so4`: the cap's cut by the hemisphere empty), the
    cap and the octant either way (`so4` and `so2`), a ball of the sphere on
    the world's axes and the cap, the wedge and the cap, a chain (the cap
    less the octant, then with a level slab: its cut two solids); declared
    `degenerate` a ball of radius `5 + 2^-30` about the centre against the
    cap (its fuse and common); declared `unsupported` (S9e.4b.3c.2) the
    hemisphere and `sphere_octant4` on one frame (the octant's equator arc
    on the hemisphere's rim circle). The reference is the constructions
    OCCT was given through S9e.3a's chained reference (two pieces `(P common
    B) op (Q common D)`, `P` and `Q` one sphere, its faces of two inputs on
    one surface), with S9e.4b.3a's checks relative to the case's size: the
    two families within 3.4e-36, each solid's closed form 4.1e-36, each
    piece's closed form 1.7e-34 (`cap_closed`, the hemisphere, the
    corners), the pair identities 5.7e-38 (the area identity left: shared
    surfaces), Monte Carlo 2.4 standard errors (100,000 points a group),
    quadrature estimates 9.2e-33, every meeting's sine at least 0.35 and
    events at least 1.3e-6 of their range apart outside the declared
    groups, each piece one solid and solid counts by rays at two
    resolutions, here with neighbouring rays' intervals joined within two
    grid spacings (the generator's `count_solids`: two pieces of one sphere
    leave slivers tapering to it, whose neighbouring rays' intervals are
    disjoint where a sliver is thinner than its slope across a spacing, and
    the chained reference's strict overlap counted a wedge less the cap as
    up to 20 solids). Every file read independently (`stored_records`): its
    faces one sphere and planes, every stored vertex within 5.9e-16 of the
    size on the construction's surfaces, a divided rim two stored vertices
    on the rim's circle. `test_one_sphere_boolean_reference.py` checks the
    cap's closed form against an integration of its slices, the chained
    reference on the hemisphere and the cap of one sphere (their common the
    cap, their fuse the hemisphere), the frames' normals off the kernel's
    whole sphere's split planes, the case list and its protocol rows. The
    generator's check is a CI group of its own (`one-sphere`, 14 minutes on
    four workers locally under Python 3.12); Python 3.9 and 3.12 write the
    same files. Corrections from the evidence, amending the refined
    decisions' plan (7): (a) `SKEW3`'s normal `(10, 11, 2)` lies in the
    split plane of the kernel's whole sphere's fourth rotation (rows `(1, 2,
    2)`, `(2, 1, -2)`), so a piece on it, split there at the second
    arrangement's seam, had its pole vertex on its own great circle within
    rounding (`Degenerate("two meetings within rounding along an arc")`):
    the bodies are on `SKEW4`, `SKEW2` and `SKEW`; (b) in `so4` and `so2`
    the wedge's pole vertex is the cap's stored pole, exactly on the world's
    axes; in turned rational frames two planes' line meets the sphere
    within rounding of a stored pole, not on it (the reference's events
    9.7e-19 of their range apart), so the octant's axis is another frame's
    and the coincidence is the DRAW trial's alone; (c) the octants of two
    frames about one centre share their corner (a vertex of both off the
    sphere, not found as an incidence on it) and are not declared; (d) a
    cap at `5/2` and a tilted chain slab left counts the reference's rays
    did not resolve, so the cap is at `3/2` and the slab level; (e) the
    ball within the resolution has no cut case (a shell `2^-30` thick).
    The capture `occt-boolean-one-sphere-preimplementation`
    (`compare_one_sphere_boolean.py`, keyed on the refusal the step removes,
    `pieces.rs`'s `OutOfDomain("faces of both inputs on one sphere
    (S9e.4b.3c)")`; the kernel's probe `unsupported` on all 22, the divided
    rims refused on import): 8 matching, 14 reviewed. OCCT's General Fuse
    on faces of both inputs on one sphere leaves 10 results invalid under
    BRepCheck_Analyzer (the hemisphere and the cap's fuse, the cap and the
    octant's fuse and cut, the ball and the cap's fuse and common, the wedge
    and the cap's fuse, the near ball's two, the declared `hemi_octant`'s
    cut and common), their solid counts the reference's and their measures
    within 1.1e-5 relative (BRepGProp's default integration; within 7.1e-12
    measured adaptively at 1e-10 and 1e-12 by a diagnostic build, the near
    ball's common 5.2e-10); the cap and the octant's common is valid, its
    volume 2.4e-5 relative off by the default integration (1.1e-11
    adaptively); and the chain is wrong: its first result (the cap less the
    octant, invalid) given to the slab, OCCT's fuse is 3.2e-2 short, its cut
    one solid of the reference's two and its common empty, consistent with
    that result missing its part above the slab's bottom plane.
    S9e.4b.3c.1's kernel next.
  * **S9e.4b.3c.1 implemented** (`curved/graph.rs`: faces on one sphere,
    the circles' crossings and `Arr::split_at`; `curved/pieces.rs` and
    `splits.rs`: a piece's sphere split at the second arrangement's seam;
    `meet.rs`, `spheres.rs`: one sphere a surface of both; `assemble.rs`:
    the new vertices' names, a stored vertex kept), as the refined
    decisions describe: two inputs' spheres one when their exact models'
    centres and radii are equal (`pieces::on_one_sphere`), their sphere
    faces then faces on one surface (`Arr::coinc`), every circle of one
    input on such a face met with every circle of the other there (a
    vertex of both, `VKey::Circles`, strictly inside both edges), a piece's
    or split zone's sphere split at the second arrangement's seam first
    (`model_of`'s `one_sphere`), and an imported sphere piece's section
    circle split at each stored vertex within the resolution of it and of
    no vertex (`VKey::Stored`, kept by the assembly), so the divided rims
    match. All 22 fixtures as declared (17 within the kernel's enclosures,
    each at most `1e-9` wide; `ball_near` refused as `Degenerate("two
    spheres within the resolution of one sphere")`, `hemi_octant` as
    `OutOfDomain("a vertex or a circle of both inputs on one sphere
    (S9e.4b.3c.2)")`), every history complete over the imported bodies'
    stored ids, results deterministic and moved rigidly, both inputs
    translated and turned keeping the reference's volumes, the divided
    rims' bodies pieces of two faces, two edges and three vertices (their
    volumes the closed forms), two whole balls of one sphere in two frames
    (their fuse and common the ball, their cut empty)
    (`tests/one_sphere_booleans.rs`, 11 s at `opt-level` 2 with debug
    assertions on a host at load 15, 11 s in release with debug assertions
    under the emulated correctly rounded `hypot`, the pieces' two earlier
    test files too). S9e.4b.3a's `one_sphere`
    (the octant and `octant_low`) and S9e.4b.3b's (the zone's half and its
    own ball), declared `unsupported` until now, are solid within the
    reference (their sets' generators and unit tests declare them so;
    `compare_imported_pieces_boolean.py` 32/13 and
    `compare_split_pieces_boolean.py` 26/22 unchanged).
    `compare_one_sphere_boolean.py` 8 matches and 14 reviewed (the 14
    captured; with the kernel, entity counts in 7: OCCT's unified results
    split their sphere faces and sections at its own points and seams where
    the kernel joins the pieces of both inputs facing one way, and the
    chain's are OCCT's wrong results'), every enclosure within the
    reference; every other comparison of `HANDOFF.md`'s table unchanged with
    0 failures; the release suite (601 tests) and the tools' unit tests
    (327) passing, the ledger unchanged. Amendments, from the
    implementation: (a) a piece's model split at the second arrangement's
    seam and degenerate there (`SKEW3`'s pieces, evidence correction (a)) is
    split at the other seams in turn, as on import (`common`'s `first`): a
    piece degenerate at every split keeps its own refusal, and a split
    alike the partner's is the second arrangement's seam conflict (two
    split great circles on one circle), retried at its next;
    (b) the circles' crossings are found before the pierces, an incidence
    of the inputs' own circles and vertices (`OutOfDomain`, S9e.4b.3c.2)
    taking precedence over a split great circle's seam conflict, which
    otherwise hid it at every seam (the survey's `so1` and `so2`); (c) a body
    of one sphere, cylinder or cone face and planes whose curved face's
    material lies outside its quadric (a groove or notch, the face's
    outward normal toward the axis or centre) is refused on import as
    S9e.4b.3c's before its model is built (`imported::piece`; the survey's
    `bcut_complex/I6` tool, a block less a half cylinder whose wall is
    tangent to two of its faces, was refused as `Degenerate("a tangency
    between the inputs (S9c)")` from its own model's arrangement; the DRAW
    survey of branch `s9-draw-6` found it). The `split` fuzz target's
    `PIECE_BOOLEANS` gives a zone's or cap's first piece also to a whole
    ball of its sphere in a turned frame (their common the piece, their fuse
    the ball, within `1e-9` of the ball's volume): of the split corpus and
    its regressions 265 inputs reach it, the common and the fuse of 260
    evaluate, 5 refused by S9's rules (a piece or result thinner than the
    resolution, a meeting at every seam tried). Replays with debug
    assertions, one process an input: the split corpus and its regressions
    (3,569 inputs) and the boolean corpus and its regressions (1,466), no
    failure, the slowest 2.4 s and 6.6 s. A trial of the DRAW survey's 38 restore cases of `so1` to `so7`
    on the Rust adapter (not the survey: nothing registered): 14 evaluate on
    both backends with every check (`bcommon_complex/B2`, `B5`, `B9`,
    `bcut_complex/C5`, `C6`, `D2`, `D3`, `bfuse_complex/B5`, `B8`,
    `bugs/modalg_2/bug413_1`, `_2`, `bugs/moddata_1/bug183_2`, `_3`, `_5`:
    `so1` and `so4`, `so4` and `so2`, and `so2` and `so6`'s common, whose
    corners within rounding at the centre leave its common clear of them);
    13 are S9e.4b.3c.2's (`OutOfDomain`: `so1` and `so2`, `so2` and `so3`,
    `so5` and `so2`); the other 11 with `so6` or `so7` `Degenerate` by S9's
    rules (two faces within the resolution of one plane, a piece or result
    thinner than the resolution); `bcut_complex/I6` and `G4` refused as
    S9e.4b.3c's (a groove; a box with a cylindrical boss), I6 before its
    model's tangency. Pending: the capture's Linux record. DRAW survey:
    that of S9e.4b.3c.1 and the switches' speed-up, below (the trial's 14
    registered, every other `so` case refused as it found).
    Campaigns at `59d0c57b`, with every switch of the boolean target on
    (`TORUS_PAIRS`, `CONE_PAIRS`, `TURNED_PARTS`, `GIVEN_MET`,
    `SPLINE_SPHERE`, `SPLINE_CONE`) and the caps' frames kept bit for bit
    (600 s each, sampled replays): `boolean` clean, 1,080 runs, the slowest
    input 13 s under AddressSanitizer at load about 6; `split` clean, 1,818
    runs, none slow.
  * **DRAW survey of S9e.4b.3c.1 and the switches' speed-up (2026-10-05,
    `UPSTREAM_TESTS.md`).** At `f4b584f7` (`s9c2-kernel` with S9e.4b.3c.1,
    the speed-up of `TORUS_PAIRS` and `TURNED_PARTS` (memos of a circle's
    resultants with a torus and of isolated roots, the shifted root
    refinement, the fields' signs, a torus meeting's jets shared by its two
    pcurves, a containment ray undecided at once beside an undecidable
    fin) and a prism's caps on its frame bit for bit (`59d0c57b`); the
    public dataset, 120 seconds a case, four at once). The 1,802
    self-contained cases of the Boolean group on both backends: every
    status, every refusal's reason and every error the last survey's
    (`d665df29`) field for field, 987 evaluating and registered, 592
    refused, 223 unsupported on both, none failing or timing out. The 1,814
    cases restoring a shape for a Boolean on the Rust adapter, and the 171
    the import reaches on native DRAW too: native DRAW's every status and
    reason the last survey's; on the Rust adapter only S9e.4b.3c.1's 38
    `so` cases and `bcut_complex/I6` move, each as its trial found, besides
    two restores that no longer time out (below). 14
    evaluate on both backends with every check (`so1` and `so4`, the
    hemisphere and the cap above `z = 5` with their rims divided:
    `bcommon_complex/B2`, `bcut_complex/C5`, `C6`, `bfuse_complex/B5`,
    `bugs/modalg_2/bug413_1`, `_2`, `bugs/moddata_1/bug183_2`, `_3`; `so4`
    and `so2`, the cap and a wedge: `bcommon_complex/B5`, `bcut_complex/D2`,
    `D3`, `bfuse_complex/B8`, `bug183_5`; `so2` and `so6`'s common,
    `bcommon_complex/B9`), 13 of them from "a piece other than its
    primitive common its planes" and `B9` from "faces of both inputs on
    one sphere". 13 are S9e.4b.3c.2's (`OutOfDomain("a vertex or a circle
    of both inputs on one sphere (S9e.4b.3c.2)")`: `so1` and `so2`, `so2`
    and `so3`, `so5` and `so2`; 5 from the former reason, 8 from the
    latter) and 11 with `so6` or `so7` are `Degenerate` by S9's rules (4
    two faces within the resolution of one plane, 4 a piece and 3 a result
    thinner than the resolution). `bcut_complex/I6`'s notched tool is now
    refused on import as a groove (S9e.4b.3c's reason for a plane piece
    other than its primitive common its planes, where the last survey found
    its own arrangement's tangency), `G4` as before. Of the 171: 41
    evaluate on both backends, 68 are refused by S9's rules, 15 are
    S9e.4b.3c's (13 S9e.4b.3c.2's, `G4` and `I6`), 36 are bodies none of
    the kernel's constructions, 6 S9e.4b.4's, 3 spline bodies (S9f) and 2
    arguments of several solids; 24 are unsupported natively too. Of the 8
    that timed out at their first restores, 6 still do; `bugs/modalg_1/
    buc60532_2` (a `SurfaceOfLinearExtrusion` the reader does not represent)
    and `bugs/modalg_6/bug23585` (`tolerance p`) ended within the 120 s, in
    113 and 102 s, as unsupported: load, not the kernel (alone, the worker
    of `d665df29` took 77 to 119 s on them and this one 78 to 116 s,
    interleaved at load 5 to 13). The volume audit (`vprops` and `sprops`
    before each `checkprops`, both backends, the 1,014 registered and the 14
    new cases): native DRAW's values the last audit's bit for bit on all
    1,014; the Rust adapter's statuses the same and 2 cases' values moved,
    `bfuse_complex/F5` and `Q2` (one pair of CTO prisms, `CTO900_pro10658a`
    and `pro10658b`, fused): 4 centre coordinates, at most 4.7e-16 relative
    (4.6e-16 of the solid's size), no volume or area. They move at
    `59d0c57b`, the caps' frames (a Rust audit with the worker of its
    parent `ff8c914f` gives the last audit's bits); S9e.4b.3c.1 and the
    speed-up move no audited value. Both still agree with native DRAW's
    printed digits; the same 35 disagreements as before. The 14 new cases'
    volume, area and centre agree with native DRAW's printed digits
    (`C6`'s cut empty on both). They are registered (`data`,
    `viewer_skipped` on both; 1,114 cases), and the ledger does not change
    (`--ledger` holds). A full contract run of the manifest holds on both
    backends with the dataset (30 seconds a case), the 14 in 0.6 to 3.3 s
    on the Rust adapter (`B9` the slowest), the slowest Boolean case 5.3 s
    (`boptuc_simple/ZK8`, on a host at load 7 to 9; `bopfuse_simple/ZK8`
    6.0 in the last survey), the restore cases 0.2 to 4.1 s
    (`bfuse_complex/N9`), the rollex 0.6 to 0.7. No case fails, crashes or
    panics; no kernel change.
  * **S9e.4b.3c.2 refined, before its code (2026-10-05).** Why the
    survey's 13 exact incidences on one sphere are refused today:
    `graph.rs`'s pass over the circles of both inputs on one sphere
    (S9e.4b.3c.1) refuses two arcs of one circle (their planes one,
    `Circ::meet_plane` none) and a crossing at an end of either where both
    are the inputs' own (`OutOfDomain("a vertex or a circle of both inputs
    on one sphere (S9e.4b.3c.2)")`), found before the pierces; past it the
    arrangement keeps each input's vertices and edges apart, so a vertex of
    one on the other's edge or vertex would be `Degenerate("a vertex of one
    input on the other's face")`, edges of both on one curve
    `Degenerate("edges of both inputs overlapping")`, an axis edge of one
    along the other's half-plane `Degenerate("an edge of one input on a
    face of the other")`, and a section along an edge of both faces a
    seam's conflict at every seam. The 13: `so1` and `so2` (5:
    `bcommon_complex/B3`, `bcut_complex/C7`, `C8`, `bfuse_complex/B6`,
    `bugs/moddata_1/bug183_4`: the wedge's equator arc on the hemisphere's
    rim circle, its base on the hemisphere's base, its corner inside that
    base, its pole the hemisphere's stored pole), `so2` and `so3` (4:
    `bcommon_complex/B4`, `bcut_complex/C9`, `D1`, `bfuse_complex/B7`: the
    corner, the pole and the axis edge of both, the bases on one plane with
    overlapping edges, each equator arc's end inside the other's arc, the
    meridian half-planes crossing along the axis) and `so5` and `so2` (4:
    `bcommon_complex/B6`, `bcut_complex/D4`, `D5`, `bfuse_complex/B9`: the
    higher wedge's axis edge inside the other's, its corner on it, the pole
    of both). The survey's files are on the world's axes: their planes'
    stored normals keep exact zeros (`gp_Dir`'s normalization and
    `Frame3::new`'s keep a zero component zero), so a meridian half-plane
    holds the world's `z` exactly and two inputs' half-planes cross exactly
    along it. Decisions. (1) *Scope.* An arrangement whose inputs have faces
    on one sphere (S9e.4b.3c.1's faces on one surface: equal exact centres
    and radii), exact incidences only: points equal in the exact models
    (rationals and quadratic surds of the stored binary64 data), never
    within a tolerance. Off one sphere S9's rules are unchanged (prisms
    touching at a vertex or along an edge stay `Degenerate`). (2) *A vertex
    of both.* Model vertices of the two inputs at one exact point are one
    arrangement vertex (`VKey::Both`), on every face of both at it; it
    continues both vertices (the tool's generated in a cut, S9a's rule).
    (3) *A vertex of one inside an edge of the other.* A model vertex of
    either strictly inside a line or circle edge of the other (exactly on
    its curve, its place strictly between the ends) splits that edge there
    and lies on its two faces; it keeps its own name. (4) *An edge of
    both.* Each input's edges split at every vertex of both on them, a part
    of B's edge with the same two ends as a part of A's, on one line or one
    circle and over the same arc (its middle strictly inside A's part), is
    A's arrangement edge, its half-edges in B's two faces with B's run
    (`Arr::shared`); it continues both edges. A face on one surface with
    the other's face has it as its own boundary already (not as the other's
    interior edge); a section of two faces along it (the half-planes of two
    wedges crossing along their axis, neither on the other's surface) is
    taken by it; and an edge of one along the other's face, on one line or
    circle with an edge of that face, is an edge of both where they
    overlap and outside the face elsewhere (else `Degenerate("an edge of one
    input on a face of the other")`). (5) *Plane faces on one plane* are
    S9c.1's faces on one surface, each holding the other's edges strictly
    inside it, their overlapping edges edges of both by (4), the kept pieces
    of both facing one way joined as before. (6) *Circles of both, pierces
    at known vertices.* Two arcs of one circle meet only at the vertices (2)
    and (3) give them; a crossing at an end of either circle's edge is such a
    vertex; an edge's pierce of a face at a vertex of that face already on
    the edge is taken. (7) *Splits alike.* A whole sphere's or a piece's
    split great circles (the second arrangement's seams, edges and vertices
    without ids) take part in (2) to (4) as the inputs' own edges do: a
    split's vertex on the other's edge or vertex is no longer a seam's
    conflict (in `so1` and `so2` a split great circle's end on the
    hemisphere's rim lies on the wedge's equator arc at every seam tried); a
    split great circle on another input's circle (one plane) stays one,
    retried. (8) *Degenerate and refused.* An incidence within the
    resolution and not exact is S9's: two faces within the resolution of
    one plane, a piece or result thinner than it, two meetings within
    rounding (the survey's `so6` and `so7` pairs, OCCT's turned copies whose
    corners are within `1e-15` of each other; any exact incidence moved by
    a turn that rounds its frames); curves tangent at a vertex of both
    (`next_on`); a result touching itself along an edge of both or at a
    vertex; spheres within the resolution of one. S9's near-plane guard
    (planes within the resolution of each other over their faces' boxes'
    overlap) stays as it is, so two inputs' faces meeting only along a line
    or at a point (two octants about one axis) are `Degenerate` with it.
    S9e.4b.3c.3's refusals (a piece not convex in its planes, two plane
    faces on one plane of one piece) unchanged. (9) *Fuzzing.* Where the
    new cases become reachable (a split piece against its own primitive on
    one frame: the stored pole and rims of both), the `split` target's
    `PIECE_BOOLEANS` takes them. (10) *Evidence first.* Bodies OCCT writes
    on the world's axes (in turned rational frames an exact incidence of two
    inputs' planes along a line holds only where both planes are stored bit
    for bit alike), each the ball common `MakeBox` on a frame whose normal
    is the axis and whose origin lies on it, reproducing `so1`, `so2`,
    `so3` and `so5` and further incidences (a half and an octant, wedges
    about another world axis, caps on one frame, a chain), declared
    `degenerate` the near ones; their reference of its own, exact (the
    chained reference's sweeps follow the world's meridians and parallels,
    which these faces lie on): the ball cut into cells by every input's
    heights along the common axis and half-planes about it, closed forms per
    cell (`generate_one_sphere_incidence_boolean_fixtures.py --check`, a CI
    group `one-sphere-incidence`,
    `test_one_sphere_incidence_boolean_reference.py`); a native capture
    keyed on the refusal it removes
    (`compare_one_sphere_incidence_boolean.py`); then the kernel and its
    tests (`tests/one_sphere_incidence_booleans.rs`). S9e.4b.3c.1's
    `hemi_octant` (the hemisphere and the octant on one turned frame)
    becomes a solid case of its set.
  * **S9e.4b.3c.2 evidence (2026-10-05).** Eleven bodies OCCT wrote
    (`boolean-one-sphere-incidence-bodies.txt`, `write` blocks of a
    sphere's row, `boolean common` and a `box` row, `BRepPrimAPI_MakeBox`
    on a frame whose normal is the axis and whose origin lies on it, written
    by `compare_one_sphere_incidence_boolean.py --write-bodies` to
    `rust/fixtures/imported/`), every one a piece of the ball of radius 5
    about `(5, 5, 4)` on the world's axes: about `z`, `incidence_hemi`, above
    the centre's parallel, its rim divided as `so1`'s; `incidence_wedge45`,
    between the half-planes at 45 and 135 degrees (`so2`'s);
    `incidence_wedge23`, from `atan2(5, 12)` a quarter turn on (`so3`'s);
    `incidence_high23`, the same above the parallel at `5/2` (`so5`'s);
    `incidence_half`, on one side of the meridian plane along `x`;
    `incidence_octant` and `incidence_back`, the quarters at 0 and 180
    degrees; the near ones `incidence_turned45` (the wedge's frame's `x`
    turned by about `2^-36`) and `incidence_lifted45` (its box `2^-40`
    above the centre); about `x`, `incidence_x45` and `incidence_x23`.
    `generate_one_sphere_incidence_boolean_fixtures.py --check`: 45 cases
    of 19 groups (29 solid, 4 empty, 12 declared `degenerate`; 34 of class
    `pieces`, 8 `sphere`, 3 `chain`): `hemi_wedge`, `wedge_hemi` (`so1` and
    `so2`), `wedge_wedge`, `wedge_back` (`so2` and `so3`), `high_wedge`,
    `wedge_high` (`so5` and `so2`), `half_octant`, `octant_half` (a vertex
    of both on the rim, the octant's corner inside the half's diameter, its
    pole inside the half's semicircle), `x_wedges`, `same_wedge` (the wedge
    and itself), `cap_cap`, `cap_back` (two of S9d.1's caps on one frame:
    the stored pole of both), `cap_wedge`, `wedge_cap`, `chain_zone` (the
    wedges' fuse, then with a zone of the sphere: its cut two solids);
    declared `degenerate`: `quadrants` (two octants about one axis touching
    along it only), `half_wedge` (the wedge's axis edge inside the half's
    meridian face, off its edges), `near_turned` and `near_lifted` (two
    faces within the resolution of one plane). The reference is exact and
    of its own: each input read off its rows as the ball cut to heights
    along the common axis and to the half-turns of its half-planes through
    the axis (each box face asserted, in the rows' rationals, to hold the
    axis, to be normal to it or to miss the ball), the ball cut into cells
    by every input's heights and half-planes, each cell kept by the set
    function at its middle, its volume, moments and sphere area in closed
    form and the faces between kept and dropped cells added (parallels'
    sectors, half-planes' pieces); checks relative to the size: the volume
    again by the divergence theorem over those faces within 4.7e-41, the
    pair identities 4.7e-41, every Monte-Carlo point's membership by the
    rows' own inequalities that of its cell (100,000 points a group, none
    differing), Monte Carlo 2.5 standard errors, the solid counts by the
    kept cells' adjacency (kept cells about the axis apart: a result
    touching itself, declared `degenerate`); each file read independently:
    its faces one sphere and planes, every stored vertex within 3.5e-17 of
    the size on the construction's surfaces, every stored plane (its normal
    the cross product of its stored axes, in rationals) holding the axis,
    normal to it or apart from the ball exactly, the divided rim two arcs.
    `test_one_sphere_incidence_boolean_reference.py` checks the cells'
    closed forms against quadrature, a wedge's measures against their
    textbook values, the hemisphere and the wedge's common and fuse, the
    octants' touching, the bodies read as axial pieces, the case list and
    its protocol rows. The generator's check is a CI group of its own
    (`one-sphere-incidence`, under a minute on four workers locally);
    Python 3.9 and 3.12 write the same files. Corrections from the evidence,
    amending the refined decisions: (a) bodies in turned rational frames
    hold an exact incidence along a line only between planes stored bit for
    bit alike, so the bodies are on the world's axes, as the survey's files
    are, and the chained reference (whose sphere families are the world's
    meridians and parallels, along these faces) gives way to the cells; (b)
    importing `incidence_hemi` panicked: S9e.4b.3c.1's `Arr::split_at`
    scaled a stored vertex's place onto a section circle without testing
    it, and the stored pole above the hemisphere's rim, exactly on the
    rim's axis on the world's axes, has none (in turned frames it rounds
    off the axis); fixed before the capture, its own commit; (c) two
    octants about one axis are `Degenerate` all three ways by S9's
    near-plane guard (their faces meet only along the axis, their boxes'
    overlap thinner than the resolution); (d) caps of one frame on one
    sphere (their stored pole of both, `cap_cap`) are refused today as this
    step's incidence, as the refined decisions' S9e.4b.3c parent said.
    The capture `occt-boolean-one-sphere-incidence-preimplementation`
    (`compare_one_sphere_incidence_boolean.py`, keyed on the refusal the
    step removes, `graph.rs`'s `OutOfDomain("a vertex or a circle of both
    inputs on one sphere (S9e.4b.3c.2)")`; the kernel's probe `unsupported`
    on the 36 solid, empty and `half_wedge` cases, every one by that
    refusal, and `refused` on the 9 near and touching ones by S9's
    near-plane guard, raised first, which `compare_boolean`'s
    `Set.refused_before_code` admits for a declared `degenerate` case): 42
    matching, 3 reviewed. OCCT leaves the hemisphere less the wedge (`so1`
    and `so2`'s `bcut_complex/C7`) invalid under BRepCheck_Analyzer, its
    measures within 3.8e-15 of the reference; of the declared near cases
    it drops the slivers thinner than its tolerance (the wedge less its
    turned copy empty where the reference's sliver is 3.0e-10 in volume;
    the hemisphere less the lifted wedge without its sliver under the
    wedge, invalid). S9e.4b.3c.2's kernel next.
  * **S9e.4b.3c.2 implemented** (`curved/graph.rs`: an arrangement of
    inputs with faces on one sphere merges their exact incidences;
    `assemble.rs`: the merged entities' names; `matched.rs`,
    `assemble.rs`: a given model's face pieces named by edges of both), as
    the refined decisions describe: model vertices of both inputs at one
    exact point one vertex (`VKey::Both`, on every face of both, continuing
    both vertices), a model vertex of either strictly inside a line or circle
    edge of the other splitting it (on that edge's faces), and, each input's
    edges split at every vertex of both on them, a part of B's edge with the
    ends, the line or circle and the arc of a part of A's that arrangement
    edge (`Arr::shared`: its half-edges in B's faces with B's run,
    continuing both edges; `Arr::model_edge` gives either input's edge of
    it); two arcs of one circle meet only at those vertices, a crossing at
    an end of either circle's edge is one of them, a pierce at a vertex of
    the face already on the edge is taken, a section of two faces along an
    edge of both is taken by it, an edge of one along the other's face on
    one line or circle with an edge of that face is an edge of both where
    they overlap and outside the face elsewhere (else `Degenerate("an edge
    of one input on a face of the other")`), and a face on one surface with
    the other's has an edge of both as its own boundary, not as the
    other's interior edge. The splits' own great circles take part alike.
    The refusal `OutOfDomain("a vertex or a circle of both inputs on one
    sphere (S9e.4b.3c.2)")` is gone. All 45 fixtures as declared (33 solid
    and empty within the kernel's enclosures, each at most `1e-9` wide; the
    12 degenerate refused with S9's reasons: two faces within the resolution
    of one plane for the near copies and the two octants about one axis, an
    edge of one input on a face of the other for the half and the wedge),
    every history complete over the imported bodies' stored ids, an edge or
    vertex of both continuing both inputs' (the two wedges' axis edges,
    corners and poles one result edge and two vertices; the higher wedge's
    axis edge inside the other's), results deterministic and moved rigidly,
    both inputs moved by exact motions (a dyadic translation, a quarter turn
    about a world axis: the stored zeros kept) keeping the reference's
    volumes and turned by a rotation that rounds their frames refused as
    `Degenerate` or within the reference (`tests/one_sphere_incidence_booleans.rs`,
    4 s in release and 6 s at `opt-level` 2 with debug assertions on a host
    at load 18, 8 s in release with debug assertions under the emulated
    correctly rounded `hypot`, where S9e.4b.3c.1's and the pieces' test
    files pass too). S9e.4b.3c.1's `hemi_octant` (the hemisphere and the
    octant on one turned frame), declared `unsupported` until now, is solid
    within its reference, its margins those of any solid case (its
    generator, unit test and `tests/one_sphere_booleans.rs` declare it so).
    `compare_one_sphere_incidence_boolean.py` 42 matches and 3 reviewed,
    every enclosure within the reference; `compare_one_sphere_boolean.py`
    8 and 14 unchanged (now all 20 solid and empty cases within the
    reference); every other comparison of `HANDOFF.md`'s table unchanged
    with 0 failures; the release suite (613 tests) and the tools' unit
    tests (333) passing, the ledger unchanged. Amendments, from the
    implementation: (a) a split's vertex on the other's edge or vertex, a
    seam's conflict in S9e.4b.3c.1's pass over the circles, is merged as the
    inputs' own are (decision (7)): in `so1` and `so2` the hemisphere's
    split great circle's end on its rim lay on the wedge's equator arc at
    every seam tried, so a first version of this step keeping those
    conflicts refused every case as a meeting at every seam tried;
    (b) a given model's face pieces named by its edges (a face of several
    stored faces, a result of several solids) take an edge of both as B's
    too (`Arr::model_edge`); (c) the `split` target's `PIECE_BOOLEANS`
    gives a zone's or cap's first piece also to the zone or cap it was split
    from, built again under other ids on its frame (a Boolean's inputs
    never share ids): their common the piece, their fuse the solid, the
    solid less the piece the rest, within `1e-9` of the solid: of the split
    corpus and its regressions 265 inputs reach it, all three evaluate for
    134, 127 are refused as two faces within the resolution of one plane and
    4 as pieces or results thinner than the resolution (S9's rules).
    Replays with debug assertions, one process an input, natively and under
    the emulated correctly rounded `hypot`: the split corpus and its
    regressions (3,569 inputs) and the boolean corpus and its regressions
    (1,466), no failure, the slowest 2.7 s and 7.3 s (2.3 s and 6.6 s
    emulated). A trial of the DRAW survey's restore cases on the Rust
    adapter and natively (not the survey: nothing registered): the 13
    targeted (`bcommon_complex/B3`, `B4`, `B6`, `bcut_complex/C7`, `C8`,
    `C9`, `D1`, `D4`, `D5`, `bfuse_complex/B6`, `B7`, `B9`,
    `bugs/moddata_1/bug183_4`) evaluate on both backends with every check;
    of the 38 `so` cases 27 evaluate (S9e.4b.3c.1's 14 and these 13), the
    11 with `so6` or `so7` `Degenerate` by S9's rules as before;
    `bcut_complex/I6` and `G4` refused as S9e.4b.3c's as before. Pending:
    the capture's Linux record. DRAW survey: that of S9e.4b.3c.2 and
    S9e.4b.3c.3a, after S9e.4b.3c.3a's implemented bullet (the 13
    registered). Campaigns at `8a55a3e6`
    (600 s each, sampled replays, every switch on): `boolean` clean, 969
    runs, the slowest input 19 s under AddressSanitizer at load about 8;
    `split` clean, 896 runs, the slowest 10 s.
  * **S9e.4b.3c.3 refined, before its code (2026-10-05).** Why the survey's
    last S9e.4b.3c cases are refused today, read off the dataset's files by
    the converter (the dataset is not committed). `bcut_complex/G4`'s part
    (`CTO900_fra50089-part`) is a box `[-100, 100] x [-50, 50] x [-100, 100]`
    fused with a boss, a cylinder of radius 52.93 along the world's `y` from
    the box's face `y = -50` (which holds the boss's circle as a hole) to a
    cap at `y = -150` (S9e.4b refined listed it among prisms with walls of two
    directions): its primitive common its planes is the cylinder through the
    whole box, which S9e.2's match refuses (`pieces::model`'s `OutOfDomain("an
    imported plane piece other than its primitive common its planes'
    half-spaces (S9e.4b.3c)")`); its tool is a cylinder of radius 37.90
    parallel to the boss from `y = -150` to `y = -50`, its caps on the boss's
    cap plane and on the box's face exactly, its wall crossing the boss's
    along two lines. `bcut_complex/I6`'s tool (`CTO902_cts20455-tool`) is a
    block `[45, 127.76] x [0, 70] x [0, 20]` less a cylinder of radius 35
    along `z` about `(45, 35)`: its wall's material outside its quadric,
    refused on import before its model (S9e.4b.3c.1's amendment (c)), and
    tangent to the block's own faces `y = 0` and `y = 70` along the two edges
    it shares with them; the case's part, an L prism, has its face `x = 80`
    tangent to the tool's wall along the generatrix at `y = 35`. `shading_132`
    (`bugs/modalg_1/buc60926`) is a frustum of radii 1 and 2 over height 2
    less the quarter between two planes through its axis, refused in its own
    first arrangement as S9d.3a's `Degenerate("a plane through a cone's
    apex")` (the planes pass through the frustum's virtual apex); its case
    fuses it with a unit ball through two of its rim's vertices. S9e.4b.3a's
    `bitten` (a ball less a box's corner) is refused as G4's part is, and the
    `boolean` fuzz target's `IMPORTED` stage refused 4 first results as
    S9e.4b.3c's (2 not their primitive common their planes, 2 with two faces
    on one plane: `hull_of`'s `OutOfDomain("an imported plane piece with two
    faces on one plane (S9e.4b.3c)")`). Decisions. (1) *Sub-steps*, by what
    unlocks the most and by the machinery each takes: **S9e.4b.3c.3a** (this
    step): a body of one sphere, cylinder or cone face and plane faces that is
    one Boolean of its primitive and the convex hull of its other planes:
    S9e.4b.3a's common, or the hull less the primitive (a hole, groove, slot,
    notch or dimple: the curved face's material outside its quadric), the
    primitive less the hull of its planes turned over (a bite: `bitten`,
    `shading_132`'s three-quarter frustum), the hull fused with the primitive
    (a boss: G4's part); two plane faces on one plane facing one way (a groove
    across a face) one plane of the hull; OCCT's vertex loops at a sphere's
    poles inside its face matched. Of the survey's cases only G4 can evaluate
    with it (I6 and `buc60926` stay `Degenerate`, (5)); the fuzz target's
    grooves, bosses and bites reach it. **S9e.4b.3c.3b**: a body of one curved
    face that needs more than one hull (a union of primitive commons and
    hulls: a primitive bitten twice, a boss, groove or bite on or by a body
    not convex in its planes), a primitive's end not normal to its axis or a
    sphere's ends (S9e.4b.3a's discs) with another form. (2) *The model.* The
    body is the given model (S9e.1, S9e.3a) of that first Boolean, as
    S9e.4b.3a's: its inputs the primitive's model and the hull leaf model
    (`pieces.rs`'s `Hull`), the hull the first input of the hull less the
    primitive (the given model's frame the primitive leaf's). The primitive: a
    whole sphere on the stored frame and radius; a cylinder or a cone on the
    curved face's stored frame over the curved face's own axial range (its
    edges' heights along the axis), past each end by a quarter of it (at least
    of the radius) but at a cap, a plane face normal to the axis at an end of
    that range whose outward normal points away from the range (the
    primitive's own end, its material inside the primitive: a boss's cap, a
    bitten cylinder's discs) or, for the hull less the primitive, into it (the
    primitive's end seen from outside: a slot's or a blind hole's floor), the
    end the cap's exact height along the axis rounded once; a cone clamped at
    its apex; S9e.4b.3a's common keeps its primitive past the body's ends and
    every plane. The hull: the planes of the plane faces but the caps, turned
    over for a bite (the hull its bite), bounded by S9e.4b.3a's cube; two
    faces on one plane facing one way are one plane of the hull, its face
    holding both stored faces (S9e.2's match takes a model face for several
    stored faces). Numbers: the engine's own, no new field or degree. (3)
    *Recognition.* After S9e.4a's construction fails, a body of S9e.4b.3a's
    shape (one solid region, one shell, one sphere, cylinder or cone face and
    plane faces) is tried as the hull less the primitive where its curved
    face's material lies outside its quadric, else as the common, the bite and
    the boss in turn: the first whose model matches the stored topology is the
    body (S9e.2's match the arbiter: every stored vertex within the resolution
    of the model's, edges through their points, faces by their edges), none
    matching `OutOfDomain("an imported plane piece other than one Boolean of
    its primitive and its planes' hull (S9e.4b.3c.3b)")`, a form's other
    refusal (S9's, raised in its own first arrangement) reported, the first in
    that order. A stored vertex loop at a pole of its sphere face's stored
    frame that the re-run lacks is left unmatched: OCCT writes one at every
    pole inside a sphere face, where the kernel's assembly closes a loop
    through one pole without a vertex at the other (S9e.4b.3a's amendment (b);
    `bitten`'s face holds its frame's south pole); a Boolean's history deletes
    it unless a result's vertex continues it. (4) *Its other queries.*
    Classification: the primitive's and the hull's locations (within the
    resolution) combined by the form (a common the lesser, a fuse the greater,
    a complement reversed); bounds S9e.4b.3a's (the primitive's between the
    body's axial ends widened by the stored edges': a body less its primitive
    lies in its hull, whose extremes are on the plane faces' edges); a rigid
    motion reads the piece off the moved stored topology again in its form; a
    result of it given to another Boolean is S9e.3a's. (5) *Degenerate and
    refused.* S9's rules unchanged in both arrangements. A tangency in a
    piece's own first arrangement is its curved face's with its own plane
    faces, `Degenerate("an imported plane piece whose curved face is tangent
    to its plane faces")`, named for the body: I6's tool, whose wall is
    tangent to two of its block's faces along the edges it shares with them,
    is a body S9 cannot arrange (a tangency is no meeting it decides), its
    reason now its own (S9e.4b.3a's survey found it as "a tangency between the
    inputs", raised in its own arrangement; S9e.4b.3c.1 refused it before its
    model as a groove). I6 stays refused either way, its part's face `x = 80`
    tangent to the tool's wall (S9's tangency between the inputs).
    `shading_132`'s three-quarter frustum is its frustum's bite by the quarter
    between two planes through its axis, `Degenerate("a plane through a cone's
    apex")` in its own arrangement as now (S9d.3a); `buc60926` stays refused.
    Two faces on one plane facing apart bound nothing: `Degenerate("an
    imported plane piece's two faces on one plane facing apart")`. Bodies no
    single form matches are S9e.4b.3c.3b's. (6) *Fuzzing.* The `boolean`
    target's `IMPORTED` stage already imports the chained stage's first
    results of one curved face and planes (S9e.4b.3a's widening) and cuts them
    by the turned box again: a box less or fused with a sphere, cylinder or
    cone, and those less a box, now reach the forms, the volume the chained
    cut's; the replay measures how many. (7) *Evidence first.* Bodies OCCT
    writes (a `write` block: one solid's rows, a `boolean` row, the other's
    rows), reproducing the survey's shapes without the dataset, every section
    a circle or a line (the reader takes no ellipse; a sphere's section off
    its frame's parallels and meridians has an uncertified pcurve):
    `form_boss` (G4's part on the world's axes: a box fused with a cylinder
    along `y` from its face to a cap), `form_scoop` (a box on `SKEW` less a
    ball centred on its top face: the top face in two faces on one plane),
    `form_slot` (a box less a cylinder along `x` ending inside it at its
    floor; a groove through a box is S9e.4a's prism), `form_dimple` and
    `form_ball_boss` (a box on a turned frame less or fused with a ball whose
    section is its parallel), `form_sink` (a box on `SKEW4` less a frustum
    through it: a conical hole), `form_bite` (a cylinder less a box between
    its caps on the world's axes: in a turned frame a plane along a cylinder's
    axis is within rounding of it, S9's refusal), declared `degenerate`
    `form_notch` (I6's tool, a box less a cylinder tangent to two of its
    faces) and `form_quarter` (`shading_132`'s three-quarter frustum on the
    world's axes); cases of each against rods, slabs, balls and boxes, the
    piece as the tool, two imported, a chain, and G4's configuration (the boss
    less a parallel cylinder whose caps lie on the boss's cap plane and the
    box's face); the reference S9e.3a's chained reference on the constructions
    OCCT was given (`(X op Y) op2 C`), with S9e.4b.3a's checks and each body's
    closed form where it has one (`generate_piece_forms_boolean_fixtures.py
    --check`, a CI group `piece-forms`,
    `test_piece_forms_boolean_reference.py`); a native capture keyed on the
    refusal it removes (`compare_piece_forms_boolean.py`, `imported.rs`'s
    refusal of a piece other than its primitive common its planes); then the
    kernel and its tests (`tests/piece_form_booleans.rs`). S9e.4b.3a's
    `bitten_box`, declared `unsupported` there, becomes a solid case of its
    set.
  * **S9e.4b.3c.3a evidence (2026-10-05).** Nine bodies OCCT wrote
    (`boolean-piece-forms-bodies.txt`, `write` blocks of one solid's rows, a
    `boolean` row and the other's rows, written by
    `compare_piece_forms_boolean.py --write-bodies` to
    `rust/fixtures/imported/form_*.brep`), every section a circle or a line:
    `form_boss` (G4's part on the world's axes: the box `[0, 10] x [0, 8] x
    [0, 10]` fused with a cylinder of radius 5/2 along `y` from its face `y =
    0` to a cap at `y = -4`), `form_scoop` (a box `12 x 6 x 6` on `SKEW4` less
    a ball of radius 9/2 centred on its top face, the ball's frame normal the
    box's `y` and its `x` the box's normal: the top face a meridian plane in
    two faces on one plane, the side faces parallels' planes, the seam above
    the body), `form_slot` (a box on the world's axes less a cylinder of
    radius 5/4 along `x`, its axis 1/2 above the top face, ending inside the
    box at its floor), `form_dimple` and `form_ball_boss` (boxes on `SKEW` and
    `SKEW2` less a ball of radius 5/2 3/2 above the top face's middle, fused
    with one 1 below it), `form_sink` (a box on `SKEW4` less a frustum of
    radii 1 and 3 along its normal through it), `form_bite` (a cylinder of
    radius 3 on the world's axes less a box across its wall between its caps),
    declared `degenerate` `form_notch` (I6's tool: a box less a cylinder of
    radius 3 along `z` about the middle of one face, tangent to the two faces
    it meets) and `form_quarter` (`shading_132`'s: the frustum of radii 2 and
    1 over 2 on the world's axes less the quarter between two planes through
    its axis). `generate_piece_forms_boolean_fixtures.py --check`: 39 cases of
    13 groups (33 solid, 6 declared `degenerate`; 33 of class `pieces`, 3
    `both`, 3 `chain`): `boss_rod` (G4's configuration: the boss and a
    cylinder of radius 2 along `y` whose caps lie on the boss's cap plane and
    on the box's face, crossing the boss's wall along two lines), `boss_slab`,
    `scoop_rod`, `slab_scoop` (the scoop the tool), `slot_rod`, `dimple_ball`,
    `ball_boss_slab`, `sink_rod`, `bite_box`, `pieces` (the boss and the slot,
    both imported), `chain_slot` (the slot less a rod, then with a `TILT`
    slab); declared `degenerate` `notch_box` and `quarter_ball` (the
    three-quarter frustum and a ball, as `buc60926`). The reference is the
    constructions OCCT was given through S9e.3a's chained reference (`(X op Y)
    op2 C`, swapped, two bodies, chained), with S9e.4b.3a's checks relative to
    the case's size: the two families within 3.4e-36, each solid's closed form
    1.4e-40, each body's closed form 3.7e-38 (G4's part, the slot's segment,
    the dimple's and the ball boss's caps on the boxes' stored axes, their
    determinant a unit within rounding), the pair identities 1.4e-41 and the
    area identity 7.0e-41 (where no two inputs share a surface), Monte Carlo
    2.6 standard errors (100,000 points a group), quadrature estimates
    9.8e-33, every meeting's sine at least 0.11 and events at least 1.7e-6 of
    their range apart outside the declared groups, solid counts by rays at two
    resolutions with S9e.4b.3c.1's join within two grid spacings (the scoop's
    side faces on `SKEW4` lie nearly along the rays, whose grazing intervals
    the strict overlap counted as further solids); the declared groups'
    checks, whose reference meets the notch's tangency and the frustum's
    planes through its apex, kept apart within 1.3e-21. Every file read
    independently (`stored_records`): its faces one sphere, cylinder or cone
    and planes, every stored vertex within 4.6e-16 of the size on the
    construction's surfaces, every body a Boolean other than a primitive
    common a box. `test_piece_forms_boolean_reference.py` checks the segments'
    and caps' closed forms against quadrature, G4's part's measures, the
    chained reference on the slot against its closed form, the frames keeping
    every section a circle or a line (the scoop's ball's axes in the rows'
    integers, the seams off the bodies), the bodies being this step's, the
    case list and its protocol rows and the files apart from every earlier
    set's. The generator's check is a CI group of its own (`piece-forms`, 8
    minutes on four workers locally); Python 3.9 and 3.12 write the same
    files. Corrections from the evidence, amending the refined decisions' plan
    (7): (a) a box less a cylinder through it along one of its axes is
    S9e.4a's prism (a profile of lines and an arc), so the groove across a
    face is the scoop's ball (the top face in two faces) and the cylinder's
    groove the slot, ending inside the box at its floor (a cap of the hull
    less the primitive); (b) a plane along a cylinder's axis in a turned frame
    is within rounding of it (S9's `Degenerate("a plane within rounding of a
    cylinder's direction")`), so the bitten cylinder is on the world's axes,
    as is the three-quarter frustum, whose planes then pass through its apex
    exactly (`shading_132`'s refusal; on `SKEW` they met the cone within
    rounding of each other instead); (c) OCCT splits a curved face its seam
    crosses, so the notch's cylinder and the frustum have their frames' `x`
    turned off the bodies; (d) the files `notch.brep` and `slot.brep` were
    S9e.4b.1's: this step's are `form_*`. The capture
    `occt-boolean-piece-forms-preimplementation`
    (`compare_piece_forms_boolean.py`, keyed on the refusal the step removes,
    `imported.rs`'s `OutOfDomain("an imported plane piece other than its
    primitive common its planes' half-spaces (S9e.4b.3c)")`; the kernel's
    probe `unsupported` on the 36 solid and notch cases, the bodies refused on
    import by that refusal, and `refused` on the three-quarter frustum's 3 by
    S9d.3a's apex, raised in its own arrangement first, which
    `compare_boolean`'s `Set.refused_before_code` admits for a declared
    `degenerate` case): every result valid, 26 matching (within 2e-8 of the
    reference), 13 reviewed: the scoop's and the slot's rods, the dimple's
    ball, the conical hole's rod and the three-quarter frustum's ball, volumes
    up to 1.8e-5 relative and centres up to 2.2e-6 of the case's size off by
    BRepGProp's default integration (a sphere's or a cone's faces met by
    planes and a cylinder), within 1.4e-9 measured adaptively at 1e-10 and
    1e-12 by a diagnostic build. G4's configuration matches. S9e.4b.3c.3a's
    kernel next.
  * **S9e.4b.3c.3a implemented** (`solid/imported.rs`'s `Form`: the forms
    tried in turn, a cylinder's or cone's primitive over its curved face's own
    axial range ending at its caps, the classification by the form;
    `curved/pieces.rs`: the form's first Boolean, the hull the first input of
    the hull less the primitive, two faces on one plane facing one way one
    plane of the hull, a tangency in the piece's own arrangement named for the
    body; `curved/given.rs` and `curved/matched.rs`: the given model's frame
    its primitive leaf's, a stored pole's vertex loop the re-run lacks left
    unmatched for an imported piece), as the refined decisions describe: a
    body of one sphere, cylinder or cone face and plane faces is tried as the
    hull less the primitive where its curved face's material lies outside its
    quadric, else as S9e.4b.3a's common, the primitive less the hull of its
    planes turned over and their fuse, the first whose given model matches its
    stored topology the body, none matching `OutOfDomain("an imported plane
    piece other than one Boolean of its primitive and its planes' hull
    (S9e.4b.3c.3b)")`. All 39 fixtures as declared (33 within the kernel's
    enclosures, each at most `1e-9` wide; the notch's 3 refused as
    `Degenerate("an imported plane piece whose curved face is tangent to its
    plane faces")`, the three-quarter frustum's 3 as `Degenerate("a plane
    through a cone's apex")`), every history complete over the imported
    bodies' stored ids, the scoop's two faces on its top plane each continuing
    apart, results deterministic and moved rigidly, both inputs translated and
    turned keeping the reference's volumes, every body in its form (its volume
    its closed form, points in and off it classified, its stored vertices on
    its boundary, `bitten` a bite through OCCT's pole vertex loop), and the
    kernel's own dimple, ball boss, rod boss and bitten rod written by its
    writer, read back and imported, their Booleans with a turned box the
    kernel's own results' (`tests/piece_form_booleans.rs`, 7.8 s in release on
    a host at load 8 to 14, 5.8 s at `opt-level` 2 with debug assertions, 6.2
    s in release with debug assertions under the emulated correctly rounded
    `hypot`, where the pieces' earlier test files pass too). S9e.4b.3a's
    `bitten_box`, declared `unsupported` until now, is solid within its
    reference (its generator declares it so, its margins unchecked: its box's
    edge runs through the ball's stored pole within rounding, the reference's
    events there 6.2e-18 of their range apart; its unit test and
    `tests/imported_piece_booleans.rs` declare it so), and S9e.4b.3c.1's notch
    (a block less a leaning rod) imports as its hull less its rod.
    `compare_piece_forms_boolean.py` 20 matches and 19 reviewed (the 13
    captured, with the kernel the boss's slab, the bitten cylinder's box and
    the dimple's fuse in entity counts: OCCT's unified results split their
    cylinders' and spheres' faces and sections at its seams and keep its pole
    edges, where the kernel's faces close over their periods), every enclosure
    within the reference; `compare_imported_pieces_boolean.py` 30 and 15 (the
    bitten ball's fuse and cut in entity counts likewise); every other
    comparison of `HANDOFF.md`'s table unchanged with 0 failures (and
    `compare_split.py` 72/56, `compare_brep.py --family spline` 10/3,
    `compare_brep_io.py` 6,835/7, `compare_step.py` 23/6 on STEP-b's SDK); the
    release suite (621 tests) and the tools' unit tests (338) passing, the
    ledger unchanged. Amendments, from the implementation: (a) a given piece's
    edge along its cylinder's section by a plane along its axis (the slot's
    rims, their bases quadratic surds) met another cylinder as
    `ComputationLimit("an irrational line against a cylinder")` (S9c's
    line-cylinder meeting took a generatrix's base rational, true of edges
    until now); where the quadratic's discriminant, in the base's field, is
    negative the line is apart from the cylinder and meets it nowhere
    (`meet.rs`; the slot's rod), elsewhere the limit stands; (b) the given
    model's frame, which places its faces' cylinders' angles, is its primitive
    leaf's where the hull is its first input; (c) the forms' refusals other
    than a mismatch are reported in the order tried, the first's (the
    three-quarter frustum's apex from the common, a notch's tangency from the
    hull less the primitive); (d) the kernel's own rods given to the roundtrip
    are on the world's axes: in a turned frame a bite's plane along a rod's
    axis, written and read back, its frame normalized again, is within
    rounding of the axis under one platform's `hypot` or the other's (S9's
    refusal; reproduced under the emulation). The `boolean` fuzz target's
    `IMPORTED` stage reaches the forms through the chained stage's first
    results of one curved face and planes (no new stage): replayed with debug
    assertions, 92 such first results of the corpus's 1,466 inputs and
    regressions reach it, 19 imported and cut within the chained cut's volume
    (17 before), 63 not written by the kernel's writer, 6 read back with
    pcurves off their edges, 4 refused as S9e.4b.3c.3b's (a ball's groove in a
    U prism, not convex in its planes; a box fused with a ball's half, a
    sphere's end with another form). Replays with debug assertions, one
    process an input, natively and under the emulated correctly rounded
    `hypot`: the boolean corpus and its regressions (1,466 inputs) and the
    split corpus and its regressions (3,569), no failure, the slowest 7.0 s
    and 2.1 s (5.2 s and 2.1 s emulated). A trial of the DRAW survey's restore
    cases on the Rust adapter and natively (not the survey: nothing
    registered): `bcut_complex/G4` evaluates on both backends with every check
    (its area 193,593.0 on the Rust adapter, native DRAW's printed 193593);
    `bcut_complex/I6` is refused for its tool alone, `Degenerate("an imported
    plane piece whose curved face is tangent to its plane faces")` (its part
    is tangent to the tool too); `bugs/modalg_1/buc60926` stays `Degenerate("a
    plane through a cone's apex")`; the 38 `so` cases as before (27
    evaluating, 11 `Degenerate`). Pending: the capture's Linux record. DRAW
    survey: that of S9e.4b.3c.2 and S9e.4b.3c.3a, below (`G4` registered,
    `I6` and `buc60926` refused as the trial found). Campaigns at
    `e8940c22` (600 s each, sampled replays,
    every switch on): `boolean` clean, 940 runs, the slowest input 17 s
    under AddressSanitizer at load 11 to 15; `split` clean, 1,635 runs, none
    slow.
  * **DRAW survey of S9e.4b.3c.2 and S9e.4b.3c.3a (2026-10-05,
    `UPSTREAM_TESTS.md`).** At `b9c7ae1b` (`s9c2-kernel` with S9e.4b.3c.2,
    its `Arr::split_at` fix and S9e.4b.3c.3a; the public dataset, 120
    seconds a case, four at once). The 1,802 self-contained cases of the
    Boolean group on both backends: every status, every refusal's reason
    and every error the last survey's (`f4b584f7`) field for field, 987
    evaluating and registered, 592 refused, 223 unsupported on both, none
    failing or timing out. The 1,814 cases restoring a shape for a Boolean
    on the Rust adapter, and the 171 the import reaches on native DRAW too:
    native DRAW's every status and reason the last survey's; on the Rust
    adapter only 15 cases move, each as the two steps' trials found. 13
    evaluate on both backends with every check, S9e.4b.3c.2's (`so1` and
    `so2`: `bcommon_complex/B3`, `bcut_complex/C7`, `C8`,
    `bfuse_complex/B6`, `bugs/moddata_1/bug183_4`; `so2` and `so3`:
    `bcommon_complex/B4`, `bcut_complex/C9`, `D1`, `bfuse_complex/B7`;
    `so5` and `so2`: `bcommon_complex/B6`, `bcut_complex/D4`, `D5`,
    `bfuse_complex/B9`), each from `OutOfDomain("a vertex or a circle of
    both inputs on one sphere (S9e.4b.3c.2)")`; `bcut_complex/G4`
    evaluates on both backends with every check (S9e.4b.3c.3a's boss, from
    `OutOfDomain("an imported plane piece other than its primitive common
    its planes' half-spaces (S9e.4b.3c)")`); `bcut_complex/I6` moves from
    that reason to `Degenerate("an imported plane piece whose curved face
    is tangent to its plane faces")`, its tool's own tangency.
    `bugs/modalg_1/buc60926` stays `Degenerate("a plane through a cone's
    apex")` and the 11 `so` cases with `so6` or `so7` `Degenerate` by S9's
    rules. Of the 171: 55 evaluate on both backends, 69 are refused by S9's
    rules, 36 are bodies none of the kernel's constructions, 6 S9e.4b.4's,
    3 spline bodies (S9f) and 2 arguments of several solids; none is
    S9e.4b.3c's any more; 24 are unsupported natively too. The same 6
    cases time out at their first restores, and `bugs/modalg_1/
    buc60532_2` and `bugs/modalg_6/bug23585` end within the 120 s again
    (104 and 100 s), unsupported as before. The volume audit (`vprops` and
    `sprops` before each `checkprops`, both backends, the 1,028 audited
    before and the 14 new cases): both backends' values the last audit's
    bit for bit on all 1,028 (no moved value: S9e.4b.3c.2 and
    S9e.4b.3c.3a move no audited value), the statuses the same, the same 35
    disagreements. The 14 new cases' volume, area and centre agree with
    native DRAW's printed digits (`C8`'s cut empty on both; `G4`'s volume
    4,812,276.20 and area 193,592.99 against native's printed 4.81228e6 and
    193593). They are registered (`data`, `viewer_skipped` on both; 1,128
    cases), and the ledger does not change (`--ledger` holds). A full
    contract run of the manifest holds on both backends with the dataset
    (30 seconds a case), the 14 in 0.6 to 1.2 s on the Rust adapter
    (`bfuse_complex/B9` the slowest), the slowest Boolean case 4.4 s
    (`bopfuse_simple/ZK8`, on a host at load 5 to 8; `boptuc_simple/ZK8`
    5.3 in the last survey), the restore cases 0.1 to 3.6 s
    (`bcommon_complex/B9`), the rollex 0.5. No case fails, crashes or
    panics; no kernel change.
  * **S9e.4b.3c.3b refined, before its code (2026-10-05).** Why it is
    refused today: S9e.4b.3c.3a tries a body of one sphere, cylinder or cone
    face and plane faces as one Boolean of its primitive and the convex hull
    of its other planes (the common, the hull less the primitive, the bite,
    the boss) and refuses one no form's model matches
    (`OutOfDomain("an imported plane piece other than one Boolean of its
    primitive and its planes' hull (S9e.4b.3c.3b)")`). Replayed with debug
    assertions, the `boolean` fuzz target's `IMPORTED` stage refuses 4 first
    results so, two bodies twice each, read off the kernel's writer's files:
    a box `[0, 1.25] x [0, 2.25] x [0, 2.25]` fused with a ball's half (a
    ball of radius 2.8125 about `(0.375, -0.125, 1.125)` and the half-space
    below its equator's plane, its disc a face of the body, the plane cutting
    the box: its primitive common its disc's half-space fused with the box,
    no single form), and a U prism on a turned frame less a ball about its
    end's corner (three of its planes through the ball's centre, the ball
    through its whole thickness: the hull less the primitive, the hull not
    convex). The survey's restore cases hold none (its last S9e.4b.3c cases
    are `bcut_complex/G4`, evaluating since S9e.4b.3c.3a, and `I6` and
    `bugs/modalg_1/buc60926`, `Degenerate` for their own tangency and apex),
    so no DRAW case is expected to move. Decisions. (1) *Sub-steps.*
    **S9e.4b.3c.3b** (this step): a body of one curved face and plane faces
    as a Boolean tree of its primitive and the convex hulls of its planes, one
    level of pockets (2): a union of the primitive's common with its ends and
    a hull (the fuzz target's ball's half fused with a box, a sphere's disc
    with another form, a cylinder's flat along its axis on a boss), a
    primitive bitten twice, a boss, groove or bite on or by a body not convex
    in its planes (the U prism's groove), a primitive's end not normal to its
    axis; deeper trees (a pocket within a pocket, a bite by a body whose own
    pockets are not convex: a tooth left between a U's arms inside the
    primitive) stay refused, now as S9e.4b.4's, with bodies of several
    curved faces. (2) *The decomposition*, deterministic, read off the stored
    topology: (a) the curved face, its primitive as S9e.4b.3c.3a's (a whole
    sphere; a cylinder or a cone over the curved face's own axial range,
    past each end by a quarter of it but at a cap, clamped at a cone's
    apex), its material inside its quadric or outside it; (b) each edge of
    two faces convex or concave: at its curve's middle point, the first
    face's outward normal against the direction into the second face from
    the edge (the second's outward normal across its loop's tangent there),
    convex where it leans away from the first face's material; within
    `1e-9` of none (the faces tangent along it) the body's own tangency, as
    (5); (c) caps as S9e.4b.3c.3a's (plane faces normal to a cylinder's or a
    cone's axis at an end of the curved face's range, facing away from it,
    or into it outside the quadric) are the primitive's ends and in no
    group; (d) two groups of the other plane faces: the primitive's (faces
    sharing an edge with the curved face convex where its material lies
    inside the quadric, concave where outside: they trim the primitive) and
    the other (sharing a concave edge, convex outside: they bound what is
    joined to it or cut by it); a face sharing edges of both kinds with the
    curved face is decomposed by neither (refused, (5)); every other plane
    face joins, in rounds over the faces in order until none joins, the
    group of a face it shares a convex edge with, else of one it shares a
    concave edge with, a face reached from both groups at once in one round
    refused; (e) each group's region: its faces connected through concave
    edges between two of its own (convexity taken in the region's material:
    for the primitive's group outside the quadric the faces turned over and
    their edges' kinds reversed) form pockets, components of at least two
    faces in the order of their first faces, each the convex hull of its
    faces' planes turned over (S9e.4b.3c.3a's bite); its other faces' planes
    bound its hull (two faces on one plane facing one way one plane, as
    S9e.4b.3c.3a's), none for a group of pockets alone; its region the hull
    less its pockets in turn; (f) the body: with its curved face's material
    inside the quadric, the primitive common the primitive's group's region,
    fused with the other group's region; outside, the other group's region
    less that; an empty group drops its Boolean. (3) *The tree.* Every
    Boolean takes two exact models: the primitive's model (S9d's, the seams
    tried in turn as now), hull leaf models (`pieces.rs`'s `Hull`, bounded
    by S9e.4b.3a's cube) and the given models of inner Booleans. The
    primitive's group's region is applied to the primitive one hull at a
    time (its hull common, then each pocket cut), the other group's region
    built from its hull less each pocket, the last Boolean the fuse or the
    cut; at most four Booleans (a deeper tree is S9e.4b.4's). An inner
    Boolean's result is its given model (S9e.1's: the re-run's assembly its
    own slots, under ids of the piece's construction operation, no stored
    topology), arranged in the next Boolean as S9e.3a's given leaves are;
    the last Boolean's given model is matched to the stored topology by
    S9e.2's match (every stored vertex within the resolution of the
    model's, edges through their points, faces by their edges), its ids the
    stored ones, so a Boolean's history is over them as now; the stored
    vertices split rims (S9e.4b.3c.1's `split_at`) in every Boolean of the
    tree; the given model's frame the leaf holding the primitive's. Numbers:
    the engine's own, no new field or degree. (4) *Recognition.* S9e.4b.3c.3a's
    forms are tried first, unchanged (their one-Boolean models); a body none
    matches is decomposed, and its tree tried where it is not one of those
    forms (one Boolean of the primitive and one hull: already tried); the
    tree's model matching the stored topology is the body (the decomposition
    kept with the piece: a rigid motion reads it off the moved stored
    topology again), else `OutOfDomain("an imported plane piece other than
    a Boolean tree of its primitive and its planes' hulls (S9e.4b.4)")`, a
    Boolean's other refusal in the tree (S9's, raised in its own
    arrangement) reported. Classification: the tree's set functions over the
    primitive's location and each hull's sides within the resolution (a
    common the lesser, a fuse the greater, a cut the first's common the
    second's reversed); bounds and mass as S9e.4b.3c.3a's. (5) *Degenerate
    and refused.* S9's rules unchanged in every Boolean of the tree; a
    tangency in any of them, or an edge of the curved face along which a
    plane face is within `1e-9` of tangent to it, is the body's own,
    `Degenerate("an imported plane piece whose curved face is tangent to its
    plane faces")` as S9e.4b.3c.3a's; two faces on one plane facing apart
    within one group as S9e.4b.3c.3a's. Refused as S9e.4b.4's: a face of
    both kinds against the curved face or reached by both groups, more than
    four Booleans, and a tree whose model does not match (a pocket within a
    pocket, a pocket not convex). (6) *Fuzzing.* The `IMPORTED` stage
    unchanged: the replay measures the bodies now imported (the 4 above
    expected). (7) *Evidence first.* Bodies OCCT writes in one Boolean (a
    `write` block; a U prism is a prism of a profile of lines), every
    section a circle or a line (a sphere's planes its frame's parallels' or
    meridians'; a cylinder's oblique section an ellipse the reader does not
    take, so a primitive's oblique end is tested on the kernel's own
    topologies, imported directly as S9e.4b.3a's split pieces were):
    `form_u_scoop` (a U prism on `SKEW4` less a ball centred on its top face
    above its notch, the ball's axis along the notch's walls' normal: the
    top face a meridian plane, the walls parallels' planes, the floor out of
    reach), `form_u_boss` (a U prism on the world's axes fused with a
    cylindrical boss from one face to a cap), `form_bites` (a cylinder on the
    world's axes less a U prism along its axis between its caps whose two
    arms bite its wall: bitten twice, by a body not convex), `form_cap_boss`
    (a box fused with a ball's zone whose disc lies below the box, the ball
    meeting only the box's planes normal to its axis), `form_cap_pocket` (a
    box less a ball's zone above a parallel's plane inside it: a dish with a
    flat floor) and `form_dee_boss` (a box fused with a cylinder along `x`
    whose axis lies on its top face, flattened by a plane along its axis
    above it); declared `degenerate` `form_u_notch` (`bcut_complex/I6`'s
    notch on a body not convex: a U prism less a cylinder tangent to two of
    its faces) and declared `unsupported` `form_tooth` (a cylinder less a U
    prism whose notch lies inside the cylinder: a tooth left between the
    arms, a pocket within a pocket, S9e.4b.4's); cases of each against rods, slabs, balls and boxes, the
    body the tool, two imported and a chain; S9e.3a's chained reference on
    the constructions OCCT was given with S9e.4b.3c.3a's checks and each
    body's closed form where it has one (`generate_piece_trees_boolean_
    fixtures.py --check`, a CI group `piece-trees`,
    `test_piece_trees_boolean_reference.py`); a native capture keyed on the
    refusal it removes (`compare_piece_trees_boolean.py`, `imported.rs`'s
    refusal of a piece other than one Boolean of its primitive and its
    planes' hull); then the kernel and its tests (`tests/piece_tree_
    booleans.rs`: every fixture within `1e-9` enclosures, the degenerate
    refused, histories complete, deterministic, moved rigidly; the fuzz
    target's two bodies, a cylinder with an oblique end fused with a box and
    a ball's half with a box built by the kernel and imported, their
    Booleans the kernel's own results').
  * **S9e.4b.3c.3b evidence (2026-10-05).** Eight bodies OCCT wrote
    (`boolean-piece-trees-bodies.txt`, `write` blocks of one solid's rows, a
    `boolean` row and the other's rows, written by
    `compare_piece_trees_boolean.py --write-bodies` to
    `rust/fixtures/imported/form_*.brep`), each one Boolean of a primitive
    and a prism of lines (a box, or a U: a body not convex in its planes),
    every section a circle or a line: `form_u_scoop` (a U prism `12 x 8 x 6`
    on `SKEW4`, its slot `[4, 12] x [3, 5]`, less a ball of radius 5/2
    centred on the top face of one arm, the ball's frame as S9e.4b.3c.3a's
    scoop's: the top face a meridian plane, the arm's walls parallels'
    planes, the groove crossing the arm into the slot, the top face in two
    faces on one plane and the arms' ends two faces on one plane whose
    stored frames OCCT rounds apart; the fuzz target's ball's groove in a U),
    `form_u_boss` (a U prism on the world's axes fused with a cylindrical
    boss of radius 3/2 from its base's top face to a cap), `form_bites` (a
    cylinder of radius 3 on the world's axes less a U prism along its axis
    between its caps, its two arms biting the wall from outside, the seam
    between them: a primitive bitten twice, by a body not convex),
    `form_cap_boss` (a box on `SKEW` fused with an upper hemisphere of radius
    3 on its frame whose disc lies 1 below the box, the ball meeting the
    box's bottom alone, a parallel's plane; the fuzz target's ball's half
    with a box), `form_cap_pocket` (a box on `SKEW2` less such a hemisphere
    whose disc lies 2 inside it: a dish with a flat floor), `form_dee_boss` (a
    box on the world's axes fused with a prism along `x` of a disc of radius
    5/2 less a segment, its flat a plane along the axis with its ends
    rational, the box's top face through the axis, the arc's seam beyond the
    flat: a cylinder's end not normal to its axis), declared `degenerate`
    `form_u_notch` (a U prism less a cylinder of radius 3 along `z` about the
    middle of its base's end face, tangent to the two faces it meets) and
    declared `unsupported` `form_tooth` (a cylinder less a U prism whose base
    and slot lie inside it, a tooth left between the arms).
    `generate_piece_trees_boolean_fixtures.py --check`: 30 cases of 10
    groups (24 solid, 3 declared `degenerate`, 3 `unsupported`; 24 of class
    `pieces`, 3 `both`, 3 `chain`): `u_scoop_rod` (an upright rod through the
    groove), `pieces` (the hemisphere's boss and the flattened boss, both
    imported), `bites_box` (a `TILT` box across one bite), `rod_cap_pocket`
    (a `TILT` rod less the dish, the body the tool), `u_boss_slab` (a `TILT`
    slab across the boss and the slot), `cap_boss_rod` (a rod along `x`
    through the disc), `chain_dee` (the flattened boss less a rod, then with
    a `TILT` slab), `dee_boss_ball` (a ball across the flat); declared
    `unsupported` `tooth_rod` (a pocket within a pocket, S9e.4b.4's) and
    `degenerate` `u_notch_rod` (the notch, its wall tangent to its own
    faces). The reference is the constructions OCCT was given through
    S9e.3a's chained reference, with S9e.4b.3c.3a's checks relative to the
    case's size: the two families within 1.6e-35, each solid's closed form
    9.7e-41, each body's closed form 4.5e-37 (the U boss, the hemisphere's
    boss and dish on their boxes' stored axes, the flattened boss), the pair
    identities 1.4e-41 and the area identity 4.8e-41 (where no two inputs
    share a surface), Monte Carlo 2.5 standard errors (100,000 points a
    group), quadrature estimates 1.2e-32, every meeting's sine at least 0.11
    and events at least 1.9e-6 of their range apart outside the declared
    groups (and `rod_cap_pocket`'s events, below), solid counts by rays at
    two resolutions with S9e.4b.3c.1's join; the declared groups' checks
    kept apart within 1.1e-21. Every file read independently
    (`stored_records`): its faces one sphere or cylinder and planes, every
    stored vertex within 4.5e-16 of the size on the construction's surfaces.
    `test_piece_trees_boolean_reference.py` checks the hemispheres' and the
    flattened boss's closed forms against quadrature, the U boss's against
    its area, the chained reference on the flattened boss against its
    closed form, the reference's U prisms against OCCT's profiles (one area,
    the slot's box past the open end and both caps, a reflex corner), the
    reference's hemispheres, the frames keeping every section a circle or a
    line and every seam off its body, the bodies' declared refusals, the
    case list and its protocol rows and the files apart from every earlier
    set's. The generator's check is a CI group of its own (`piece-trees`,
    15 minutes on four workers locally, 11 under Python 3.12); Python 3.9 and 3.12 write the same
    files. Corrections from the evidence, amending the refined decisions'
    plan (7): (a) the reference models convex profiles and whole spheres
    only, so it takes each U prism as the box of its outer planes less its
    slot's box (past the open end and both caps) and each hemisphere as its
    ball common a cylinder from its equator's plane (twice its radius wide,
    four radii high): the same sets; the dish's cylinder's base, through the
    ball's rounded centre, lies within rounding of a family line of the
    dish's box's walls, so `rod_cap_pocket`'s events are found twice
    rounding apart there and its spacing is not checked (its sines and gaps
    are); (b) the reference's cost grows with a case's faces (the first set,
    with boxes and slabs as partners and a box less a slot in every U, ran
    over half an hour on four workers), so the partners are rods where they
    can be, the chain the flattened boss's and `pieces` the hemisphere's boss
    with the flattened boss; a rod across the bitten cylinder meets its
    bites' lines along its axis, irrational lines
    (`ComputationLimit("an irrational line against a cylinder")`,
    S9e.4b.3c.3a's), so its partner stays a box; (c) the U boss is 11 long:
    at 10 the arms' end plane was a face of S9e.4b.3c.3a's cube about its
    bite form's primitive, facing it, so the current import refused it as
    `Degenerate("an imported plane piece's two faces on one plane facing
    apart")` rather than as this step's; (d) partners moved: the flattened
    boss's ball off its box's face and its cap's plane (tangent), the
    hemisphere's rod from the box's normal (a plane along a cylinder's axis
    within rounding, S9's refusal) to the world's `x`, the U boss's slab
    thicker and the chain's slab higher (solids within two grid spacings at
    the coarse resolution), the bitten cylinder's box and the flattened boss
    off round coordinates (events of different features a rounding apart).
    The capture `occt-boolean-piece-trees-preimplementation`
    (`compare_piece_trees_boolean.py`, keyed on the refusal the step removes,
    `imported.rs`'s `OutOfDomain("an imported plane piece other than one
    Boolean of its primitive and its planes' hull (S9e.4b.3c.3b)")`; the
    kernel's probe `unsupported` on all 30, every body refused on import by
    that refusal): every result valid, 27 matching (within 2e-8 of the
    reference), 3 reviewed: the U's groove's, the flattened boss's ball's
    and the chain's commons, volumes up to 4.0e-7 relative and centres up to
    5.8e-8 of the case's size off by BRepGProp's default integration (a
    sphere's or a cylinder's faces met by planes and other quadrics), within
    4.8e-9 measured adaptively at 1e-10 and 1e-12 by a diagnostic build.
    S9e.4b.3c.3b's kernel next.
  * **S9e.4b.3c.3b implemented** (`solid/imported.rs`'s `Tree`: the bends of
    the stored edges, the faces' groups, each group's hull and pockets, the
    tree, its classification; `curved/pieces.rs`: the tree's Booleans in turn,
    an inner one's assembly named under an operation of its own and given to
    the next as its given model; `curved/given.rs`: `built_on`, a given model
    on a stored topology or an inner assembly, its frame its first leaf's
    holding the primitive; `topology.rs`: an inner assembly's unchecked
    topology), as the refined decisions describe: a body of one sphere,
    cylinder or cone face and plane faces that no form of S9e.4b.3c.3a's
    matches is decomposed into a Boolean tree of its primitive and the convex
    hulls of its planes, the first tree whose given model matches its stored
    topology the body, none matching `OutOfDomain("an imported plane piece
    other than a Boolean tree of its primitive and its planes' hulls
    (S9e.4b.4)")` (the forms' mismatch too). All 30 fixtures as declared (24
    within the kernel's enclosures, each at most `1e-9` wide; the U's notch's
    3 refused as `Degenerate("an imported plane piece whose curved face is
    tangent to its plane faces")`, read off its edges' bends before any tree;
    the tooth's 3 as S9e.4b.4's), every history complete over the imported
    bodies' stored ids, the U's two arm ends on one plane each continuing
    apart, results deterministic and moved rigidly, both inputs translated and
    turned keeping the reference's volumes, every body its tree (its volume
    its closed form, points in and off it classified, its stored vertices on
    its boundary), and the kernel's own bodies of trees imported, their
    Booleans with a turned box the kernel's own results': the `boolean` fuzz
    target's two (a box fused with a ball's half whose disc plane cuts it, a U
    prism less a ball about its end's corner) and a rod bitten twice by a U,
    written by its writer and read back, and a cylinder's piece of an oblique
    plane fused with a box (its end an ellipse, which the reader does not
    take), its topology imported directly (`tests/piece_tree_booleans.rs`,
    8.1 s in release, 9.0 s at `opt-level` 2 with debug assertions, 8.7 s in
    release with debug assertions under the emulated correctly rounded
    `hypot`, where the pieces' earlier test files pass too). Amendments, from
    the implementation: (a) the groups: an edge of the kind no seam between
    the groups takes (convex, but concave for the primitive's group outside
    the quadric: within a hull, or between a hull and its pocket) keeps its
    faces in one group, every other face joining through such edges; a
    component of faces met only across the other kind (a pocket's own edges,
    or the seam where the groups are fused or cut) is tried in the other group
    than the faces it meets and in the same, in turn (the same first where it
    meets the other group or a face of a pocket; at most three such
    components), the match the arbiter. The refined decisions' rounds through
    convex then concave edges put the fuzz target's box, whose faces meet only
    the ball's half's disc, with the disc; (b) inside the quadric the union is
    tried two ways, the primitive common its group's region fused with the
    other's, and the primitive common both regions' union, the latter first
    where no face of the other group meets the curved face: the fuzz target's
    box lies within the ball, its bottom face hidden, so its hull is open
    below and only the primitive bounds it; (c) planes of faces within the
    resolution of one plane facing one way (a U prism's two arm ends, whose
    stored frames OCCT rounds apart) take the first's frame, so the hull takes
    them as one plane (S9e.4b.3c.3a's exact test of one plane missed them, and
    the U's groove was a face too many); (d) only the first tree's refusal
    other than a mismatch is reported, a later choice's being a tree the body
    is not (a later tree of the tooth's met S9's `Degenerate("a vertex of one
    input on the other's face")`); (e) the tree's hulls are bounded by a cube about the
    primitive whose half side reaches past every stored vertex, twice and
    `3/8` more (S9e.4b.3c.3a's cube about the primitive alone met the U boss's
    arms' end plane), each hull's and inner Boolean's entities under an
    operation of their own. `compare_piece_trees_boolean.py` 22 matches and 8
    reviewed (the 3 captured and, with the kernel, entity counts: the kernel's
    meetings of a cylinder with a ball split at their turning points into
    S9d.2's graphs where OCCT's are one approximated edge, OCCT's unified
    result keeping a boss's cylinder's seam edge), every enclosure within the
    reference; every other comparison of `HANDOFF.md`'s table unchanged with 0
    failures (and `compare_split.py` 72/56, `compare_brep.py --family spline`
    10/3); the release suite (628 tests) and the tools' unit tests (344)
    passing, the ledger unchanged. The `boolean` fuzz target's `IMPORTED`
    stage: replayed with debug assertions, 92 first results of the corpus's
    1,467 inputs and regressions reach it, 23 imported and cut within the
    chained cut's volume (the 4 refused as this step's among them: the box
    with a ball's half and the U's groove, twice each), 63 not written by the
    kernel's writer, 6 read back with pcurves off their edges, none refused.
    Replays with debug assertions, one process an input, natively and under
    the emulated correctly rounded `hypot`: the boolean corpus and its
    regressions (1,467 inputs) and the split corpus and its regressions
    (3,570), the slowest 6.6 s and 3.2 s (6.3 s and 2.7 s emulated), no
    failure but two inputs another campaign added to the local corpus after
    this branch's base, which fail on that base alike (a boolean input's
    `coordinates cannot resolve the requested linear tolerance`, natively and
    emulated; a split input's `a piece below its plane`, emulated only). A
    trial of the DRAW survey's 171 restore cases the import reaches, on the
    Rust adapter and natively (not the survey: nothing registered): no case
    moves with this step, every status and reason S9e.4b.3c.3a's trial found
    (55 evaluating on both backends, S9e.4b.3c.2's 13 and `G4` among them;
    `I6` refused for its tool's own tangency, `buc60926` at its frustum's
    apex), none refused as S9e.4b.3c's; native DRAW's statuses as before.
    Pending: the capture's Linux record. DRAW survey: that of S9e.4b.3c.3b,
    S9e.4b.4a and the scheduled replay's fixes, after S9e.4b.4a's
    implemented bullet (no case moving with this step, as its trial
    found). Campaigns at `699b9b85` (600 s each, sampled replays, every
    switch on): `boolean` clean, 978 runs, the slowest input 22 s under
    AddressSanitizer; `split` clean, 1,241 runs, the slowest 16 s.
  * **S9e.4b.4 refined, before its code (2026-10-05).** Why each class is
    refused today, from S9e.4b.3c.3b's trial of the 171 restore cases the
    import reaches and the files they restore (read by the converter; the
    dataset is not committed). (a) *Joints of two circles* (6 cases, each
    `OutOfDomain("an imported prism's arcs of two circles meeting at a joint
    (S9e.4b.4)")`, S9e.4b.1's `snapped::points` in the Boolean's exact model):
    `bcut_complex/E8` and `bfuse_complex/D5` (`CTO900_jap60038-part`,
    `CTO900_ksi0014a`: the common of four discs of radius 80 about `(0, +-50)`
    and `(+-50, 0)` on the world's axes, its corners the circles' irrational
    common points; the tool a cylinder of radius 80 on one of those circles,
    below the part's bottom face), `bcut_complex/P4` (`CTO909_tool_2`, a
    profile of lines and fillets whose arcs of radii 2 and 0.47 meet tangent
    from either side, inside the box `CTO909_part_1`), `bfuse_complex/E1`
    (`CTO900_fra11018a`, an arc of radius 10 tangent inside one of radius 26.6
    at their joint) and `bugs/modalg_2/bug4993_1`, `_2` (`OCC4993-s1` and
    `-s2`, arcs of radii 32 and 4.5 crossing at 0.011 rad): S9e.4a's prism,
    whose stored joint, rounded once in its cap's frame, lies off both circles
    (by up to 9.2e-14 in these files), and whose circles' common point is a
    quadratic surd that a profile segment's rational ends do not hold. (b)
    *Bodies of several curved faces* (10): `bcut_complex/G9` and
    `bugs/modalg_2/bug417` (`cts21128c`, a frustum of half-angle 0.145 rad
    from `z = 0` common a ball of radius 58.9 about `(0, 0, -45.7)`;
    `cts21128d`, a rod of radius 2.95 between a plane and the same sphere: a
    pin drilled out of the dome, both inputs' faces on one sphere), and
    `bugs/modalg_2/bug476_1` to `_8` (`OCC485a`, a turned part of two
    cylinders, three cones and two tori, smoothly joined, whose stored axes
    lean 1.1e-8 apart and whose origins drift up to 4.4e-7, two of its cones
    of half-angles within 1.3e-6 of a right angle; against `OCC485b`, a
    triangular prism). (c) *Walls of two directions* (3): `bfuse_complex/E5`
    (`CTO900_pro9476-part`, a stepped shaft: two coaxial cylinders along `x`
    of radii 118.4 and 419.4 between three discs), `bugs/modalg_6/bug28773`
    (`bug28773_2`, a stepped shaft of five coaxial cylinders, against `_1`, of
    two) and `bfuse_complex/K1` (`CTO904_cts20370-part`, a box with its four
    vertical edges rounded and a cylindrical boss of radius 40 along `y`);
    `bcut_complex/G4`'s part, listed among these in S9e.4b refined, is
    S9e.4b.3c.3a's boss. (d) The other 23 bodies none of the kernel's
    constructions are the hollow sphere's, in Booleans with wires (no solid
    argument; unsupported natively too). (e) From the fixtures and the fuzz
    target: an imported polyhedron against curved faces (one first result of
    the `boolean` target's replay) or with a cavity (S9e.4b.2's `hollow`), and
    a pocket within a pocket (S9e.4b.3c.3b's `tooth`). Decisions. (1)
    *Sub-steps*, ordered by the cases each unlocks and the machinery each
    takes: **S9e.4b.4a** (this step): an imported prism whose arcs of two
    circles meet at a joint, (a), on S9e.4a's construction with each such arc
    taken through its two ends (2): the prism model alone, no new field or
    degree (the S9e text's plan, a joint the circles' common point in a
    quadratic field, would carry surds through every profile test and every
    pierce of the vertical edge there, nested under a partner cylinder's);
    **S9e.4b.4b**: bodies of several primitives' curved faces, (c) and G9's
    pair: S9e.4b.3c.3b's tree with a leaf per curved face's primitive (a
    stepped shaft the fuse of coaxial prisms, a boss along another axis the
    fuse of two prisms, the dome a frustum common a ball, the pin a rod common
    that ball), its faces grouped by their primitives' surfaces, with
    S9e.4b.3c.1's faces on one sphere between the inputs; **S9e.4b.4c**:
    deeper trees (a pocket within a pocket) and an imported polyhedron against
    curved faces or with a cavity (S9b.2's stored model of exactly planar
    triangles as a leaf of the curved engine, membership by parity);
    **S9e.4b.4d**: turned bodies of smooth joins and nearly degenerate
    surfaces (`OCC485a`): the S9e text holds stored faces tangent along an
    edge `Degenerate` until such edges keep exact data, and a cone within
    1.3e-6 of a disc is no exact model's, so they stay refused until a
    decision of their own. The reader's version-3 header check stays open on
    the import track (none of the 171 restores such a file). (2) *S9e.4b.4a's
    representation.* A *crossing joint*: arcs of two circles (centres or radii
    different in the cap's frame) meeting at a point that, rounded once, lies
    off either exactly. An arc ending at a crossing joint whose circle no
    other arc of its path shares is *taken through its ends*: in the exact
    model its circle is the one through its two model ends `a` (where it
    starts) and `b`, of centre `a + rho e` and radius `rho`, `e` the rational
    unit vector nearest the stored centre's direction from `a` (`(+-(1 - s^2),
    2 s) / (1 + s^2)`, the sign the direction's `x`'s, `s` the binary64
    rounding of `d_y / (r + |d_x|)` computed exactly: S9e.4b.1's half-angle
    tangent) and `rho = |b - a|^2 / (2 (b - a) . e)` exactly (a chord's
    projection on the radius at its end is half its square over the radius).
    Every circle through two rational points with a rational centre and radius
    has that form, and `(b - a) . e > 0` for any chord seen from its end, so
    `rho` exists and is positive. A crossing joint of such an arc with a line
    or with another arc taken through its ends is its rounded point (as two
    lines' joint); with an arc kept (an arc meeting only lines, or one of
    several arcs on one circle) that circle's point (S9e.4b.1's `onto`), the
    same for both. Numbers: the model's centre and radius rationals of a few
    hundred bits (`e`'s 106-bit terms, `rho` over their products), its
    vertices binary64 or S9e.4b.1's points: no new field or degree. The circle
    moves by at most about `(|a - p| + |b - q|) / sin(theta / 2)` plus `r
    2^-52` for an arc of angle `theta` whose stored ends `p` and `q` lie off
    its stored circle: checked exactly, its centre and radius within the
    body's resolution of the stored ones, else `Degenerate("an imported
    prism's arc too short to take through its ends")`. A joint tangent in the
    body OCCT was given (a fillet chain, `P4`'s and `E1`'s) is tangent in no
    rounded data: the two circles cross at the joint within rounding of
    tangency and again within rounding of it on one circle's continuation,
    outside that arc, and the arrangement decides the joint by its exact turn,
    as S9e.4b.1's line tangent to its arc. The flag asking for it is
    S9e.4b.1's (`Profile::rounded_arcs`); the code stays `snapped.rs`'s, its
    refusal going. (3) *Degenerate and refused.* S9's rules unchanged. A
    partner's surface on an arc's stored circle (`D5`'s tool, a cylinder on
    the part's circle below its face) lies within rounding of the circle taken
    through the arc's ends, not on it: two surfaces within the resolution of
    one, which the arrangement refuses where it meets them (the declared case
    below); `D5` is expected refused. The arc too short, as (2). Refused: a
    crossing joint of two arcs whose circles each hold other arcs of the path
    (a lens whose arcs OCCT split, none among the surveyed files),
    `OutOfDomain("an imported prism's joint of two circles each holding
    several arcs (S9e.4b.4)")`; a kernel profile's arc ending off its circle
    stays refused by design (S9c). (4) *Its other queries* are S9e.4b.1's:
    classification by the construction, mass by the stored topology, a rigid
    motion rebuilding the construction whose arcs are taken through their ends
    again in the moved frame, a result given to another Boolean re-running the
    same model. (5) *Fuzzing.* No profile of the `boolean` target has a joint
    of two circles (S9e.4b.1's replay met none), so the `IMPORTED` stage gains
    a lens (two arcs of circles about `(+-3k, 0)` of radius `5k` meeting at
    `(0, +-4k)`, `k` dyadic) in the object's frame given the chosen operation
    with the tool, then written, read back, imported and given it again, both
    volumes equal within `1e-9` (`JOINTS`, by the chained byte's next bit; in
    the tilted frame the joints round off both circles once read back). (6)
    *Evidence first.* Bodies OCCT writes (`MakePrism` of profiles in turned
    frames, every joint of two circles a rational common point in the
    profile's frame, so the construction is exact and the stored joints its
    roundings): `quad` (E8's and D5's shape, four discs' common), `ogee`
    (P4's, two arcs tangent from either side between lines), `cam` (E1's, an
    arc tangent inside another), `blade` (bug4993's, two circles crossing at a
    small angle) and, declared `unsupported`, `split_lens` (a lens each of
    whose arcs is split in two); cases against boxes, a rod, a ball and a
    `TILT` slab, as object and as tool, two imported, a chain, and declared
    `degenerate` a rod on one of the quad's circles below its bottom face
    (D5's configuration); the reference the constructions OCCT was given
    through S9e.3a's chained reference with S9e.4a's checks, each file's
    stored vertices on the construction's surfaces and its crossing joints off
    both circles in the construction's frame once rounded
    (`generate_imported_joints_boolean_fixtures.py --check`, a CI group
    `imported-joints`, `test_imported_joints_boolean_reference.py`); a native
    capture keyed on the refusal it removes
    (`compare_imported_joints_boolean.py`, `snapped.rs`'s "an imported prism's
    arcs of two circles meeting at a joint"); then the kernel and its tests
    (`tests/imported_joint_booleans.rs`: every fixture within `1e-9`
    enclosures, the degenerate refused, histories over the stored ids,
    deterministic, moved rigidly), S9e.4b.1's `lens` cases solid with it, and
    the DRAW trial of the 6.
  * **S9e.4b.4a evidence (2026-10-05).** Five bodies OCCT wrote
    (`boolean-imported-joints-bodies.txt`, `MakePrism` of profiles whose world
    coordinates are their turned frames' roundings, written by
    `compare_imported_joints_boolean.py --write-bodies` to
    `rust/fixtures/imported/`), every joint of two circles a rational common
    point of both in the profile's frame: `quad`, the common of four discs of
    radius 5 about `(+-1, 0)` and `(0, +-1)`, its corners `(+-3, +-3)` (`E8`'s
    and `D5`'s part), in the `R125` frame; `arch`, a base from `(-2, 0)` to
    `(2, 0)` under two arcs of radius 5 about `(-+3, 0)` meeting at `(0, 4)`
    (a crossing joint whose arcs end at a line too), in the `TURN30` frame;
    `cam`, an arc of radius 3 tangent inside one of radius 5 at `(5, 0)`
    (`E1`'s), in the `R125` frame; `blade`, arcs of radii 25/2 and 13/2
    crossing at 0.11 rad at `(12, 7/2)` (`bug4993`'s), in the `TURN30` frame;
    and `split_lens`, S9e.4b.1's lens with each arc split in two, in the
    `TURN30` frame. `generate_imported_joints_boolean_fixtures.py --check`: 39
    cases of 13 groups (33 solid, 3 declared `degenerate`, 3 `unsupported`; 33
    of class `prism`, 3 `both`, 3 `chain`): the quad against a box across its
    corner, as object and as tool, a rod through its arc wall and a ball about
    its corner's vertical edge; the arch against a box across its apex and a
    `TILT` slab (two solids cut); the cam against a box across its tangent
    joint and a rod through its small arc's wall; the blade against a box
    across its crossing; the quad and the arch both imported; the quad less
    the rod, then with the box; declared `degenerate` `quad_seat`, a rod in
    the quad's frame on the circle of its arc about `(-1, 0)` below its bottom
    face (`D5`'s configuration); declared `unsupported` the split lens against
    a box (a joint of two circles each holding several arcs, S9e.4b.4's). The
    reference is the constructions OCCT was given through S9e.3a's chained
    reference with S9e.4a's checks (`generate_imported_boolean_fixtures.
    evaluate_chain`), relative to the case's size: the two families within
    2.1e-36, each solid's closed form 5.7e-41, the pair identities 3.4e-42 and
    the area identity 2.6e-41, Monte Carlo 3.3 standard errors, quadrature
    estimates 8.5e-33, solid counts by rays at two resolutions, every
    meeting's sine at least 0.47 and events at least 1.6e-3 of their range
    apart outside the declared groups. Every file read independently
    (`stored_records`): its faces' kinds the construction's, every stored
    vertex within 1.9e-15 of the size on one of the construction's walls and
    one of its caps (OCCT's 15 digits; a vertex on a cap's plane exactly no
    longer hides its wall's distance in one minimum), and every joint of two
    circles, its stored vertices' local coordinates in the construction's
    frame rounded once, off one of its circles at least: all 18 (the quad's 8,
    the arch's, the cam's and the blade's 2, the split lens's 4).
    `test_imported_joints_boolean_reference.py` checks every construction's
    arcs ending on their circles and its joints on both exactly, the joints'
    kinds (the cam's circles touching inside, the arch's, the blade's and the
    quad's crossing at `asin(24/25)`, `asin(36/325)` and `asin(7/25)`), the
    chained reference on each body alone against its profile's closed form (a
    polygon and its arcs' segments, times the stored axes' determinant), the
    quad in common with the box against a direct quadrature of its sections
    within 1e-12, the off-joint test, and the case list and its protocol rows.
    The generator's check is a CI group of its own (`imported-joints`, 4.8
    minutes on four workers locally under Python 3.9, 3.2 under 3.12); Python
    3.9 and 3.12 write the same files. Corrections from the evidence, amending
    the refined decisions' plan (6): (a) the chained reference takes convex
    profiles only (counter-clockwise arcs), so the S curve of two arcs
    touching from either side (`P4`'s fillets, `ogee`) is no body of this set:
    the kernel's tests write such a body of its own and read it back, and the
    DRAW trial reads `P4`'s; the `arch` takes its place; (b) the arch's slab
    moved off its base corner (a slab plane through a stored vertex: events
    3.7e-17 of their range apart); (c) a trial of the cases on a draft of the
    step's code (not committed) put every solid case within the reference but
    evaluated the declared `quad_seat`: the rod's circle and the quad's arc's
    circle taken through its ends lie a rounding apart, and the walls leave a
    sliver between them at the shared cap plane that none of the curved
    engine's rules finds (no two vertices within the resolution), where S9a's
    one frame refuses such a pair as "a boundary within the resolution of the
    other profile's" and S9e.4b.3c.1's spheres as two spheres within the
    resolution of one; the step's kernel takes that rule for parallel circular
    cylinders (refined decisions (3)), which refuses `quad_seat` and, in the
    same draft, DRAW's `D5` (otherwise left to the validator's undecided loop
    winding, a `ComputationLimit`). The capture
    `occt-boolean-imported-joints-preimplementation`
    (`compare_imported_joints_boolean.py`, keyed on the refusal the step
    removes, `snapped.rs`'s `OutOfDomain("an imported prism's arcs of two
    circles meeting at a joint (S9e.4b.4)")`; the kernel's probe `unsupported`
    on all 39, every body refused by it in the Boolean's exact model): every
    result valid, 36 matching (volumes within 2.3e-10, areas 8.8e-11, centres
    3.2e-11 of the size: OCCT joins the declared `quad_seat`'s rod and arc
    walls, its fuse one solid with the reference's measures), 3 reviewed: the
    quad's ball cases, BRepGProp's default integration on a sphere's faces met
    by cylinders (up to 7.0e-7 relative; within 2.3e-9 measured adaptively at
    1e-10 and 1e-12 by a diagnostic build), in
    `occt-boolean-imported-joints-divergences.json`; 5 results' counts change
    when unified. S9e.4b.4a's kernel next.
  * **S9e.4b.4a implemented** (`solid/boolean/curved/snapped.rs`'s `path`,
    `model.rs` taking each arc's circle from it, `meet.rs`'s `cyl_pair`), as
    the refined decisions describe: an arc of an imported prism ending at a
    joint where arcs of two circles meet off either, whose circle no other arc
    of its path shares, is taken through its two ends in the exact model (its
    centre `a + rho e` from its start, `e` the rational unit vector nearest
    the stored centre's direction by its half-angle tangent rounded once, `rho
    = |b - a|^2 / (2 (b - a) . e)` exactly), the joint its rounded point or a
    kept circle's, the circle within the resolution of the stored one. All 39
    fixtures as declared (33 within the kernel's enclosures, each at most
    `1e-9` wide; `quad_seat`'s 3 refused as `Degenerate("two cylinders within
    the resolution of one cylinder")`; the split lens's 3 `OutOfDomain` as
    S9e.4b.4's), every history complete over the imported bodies' stored ids,
    results deterministic and moved rigidly, both inputs translated (and
    turned where no face of one input is exactly parallel to the other's
    cylinder) keeping the reference's volumes, every body its construction
    (its closed-form volume, points classified, its stored vertices on its
    boundary), and the kernel's own prisms whose arcs of two circles meet at
    rational points in turned frames (four discs' common, a lens in the tilted
    frame, an S curve of two arcs touching from either side, an arc tangent
    inside another) written by its writer, read back and imported, each
    Boolean with a box and a rod its own result's volume
    (`tests/imported_joint_booleans.rs`, 10.9 s at `opt-level` 2 with debug
    assertions and 9.8 s in release on a host at load 23 to 29, 7.3 s in
    release with debug assertions under the emulated correctly rounded
    `hypot`, where S9e.4b.1's and S9e.4a's test files pass too). S9e.4b.1's
    `lens` cases are solid, declared so (its generator's group, its tests).
    `compare_imported_joints_boolean.py` 36 matches and 3 reviewed (the native
    measures and, with the kernel, entity counts: a ball's meetings and sphere
    face split at the kernel's own points and at OCCT's), every enclosure
    within the reference; `compare_imported_arcs_boolean.py` 27 and 9, the
    kernel within the reference on all 30 solid cases; every other comparison
    of `HANDOFF.md`'s table unchanged with 0 failures (and `compare_split.py`
    72/56, `compare_brep.py --family spline` 10/3, `compare_brep_io.py`
    6,835/7, `compare_step.py` 23/6 on STEP-b's SDK); the release suite (640
    tests) and the tools' unit tests (350) passing, the ledger unchanged.
    Amendments, from the implementation: (a) the refined decisions' (3)
    expected the arrangement to refuse a partner's cylinder on an arc's stored
    circle, but the curved engine had no rule for it (S9a's one frame refuses
    two boundaries within the resolution, S9e.4b.3c.1 two spheres within the
    resolution of one): two parallel circular cylinders within the resolution
    of one and not one, their faces' boxes meeting, are now `Degenerate("two
    cylinders within the resolution of one cylinder")` (`cyl_pair`), refusing
    `quad_seat` and DRAW's `D5` where a sliver between the walls or the
    validator's undecided winding was left; no case of the comparisons, the
    suite or the replays moves with it; (b) `onto` keeps S9e.4b.1's points bit
    for bit, its rational unit vector now `unit`'s, shared with `through`.
    Open: where a partner holds an arc's stored circle exactly (`E8`'s and
    `D5`'s tools, on the world's axes), keeping that circle and taking the
    other arc at each of its joints through its ends would decide the two
    walls on one surface; the choice needs the partner in the prism's model,
    left open (the trial's three such cases refused by (a)). The `boolean`
    fuzz target's `JOINTS` stage (a lens of two arcs of circles about `(+-3k,
    0)` of radius `5k` in the object's frame, by the chained byte's next bit,
    given the chosen operation, then written, read back, imported and given it
    again): replaying the corpus and the regressions with debug assertions
    (1,468 inputs), 480 reach it, 298 evaluate both ways with equal volumes
    (96 in the tilted frame, where the joints round off both circles once read
    back), the others refused by S9's rules on one side or both (most of them
    the reimported frames, normalized again, within rounding of the tool's
    directions); no failure, the slowest input 10.9 s, and under the emulated
    `hypot` 15.3 s on a host at load 20 to 29; the split corpus and its
    regressions (3,571 inputs) clean natively and emulated, the slowest 3.3 s
    and 1.6 s. A trial of the DRAW survey's 171 restore cases the import
    reaches, on the Rust adapter and natively (not the survey: nothing
    registered): only the 6 of two circles at a joint move, each from "an
    imported prism's arcs of two circles meeting at a joint (S9e.4b.4)":
    `bcut_complex/P4` (a fillet chain's tool, arcs touching from either side,
    inside its box) evaluates on both backends with every check;
    `bcut_complex/E8` and `bfuse_complex/D5` (the four discs' part against a
    cylinder on one of its circles below it) and `E1` are refused as two
    cylinders within the resolution of one cylinder (each tool on a stored
    circle of the part's arcs taken through their ends; `E1`'s tangency
    besides), `bugs/modalg_2/bug4993_1` and `_2` as two faces within the
    resolution of one plane; native DRAW's statuses as before. Of the 171, 56
    evaluate on both backends. Pending: the capture's Linux record (the
    campaigns below). DRAW survey: that of S9e.4b.3c.3b, S9e.4b.4a and the
    scheduled replay's fixes, below (`P4` registered; the parallel
    cylinders' rule moving no registered case).
    Campaigns at `d1869f2f` (600 s each, sampled replays, every switch on):
    `boolean` clean, 1,025 runs, the slowest input 19 s under
    AddressSanitizer (its startup replay 1,601 s; at `64673ebd` it overran
    the hour's budget at load 13 to 16, no input failing); `split` clean,
    1,972 runs, none slow.
  * **DRAW survey of S9e.4b.3c.3b, S9e.4b.4a and the scheduled replay's
    fixes (2026-10-05, `UPSTREAM_TESTS.md`).** At `16121052`
    (`s9c2-kernel` with the fixes of CI's scheduled full replay,
    S9e.4b.3c.3b and S9e.4b.4a; the public dataset, 120 seconds a case,
    four at once, on a host at load 10 to 64). The 1,802 self-contained
    cases of the Boolean group on both backends: every status, every
    refusal's reason and every error the last survey's (`b9c7ae1b`) field
    for field, 987 evaluating and registered, 592 refused, 223 unsupported
    on both, none failing or timing out. The 1,814 cases restoring a shape
    for a Boolean on the Rust adapter, and the 171 the import reaches on
    native DRAW too: native DRAW's every status and reason the last
    survey's; on the Rust adapter only S9e.4b.4a's 6 cases move, each from
    `OutOfDomain("an imported prism's arcs of two circles meeting at a
    joint (S9e.4b.4)")` as its trial found: `bcut_complex/P4` evaluates on
    both backends with every check; `bcut_complex/E8`, `bfuse_complex/D5`
    and `E1` are `Degenerate("two cylinders within the resolution of one
    cylinder")`; `bugs/modalg_2/bug4993_1` and `_2` `Degenerate("two faces
    within the resolution of one plane")`. S9e.4b.3c.3b and the scheduled
    replay's fixes move no case. The parallel cylinders' rule
    (S9e.4b.4a's amendment (a), `meet.rs`'s `cyl_pair`): no registered
    case and no self-contained case's status or reason moves with it, and
    of all the surveyed cases it refuses only `E8`, `D5` and `E1`, none
    registered and none evaluating before. Replayed at `16121052` with
    the rule disabled (a scratch build, no kernel change): `E8` evaluates
    with native DRAW's volume and area to its printed digits (4,697,379.48
    and 173,509.85, native's 4.69738e6 and 173510; 3.9 s), so the rule
    costs it a result, the open item's case (its tool holds a stored
    circle of the part exactly); `D5` is `ComputationLimit("a face's loop
    winding the validator leaves undecided")` and `E1` `Degenerate("a
    tangency between the inputs (S9c)")`, so the rule only renames their
    refusals; `P4` evaluates either way. Of the 171: 56 evaluate on both
    backends, 74 are refused by S9's rules, 36 are bodies none of the
    kernel's constructions, 3 spline bodies (S9f) and 2 arguments of
    several solids; none is S9e.4b.4's any more; 24 are unsupported
    natively too. Six cases time out at their first restores again, the
    load moving two of the reader's slowest: `bugs/modalg_1/buc60532_2` and
    `bugs/modalg_6/bug23585` time out where the last survey's ended in 104
    and 100 s, `buc60532` and `buc60532_1` end in 114 s where they timed
    out; alone at load 6 to 10 the four end in 88 to 117 s, unsupported as
    before (the reader takes no `SurfaceOfLinearExtrusion`; `tolerance`).
    The volume audit (`vprops` and `sprops` before each `checkprops`, both
    backends, the 1,042 audited before and `P4`): both backends' values the
    last audit's bit for bit on all 1,042 (no moved value), the statuses
    the same, the same 35 disagreements. `P4`'s volume 88,171,481.37, area
    1,642,145.80 and centre agree with native DRAW's printed digits
    (8.81715e7, 1.64215e6). It is registered (`data`, `viewer_skipped` on
    both; 1,129 cases), and the ledger does not change (`--ledger` holds).
    A full contract run of the manifest holds on both backends with the
    dataset (30 seconds a case), `P4` in 3.0 s on the Rust adapter, the
    slowest Boolean case 5.8 s (`bopcut_simple/ZK8`, on a host at load 10
    to 45; `bopfuse_simple/ZK8` 4.4 in the last survey), the restore
    cases 0.2 to 5.1 s (`bcommon_complex/B9`), the rollex 0.6 to 0.7. No
    case fails, crashes or panics; no kernel change.
  * **S9e.4b.4b refined, before its code (2026-10-05).** Why S9e.4b.4's
    bodies of several curved faces are refused today, read off the survey's
    files by the converter (the dataset is not committed): the import takes
    S9e.4a's constructions and, for one sphere, cylinder or cone face among
    planes, S9e.4b.3's pieces; a body whose curved faces lie on two
    surfaces or more is neither, `OutOfDomain("an imported solid other than
    a prism, a sphere, a cone or a torus (S9e.4b)")` (`imported.rs`'s
    `general`). `bfuse_complex/E5`'s part (`CTO900_pro9476-part`) is a
    stepped shaft along the world's `x`: a cylinder of radius 419.41 from
    its disc `x = 0` to a disc `x = 138.37` holding a ring, and a coaxial
    one of radius 118.38 from that ring to a disc `x = 595.65` (rings, no
    vertex), fused with a box (`-tool`, S9e.4a's prism) whose face lies on
    the step's plane. `bugs/modalg_6/bug28773`'s `_2` is six coaxial
    cylinders along `x` and seven discs: a cylinder of radius 1 over `[1,
    2]`, and over `[0, 1]` a tube of radii 0.025 and 0.03, a tube of radii
    0.01 and 0.02 and a rod of radius 0.006 joined to it (three of its
    cylinders holes, its discs at `x = 1` stored on frames OCCT rounds
    1e-17 apart); its `_1`, cut from it, is a tube of radii 0.006 and 0.01
    over `[-0.01, 1]` lying in `_2`'s gap, on its cylinders (S9e.4a's prism
    on its bottom disc's stored frame, whose normal leans 2.2e-33 off the
    walls' `x`). `bfuse_complex/K1`'s part (`CTO904_cts20370-part`) is the
    box `[0, 200]^3` with its four edges along `z` rounded by cylinders of
    radius 50 tangent to its faces, less a cylinder of radius 40 along `y`
    through it (its tool a rod of radius 40 along `x`, its cap on the box's
    face `x = 200`). `bcut_complex/G9` and `bugs/modalg_2/bug417`
    (`cts21128c` and `d`): a frustum of half-angle 0.145 rad narrowing up
    from its disc `z = 0` of radius 14.8 common a ball of radius 58.9 about
    `(0, 0, -45.7)` on its axis (a dome, its pole a vertex loop), less a
    rod of radius 2.95 from the plane `z = 0` common the same ball (a pin,
    its disc inside the dome's, its rim on the dome's sphere split at its
    cylinder's and its sphere's seams by two stored vertices). Decisions.
    (1) *Sub-steps*, by what each unlocks and the machinery each takes:
    **S9e.4b.4b.1** (this step): a body of several sphere, cylinder and
    cone faces every plane face of which is an end of their primitives
    (S9e.4b.3c.3a's caps), a Boolean chain of its primitives (2, 3): E5's
    shaft, `bug28773`'s comb and G9's dome and pin; **S9e.4b.4b.2**: such a
    body with other plane faces (a primitive's flat, a box with a hole and
    fillets) and a prism among its leaves (walls along one axis between
    two caps, S9e.4a's construction of those faces, its arcs tangent to its
    lines inside it): K1's part. (2) *Primitives.* The curved faces on one
    stored surface (equal stored data) are one primitive's. A sphere's is
    the whole ball on its stored frame; a cylinder's or a cone's spans its
    faces' axial range on its frame, past each end by a quarter of it (at
    least of its radius) but at a cap (a plane face normal to its axis at
    that end whose outward normal points away from the range where its
    material lies inside its quadric, into it where outside: a blind
    hole's floor), its height there rounded once, a cone clamped at its
    apex (S9e.4b.3c.3a's). Coaxial primitives share one axis: a cylinder or
    a cone whose stored axis lies on an earlier one's in the chain's order
    within the resolution (parallel within `1e-9`, its stored origin within
    the resolution of that axis; a sphere's axis its stored frame's normal
    through its centre) is built on that primitive's frame (`Frame3::at`
    where it moves along it), its heights read in it, so `bug28773`'s discs
    at `x = 1` are its caps at one height and G9's pin's rod lies on its
    ball's axis exactly (its stored origin 8.9e-16 off it, where the rod's
    meeting with the ball would be no circle). Every plane face must be a
    cap of some primitive (else S9e.4b.4b.2's). (3) *The chain.* The
    primitives in the order of their reach from their axes or centres (a
    sphere's radius, a cylinder's, a cone's widest end over its faces'
    range), widest first, ties in face order; the first's material inside
    its quadric; each next one cut where its material lies outside its
    quadric, else in common where its faces meet the earlier ones' curved
    faces along convex edges only (S9e.4b.3c.3b's bends: the dome's cone and
    sphere, the pin's rod and sphere), fused where along concave edges only
    or none (joined through the caps: E5's thin cylinder through the wide
    one's ring). E5 is one fuse, G9's dome and pin one common each,
    `bug28773`'s comb five Booleans (the radius-1 cylinder fused with the
    0.03 tube's, less the 0.025 bore, fused with the 0.02, less the 0.01,
    fused with the rod). Each Boolean is S9e.4b.3c.3b's tree's (an inner
    one's assembly its given model, the stored vertices splitting rims in
    every one), the last one's given model matched to the stored topology
    by S9e.2's match, its ids the stored ones. Numbers: the engine's own,
    no new field or degree. A draft of the step's code (not committed)
    imports all five bodies but K1's and needs three changes of the
    arrangement, each the existing rule carried to a meeting of two curved
    faces: (a) a stored vertex within the resolution of a cylinder's or a
    cone's meeting with a sphere splits it (S9e.4b.3c.1's `split_at`, until
    now a sphere's circles only), at the carrier's rational unit direction
    nearest the vertex's (its half-angle tangent rounded once), a point of
    one quadratic field: the pin's rim, which OCCT splits at its seams'
    vertices; (b) such a meeting matched to a stored circle (OCCT's
    coaxial rim) whose frame turns against the carrier's angle is flipped,
    as conics and a sphere's circles are (S9e.2's `flip`); (c) an edge of a
    given meeting of two curved faces against a face of the partner on one
    surface with one of the edge's faces lies on that face's surface and
    is taken with the faces on one surface (as `Along`), not met again by
    S9e.3b's three surfaces, whose fibres there are not separated (the
    pin's rim on the dome's sphere: `ComputationLimit("a given meeting's
    points not separated by a projection (S9e.3b)")` without it). (4) *Its
    other queries.* Classification by the chain's set functions over each
    primitive's location within the resolution (a common the lesser, a
    fuse the greater, a cut reversed); bounds every primitive's between the
    body's axial ends widened by the stored edges'; mass by the stored
    topology; a rigid motion reads the chain off the moved stored topology
    again; a result given to another Boolean is S9e.3a's. (5) *Degenerate
    and refused.* S9's rules unchanged in every Boolean of the chain. A
    body's own faces tangent along an edge (a bend within `1e-6` of flat:
    a hemisphere's end on a rod of its radius, a fillet on its faces) are
    `Degenerate("an imported body of several primitives whose faces are
    tangent along an edge")`, named for the body (the S9e text holds stored
    faces tangent along an edge degenerate), after the plane faces' check,
    so K1's part is S9e.4b.4b.2's. Refused: plane faces other than ends,
    `OutOfDomain("an imported body of several primitives with plane faces
    other than their ends (S9e.4b.4b.2)")`; faces of one surface whose
    material lies on different sides, a widest primitive whose material
    lies outside it, a primitive meeting the earlier ones along edges of
    both kinds and a chain whose model does not match the stored topology,
    `OutOfDomain("an imported body of several primitives other than a
    Boolean chain of them (S9e.4b.4c)")`; a torus among the faces stays
    S9e.4b's. The general refusal becomes `OutOfDomain("an imported solid
    other than a prism, a sphere, a cone, a torus or a Boolean of its
    primitives (S9e.4b)")`. Expected in the DRAW trial: E5, G9 and `bug417`
    evaluating (the draft's: G9's area 2,814.0138, DRAW's printed
    2814.01); `bug28773` refused by S9's rules as `Degenerate("two
    cylinders' axes within rounding of parallel")` (`_1`'s prism on its
    leaning disc frame against `_2`'s primitives on its walls' axis); K1
    as S9e.4b.4b.2's. (6) *Fuzzing.* The `boolean` target's `IMPORTED`
    stage imports the chained stage's first results of no or one curved
    face; it takes those of several sphere, cylinder and cone faces too
    (the kernel's writer writes circles alone, so coaxial meetings and plane
    sections), cut by the turned box again, the volume the chained cut's;
    the replay measures how many. (7) *Evidence first.* Bodies OCCT writes
    in one Boolean of two primitives (a `write` block), every section a
    circle or a line (coaxial meetings, planes normal to the axes): `shaft`
    (E5's: a cylinder fused with a coaxial narrower one from inside it to a
    cap, on a turned frame), `cup` (a cylinder less a coaxial bore ending
    inside it, its floor a cap), `dome` (G9's `c`: a frustum common a ball
    on its axis, meeting it at a rational circle), `pin` (G9's `d`: a rod
    from the dome's base plane common that ball), `bead` (a ball less a
    coaxial rod through it: no plane face), `knob` (a ball fused with a
    coaxial rod from inside it to a cap); declared `degenerate` `capsule`
    (a rod fused with a ball of its radius about its end's centre: tangent
    along the rim) and declared `unsupported` `rounded` (K1's part, a prism
    of a rounded square less a cylinder across it: S9e.4b.4b.2's); cases of
    each against boxes, slabs, balls and rods, the body the tool, the dome
    and the pin both imported (G9's configuration: one sphere, the pin's
    disc on the dome's base plane) and a chain; the reference S9e.3a's
    chained reference on the constructions OCCT was given with S9e.4b.3c.3a's
    checks and each body's closed form
    (`generate_primitive_chains_boolean_fixtures.py --check`, a CI group
    `primitive-chains`, `test_primitive_chains_boolean_reference.py`); a
    native capture keyed on the refusal the step removes
    (`compare_primitive_chains_boolean.py`, `imported.rs`'s general
    refusal's present text); then the kernel and its tests
    (`tests/primitive_chain_booleans.rs`: every fixture within `1e-9`
    enclosures, the degenerate refused, histories over the stored ids,
    deterministic, moved rigidly; the kernel's own stepped shaft and bead
    written, read back and imported) and the DRAW trial of the five.
  * **S9e.4b.4b.1 evidence (2026-10-06).** Eight bodies OCCT wrote
    (`boolean-primitive-chains-bodies.txt`, `write` blocks of one
    primitive's rows, a `boolean` row and another's, written by
    `compare_primitive_chains_boolean.py --write-bodies` to
    `rust/fixtures/imported/chain_*.brep`), each one Boolean of two coaxial
    primitives, every section a circle: `shaft` (E5's: a cylinder of radius
    3 over `[0, 4]` on `SKEW` fused with a coaxial one of radius 3/2 over
    `[2, 9]`), `cup` (a cylinder of radius 3 over `[0, 5]` on `SKEW2` less a
    coaxial bore of radius 2 over `[1, 7]`), `dome` (G9's `c`: a frustum of
    radii 7/2 and 5/2 over 4 on the world's axes common a ball of radius 5
    about `(0, 0, -2)`, meeting at radius 3, 2 above the base), `pin` (G9's
    `d`: a rod of radius 1.4 over `[0, 6]` common that ball, its disc on the
    dome's base plane), `bead` (a ball of radius 5 on `SKEW4` less a coaxial
    rod of radius 3 through it: no plane face), `knob` (a ball of radius 2
    about `(1, 1, 1)` fused with a coaxial rod of radius 1.2 from its centre
    to a cap 5 above it), declared `degenerate` `capsule` (a rod of radius 2
    fused with a ball of its radius about its top's centre: tangent along
    the rim) and declared `unsupported` `rounded` (K1's part: a prism of a
    square of side 8 its corners rounded by arcs of radius 2, less a
    cylinder of radius 3/2 along `y` through it).
    `generate_primitive_chains_boolean_fixtures.py --check`: 24 cases of 8
    groups (21 solid, 3 `unsupported`; 18 of class `pieces`, 3 `both`, 3
    `chain`): `cup_ball` (a ball centred on the cup's rim between its wall
    and its bore), `shaft_box` (a `TILT` box across the step), `dome_pin`
    (G9's configuration, both imported: one sphere, the pin's disc inside the
    dome's on one plane), `chain_dome` (the dome less a `TILT` box biting its
    side, then with a `TILT` slab), `dome_slab`, `bead_slab` (`TILT` slabs
    through them, their cuts two solids), `knob_box` (a `TILT` box less the
    knob: the body the tool); declared `unsupported` `rounded_rod`. The
    reference is the constructions OCCT was given through S9e.3a's chained
    reference with S9e.4b.3c.3a's checks, relative to the case's size: the
    two families within 1.8e-35, each solid's closed form 3.6e-40, each
    body's closed form 4.3e-42 (the shaft's and the cup's sections in their
    chart times the stored axes' determinant; the dome, the pin and the knob
    on the world's axes), the pair identities 1.4e-41 and the area identity
    2.6e-40, Monte Carlo 2.7 standard errors (100,000 points a group),
    quadrature estimates 1.3e-32, every meeting's sine at least 0.29 and
    events at least 1.0e-4 of their range apart outside the declared group
    (the bead's events but its spacing, below), solid counts by rays at two
    resolutions with S9e.4b.3c.1's join; the declared group's checks kept
    apart within 1.1e-35. Every file read independently (`stored_records`):
    its curved faces' kinds its primitives', every stored vertex within
    2.3e-14 of the size on the construction's surfaces.
    `test_primitive_chains_boolean_reference.py` checks the dome's, the
    pin's, the knob's, the shaft's and the cup's closed forms against
    quadratures of their sections, the chained reference on the cup and the
    dome against theirs, every body's two primitives on one axis exactly in
    the rows (every meeting a circle: the dome's at radius 3, the bead's at
    `+-4`), the pin's disc on the dome's base plane, the capsule's tangency
    and the rounded square's arcs tangent to its lines, the case list and
    its protocol rows and the files apart from every earlier set's. The
    generator's check is a CI group of its own (`primitive-chains`, 10.9
    minutes on four workers locally under Python 3.9, 7.2 under 3.12);
    Python 3.9 and 3.12 write the same files. Corrections from the evidence,
    amending the refined decisions' plan (7): (a) the capsule is given to
    no case: its own tangency along its rim kept the reference's sweeps over
    twenty minutes of one worker for one group, which CI's job cannot
    afford; the kernel's tests refuse it on import (its file is checked
    with the others); (b) the dome's frustum first chosen, of radii 4 and 2,
    put the open end of its primitive in a draft of the step's code (a
    quarter of its face's range past it) on the plane tangent to the ball at
    its pole, `Degenerate` in the dome's own common; and a primitive in
    common or cut needs to reach past the other primitives anyway (the
    dome's cone past the ball's top, a bore past the ball it goes through):
    the first primitive and those in common or cut reach past the other
    primitives' bounds along their axes as well (amending (2)), and the
    set's frustum is of radii 7/2 and 5/2; (c) a stored circle of a
    coaxial meeting is also a ring (the dome's rim): such a ring turns with
    the meeting's carrier (`ring_about`, as a cone's section), amending
    (3)'s (b); (d) the `TILT` slabs through the dome and the bead cut each
    in two, their cuts two solids; the chain is the dome's (a shaft less a
    rod parallel to its axis, then with a slab, kept one worker over thirty
    minutes); (e) the bead's bore, a cylinder on `SKEW4`'s stored axes,
    meets its ball in its rims, where the reference finds each rim as both
    surfaces' event a rounding apart, so its spacing is not checked (its
    sines and gaps are). The capture
    `occt-boolean-primitive-chains-preimplementation`
    (`compare_primitive_chains_boolean.py`, keyed on the refusal the step
    removes: `imported.rs`'s general `OutOfDomain("an imported solid other
    than a prism, a sphere, a cone or a torus (S9e.4b)")`, whose text the
    step changes; the kernel's probe `unsupported` on all 24, every body
    refused on import by it): every result valid, 20 matching (within 2e-8
    of the reference), 4 reviewed, in
    `occt-boolean-primitive-chains-divergences.json`: the cup's three with a
    ball across its rim, whose faces OCCT bounds by its approximated quartics
    (its edges' tolerance up to 8.0e-6, 10 of its 29 edges B-splines),
    volumes up to 1.4e-4 relative by BRepGProp's default integration and,
    measured adaptively at 1e-10 and 1e-12 by a diagnostic build, the fuse
    and the cut within 6.9e-8 and 8.8e-8 (no closer than those curves'
    tolerance; the native cup alone measures its 29 pi to 1e-15), the
    common within 4.0e-9; and the knob's common with a box, its volume 2.2e-8
    off by the default integration, 2.3e-11 adaptively. S9e.4b.4b.1's kernel
    next.
  * **S9e.4b.4b.1 implemented** (`solid/imported.rs`'s `primitives_piece`:
    the primitives of a body of several sphere, cylinder and cone faces,
    their Booleans and frames, `Piece`'s other primitives and
    `Tree::Primitive(i)`, the classification by the chain; `curved/pieces.rs`:
    every primitive in the tree's Booleans, each at a seam of its own, a
    body's sphere among its primitives for one-sphere pairs, its own tangency
    and mismatch named for it; `curved/graph.rs`, `curved/given.rs`,
    `curved/assemble.rs` and `curved/procedural.rs`: the three changes of the
    arrangement), as the refined decisions describe with the evidence's
    amendments: a body of several sphere, cylinder and cone faces whose
    plane faces are all ends of their primitives is a Boolean chain of
    those primitives, widest first, each next one cut where its material
    lies outside its quadric, in common along convex edges with the earlier
    ones' curved faces, fused otherwise; the first primitive and those in
    common or cut reach past the other primitives' bounds along their axes,
    those fused a quarter of their faces' range past their open ends;
    coaxial ones on one frame. All 24 fixtures as declared (21 within the
    kernel's enclosures, each at most `1e-9` wide; the rounded box's 3
    `OutOfDomain("an imported body of several primitives with plane faces
    other than their ends (S9e.4b.4b.2)")`), every history complete over
    the imported bodies' stored ids, results deterministic and moved
    rigidly, both inputs translated and turned keeping the reference's
    volumes, every body its chain (its closed-form volume, points in and off
    it classified, its stored vertices on its boundary), the capsule refused
    on import as `Degenerate("an imported body of several primitives whose
    faces are tangent along an edge")`, and the kernel's own stepped shaft,
    cup and comb (a disc fused with a tube from below, its bore a third
    primitive) written by its writer, read back and imported, each Boolean
    with a turned rod its own result's volume
    (`tests/primitive_chain_booleans.rs`, 7.9 s in release on a host at load
    10 to 18, 6.6 s at `opt-level` 2 with debug assertions, 6.5 s in release
    with debug assertions under the emulated correctly rounded `hypot`,
    where the pieces' and the joints' earlier test files pass too).
    `compare_primitive_chains_boolean.py` 14 matches and 10 reviewed (the 4
    captured and, with the kernel, entity counts: the kernel's meetings of
    the cup's two cylinders with the ball split at their turning points into
    S9d.2's graphs where OCCT's are approximated edges, OCCT's unified
    results keeping its closed surfaces' seams), every enclosure within the
    reference; every comparison of `HANDOFF.md`'s table unchanged with 0
    failures (and `compare_split.py` 72/56, `compare_brep.py --family
    spline` 10/3, `compare_brep_io.py` 6,835/7, `compare_step.py` 23/6 on
    STEP-b's SDK); the release suite (646 tests) and the tools' unit tests
    (355) passing, the ledger unchanged. Amendments, from the
    implementation: (a) the kernel's writer writes circles alone, and a
    Boolean on a given result of coaxial cylinders stores its rings as
    ellipses of equal axes, so the kernel's comb is two Booleans (a disc
    fused with a tube prism), not `bug28773`'s five; (b) `split_at`'s
    meetings are placed at the rational unit direction nearest the stored
    vertex's by the half-angle tangent rounded once (`unit_near`), a binary64
    division correctly rounded on every host. The `boolean` fuzz target's
    `IMPORTED` stage takes the chained stage's first results of several
    sphere, cylinder and cone faces too: replayed with debug assertions, 80
    such first results of the corpus's 1,468 inputs and regressions reach
    it, 26 imported and cut within the chained cut's volume, 51 not written
    by the kernel's writer (their meetings no circles), 3 refused as
    S9e.4b.4b.2's (a plane face other than the primitives' ends). Replays
    with debug assertions, one process an input, natively and under the
    emulated correctly rounded `hypot`: the boolean corpus and its
    regressions (1,468 inputs) and the split corpus and its regressions
    (3,571), no failure, the slowest 9.1 s and 2.0 s (6.4 s and 1.5 s
    emulated) on a host at load 10 to 18. A trial of the targeted DRAW
    restore cases on the Rust adapter and natively (not the survey: nothing
    registered): `bfuse_complex/E5` (its area 1,831,656.94, native DRAW's
    printed 1.83166e6), `bcut_complex/G9` and `bugs/modalg_2/bug417` (their
    area 2,814.0138, native's 2814.01) evaluate on both backends with every
    check; `bugs/modalg_6/bug28773` is `Degenerate("two cylinders' axes
    within rounding of parallel")` (its tube `_1` S9e.4a's prism on its
    bottom disc's stored frame, leaning 2.2e-33 off `_2`'s primitives on the
    walls' axis), `bfuse_complex/K1` S9e.4b.4b.2's (plane faces other than
    its primitives' ends); `bugs/modalg_2/bug476_1` to `_8` (`OCC485a`, its
    tori among its faces) keep the general refusal, under its new text.
    Open: S9e.4b.4b.2 (K1's rounded box less a cylinder: hulls with several
    primitives, prism leaves); `bug28773`'s tube would need S9e.4a's prisms
    on their walls' axes. Pending: its capture's Linux record. DRAW survey:
    below (E5, G9 and `bug417` registered; nothing else moving but the
    general refusal's text). Campaigns at `d1869f2f` (600 s each, sampled replays, every
    switch on): `boolean` clean, 1,025 runs, the slowest input 19 s under
    AddressSanitizer (its startup replay 1,601 s; at `64673ebd` it overran
    the hour's budget at load 13 to 16, no input failing); `split` clean,
    1,972 runs, none slow.
  * **DRAW survey of S9e.4b.4b.1 (2026-10-06, `UPSTREAM_TESTS.md`).** At
    `c62470a3` (`s9c2-kernel` with S9e.4b.4b.1 and its campaigns' record;
    the public dataset, 120 seconds a case, four at once, on a host at load
    3 to 15). The 1,802 self-contained cases of the Boolean group on both
    backends: every status, every refusal's reason and every error the last
    survey's (`16121052`) field for field, 987 evaluating and registered,
    592 refused, 223 unsupported on both, none failing or timing out. The
    1,814 cases restoring a shape for a Boolean on the Rust adapter, and
    the 171 the import reaches on native DRAW too: native DRAW's every
    status and reason the last survey's; on the Rust adapter only the 36
    bodies refused by `imported.rs`'s general refusal move, as S9e.4b.4b.1's
    trial found: `bfuse_complex/E5` (the stepped shaft fused with a box on
    its step's plane), `bcut_complex/G9` and `bugs/modalg_2/bug417` (the
    dome less the pin) evaluate on both backends with every check;
    `bugs/modalg_6/bug28773` is `Degenerate("two cylinders' axes within
    rounding of parallel")`, `bfuse_complex/K1` `OutOfDomain("an imported
    body of several primitives with plane faces other than their ends
    (S9e.4b.4b.2)")`; the other 31 keep the general refusal under its new
    text, `OutOfDomain("an imported solid other than a prism, a sphere, a
    cone, a torus or a Boolean of its primitives (S9e.4b)")`:
    `bugs/modalg_2/bug476_1` to `_8` (`OCC485a`, tori among its faces) and
    23 cases of the `_2d` grids (`bcommon_2d/M6` to `N2`, `bcut_2d/N1` to
    `N5`, `bopcommon_2d/M6` to `N2`, `boptuc_2d/M3` to `M8`: a solid of
    several curved faces, `case_8_solid_repaired`, against wires, which the
    trial did not name). Of the 171: 59 evaluate on both backends, 75 are
    refused by S9's rules, 31 are bodies none of the kernel's constructions,
    1 is S9e.4b.4b.2's, 3 spline bodies (S9f) and 2 arguments of several
    solids; 24 are unsupported natively too. Four cases time out at their
    first restores, as in the last three surveys
    (`bugs/modalg_1/buc60531_1`, `_2`, `bugs/modalg_5/bug23849_1`, `_3`);
    the reader's slowest besides end within the 120 s at this load:
    `buc60532_2` and `bugs/modalg_6/bug23585`, which timed out at the last
    survey's, in 97 and 94 s, `buc60532` and `buc60532_1` in 114 s,
    unsupported as before (the reader takes no
    `SurfaceOfLinearExtrusion`; `tolerance`). No case fails, crashes or
    panics on either backend (the 276 restore cases `failed` on the Rust
    adapter are the reader's validator rejecting their files, as before).
    The volume audit (`vprops` and `sprops` before each `checkprops`, both
    backends, the 1,043 audited before and the three): both backends'
    values the last audit's bit for bit on all 1,043 (no moved value, the
    largest relative move 0), the statuses the same, the same 35
    disagreements. E5's volume 97,940,759.17, area 1,831,656.94 and centre
    agree with native DRAW's printed digits (9.79408e7, 1.83166e6); G9's
    and `bug417`'s volume 8,941.5371, area 2,814.0138 and centre with
    native's (8941.54, 2814.01). They are registered (`data`,
    `viewer_skipped` on both; 1,132 cases), and the ledger does not change
    (`--ledger` holds). A full contract run of the manifest holds on both
    backends with the dataset (30 seconds a case), E5, G9 and `bug417` in
    0.7, 1.4 and 1.3 s on the Rust adapter, the slowest Boolean case 3.9 s
    (`bopfuse_simple/ZK8`, on a host at load 4 to 15; `bopcut_simple/ZK8`
    5.8 in the last survey), the restore cases 0.1 to 3.3 s
    (`bcommon_complex/B9`), the rollex 0.4. No kernel change.
  * **S9e.4b.4b.2 refined, before its code (2026-10-06).** Why it is
    refused today: S9e.4b.4b.1's chain takes the curved faces on each stored
    surface as one primitive and refuses a body whose plane faces are not
    all ends of those primitives, before its own tangency check
    (`OutOfDomain("an imported body of several primitives with plane faces
    other than their ends (S9e.4b.4b.2)")`, `imported.rs`'s
    `primitives_piece`). `bfuse_complex/K1`'s part (`CTO904_cts20370-part`,
    read off its file by the converter; the dataset is not committed) has 11
    faces on 11 stored surfaces: the box `[0, 200]^3`'s six planes, its four
    edges along `z` rounded by quarter cylinders of radius 50 (each its own
    surface, about `(50, 50)`, `(50, 150)`, `(150, 50)` and `(150, 150)`,
    their 8 edges with the side faces stored `G1`, tangent), the faces `z =
    0` and `z = 200` each one loop of four lines and four arcs, less a
    cylinder of radius 40 along `-y` about `(100, *, 100)` through the faces
    `y = 0` and `y = 200` (a hole in each). As a chain its fillets are four
    primitives whose side faces are no primitive's ends and whose tangent
    edges are the body's own tangency; it is S9e.4a's prism of a rounded
    square (its fillets the profile's arcs, tangent to its lines inside it:
    S9e.4b.1's joints) less the bore. Its tool (`cts20370-tool`) is a rod of
    radius 40 along `x` about `(*, 100, 100)` from `x = -10` to its cap on
    the face `x = 200`: two cylinders of equal radii whose axes cross, tangent
    at `(100, 100, 60)` and `(100, 100, 140)`, so the case stays refused by
    S9's rules (a draft of this step's code imports the part, its volume
    6,565,486.678 its closed form, and refuses the fuse as `Degenerate("a
    tangency between the inputs (S9c)")`). S9e.4b.4b.1's trial and the
    surveys list no other restore case refused so; the `boolean` fuzz
    target's `IMPORTED` stage refuses 3 first results so (S9e.4b.4b.1's
    replay). Decisions. (1) *Sub-steps*, by what each unlocks:
    **S9e.4b.4b.2a** (this step): a body of several sphere, cylinder and
    cone faces one group of whose faces is a *prism leaf* with both its caps
    (2), the chain led by it (3): K1's part, a rounded or stadium plate with a
    bore across it, a boss, a dimple, a dome or a conical pocket;
    **S9e.4b.4b.2b**: other plane faces (a primitive's flat or a hull of
    planes among several primitives, S9e.4b.3c.3b's groups with several
    primitives), a prism with one cap (a prism boss on a primitive) or with
    both caps cut at their rims by another primitive (its profile read off
    its walls), and two prism leaves. (2) *The prism leaf*, read off the
    stored topology: two plane faces facing apart (their outward normals
    parallel within `1e-9`; the pairs in face order, each face as the bottom
    in turn), the bottom's loops each of whose edges' other faces is a *wall*
    along the bottom's inward normal `n` (a plane whose outward normal is
    perpendicular to `n` within `1e-9`, or a cylinder whose stored axis is
    parallel to it) sharing an edge with the top cap's faces, its outer loop
    among them (else no leaf: a prism on another primitive, or cut at its
    rim, is S9e.4b.4b.2b's); its faces the caps, those walls and every face on
    one of their surfaces (equal stored surfaces, or planes within the
    resolution of one plane facing one way: a wall or cap another primitive
    splits, its stored frames rounded apart). Its solid is S9e.4a's prism on
    those loops of the bottom alone (the bottom's stored frame, turned over to
    `n` exactly by `Frame3::flipped`; its profile their stored vertices and
    arcs rounded once into it, the arcs taken onto their circles by S9e.4b.1
    and S9e.4b.4a; its height the top cap's), the bottom's other loops (a
    primitive's crossing of the cap) not its. Its walls' edges with each
    other may be tangent (S9e.4b.1's joints of a line and an arc, as the
    construction's own). (3) *The chain.* The leaf first (its material
    inside), then the other surfaces' primitives as S9e.4b.4b.1's: widest
    first, each cut where its material lies outside its quadric, in common
    where its faces meet the earlier ones (the leaf's among them) along convex
    edges only, fused otherwise; a cylinder or a cone whose stored axis is
    parallel to the leaf's within `1e-9` is built on the leaf's frame's axes
    bit for bit at its own stored origin (`Frame3::at`: a boss or a pocket
    along the plate's axis, no axes within rounding of parallel between the
    chain's own Booleans), coaxial primitives on one frame as S9e.4b.4b.1's;
    the cut and common primitives reach past the leaf's bounds as past the
    other primitives'. Every plane face must be the leaf's or a primitive's
    end. Each Boolean is S9e.4b.3c.3b's tree's (an inner one's assembly its
    given model, the stored vertices splitting rims in every one), the last
    one's given model matched to the stored topology by S9e.2's match, its
    ids the stored ones. A leaf is tried only where S9e.4b.4b.1's chain of
    the primitives alone refuses for plane faces or for the body's own
    tangency (every body it imports keeps its chain), each leaf in turn, the
    first whose model matches the stored topology the body (`Piece`'s
    choice, the leaf read off the moved stored topology again by a rigid
    motion). Numbers: the engine's own and S9e.4a's prism's, no new field or
    degree. (4) *Its other queries* are S9e.4b.4b.1's: classification by the
    chain's set functions over each leaf's location (the prism's within the
    resolution), bounds every leaf's widened by the stored edges', mass by
    the stored topology, a result given to another Boolean S9e.3a's. (5)
    *Degenerate and refused.* S9's rules unchanged in every Boolean of the
    chain. A body's own faces tangent along an edge other than between two of
    the leaf's faces (a fillet on a primitive's faces, a primitive tangent to
    a cap) are `Degenerate("an imported body of several primitives whose
    faces are tangent along an edge")`, named for the body as S9e.4b.4b.1's
    (also where no leaf takes a body's tangent faces). Refused: plane faces
    other than the leaf's and the primitives' ends, `OutOfDomain("an
    imported body of several primitives with plane faces other than their
    ends or a prism's (S9e.4b.4b.2b)")`, replacing S9e.4b.4b.1's text; a leaf
    whose chain's model does not match the stored topology as S9e.4b.4b.1's
    (`OutOfDomain(... other than a Boolean chain of them (S9e.4b.4c))`). Of
    the leaves tried the first refusal their Booleans give is reported, else
    the body's own tangency, else S9e.4b.4b.2b's. (6) *`bug28773` and the
    walls' axis.* S9e.4a's prisms stay on their bottom cap's stored frame.
    A draft building them on their first wall cylinder's stored axes where
    those are not the cap frame's normal bit for bit (`_1`'s tube, its disc
    frame leaning 2.2e-33 off its walls' `x`; `_2`'s primitives on that axis
    since S9e.4b.4b.1) takes the case past its axes within rounding of
    parallel to `Degenerate("two faces within the resolution of one plane")`:
    `_1`'s top cap at `x = 1` and `_2`'s discs there lie 1e-17 apart, each
    height rounded once in its own frame, a flush contact S9's rules refuse.
    The change would unlock nothing and stays undone. (7) *Fuzzing.* The
    `IMPORTED` stage unchanged: its replay measures how many of its 3
    refusals import now. (8) *Evidence first.* Bodies OCCT writes in one
    Boolean of a prism of lines and arcs and a primitive (a `write` block),
    every section a circle or a line: `bore` (K1's part, a square its
    corners rounded on a turned frame less a cylinder across it through two
    walls), `boss` (a stadium plate fused with a cylinder along its axis from
    inside it to a cap), `dimple` (a rounded plate less a ball centred above
    its top cap), `pocket` (a stadium plate less a frustum along its axis
    narrowing to a floor inside it: a conical pocket), `dome` (a rounded
    plate fused with a ball about a point inside it, meeting its top cap
    only); declared `degenerate` `groove` (a rounded plate less a rod across
    it tangent to its bottom cap from inside: the body's own tangency, given
    to no case) and declared `unsupported` `notch` (a rounded plate less a
    blind rod across its rim through both caps: S9e.4b.4b.2b's); cases of
    each against boxes, slabs, balls and rods, the body the tool, two
    imported and a chain, and declared `degenerate` K1's configuration (the
    bore and a rod of its radius whose axis crosses its axis: two cylinders
    tangent at two points); S9e.4b.4b.1's `rounded` (K1's on the world's
    axes) and its `rounded_rod` cases declared solid with the kernel; the
    reference S9e.3a's chained reference on the constructions OCCT was given
    with S9e.4b.3c.3a's checks and each body's closed form
    (`generate_prism_leaves_boolean_fixtures.py --check`, a CI group
    `prism-leaves`, `test_prism_leaves_boolean_reference.py`); a native
    capture keyed on the refusal the step removes
    (`compare_prism_leaves_boolean.py`, `imported.rs`'s present text "with
    plane faces other than their ends (S9e.4b.4b.2)"); then the kernel and
    its tests (`tests/prism_leaf_booleans.rs`: every fixture within `1e-9`
    enclosures, the degenerate refused, histories over the stored ids,
    deterministic, moved rigidly; the kernel's own rounded plate with a bore
    and with a boss written, read back and imported) and the DRAW trial of
    K1 (expected refused as above).
  * **S9e.4b.4b.2a evidence (2026-10-06).** Seven bodies OCCT wrote
    (`boolean-prism-leaves-bodies.txt`, `write` blocks of a prism's rows, a
    `boolean` row and a primitive's, written by
    `compare_prism_leaves_boolean.py --write-bodies` to
    `rust/fixtures/imported/leaf_*.brep`), each one Boolean of a prism of
    lines and arcs (a square of side 9 its corners rounded by arcs of radius
    2, tangent to its lines; or a stadium of half discs of radius 5/2 about
    `(0, 0)` and `(6, 0)`) and a primitive, every section a circle or a line
    (lines and circles alone in every file): `bore` (K1's part: the rounded
    square over `[0, 6]` on the world's axes less a cylinder of radius 3/2
    along `y` about `(9/2, *, 3)` through its walls `y = 0` and `y = 9`),
    `boss` (the stadium over `[0, 2]` on `SKEW2` fused with a cylinder of
    radius 5/4 along its axis from inside it to a cap 3 above it), `dimple`
    (the rounded square over `[0, 3]` on the world's axes less a ball of
    radius 5/2 meeting its top cap at radius 2), `pocket` (the stadium over
    `[0, 3]` on `SKEW` across the boss's plate, less a frustum along its axis
    from radius 1/2 at height 1: a conical pocket with a flat floor), `dome`
    (the rounded square over `[0, 3]` on `SKEW2` fused with a ball of radius
    2 about a point 1/2 below its top cap, meeting it alone), declared
    `degenerate` `post` (the rounded square fused with a cylinder of radius 2
    about its corner's arc's centre from inside it to above it: the fillet's
    face and the post's on one surface along their arc) and declared
    `unsupported` `notch` (the rounded square over `[0, 2]` less a blind rod
    across its rim through both caps: no cap's loop the prism's, S9e.4b.4b.2b's).
    `generate_prism_leaves_boolean_fixtures.py --check`: 24 cases of 8 groups
    (21 solid, 3 declared `degenerate`; 18 of class `pieces`, 3 `both`, 3
    `chain`): `bore_rod` (a `TILT` rod across the bore's plate), `boss_box`
    (a `TILT` box across the boss and the plate's rim), `dimple_ball` (a ball
    across the dimple's rim), `pocket_rod` (a `TILT` rod through the pocket's
    floor), `dome_box` (a `TILT` box less the dome: the body the tool),
    `boss_pocket` (both imported, their plates crossing), `chain_dimple` (the
    dimple less a `TILT` rod through its rounded corner, then with a `TILT`
    slab: its cut two solids) and declared `degenerate` `bore_touch` (a rod
    of radius 1 along `x` tangent to the bore from above at one point:
    `Degenerate("a tangency between the inputs (S9c)")`). The reference is
    the constructions OCCT was given through S9e.3a's chained reference with
    S9e.4b.3c.3a's checks, relative to the case's size: the two families
    within 3.8e-36, each solid's closed form 1.5e-36, each body's closed
    form 4.0e-42 (its prism's section in its chart times the stored axes'
    determinant, less or plus its primitive's part: the bore's part between
    the walls by its own chart, the dome's ball beyond the cap at its world
    distance from the centre, the pocket's frustum from its base's height in
    the plate's chart), the pair identities 6.8e-42 and the area identity
    1.0e-40, Monte Carlo 2.5 standard errors (100,000 points a group),
    quadrature estimates 1.1e-32, every meeting's sine at least 0.11 and
    events at least 4.5e-5 of their range apart outside the declared group
    (the chain's spacing below), solid counts by rays at two resolutions;
    the declared group's checks kept apart within 6.4e-36. Every file read
    independently (`stored_records`): its curved faces' kinds its prism's
    arcs' and its primitive's, every stored vertex within 6.7e-16 of the size
    on the construction's surfaces. `test_prism_leaves_boolean_reference.py`
    checks every body's closed form against quadratures of its sections, the
    chained reference on the dimple and the pocket against theirs, every
    profile's arcs tangent to the lines they join exactly, every primitive's
    place (the bore within the flats of the walls it crosses and clear of the
    caps, the boss and the pocket on their plates' axes, the balls meeting the
    top caps alone within the flats), the post on its fillet's circle, the
    notch through both caps, the touching rod's axis the radii's sum above
    the bore's, the case list and its protocol rows and the files apart from
    every earlier set's. The generator's check is a CI group of its own
    (`prism-leaves`, 12 minutes on four workers locally under Python 3.9,
    8.3 under 3.12); Python 3.9 and 3.12 write the same files. Corrections
    from the evidence, amending the refined decisions' plan (8): (a) K1's
    part on a turned frame is degenerate by S9's rules: the bore's axis lies
    within rounding of the walls `x = 0` and `x = 9` it runs along
    (`Degenerate("a plane within rounding of a cylinder's direction")` on
    the kernel's own plates on turned frames), so the bore is on the world's
    axes, as K1's; (b) a rounded square whose side is four
    times its arcs' radius (K1's own proportions, S9e.4b.4b.1's `rounded`) on
    a turned frame made its Booleans with a cross rod and with a ball
    `Degenerate("two cylinders' section within the resolution of a node")`
    (the kernel's own plates too: S9's rule on the arcs' whole cylinders), so
    the set's squares are of side 9; (c) the first degenerate body, a groove
    tangent to the bottom cap from inside, OCCT wrote with a non-manifold
    edge the converter rejects: the post replaces it; (d) K1's own
    configuration (a rod of the bore's radius whose axis crosses it, tangent
    at two points) the kernel does not refuse alike in every Boolean (a
    kernel prism with a round hole against that rod: the fuse
    `InvalidTopology("non_manifold_vertex")`, the cut evaluating, the common
    `Degenerate("solids touching at a vertex")`, before this step's code),
    so the declared case is a rod tangent at one point, K1's own left to the
    DRAW trial and the engine's behaviour to a task of its own; (e) the
    reference's cost grows with a case's faces: a first set with the bore and
    the boss both imported, the bore's chain with a box and a slab and the
    notch against a box ran over twenty minutes of one worker a group, so
    the pair is the boss and the pocket, the chain the dimple's with a rod and
    the notch and the post are given to no case (the kernel's tests refuse
    them on import); (f) the chain's slab crosses the dimple's rim circle,
    where its top cap and its ball meet the slab's face at one point, which
    the reference finds as several events of the slab's chords within 1e-21
    of each other, so that group's spacing is not checked (its sines and gaps
    are). The capture `occt-boolean-prism-leaves-preimplementation`
    (`compare_prism_leaves_boolean.py`, keyed on the refusal the step
    removes: `imported.rs`'s `OutOfDomain("an imported body of several
    primitives with plane faces other than their ends (S9e.4b.4b.2)")`, whose
    text the step changes; the kernel's probe `unsupported` on all 24, every
    body refused on import by it): every result valid with the reference's
    solid counts, 19 matching (within 2e-8 of the reference), 5 reviewed, in
    `occt-boolean-prism-leaves-divergences.json`: the bore's three with a
    `TILT` rod (its cylinder met off its axis in a quartic; volumes up to
    1.7e-5 relative by BRepGProp's default integration), the boss and the
    pocket's common and the dimple's chain's cut (2.8e-8 and 3.3e-8), each
    within 2.4e-9 measured adaptively at 1e-10 and 1e-12 by a diagnostic
    build (its edges' tolerance at most 1.5e-7). S9e.4b.4b.2a's kernel next.
  * **S9e.4b.4b.2a implemented** (`solid/imported.rs`'s `leaves`,
    `chain_piece` and `joint`, `prism_on` read off chosen loops of a cap,
    `Piece`'s choice the leaf's; `curved/pieces.rs` unchanged but for its
    text: the leaf is `Tree::Primitive(0)`, S9e.4a's prism's model), as the
    refined decisions describe with the evidence's amendments: a body of
    several sphere, cylinder and cone faces whose plane faces are not all
    its primitives' ends, or whose faces are tangent along an edge, is
    S9e.4b.4b.1's chain led by a prism leaf, each leaf tried in turn and the
    first whose model matches the stored topology the body. All 24 fixtures
    as declared (21 within the kernel's enclosures, each at most `1e-9`
    wide; the rod tangent to the bore refused as `Degenerate("a tangency
    between the inputs (S9c)")`), every history complete over the imported
    bodies' stored ids, results deterministic and moved rigidly, both inputs
    translated and turned keeping the reference's volumes, every body its
    chain (its closed-form volume, points in and off it classified, its
    stored vertices on its boundary), the post refused on import as
    `Degenerate("an imported body of several primitives whose faces are
    tangent along an edge")` and the notch as `OutOfDomain("an imported body
    of several primitives with plane faces other than their ends or a
    prism's (S9e.4b.4b.2b)")`, and the kernel's own boss, dimple and dome
    plates (a rounded square of side 9 with a rod along its axis, less and
    fused with a ball) written by its writer, read back and imported, each
    Boolean with a turned rod its own result's volume
    (`tests/prism_leaf_booleans.rs`, 12.8 s in release on a host at load 10
    to 13, 18.0 s at `opt-level` 2 with debug assertions, 13.0 s in release
    with debug assertions under the emulated correctly rounded `hypot`, where
    the primitive chains', piece trees' and piece forms' test files pass
    too). S9e.4b.4b.1's `rounded` (K1's on the world's axes), declared
    `unsupported` until now, is solid within its reference
    (`generate_primitive_chains_boolean_fixtures.py` declares it so, its
    margins checked: sines at least 0.29, events 1.0e-4 apart; its unit test
    and `tests/primitive_chain_booleans.rs` too, the body its closed form).
    `compare_prism_leaves_boolean.py` 15 matches and 9 reviewed (the 5
    captured and, with the kernel, entity counts: the kernel's meeting of the
    bore with the crossing rod split at its turning points into S9d.2's
    graphs where OCCT's is an approximated edge, OCCT's unified results
    keeping the boss's, the rod's and the frustum's seams and a vertex on a
    closed edge), every enclosure within the reference;
    `compare_primitive_chains_boolean.py` 14 and 10, the kernel within the
    reference on all 24; every other comparison of `HANDOFF.md`'s table
    unchanged with 0 failures; the release suite (652 tests) and the tools'
    unit tests (361) passing, the ledger unchanged. Amendments, from
    the implementation: (a) a prism's walls' joints are exempt from the
    body's own tangency only as such (two walls along one axis tangent along
    a line along it, `joint`), checked before the plane faces in the chain of
    the primitives alone too, so a primitive on a fillet's surface (the post)
    is the body's own tangency and a body with a prism no leaf takes (the
    notch) S9e.4b.4b.2b's; (b) a leaf's walls must reach its top cap, so a
    prism on another primitive (a rounded boss on a disc, whose disc's bottom
    would be taken for the boss's top) is S9e.4b.4b.2b's rather than a chain
    the match rejects; (c) S9e.4b.4b.1's refusal of plane faces other than
    the primitives' ends is gone, the chain of the primitives alone refusing
    with S9e.4b.4b.2b's text before the leaves are tried. The `boolean` fuzz
    target's `IMPORTED` stage takes such first results already (S9e.4b.4b.1's
    stage takes every first result of sphere, cylinder and cone faces and
    planes): replayed with debug assertions, the same 80 first results of
    several curved faces of the corpus's 1,468 inputs and regressions reach
    it, 26 imported and cut within the chained cut's volume, 51 not written by
    the kernel's writer, the 3 refused as S9e.4b.4b.2's refused as
    S9e.4b.4b.2b's (two of prisms of lines and arcs fused at two heights,
    their caps on two levels: two leaves; one of a cone and a ball with a
    plane face no primitive's end: a hull among several primitives); under
    the emulated
    `hypot` 89 reach it (32 imported, 54 not written, the same 3 refused).
    Replays with debug assertions, one process an input, natively and under
    the emulated correctly rounded `hypot`: the boolean corpus and its
    regressions (1,468 inputs) and the split corpus and its regressions
    (3,571), no failure, the slowest 6.0 s and 1.2 s (5.5 s and 1.2 s
    emulated) on a host at load 6 to 12. A trial of the targeted DRAW restore
    cases on the Rust adapter and natively (not the survey: nothing
    registered): `bfuse_complex/K1`'s part (`CTO904_cts20370-part`) imports
    as its chain led by its rounded box, and the fuse is
    `Degenerate("a tangency between the inputs (S9c)")` (its tool a rod of
    the bore's radius whose axis crosses it, tangent at two points; native
    DRAW evaluates it); `bugs/modalg_6/bug28773` stays `Degenerate("two
    cylinders' axes within rounding of parallel")` (refined decision (6));
    `bfuse_complex/E5`, `bcut_complex/G9` and `bugs/modalg_2/bug417`
    evaluate on both backends with every check as before. Open:
    S9e.4b.4b.2b (other plane faces: a primitive's flat, a hull of planes
    among several primitives, a prism on a primitive or cut at both caps'
    rims, two leaves: the fuzz target's 3); the kernel's handling of two
    cylinders of equal radii whose axes cross at right angles (K1's own
    configuration, tangent at two points: on kernel-built prisms the fuse
    `InvalidTopology("non_manifold_vertex")`, the cut evaluating, the common
    `Degenerate("solids touching at a vertex")`, the imported K1 refused as
    a tangency), a task of its own (fixed since: a result touching itself at
    a vertex, `9c06cdf3`). Pending: its capture's Linux record. DRAW
    survey: that of S9e.4b.4b.2a, S9e.4b.4b.2b.1 and the cylinder pairs'
    fixes, below (K1 refused as a tangency between the inputs, as its trial
    found; nothing registered). Campaigns at `3aba844c`, with the fix of a hole and a rod of its
    radius crossing it (600 s each, sampled replays, every switch on):
    `boolean` clean, 892 runs, the slowest input 31 s under
    AddressSanitizer, its startup replay 2,200 s; `split` clean, 1,571 runs,
    none slow.
  * **S9e.4b.4b.2b refined, before its code (2026-10-06).** Why it is
    refused today: S9e.4b.4b.2a's chain, of the primitives alone or led by a
    prism leaf, takes every plane face as a primitive's end or a leaf's face
    and refuses any other (`OutOfDomain("an imported body of several
    primitives with plane faces other than their ends or a prism's
    (S9e.4b.4b.2b)")`, `imported.rs`'s `chain_piece`). The `boolean` fuzz
    target's `IMPORTED` stage refuses 3 first results so (S9e.4b.4b.2a's
    replay, again at `298fcf3a`: of the corpus's 1,470 inputs and
    regressions none other), read off the kernel's writer's files: (a) a
    ball of radius 2.8125 about `(0.375, -0.125, 1.125)` below its equator's
    plane `z = 1.125` fused with a frustum along `z` about the origin
    narrowing from radius 0.9375 at `z = 0` to its disc at `z = 2.25`, its
    base hidden in the ball: the plane face an annulus between the ball's
    rim and the frustum's circle, convex along the rim, concave along the
    circle (a primitive's flat; four faces, no leaf); (b) a regular pentagon
    of circumradius 1 about the origin over `[0, 1/4]` fused with a stadium
    (half discs of radius 9/4 about `(-3/4, -5/4)` and `(1/2, -5/4)`) over
    `[1/4, 3/4]`, flush: the stadium a leaf (its bottom's outer loop's walls
    reaching its top), the pentagon's walls and bottom other plane faces, its
    top hidden in the stadium (a prism of one cap); (c) the square `[-4,
    4]^2` with a hole of radius 2 about the origin over `[0, 3/2]` less the
    same shape moved by `(9/8, 9/8)` over `[3/8, 9/8]`: the plate a leaf
    (with its hole), the moved copy's floor, roof and two walls other plane
    faces joined along concave edges (a pocket), the moved hole's cylinder a
    crescent standing in it (its faces meeting the pocket's along concave
    edges, the plate's hole along two lines). The survey's restore cases
    hold none (S9e.4b.4b.2a's trial: K1 imports; the rest of the 171 none of
    these), so no DRAW case is expected to move. Decisions. (1) *Sub-steps*,
    by what each unlocks: **S9e.4b.4b.2b.1** (this step): the three (2):
    flats of one part, prisms of one cap, pockets with their teeth;
    **S9e.4b.4b.2b.2**: the rest of S9e.4b.4b.2b's list: a hull of planes
    across several parts (a plane face convex along two primitives' faces),
    a convex group of planes meeting the parts along concave edges alone (a
    polyhedral boss no prism of one cap takes), a prism with both caps cut at
    their rims (S9e.4b.4b.2a's notch, its profile read off its walls), a
    pocket meeting a part along edges of both kinds or holding a prism (a
    square post in a slot), a group of planes of both kinds, two leaves with
    both caps each. (2) *The parts*, read off the stored topology, tried only
    after S9e.4b.4b.2a's chains (of the primitives alone, then each leaf)
    refuse for other plane faces: for each of those chains in turn, the
    leaves' first and the primitives' alone last, (a) *prisms of one cap*: in
    face order, a plane face no cap of a leaf (S9e.4b.4b.2a's pairs) and no
    face of the chain's leaf or an earlier such prism, whose outer loop's
    every edge's other face is a wall along its inward normal (S9e.4b.4b.2a's
    walls, none of the leaf's): its faces it, those walls and every face on
    their surfaces; S9e.4a's prism on that loop alone, past its walls' far
    end (their range along the normal) by a quarter of it (its other cap
    hidden in the part it stands on, as a fused primitive's open end), on the
    face's stored frame or, where the walls' stored cylinders' axes point
    against the inward normal (a boss read off its top cap), on that frame
    turned over and extruded back from the cap, so its walls' frames and the
    stored ones point one way; a prism parallel to the chain's first prism
    within `1e-9` on that prism's axes bit for bit at its own stored origin
    (`Frame3::at`, as S9e.4b.4b.2a's primitives along the leaf); (b) the
    primitives as S9e.4b.4b.2a's (each surface's, every one's caps its ends);
    (c) the other plane faces (no end, no prism's), in components through
    their edges with each other: a component whose edges among its faces are
    all convex (or of one face) is a *flat* of the one part (a prism or a
    primitive with its ends) it meets along convex edges, that part common
    the convex hull of its planes (S9e.4b.3c.3b's hull leaf; the others it
    meets along concave edges, fused to it); a component whose edges among
    its faces are all concave is a *pocket*: the convex hull of its planes
    turned over (S9e.4b.3c.3b's pocket), less its *teeth*, the primitives
    it meets along concave edges alone (none of them in the chain), cut from
    the chain; it may meet the other parts along convex edges alone (its
    rims). A plane of a hull parallel to a prism's or a cylinder's or a
    cone's axis of the chain within `1e-9` lies on that axis's axes bit for
    bit at its stored origin (a pocket's floor normal to the plate's walls
    exactly). (3) *The tree.* The prisms first (the leaf, then the prisms of
    one cap), then the primitives but the teeth as S9e.4b.4b.2a's (widest
    first; each next one cut where its material lies outside its quadric, in
    common where its faces and its flats' meet the earlier parts' along
    convex edges only, fused otherwise), each part common its flats' hull,
    then each pocket (its hull less each tooth in turn) cut from the chain;
    each Boolean S9e.4b.3c.3b's tree's (an inner one's assembly its given
    model, the stored vertices splitting rims in every one), the last one's
    given model matched to the stored topology by S9e.2's match, its ids the
    stored ones; the first chain whose model matches the body (`Piece`'s
    choice, read off the moved stored topology again by a rigid motion). The
    fuzz target's three: (a) the ball common its flat's half-space, fused
    with the frustum; (b) the stadium fused with the pentagon past its top by
    a sixteenth; (c) the plate less the pocket less the crescent's cylinder.
    Numbers: the engine's own and S9e.4a's prism's, no new field or degree.
    (4) *Its other queries* are S9e.4b.4b.2a's: classification by the tree's
    set functions (a hull's sides within the resolution), bounds every part's
    widened by the stored edges', mass by the stored topology, a result given
    to another Boolean S9e.3a's. (5) *Degenerate and refused.* S9's rules
    unchanged in every Boolean of the tree. A body's own faces tangent along
    an edge are S9e.4b.4b.2a's `Degenerate`, a prism's walls' joints exempt
    within one prism, two faces of one stored surface (a face OCCT split at
    its surface's seam) no tangency; a chain with the other plane faces whose
    prisms part such a joint is one the body is not (the next is tried).
    Refused as S9e.4b.4b.2b.2's: a component of both kinds, a flat meeting no
    part or several along convex edges, a pocket meeting a part along both
    kinds or meeting a prism along concave edges, a tooth with a flat,
    `OutOfDomain("an imported body of several primitives with plane faces
    other than their ends, a prism's, a flat or a pocket (S9e.4b.4b.2b.2)")`,
    replacing S9e.4b.4b.2a's text (the notch now refused so); a chain whose
    model does not match as S9e.4b.4b.1's (`OutOfDomain(... other than a
    Boolean chain of them (S9e.4b.4c)")`); of the chains tried the first
    refusal their Booleans give is reported, else the body's own tangency,
    else S9e.4b.4b.2b.2's. (6) *Fuzzing.* The `IMPORTED` stage unchanged: its
    replay measures how many of its 3 refusals import now (a draft of the
    step's code, not committed, imports all three, each cut by the turned
    box within the chained cut's volume, and no other import moves). (7)
    *Evidence first.* Bodies OCCT writes in one Boolean (a `write` block),
    every section a circle or a line: `flat` ((a): a ball's half, a sphere's
    row from its south pole to its equator, fused with a frustum off its
    axis on a turned frame), `stack` ((b): a hexagon of rational corners
    under a stadium, flush, on a turned frame), `cake` (a stadium boss read
    off its top cap, from inside a box to its cap, on a turned frame: a prism
    of one cap with arcs), `slot` ((c): a plate with a hole less its copy's
    band on the world's axes, the copy turned half a turn so each hole's
    seam lies inside the other's disc and OCCT splits no face at a seam);
    cases of each against boxes, rods and balls, the body the tool, two
    imported and a chain, and a declared `degenerate` ball resting on the
    flat's frustum's top; S9e.4b.4b.2a's notch refused with the new text by
    the kernel's tests; the reference S9e.3a's chained reference on the
    constructions OCCT was given (a ball's half its ball common a cylinder
    below its equator, a plate with a hole its box less the hole's cylinder)
    with S9e.4b.3c.3a's checks and each body's closed form
    (`generate_plane_parts_boolean_fixtures.py --check`, a CI group
    `plane-parts`, `test_plane_parts_boolean_reference.py`); a native capture
    keyed on the refusal the step removes (`compare_plane_parts_boolean.py`,
    `imported.rs`'s present text ending "prism's (S9e.4b.4b.2b)"); then the
    kernel and its tests (`tests/plane_part_booleans.rs`: every fixture
    within `1e-9` enclosures, the degenerate refused, histories over the
    stored ids, deterministic, moved rigidly; the kernel's own cake and stack
    written, read back and imported, their Booleans the kernel's own
    results'; the fuzz target's three bodies rebuilt by the kernel and
    imported) and a DRAW trial of the restore cases the import reaches.
  * **S9e.4b.4b.2b.1 evidence (2026-10-06).** Four bodies OCCT wrote
    (`boolean-plane-parts-bodies.txt`, `write` blocks of one solid's rows, a
    `boolean` row and the other's, written by
    `compare_plane_parts_boolean.py --write-bodies` to
    `rust/fixtures/imported/part_*.brep`), each one Boolean, every section a
    circle or a line (lines and circles alone in every file): `flat` (the
    fuzz target's (a): a ball of radius 3 on `SKEW` below its equator's plane,
    a sphere's row from its south pole to its equator, fused with a frustum
    on its frame's axes about `(3/4, -1/2)` from radius 5/4 at height -1,
    inside the ball, to radius 3/4 at 5/2), `stack` ((b): a hexagon of
    rational corners over `[0, 1]` on `SKEW2` fused with a stadium of half
    discs of radius 5/2 about `(0, 0)` and `(6, 0)` over `[1, 3]`, flush, its
    footprint holding the hexagon's), `cake` (a box `9 x 7 x 2` on `SKEW4`
    fused with a stadium boss of half discs of radius 3/2 about `(5/2, 7/2)`
    and `(13/2, 7/2)` from height 3/2 inside it to its cap at 4: a prism of
    one cap with arcs, read off its top) and `slot` ((c): the square `[0,
    8]^2` with a hole of radius 2 about `(4, 4)` over `[0, 3]` on the world's
    axes less the square `[3/4, 35/4] x [5/4, 37/4]` with a hole of radius
    3/2 about `(11/2, 4)` over `[1, 2]` on those axes turned half a turn: the
    pocket's floor, roof and two walls, and a crescent of the second hole's
    disc standing in it, its cylinder meeting the first hole's in two lines).
    `generate_plane_parts_boolean_fixtures.py --check`: 24 cases of 8 groups
    (21 solid, 3 declared `degenerate`; 18 of class `pieces`, 3 `both`, 3
    `chain`): `slot_rod` (an upright rod through the crescent and the plate's
    hole), `chain_cake` (the cake less a `TILT` rod through the boss, then
    with a `TILT` box), `cake_stack` (both imported, the cake's box across the
    stadium), `flat_box` (a `TILT` box across the disc and the frustum),
    `stack_rod` (a `TILT` rod through both levels), `cake_ball` (a ball across
    the boss's rim), `rod_flat` (a `TILT` rod less the flat: the body the
    tool) and declared `degenerate` `flat_touch` (a ball resting on the
    frustum's top disc at its centre: `Degenerate("a plane crossing a sphere
    within the resolution of tangency (S9d.1)")`). The reference is the
    constructions OCCT was given through S9e.3a's chained reference with
    S9e.4b.3c.3a's checks (the ball's half its ball common a cylinder below
    its equator's plane, each plate with a hole its box less its hole's
    cylinder: the reference models convex profiles), relative to the case's
    size: the two families within 1.4e-31 (the cake's groups; the others'
    within 5.1e-36), each solid's closed form 2.4e-36, each body's closed form
    6.2e-33 (the cake's; the others' 2.0e-37: the flat's ball's half and its
    frustum above the disc in the frustum's chart, the stack's and the
    cake's sections times their axes' determinant, the cake's boss above the
    box by its frame's origin's height in the box's chart, the slot's plate
    less its band's rectangle less both discs' union), the pair identities
    1.4e-42 and the area identity 2.0e-41, Monte Carlo 2.2 standard errors
    (100,000 points a group), quadrature estimates 1.5e-32, every meeting's
    sine at least 0.17 and events at least 8.0e-6 of their range apart
    outside the declared group, solid counts by rays at two resolutions; the
    declared group's checks kept apart within 1.3e-37. Every file read
    independently (`stored_records`): its curved faces' kinds its
    constructions', every stored vertex within 4.7e-16 of the size on the
    constructions' surfaces. `test_plane_parts_boolean_reference.py` checks
    every body's closed form against quadratures of its sections (the slot's
    by rows, its chords parting where the circles cross), the chained
    reference on the stack and the cake against theirs, every part's place
    (the frustum's base inside the ball, the hexagon inside the stadium's
    footprint, the boss inside the box's from inside it, the crescent's
    circle crossing the plate's hole and each hole's seam inside the other's
    disc), the declared case, the case list and its protocol rows and the
    files apart from every earlier set's. The generator's check is a CI group
    of its own (`plane-parts`, 18.5 minutes on four workers locally under
    Python 3.9, 12.7 under 3.12); Python 3.9 and 3.12 write the same
    files. Corrections from the evidence, amending the refined decisions'
    plan (7): (a) the slot's copy is turned half a turn about `z`: with both
    holes' seams toward `+x` OCCT split the crescent's wall at its seam into
    two faces on one surface (a draft's import found them no tangency, then
    its model, of one face there, did not match) or, the copy moved, wrote the
    plate's hole's wall with a seam edge the converter rejects (`SeamEdge`),
    so each hole's seam lies inside the other's disc; (b) the reference's cost
    grows with a case's faces: a first set with the flat and the stack both
    imported, the stack's chain with a slab and the slot with a `TILT` box ran
    18 to 26 minutes of one worker a group, so the pair is the cake and the
    stack, the chain the cake's with a box and the slot's partner an upright
    rod; (c) a `TILT` box less the flat put two of the reference's events
    within 6e-18 of each other wherever it was placed, so the body is the
    tool of a rod; (d) the declared ball resting on the frustum's top is S9's
    plane tangent to a sphere, `Degenerate("a plane crossing a sphere within
    the resolution of tangency (S9d.1)")`. The capture
    `occt-boolean-plane-parts-preimplementation`
    (`compare_plane_parts_boolean.py`, keyed on the refusal the step removes:
    `imported.rs`'s `OutOfDomain("an imported body of several primitives with
    plane faces other than their ends or a prism's (S9e.4b.4b.2b)")`, whose
    text the step changes; the kernel's probe `unsupported` on all 24 before
    the code, every body refused on import by it but the cake, which
    S9e.4b.4b.2a refuses as `Degenerate("an imported body of several
    primitives whose faces are tangent along an edge")`, its boss's arcs read
    as primitives tangent to its line walls, so the comparison counts the
    cake's refusal `unsupported` before the code): every result valid, 19
    matching (within 2e-8 of the reference), 5 reviewed, in
    `occt-boolean-plane-parts-divergences.json`: the cake's common with the
    ball and the rod's three with the flat, volumes up to 2.8e-7 relative by
    BRepGProp's default integration, within 9.7e-9 measured adaptively at
    1e-10 and 1e-12 by a diagnostic build but the rod's common, 3.3e-8 (no
    closer than its approximated meetings: 10 of its 24 edges B-splines, their
    tolerance 1e-7); and the touching ball's fuse, which OCCT keeps as two
    solids touching at a point where the reference's rays count one, their
    totals the reference's within 4.8e-15. S9e.4b.4b.2b.1's kernel next.
  * **S9e.4b.4b.2b.1 implemented** (`solid/imported.rs`'s `chain_piece`
    given the body's leaves, `bosses`, `plane_along`, `prism_on`'s prism of
    one cap; `curved/pieces.rs` unchanged: the flats' and pockets' hulls are
    S9e.4b.3c.3b's `Tree::Hull` leaves), as the refined decisions describe:
    after S9e.4b.4b.2a's chains (the primitives alone, then each leaf), each
    chain again with the other plane faces, the leaves' first and the
    primitives' alone last (`primitives_choices`, `Piece`'s choice read off
    the moved stored topology again by a rigid motion): prisms of one cap
    fused after the leaf, flats making their part common their hull, pockets
    less their teeth cut from the chain, the first whose model matches the
    stored topology the body. All 24 fixtures as declared (21 within the
    kernel's enclosures, each at most `1e-9` wide; the ball resting on the
    frustum refused as `Degenerate("a plane crossing a sphere within the
    resolution of tangency (S9d.1)")`), every history complete over the
    imported bodies' stored ids, results deterministic and moved rigidly,
    both inputs translated and turned keeping the reference's volumes, every
    body its chain (its closed-form volume, points in and off it classified,
    its stored vertices on its boundary), S9e.4b.4b.2a's notch refused as
    `OutOfDomain("an imported body of several primitives with plane faces
    other than their ends, a prism's, a flat or a pocket (S9e.4b.4b.2b.2)")`,
    and the kernel's own cake and stack and the `boolean` fuzz target's three
    bodies (a ball's half fused with a frustum on its disc, a pentagon under a
    stadium, a plate with a hole less its moved copy's band) built by the
    kernel, written by its writer, read back and imported, their volumes, and
    their Booleans with a turned box the kernel's own results' but the
    stack's (`tests/plane_part_booleans.rs`, 13.1 s in release on a host at
    load 10 to 12, 18.3 s at `opt-level` 2 with debug assertions, 14.5 to
    18.4 s in release with debug assertions and overflow checks under the
    emulated correctly rounded `hypot`, where the prism leaves', primitive
    chains', piece trees' and piece forms' test files pass too).
    `compare_plane_parts_boolean.py` 15 matches and 9 reviewed (the 5
    captured and, with the kernel, entity counts: the faces agree, OCCT's
    unified results keeping seam edges of their closed surfaces or dividing
    a meeting into edges otherwise, the kernel's meetings of the ball with
    the boss's arc walls split at their turning points), every enclosure
    within the reference; every other comparison of `HANDOFF.md`'s table
    unchanged with 0 failures (and `compare_split.py` 72/56, `compare_brep.py
    --family spline` 10/3, `compare_brep_io.py` 6,835/7, `compare_step.py`
    23/6 on STEP-b's SDK); the release suite (661 tests) and the tools' unit
    tests (366) passing, the ledger unchanged. Amendments, from the
    implementation: (a) the cake, which S9e.4b.4b.2a refused as the body's
    own tangency (its boss's arcs read as primitives tangent to its line
    walls), now evaluates: a joint of one prism of one cap is the
    construction's own; (b) a hull plane's axis may be a sphere's frame's
    normal too (each primitive's frame is among the chain's axes, as
    S9e.4b.4b.1's coaxial test takes them); (c) the kernel's own stack (a
    hexagon prism fused with a stadium prism on one frame) given to a Boolean
    is refused as `Degenerate("an edge of one input on a face of the other")`
    with boxes, rods and balls alike, flush or with the hexagon reaching into
    the stadium, where its import evaluates them: the kernel's tests take its
    import alone, and the given model of such a fuse is left open; (d) OCCT's
    cake with a ball wholly inside it (a cavity) is `ComputationLimit("a
    cavity's containment the validator's rays leave undecided")` where the
    kernel's own cake evaluates, given to no case. The `boolean` fuzz
    target's `IMPORTED` stage (unchanged): replayed with debug assertions, the
    same 80 first results of several curved faces of the corpus's 1,470
    inputs and regressions reach it, 29 imported and cut within the chained
    cut's volume (the 3 refused as S9e.4b.4b.2b's among them), 51 not written
    by the kernel's writer, none refused; under the emulated `hypot` 89 (35
    imported, 54 not written). Replays with debug assertions, one process an
    input, natively and under the emulated correctly rounded `hypot`: the
    boolean corpus and its regressions (1,470 inputs) and the split corpus
    and its regressions (3,571), no failure, the slowest 5.8 s and 1.7 s (5.3
    s and 1.2 s emulated) on a host at load 6 to 12. A trial of the
    survey's 171 restore cases the import reaches, on the Rust adapter and
    natively (not the survey: nothing registered): no case moves with this
    step, 59 evaluating on both backends as before, `bfuse_complex/K1`
    refused as a tangency between the inputs as S9e.4b.4b.2a's trial found,
    `bugs/modalg_6/bug28773` as two cylinders' axes within rounding of
    parallel, none refused as S9e.4b.4b.2b's; native DRAW's statuses as
    before (147 evaluating, 24 unsupported). The trial found a draft's
    mistake, fixed before the commit: counting a primitive's caps among its
    faces for its bends with the earlier parts (S9e.4b.4b.1 counts its curved
    faces alone) changed `bug28773`'s comb's chain, which then did not match
    (`OutOfDomain(... other than a Boolean chain of them (S9e.4b.4c))`); a
    primitive's faces there are its curved faces and its flats. Open:
    S9e.4b.4b.2b.2 (hulls across several parts, polyhedral bosses, the
    notch's prism cut at both caps' rims, mixed pockets and components, two
    leaves of two caps each); the kernel's own stack given to a Boolean (c).
    Pending: its capture's Linux record. DRAW survey: below (no case moving
    with this step, as its trial found). Campaigns at `c3ce4d42` (600 s
    each, sampled replays, every switch on): `boolean` clean, 930 runs, the
    slowest input 24 s under AddressSanitizer; `split` clean, 1,519 runs,
    none slow.
  * **DRAW survey of S9e.4b.4b.2a, S9e.4b.4b.2b.1 and the cylinder pairs'
    fixes (2026-10-09, `UPSTREAM_TESTS.md`).** At `826346b7` (`s9c2-kernel`
    with S9e.4b.4b.2a, a result touching itself at a vertex refused
    (`9c06cdf3`), the near node's wider margin for crossing cylinders
    (`8dd43713`), S9e.4b.4b.2b.1 and their campaigns' records; the public
    dataset, 120 seconds a case, four at once, on a host at load 8 to 18).
    The two fixes are general rules of the curved engine; neither moves any
    case. The 1,802 self-contained cases of the Boolean group on both
    backends: every status, every refusal's reason and every error the last
    survey's (`c62470a3`) field for field, 987 evaluating and registered,
    592 refused, 223 unsupported on both, none failing or timing out; no
    case is refused as `Degenerate("a result touching itself at a
    vertex")`, and the near node (`two cylinders' section within the
    resolution of a node`) refuses the same 8 (`ZD9` and `ZE2` of the four
    `bop*_simple` grids). The 1,814 cases restoring a shape for a Boolean
    on the Rust adapter, and the 171 the import reaches on native DRAW too:
    native DRAW's every status and reason the last survey's; on the Rust
    adapter one case moves, `bfuse_complex/K1`, from `OutOfDomain("an
    imported body of several primitives with plane faces other than their
    ends (S9e.4b.4b.2)")` to `Degenerate("a tangency between the inputs
    (S9c)")`, as S9e.4b.4b.2a's trial found (its part imports as its
    rounded box less its bore; its tool a rod of the bore's radius whose
    axis crosses it, tangent at two points; native DRAW evaluates it). A
    replay at `1cd7bb3a` (S9e.4b.4b.2a merged, before either fix) gives
    K1 the same refusal, and the same statuses and reasons to the 28 cases
    refused as a near node or as a result or face touching itself
    (`bug29807_b1`'s near node and `bcut_complex/J3` among them): the
    rules are not the cause of any move, and S9e.4b.4b.2b.1 moves none, as
    its trial found. Of the 171: 59 evaluate on both backends, 76 are
    refused by S9's rules, 31 are bodies none of the kernel's constructions
    (the general refusal), 3 spline bodies (S9f) and 2 arguments of several
    solids, none S9e.4b.4b.2b's; 24 are unsupported natively too. Eight
    cases time out at their first restores: the four of the last four
    surveys (`bugs/modalg_1/buc60531_1`, `_2`, `bugs/modalg_5/bug23849_1`,
    `_3`) and, at this run's load, the reader's slowest besides
    (`buc60532`, `_1`, `_2` and `bugs/modalg_6/bug23585`), which end alone
    in 68 to 73 s (on a host at load 4 to 8), unsupported as before (the
    reader takes no `SurfaceOfLinearExtrusion`; `tolerance`). No case
    fails, crashes or panics on either backend (the 276 restore cases
    `failed` on the Rust adapter are the reader's validator rejecting their
    files, as before). The volume audit (`vprops` and `sprops` before each
    `checkprops`, both backends, the 1,046 audited before): both backends'
    values the last audit's bit for bit on all 1,046 (no moved value, the
    largest relative move 0), the statuses the same, the same 35
    disagreements. No case newly evaluates, so none is registered (1,132
    cases), and the ledger does not change (`--ledger` holds). A full
    contract run of the manifest holds on both backends with the dataset
    (30 seconds a case), the slowest Boolean case 9.8 s
    (`bopfuse_simple/ZK8`, on a host at load 7 to 21; 3.9 in the last
    survey), the restore cases 0.2 to 7.0 s (`bcommon_complex/B9`), E5, G9
    and `bug417` 0.9, 1.9 and 2.0 s, the rollex 0.7 to 1.8. No kernel
    change.
  * **CI's scheduled full replay at `c3ce4d42` (2026-10-07 to 10-09).**
    Two jobs of the scheduled "Rust geometry fuzzing" runs failed, no
    crash among them. (1) `Fuzz / split` exited 124 on all three runs
    (`37629390970`, `37785997266`, `37936642519`): its full replay of
    CI's 1,749 inputs had not finished its startup in the hour, the
    artifacts only slow units (`66b4f249`, `7bd1734f`, `a8a7bafb`,
    `27a0b6f9`, `cbf5efb7`). The full replay had been near the hour
    (1,662 inputs, 3,523 s of startup at `d1869f2f`, green), and
    S9e.4b.3b and S9e.4b.3c's `PIECE_BOOLEANS` stage made many inputs
    slower. Fix: `split` joins `REPLAY_SHARDS` with eight shards (the
    parallel track "The boolean target's full replay", below; the replay
    job's matrix now lists shards up to eight and excludes `boolean`'s and
    `degree_elevation`'s past their four). Replayed on the development
    Mac under the sanitizer with the shard jobs' limits, CI's 1,749
    inputs took 250 to 413 s of CPU a shard (2,729 s in all, every shard
    passing, the slowest input 41 s at load 55, within the 60-second
    limit) and the local corpus's first shard of 448 (of 3,551) 583 s: at
    CI's 2.6 times, 650 to 1,100 s a shard of CI's corpus and about
    1,500 s at the local corpus's size, plus a build, where one process
    would take about 7,100 s. (2) `Fuzz / step` exited 1 on the last two
    runs (job `113840241078` of `37936642519`, and `37785997266`):
    `timeout-f692f018…` over the 20-second limit in the full replay. The
    input is not a seed: the 10-07 campaign's mutation found it (a
    15-second slow unit there) and CI kept it in the corpus. It mutates
    the `bspline_trimmed` fixture: the trimmed face's corner `#70` moved
    off both its edges' curves (`z = 0` for `0.75`) and the trim pcurve's
    end `#85` a hair outside the surface (`u = 1.000000001`); the import
    is an invalid body (`vertex_off_curve` twice, `uncertified_uv_gap`,
    `enclosure_exceeds_resolution`). Not Linux's rounding and not a hang,
    a real cost: 1.5 to 2.2 s in release on the Mac and 65 s under the
    sanitizer, deterministic and finishing. Profiled, three quarters of it
    was `step::spline::locate`: locating a vertex inside an edge (the
    trimmed face's corners are inside its boundary curves and line
    pcurves) samples about 65 points and refines with 100 golden-section
    steps, and every point was the kernel's correctly rounded evaluation,
    BigRational de Boor on a degree-9 curve and on the bicubic surface
    through the pcurve; the unmutated fixture itself took 0.47 s for it.
    The search is a claim the validator certifies, so it needs no exact
    points: it now evaluates in binary64 de Boor on the homogeneous poles
    (periodic vectors extended as the kernel's, no libm, deterministic),
    the ends' tests within the tolerance still on correctly rounded
    points; the validator's checks are unchanged. The input takes 0.56 s
    in release and 4.4 s under the sanitizer (the validator's certified
    enclosures the rest; 11.5 s at load 25 with three sanitizer replays
    beside it), the companion slow unit `dd91a75e` (CI's 10 s) 0.55 and
    2.4 s, the fixture 0.04 s. Of 3,122 inputs (the fixtures, the
    regressions, CI's 2,506 and the local 743), 41 imports differ, each
    in located parameters alone (at most `1.4e-8` relative, where a vertex
    off its curve leaves the distance's minimum flat; the two valid
    bodies among them in the last bit), none in a result or an issue; the
    fixtures' imports are identical. The glibc emulation patch was not
    applied (not permitted in this session): the search's hot path calls
    no `hypot` and the cost reproduces on macOS. Regression:
    `tests/step.rs`, `a_vertex_off_its_spline_edges_is_located_and_refused`
    (the input's text, its issues and determinism), and in
    `step/spline.rs` the binary64 points against the kernel's (clamped,
    unclamped, periodic and rational curves, a rational surface) and the
    search taking two exact points; both inputs in
    `fuzz/regressions/step` with a README entry. Checks: fmt, clippy, the
    1.85 check, the release suite (664 tests), the tools' tests (366),
    `compare_step.py` with STEP-b's SDK (23 / 6, 0 failures), the step
    corpora and regressions replayed with debug assertions (3,122 inputs,
    no failure, the slowest 0.9 s).
  * **S9e.4b.4c refined, before its code (2026-10-06).** Why each class is
    refused today. (a) *An imported polyhedron against curved faces*:
    S9e.4b.2 decides an imported body of plane faces and line edges that is
    no prism on its stored vertices (S9b.2's stored model) in the polyhedral
    engine alone, and refuses it against a solid with curved faces or edges
    (`polyhedra/imported.rs`'s `involved`, `OutOfDomain("an imported
    polyhedron against curved faces (S9e.4b.4)")`) and a result of one given
    to a Boolean of curved faces (`curved/given.rs`, `OutOfDomain("an
    imported polyhedron's result given to a Boolean of curved faces
    (S9e.4b.4)")`): the curved engine's leaves are constructions. (b) *A
    polyhedron with a cavity*: the import takes one solid region of one
    shell (`imported.rs`, `OutOfDomain("an imported polyhedron with a cavity
    or several shells (S9e.4b.4)")`): S9e.4b.2's `hollow` and STEP-b's
    `box_void` (a `BREP_WITH_VOIDS`). (c) *Deeper trees*: S9e.4b.3c.3b's
    `tooth` (a pocket within a pocket, `OutOfDomain("an imported plane piece
    other than a Boolean tree of its primitive and its planes' hulls
    (S9e.4b.4)")`) and S9e.4b.4b's chains no model matches (`OutOfDomain(...
    other than a Boolean chain of them (S9e.4b.4c))`). Where they occur:
    none of the DRAW survey's 171 restore cases the import reaches is
    refused for any of them (the survey at `c62470a3` and the trials of
    S9e.4b.4b.2a and S9e.4b.4b.2b.1 since: 59 evaluating on both backends,
    75 refused by S9's rules, 31 by the general refusal, `K1` a tangency,
    `bug28773` axes within rounding of parallel, 3 spline bodies and 2
    arguments of several solids), so no DRAW case is expected to move; in
    the `boolean` fuzz target's replay at `c3ce4d42` (1,470 inputs and
    regressions, an instrumented build not committed) the `IMPORTED` stage
    gives 37 planar first results to a round or ball partner, one refused
    as (a), and imports no first result with a cavity; nothing reaches (c).
    Decisions. (1) *Sub-steps*, by what each unlocks and the machinery each
    takes: **S9e.4b.4c.1** (this step): (a) and (b) on one representation,
    an imported polyhedron's stored model of exactly planar triangles as a
    leaf of the curved engine (2), its cavities among its triangles (5); it
    unlocks every imported polyhedral part (a `.brep` or STEP body of planes
    no prism) drilled, bored or rounded by a construction with curved faces,
    and the fuzz target's polyhedra (6); **S9e.4b.4c.2**: (c), the trees'
    recognition (a pocket's own pockets, a component of both kinds), on the
    engine's own models, with no restore case and no fuzz input behind it.
    (2) *The polyhedron as a leaf of the curved engine*
    (`curved/meshes.rs`). S9e.4b.2's stored model (each stored face the
    polygon of its stored vertices, its ears clipped into triangles of its
    own vertices, each exactly planar) is the model: each triangle a model
    face on its exact plane (`Surf::Plane` through its first corner, its
    normal twice its vector area: a hull's kind of face), named by its stored
    face (its id, its stored plane, the sense taking that plane's normal to
    the outward one); each triangle's edge a model edge, the line between
    two stored vertices, named by the stored edge with those ends, else none
    (a diagonal inside a face); the stored vertices its vertices. A face
    whose clipping fails (S9b.2's trapezoids zipped, corners off the stored
    vertices) is `ComputationLimit("an imported polyhedron's face not cut
    into triangles of its own vertices")`, an edge of other than two
    triangles `InvalidTopology`. A point's membership is the parity of an
    exact ray's crossings with the triangles (S9b.2's `parity` on points of
    the arrangement's fields: each test a linear form of the point with
    rational coefficients, a crossing's side of an edge `(n . r) (w . x) -
    (w . r) (n . x)` over the point's offset `x`, `w` the edge's inward
    normal in the plane), retried along other rays where one meets an edge;
    a point on a triangle `On`. Pushed along directions (`p + e d1 + e^2 d2
    ...`), the triangles holding the point decide it by their wedges there
    (the triangle's plane through it bounded by the rays of its edges at it:
    none inside, one on an edge, two at a corner): the pushed point on a
    wedge `On`, else the parity of the crossings of a ray from it, the
    wedges' by the pushes (each test's first nonzero sign) and the other
    triangles' by the ray from the point itself. A point of a face's plane
    lies in it by its triangle's edges. Numbers: the stored binary64 points
    and their planes' normals (products of their differences); a triangle's
    meeting with a quadric is any plane's (S9c's conics, S9d's circles and
    sections), no new field or degree. (3) *Faces and names in the result.*
    The pieces of every triangle of one stored face are pieces of one input
    face (`assemble.rs` joins pieces of one input face's id kept the same
    way), so the result's face of a stored face is one face across its
    diagonals (each used both ways, dropped), on the stored plane; a section
    crossing a diagonal keeps a vertex there between its conics on the two
    triangles' planes (within rounding of each other). The history is over
    the stored ids directly, as S9e.4b.2's. (4) *Given results.* A result
    with an imported polyhedron given to a Boolean of curved faces re-runs
    its construction's arrangement, the polyhedron its stored triangles
    (S9e.1's direct slots for a curved result, S9e.2's match for S9b's
    polyhedral result); `given.rs`'s refusal goes. (5) *Cavities.* An
    imported body of plane faces and line edges of one solid region of
    several shells (an outer shell and its voids') is a polyhedron on its
    stored vertices, its triangles every shell's and its membership the
    parity over all of them, in both engines; several solid regions stay
    refused (`OutOfDomain("an imported polyhedron of several solids
    (S9e.4b)")`). The results: both assemblies put all of a solid's cavities
    into one shell and one void region, which a cavity split by the partner
    or two kept leave a disconnected shell (a draft of the step's code: the
    hollow box fused with a slab through its cavity, `InvalidTopology
    ("disconnected_shell")`, the kernel's own hollow box alike); each cavity
    becomes a shell of the solid and a void region of its own, its void's
    provenance as the one cavity's. In the curved assembly a shell of one
    input's faces alone is preset (the tool's in a cut a cavity, else a
    solid); where an input may hold a cavity (an imported polyhedron of
    several shells, a given solid with a cavity) every shell is tried by its
    orientation instead (built alone, its flux turned inward a cavity), as
    shells of both inputs' faces are already: the hollow box's kept cavity
    and a hollow tool's void turned into a solid (`InvalidTopology
    ("shell_orientation")` in the draft before, the kernel's own hollow box
    alike). S9b's containment of a cavity in an outer shell takes a ray's
    start on a fragment's plane outside the fragment as no crossing (a
    cavity's corner on the plane of the outer shell's face elsewhere left
    every ray undecided: `InvalidTopology("a cavity outside every shell")`).
    (6) *Fuzzing.* The `IMPORTED` stage's chained cut (the first result less
    the turned box, a polyhedron with the box's hole in it where its faces
    are planes) written, read back and imported, and the kernel's own cut,
    each less a ball about the box's axis crossing the hole's walls, by the
    chained byte's next bit (`MESHES`): where both evaluate, one volume. (7)
    *Degenerate and refused.* S9's rules unchanged against the triangles'
    planes (a partner's face within the resolution of a stored face's
    triangle and not on it, a tangency, a result touching itself); a face
    folded by its stored vertices as S9e.4b.2's. Staying refused: a cavity
    among several solids (S9c, by design: the hollow box fused with a ball
    inside its cavity); a polyhedron against spline walls (S9f's); a kept
    cavity whose containment the validator's rays leave undecided, beside a
    plane face holding a section's projected pcurve on every ray
    (`ComputationLimit`, the kernel's own cavity bodies alike: the hollow box
    fused with a rod through its cavity, its void a ring), the validator's
    track's. (8) *Evidence first.* Bodies OCCT wrote: S9e.4b.2's (their
    files read again: the pyramid, the frustum of a pyramid, the slanted
    wedge, the tetrahedron, the octahedron, the notched box and the hollow
    box), S9e.4a's ball and this step's `cavity` (a box in the `TURN30` frame
    less a box in the `TILT` frame inside it: each face's corners rounded,
    each face two triangles); cases against rods, balls and a frustum, the
    body the tool, two imported, two chains (S9b's polyhedral result and a
    curved one given), the cavities' ball and slab, declared `degenerate` a
    ball tangent to the tetrahedron's base and declared `unsupported` the
    ball inside the hollow box's cavity; the reference S9e.3a's chained
    reference with convex hulls of exact points
    (`polyhedra_curved_boolean_reference.py`; a body not convex a Boolean of
    hulls and boxes), each body's closed form, `generate_polyhedra_curved_
    boolean_fixtures.py --check` (a CI group `polyhedra-curved`),
    `test_polyhedra_curved_boolean_reference.py`; a native capture keyed on
    the refusal the step removes (`compare_polyhedra_curved_boolean.py`,
    `polyhedra/imported.rs`'s "an imported polyhedron against curved faces
    (S9e.4b.4)"); then the kernel and its tests
    (`tests/polyhedra_curved_booleans.rs`: every fixture within `1e-9`
    enclosures, the degenerate refused, histories over the stored ids,
    deterministic, moved rigidly, the kernel's own polyhedra written, read
    back and imported), S9e.4b.2's `hollow_slab` declared solid with it and
    STEP-b's `box_void` imported, and the DRAW trial of the restore cases
    the import reaches.
  * **S9e.4b.4c.1 evidence (2026-10-09).** Bodies OCCT wrote: S9e.4b.2's
    `pyramid`, `truncated`, `wedge`, `tetra`, `octa`, `notched` and
    `hollow` (their files read again), S9e.4a's `ball`, and this step's
    `cavity` (`boolean-polyhedra-curved-bodies.txt`, a `write` block of two
    `box` rows about a `boolean cut`, written by
    `compare_polyhedra_curved_boolean.py --write-bodies` to
    `rust/fixtures/imported/poly_cavity.brep`: a box `10 x 8 x 6` in the
    `TURN30` frame at `(1, 2, 0)` less a box `3 x 2.5 x 2` in the `TILT` frame
    at `(7/4, 51/8, 3)` inside it, its inner corners at least 1.35 inside the
    outer box in its chart). `generate_polyhedra_curved_boolean_fixtures.py
    --check`: 34 cases of 12 groups (30 solid, 3 `degenerate`, 1
    `unsupported`; 18 of class `curved`, 3 `both`, 6 `chain`, 7 `cavity`):
    `pyramid_rod` (an upright rod through the pyramid's tilted base and its
    lower face, crossing the base's diagonal), `wedge_ball` (a ball across
    the wedge's slanted top and its back face), `octa_cone` (an upright
    frustum through two of the octahedron's faces), `rod_truncated` (a `TILT`
    rod less the truncated pyramid: the body the tool, its cut two solids),
    `notched_ball` (a ball in the notched box's notch), `tetra_ball` (the
    tetrahedron and S9e.4a's imported ball), `chain_pyramid` (the pyramid
    less a box across its corner, S9b's polyhedral result, then with a ball),
    `chain_wedge` (the wedge less an upright rod, a curved result, then with
    a `TILT` box, their common two solids), `hollow_ball` (a ball through the
    hollow box's wall `x = 10` into its cavity across its wall `x = 7`: the
    fuse keeps the cavity, the cut opens it), `cavity_slab` (a `TILT` slab
    through the turned cavity parallel to its faces: the fuse's cavity in
    two, the cut two solids), declared `degenerate` `tetra_touch` (a ball of
    radius 3/2 resting on the tetrahedron's base at `(4, 9/2, -1/2)`, its
    corners' heights `-1/2` exactly: `Degenerate("a tangency between the
    inputs (S9c)")`) and declared `unsupported` `hollow_inner_fuse` (a ball
    inside the hollow box's cavity: `OutOfDomain("a cavity among several
    solids (S9c)")`). The reference (`polyhedra_curved_boolean_reference.py`)
    is S9e.3a's chained reference with a convex hull of exact points as an
    input (its faces' planes, their chords' two families, its exact closed
    form, its binary64 rays), each body the construction OCCT was given (a
    wedge's corners on its frame's stored axes, a polyhedron's binary64
    points, a Boolean of boxes its two prisms), with S9e.4b.3c.3a's checks
    relative to the case's size: the two families within 4.2e-36, each
    solid's closed form 8.7e-41, each body's closed form 1.0e-41 (the
    hulls' exact volumes, the hollow box's 936, the cavity's boxes times
    their axes' determinants), the pair identities 1.3e-41 and the area
    identity 8.9e-41, Monte Carlo 2.7 standard errors (100,000 points a
    group), quadrature estimates 1.0e-32, every meeting's sine at least 0.27
    and events at least 8.0e-6 of their range apart outside the declared
    groups (the notched ball's spacing below), solid counts by rays at two
    resolutions; the declared groups' checks kept apart within 4.7e-41.
    Every file read independently (`stored_records`): its faces planes (the
    ball's a sphere), every stored vertex within 2.2e-16 of the size on its
    construction's planes. `test_polyhedra_curved_boolean_reference.py`
    checks the hull input (a cube's and a skew tetrahedron's faces and exact
    volumes, their sweep against their closed form, their binary64 rays), each
    body the construction S9e.4b.2 or S9e.4a gave OCCT (the hulls' points
    and volumes S9e.4b.2's, the Booleans of boxes its rows), the cavity's
    inner box inside its outer box and its volume, the slab between the
    cavity's faces and past the outer box's corners, the partners' places
    (the hollow box's ball crossing its two walls within those faces and
    tangent to no plane of either box, the inner ball inside the cavity, the
    touching ball's foot on the tetrahedron's base inside it), the case list
    and its protocol rows and the new body's file apart from every earlier
    set's. The generator's check is a CI group of its own
    (`polyhedra-curved`, 9.7 minutes on four workers locally under Python
    3.9, 6.1 under 3.12); Python 3.9 and 3.12 write the same files.
    Corrections from the evidence, amending the refined decisions' plan (8),
    each found by the reference or by a draft of the step's code (not
    committed) run on the cases: (a) the hollow box's ball first lay tangent
    to the cavity's planes `z = 7` and `y = 3` outside their faces (the
    draft `Degenerate("a tangency between the inputs (S9c)")`), and the
    chain's `TILT` box was tangent to the wedge's bore, so both moved; (b) a
    ball resting on the octahedron's face `x + y + z = 18` (its radius the
    distance `sqrt(3) / 2` rounded down) misses the plane by 4.6e-17, which
    S9d.1's rule (a plane crossing a sphere within the resolution of
    tangency) does not refuse: the draft fused them into two solids; the
    declared case is a ball tangent to the tetrahedron's base exactly; (c) a
    rod through the cavity leaves a void ring in the fuse whose containment
    the validator's rays leave undecided (the kernel's own hollow box
    alike), so the cavities' partners are a ball through one wall and a slab;
    (d) a slab reaching far past the turned cavity's box made the reference's
    ray count join the cut's two solids (its spacing past the gap between
    them), so its sides lie just past the box's corners; (e) the notched
    box's corners where its notch's wall meets the box's edges are events
    of two of the box's planes on that wall at one height, found twice
    2.6e-16 apart, so `notched_ball`'s spacing is not checked (its sines and
    gaps are). Every solid case of the draft lay within the reference's
    volume. The capture `occt-boolean-polyhedra-curved-preimplementation`
    (`compare_polyhedra_curved_boolean.py`, keyed on the refusal the step
    removes: `polyhedra/imported.rs`'s `OutOfDomain("an imported polyhedron
    against curved faces (S9e.4b.4)")`; the kernel's probe `unsupported` on
    all 34 before the code, every case refused as S9e.4b.4's: against curved
    faces, a result given, or a cavity on import): every result valid with
    the reference's solid counts, 32 matching (within 2e-8 of the
    reference), 2 reviewed, in `occt-boolean-polyhedra-curved-
    divergences.json`: the wedge's common with its ball (volume 2.1e-8
    relative by BRepGProp's default integration, 3.2e-10 adaptively at 1e-10
    and 1e-12 by a diagnostic build) and the tetrahedron's common with the
    imported ball (4.4e-8, 4.0e-8 adaptively: no closer than its edges'
    tolerance, 1.3e-7, none a B-spline; its fuse and cut within the
    comparison's limit, 2.6e-9 and 8.0e-9 adaptively). S9e.4b.4c.1's kernel
    next.
  * **S9e.4b.4c.1 implemented** (`solid/boolean/curved/meshes.rs`: the
    polyhedron's stored triangles as a leaf model, `Mesh::member` and
    `in_face`; `polyhedra/imported.rs`'s `triangles` and `involved`;
    `curved/mod.rs`'s `applies` and `model_of`; `model.rs`'s `mesh` and
    `may_hold_cavities`; `given.rs` without its refusal and with the given
    solid's `cavities`; `assemble.rs` and `polyhedra.rs`: each cavity a shell
    and a void region of its own, every shell's orientation tried where an
    input holds a cavity, S9b's containment of a cavity's point on a far
    fragment's plane; `imported.rs`: one solid region of several shells a
    polyhedron), as the refined decisions describe. All 34 fixtures as
    declared (30 within the kernel's enclosures, each at most `1e-9` wide;
    the touching ball `Degenerate("a tangency between the inputs (S9c)")`,
    the hollow box's inner ball's fuse `OutOfDomain("a cavity among several
    solids (S9c)")`), every history complete over the imported bodies'
    stored ids, results deterministic and moved rigidly, both inputs
    translated and turned keeping the reference's volumes, every body its
    stored triangles (its volume its construction's, points in its material
    inside and in its cavity outside, its stored vertices on its boundary),
    the kernel's own notched and hollow boxes written by its writer, read
    back and imported, their Booleans with a leaning rod and a ball the
    kernel's own results' (15 of 16 evaluating, the hollow box's fuse with
    the rod through its cavity refused alike: its void a ring the
    validator's rays leave undecided), and the hollow box's cavity split in
    two by a slab and that result given to a ball, each cavity a shell of
    its own (`tests/polyhedra_curved_booleans.rs`, 6.8 s in release on a host
    at load 12 to 37, 7.5 s at `opt-level` 2 with debug assertions, 9.0 s in
    release with debug assertions and overflow checks under the emulated
    correctly rounded `hypot`, where the imported polyhedra's, imported
    bodies', plane parts', prism leaves', primitive chains' and piece trees'
    test files pass too). S9e.4b.2's `hollow_slab`, declared `unsupported`
    until now, is solid within its reference
    (`generate_imported_polyhedra_boolean_fixtures.py` declares it so, its
    margins checked; `tests/imported_polyhedra_booleans.rs`'s hollow box
    imports and a pyramid against a cylinder and a result of it given with
    a cylinder evaluate with the pair identities), and STEP-b's `box_void`
    imports, its cut and common by a box its volume
    (`tests/imported_booleans.rs`). `compare_polyhedra_curved_boolean.py` 18
    matches and 16 reviewed (the 2 captured and, with the kernel, entity
    counts: the kernel's meetings with a stored face of corners rounded in a
    turned frame cross its triangles' diagonals at vertices where OCCT's one
    edge crosses the face, OCCT's unified results keeping its sphere's and
    the bore's seams), every enclosure within the reference;
    `compare_imported_polyhedra_boolean.py` 47 and 1, the kernel within the
    reference on all 47 solid and empty cases; every other comparison of
    `HANDOFF.md`'s table unchanged with 0 failures; the release suite (668
    tests) and the tools' unit tests (371) passing, the ledger unchanged. The
    `boolean` fuzz target's `MESHES` stage (the chained cut's first solid of
    plane faces written, read back and imported, and the kernel's own cut,
    each less a ball about the turned box's axis, by the chained byte's next
    bit): replayed (an instrumented build, not committed), 15 of the
    corpus's 1,470 inputs and regressions reach it, 11 evaluating both ways
    with one volume, 2 refused alike as a tangency and 2 where the kernel's
    own polyhedral cut given is `ComputationLimit("a given result rebuilt
    differently")` (S9e.2's match) and the import a tangency; the `IMPORTED`
    stage's planar first result cut by a round partner that S9e.4b.2
    refused evaluates (37 such cuts, none refused). Replays with debug
    assertions, one process an input, natively and under the emulated
    correctly rounded `hypot` (with overflow checks): the boolean corpus and
    its regressions (1,470 inputs) and the split corpus and its regressions
    (3,571), no failure, the slowest 15.0 s and 4.5 s natively on a host at
    load 30 to 80 (15.0 s and 97 s emulated at that load: the emulation's
    `hypot` exact by big rationals, a split input 33 s emulated where 0.17 s
    natively at load 7). A trial of the survey's 171 restore cases the
    import reaches, on the Rust adapter and natively (not the survey: nothing
    registered): no case moves with this step, 59 evaluating on both
    backends as before, `bfuse_complex/K1` refused as a tangency between the
    inputs and `bugs/modalg_6/bug28773` as two cylinders' axes within
    rounding of parallel as S9e.4b.4b.2a's trial found, none reaching an
    imported polyhedron against curved faces or a cavity; native DRAW's
    statuses as before (147 evaluating, 24 unsupported). Open:
    S9e.4b.4c.2 (deeper trees: a pocket within a pocket, chains no model
    matches); a kept cavity's containment the validator's rays leave
    undecided beside a plane face holding a section's projected pcurve (the
    kernel's own cavity bodies alike), the validator's track's; a ball within
    the resolution of a plane but missing it (S9d.1's rule refuses only a
    crossing within it), found by the evidence's draft. Pending: its DRAW
    survey, its capture's Linux record and its campaigns.
  * **S9f.2b.2 refined, before its code (2026-10-03).** Why it is refused
    today: `spline_crossing::section` refuses a turning point of a spline
    wall's meeting with a crossing cylinder inside both faces
    (`OutOfDomain("a spline wall's meeting with a cylinder turning back
    inside the faces (S9f.2b.2)")`): S9f.2b.1's pieces are graphs over the
    run parameter `tau`, `w = (-B +- sqrt(D)) / A`, whose slope is
    unbounded at a root of `D`, so no piece of theirs reaches a turning
    point; and `conic_wall` refuses a cylinder's cap circle meeting a wall
    in a cap plane holding the wall's axis direction (`OutOfDomain("a
    cylinder's cap circle meeting a spline wall in a plane along the wall's
    axis (S9f.2b.2)")`): that plane meets the wall in generatrices at the
    roots of its trace's equation (degree `p`, `Q(alpha)`), the circle
    meets each where a quadratic over `Q(alpha)` vanishes, its points in
    the tower `Q(alpha)(sqrt(delta))`, `delta` in `Q(alpha)`, which the
    engine's numbers (`Qd`: `a + b sqrt(d)`, `d` rational) do not hold.
    Decisions.
    (1) *Representation: a second graph kind of the same curve.*
    `Crv::WallMeet` and `Curve3::WallMeet` are extended, no new curve: a
    piece about a turning point is a graph over the spline prism's height
    `w` (the wall's `v`), its run parameter the one root of `F(tau, w) = A
    w^2 + 2 B(tau) w + C(tau)` in a rational window `(t0, t1)` of the run
    inside one Bézier arc (one knot span of the wall), over a range of
    `w`. Over the height rather than over the cylinder's angle: at a
    rational height `F` is a polynomial of degree `2 p` in the arc's
    parameter (exact Sturm counts, no trigonometric chart), and the wall's
    `v` is already the curve's. In the engine `WallMeetCrv` gains `window:
    Option<[R; 2]>` (`None`: S9f.2b.1's graph over `tau` on its branch;
    `Some`: over `w`, the branch unused), its `range` then two heights,
    placed by the height (`Pos::T(w)`); a point lies on it when its profile
    point lies on the segment at a run parameter strictly inside the
    window, it lies on the cylinder exactly and its height lies in the
    range; its point at a rational height is the window's root, in a field
    `Q(beta)` of degree at most `2 p` (S9f.1's generators, kept by
    polynomial); its tangent `F_w S' - F_tau n` turned to rising `w`
    (`F_tau` never vanishes on it, (3)). In the topology `Curve3::WallMeet`
    gains `window: Option<[f64; 2]>` (the wall's `u`): with it `v = start +
    sweep f` and `u` the root of `a t^2 + 2 b t + c` (`t = v - v0`) in the
    window, by bisection in binary64; certified, `u` over a base by
    interval Newton inside the window and its Taylor coefficients by the
    implicit function theorem term by term (as `Toric`'s), on the window's
    span's exact Bernstein polynomials (`wall_meet::Span`), for every number
    of the integrands (`quadrature::Num`: enclosures, series, jets, the
    trait reading and setting a coefficient); on its own wall its pcurve is
    its own `(u, v)` (`Projection::own_wall`); its integrals are one piece
    on that span. Validity: a window inside one knot span of the wall's
    `u` domain whose ends give `a t^2 + 2 b t + c` opposite signs at sampled
    fractions, the rest as S9f.2b.1's. Rigid motion, the writer's and the
    importer's refusals, history, `curve_curve` and `curve_surface` as
    S9f.2b.1's (the window is in the wall's own `u`, which a motion keeps).
    (2) *Switches.* Each turning point `tau*` strictly inside both faces
    (S9f.2b.1's rules for one on a boundary, at a knot or repeated stand)
    gets one piece over the height: on the side of `tau*` where `D > 0` a
    rational `tau_s` where the steeper branch's slope `|dw / ds|` (`s` the
    profile's arc length, binary64) has fallen to two, at most a third of
    the way to the next root of `D` or the segment's end and inside the
    arc; the window from a rational `t0` beyond `tau_s` (where the slope is
    one) to a rational `t1` past `tau*` on the side `D < 0` (as far, inside
    the arc and short of the next root of `D`). The piece runs over
    `[w-(tau_s), w+(tau_s)]`, heights in `Q(sqrt(D(tau_s)))`, and S9f.2b.1's
    graphs over `tau` on both branches end at `tau_s` instead of `tau*`.
    The two switch points `(tau_s, w+-)` have a rational profile point and
    a surd height; inside both faces they are vertices of the arrangement
    (`CylPair::Mixed`'s switches, as S9d.2b's: a spline wall's meeting with
    a crossing cylinder is now found once per pair of faces in the pairs'
    pass, as a torus's), on a face's boundary a seam.
    (3) *Exact verification* of a piece over the height before it is kept
    (on failure the distances halved, up to twenty times, then
    `ComputationLimit`): (i) at `t0` the roots of `F` in `w` (surds, or
    none) lie outside the range, and at `t1` `D < 0` (exact signs of two
    rational surds, `tower_sign`); (ii) at a rational height inside the
    range (`-B(tau_s) / A`) exactly one root of `F` lies in the window
    (isolation, degree `2 p`); (iii) `H = A C'^2 - 4 B B' C' + 4 C B'^2`
    (degree `4 p - 2` in the arc's parameter) has no root in the closed
    window: a point with `F = F_tau = 0` has `H = 0` (`w = -C' / (2 B')`
    where `B' != 0`; `C' = 0` where `B' = 0`, and there `H = A C'^2`), so
    every root in the window is simple at every height. By (i) to (iii)
    the window holds one root at every height of the range: the piece is a
    graph, its point the window's root.
    (4) *Tower fields by a primitive element.* A cap circle of the cylinder
    prism in a plane holding the wall's axis direction (exactly
    perpendicular axes) meets the wall where its projection along the axis
    meets the profile: each arc's implicit equation `f(u, v) = 0` (degree
    `p`, S9f.2a's `Implicit`) at the circle's point in a chart of its
    half-angle tangent `t`, times `(1 + t^2)^p`, a polynomial of degree `2
    p` in `t`; each real root gives the point and its angle on the circle in
    one field `Q(beta)`, `beta = t` of degree at most `2 p`: the tower's
    primitive element instead of the tower, no new number. A root whose
    point is off the arc (another branch of the implicit curve, `locate`)
    is dropped, a point found on two arcs (at a knot) kept once; a chart
    whose antipode is a crossing turns by S9d.2c's rotations; a repeated
    root (the circle tangent to the generatrix: a turning point on the
    cap's rim) is `Degenerate`. The vertex is found again on the
    generatrix's line, on the meeting and on the circle by the engine's
    exact tests across fields (S9c.2b.2's enclosures).
    (5) *Degenerate*: S9f.2b.1's (7): the cylinder tangent to the wall; a
    turning point at an interior knot inside both faces' closures; a
    turning point on a face's boundary (a loop's on a cap's rim or a
    seam); a vertex's polynomial with a multiple root (a cap circle
    tangent to the wall); and the engine's rules (a vertex of one input on
    the other's face, crossings within the resolution). A piece over the
    height that does not verify is a `ComputationLimit`, never taken for a
    contact.
    (6) *Stays refused*: spline walls against spline walls on crossing axes
    (by design), spheres and cones (S9f.3), tori, given results with spline
    walls, rational and periodic profile splines. The fields' degrees stay
    at most `2 p <= 14` (`H`, of degree `4 p - 2 <= 26`, is only counted);
    no new limit.
    (7) *Evidence first*: S9f.2b's fixtures extended in
    `generate_spline_crossing_boolean_fixtures.py` (a pair is S9f.2b.2's
    when a turning point or a tower point lies inside its faces, the
    reference's `wall_cylinder_events` on a cap plane holding the axis
    counting the latter): towers alone (a rod ending inside the dome on its
    arch, a wide rod ending inside the lens across its two cubics), a loop
    cut by a tower's cap circle (the bulge), a loop cut by the spline
    prism's top cap (the bulge), a loop cut by a tilted rod's cap circle
    (the lens), and a declared degenerate pair whose loop's turning point
    lies on the rod's cap rim (`cap_turn`); the reference's checks as
    S9f.2b's, the same files under Python 3.9 and 3.12; a recapture of the
    whole set (`compare_spline_crossing_boolean.py` keyed on S9f.2b.2's two
    refusals in `spline_crossing.rs`: every S9f.2b.2 case `unsupported`
    while they stand, none after); then the kernel and its tests, and the
    `boolean` fuzz target's crossing variants reaching loops and towers
    through `SPLINE_CROSSING` (a switch of their own if they are slow under
    AddressSanitizer).
    * **S9f.2b.2 evidence (2026-10-03), before its code.**
      `generate_spline_crossing_boolean_fixtures.py --check` writes 51
      cases (16 fuses, 17 cuts, 18 commons; 45 solid, 6 degenerate; 28 of
      S9f.2b.1, 23 of S9f.2b.2), 17 of them new: `dome_cap_tower` (a rod
      along `x` ending at `x = 5/2` inside the dome, its cap circle on the
      arch's generatrix at `tau = 3/8`, no turning point), `lens_cap_tower`
      (a rod of radius 2.5 ending inside the lens, its cap circle on both
      cubics' generatrices: degree six), `bulge_cap_loop` (the bulge's loop
      cut by a tower cap circle at `x = 10.95`, 0.061 beyond the turning
      points), `bulge_side_cap_loop` (a loop cut by the bulge's top cap,
      its turning points 0.5 below it), `lens_tilt_cap_loop` (the tilted
      rod's loop on the lens cut by its cap circle, 0.30 beyond the turning
      points), and the declared degenerate `cap_turn` (a rod whose loop on
      the dome turns back on its cap's rim, `y = 15/8` at `x = 5/2`). A
      pair is declared S9f.2b.2's exactly when a turning point lies inside
      both faces or a tower point inside the wall's heights (`towers`).
      Checks as S9f.2b's: the divergence theorem within 2.6e-41, inclusion
      and exclusion 2.0e-41, the area identity 9.4e-41, faces' classes
      4.5e-41, the nine perpendicular pairs' commons as the product of
      chords 4.4e-42; margins outside the declared pairs at least 0.006
      (`bulge_cap_loop`'s cap circle crossing the bulge near its apex at a
      shallow angle), the declared pairs' zero; Python 3.9 and 3.12 the same
      files; the reference test checks the dome's tower points in closed
      form. The recapture (`compare_spline_crossing_boolean.py`, keyed on
      S9f.2b.2's two refusals in `spline_crossing.rs`): the 34 earlier rows
      to the bit, the 17 new done and valid with the reference's counts, 4
      matching, 13 reviewed (BRepGProp's default up to 2.0e-6 in volume;
      adaptive BRepGProp and Green's theorem over OCCT's faces, the better
      within 2.1e-9 in volume, 3.7e-9 in area, 1.8e-10 in the centre), five
      more solids' counts changing when unified; 17 matches and 34 reviewed
      in all, S9f.2b.2's 23 `unsupported` until its code, failures after.
      No correction to the decisions from the evidence. S9f.2b.2's kernel
      next.
  * **S9f.2b.2 implemented** (`solid/boolean/curved/spline_crossing.rs`'s
    `meeting`, `height_piece` and `tower_points`; `graph.rs`'s pairs pass,
    `meet.rs`'s sections, `assemble.rs`'s curves, `spline_parallel.rs`'s
    `Implicit::homogeneous`; the topology's `Curve3::WallMeet::window` with
    `topology/validate/wall_meet.rs`'s `eval_height`, `quadrature.rs`'s
    `HeightPiece` and its `Num` numbers' coefficients, `spline_flux.rs`'s
    `tensor_jets`, the validity rule): loops and towers as the refined
    decisions describe, S9f.2b.2's two `OutOfDomain` refusals gone. A spline
    wall's meeting with a crossing cylinder is found once per pair of faces
    (`CylPair::Mixed`: its pieces and switches); each turning point strictly
    inside both faces gets its graph over the height, verified exactly ((i)
    to (iii)) before it is kept, the graphs over the run on both branches
    ending at its switch; a cap circle in a plane holding the wall's axis
    meets each arc where the arc's implicit equation vanishes along the
    circle's half-angle tangent, its points in that one field. All 51
    fixtures as the reference (the 45 results within the kernel's enclosures
    of volume, area and centre, each at most `1e-9` wide; the 6 degenerate
    refused, `cap_turn` as a turning point on a face's boundary); every
    history complete, results deterministic and moved rigidly
    (`tests/spline_crossing_booleans.rs`, 8 tests, 43 s in release and 50 s
    at `opt-level` 2 with debug assertions: also the cap circle on the
    dome's wall along its axis evaluated, its tower vertex at `(2, 2,
    sqrt(3) - 1)`, and the loops' graphs over `v` each in one knot span,
    bracketing their roots; the module's 5: a graph over `v` on the wall and
    the rod with jets over points and ranges enclosing its points and
    slopes, a range's coefficients holding the point's and narrow, one piece
    on its window's span). `compare_spline_crossing_boolean.py` 14 matches
    and 37 reviewed, the kernel within the reference on all 45 results and
    none `unsupported` (the 15 loop cases' reviews name their counts: the
    kernel cuts each loop at its four switches and where the meeting crosses
    its rod's faces' boundaries in the arrangement, OCCT at its rod's seam
    touching the loop at a turning point and at its intersection edges' own
    splits); every other comparison unchanged, `compare_split.py` 72/56,
    `compare_brep.py --family spline` 10/3, `compare_brep_io.py` 6835/7.
    Amendments to the decisions, from the implementation: (a) the switch
    lies a quarter of the way from the turning point to where the gentler
    branch's slope over the profile's arc length has fallen to one (about
    two there: the slope grows as the distance's inverse square root), the
    window's near end at that point: half the way, the slowest crossing fuzz
    variant took a fifth to a third longer (a graph over the height costs
    more per point than one over the run, the run's more per piece the
    nearer it ends to the turning point); (b) a graph over the height's
    series about a point comes from Newton's steps on the series (thirteen
    coefficients in four steps), over a range coefficient by coefficient:
    Newton's series quotient over a range overestimated by `10^6` (the
    tests' wall at a quarter of the piece), the recurrence about a point
    widened the thirteenth coefficient to `±1.2` by the rounding of its
    twelve divisions; (c) the wall's Green integrals evaluate their ten
    moment tensors (degree up to about `5 p` in `ū`) through the Bernstein
    bases' jets made once per degree, one product of jets per row instead of
    de Casteljau's two hundred: about points and over ranges the integrals
    needed no more pieces (somewhat fewer), and the slowest loop variant
    took a third less (S9f.2b.1's meetings gain as much); (d) a tower
    polynomial's root whose point lies off the segment (another branch of an
    arc's implicit curve) is dropped before its repetition is asked, so only
    the segment's own tangencies are `Degenerate`; (e) a cylinder's halves
    each have their pieces: a turning point inside one half lies outside the
    other, whose graphs over the run end there as S9f.2b.1's. The `boolean`
    fuzz target's crossing variants (`SPLINE_CROSSING`, on) now reach loops
    and towers. Replays with debug assertions: the corpus (1,430 inputs) and
    the 24 regressions, none failing; the 1,431 S9f.2b variants of S9f.2b.1,
    none failing, 3,646 of their operations evaluating (3,232 before; the
    rest refused as documented, none naming S9f.2b.2), median 0.45 s and the
    slowest 6.8 s on a host at load 30 to 48; the 150 of them S9f.2b.1
    refused as S9f.2b.2's (loops and towers): 414 operations evaluating,
    median 0.85 s, ninth decile 3.7 s, the slowest 6.8 s (32 s under
    AddressSanitizer on this host at load 6 to 15, where the corpus's
    slowest input took 36 s). Every check of `HANDOFF.md`'s "Verification"
    holds: the suite, all 25 comparisons unchanged but this one with 0
    failures, the 13 generators importing `curved_boolean_reference.py`
    with `--check`, unittest (287 tests) and the ledger. DRAW survey: that
    of S9f.2b.2 and S9e.4b.1 (no case reaches its loops or towers).
    Campaign: the boolean campaign at `7199e06a`, S9e.4a's and
    S9f.2b.2's together (600 s, a sampled replay, `IMPORTED` and
    `SPLINE_CROSSING` on) clean, 845 runs, the slowest input 47 s under
    AddressSanitizer at load 12 to 17 (`17e131e3`, an existing corpus input,
    the torus against prisms). A replayed boolean variant with the shipped
    switches then panicked in the validator (`a conic or section
    evaluates`, at `7199e06a` and `93e6fcd0` alike): a holed square in the
    tilted frame, the fuse's first solid with a torus band (the prism
    itself), given to the `GIVEN_ROUND` cylinder, whose frame's normal is
    the tilted frame's normalized again, an ulp off it. The two cylinders'
    models crossed within rounding of parallel, S9c.2b.1's quartic gave a
    `Meet` of an ulp's sweep of the carrier's angle, and on the stored
    axes, exactly parallel, it has no point (latent since S9c.2b.1 and
    S9e.2's `GIVEN_ROUND`; the band is incidental, the holed prism alone
    panics too). Two cylinders whose axes are within `10^-12` of parallel
    without being parallel are now `Degenerate` in `meet::cyl_pair`, as a
    plane within rounding of a cylinder's direction is, unless certainly
    apart within their faces' bounds (the other's section beyond the first
    circle, within it or holding it by a certified margin over its axis's
    drift there), where they evaluate as before
    (`fuzz/regressions/boolean/crash-3e989fac80f9cef1f1a6abe1d874ec76855724c3.bin`,
    `tests/turned_booleans.rs`). Replays with debug assertions: the corpus
    (1,430 inputs), the 25 regressions, 477 variants beside the input and
    363 single-byte mutations of it, none failing; every comparison
    unchanged with 0 failures. An audit of the other pair kinds for the
    same class (two surfaces whose axes, centres or planes lie within
    rounding of a special relation without holding it) then ran each
    candidate through fuse, cut and common in the kernel, on frames built
    from another's normal normalized again or turned by an ulp (about 250
    chosen cases and 259 random ones in three tilted frames), and 2,304
    boolean variants with the shipped switches. One kind failed: a spline
    wall against a cylinder (S9f.2b). `spline_parallel::map2` takes only
    exactly parallel axes, so axes an ulp apart went to S9f.2b.1's
    crossing meeting, whose `A` (the `w^2` coefficient along the rulings)
    is then within rounding of zero and whose meeting within the faces is a
    sliver of the run: usually refused (`a wall's meeting's piece within
    rounding`), but 16 of the variants failed at `b903f3fa`
    (`degenerate_curve` from the validator, two `PrecisionLoss` where a
    projection's lift was left unpinned): a square with the lens hole, the
    cut's first result, given to the `GIVEN_ROUND` cylinder. Such a pair is
    now `Degenerate` in `spline_crossing::meeting_with` unless certainly
    apart within the faces' bounds (on every Bézier arc `|X|^2 - (r +-
    m)^2` of one sign at the overlap's middle height, `m` a rational bound
    of the ruling's drift over half the overlap's heights, by exact root
    isolation), where each operation is the same cylinder's on the wall's
    own frame (`fuzz/regressions/boolean/crash-a3fb9da3e3d83ef1cc3b272a30e6716eed8145b2.bin`,
    `tests/spline_crossing_booleans.rs`). Two whole tori of equal radii
    about centres within rounding, their axes within rounding of parallel
    or one normal's axes turned, gave `PrecisionLoss` (not a panic; two
    tori are off in the target, `TORUS_PAIRS`): they are now one surface,
    `Degenerate` (`torus_curved::torus_torus`, `tests/torus_curved_booleans.rs`).
    Found safe, each with its near cases evaluating or refused as
    documented: a spline wall against a sphere, cap or zone (a sphere's `A`
    is `|n|^2`; a rim's plane nearly normal to the wall's axis meets it in
    a crease, one nearly holding it is refused before); a cylinder against
    a cone nearly coaxial or parallel (the cylinder the carrier, its `A`
    within rounding of minus the cone's slope squared, rings over the
    turn; a cylinder within rounding of a ruling refused by `a_bound`) and
    two cones (unequal slopes: `A` within rounding of a nonzero constant;
    equal slopes: `A` within rounding of zero with simple roots, S9d.3b.2's
    split there, the far branch past `10^16` and the near one cancel-free
    in `meet_jet`); spheres centred within rounding of a cylinder's or a
    cone's axis (rings; a cylinder's equal radius within the resolution of
    a node); a torus against a cylinder, cone, sphere or torus nearly
    coaxial (`G` depending on `u` by rounding only, traced and verified
    exactly); planes nearly normal to a cylinder's, cone's, torus's or
    spline wall's axis, nearly through a cone's apex (refused, `a plane
    within the resolution of a cone's apex` or `a piece thinner than the
    resolution`, or the cone's own volumes) or nearly holding a torus's
    axis. The near coaxial results are the coaxial ones' (a cylinder in a
    frustum 4.97419, a torus's inner half in its major radius's cylinder
    9.33616 by Pappus, nested tori 12.33701). Of the `expect`s on a
    geometric evaluation, the validator's `curve_at` (`a conic or section
    evaluates`) is the one a near-degenerate edge reaches unguarded:
    `measure` (`TopologyParts::with_measured_enclosures`, before
    validation) takes every edge's end gaps, and the checks after it test
    only `curve_valid`'s parameters (a `Meet` of axes within `1e-6` of
    parallel is `degenerate_curve`, but a `Meet`, `Rise`, `Toric` or
    `WallMeet` that fails to evaluate elsewhere is not caught), so the
    kernel's refusals above keep it unreached; the engine's `recip`s off an
    apex or an axis (`cone.rs`, `procedural.rs`, `torus.rs`) are exact and
    guarded by exact refusals (a plane through the apex, the apex on the
    other's surface; a ring torus's points off its axis). Replays with
    debug assertions: the corpus (1,430 inputs), the 26 regressions, the
    2,304 variants and 605 single-byte mutations of the new input, none
    failing; every comparison unchanged with 0 failures.
  * **The validator's curve points guarded (2026-10-04).** After the
    near-parallel audit above, `validate::curve_at` returns `None` where a
    meeting, a section or a conic has no point in a tier
    (`projection::conic_point`'s `None`, until now `expect`ed: `a conic or
    section evaluates`), and each caller reports it instead of panicking.
    `measure` (`TopologyParts::with_measured_enclosures`, a builder's
    `measure_enclosures`) leaves that end's vertex without a bound
    (`root_bound_of`: `None` when neither tier evaluates), so the check
    reports `enclosure_missing` and a builder fails `Unrepresentable`.
    `check`'s `curve_ok` also requires the curve to evaluate, in either
    tier, at each end that has a vertex (`ends_evaluate`), else
    `degenerate_curve` on the edge, whose vertex and deviation checks are
    then skipped as for any degenerate curve. A `deviation` sample that
    does not evaluate leaves the tier undecided (`Unknown`, then the exact
    tier, then `uncertified_pcurve_off_edge`). No new issue kind: a curve
    whose parameters pass but which has no point at an end is degenerate,
    as a `Meet` of axes within `1e-6` of parallel already is. Where every
    point evaluates, which is every result so far (a failure panicked),
    nothing changes: the same bounds, verdicts, volumes and enclosures.
    The other `expect`s and `unwrap`s outside tests in
    `topology/validate.rs` and `validate/*.rs` are not evaluations that
    can fail on parameters `curve_valid` accepts: a spline's exact point
    inside its span's range (`SplineSpan` keeps the range in the domain),
    `arc_of` on the arcs its match leaves, exact conversions after
    `finite`, nonempty lists and bases; a projection pcurve's point
    already widens where it does not evaluate. A unit test
    (`curves_without_a_point_are_degenerate_not_a_panic`) measures and
    checks a wire edge on 7567a04b's `Meet` (an ulp's sweep on exactly
    parallel axes, which `curve_valid` refuses but `measure` panicked on)
    and, past `curve_valid`, on a `Toric` whose cylinder lies far off the
    torus and a `Rise` from a cone's apex: no bound on the start vertex,
    `degenerate_curve` and `enclosure_missing`, and `from_parts` refusing
    with the same issues (with the `expect` restored, the test panics).
    Checks: fmt, clippy, the 1.85 check, the release suite (551 tests);
    unchanged with 0 failures `compare_brep.py` (prisms 44/8, splines
    10/3, sheets 17/1), `compare_brep_io.py` 6835/7,
    `compare_turned_boolean.py` 2/13, `compare_given_curved_boolean.py`
    25/23 and `compare_spline_crossing_boolean.py` 14/37.
  * **A loop through one pole lifted at it (2026-10-04).** The merge of
    S9e.4b.3a's kernel with `GIVEN_MET` switched on (`86b1d834`) failed
    the boolean target's debug-assertion replay on corpus input
    `21928984` (`unexpected error invalid topology: uv_gap`), which passed
    at either parent (`GIVEN_MET` off at the one with S9e.4b.3a, whose
    kernel fails the same operations called directly). A bar in the tilted
    frame less a torus band; the cut's first result given to the
    `GIVEN_BALL` sphere about the middle of its first torus section, which
    lies on the bar's wall, so the wall holds the sphere's axis and its
    section is a meridian through the north pole inside the bar: the
    sphere face's loop runs up the meridian, through the pole, down the
    other side and back along the bar's bottom circle, its pcurves turning
    half a turn at the pole and the loop once in all. S9e.4b.3a's
    amendment (b) set such a loop's winding to none without lifting its
    pcurves, so its closing fin was a turn off its first unless the loop
    happened to start at the pole (as the pieces' fixtures' do), and the
    chained cut failed `uv_gap`; before it the face closed at a pole
    vertex loop beside the loop's own pole vertex. Not the validator's
    change and not the rounding: a box `4 x 2 x 2` with a ball of radius
    1.25 about a point of its wall `y = 0`, its axis along `z`, fails
    alike in 76 of 80 frames and operations at `86b1d834` (none at
    S9e.4b.3a's parent). The rule now takes only the pole the face would
    close at (`(total == 1) == forward`, as the pole vertex loop's), finds
    the loop's fin starting there and lifts that fin and the ones after it
    by the turn (`assemble::shift_u`, the loops' lift shared), so the half
    turn at the pole is the face's and the loop closes on its first fin;
    a loop starting at the pole is unchanged. The face is the one the pole
    vertex loop gave (its sector at the pole the same), one vertex fewer:
    the chained cut's volumes and areas those at `d82ffe24` (the merge's
    other parent, without S9e.4b.3a) within `3e-15`. A loop through the other pole only (the face holding the pole
    it closes at) keeps the pole vertex loop, as before S9e.4b.3a.
    `tests/sphere_booleans.rs`
    (`a_meridian_loop_through_one_pole_closes_on_its_first_fin`: the box
    and the ball about either pole, five reference directions, the four
    operations' volumes the closed forms' and the histories complete) and
    `tests/given_met_booleans.rs`
    (`a_ball_about_a_section_on_a_wall_through_its_axis`: the input's
    chain, the three operations' identities), and the input kept
    (`fuzz/regressions/boolean/replay-21928984a54f016cf7e3e987dcdd3f07df37c15c.bin`).
    Checks: fmt, clippy, the 1.85 check, the release suite (567 tests);
    replays with debug assertions of the corpus (1,432 inputs), the 27
    regressions and 540 single-byte mutations of the input (every value of
    the chained byte among them), none failing; the 30 boolean comparisons
    of `HANDOFF.md`'s table unchanged with 0 failures
    (`compare_imported_pieces_boolean.py` 32/13, `compare_sphere_boolean.py`
    30/0, `compare_spheres_boolean.py` 12/21, `compare_given_met_boolean.py`
    8/42 and the rest).
  * **A torus cap's round ends lifted along their branches (2026-10-04).**
    The `split` campaign at `6d26886d` found a whole torus (major 1.625,
    minor 0.875, the tilted frame) whose oblique split failed
    `InvalidTopology("an open loop in a torus piece")` in S8d.3's own
    split, before any piece reached S9e.4b.3b's Booleans. The plane cuts
    it in one contractible loop (a cap) reaching round the tube's inner
    side. `spiric::cap_edges` lifted both ends of each round end (the
    graphs over `v` through the turning points) to the turn nearest the
    turning point's `v`, assuming each half of a round end sweeps less
    than half a turn; here the `+` branch's half sweeps 3.42 radians, so
    the end ran the other way round, off the section, both round ends'
    sides came out reversed and the cap's loop did not close. Latent since
    S8d.3 (`d7049aac`, where the kernel call fails alike; `spiric.rs`
    unchanged since): the input fails alike
    at `788f8861`, `bb930aa0` (the frames' change), `a32d256f` (the
    arithmetic and charts) and `b99799ec`; not a regression. Each branch's
    `v` (`atan2` of the trace, `W > 0` on the cap, plus or minus `acos`)
    is continuous over the cap's `u`, so the ends now keep their branches'
    own values and the `-` end is lifted by the turn between the branches
    at the turning point (none where both are `psi`, one where they are
    `psi +- pi`). Where each half sweeps less than half a turn this is the
    old lift, so every fixture's pieces are unchanged. The cap and the rest
    split, the piece below the plane within 1e-6 of a midpoint quadrature
    over the tube with the radial integral exact (`tests/split.rs`,
    `a_long_torus_cap_splits_round_its_tube`, its history checked; the
    input kept,
    `fuzz/regressions/split/crash-366b6402fd50dddd7b52db2fd3e3de9dc02eed13.bin`).
    Checks: fmt, clippy, the 1.85 check, the release suite; replays with
    debug assertions of the split corpus (3,548 inputs), the 18
    regressions, 279 single-byte mutations of the input and 400 random
    whole tori cut by the target's oblique planes (2 of them failing alike
    at `a32d256f`), none failing;
    `compare_split.py` 72/56, `compare_split_pieces_boolean.py` 26/22 and
    the torus comparisons (11/24, 15/14, 15/29, 16/37) unchanged with 0
    failures.
  * **A plane within rounding of a cone's apex refused before its loops
    (2026-10-04).** The `split` campaign at `770bcdbc` found a frustum
    (radii 3.5 and 0.25 over 1.75, upright) split through its virtual
    apex's height `49/26` rounded by an oblique plane, whose piece below
    failed every Boolean of S9e.4b.3b's stage, the holding box's common
    too, as `InvalidTopology("an open Boolean of arcs in any position")`.
    The split takes the plane within the resolution of the apex and cuts
    rulings (S8d.2); the piece's model is the frustum common the plane's
    exact half-space, which misses the apex by about `1e-16`, so the
    section is a hyperbola within rounding of the rulings. Its edges'
    binary64 samples, a graph over the cone's angle (`rho = G / (mu + k
    (alpha cos + beta sin))`, both terms within rounding of zero) whose
    arms turn by less than an ulp, ran off the section: both of the cut
    face's loops round the frustum's section came out clockwise, both were
    nested as holes in the face's outer piece, and the assembly was open.
    The piece above reached `cone_curve3`, which refuses a plane whose
    section rounds to S8d.2's rulings (`Degenerate("a plane within the
    resolution of a cone's apex")`, the near-parallel audit's refusal), and
    a whole frustum against a slab with a face on that plane fails as the
    piece below did, at `770bcdbc` and at `f87d150d`: S9d.3a's plane
    section, reached through S9e.4b.3b, and not the coarse polygons
    `c82af2ec` samples again (these samples are off the section, not too
    few). The arrangement now refuses such a section before any edge of it
    enters (`ConeSec::rounds_as_rulings`, `cone_conic`'s own test on the
    cone model its edges round on, in `graph.rs` where a section's edges
    are made), with the refusal `cone_curve3` gave; a section whose edges
    lie outside both faces still meets nothing. Not the validator's change:
    every result that reached `cone_curve3` was refused there already.
    `tests/cone_booleans.rs`
    (`a_plane_within_rounding_of_a_virtual_apex_is_degenerate`: the
    frustum against the slab refused in all three operations, the slab
    moved `1/64` along its normal evaluating with the pair identities) and
    `tests/split_piece_booleans.rs`
    (`pieces_by_a_plane_within_rounding_of_a_virtual_apex`: both pieces
    refused against the holding box and the turned box), and the input kept
    (`fuzz/regressions/split/crash-baa7bc2d89a6658d6e58e7e931d5f60eb90a4265.bin`).
    Checks: fmt, clippy, the 1.85 check, the release suite (593 tests);
    replays with debug assertions of the split corpus and its regressions
    (3,568 inputs), the boolean corpus and its regressions (1,463 inputs),
    645 single-byte mutations of the input (every value of the mode byte
    among them; 337 failed alike before) and 500 random cones, zones and
    caps cut by the target's planes, none failing; the 31 boolean
    comparisons of `HANDOFF.md`'s table unchanged with 0 failures
    (`compare_imported_pieces_boolean.py` 32/13,
    `compare_split_pieces_boolean.py` 26/22, `compare_cone_boolean.py`
    25/5 and the rest), and `compare_split.py` 72/56.
  * **S9f.3 refined, before its code (2026-10-03).** Why it is refused
    today: `curved::spline_pairs` refuses a spline prism against a sphere
    or a cone ("a spline prism against a sphere or a cone (S9f.3)");
    behind that refusal `meet::section` and `meet::edge_surface` have no
    spline wall against a sphere or a cone face, no spline cap edge or
    crease against them and no sphere's circle or cone's rim against a
    spline wall (`spline_curved`: "a spline wall against a curved face in
    any position (S9f.3)"), the cone pairs' pass would take a spline wall
    for a quadric (`cones::cone_pair`), and the topology's
    `Curve3::WallMeet` holds a cylinder only. Decisions.
    (1) *Sub-steps.* **S9f.3a** (this step): a sphere, a cap or a zone
    (`Solid::sphere_with`) against a spline prism in any position, either
    the object; **S9f.3b**: a cone or frustum against a spline prism. They
    differ where the cone's radius term enters: along a spline wall's
    ruling a cone's function is `A w^2 + 2 B w + C` with `A = sum q_i^2 -
    k^2 q_h^2` of either sign (a ruling steeper than the cone's generatrices
    meets both nappes, one point on each and no turning point but where the
    ruling passes the apex; a ruling along a generatrix direction, `A = 0`,
    meets it once), the other nappe's points lie on the quadric but not on
    the face, the apex is a vertex on the meeting when the wall holds it,
    and the topology's cone is stored by its half angle, whose tangent is
    not rational (the certified evaluation's polynomials split by powers
    of `tan(a)`, enclosed per tier, where the sphere's and cylinder's are
    exact). S9f.3b refines those before its own evidence; this bullet's
    remaining decisions are S9f.3a's.
    (2) *The meeting.* Along the spline wall's ruling at the run parameter
    `tau`, `X = o + S_x(tau) x + S_y(tau) y + w n` (the spline prism's exact
    model), the sphere's function in the world's rows (`procedural::Other`
    of `other_sphere`: `|X - c|^2 - r^2`, its stored centre and radius) is
    `F = A w^2 + 2 B(tau) w + C(tau)` with `A = n . n` (positive: every
    ruling, so no condition on the axes), `B = n . (P - c)` of degree `p`
    and `C = |P - c|^2 - r^2` of degree `2 p` on each Bézier arc (`P = o +
    S_x x + S_y y`), and `D = B^2 - A C = A r^2 - |(P - c) x n|^2` (Lagrange)
    of degree `2 p`. These are S9f.2b's polynomials with three rows instead
    of a cylinder's two: S9f.2b.1's branches over `tau` between the roots
    of `D`, S9f.2b.2's graphs over the height about each turning point
    strictly inside both faces (the window's root, verified exactly by its
    (i) to (iii), switched at rational run parameters, the switches
    vertices) and the classification of turning points carry over
    unchanged (`spline_crossing::meeting` takes the partner's rows, the
    sphere's face and labels). A sphere smaller than the wall's height
    straddling it meets it in a loop with two turning points inside both
    faces (graphs over the height); one crossing a cap or larger than the
    wall in branches over `tau`, the turning points outside a face. The
    meeting is found once per pair of faces (`CylPair::Mixed`), for each
    hemisphere (S9d.1's split): a turning point inside one lies outside
    the other, whose branches over `tau` end there (S9f.2b.2's (e)).
    (3) *Representation.* In the engine `Crv::WallMeet` unchanged, its
    `other` the sphere's three rows, its carrier's partner a hemisphere
    (`FaceKind::Half`). In the topology `Curve3::WallMeet` gains
    `other_sphere: bool` (as `Curve3::Meet`'s): with it the other surface
    is the sphere of `other`'s origin (its stored frame, the centre) and
    `other_radius`, its function along the ruling `|L(u) + t M(u) -
    o|^2 - r^2`, and `topology/validate/wall_meet.rs`'s spans take the
    partner's rows (a cylinder's `x2` and `y2`, a sphere's three world axes,
    exact): `a = sum m_i^2`, `b = sum w_i m_i`, `c = sum w_i^2 - r^2` and `d
    = a r^2 - sum_{i<j} (w_i m_j - w_j m_i)^2` (Lagrange's identity, as the
    cylinder's two-row form), each made exactly once per wall and sphere as
    Bernstein polynomials, the graphs over `v` with `g = sum X_i^2 - r^2`.
    On the sphere face its pcurve is a `Projection` by the sphere's inverse,
    on its own wall its own `(u, v)`. Validity: S9f.2b's, the ruling's
    crossing of the cylinder's axis asked only of a cylinder. Rigid
    motion (the sphere's frame moves), the writer's and importer's
    refusals, history, `curve_curve` and `curve_surface` as S9f.2b's.
    (4) *Vertices, degrees and fields.* A spline cap edge or crease against
    the sphere: `F` along it, degree `2 p` (S9f.2b's `wallcrv_cyl` with the
    sphere's rows); the spline prism's vertical edges against the sphere:
    quadratic surds (S9d.1's `line_sphere`); the sphere's circles (`Circ`:
    the rims of a cap or zone and the split's great circle, each on its
    own sphere) against a spline wall: (a) in a plane not holding the
    wall's axis direction, its crease on the wall `w = h0 + h1 S_x + h2
    S_y` in `F`, degree `2 p`, each root a point in `Q(alpha)` and its place
    `(dx, dy)` on the circle read off rationally (S9f.2b's `conic_wall`
    with the circle's own sphere); (b) in a plane holding the wall's axis
    direction (a cap's or zone's split great circle when the sphere's axis
    is parallel to the wall's, a rim when perpendicular): the plane's
    generatrices on the wall meet the circle in the tower
    `Q(alpha)(sqrt(delta))`, and a circle of a surd radius has no rational
    chart (S9f.2b.2's half-angle tangent), so a primitive element by
    elimination instead: each arc's implicit equation `f(u, v)` (degree
    `p`, S9f.2a's `Implicit`) at the circle's projection `l0 + dx lx + dy
    ly`, reduced by the circle's equation `dx^2 |x|^2 + dy^2 |y|^2 = r2` to
    `E(dx) + dy O(dx)`, gives `R = E^2 - (r2 - |x|^2 dx^2) / |y|^2 O^2` of
    degree at most `2 p` in `dx`, and at each real root `dy = -E / O` in
    the same field `Q(dx)`; where `O` vanishes at a root on the segment
    (two of the circle's points share `dx`: a great circle whose `x` is the
    axis, its generatrices' two points) the roles of `dx` and `dy` are
    swapped, and where both vanish `ComputationLimit`. A root whose point
    lies off the segment (another branch of the implicit curve) is dropped,
    one found on two arcs (a knot) kept once, a repeated root on the
    segment (the circle tangent to a generatrix: a turning point on the
    circle) `Degenerate`. The sphere's poles are its vertices (on a spline
    wall, the engine's vertex-on-face rule). Every vertex on a wall lies in
    `Q(alpha)` of degree at most `2 p` or in `Q(sqrt(d))`; turning points
    are roots of `D`, degree `2 p`; `H` of S9f.2b.2's (iii), degree `4 p -
    2`, is only counted; no new limit.
    (5) *Degenerate*: the sphere tangent to the wall (a root of `D` of
    multiplicity above one inside both faces' closures: an isolated point
    or a node of the meeting, "a sphere tangent to a spline wall"); a
    turning point at an interior knot inside both faces' closures; a
    turning point on a face's boundary (a cap's edge, a rim, the
    segment's end; on the hemispheres' split, which is no edge of the
    input, the split is tried at another seam instead); a vertex's
    polynomial with a multiple root (an edge tangent to the other's face);
    and the engine's rules (a vertex of one input on the other's face,
    crossings within the resolution). A meeting within rounding of tangency
    validates or is `PrecisionLoss`, never taken for one; fixtures keep a
    margin.
    (6) *Stays refused*: cones and frustums (S9f.3b's), tori (by design,
    "S9f refined"), spline walls against spline walls on crossing axes,
    given results with spline walls or made from spheres against a spline
    prism (S9e's general faces), rational and periodic profile splines.
    (7) *Evidence first*: an independent reference of its own,
    `spline_sphere_boolean_reference.py` (the spline prism's exact model
    and profile parsing shared with `curved_boolean_reference.py`, nothing
    else): the pair sliced by planes, each slice's sections exact 2D
    regions (the profile cut by the prism's caps' lines; the sphere's
    circle projected along the prism's axis, cut by a cap's or zone's
    planes' lines), their Boolean pieces by Green's theorem on the cut
    boundaries, the breakpoints every vertex's slice and every edge's
    extreme ones (roots of exact polynomials: the meeting's by a resultant
    of the quadratic `F` and its tangent condition), each operation's
    volume and moments two ways (slicing along the prism's axis and along
    an oblique direction), the faces' classes (the walls swept along their
    generatrices, the sphere's by Archimedes' area element, the planar
    faces in their planes), closed forms of each input, inclusion and
    exclusion, the area identity, a Monte Carlo estimate, solids by the
    slices' union-find, and margins as S9f.2b's (turning points from faces
    and knots, crossings' sines at caps and circles, the discriminant's
    critical values, vertices from faces); `generate_spline_sphere_boolean_fixtures.py
    --check` (the same files under Python 3.9 and 3.12) across spheres
    straddling the bulge's, dome's, lens's and `knot`'s walls (loops),
    spheres larger than the wall or crossing a cap (branches), hemispheres
    on the prism's axis (their split great circles' towers) and on its
    side (their rims'), tilted prisms, either input the object, and
    declared degenerate pairs (a sphere touching the dome's arch, a
    turning point at `knot`'s knot, a loop turning back on a cap's edge);
    `test_spline_sphere_boolean_reference.py`; the generator in a CI group
    of its own (`spline-sphere`); a native capture
    `occt-boolean-spline-sphere-preimplementation` before
    `solid/boolean/curved/spline_sphere.rs` exists
    (`compare_spline_sphere_boolean.py` keyed on it); then the kernel in
    that file (the sphere's meeting, its cap edges' and circles' vertices)
    and `spline_crossing.rs` (the partner's rows), its tests, and the
    `boolean` fuzz target's spline prisms against its sphere, cap and zone
    tools (`SPLINE_SPHERE`).
    * **S9f.3a evidence (2026-10-03), before its code.**
      `spline_sphere_boolean_reference.py`, an independent reference of its
      own: the pair sliced by planes, each slice's sections exact regions
      of the profile's plane (the profile cut by the caps' lines; the
      sphere's circle projected along the prism's axis, an ellipse
      parameterised by the circle's own angle, cut by a hemisphere's line),
      their Boolean pieces by Green's theorem on the cut boundaries in
      closed form; breakpoints at every vertex's slice and every edge's
      extremes (the meeting's from `F` and its tangent condition, linear in
      `w`, eliminated exactly); the walls swept along their generatrices,
      the sphere's face by Archimedes' area element, a hemisphere's disc in
      its plane (by the profile's chords where it holds the axis); solids by
      the slices' union-find. `generate_spline_sphere_boolean_fixtures.py
      --check` writes 33 cases (9 fuses, 12 cuts, 12 commons; 27 solid, 6
      degenerate): loops (spheres straddling the bulge's wall, the lens's
      lower cubic, `knot`'s second span in `TILT`, the blob's wall under its
      top cap, and a sphere object across the wave's knot), branches over
      the run (a large sphere over the dome, its top cap inside it; a
      sphere under the blob in `TILT`, arches from its bottom cap's edge
      back to it), hemispheres (`hemi_bulge` through the bulge's bottom cap
      on its axis, its split great circle on the wall in a tower field;
      `side_hemi_bulge` on its side, its rim's plane holding the wall's
      axis, one turning point inside the hemisphere and one beyond its
      rim), and `degenerate` a sphere touching the dome's apex, a turning
      point at `knot`'s knot and a loop turning back on the bulge's top
      cap's edge. Checks: the two slicings (along the caps' normal and
      along `(2, -3, 5)`) within 1.8e-41 of the size, both inputs' closed
      forms 5.7e-42, faces' classes 2.8e-41, the area identity 8.0e-41,
      Monte Carlo (100,000 points a pair) within 2.8 standard errors;
      margins outside the declared pairs at least 0.023 (a vertex of the
      prism from the sphere; crossings at caps and the rim 0.071, turning
      points 0.27 outside a face and loops 0.34 inside both), the declared
      pairs' zero; Python 3.9 and 3.12 the same files;
      `test_spline_sphere_boolean_reference.py` (a sphere across a straight
      spline wall: its common, moments, areas and turning points in closed
      form; a rim's tower points in closed form); the CI group
      `spline-sphere`. The capture `occt-boolean-spline-sphere-preimplementation`
      (`compare_spline_sphere_boolean.py`, keyed on
      `solid/boolean/curved/spline_sphere.rs`; the probe `unsupported` on all
      33, refused by `spline_pairs`): every result valid with the
      reference's solids, 4 matches, 29 reviewed (BRepGProp's default
      integration up to 4.1e-5, the wave's up to 1.0e-3; a diagnostic
      build's adaptive BRepGProp and Green's theorem over OCCT's own faces
      and pcurves, the better within 6.0e-9 in volume, 7.8e-9 in area and
      8.5e-9 in the centre but on the wave's cut and common, 9.3e-7 within
      its edges' own tolerance of 3.7e-5); two solids' counts change when
      unified. No correction to the decisions from the evidence. S9f.3a's
      kernel next.
  * **S9f.3a implemented** (`solid/boolean/curved/spline_sphere.rs`;
    `spline_crossing.rs`'s `meeting_with` over the partner's rows and
    labels, `wallcrv_quadric`; `mod.rs`'s `spline_pairs`, `graph.rs`'s pairs
    pass, `meet.rs`'s sections and edge meetings, `assemble.rs`'s curves,
    `spline_parallel.rs`'s `Implicit::terms`; the topology's
    `Curve3::WallMeet::other_sphere` with `topology/validate/wall_meet.rs`'s
    spans over the other surface's rows and the validity rule): spline
    prisms against spheres, caps and zones in any position as the refined
    decisions describe, cones refused as S9f.3b's ("a spline prism against
    a cone (S9f.3b)"). The meeting is S9f.2b's over the sphere's three rows,
    found once per pair of a spline wall and a hemisphere; a sphere's circle
    meets a wall along its crease or, in a plane holding the wall's axis,
    at the roots of `E^2 - rho O^2` with `dy = -E / O`, the roles of `dx` and
    `dy` swapped where `O` vanishes (the hemispheres' great circles on the
    prism's axis, whose `x` is that axis, take the swapped order). All 33
    fixtures as the reference (the 27 results within the kernel's
    enclosures of volume, area and centre, each at most `1e-9` wide; the 6
    degenerate refused with the decisions' reasons); every history
    complete, results deterministic and moved rigidly
    (`tests/spline_sphere_booleans.rs`, 7 tests, 38 to 51 s in release on
    a host at load 15 to 50, 42 s at `opt-level` 2 with debug assertions: also the loops' graphs over `v` each in one knot span, a
    rim's tower vertices on the bulge, zones and caps in turned frames by
    inclusion and exclusion, a cone refused; the module's sixth test: a
    graph over `u` and one over `v` on a sphere with jets over points and
    ranges enclosing its points and slopes). `compare_spline_sphere_boolean.py`
    4 matches and 29 reviewed, the kernel within the reference on all 27
    results and none `unsupported` (the 18 loop cases' reviews name their
    counts: the kernel cuts each loop at its switches and at the
    hemispheres' split, OCCT's edges follow its own seams); every other
    comparison unchanged, `compare_split.py` 72/56, `compare_brep.py
    --family spline` 10/3, `compare_brep_io.py` 6835/7, `compare_step.py`
    23/6. Amendments to the decisions, from the implementation: (a) a
    sphere's loop switches a sixty-fourth of the way from its turning point
    to where the gentler branch's slope has fallen to one (a cylinder's
    stays a quarter): its graph over the height spans the sphere's section
    there, and at a quarter the slowest fuzz variants (spheres larger than
    a thin prism) spent ten times longer in the certified integrals over
    it, at a sixty-fourth and a two-hundred-fifty-sixth alike, the
    fixtures' tests no slower; (b) a turning point inside both faces but
    within the resolution of a prism wall's cap plane or a sphere's rim
    plane is `Degenerate` like one on it (a sphere about a turned prism's
    frame origin: a loop turning back `1e-17` above the cap, its graphs
    over the height cut there); (c) the binary64 root of a graph over `v`
    (`WallMeet::root_at`) builds the wall's pole rows once instead of at
    each of its bisection's steps; (d) the cone pairs' pass never sees a
    spline wall (`spline_pairs` refuses a cone first). The `boolean` fuzz
    target decodes a spline prism object against its sphere, cap or zone
    tool through S9f.3a (`SPLINE_SPHERE`), switched off for time: of 477
    variants of every third corpus input (the spline byte from 192 with its
    low bit, the flags' bits 5 and 6 clear) replayed with debug assertions,
    none failing, the median took 0.55 s, the ninth decile 1.5 s and the
    slowest 37 s (`d7599dbe`, a lens hole under a leaning sphere whose
    loops turn back near the lens's tips: the graphs over the run halved
    toward the turning points, sweep after sweep), 472 s under
    AddressSanitizer, another 297 s; the corpus holds `d7599dbe` itself.
    Replays with debug assertions of the corpus (1,430 inputs) and the 24
    regressions, none failing, with the switch off (the slowest 12.6 s) and
    on (48 s, `d7599dbe`). Every check of `HANDOFF.md`'s "Verification"
    holds: fmt, clippy, the 1.85 check, the suite, the table's 26 comparisons
    and the four beside it with 0 failures, unittest (296 tests) and the
    ledger. DRAW survey: that of S9e.4b.2, S9f.3a and S9f.3b, below (no
    case reaches it). Pending: the speed of the certified integrals beside
    a loop's turning points before it can be switched on.
    Campaign: the boolean campaign at `788f8861`, with S9e.4b.1, S9e.4b.2,
    S9f.3a, S9f.3b and the near-parallel guards (600 s, a sampled replay,
    `SPLINE_SPHERE` and `SPLINE_CONE` off) clean, 909 runs, the slowest
    input 45 s under AddressSanitizer at load 3 to 7 (`e36969f1`, an
    existing corpus input).
  * **S9f.3b refined, before its code (2026-10-04).** Why it is refused
    today: `curved::spline_pairs` refuses a spline prism against a cone or
    frustum ("a spline prism against a cone (S9f.3b)"); behind it
    `meet::section` and `meet::edge_surface` have no spline wall against a
    cone's wall, no spline cap edge or crease against it and no cone's rim
    against a spline wall (`spline_curved`), `spline_crossing`'s quadratic
    has no radius term (its rows, coefficients, slope and tangent are sums
    of squares less `r^2`; `meeting_with` asserts `t = 0`), the cone pairs'
    pass would take a spline wall for a ruled quadric (`cones::cone_pair`,
    never reached since S9f.3a's amendment (d)), and `Curve3::WallMeet`
    holds a cylinder or a sphere. Decisions, refining "S9f.3 refined" (1)
    for the cone; its (2) to (7) hold where nothing below differs.
    (1) *The meeting.* Along a spline wall's ruling `X = P(tau) + w n` the
    cone's function on its exact model (`procedural::other_cone`: `u^2 +
    v^2 - (b + k w_c)^2` in the cone frame's exact rows, `k` its rational
    slope) is `F = A w^2 + 2 B(tau) w + C(tau)` with `A = q_u^2 + q_v^2 -
    k^2 q_w^2` (`q` the prism's axis in the cone's rows), `B = q_u P_u + q_v
    P_v - k q_w rho(P)` of degree `p` and `C = P_u^2 + P_v^2 - rho(P)^2` of
    degree `2 p` (`rho(P) = b + k P_w` the radius term at the ruling's
    foot): S9f.2b's quadratic with a radius row of negative sign. `A` is
    one constant for the whole pair (every ruling is parallel to the
    prism's axis) and its sign is the case: (a) `A > 0`, the prism's axis
    farther from the cone's than its half angle: each ruling meets the
    quadric twice on one nappe or not at all, its turning points (roots of
    `D = B^2 - A C`) where it touches the cone, and S9f.3a's meeting
    carries over unchanged (branches over the run, loops' graphs over the
    height about turning points inside both faces, their switches; a
    turning point within the resolution of a rim's plane `Degenerate` as a
    sphere's, S9f.3a's amendment (b)); (b) `A < 0`, within the half angle
    (the prism's axis along the cone's among them): each ruling meets the
    double cone once on each nappe, so `D > 0` but where the ruling passes
    the apex, and the plus and minus branches run over the whole run, one
    on each nappe; the other nappe's lies beyond the apex, outside the
    cone's face (its heights: a cone's radii are nonnegative, so the apex,
    real or virtual, lies at or beyond an end), and the arrangement drops
    it as it drops any piece outside a face; (c) `A = 0`, the prism's axis
    exactly along a generatrix direction: one finite root per ruling,
    running to infinity where `B` vanishes, `Degenerate("a spline wall
    along a cone's ruling")`, as S9d.3c's cylinder along a cone's ruling
    (the limit of the rounding cases, whose second branch comes in from
    either side), checked first. (d) The apex, a cone's or a frustum's
    virtual one, on a spline wall's surface (the segment's wall, at any
    height) is `Degenerate("a cone's apex on the other input's surface")`,
    `cone_pair`'s rule: every ruling's two points meet there (a node of the
    meeting or an isolated point, `D`'s double root), checked before the
    turning points. The apex of a cone (a vertex of its model) inside the
    prism, outside it or on a face other than a spline wall is the
    engine's as before.
    (2) *Representation.* In the engine `Crv::WallMeet` unchanged, its
    `other` the cone's `Other` (`t = k`); `spline_crossing`'s rows carry
    their signs (a cylinder's two and a sphere's three positive less
    `r^2`, a cone's two positive and its radius row negative), so `A`,
    `B`, `C`, a point's `F_w / 2`, the tangent, the binary64 views and the
    graphs over the height are one code for the three quadrics; `Partner`
    gains `Cone` (its reasons: "a cone tangent to a spline wall", "a spline
    wall's meeting with a cone turning back at a knot" and "... on a face's
    boundary"), its loops' switch a sixty-fourth of the way as a sphere's
    until timing says otherwise. In the topology `Curve3::WallMeet` gains
    `other_half_angle` (as `Curve3::Meet`'s and `Toric`'s): nonzero, the
    other surface is the cone `|(w . x2, w . y2)| = other_radius + (w . n2)
    tan a2` of `other`'s frame (its stored base frame, `other_radius` its
    bottom radius, zero at an apex) and its function along the ruling `a
    t^2 + 2 b t + c` with `a = m_x^2 + m_y^2 - tan^2 m_n^2`, `b = w_x m_x +
    w_y m_y - (R + tan w_n) tan m_n`, `c = w_x^2 + w_y^2 - (R + tan w_n)^2`
    and `d = b^2 - a c = sum_i (R m_i + tan (w_n m_i - m_n w_i))^2 - (w_x
    m_y - w_y m_x)^2` (`i` over `x`, `y`; Lagrange's identity with the
    radius row). The stored half angle's tangent is not rational, so
    `wall_meet.rs`'s spans keep each of `a`, `b`, `c`, `d` as exact
    Bernstein polynomials per power of `tan` (`a0 + tan^2 a2`, `b0 + tan b1
    + tan^2 b2`, ...), combined at evaluation with `tan` enclosed in the
    tier (`cos_sin`, then a quotient, as `projection::meet_jet`'s), the
    graphs over `v`'s `g` with the radius row's square subtracted. The
    curve's evaluation and its `(-b + s sqrt(d)) / a` or `c / (-b - s
    sqrt(d))` choice are unchanged (`a` of either sign; at `a` near zero the
    second form, as `Meet`'s). On the cone face its pcurve is a
    `Projection` by the cone's inverse. Validity: a cone's half angle
    nonzero and under a right angle, its radius positive or zero (an apex),
    never both a sphere and a cone, the axis crossing asked of a cylinder
    only. Rigid motion, the writer's and importer's refusals, history,
    `curve_curve` and `curve_surface` as S9f.2b's.
    (3) *Vertices.* A spline cap edge or crease against the cone: `F` along
    it (degree `2 p`, `wallcrv_quadric` with the radius row); the spline
    prism's vertical edges: S9d.3a's `line_cone`; a cone's rim (a circle of
    rational radius `b` or `t` in its end plane) against a spline wall:
    along its plane's crease on the wall in the rim's own elliptic
    cylinder (the conic `c + a cos + b sin`'s, degree `2 p`, placed by its
    angle), or, in a plane holding the wall's axis direction (the cone's
    axis across the prism's), S9f.2b.2's half-angle chart (`tower_points`:
    the rim's radius is rational, so its chart is, and no elimination is
    needed); the cone's end discs against spline walls by S9f.1's creases
    and generatrices. Every vertex on a wall lies in `Q(alpha)` of degree
    at most `2 p` or in `Q(sqrt(d))`; no new limit.
    (4) *Degenerate*: (1)(c) and (1)(d); the cone tangent to a spline wall
    (a multiple root of `D` inside both faces' closures: "a cone tangent to
    a spline wall"); a turning point at an interior knot or on a face's
    boundary (a cap's edge, a rim, the segment's end) as S9f.3a's with the
    cone's reasons; a vertex's polynomial with a multiple root; and the
    engine's rules (a vertex of one input on the other's face, the rims'
    seams tried again). A meeting within rounding of these validates or is
    `PrecisionLoss`; fixtures keep a margin.
    (5) *Stays refused*: tori (by design), spline walls against spline
    walls on crossing axes, given results with spline walls or made from
    cones against a spline prism (S9e's general faces), rational and
    periodic profile splines.
    (6) *Evidence first*: `spline_cone_boolean_reference.py`, an extension
    of `spline_sphere_boolean_reference.py` (its prism, profile elements,
    Green's integrals of the pieces and solids reused; the cone its own):
    the cone or frustum on an exact frame (axes along the world's, either
    sense; any relative position comes from the prism's frame), the pair
    sliced along two directions that cut the cone in ellipses (its axis
    leaned toward the prism's), each slice's cone section the quadric's
    restriction to the slice (an ellipse about its centre on its principal
    axes) projected along the prism's axis and cut by the end planes'
    lines; the cone's wall area by its element `sqrt(1 + k^2) |G(s)| /
    D(theta)^2 dtheta ds` (`r = G / D` along the slice), in closed form
    over the section's arcs inside the prism; the prism's walls along
    their rulings (`F <= 0` within the end planes' slab, either sign of
    `A`); the caps by chords in their planes against the cone's conic of
    any type; the end discs as the sphere's; volumes and moments two ways,
    the inputs' closed forms, inclusion and exclusion, the area identity,
    Monte Carlo, solids by the slices' union-find; margins as S9f.3a's plus
    `|A|` and the apex's distance from the spline walls' surfaces.
    `generate_spline_cone_boolean_fixtures.py --check` (the same files
    under Python 3.9 and 3.12) across cones and frustums in either
    relation of the axes (`A < 0`: coaxial and leaning, both nappes met by
    every ruling, the apex inside the prism and outside it; `A > 0`: a
    leaning cone's loop through a wall, branches over the run, a cone
    across the prism's axis whose rims' planes hold it and whose cap
    sections are hyperbolas), either input the object, and declared
    degenerate pairs (the apex on a wall, `A = 0`, a cone touching a
    wall); `test_spline_cone_boolean_reference.py`; the CI group
    `spline-cone`; a native capture
    `occt-boolean-spline-cone-preimplementation` before
    `solid/boolean/curved/spline_cone.rs` exists
    (`compare_spline_cone_boolean.py` keyed on it); then the kernel
    (`spline_cone.rs` with `spline_crossing.rs`'s signed rows), its tests,
    and the `boolean` fuzz target's spline prism object against its cone
    tool (`SPLINE_CONE`).
    * **S9f.3b evidence (2026-10-04), before its code.**
      `spline_cone_boolean_reference.py`, an extension of S9f.3a's: the
      cone on an exact frame, the pair sliced along two directions that cut
      it in ellipses (its axis leaned toward the prism's), each slice's cone
      section its quadric's restriction to the slice on its principal axes
      projected along the prism's axis and cut by the end planes' lines,
      the Boolean pieces by Green's theorem on the cut boundaries; the
      cone's wall by its area element `sqrt(1 + k^2) |G(s)| / D(theta)^2`
      in closed form over the section's arcs inside the prism (an `atan`
      and a rational term, unwrapped across turns); the prism's walls along
      their rulings (`F <= 0` within the slab, either sign of `A`); the caps
      by chords against the cone's conic of any type; the end discs as the
      hemisphere's; breakpoints at every vertex's slice and every edge's
      extremes (the cap planes' conics' on the line where the quadric's
      gradient lies in the span of the cap's normal and the slicing
      direction, the meeting's by S9f.3a's elimination). Solids by
      intervals on chords of slices across the prism's axis, joined where
      they overlap (amending (6)'s "solids by the slices' union-find": a
      slice leaning from the prism's axis crosses a thin layer under a cap
      as a moving sliver whose sections never overlap in the profile's
      plane, and a cone across the prism's axis cuts the caps in
      hyperbolas, so no slicing of S9f.3a's kind serves).
      `generate_spline_cone_boolean_fixtures.py --check` writes 27 cases (7
      fuses, 10 cuts, 10 commons; 21 solid, 6 degenerate): `A < 0` a
      frustum on the bulge's axis across both caps (`bulge_frustum`), a
      cone with its apex inside the lens prism (`lens_apex`: the other
      nappe meets the lens's walls inside the prism's heights, outside the
      cone's face), a cone hanging over the blob with its apex below the
      prism (`blob_down`), a frustum object across the wave's middle span
      (`cone_wave`); `A > 0` a thin frustum across the dome's axis piercing
      its arch (`dome_pierce`, a loop), one whose bottom rim's plane holds
      the dome's axis (`dome_side_cone`: the rim's points on the arch's
      generatrices in a tower field, the top cap cutting the cone in a
      hyperbola), a frustum on `z` against `knot` in `TILT`
      (`knot_tilt_cone`, a loop); and `degenerate` the cone's apex on the
      dome's arch (`apex_wall`), the bulge in `TILT` against a cone whose
      slope is the stored axis's ratio (`ruling_tilt`, `A = 0` exactly), a
      cone across the dome's axis touching its arch (`dome_touch`). Checks:
      the two slicings within 1.1e-39 of the size, both inputs' closed
      forms 1.0e-39, faces' classes 1.4e-40 (the cone's wall by its
      element against `pi (b + t)` times its slant), the area identity
      9.2e-41, Monte Carlo (100,000 points a pair) within 2.1 standard
      errors; margins outside the declared pairs at least 0.0156 (`|A|` of
      the coaxial frustum; the apex 0.25 from the walls, crossings 0.16,
      loops 0.16 inside both faces), the declared pairs' zero; Python 3.9
      and 3.12 the same files; `test_spline_cone_boolean_reference.py` (a
      frustum halved by a straight spline wall in closed form, the area
      element over a turn, `A`'s sign, a rim's tower points in closed
      form); the CI group `spline-cone`. The capture
      `occt-boolean-spline-cone-preimplementation`
      (`compare_spline_cone_boolean.py`, keyed on
      `solid/boolean/curved/spline_cone.rs`; the probe `unsupported` on all
      27, refused by `spline_pairs`): every result valid with the
      reference's solids, 4 matches, 23 reviewed (BRepGProp's default
      integration up to 5.5e-6 off on the quadratic and cubic walls, 7.4e-5
      on `knot` in `TILT`, the wave's 1.0e-3; a diagnostic build's adaptive
      BRepGProp and Green's theorem over OCCT's own faces and pcurves, the
      better within 3.2e-8 in volume, 7.0e-9 in area and 9.0e-9 in the
      centre); three solids' counts change when unified. No other correction
      to the decisions from the evidence. S9f.3b's kernel next.
  * **S9f.3b implemented** (`solid/boolean/curved/spline_cone.rs`;
    `spline_crossing.rs`'s signed rows (`terms`), `conic_crease` and
    `tower_points` over a partner, `inner_gap`'s cone; `mod.rs`'s
    `spline_pairs`, `graph.rs`'s pairs pass (the cone pairs' pass leaving
    spline walls to it), `meet.rs`'s sections and edge meetings,
    `assemble.rs`'s curves; the topology's `Curve3::WallMeet::other_half_angle`
    with `WallMeet::radius_row`, `topology/validate/wall_meet.rs`'s
    polynomials per power of the half angle's tangent and the validity
    rule): spline prisms against cones and frustums in any position as the
    refined decisions describe. The meeting is S9f.2b's over the cone's
    signed rows, found once per pair of a spline wall and the cone's wall
    after the decisions' two refusals (`A = 0`; the apex, real or virtual,
    on the wall's surface); `A > 0` gives S9f.3a's loops and branches (the
    loops switched a sixty-fourth of the way, as a sphere's), `A < 0` the
    two branches over the whole run on opposite nappes, the other nappe's
    dropped as outside the cone's face. A rim meets a wall along its
    plane's crease in its own elliptic cylinder (the rows dual to the
    conic's axes) or, in a plane holding the wall's axis, at S9f.2b.2's
    points in its half-angle chart. All 27 fixtures as the reference (the
    21 results within the kernel's enclosures of volume, area and centre,
    each at most `1e-9` wide; the 6 degenerate refused with the decisions'
    reasons); every history complete, results deterministic and moved
    rigidly (`tests/spline_cone_booleans.rs`, 6 tests, 20 s in release on a
    host at load 10 to 17, 23 s at `opt-level` 2 with debug assertions:
    also the loops' graphs over `v` each in one knot span and the coaxial
    frustum's branch on its own nappe, a rim's tower vertices on the dome,
    a leaning frustum (`A > 0`) and a tilted one (`A < 0`) in turned frames
    by inclusion and exclusion; `wall_meet.rs`'s seventh test: cones
    across and along a wall, their points on both surfaces, jets enclosing
    them, `a`'s sign). `compare_spline_cone_boolean.py` 4 matches and 23
    reviewed, the kernel within the reference on all 21 results and none
    `unsupported` (the 6 loop cases' reviews name their counts: the
    kernel's loops cut at their switches; in `knot_tilt_cone_fuse` the
    cone's wall outside the prism, two regions, two faces of OCCT's and one
    of the kernel's, its cone wall one face traced in projection); every
    other spline comparison unchanged with 0 failures
    (`compare_spline_sphere_boolean.py` 4/29, `compare_spline_crossing_boolean.py`
    14/37, `compare_spline_parallel_boolean.py` 31/12,
    `compare_spline_any_boolean.py` 22/16). Amendments to the decisions,
    from the implementation: (a) a turning point inside both faces but
    within the resolution of a cone's rim's plane is `Degenerate`, as a
    sphere's rim's (S9f.3a's amendment (b), `inner_gap`); (b) a sphere's
    and a cylinder's binary64 curve and certified polynomials are as
    before bit for bit (the radius row's terms only where there is one).
    The `boolean` fuzz target decodes a spline prism object against its
    cone or frustum tool through S9f.3b (`SPLINE_CONE`), switched off for
    time: of 405 variants of every third corpus input (the spline byte in
    `161 + 4 m`, the flags' bits 5 and 6 clear) replayed with debug
    assertions on a host at load 26, none failing, the median took 0.47 s,
    the ninth decile 2.0 s and the slowest 12.6 s, and under
    AddressSanitizer the three slowest 97, 76 and 49 s (a cone on its side
    or tilted across a lens hole, its loops: the certified integrals beside
    the turning points, as `SPLINE_SPHERE`'s). Replays with debug
    assertions of the corpus (1,430 inputs) and the 25 regressions, none
    failing, with the switch off (the slowest 18 s) and on (16 s). Every
    check asked holds: fmt, clippy, the 1.85 check, the suite, the five
    spline comparisons with 0 failures, the generator's `--check` under
    Python 3.9 and 3.12 and its unit tests. DRAW survey: that of S9e.4b.2,
    S9f.3a and S9f.3b, below (no case reaches it). Pending: the Linux
    record of the capture, and the speed of the loops' certified
    integrals before `SPLINE_SPHERE` and `SPLINE_CONE` can be switched on.
    Campaign: the boolean campaign at `788f8861`, with S9e.4b.1, S9e.4b.2,
    S9f.3a, S9f.3b and the near-parallel guards (600 s, a sampled replay,
    `SPLINE_SPHERE` and `SPLINE_CONE` off) clean, 909 runs, the slowest
    input 45 s under AddressSanitizer at load 3 to 7 (`e36969f1`, an
    existing corpus input).
  * **DRAW survey of S9e.4b.2, S9f.3a and S9f.3b (2026-10-04,
    `UPSTREAM_TESTS.md`).** At `12b6c176` (`s9c2-kernel` with S9e.4b.2,
    S9f.3a, S9f.3b and the near-parallel guards merged; the public dataset,
    120 seconds a case, four at once). The 1,802 self-contained cases of
    the Boolean group on both backends: every status and every refusal's
    reason the last survey's (`93e6fcd0`) field for field, 987 evaluating
    and registered, 592 refused, 223 unsupported on both, none failing or
    timing out, the sentinels refused as before, `bopfuse_simple/ZP6` and
    the `gdml_public` tori too. The 1,814 cases restoring a shape for a
    Boolean on the Rust adapter, and the 171 the import reaches on native
    DRAW too: only S9e.4b.2's 7 polyhedra move, as its trial found: 4
    evaluate on both backends with every check (`bugs/modalg_1/buc60803`,
    `bug102_1`, `bug102_2`: two frustums of a pyramid, the second on the
    first's top, fused; `bopfuse_complex/K5`: a frustum on a box's top), 2
    are refused by the adapter as the next Boolean's argument
    (`bugs/modalg_2/bug578_1`, `_2`: the frustums' bases 6.6e-7 to 2.0e-6
    apart, their fuse two solids) and 1 by S9's rules (`bfuse_complex/D9`:
    a corner the two files share stored 1e-13 apart, a face thinner than
    the resolution). Of the 171: 27 evaluate on both backends, 56 are
    refused by S9's rules, 77 are bodies none of the kernel's constructions,
    6 S9e.4b.4's, 3 spline bodies (S9f) and 2 arguments of several solids;
    24 are unsupported natively too (`checksection`, `bopargcheck`). No other
    status or reason moves; the same 8 time out at their first restores.
    S9f.3a's and S9f.3b's spline walls reach no case: the Boolean group's
    spline solids are the 96 cases converting boxes by `nurbsconvert` (no
    sphere or cone among them), which neither host runs, and the 3
    restored spline bodies (`bcut_complex/L9`, `O1`, `bfuse_complex/N7`)
    are refused as imported spline faces (S9f) before any wall meets a
    partner. The near-parallel guards (two cylinders, a spline wall and a
    cylinder, two tori within rounding of parallel or of one surface)
    refuse no case: none of their reasons appears, and no registered case's
    status or value moves. The volume audit (`vprops` and `sprops` before
    each `checkprops`): the 4 new cases' volume, area and centre native
    DRAW's to its printed digits (`buc60803`'s 805,475,092.02 and
    18,606,382.71; `K5`'s area 448,769.50, its `checkprops -s` 448,769);
    the 1,010 registered cases' values the last audit's (`93e6fcd0`'s
    1,003 and 7) bit for bit on both backends; the same 35 disagreements,
    native off in each. The 4 are registered (1,100 cases; the contract
    holds for them on both backends within 30 seconds, Rust 0.4 to 1.1 s);
    the ledger does not change (no `checknbshapes` among them). A full
    contract run of the manifest holds on both backends with the dataset
    (30 seconds a case), the slowest Boolean case 16.9 seconds
    (`bopcommon_simple/ZK8`, on a host at load 5 to 10; `boptuc_simple/ZK8`
    16.2 in the last survey, 13.3 here), the restore cases 0.1 to 3.9 s,
    the rollex 2.1 to 3.0. No case fails, crashes or panics; no kernel
    change.
  * **CI on S9e.4b.1 (2026-10-04).** The "Minimum Rust 1.85" job (debug
    assertions and overflow checks on Linux) failed at `788f8861` in
    `tests/imported_arc_booleans.rs`: five tests panicked at the Boolean's
    debug history check, two `split_support_differs` each. It is (a) and
    (j) of "CI on `rust-kernel`" above again, glibc's correctly rounded
    `hypot`: a prism built its walls' frames by `Frame3::new` on its own
    frame's axes, normalizing the normal again, which is not idempotent and
    not the same on every host (the slot's tilted axis normalized once,
    `(0, 0.5999999999999999, 0.8)`, stays under macOS's `hypot` and turns
    to `(0, 0.6, 0.8000000000000002)` under glibc's, which macOS's turns
    back). An imported prism's construction (S9e.4a) takes its bottom
    cap's stored frame, normalized once by the converter as each stored
    cylinder's axis is, so on Linux its walls, and the Boolean's split
    walls on them, lay an ulp off the stored cylinders, and the history
    check asks a split cylinder for an axis exactly parallel to its
    source's. Reproduced on macOS with `Vec3::length` made correctly
    rounded (exactly, by big integers) in a local patch: the same five
    tests fail with the same entity ids. Fix: the walls derived from a
    prism's frame take its axes bit for bit (`Frame3::at`): `Topology`'s
    prism's arc and circle walls, S8a.2's plane split's walls
    (`level_frame`) and S9a.2's stacked Booleans' cylinders, which must
    agree with the prism's (the prism's walls alone, under the emulated
    rounding, broke three of `tests/split.rs`'s histories, the split's
    pieces normalized again). `history::check` is unchanged. Regressions,
    each failing on macOS without the fix, both frames run on any host:
    `snapped::tests::split_walls_keep_the_stored_axes_on_either_platforms_frames`
    (the slot with its caps' and cylinders' stored frames macOS's and
    glibc's normalized axes bit for bit, fused with and cut by a box) and
    `oblique::tests::split_walls_keep_their_prisms_axes_on_either_platforms_frames`
    (a stadium prism on either frame split by an oblique plane and fused
    with a box across its middle in its frame). Under the emulated
    rounding with debug assertions the whole suite passes but
    `properties_baseline`, whose committed baseline is macOS's (CI
    regenerates it per host); its rows are bit for bit those before the
    fix there, and the committed baseline holds here. S9e.4a's and
    S9e.4b.2's paths had no such failure: `imported_booleans.rs` and
    `imported_polyhedra_booleans.rs` pass under the emulated rounding with
    debug assertions before and after the fix (the polyhedra have no
    cylinder). Checks: fmt, clippy, the 1.85 check, the 1.85 suite with
    debug assertions at `opt-level` 2, the release suite;
    `compare_imported_boolean.py` 54/15, `compare_imported_arcs_boolean.py`
    27/9 and `compare_imported_polyhedra_boolean.py` 47/1, none failing.
  * **CI on the boolean target's imported stage (2026-10-05).** The
    per-push "Fuzz / boolean" job on Linux crashed at `dd51af05` replaying
    `regressions/boolean/replay-26c72abf…` (a corpus input): two
    `split_support_differs` at the Boolean's debug history check. It passed
    on macOS and failed alike at `ff8c914f` under the emulated rounding of
    the note above (`Vec3::length` correctly rounded in a local patch),
    with the same ids. It is the note above's kind again, at the caps: a
    stadium prism in the tilted frame is written, read back, imported
    (S9e.4a's stage, `IMPORTED`) and cut by the target's cone, which
    crosses its wall about the origin. Since `e961e477` a prism's walls
    take its frame's axes bit for bit, but `Topology`'s prism still built
    its caps by `Frame3::new` on its normal, normalizing it again:
    `(0, 0.6, 0.8)` to `(0, 0.5999999999999999, 0.8)` on every host, so the
    caps were written an ulp off the walls' axis. The import normalizes
    each stored frame once more: glibc's `hypot` turns the caps' to
    `(0, 0.6, 0.8000000000000002)` and keeps the walls' at
    `(0, 0.5999999999999999, 0.8)`, and the imported prism's construction,
    on its bottom cap's frame, split the walls an ulp off the stored
    cylinders. macOS's `hypot` keeps that axis, but turns `(0, 2, 3)`'s
    normalization again where glibc's keeps it, so the same failure was
    there for such frames on macOS. Fix: a prism's caps take its axes bit
    for bit (`Frame3::at`, the bottom `flipped`), so its caps and walls are
    written on one axis up to sign and imported on one, whatever the
    platform's `hypot`. `history::check` is unchanged. Regression, failing
    on macOS without the fix on its second frame and under the emulated
    rounding on its first, both frames run on any host:
    `imported::tests::imported_caps_and_walls_share_their_axis_on_either_platforms_frames`
    (the input's stadium on frames whose normals are `(0, 3, 4)`'s and
    `(0, 2, 3)`'s normalizations bit for bit, written, read back, imported,
    fused with, cut by and in common with the cone, and its caps' normals
    and walls' axes one bit pattern up to sign); `fuzz/regressions/README.md`
    describes the input. Sweep under the emulated rounding with debug
    assertions and overflow checks, one process per input: before the fix 7
    distinct boolean corpus inputs failed alike (12 files with the
    regression and duplicates), and 3,103 of the input's 4,080 single-byte
    mutations; after it none of the 1,466 boolean corpus inputs and
    regressions, the 3,569 split ones or the 4,080 mutations fails, and the
    kernel suite passes but `properties_baseline` (macOS's baseline, as
    above). Natively the suite passes with `properties_baseline` too: the
    caps' frames move no mass property. Checks: fmt, clippy, the 1.85
    check, the release suite (606 tests); every comparison of
    `HANDOFF.md`'s table with its matches and reviews unchanged, none
    failing, and `compare_split.py` 72/56, `compare_brep.py --family spline`
    10/3, `compare_brep_io.py` 6,835/7 and `compare_step.py` 23/6; the
    native replays with debug assertions of both corpora and their
    regressions, none failing. In the DRAW survey of S9e.4b.3c.1 and the
    switches' speed-up the fix moves 4 centre coordinates of two
    registered restore cases (`bfuse_complex/F5`, `Q2`) within rounding,
    at most 4.6e-16 of the solid's size, still native DRAW's to its
    printed digits; no status.
  * **CI's scheduled full replay at `59d0c57b` (2026-10-05).** The
    scheduled "Rust geometry fuzzing" run replayed the full corpora under
    AddressSanitizer on Linux and found two crashes while mutating. (1)
    `boolean/crash-43d93718…`: `unexpected error coordinates cannot
    resolve the requested linear tolerance`, `Error::PrecisionLoss`, which
    a Boolean does not document. It fails natively at `973b0bb4`. The
    cut of a holed prism by the sphere tool's cap is given to the
    `GIVEN_BALL` sphere about the middle of its first meeting of two
    curved faces (`GIVEN_MET`), the hole's wall with the cap's sphere; that
    point lies on the hole's cylinder within rounding and the ball's axis
    is parallel to the hole's, so the ball meets the wall in a loop through
    its north pole within about `1e-16`, inside the wall. Exactly through
    the pole the graph's pole pass refuses it (S9d.1's `OutOfDomain`, a
    section through a pole off its meridians); within rounding the pass
    (exact for curves other than circles) misses it, and the loop's
    projection pcurve on the ball turns half a turn at the pole, which
    `loop_fins`' 16, 64 and 256 anchors cannot pin (S9d.4b.2's addition):
    `PrecisionLoss`. A ball centred `1e-4` off the cylinder fails alike,
    `1e-3` off evaluates: the anchors are a computation budget, not the
    coordinates' resolution. Fix: a projection whose lift 256 anchors
    leave unpinned is `ComputationLimit` (`assemble::loop_fins`), the
    exact case unchanged. Regression: `tests/spheres_booleans.rs`,
    `a_section_within_rounding_of_a_spheres_pole_is_a_computation_limit`
    (the holed prism's direct form: a ball on the hole's cylinder exactly,
    `OutOfDomain`; an ulp off and the fuzz input's centre,
    `ComputationLimit` for fuse, cut and common; `0.0125` off, all three
    evaluating with the pair identities). (2) `split/crash-93175910…`:
    `4.4e-16: a piece below its plane`, the target's check of a body
    piece's centre. It passes natively and fails under the emulated
    glibc rounding at `973b0bb4`. An S8e closed wire (a square `17.5`
    across in the tilted frame) is split by a plane through a point of its
    plane along three times its frame's normal, which `Frame3::new`
    normalizes again: an ulp's tilt, whose trace `a u + b v + d = 0` misses
    the square on macOS and crosses it under glibc's `hypot`, cutting the
    wire into two open wires each within about `1e-15` of the plane, on
    its side only exactly (the target's binary64 centre test cannot tell).
    The kernel was the one to fix, not the assertion: the pieces' sides
    were decided below the resolution, as a prism's cap within it is
    refused (`a vertex within the resolution of the plane`, `a split
    within binary64 of a cap`), and flipped with an ulp of the input. Fix:
    `Body::split_by_plane` refuses a plane within the resolution of the
    body's plane over the whole body (`|d| + |(a, b)| R` within it, `R`
    the profile's reach from the frame's origin over its points, arcs'
    and circles' extents and splines' poles) that would split it,
    `Degenerate`, sheets and wires alike; an exactly parallel plane and a
    body on one side are returned as before. Regression:
    `tests/sheet_splits.rs`,
    `a_plane_within_the_resolution_of_the_bodys_plane_is_degenerate` (the
    input's square as a wire and a sheet, normals an ulp or two off the
    frame's through its centre refused or returning the body, never
    split; `1e-3` off the plane the body whole below; a `1e-6` tilt still
    splitting it in two with its length or area kept). Both inputs are in
    `fuzz/regressions/` with README entries. Sweeps with debug assertions,
    one process per input: of the inputs' single-byte mutations under the
    emulated rounding 3,108 of 4,335 (boolean) and 7,086 of 8,670 (split)
    failed alike before the fix, none after, natively or emulated; the
    corpora and regressions (1,467 boolean and 3,570 split inputs) replay
    without a failure natively and emulated. Checks: fmt, clippy, the 1.85
    check, the release suite (623 tests); every comparison of
    `HANDOFF.md`'s table with its matches and reviews unchanged, none
    failing, and `compare_split.py` 72/56.
    Campaigns at `f614b6ce` (600 s each, sampled replays, every switch on):
    `boolean` clean, 841 runs, the slowest input 18 s under
    AddressSanitizer; `split` clean, 1,474 runs, the slowest 18 s, at load
    about 10.
  * **Equal cylinders crossing at right angles, a result touching itself
    (2026-10-06).** S9e.4b.4b.2's evidence found K1's configuration on
    kernel-built prisms decided three ways: a prism with a round hole and a
    rod of the hole's radius whose axis crosses the hole's at right angles
    (exact frames) fused as `InvalidTopology("non_manifold_vertex")`, cut
    evaluating, common `Degenerate("solids touching at a vertex")`. The
    pair is S9c.1's, decided exactly: equal circular cylinders with
    crossing axes meet in two ellipses, which cross at the two points where
    the walls are tangent (`meet.rs`'s `cyl_pair`, `CylPair::Crossing`;
    S9c's refusal of tangencies between the inputs covers tangencies that
    leave no such curves, not these), and whether a result keeps both sides
    of a wall there is its own: as declared for the Steinmetz fixtures
    (fuse and common evaluating, the cut two pieces touching at those
    points, `degenerate`). With a hole the sides swap: the cut keeps the
    rod's wall's outside and is a manifold solid, the common is the rod's
    two halves either side of the hole touching there (two shells meeting
    at a vertex, refused as before), and the fuse is one solid whose void
    is the hole's two halves either side of the rod, touching at those
    points: one shell meeting itself at a vertex. `assemble.rs` refused
    shells meeting each other at a vertex, not a shell meeting itself, so
    the fuse was built and the validator found the pinch. Root cause and
    fix: the assembly now checks every shell before building it, linking
    the arrangement's edges at each vertex where a face's loop runs from
    one to the next (the validator's link test, on the exact arrangement)
    and refusing a vertex whose edges fall into more than one fan as
    `Degenerate("a result touching itself at a vertex")`, S9a's rule and
    reason; the validator is unchanged. The three operations are now each
    decided exactly as before, the fuse refused for its own result. Other
    crossings of equal radii were consistent already: in turned frames
    (the tilted frame's rod across its hole) S9c.2b.1's near node
    (`two cylinders' section within the resolution of a node`), at an
    oblique angle (axes `z` and `(0, 3, 4)`) a tangency (`a tangency
    between the inputs (S9c)`), all three operations alike, holes or rods.
    Radii `1` and `1 + 2^-k` evaluate in all three frames for `k` = 10, 30
    and 40 with the pair identities (rods and holes); nearer, outside this
    fix: in exact frames S9c.2a has no near-node rule, and at `k` = 50 and
    52 the validator leaves a loop's winding undecided
    (`ComputationLimit`, rods and at 52 holes too); in turned frames `k` =
    50 and 52 are near nodes, and at `k` = 46 (the two rings a few
    resolutions apart where the walls nearly touch) the holes evaluate,
    the tilted rods are `ComputationLimit` but the oblique rods fail
    validation (`InvalidTopology("enclosure_exceeds_resolution")`), the
    near-node margin's own question, left open (the fuzz target's radii,
    eighths, cannot reach it). The boolean fuzz target reaches the fixed
    case: a square with a round hole and the stadium tool of its radius
    stood on its side through it (`pick` 128, `dx = dy = 0`), its first
    arc's axis crossing the hole's, crashed with `unexpected error invalid
    topology: non_manifold_vertex` before the fix; that hand-written input
    is kept (`fuzz/regressions/boolean/replay-47992d6f….bin`). Regressions:
    `tests/curved_booleans.rs`,
    `a_hole_and_a_rod_of_its_radius_crossing_it_touch_at_two_points` (the
    fuse touching itself, the common two solids touching, the cut its
    closed form `144 - 10 pi + 16 / 3`, either way round: the rod less the
    slab three solids, `2 pi + 16 / 3`; the tilted frame's near node and the
    oblique tangency in every operation; radius `1 + 2^-10` evaluating in
    all three frames with the pair identities) and
    `a_stadium_on_its_side_through_an_equal_hole` (the fuzz input's
    geometry: fuse and common each one solid touching itself, the cut its
    closed form). Checks: fmt, clippy, the 1.85 check, the release suite
    (648 tests); every comparison of `HANDOFF.md`'s table with its matches
    and reviews unchanged, none failing; replays with debug assertions of
    the boolean corpus and its regressions (1,469 inputs) and the split
    corpus and its regressions (3,571), none failing; of the input's 3,570
    single-byte mutations 1,527 failed alike before (every one
    `non_manifold_vertex`) and none after; the new tests under the emulated
    glibc `hypot` with debug assertions and overflow checks (8.6 s). DRAW
    survey: that of S9e.4b.4b.2a, S9e.4b.4b.2b.1 and the cylinder pairs'
    fixes at `826346b7` (no case refused so, none moving with the rule).
    Superseded in part (2026-10-09) by "Equal cylinders with crossing axes"
    (branch `steinmetz-tangency`, the user's decision: contact only): where
    the faces' outward normals are opposite at those points (a hole and a
    rod of its radius, the stadium through the hole) the inputs touch, a
    tangency between the inputs in every operation, so the cut no longer
    evaluates and the fuse and common are refused as that tangency before
    the assembly; the two tests above assert it. The assembly's check of a
    shell touching itself at a vertex stays, with that note's fallback in
    `polyhedra.rs`, for any other path.
  * **Crossing cylinders near a node, enclosures past the resolution
    (2026-10-06).** The open item the previous note left: radii `1` and `1 +
    2^-k` crossing, between the near-node rule's margin and exact
    evaluation. Mapped directly over `k` = 36 to 56, the exact, tilted and
    oblique frames, rods and the slab's hole, either radius the larger, and
    the slab also standing on the plane through the near nodes (its bottom
    face's vertices there), fuse, cut and common alike in every row:
    `InvalidTopology` in the exact frames with the radii swapped (the `1 +
    2^-k` cylinder on `z`; rods and holes) for `k` = 43 to 52
    (`enclosure_exceeds_resolution` to 45, with it
    `uncertified_loop_winding` from 46), the oblique rods (either order) for
    44 to 49, and the standing slab with the larger hole for 45 to 49
    (exact), 45 to 47 (tilted) and 47 to 49 (oblique); elsewhere in the band
    results evaluated, the tilted rods were `ComputationLimit` at 46 and 47,
    S9c.2b.1's near node began at 48 (tilted) or 50 (oblique), and the exact
    frames had no near-node rule (S9c.2a), unswapped rods there
    `ComputationLimit` from 50. From 53 the radii are equal in binary64, the
    previous note's cases. Root cause: a stored meeting (`Curve3::Meet`) is
    the root `w = (-B + s sqrt(D)) / A` of a ruling's quadratic, its
    coefficients rounded from the stored frames, so `D` is uncertain by
    about `eps A L^2` (`L^2 = 2 rho^2 + r^2`, `rho` the reach from the
    carrier's stored origin to the other's axis plus its radius, `r` the
    other's radius: `51` for the exact frames' rods from `-4`), and near an
    extremum of `D` the height is uncertain by about `eps L^2 / sqrt(D)`.
    The validator encloses what lies there: a ring edge's closing point (the
    carrier's angle `0`, which in these frames is a near node: the face's
    loop gap) or a vertex on a plane through the node. Measured, those
    enclosures are `kappa eps L^2 / sqrt(D)` with `kappa` 1.3 to 5.0 over
    the failing rows (doubling every two `k`: the exact frames' swapped
    rods' ring `3.6e-8` at `k` = 40, `1.03e-7` at 43). The near-node margin,
    `|D| < (A res / 2)^2` (the branches within the resolution), takes no
    account of the image's rounding: at `k` = 46 the branches are a few
    resolutions apart while the image's enclosures are three times the
    resolution; rings closing elsewhere (the exact frames' unswapped rods,
    closing a quarter turn from the node) passed. Fix:
    `turned::conditioned_node` refuses an extremum of either cylinder's `D`
    with `|D| < max((A res / 2)^2, (8 eps L^2 / res)^2)` as the same near
    node (`two cylinders' section within the resolution of a node`), in
    turned frames after `near_node` and now in exact frames too
    (`procedural::perpendicular`, after its exact tangencies), `L^2` from
    the faces' stored surfaces (`meet::reaches`, given to `cyl_pair`). Each
    chart at an axis point is read over the quarter turns either side of it:
    an extremum at a chart's antipode shows as critical points of its
    polynomial far out (`|t|` about `10^6` here) where `|D|` is a fraction
    of its extremum's, which refused `k` = 40 on a first try. Eight times
    the margin keeps the largest measured `kappa` 1.6 times inside it and
    `k` = 40 1.5 times outside it (exact frames; 1.9 oblique). The validator
    is unchanged. After the fix the same map has no `InvalidTopology`: every
    row evaluates through `k` = 41 as before (the pre-existing refusals of
    the standing rod with the larger radius, `Unrepresentable` and
    `ComputationLimit`, unchanged) and is the near node from 42 to 52 in
    every frame; two rods standing on the nodes' plane (stored origins
    there, `L^2 = 3`) still evaluate validly to 47. The boolean fuzz
    target's radii (eighths) still cannot reach the band. Regression:
    `tests/curved_booleans.rs`,
    `crossing_cylinders_near_a_node_are_refused_or_valid` (`k` = 42 to 52,
    the three frames, rods, the slab and the standing slab, either radius
    the larger, every operation the near node; at `k` = 40 the exact frames'
    swapped rods, the oblique rods and the exact standing slab evaluating
    with the pair identities); before the fix it fails at `k` = 42. Other
    quadric pairs' near nodes (`spheres.rs`, `cones*.rs`, `torus*.rs`) keep
    `near_node`'s margin alone, unexamined here. Checks: fmt, clippy, the
    1.85 check, the release suite (655 tests); every comparison of
    `HANDOFF.md`'s table with its matches and reviews unchanged, none
    failing, and `compare_split.py` 72/56; replays with debug assertions of
    the boolean corpus and its regressions (1,470 inputs) and the split
    corpus and its regressions (3,571), none failing;
    `tests/curved_booleans.rs` under the emulated glibc `hypot` with debug
    assertions and overflow checks, the new test 6.6 s.
    Campaigns at `298fcf3a` (600 s each, sampled replays, every switch on):
    `boolean` clean, 929 runs, the slowest input 30 s under
    AddressSanitizer; `split` clean, 1,353 runs, the slowest 11 s. DRAW
    survey: that of S9e.4b.4b.2a, S9e.4b.4b.2b.1 and the cylinder pairs'
    fixes at `826346b7` (the same 9 cases the near node, none moving with
    the rule).
  * **Where S9 stands (2026-09-30, paused).** Done and pushed: S9a to S9d
    (every sub-step with its DRAW survey and a clean campaign), S9e.1
    (campaign clean at `51c08edf`) and S9e.2 (`8e060c67`), S9f's decisions
    and S9f.1's evidence (38 cases captured before its kernel code).
    Since the pause (2026-10-01): S9e.2's campaign clean at `6c221525`
    and the DRAW survey of S9e.1 and S9e.2 (no case changes but the
    rollex, already registered), so S9e.1 and S9e.2 are done; `HANDOFF.md`
    summarizes this state. Open, in order: S9e.3 (results with spheres,
    cones or tori, procedural edges, deeper
    chains, `G9` and `H3`): its evidence first, a local WIP branch
    `s9e3-wip` began kernel code before any evidence and is to be restarted
    evidence first; S9f.1's kernel (spline walls in the curved engine and
    the R4 lifting fix found by its evidence: a C1 knot of multiplicity
    `p` in a turned frame fails at the extrusion, `edge_not_c1`), a local
    WIP branch `s9f1-kernel` unfinished and unverified; then S9e.4, S9f.2a,
    S9f.2b, S9f.3; then S9's acceptance (U6: CI green at the revision, the
    schedule run's sharded full replays, a clean campaign).

### Parallel tracks

* **The boolean target's full replay under AddressSanitizer (S9
  acceptance), done: the replay split across schedule jobs.** The local
  corpus's full replay (1,419 inputs) exceeds the 3,600 s startup hour
  under the sanitizer (S9d.3a's note). CI's corpus is the cached one, 356
  inputs (seeds, regressions and the few inputs campaigns kept): its only
  scheduled full replay (`fc695afd`) stopped on crash `6d1fd061` after
  1,667 s (a regression since), and per-push runs, which replay it whole
  without a manifest, took 1,242 s of startup at `85104dc3`, 1,542 s at
  `b96730a7`, 2,074 s at `f4ea7a2a` and 2,108 s at `2efa1ec7` (a 275 s
  build; the slowest input `d0a3de29`, 42 s, within the 60 s limit), the
  same inputs slower with each curved family. The schedule and manual
  campaigns now replay `REPLAY_SHARDS` targets (`boolean`, four shards) in
  jobs of their own (`FUZZING.md`): one published snapshot of the corpus,
  each input replayed once with `-runs=0` in the shard its contents hash
  to, under the target's limits and a startup hour per shard, and a check
  job that passes only if the shards' union is the snapshot, each input
  once, every shard green, and then writes the manifest. The campaign job
  of the same run replays regressions and a seeded sample before its
  mutation. Replayed locally in the four shards at `2efa1ec7` over CI's
  356 inputs: 80, 96, 96 and 84 inputs in 106, 216, 196 and 121 s (639 s
  against CI's 1,830 s in one process), the slowest 15 s, the check clean
  and the manifest written. What remains: the next schedule
  run's green shards and check (acceptance evidence); the Sunday
  minimisation still merges the whole corpus in one process under the
  hour (CI's 356 inputs fit, the local 1,419 did not); the sampled
  campaign keeps no new inputs, so CI's boolean corpus grows by seeds and
  regressions, as it did in practice before. `degree_elevation` overran
  the hour once on the schedule (471 inputs, `7bf4e3a1`; 2,666 s the next
  day), a candidate for `REPLAY_SHARDS` if it recurs.

  The first scheduled run after the split (`36716623883`, at `428349e8`)
  replayed the boolean corpus completely in its four shards, and the check
  job wrote its manifest; `degree_elevation`'s full replay reached 256 of
  its 485 inputs in the hour there, so it is sharded the same way now.

  `split` is sharded too since its full replay overran the hour on the
  schedule runs at `c3ce4d42` (CI's 1,749 inputs; "CI's scheduled full
  replay at `c3ce4d42`" under S9): eight shards, the replay job's matrix
  listing shards up to the largest count and excluding each target's
  past its own (`test_fuzz_runner.py` holds the matrix less its
  exclusions equal to `REPLAY_SHARDS`). Measured locally, a shard of CI's
  corpus takes about 650 to 1,100 s on CI and one of the local corpus's
  size (3,551 inputs) about 1,500 s. As for `boolean`, its scheduled
  campaign now keeps no new inputs, and the Sunday minimisation still
  merges the whole corpus in one process under the hour, which CI's
  1,749 inputs would exceed (about 7,100 s on CI): open.
* **Certified integrals along procedural meetings (S9d.4b.2b), done.** The
  validator's and mass's certified integrals along `Curve3::Toric` meetings
  on two tori took tens of seconds per Boolean under AddressSanitizer, and
  the fuzz target's torus pairs were off (`TORUS_PAIRS`). Profiles (macOS
  `sample` of the release replay): a torus band in the tilted frame against
  a prism spent 64% in the integrals (the torus face's fourteen moments
  along its projections 41%, each evaluation a hundred and more series
  products; the validator's two runs, the result's and its moved copy's,
  24%, two thirds of it the projections' own jets) and 25% in the exact
  arrangement; a whole torus in a turned frame against a prism with arcs
  67% in the arrangement (fields of degree eight: rational gcds, signs and
  binary64 views of `Q(alpha)`) and 30% in the integrals; two tori 79% in
  the integrals, on pieces of `1/64` to `1/2048` of an edge: over a piece
  the interval jets of a meeting's angles blew up (the thirteenth
  coefficient `1e13` to `1e44` over a piece against `1e-8` to `1e24` at its
  middle), the interval Newton for its other angle failed over `1/32` of an
  edge (the natural extension of `G_s`, `-14.6` at the middle, `[-40,
  10.8]` over it), and the edges' ends lie near their meetings' turning
  points. Changes (`MATHEMATICS.md`, certified quadrature and the torus
  meetings' jets; `BOOLEAN.md`, S9c.1): (a) a binary64 series product's
  coefficients summed with one error bound rather than an exactly signed
  rounding per operation; (b) the projections' jets kept per piece by the
  projection's content, shared by the areas, the fluxes and the moments; (c)
  sign decisions' integrals at the width `1e-6` first, `1e-12` only when
  that leaves them undecided; (d) a sphere's or torus's fourteen moments as
  scalar combinations of the integrals of a trigonometric basis, one set of
  jets for all; (e) `G_s` (Newton's divisor and the recurrence's) and the
  first coefficient's `G_f` also in mean-value forms about the base's
  middle, the narrower kept; (f) the arrangement shared by fuse, cut and
  common of one pair (kept by the inputs' content; the operation only keeps
  pieces), each vertex's binary64 view computed once for the meeting's
  pieces, signs across two fields and at a generator's isolator tried by
  binary64 enclosures from the rationals' leading bits first, and a field
  number's binary64 view (`K::value_near`) by an integer Horner reduced once.
  Tried and dropped: other orders (6 to 24) and each piece's Taylor degree
  chosen by its least remainder (no fewer evaluations). Results (release
  with debug assertions, instructions retired as the load-free measure, a
  host at load 10 to 20; AddressSanitizer's targets run side by side):
  the band against the prism (two corpus inputs) 59.1 to 18.2 G, 5.3 to
  1.6 s, 63 to 17 s under AddressSanitizer; the turned torus 179 to 50 G,
  14.4 to 5.6 s, 213 (63 on a quieter host) to 30 s, still 66% the
  arrangement; 44 intersecting two-tori variants of 60 corpus inputs 30.7
  to 12.0 G at the median, 242 to 53 G at the ninth decile and 593 to 71 G
  at the slowest (0.6, 2.75 and 3.7 s), ten of them under AddressSanitizer
  3 to 51 s (26 to 306 s before on the same host). The torus pairs stay
  off: within the target's 60 s here, but at the Linux runners' 2.6-fold
  the slowest tenth would not be, so CI's fuzz runs would time out on them
  now and then; the degree-eight arrangement's arithmetic (below) is
  next. Enclosures: every
  curved Boolean fixture's volume and area within `6.1e-11` relative
  (`1.7e-11` before), medians 1.0 to 4.2 times as wide per set, at most 24
  times, none apart from its former one; all fixture tests, every
  comparison (counts unchanged, 0 failures), `compare_brep.py` (and its
  spline and sheet families) and `compare_brep_io.py` unchanged. The
  corpus and the regressions (1,435 inputs) replay with debug assertions
  without a failure, 1,080 s to 732 s in all and the slowest 21.4 s to 7.1
  s (three processes each, side by side). The band's input is kept as a
  regression (`fuzz/regressions/README.md`); the turned torus's is not (30
  s under AddressSanitizer, over the Linux runners' 60 s at their
  2.6-fold). Open: the arrangement's exact arithmetic for fields of degree
  eight (a rational interval Horner per binary64 view, `Qd::to_f64`), and
  the moved result's second validation.
* **Certified integrals beside spline walls' loops (S9f.3a, S9f.3b), done:
  `SPLINE_SPHERE` and `SPLINE_CONE` on.** The fuzz target's spline prisms
  against its sphere and cone tools were off for time: the corpus's own
  `d7599dbe` took 472 s under AddressSanitizer, a sphere variant
  (`e2fbf381`) 297 s, three cone variants 97, 76 and 49 s. Profiles (macOS
  `sample` of the replay with debug assertions and of the sanitizer's
  target; instructions retired as the load-free measure): `d7599dbe`, 472
  G, spent 83% in its results' mass. 70% was the quadrature of the spline
  walls' four `|N|` integrals over the meetings' pieces
  (`quadrature::spline_face`'s sweeps: 93 tensor-rule integrals examined
  9,747 boxes, 4,594 of them with a series undefined over the box, since
  near a turning point the discriminant `d` is far below its terms, and its
  factors' products over a range, combined, overestimated it by their own
  size, so its square root failed until a piece was a hundred times
  shorter than it needed; such a box was halved across its less divided
  side, nearly half the time `σ`, along which the integrand is linear),
  and a quarter of those sweeps was rational arithmetic made again at every
  evaluation (the pieces' affine maps, their brackets and gcds). 27% was
  the jets along the meetings' pieces (`jet::integrate_many` through
  `projection::wall_pieces`: the wall's ten Green moment tensors, the
  sphere's moments, the validator's areas and fluxes), each integral
  computing the pieces' wall jets again (72% of them repeats). Under
  AddressSanitizer the same integrals cost most in allocations (jets of
  order zero at every de Casteljau step of every Taylor coefficient, a jet
  per term of a cone's moment integrands, de Boor's points at every node of
  a graph over `v`), in the comparison tracing libFuzzer's instrumentation
  adds to every integer test and to every shadow check, and in rational
  arithmetic outside the integrals (the validator's exact spline points;
  the strips, where Green's exact route refused a wall for a segment along
  a ruling an ulp past its domain). Changes (`MATHEMATICS.md`, the
  certified evaluation of `WallMeet` and certified quadrature's (e) and
  (f)): (a) every polynomial of a spline wall's meeting taken over a range
  in centred form (`wall_meet::Centre`: the Taylor coefficients at the
  range's binary64 middle, enclosed at that point, shifted over the range
  exactly; a cone's tiers combined at the point first), `a`, `b`, `c` and
  `d` from the polynomials over a range too, the Taylor coefficients about
  a point in the tier (not as jets of order zero), and a graph over `v`'s
  rows expanded once about the root's base and solved for on series cut to
  `k + 1` terms; (b) the wall pieces' jets kept in the projections' memo by
  the knot span too; (c) constants made once instead of at every
  evaluation: a wall piece's affine maps (`quadrature::Placed`), the Green
  tensors' patch maps, the binary64 Gauss–Legendre rule, de Boor's points
  on the stack; (d) a bisected piece's halves' remainders bounded by its
  own next coefficients where they suffice (`jet::integrate_many`, and the
  tensor rule's boxes likewise); (e) the binary64 tier's rounding errors
  kept as binary64 values and its products' corners decided by
  floating-point tests (the same bounds bit for bit, no integer comparison
  to trace); (f) Green's exact route integrates nothing along a constant-`u`
  piece, a wall's moment tensors take one product of jets per column, and a
  cone's or cylinder's moment integrands are summed in place. Tried and
  dropped: the tiers centred one by one (the cone variants slower: their
  terms no longer cancelled, 710 boxes against 250) and the polynomials
  over a range by de Casteljau's range form (over ten minutes on a cone
  variant). Results, instructions retired with debug assertions
  (`d7599dbe`, `e2fbf381` and `ab0e5716` spheres, `2dddfe21` and
  `bb24e50b` cones): 472, 200, 107, 87.7 and 74.4 G before; 58, 47, 36, 87
  and 48 with (a) but the tiers and (b); 58, 47, 36, 42 and 43 with the
  tiers combined; 52, 41, 34, 41 and 40 with (c)'s maps and rule; 46, 36,
  28, 36 and 34 with the tensors' columns and the point coefficients in
  the tier; 43, 33, 22, 33 and
  31 with (e), the graph over `v`'s expansions, de Boor's points and the
  constant-`u` pieces; 37, 30, 21, 29 and 27 with (d); 36.8, 29.2, 21.1,
  27.1 and 25.5 with the cones' integrands: 12.8, 6.9, 5.1, 3.2 and 2.9
  times fewer. The 879 variants of every third corpus input (S9f.3a's 477
  spheres, S9f.3b's 402 cones), none failing: spheres 5.5 to 4.9 G at the
  median, 14.2 to 11.1 G at the ninth decile and 472 to 36.8 G at the
  heaviest, cones 3.5 to 3.3, 11.7 to 8.6 and 87.7 to 27.1 G; none
  heavier than before. Under AddressSanitizer with the target's own
  options (the allocator feature, a 64 MB quarantine, five-frame stacks)
  on a host at load 8 to 22, which made the corpus's slowest input
  `e36969f1`, unrelated to the switches, take 87 to 138 s (45 s in the last
  campaign at load 3 to 7): `d7599dbe` 24 to 38 s (472 s before),
  `2dddfe21` 23 to 33 s (97 s), `bb24e50b` 21 s (76 s), the 25 heaviest
  variants 16 to 33 s, and the six heaviest of the corpus's 215 inputs
  reaching either switch 19 to 38 s (`d7599dbe` and `85659e3a` 422 G, then
  four of 337 to 391 G, against `e36969f1`'s 1,100 G). Both switches are
  on. Enclosures: every fixture test passes within its `1e-9` enclosures;
  the release suite (548 tests), fmt, clippy and the 1.85 check pass; every
  spline comparison unchanged with 0 failures (`compare_spline_sphere_boolean.py`
  4/29, `compare_spline_cone_boolean.py` 4/23,
  `compare_spline_crossing_boolean.py` 14/37,
  `compare_spline_parallel_boolean.py` 31/12, `compare_spline_any_boolean.py`
  22/16, `compare_boolean.py --splines` 33/13). The corpus (1,432 inputs)
  and the 26 regressions replay with debug assertions and both switches
  on without a failure, the slowest 19 s on a host at load 20 to 27 (91 G,
  inputs unrelated to the switches). Campaign: the boolean campaign at
  `cfeab65d` with both switches on (600 s, a sampled replay) clean, 803
  runs, the slowest input 52 s under AddressSanitizer at load 5 to 8
  (`e36969f1`, an existing corpus input reaching neither switch), the next
  28 s (`a46b4614`). Open: the Linux runners' 2.6 times puts the heaviest of
  these inputs near the target's 60 s, as the corpus's slowest inputs
  already are; the moved result's second validation and the validator's
  exact spline points are the next levers.
* **The degree-eight arrangement's arithmetic (S9e.3b, S9d.4b.2), done:
  `GIVEN_MET` on.** The fuzz target's chained partners centred on a given
  result's meeting of two curved faces were off for time: the corpus's
  slowest chained operations reaching one took 60 to 71 s an input under
  AddressSanitizer, 100 to 230 s on a host at load 7 to 31 (`41833cc4`,
  `487e8cac`, `972c8678`, `0fbd6914`, `442d4917`: 122.8, 112.5, 106.8, 97.2
  and 92.1 G instructions with debug assertions, 12.3, 23.8, 13.2, 4.7 and
  7.4 G with the switch off). Profiles (macOS `sample` of the replay with
  debug assertions and of the sanitizer's target; instructions retired as
  the load-free measure): `0fbd6914` spent 32% in the exact signs of a
  field's numbers (`Gen::sign_of`), 1,126 Sturm-Tarski queries every one of
  them nonzero: the binary64 filter's interval held zero where an element's
  coefficients of 600 to 1,500 bits cancel, and some isolators were still
  `1e66` wide after their 96 bisections (an eliminant's root isolated over
  its root bound); 29% of it hashing the kept signs' rational keys
  (`num_rational` hashes a rational by its continued fraction, a division
  per partial quotient); 29% a torus meeting (`torus_far`, S9d.4b.2's
  discriminants of degree 24 and their roots). `41833cc4`: the signs 39%,
  binary64 views of field numbers 13% (`Qd::to_f64`: a rational interval
  Horner, then the value at isolators bisected up to 1,280 times).
  `487e8cac`: a cone carrier's place on its meeting 49% (an inverse in the
  field by Euclid's algorithm over the rationals at every test of a vertex),
  the rationals' gcds 53% of all. `972c8678`: the arrangement's direction
  tests (`cross2`, the exact products of a dot product in the field) 32%.
  Under AddressSanitizer the same, and `num_rational`'s reductions by
  Stein's binary gcd, whose every comparison libFuzzer's instrumentation
  traces, 19%; and libFuzzer runs again any input whose run ends with more
  allocations than frees (the kernel's kept arrangements and meetings), so a
  wall time is about two runs. Changes, every decision exact as before: (a)
  each generator's root as a dyadic point `x / 2^b` within `2^-b` of it
  (Newton's iteration in integers from the isolator narrowed to `2^-24`,
  certified by the polynomial's opposite signs at `(x -+ 1) / 2^b` inside
  the isolator), and an element's sign from its exact value at that point
  with the absolute coefficients' derivative bounding the rest, at three
  precisions from its size, before Sturm-Tarski (a zero still goes the exact
  way); the binary64 filters' interval about that point, or the isolator
  narrowed by its width's bits too; the kept signs keyed by integers; (b)
  inverses kept per generator, and found by the extended Euclid of integer
  pseudo-divisions; (c) binary64 views rounded to nearest from dyadic
  enclosures at 96, 384 and 1,536 bits after the point (both ends rounding
  alike), the isolator's view where they do not (a rational root, zero); (d)
  a field's numbers as integer numerators over one denominator, a product
  reduced by its content once (the old way a gcd per coefficient and
  operation), a sum not at all (compared by value, keyed by its one form),
  in one vector so an element is no larger beside a rational; (e) a dot
  product's sign in a field (the arrangement's `cross2`, `dot2`), a surd's
  and a tower's, from the dyadic enclosures before the exact products; (f) a
  given meeting's fibre of two quadratics by their combination without the
  leading terms and the point homogenized at its root (one inverse where
  Euclid's algorithm over the field took four), and a cone carrier's `on`
  test on the direction its radius scales (no inverse); (g) rational
  arithmetic in lowest terms by the crate's Lehmer gcd where
  `num_rational`'s operators reduced by Stein's (S9d.4b.2's charts, forms,
  polynomial sums and remainders, the nearest binary64 of a rational, an
  isolator's ends after bisection), products of polynomials in integers over
  their denominators, and a quartic's discriminant in integers over one
  common denominator (a positive multiple: its roots); (h) coprime integer
  polynomials certified modulo `2^61 - 1` before the subresultant chain, and
  a polynomial's gcd with its derivative where it is a chart's `(1 + t^2)^m`
  times a part certainly square-free that way (`(1 + t^2)^(m - 1)`, the
  chain's). Every sign, inverse, gcd and point is the one the exact path
  gives; a binary64 view is now the number rounded to nearest (before, the
  middle of a narrow interval, or a value at a bisected isolator: one ulp
  off in 30 of 1,458 views on `972c8678`, each checked against the isolator
  narrowed by 400 bisections more). Tried and dropped: integer products of a
  one-term operand or of denominators far apart (slower: term by term
  there), the isolators' bisections jumped by Newton's iteration to the same
  intervals (1 to 3% fewer instructions, not worth its code). Results,
  instructions retired with debug assertions (`41833cc4`, `487e8cac`,
  `972c8678`, `0fbd6914`, `442d4917`, and `e36969f1`, the corpus's slowest,
  which the switch reaches too): 122.8, 112.5, 106.8, 97.2, 92.1 and 112.1 G
  before; 69.2, 109.6, 73.2, 60.9, 58.7 and 109.7 with (a)'s signs; 66.3,
  65.8, 63.3, 60.8, 58.7 and 109.6 with the inverses kept; 50.1, 58.9, 52.1,
  59.2, 56.6 and 102.3 with (c); 48.2, 60.7, 51.0, 51.5, 49.1 and 88.5 with
  the integer products and (f); 42.9, 51.4, 43.5, 45.2, 42.5 and 82.5 with
  (d)'s integer form; 38.6, 51.0, 40.3, 36.1, 34.3 and 56.2 with (g); 37.0,
  47.1, 36.2, 33.4, 31.9 and 55.2 with (e) and the integer inverses; 35.8,
  47.0, 35.2, 31.9, 31.1 and 53.8 with (a)'s Newton start; 34.0, 47.0, 33.3,
  27.4, 26.7 and 47.7 with (h); 30.0, 42.0, 29.6, 24.0, 24.0 and 43.6 with
  the sums unreduced: 4.1, 2.7, 3.6, 4.1, 3.8 and 2.6 times fewer. The
  corpus (1,432 inputs) and the 26 regressions replay with debug assertions
  and the switch on without a failure, 9.71 T instructions to 7.22 T in all,
  the slowest 43.6 G (`e36969f1`; 112.1 G before, 90.8 G with the switch
  off). Under AddressSanitizer with the target's own options (wall and CPU,
  user plus system, time), on a host at load 25 to 50 (other agents' builds
  and replays): the inputs the switch reaches 33 to 62 s wall and 28 to 40 s
  CPU (`41833cc4`, `487e8cac`, `1745232e` 62, 51 and 59 s wall, 32, 34 and
  40 s CPU; 186, 142 and 120 s wall, 158, 125 and 105 s CPU before),
  `e36969f1` 87 and 54 s (225 and 192 s before); at an intermediate commit
  (`32297b78`, 10 to 20% heavier than the last) on a host at load 6 to 7, 21
  to 35 s wall, `e36969f1` 50 s (52 s in the last campaign, with the switch
  off). Each of these is two runs: libFuzzer runs an input again, its leak
  detection on, when a run ends with more allocations than frees (the
  kernel's kept arrangements and meetings): `1745232e` 20 s CPU with
  `-detect_leaks=0` against 40 s; the target's per-run limit sees one. The
  "points not separated by a projection" limits: the corpus's 8 were two
  cylinders with parallel axes, meeting in two rulings a sphere crosses at
  equal heights (each point shares its ruling with one and its height with
  another, so no ruling's angle or height separates them); `triple.rs` now
  projects along each ruled surface's place sheared by its angle's chart
  (the place plus `c` times the chart's variable, `c` = 1, -1, 2 and 1/2),
  after the projections before: 4 of their chained operations evaluate, 4
  are refused as degenerate (thin pieces, two meetings within rounding),
  none as a limit (`triple.rs`'s test of the four points `(3, +-4, +-3)`).
  Checks: fmt, clippy, the 1.85 check, the release suite (555 tests), the
  fuzz crate's fmt and check; every comparison unchanged with 0 failures
  (`compare_given_met_boolean.py` 8/42, `compare_given_curved_boolean.py`
  25/23, `compare_given_boolean.py` 36/0, `compare_chained_boolean.py` 24/6,
  `compare_torus_boolean.py` 11/24, `compare_torus_segment_boolean.py`
  15/14, `compare_torus_curved_boolean.py` 15/29,
  `compare_torus_parts_boolean.py` 16/37, `compare_cones_loops_boolean.py`
  5/26, `compare_spheres_turned_boolean.py` 0/18). Open: a campaign with the
  switch on; the leak-check rerun (clearing the kernel's kept arrangements
  at an input's end would halve the sanitizer's wall time on curved inputs);
  the certified integrals along a cone carrier's meetings (`487e8cac`: 58%
  of its instructions in `projection::along`); the torus meetings' root
  sampling (`e36969f1`'s `all_roots`, 72 bisections a root at each of 17 to
  49 samples).

  Campaign: the boolean campaign at `a32d256f` with `GIVEN_MET`,
  `SPLINE_SPHERE` and `SPLINE_CONE` on, after this track, the slowest
  inputs' and the meridian loop's fix (600 s, a sampled replay) clean, 895
  runs, the slowest input 31 s under AddressSanitizer at load 5 to 9
  (`c4d14627`, a polyhedral chained cut; `e36969f1` no longer among the five
  slowest).
* **The boolean target's slowest inputs (S9's campaign), done.** The 600 s
  boolean campaigns were clean but the corpus's slowest input, `e36969f1`,
  took 45, 52 and 60 s under AddressSanitizer at the target's 60 s limit
  (load 3 to 8), and the Linux runners take about 2.6 times this host's;
  next `62cf4d8f`, `a46b4614`, `17e131e3` and `c4d14627`. Profiles (macOS
  `sample` of the replay with debug assertions and of the sanitizer's
  target; instructions retired as the load-free measure; the corpus decoded
  stage by stage): `e36969f1` (90.9 G, `a67d82c5` the same Boolean) is a
  stadium prism against the outer half of a torus in the tilted frame, its
  two arcs' cylinders met by the torus (S9d.4b.2a's `torus_far`, 71% of its
  instructions); 54% was the S9e.4a stage, whose imported construction
  (its cylinders' centres read back, with more bits) arranges the same
  pair again, a third heavier. Within a torus meeting: the windows of its
  pieces 42% (`window`'s 17 samples a run, 136 a meeting, each root's angle
  from its isolator bisected 72 times, then turned through the chart and
  rounded: `Chart::at`, `rational_f64`, all in rationals reduced at every
  operation), the critical values 24% (the quartic's discriminant of
  degree 24, a product of rational polynomials, 16%, and its square-free
  roots), the pieces' verification 12% (`clear_over`'s Sturm chains in
  rationals); the rationals' gcds a third of everything. `62cf4d8f`,
  `17e131e3`, `f0987e7a`, `1745232e` and `01b580b7`: prisms against torus
  bands, wedges and halves in the tilted or turned frames, the same meetings
  and the arrangement's tests in fields (`cross2`, `Qd::sign`, `Qd::to_f64`:
  the degree-eight arrangement's arithmetic above); `c4d14627`, a polyhedral
  chained cut, 16 s with the sanitizer against 1.9 s without (its rational
  vectors' allocations and gcds), unchanged. `a46b4614` is not in this
  corpus. Under the sanitizer each of these runs twice in a replay (the
  leak check's rerun, above); the per-input limit sees one run, so the
  timings here are one run (`-detect_leaks=0`). Changes, every number the
  same (each result's text, its history's and its mass properties' hashed
  over the corpus and the regressions, identical before and after each
  step): (a) `rational_f64` by one integer division in the comfortably
  normal range (the quotient's 55 or 56 leading bits and the remainder: the
  nearest binary64, a tie the lower, as the bounds' comparison gives it);
  (b) `Chart::at`, `Form::poly` and a torus meeting's form at a rational
  angle (`Bi::at`) in integers over one denominator, each coefficient
  reduced once; (c) the quartic's discriminant over its coefficients'
  common denominator (with the arithmetic above, a positive multiple); (d)
  `roots_repeated`'s square-free part by an exact integer quotient (Gauss's
  lemma: both primitive), the gcd as above (modulo a prime, a chart's
  factor); (e) `clear_over`'s Sturm chain in integers at its rational ends
  (positively scaled pseudo-remainders, contents removed: each member a
  positive multiple, the same sign changes); (f) what the arrangement asks
  again kept: a half-edge's travel at its start across the faces traced at
  a vertex, a vertex's binary64 view and an edge's binary64 samples (each
  face along it and each operation's assembly asked again); (g) the mass
  moments' integrals along a projection kept by their content (its and its
  integrands' `Debug` texts): a pair's other results and an imported
  input's own construction ask one face's integrals along one pcurve again.
  Results, instructions retired with debug assertions: with `GIVEN_MET` off
  (as before the arithmetic above), the corpus and regressions 8,145 G to
  6,882 G with (a) to (e) but `Bi::at`, and 6,569 G with it, (f) and a
  modular gcd (since replaced by the arithmetic's own);
  `e36969f1` 90.9 to 50.6 and 46.9 G, `62cf4d8f` 91.1 to 76.5 and 68.7,
  `f0987e7a` 53.5 to 34.0 and 32.2, `17e131e3` 52.0 to 43.5 and 41.2. Merged
  with the arithmetic above (`86b1d834`, the switch on: 7,307 G, `e36969f1`
  43.8 G), its own quartic and gcd kept: 6,418 G, and 6,312 G with (g);
  `e36969f1` 36.4 G, `487e8cac` 42.1 to 39.3, `1b405929` 31.2 to 27.3, the
  slowest now `487e8cac` (its results' mass along a cone carrier's
  meetings). Under AddressSanitizer with the target's options, one run each,
  the original target and the merged one side by side on a host at load 10
  to 15 (before the arithmetic above and these changes; after both, the
  switch on): `e36969f1` 67.2 s and 1,051 G to 21.0 s and 275 G, `62cf4d8f`
  64.5 s to 17.7 s, `17e131e3` 36.9 to 17.8, `f0987e7a` 72.3 to 10.6,
  `1745232e` 33.2 to 17.1, `85659e3a` 19.9 to 18.8; `487e8cac`, `41833cc4`
  and `972c8678`, which the switch reaches now, 17.6, 15.0 and 16.5 s
  (10.1, 7.4 and 8.5 s off); every one under 23 s, and with (g) at the
  last commit `e36969f1` 20.4 s (267 G), `487e8cac` 16.5 s; `c4d14627`
  (251 G either way) 16 s at load 8, 26 s at load 25. Checks: fmt, clippy
  (release, all targets), the 1.85 check, the release suite (569 tests),
  the fuzz crate's fmt and check; every comparison unchanged with 0
  failures (the table's 31 boolean comparisons, `compare_split.py` 72/56,
  `compare_brep.py --family spline` 10/3, `compare_brep_io.py` 6,835/7);
  the corpus and regressions replay with debug assertions (the slowest
  39.3 G, `487e8cac`). Open:
  the root sampling's 72 bisections a root (now 15% of `e36969f1`, in
  `refine_for_signs`, which must give the same isolator); the leak-check
  rerun; `487e8cac`'s certified integrals along the cone carrier's meetings.
  The corpus's `21928984` panics with `InvalidTopology` (`uv_gap`) at
  `86b1d834` and here alike (`GIVEN_MET` on, S9e.4b.3a's stage and the
  validator together; fixed on its own track), so the replays above count
  1,457 inputs.

  Campaign: the boolean campaign at `a32d256f` with `GIVEN_MET`,
  `SPLINE_SPHERE` and `SPLINE_CONE` on, after this track, the slowest
  inputs' and the meridian loop's fix (600 s, a sampled replay) clean, 895
  runs, the slowest input 31 s under AddressSanitizer at load 5 to 9
  (`c4d14627`, a polyhedral chained cut; `e36969f1` no longer among the five
  slowest).
* **The boolean target's off switches (S9d.3c, S9d.4b.2b, S9d.4c), done:
  `CONE_PAIRS` on; `TORUS_PAIRS` and `TURNED_PARTS` stay off.** The three
  switches were off for time (their notes above and S9d's): two tori's
  integrals, two cones' integrals, and a cap's or part's crossings with a
  torus in the tilted and turned frames. Each was measured on: the corpus
  and its regressions replayed with debug assertions and the switch on in
  the source, and variants that reach it, every corpus input moved onto the
  switch's configuration (the object's bits and the tool's byte; three
  thirds of the corpus, 1,203 distinct cone-pair, 1,203 two-tori and 752
  turned-part variants), instructions retired as the load-free measure;
  the slowest under AddressSanitizer with the target's options, one run
  each (`-detect_leaks=0`). Failures, all `CONE_PAIRS`'s, each a kernel bug
  fixed with a test and a regression input (`fuzz/regressions/README.md`):
  (1) corpus input `6fab9d41` and 28 of the first 459 cone-pair variants
  failed `vertex_off_curve` in the chained stage, and one more was refused
  as a plane at the wrong cone's apex: a given edge of one cone's wall and
  the other cone's end disc has both cones among its faces' views, and its
  section was rounded on the first with a funnel, the disc's, whose frame
  turned the section's plane into another (a hyperbola off the vertices);
  it is rounded on the cone it lies on now (`ConeSec::lies_on`, `be38e1c2`).
  (2) Two coaxial cones whose ring lies `0.004` below the frustum's top rim
  failed their fuse (`an open Boolean of arcs in any position`): a wall's
  pieces are nested by their loops' binary64 polygons, 24 sides a circle,
  and the rim's and the ring's crossed; a hole's point near another loop's
  polygon now has the loops sampled again, four times as finely, to choose
  among the pieces that may hold it (a hole none holds is refused as
  before: finer polygons nested one in a `TURNED_PARTS` variant and its
  result was open; `c82af2ec`, `ffd51b85`). (3) Two coaxial cones' ring
  within rounding of the turned box's floor, given to the box, crossed its
  walls at the floor's edges and the floor not at all (`a cone's wall
  winding without an apex`): a given closed meeting all round within the
  resolution of a partner's plane is refused as a tangency now
  (`triple::meet`). Profiles (macOS `sample` of the replay with debug
  assertions): the slowest turned part, `893d996d` (a cap against a wedge in
  the tilted frame, 329 G), spent 80% isolating the roots of its circles'
  crossings with the torus (`torus_parts::circ_torus`): the Sturm chains'
  exact Horner at every bisection's midpoint (two thirds of all in bigint
  products) and the midpoints' rational gcds; then the fields' products and
  inverses, one inverse per crossing of one resultant. The slowest cone
  pair, `2f5a979d` (two coaxial cones in the tilted frame, 159 G), spent 95%
  in the certified integrals along their ring: 3,000 pieces a turn, its
  quadratic's coefficients' jets over a piece wide (`cos^2 + sin^2` and the
  discriminant `b^2 - a c` as natural extensions, their terms cancelling;
  the square root's recurrence widened the top coefficient to `1e17`); the
  next, `27c68428` (cones crossing on their side), 2,600 pieces an integral
  where `d`'s interval over a piece was `[6.5, 20.8]`. Two tori, `73b90a76`
  and `86e3875f`: half in the integrals along their `Toric` meetings, whose
  jets' widths grew twentyfold an order even at a point (`+-3.4e4` at order
  13: each step's corrections counted its provisional terms twice), a
  quarter of that the moved result's second validation, the rest the
  degree-eight arrangement. Changes: (a) an integer polynomial's sign at a
  rational by a binary64 Horner with unbounded exponents first (each
  coefficient's top 128 bits rounded once; the value beyond `(6 d + 8) 2u`
  times the terms' absolute sum, Higham's bound, has the exact sign; the
  exact Horner otherwise; `IntPolynomial::sign_filter`); (b) an isolation
  bisected over one denominator, the same midpoints reduced only where a
  root is published; (c) `turned::roots`' gcd given to the isolation, and
  `circ_torus`'s resultant's roots isolated once; (d) what a generator's
  polynomial gives every root alike, `x^j mod poly` and the inverses of the
  elements prime to it, shared by its roots' generators (`num::Shared`);
  (e) a `Meet`'s quadratic as forms of degree two in the angle, their
  constants combined before any jet (the narrower coefficient by coefficient
  beside the plain products), and its discriminant's constant term also in
  its mean-value form about the piece's middle; (f) a `Toric`'s recurrence
  with each step's outputs linear forms in its provisional terms, each
  taken once (`+-23` at order 13 now). Tried and dropped: a transversality
  test at a given meeting's crossings for (3) (the coaxial ring crossed
  the box's walls, not its floor), and the finer polygons for a hole no
  piece holds (2). Results, instructions retired with debug assertions:
  `893d996d` 329 G to 114 G with (a) and 72 G with (b) and (c), 56 G at the
  end; `5467edc4` (a cap against a band) 109 G to 80 G with (d);
  `2f5a979d` 159 G to 14 G and `27c68428` 90 G to 18 G
  with (e); `73b90a76` 80 G to 58 G and `86e3875f` 93 G to 82 G with (d)
  and (f); the corpus's `487e8cac` (a cone carrier's meetings, the last
  track's open lever) 39 G to 20 G. The variants, none failing: cone pairs
  3.3 G at the median, 12.5 G at the ninth decile and 48 G at the slowest;
  two tori 3.8, 25 and 82 G; turned parts 16, 36 and 80 G (20, 52 and 329 G
  before over the first 300). Under AddressSanitizer at load 9 to 14, one run
  each: the slowest cone pairs 28 s (`99d30aca`, 48 G), then 21 to 23 s;
  the slowest two tori 40 and 36 s (`86e3875f`, `864b098c`), then 29 and 27
  s; the slowest turned parts 43, 43, 39 and 38 s; the corpus's slowest
  `d7599dbe` 19 s, `85659e3a` 18 s, `e36969f1` 17.5 s (20 s in the last
  track), `487e8cac` 9.7 s (16.5 s). `CONE_PAIRS` is on; `TORUS_PAIRS` and
  `TURNED_PARTS` stay off, above the 35 s the switches on keep to here.
  Exactness: (a) to (d) give every result, history and refusal of the
  corpus, its regressions and 1,218 variants as before (each result's
  written text, mass properties and history hashed), every isolator the
  same, but for the fixes' own: the 29 failures, a refusal at the wrong
  cone's apex now evaluating, and one refusal's message. (e) and (f) are
  the same functions: of the 9,809 results 1,660 report other volumes or
  areas, within their enclosures
  (3.6e-10 relative at most, a volume of `7e-5`), 2,034 enclosures narrower
  and 252 wider (5.8 times at most), every one within the integrator's
  width. Checks: fmt, clippy (release, all targets), the 1.85 check, the
  fuzz crate's fmt and check, the release suite (578 tests); every
  comparison unchanged with 0 failures (the table's 30 boolean comparisons,
  `compare_imported_pieces_boolean.py` 32/13 among them); the corpus (1,433
  inputs) and the 30 regressions replay with debug assertions and the
  switches as left without a failure, 6.1 T instructions in all, the slowest
  37.8 G (`734bc584`, a regression of this track), and with all three on
  too, as do the 3,158 variants. Open: `TORUS_PAIRS`' slowest, the moved
  result's second validation (a quarter of it) and the toric jets' cost a
  piece (interval Newton and the recurrence, now near their natural piece
  count); `TURNED_PARTS`', the Sturm chains of a circle's three resultants
  (the count's with the torus offset either way) and the fields' products
  over the rounded frames' inverses.
  Campaign: the boolean campaign at `cc7ea7ff` with `CONE_PAIRS` on (600 s,
  a sampled replay) clean, 996 runs, the slowest input 13 s under
  AddressSanitizer at load 8 to 10.
* **The boolean target's last two off switches (S9d.4b.2b, S9d.4c), done:
  `TORUS_PAIRS` and `TURNED_PARTS` on.** The off switches' track left
  them off, their slowest variants 40 and 43 s a run under
  AddressSanitizer. Measured as there: the 1,203 two-tori and 752
  turned-part variants (every corpus input moved onto the switch's
  configuration), instructions retired with debug assertions as the
  load-free measure, and the slowest one run each under AddressSanitizer
  with the target's options (`-detect_leaks=0`), the original target and
  the changed one side by side; every result's written text, mass
  properties and history hashed over the corpus, its 30 regressions and the
  variants with both switches on (3,421 inputs), identical before and
  after each change. Profiles (macOS `sample` of the replay with debug
  assertions and of the sanitizer's target): the slowest turned part,
  `5467edc4` (a cap against a part in a turned frame, 80 G; 595 G under the
  sanitizer), spent a third in its cap circles' crossings with the torus
  (`circ_torus`): each circle met each face on the part's torus, for each
  seam tried, each operation and the imported input again, so eight times
  the same resultants and their Sturm chains (32 calls, 4 different);
  `64c75743` 18% in a cavity's containment rays, whose exact tier
  evaluated a torus section's ends in rational series for a loop with a
  projection's fin, which the boundary test never decides; then the
  fields' signs: a torus's seam sides at every vertex by exact products in
  `Q(alpha)` (`Ring::sides`), two fields' dot products by rational
  interval enclosures asked again and again, and on `b48d4be3` (the
  slowest after the first changes) 244 Sturm-Tarski queries of elements
  with coefficients of 300 to 25,000 bits, every one nonzero: `tight`'s
  Newton iteration never certified a dyadic point from its starts (complex
  roots beside the real one put them outside its convergence), so the
  dyadic filter had none, and it was tried again at every query (4,000
  times). The slowest two tori, `86e3875f` (82 G; 519 G), spent a third in
  the integrals along their `Toric` meetings, each edge's jets computed for
  both its faces' pcurves (they run opposite ways), and a quarter in the
  moved result's second validation; under the sanitizer two thirds in
  their meeting's arrangement, up to half of that (`864b098c`) the root
  sampling's bisections (`refine_for_signs`, 72 a root: past the binary64 filter's
  46 bits every midpoint's sign by the exact Horner on growing numbers,
  and the rational root test and the ends' reductions every 16 steps).
  Changes, every sign, isolator and decision the one before: (a) a
  polynomial's isolated roots kept by its integers (`turned::roots`), and
  a circle's reduced quartics and resultants with a torus by its and the
  torus's numbers, the rotation and the offset (`circ_torus`; the fields'
  generators made afresh, so no two crossings share one that did not);
  (b) a constant gcd's quotient the primitive part, without the rational
  division (`isolate_with_gcd`); (c) a containment ray's hit undecided at
  once where a loop holds a fin the boundary test does not decide (a
  projection, a spline, a sinusoid, an arc on a periodic surface): it was
  never clear, and `face_hits` took either answer alike; (d) a torus
  meeting's jets kept in the binary64 tier by its numbers and the
  variable's base (`toric_jet`), a reversed use's `B - s` given the jets in
  `B + s` with the odd coefficients negated (binary64 intervals round both
  ways alike, so the recurrences give the same bits with every odd term
  negated, zeros' signs aside, which no bound takes: checked on every
  two-tori variant and in a test); (e) an isolator bisected, once the
  binary64 filter leaves a midpoint undecided, by the polynomial taken onto
  the isolator then (`Shifted`: `D^d p((a + y w) / D)`, the same signs, by
  the filter again), and where a small prime shows the defining polynomial
  free of rational roots no rational root test and the ends reduced at the
  end alone (`RootPolynomial::rational_free`; tests against the plain
  bisection); (f) a torus's seam sides from binary64 enclosures of the
  point's local coordinates before the exact products (`Ring::sides`); (g)
  two fields' dot products' enclosure signs kept by their numbers, the
  generators kept alive with their keys (`approx_dot_sign`); (h) `tight`'s
  failed precisions kept, and where it has no point one for signs alone
  from isolators narrowed 192 to 1,536 bisections more (`sign_point`),
  `tight` asked only what it was before, so its points, and the binary64
  views from them, are the ones they were (a first version asked it at
  the sign filter's precisions, found points it had not and moved two
  variants' views: caught by the hashes, kept apart; test). Tried and
  dropped: the projections' memo enlarged (1 to 2%). Not done: the moved
  result's validity carried over. Its geometry is the motion's rounded
  image, not the exact one, so the original's certificates hold of it only
  with each certified integral's margin against the rounding's effect on
  it, which nothing bounds yet; the second validation is cheaper through
  (d) instead (the moved jets' two pcurves share them too). Results,
  instructions retired with debug assertions, before and after: turned
  parts `5467edc4` 80.2 to 43.7 G, `93091c44` 80.1 to 31.7, `64c75743` 71.8
  to 39.6, `fabb98f8` 70.9 to 32.8, `b48d4be3` 66.7 to 34.0, `d32da787` 62.6
  to 43.8 (now the slowest); two tori `86e3875f` 81.8 to 54.6, `f799cd7f`
  64.8 to 45.8, `864b098c` 62.9 to 37.9, `9b1c1192` 59.2 to 41.5. By step on
  the slowest: (a) and (b) the turned parts 80 to 66, 80 to 56 and 71 to 56
  G (the resultants a further 10 to 18%), (c) `64c75743` 66 to 54, (d)
  `86e3875f` 79 to 68 and `f799cd7f` 63 to 52, (e) `864b098c` 58 to 43 and
  `86e3875f` 68 to 60, (f) 60 to 55, (g) `5467edc4` 56 to 49, (h)
  `b48d4be3` 53 to 34. The variants: two tori 3.6 G at the median, 18.9 G
  at the ninth decile and 54.6 G at the slowest (3.8, 25 and 82 G before),
  8.4 T in all (11.0 T); turned parts 10.3, 23.3 and 43.8 G (16, 36 and 80
  G), 9.2 T (14.0 T). Under AddressSanitizer, one run each side by side
  at load 6 to 16: `86e3875f` 30.1 to 19.8 s (519 to 343 G), `864b098c`
  27.0 to 16.1 s, `9b1c1192` 23.8 to 16.1 s, `5467edc4` 53.1 to 33.4 s wall
  at load 14 to 16 (42.6 to 24.4 s of CPU), `93091c44` 37.0 to 12.0 s,
  `d32da787` 26.1 to 18.1 s, `a4843568` 29.6 to 16.6 s; the 40 heaviest
  variants of each kind 8 to 22 s (two tori) and 13 to 22 s (turned
  parts) of CPU at load 7 to 35, where the corpus's slowest inputs took 15
  s (`d7599dbe`, `85659e3a`) and 11 s (`e36969f1`). Both switches are on.
  Checks: fmt, clippy (release, all targets), the 1.85 check, the fuzz
  crate's fmt and check, the release suite (597 tests); every
  comparison unchanged with 0 failures (the table's 31 boolean comparisons,
  `compare_imported_pieces_boolean.py` 32/13 and
  `compare_split_pieces_boolean.py` 26/22 among them, `compare_split.py`
  72/56); the boolean corpus (1,436 inputs) and its 30 regressions replay
  with debug assertions and the switches on without a failure, 6.04 T
  instructions in all (6.33 T before with both on), the slowest 35.1 G
  (`85659e3a`, `d7599dbe`; 46.9 G before, `ac19cff5`, a turned part), and
  the split corpus (3,550) and its 19 regressions too; the torus test
  files (31 tests) and the library's unit tests pass under Linux's
  correctly rounded `hypot` emulated. Open: a campaign with both switches
  on; the arrangement's direction tests at a vertex (`graph::next_on_kept`
  and `model::angle_cmp`: the exact coordinates `qqdot` and `qcross`
  build before any enclosure, a tenth to a sixth of these inputs under
  the sanitizer; S9e.4b.3c's files, left to that track), a torus
  meeting's pieces (`pieces_of`'s windows and exact points), and the
  moved result's validity carried over.
  Campaigns at `59d0c57b`, with every switch of the boolean target on
  (`TORUS_PAIRS`, `CONE_PAIRS`, `TURNED_PARTS`, `GIVEN_MET`,
  `SPLINE_SPHERE`, `SPLINE_CONE`) and the caps' frames kept bit for bit (600
  s each, sampled replays): `boolean` clean, 1,080 runs, the slowest input
  13 s under AddressSanitizer at load about 6; `split` clean, 1,818 runs,
  none slow. DRAW survey: that of S9e.4b.3c.1 and the switches' speed-up
  (S9e.4b.3c's section), at `f4b584f7`: the speed-up moves no status,
  reason or audited value.
* **The DRAW adapter's curved primitives (S9d.4b.2b), done in S9d.4b's
  survey.** It built every `ptorus`, `psphere` and `pcone` with the same
  ids, so two tori (and two spheres or cones) were refused as solids
  sharing ids. Each is built under an operation of its own now, and a copy
  of one sharing ids with the other argument is built again from its
  constructor's numbers in its frame: 20 upstream cases of two cones
  evaluate, the other 41 are the kernel's `Degenerate`, and the bridge
  tests two spheres and a torus with its moved copy.
* **Equal cylinders with crossing axes (S9c.1, found against
  `bfuse_complex/K1`'s part bore and tool rod), done (2026-10-09).** A
  plate's round hole of radius 1.5 along `y` against a rod of the same
  radius along `x` whose axis crosses the hole's at a right angle (exact
  frames, the plate built or a box already cut by a rod along `y`) gave
  three answers: the fuse `InvalidTopology("non_manifold_vertex")` (one
  shell touching itself, which the curved and polyhedral results' builder
  passed on as invalid), the cut one solid and the common `Degenerate`
  (solids touching at a vertex), where a rod of radius 1 tangent to the
  hole at one point is a tangency between the inputs for all three. Why:
  S9c.1's two ellipses of equal circular cylinders with meeting axes cross
  at the two ends of the cylinders' common extent across the axes' common
  perpendicular, where the surfaces are tangent (both normals along it),
  and `curved/graph.rs` took those points as vertices whatever the faces'
  sides. Refusing every such point (S9c.2a's rule for a single node, and
  S9c.2b's for the pair in turned frames) would also refuse the whole
  rods' fuse and common, which do not touch themselves: S9c.1's
  `steinmetz_fuse` and `steinmetz_common`, and the registered upstream
  cases `bfuse_complex/J5`, `bopfuse_simple/ZD8` and `ZE1` and
  `bopcommon_simple/ZD8` and `ZE1`. Decision (the user's, 2026-10-09):
  contact only. A crossing point inside both faces where the faces'
  outward normals are opposite (the inputs touching each other there: a
  bore against a rod of its radius across it) is a tangency between the
  inputs, `Degenerate` for every operation; where they are alike (whole
  rods, overlapping there) it stays a vertex, the fuse and common
  evaluating and the cut touching itself as before (the exact sign of
  the normals' dot product, `normal_at` on both models, given and imported
  ones too). Also, S9a's rule now holds for the curved and polyhedral
  results (`polyhedra.rs`'s `solid`): a result whose validation finds a
  non-manifold vertex is `Degenerate` ("a result touching itself at a
  vertex"), never `InvalidTopology`.
  `tests/curved_booleans.rs`'s `equal_crossing_cylinders_are_a_tangency`:
  the plate built and as a cut box against the crossing rod and the
  tangent rod, each operation either way round, the tangency; the bore as
  a whole rod against the crossing rod, the Steinmetz fuse and common
  (`16 r^3 / 3`) and the cut `Degenerate`. S9c.1's fixtures, its
  comparison (42 matches, 2 reviewed) and the five upstream cases
  unchanged.

* **The Boolean captures' Linux records, open: CI has no observations to
  take them from.** S9e.4b and S9f.3's notes leave each capture's Linux
  record pending until CI's run, but the B-rep job runs only
  `compare_boolean.py`'s default set (S9a.1's two prisms): the
  `source-pinned-brep-results` artifacts of the green runs at `d1869f2f`
  (37438561216, push; 37454728856, schedule) and `64673ebd`
  (37415799498) hold `boolean-oracle/native-observed.txt` alone, its 45
  rows equal to that capture's `platform-linux/native.txt`. No other
  Boolean comparison has ever run in CI, and none of the other 38
  `occt-boolean-*-preimplementation` captures (S9a.2's splines, S9c to
  S9f, plane-parts not yet pushed) has a `platform-linux/` record; their
  comparisons are run on macOS only, where `platform_record` reads the
  capture itself. No record was added. Taking them needs the B-rep job to
  run each set's compare script on its pinned SDK and upload each set's
  output directory (the kernel probe's time counted against the job's
  limit), then one green run's
  observations reviewed against the macOS captures as the other
  `platform-linux/` records were.

* **CI budget (U6).** Per-push fuzz runs replay a bounded sample plus every
  regression and new seed; the daily schedule replays everything; the heavy
  exact targets run on schedule only; Windows runs a smoke subset per push
  and the full suite nightly. The gates do not change: acceptance still
  needs the full replay, on the schedule run of the accepted revision.
* **Higher-order certified quadrature** for spline mass properties, to bring
  F8 from 5% to the enclosure widths of the analytic family.

  F8 quadrature decisions, recorded before its code (2026-09-28):

  * **Method.** Gauss–Legendre with `n = 8` nodes per direction on
    subdivided pieces, with the rule's exact remainder bounded from
    enclosed derivatives. On `[a, b]`, `∫ f = (b - a)/2 Σ w_i f(x_i) + K_n
    (b - a)^(2n+1) f_(2n)(ξ)` with `K_n = (n!)^4 / ((2n + 1) ((2n)!)^2)`
    and `f_(2n) = f^(2n)/(2n)!` the Taylor coefficient. The nodes are
    enclosed by exact sign changes of the Legendre polynomial `P_8` and the
    weights `2/((1 - x^2) P_8'(x)^2)` over those brackets (rational
    intervals, once per thread). `f_(2n)` over the piece comes from
    interval Taylor arithmetic: the recurrences of `+`, `×`, `÷`, `√` and
    `cos`/`sin` applied to the inputs' series at `X + ε`, `X` the piece, are
    inclusion isotone, so they contain the coefficient at every point of
    `X`. On a box, the tensor rule's remainder is `K_n Δx^(2n+1) Δy
    f_(2n,0)(X × Y) + K_n Δx Δy^(2n+1) f_(0,2n)(X × Y)` (the rule in `x`
    applied to `∫ f dy`, then the rule in `y` at each node in `x`, whose
    weights are positive and sum to `Δx`), so only univariate series along
    each axis with the other held as an interval are needed. Why not Taylor
    models (the integral of the expansion at the centre plus a Lagrange
    term): on a piece of half width `h`, for the same sixteenth
    coefficient, theirs is `2 h^17/17` times it and Gauss's `K_8 (2h)^17`,
    about `π h (h/2)^16`, 2,500 times smaller; they need a second series at
    the centre instead of eight point values, and in two dimensions
    bivariate series. The node values are cheap and tight. `n = 8` puts
    the remainder at the sixteenth coefficient: for an integrand analytic
    within `ρ` of the piece, `(h/2ρ)^16` falls below `2^-40` at `h/ρ`
    about `0.35`, a few halvings on the fixtures.
  * **What it integrates.** (a) `∫ M/W^k` of Bernstein polynomials along a
    rational spline pcurve piece (planes, and rational pcurves on
    nonrational patches). (b) `-∫ F(u(τ), v(τ)) u'(τ) dτ` along spline
    pcurves on cylinders, cones, spheres and tori, `F` the closed-form
    antiderivative in `v` (polynomials in `v` and `cos`/`sin` of `u` and of
    multiples of `v`). (c) Every term of a rational spline surface and the
    four `|N|` terms of a nonrational one (the other ten stay exact): by
    Green's theorem with `G` from the domain's start through the patches,
    a boundary piece in the patch `[u0, u1] × [v0, v1]` contributes, in its
    local coordinates, `J = -∫_0^1 ∫_0^1 ū'(τ) v̄(τ) f̄(ū(τ), σ v̄(τ)) dσ
    dτ` plus, for every patch below it in its column, the same `J` of the
    line from `ū(0)` to `ū(1)` at `v̄ = 1`. The mass integrands are
    homogeneous of degree one in `N`, so with `N̄ = S_ū × S_v̄` the change
    to local coordinates cancels exactly. `f̄` is analytic on each patch
    (its weights positive), so the rule's smoothness holds; a piece must
    lie in one patch (lines are split exactly at knot lines, spline pieces
    are located by their control points, chords by certain ends).
  * **Adaptivity.** A piece (or box) is accepted when every integrand's
    remainder is at most `2^-40` of its share (by length or area) of the
    piece's absolute integral, taken from the node values of the first
    rule, or at most four times the rounding width of its node sum;
    otherwise it is halved across the direction with the larger remainder,
    at most twelve times per direction. A box still too wide there is
    accepted as it is: always sound, only wider, and the widths are
    reported.
  * **Fallback and scope.** Where the rule cannot run (a piece across a
    patch boundary, a weight or `|N|` not certainly positive, a periodic
    spline surface) the first-order strips and Green integrals of S4d stay.
    Only mass properties use the rule (`mass_enclosure`, face areas and
    centres, sheet measures); validation's orientation fluxes and loop
    areas need signs only and keep their routes, so validation costs
    nothing more.
  * **Target widths.** Every enclosure of every spline mass fixture within
    `1e-11` relative of its quantity's scale (volume, area; the centroid
    by `V^(1/3)`, the inertia by `V^(5/3)`), asserted in tests, where S4d
    gave `2e-3`–`5e-3` on nonrational walls' areas and `5e-2`–`0.6` on
    the rational corner. The reference's quadrature (mpmath at 40 digits)
    is checked against itself on halved pieces to `1e-25` of each row's
    scale, and against closed forms where they exist, so containment is
    asserted with a slack of binary64 rounding, not `1e-12`.
  * **Evidence added.** Two spline models whose routes the fixtures do not
    yet exercise, captured natively before the kernel change:
    `spline_stadium_parallel` (the stadium with a spline pcurve along a
    cylinder's parallel: route (b) with `du ≠ 0`) and
    `spline_bulge_split_wall` (the bulge with a knot in its wall's `v`
    direction: route (c) with a column below).

  Amended during the implementation (2026-09-28):

  * **Acceptance** at `2^-44` rather than `2^-40`, with the floor taken
    from the first rule (four times its rounding width, scaled by the
    piece's relative size), so that the remainders are formed first and the
    node sum only for accepted pieces. A series undefined over a whole
    piece (an enclosure too wide to exclude a zero, as `|N|^2` over the
    whole rational wall) sends the piece to its halves instead of failing.
    The widths on the fixtures come from rounding: `2^-40` and `2^-44` gave
    the same enclosures.
  * **Nodes and weights** are bisected exactly down to `2^-100`, so that
    they lift into binary64 as adjacent values: brackets of `2^-52` had
    left node sums `1e-13` wide.
  * **Integrands are written once** over a small `Num` trait: node values
    evaluate on plain enclosures, remainders on series (the rational wall
    took twice as long with every node value a one-term series).
  * **Rational pcurves on planes** go through route (b), `F(u(τ), v(τ))
    u'(τ)` from the pcurve's own coordinates translated exactly to the
    face's reference point, not through (a): the high-degree Bernstein
    `M` of `∫ M/W^k` evaluates in binary64 with node values up to `6e-12`
    wide relative to themselves on the rounded corner's caps. (a) remains
    for rational pcurves on nonrational patches.
  * **The reference** agrees with itself on halved pieces within `1e-20` of
    each property's scale (the bulge's area differs by `1.9e-22` relative,
    the rest by `1e-30` or less), not `1e-25`; the rounded corner's weight
    is the binary64 value of `√2/2`, so only the bulges have exact closed
    forms (`20/3`, `40/3 + 8 + √2 + asinh 1`), checked to the same bound.
  * **Evidence**: `spline_rounded_corner_split_wall` beside the bulge's
    (every term of a rational wall with a column below). The stadium is
    narrow (half width `0.25`) at tolerance `1e-6`: M5's measured
    enclosure of an arc edge's use by a spline pcurve along a parallel is
    the Taylor bound on at most `2^8` pieces (`Taylor::upper_bound`),
    about `r (π/256)^3/4`, `4.6e-7` for a unit half turn, above `1e-7` and
    above the bound the reference declares. Validation itself certifies
    such a use at `1e-7`; the measurement's cap is recorded below.
  * **Work** is bounded: at most 2,048 pieces or boxes per integral, and a
    series still undefined after twelve halvings in all is a singularity,
    not a wide enclosure; either hands the integral to the first-order
    route.
  * **Target**: asserted at `1e-12` of each property's scale, since that is
    what the rule reaches with room to spare.
  * **Sheet bridge**: the bulge's lone spline wall (`sheet_spline_wall`) is
    now enclosed to `1e-15`, and OCCT's area of it is `1.8e-8` relative
    above `√2 + asinh 1` (its centre `4.6e-9` off), beyond S6's `1e-9`
    allowance; that row allows `2e-8`, as S4d allows `1e-8` for the whole
    bulge.

  F8 status (2026-09-28): implemented (decisions `b960b1a4`, evidence
  `3ad6b53f`, kernel `eca706b2`). Relative widths before → after, as
  volume / area / centroid (by `V^(1/3)`) / inertia (by `V^(5/3)`):

  | Case | S4d first order | F8 quadrature |
  | --- | --- | --- |
  | `spline_bulge` | `6e-15` / `2.2e-3` / `2e-14` / `1.3e-13` | `6.4e-15` / `3.6e-15` / `2.0e-14` / `1.3e-13` |
  | `spline_cubic_bulge` | `8e-15` / `5.4e-3` / `2e-14` / `1.6e-13` | `7.6e-15` / `3.6e-15` / `2.4e-14` / `1.6e-13` |
  | `spline_bulge_hole_inside` | `7e-15` / `2.1e-3` / `2e-14` / `1.4e-13` | `7.4e-15` / `4.9e-15` / `2.2e-14` / `1.4e-13` |
  | `spline_bulge_split_wall` | `6e-15` / `2.2e-3` / — / — | `6.4e-15` / `3.8e-15` / `2.0e-14` / `1.3e-13` |
  | `spline_rounded_corner` | `5.2e-2` / `1.4e-2` / `0.10` / `0.62` | `7.2e-15` / `5.6e-15` / `1.3e-14` / `7.4e-14` |
  | `spline_rounded_corner_far` | `5.2e-2` / `1.4e-2` / `0.10` / `0.62` | `7.2e-15` / `5.6e-15` / `5.0e-13` / `7.4e-14` |
  | `spline_rounded_corner_split_wall` | `5.2e-2` / `1.3e-2` / — / — | `7.2e-15` / `5.6e-15` / `1.3e-14` / `7.4e-14` |
  | `spline_stadium_parallel` | `2.2e-2` / `3e-15` / — / — | `4.1e-15` / `2.1e-15` / `1.5e-14` / `1.1e-13` |
  | `spline_face_c1` | `7e-15` / `7e-15` / `9e-15` / `4.5e-14` | `7.4e-15` / `3.6e-15` / `8.7e-15` / `4.5e-14` |

  The eight plane and edge cases stay at `1.5e-15` / `1.6e-15` /
  `2.6e-15` / `1.6e-14` and the two stadiums at `3e-15` (exact before and
  after). The far corner's centroid is the rounding of world coordinates
  near 7.5 against a body `1.8e-3` across. (The S4d column for the three new
  models is from the native bridge's absolute widths.) The native bridge
  (`--family spline`, 10 matches, 3 reviewed, no failures) and the sheet
  bridge pass; `spline_parallels_integrate_as_their_lines` exercises the
  cone, sphere and torus routes. Performance, the mean of ten
  `mass_enclosure` calls in release on this machine, first order → F8:
  `spline_bulge` 114 → 56 ms, `spline_cubic_bulge` 204 → 130,
  `spline_face_c1` 286 → 200 (the exact moments of nonrational walls now
  dominate), `spline_rounded_corner` 152 → 65, its knotted wall 150 → 111,
  `spline_stadium_parallel` 10.8 → 11.7, the plane and stadium cases
  unchanged (1.3–5.5). Validation is unchanged: it does not call the rule.
  A series still undefined after twelve halvings in all, or more than
  2,048 pieces or boxes in one integral, means a singularity (a degenerate
  patch edge where `|N|` vanishes) and hands the integral to the first-order
  route at once; `√τ` on the square gives up in 0.15 s in release.

  Clean local 600-second campaign of `brep_validation` at `6bb3d280`
  (AddressSanitizer, standard 20-second/2 GiB limits, allocator purge,
  fresh seed corpus of 601 inputs): 1,839 mutation executions after 312 s
  of replay (2,483 in all), 19,922 edges, peak RSS 1,307 MB, no artifact.
  Mutation 32 now also runs the rational corner, the knotted walls and the
  stadium's spline parallel and checks the area, so the target is slower
  per input than S4d's campaign (4,266 mutation executions at `ac92ec89`).

  Left open: M5's measured enclosure of an arc edge's use by a spline
  pcurve along a parallel stops at `2^8` Taylor pieces (above) and could
  refine towards the resolution; route (a) along rational pcurves on
  nonrational patches has no fixture; periodic spline surfaces stay
  unsupported (S4).

  S6's interop bridge after F8 (2026-09-29): from `dc5fd1a3` the tight
  enclosures excluded OCCT's default measures of seven free spline faces of
  `hammer.brep` (records 158, 167, 200, 212, 225, 241, 243) and
  `compare_brep_io.py` failed. None is a kernel error: an independent
  measure (`free_face_measure.py`, Green's theorem over the face's surface
  and pcurves at 40 digits, 12 and 24 nodes agreeing to 20 digits) lies
  inside every enclosure. On six, BRepGProp's default integration errs by
  `1.3e-7` to `7.2e-6` relative in area, and its adaptive integration (Eps
  `1e-12`, the second probe `occt_free_shape_adaptive_oracle.cpp`) lies
  within `MEASURE_BOUND` (243's by `2.7e-11`, OCCT's own error). On 225, a
  plane bilinear patch, the wire does not close in UV: pcurves 122 and 193
  miss by `(5e-11, -2.5e-10)` at their shared vertex (`2.5e-7` in space).
  The kernel closes its loops with the UV chord; BRepGProp sums Green's
  integrals of the open pcurves from the face's UMin, which closes the gap
  through `u = UMin` and drops a sliver `1.25` long, `3.1e-4` of area
  (`2.0e-9` relative). That open sum reproduces OCCT's default and adaptive
  values to `4e-14`; the chord-closed one is the kernel's. The seven are
  reviewed differences (`occt-brep-io-divergences.json`), fingerprinted by
  the file and both native rows, each with an independent measure the
  enclosure must contain.
* **Tessellation**, deflection-controlled and watertight, once S5 lands;
  the application needs it for display and it needs nothing from S7–S9.

  Decisions for tessellation, recorded before its code (2026-09-27):

  * **Contract.** `tessellation::tessellate(&Topology, Parameters)` (and
    `Solid::tessellate`, `Body::tessellate`) returns a `Mesh`: nodes,
    triangles grouped by face, and one polyline per edge. Every edge is
    discretized once, and every face that uses it takes exactly its
    polyline's nodes and segments as boundary, so neighbouring faces share
    nodes and the mesh is watertight by construction: no gaps, no
    T-junctions, no welding by distance. A vertex is one node at its stored
    position. Every triangle's normal (right-hand rule) leaves the solid
    region behind its face; a face between void regions follows its
    oriented normal. So a closed shell gives a closed, consistently oriented
    2-manifold mesh (every mesh edge in exactly two triangles, in opposite
    directions) whose Euler characteristic is the shell's; an open sheet a
    manifold with its boundary polylines; wire edges polylines only; an
    acorn one node. These properties are checked on every result; a failed
    check is an error, never a returned mesh.
  * **Deflection.** Two absolute bounds, as `IMeshTools_Parameters::
    Deflection` and `Angle` of OCCT's BRepMesh: a linear deflection `δ > 0`
    and an angle `0 < θ <= π/2`. Each triangle comes with the explicit map
    `Σ λ_i X_i ↦ S(Σ λ_i p_i)` onto its surface, `X_i` its nodes and `p_i`
    their parameter points on the universal cover, and each edge segment with
    the map to its curve's arc between the same fractions. The certified
    bound is that no point moves farther than `δ` under these maps (a
    parametric, Fréchet-type bound): every mesh point lies within `δ` of its
    face's surface, as BRepMesh measures deflection, and every polyline point
    within `δ` of its edge. The parameter triangles cover the face's domain
    up to the slivers between its pcurves and their chords (none along a
    line pcurve; on a plane no wider than the edge's own deflection), where
    the map lands on the surface just outside the face, so every mesh point
    lies within `2δ` of the face itself. The surface normal turns by at most
    `θ` over each triangle's parameter triangle, and an edge's tangent by at
    most `θ` over each segment. The mesh reports each triangle's and
    segment's bounds and their maxima; both maxima are at most the request.
  * **How the bounds are certified.** For a parameter triangle with extents
    `U`, `V` and bounds `a >= |S_uu|`, `b >= |S_uv|`, `c >= |S_vv|` over its
    box, linear interpolation deviates by at most `(a U^2 + 2 b U V + c
    V^2) / 8` (Taylor with integral remainder at the interpolated point, and
    the variance of a distribution on an interval of length `U` is at most
    `U^2 / 4`); the triangle's bound adds its nodes' gaps `|X_i - S(p_i)|`,
    the vertex and edge points against the face's surface at their pcurve
    points. `a`, `b`, `c` are closed forms per surface (plane 0; cylinder
    `r`, 0, 0; cone `|R + v sin α|`, `|sin α|`, 0; sphere `R |cos v|`,
    `R |sin v|`, `R`; torus `R + r cos v`, `r |sin v|`, `r`) times a bound on
    the stored frame's spectral norm, and everything is evaluated in the
    `Fast` outward-rounded tier of `certified.rs`, trigonometry included. An
    arc segment deviates by at most `r φ^2 / 8` for its sweep `φ`, a line
    segment by its end gaps. The normal turn is at most `n_u U + n_v V` with
    `|N_u|`, `|N_v|` in closed form (cylinder 1, 0; cone `cos α`, 0; sphere
    and torus `|cos v|`, 1), with a `1e-12` relative allowance for the stored
    frames' departure from orthonormality. The derivation goes in
    `MATHEMATICS.md`.
  * **Algorithm.** Edges first: each edge's segment count, uniform in its
    fraction, is the least meeting its curve's `δ` and `θ` and, for each
    curved face using it, the thin-triangle condition at `δ/2` and `θ/2`
    (a triangle on a boundary segment cannot have smaller extents than the
    segment, so the segment must leave the face room). Then each face in a
    planar chart of its domain on the cover: a plane in its frame
    coordinates; a periodic face whose loops do not wind in the sinusoidal
    chart `((u - u_c) |S_u|(v), s(v))`, `s` the arc length along `v`, which
    collapses a pole to a point; a face wound in `u` in an annulus chart
    `P(v) (cos u, sin u)`, `P` exponential in `s(v)` with its range capped,
    or the distance along the meridian from the pole when a vertex loop
    closes the band; a torus face wound in `v` the same with `u` and `v`
    exchanged. A constrained Delaunay triangulation of the chart polygons
    with exact orientation predicates (and a filtered in-circle test that
    never flips on doubt) gives the triangles; a boundary chord crossing
    another, a node on another's segment, or a loop on the wrong side of its
    chords (nesting changed by coarse chords) doubles the counts of that
    face's curved edges and starts again, within a budget. Curved faces are
    refined by Steiner points (the midpoint of the longest unconstrained
    edge, else the centroid, both strictly inside the domain) until every
    triangle's certified bounds hold, its lifted parameter triangle has the
    chart's orientation and its flat normal agrees with the surface's.
    Planar faces need no interior points. A face without loops (the whole
    sphere, the whole torus) is a structured grid in `(u, v)` refined until
    every triangle's bounds hold. Everything is deterministic: fixed
    orders, no hashing, no threads, a pseudo-random walk with a fixed seed.
  * **Seams, poles and windings.** No seam is meshed: the charts are
    homeomorphisms of the face's domain, so there is nothing to weld (OCCT
    duplicates seam nodes and joins them through its seam edges). A pole is
    one node (the vertex of its vertex loop, or of a loop passing through
    it); a triangle's parameter point at a pole takes the mean `u` of its
    other two, since every `u` maps there. A loop wound more than once, or
    in both directions, and a vertex loop that is not a pole are
    `OutOfDomain`.
  * **Degenerate and limit cases.** Non-finite or non-positive deflection,
    or an angle outside `(0, π/2]`: `NonFinite` or `OutOfDomain`. A node gap
    above `δ/2` (a pcurve far from its edge, possible only on imported or
    invalid bodies): `InvalidTopology`. Boundary chords that still cross
    after the refinement budget, or crossing line edges: `InvalidTopology`.
    More than 4,000,000 nodes: `ComputationLimit`. Tessellation does not
    validate the body, but every bound it reports is certified.
  * **Evidence first.** An independent reference
    (`tessellation_reference.py`, mpmath) derives each fixture body's
    boundary from the case alone (planar regions with arcs, cylinder
    patches, the meridians of cones, spheres and tori, exact distances from a
    point to each), its exact area, volume and Euler characteristic, and
    checks a mesh: closedness and opposite orientations, the Euler
    characteristic, nodes on the boundary, sampled deflection from each
    face's whole surface against the request and against each triangle's
    reported bound, sampled distance from the boundary within `2δ`, edge
    polylines within `δ`, normals leaving the solid, and volume within `δ`
    times the areas. A native capture of `BRepMesh_IncrementalMesh` on
    the same bodies, read from the `.brep` text the kernel's writer already
    produces, records node and triangle counts, OCCT's own deflection, the
    deflection measured by sampling, watertightness after joining nodes
    through OCCT's edge polygons, orientation, area and volume, before any
    kernel tessellation code exists. Differences (OCCT beyond the requested
    deflection, gaps) are reviewed with fingerprints.
  * **Split.** T-a (now): faces on planes, cylinders, cones, spheres and
    tori with line, circle and arc edges: profiles, prisms, the S3
    primitives, S6 sheets and wires, and imported bodies of those. T-b:
    spline edges, pcurves and faces, with second-derivative bounds from the
    second differences of the control net over each Bézier piece (rational
    ones through the weights' lower bound); until then spline geometry is
    `OutOfDomain`. T-c: procedural intersection edges (D13) once faces
    carry them, through their certified parameterisation. Later, not
    decided: relative deflection, a minimum-angle quality guarantee,
    parallel faces.

  T-b decisions, recorded before its code (2026-09-28):

  * **Scope.** Spline edges (`Curve3::BSpline`, rational or not, periodic
    or not, over any range, reversed spans and ring edges over a full
    period included), spline pcurves on every surface T-a meshes and on
    spline surfaces, and faces on nonperiodic spline surfaces, rational or
    not, with the other T-a kinds beside them. Out of scope, `OutOfDomain`:
    a spline surface periodic in either parameter (validation, mass
    properties and the reader do not certify or take one either; S4); a
    spline face without loops or with a vertex loop (no poles); a spline
    edge or surface that is not C1 inside its range, by R4's exact test
    (the bounds below need C1 across knots; pcurves enter no bound, only
    their nodes' gaps, so they are not tested). Where a spline surface's
    normal vanishes in a face (a collapsed row of poles) or an edge's
    derivative vanishes, no turn can be certified: `ComputationLimit`,
    after a bounded refinement.
  * **Contract.** T-a's, unchanged. A spline edge's node at fraction `k/n`
    carries the parameter `u_k = SplineSpan::parameter(k/n)` (binary64,
    monotone, the range's ends exact), and a segment's certified map is
    onto the curve's arc between `u_k` and `u_{k+1}`, affine in the
    parameter; a spline face's parameter points are in the surface's own
    `(u, v)`.
  * **Bounds.** From the exact Bézier pieces of the span's range
    (`bezier_arcs_in`) and patches of the surface (`bezier_patches`),
    enclosed in the `Fast` tier and halved by de Casteljau (curves three
    times, surfaces twice per direction, less when a surface has many
    patches) into cells. On a cell of parameter length `L` (or `L_u`,
    `L_v`), degree `p` (`q`), homogeneous controls `(w_i P_i, w_i)`,
    `Q_i = P_i - c` for its first control `c`: `R = max |Q_i|`,
    `ω = min w_i` (the weights' lower bound, positive), the numerator
    `A = w (C - c)` and `w` bounded by their control nets' first and
    second differences (`|A'| <= p/L max |Δ(wQ)_i|`, `|A''| <= p(p-1)/L²
    max |Δ²(wQ)_i|`, mixed differences for `A_uv`), then
    `|C'| <= D1 = (|A'| + |w'| R)/ω` and `|C''| <= D2 = (|A''| + 2 |w'| D1
    + |w''| R)/ω`, and for a surface `|S_uu|`, `|S_uv|`, `|S_vv|` the same
    way from `A = w (S - c)` (`MATHEMATICS.md`). A nonrational cell reduces
    to `p(p-1)/L² max |Δ² P_i|`. A segment of parameter length `h` deviates
    by at most `D2 h² / 8` over the cells it meets, plus its end nodes'
    gaps, and its tangent turns by at most `D2 h / s`, `s` the enclosed
    `|C'|` at its middle less `D2 h / 2` (no bound when not positive). A
    triangle deviates by T-a's `(a U² + 2 b U V + c V²) / 8` with `a`, `b`,
    `c` over the cells meeting its parameter box, plus its nodes' gaps; its
    normal turns by at most `(M_u U + M_v V) / m`, `M_u = a D_v + D_u b >=
    |∂_u (S_u × S_v)|`, `M_v = b D_v + D_u c`, `m` the enclosed
    `|S_u × S_v|` at the box's centre less `(M_u U + M_v V) / 2`. T-a's
    tangential correction and fan map rely on analytic structure and are
    not used on splines. Node positions and gaps are evaluated by de
    Casteljau on the cells in the same tier.
  * **Algorithm.** A spline edge starts from the least uniform count its
    largest `D2` allows (0.9 of the request, as T-a), then every segment's
    certified bound is checked and any excess doubles the count. Every edge
    with a fin on a spline face (the line edges of a spline wall too) also
    needs, for each boundary segment, T-a's thin-triangle condition over
    the cells of its chord's box (0.45 of the request); a failure doubles.
    At most twelve doublings, then `ComputationLimit`. A spline pcurve on
    an analytic surface gives T-a's thin-triangle condition its extents
    from `D1` in `u` and `v` and its `v` range from its control hull. A
    spline face is meshed in the affine chart `(g_u (u - u_c), g_v (v -
    v_c))`, `g_u`, `g_v` the lengths of `S_u`, `S_v` at its region's
    centre, and refined as T-a's curved faces; a triangle still failing
    when its parameter box is within `2^-30` of the face's parameter span
    in both directions is `ComputationLimit`. Analytic-only topologies mesh
    exactly as in T-a.
  * **Evidence first.** `tessellation_reference.py` gains planar B-spline
    pieces (exact Bézier pieces in `Fraction`s, distances by projection in
    binary64 from sampled starts, crossing parity from each piece's roots,
    exact area and length by mpmath), so prisms with spline profiles are
    `Prism`s, plus a dome (a box under a spline graph) and spline sheets
    with holes, measured by projection onto their surfaces and edges.
    `generate_tessellation_fixtures.py` adds `tessellation-spline-cases.txt`
    and `-expected.tsv`: twelve bodies in the B-rep line protocol
    (`generate_brep_fixtures.prism`: quadratic, cubic and rational spline
    walls, a wave with a spline and a circular hole, a sharp rational
    corner, a rigidly moved cubic bulge, degree-1 spline edges and pcurves
    on a stadium's cylinder; a plane face in a periodic spline ring; the
    dome; a trimmed bicubic sheet whose hole edges are iso-curve ranges; a
    rational biquadratic sheet) at T-a's two settings, the T-a files
    unchanged; `--check` in CI. The native capture
    (`occt-spline-tessellation-preimplementation`, the same probe on the
    kernel writer's `.brep` texts, `brep_io_probe parts`) comes before any
    kernel code, its `capture.json` recording
    `rust_spline_tessellation_exists: false` (no
    `tessellation/spline.rs`). The kernel's meshes must pass every check of
    the reference; OCCT's differences get fingerprinted reviews. The kernel
    tests also mesh every valid spline case of `brep-cases.txt` and the
    certified `data/occ` solids with spline geometry; the `tessellation`
    fuzz target gains spline prisms and sheets.
* **STEP import**, after S6, reusing the converter architecture, with the
  OCCT STEP reader as the native oracle.

  STEP import decisions, recorded before its code (2026-09-28, by the
  implementing agent of the parallel track; none is a user decision):

  * **Architecture.** `step::read` parses a Part 21 exchange structure into
    instances by entity number; `step::import` translates its bodies into
    OCCT's shape structure, the `occt_brep::Document` the `.brep` reader
    produces (as `StepToTopoDS` translates into `TopoDS`), and converts that
    with the existing `occt_brep::import`: seam merging, poles, windings,
    orientations, enclosures and validation are the `.brep` converter's,
    unchanged. What the STEP layer adds is the entity graph, units, the
    uncertainty, the degenerated edges OCCT's reader adds at poles, and the
    pcurves STEP does not carry. Every body passes `Topology::from_parts`
    or is an import failure with the validator's issues, as a `.brep` body.
  * **Schemas.** AP203 (`CONFIG_CONTROL_DESIGN` and the second edition's
    `AP203_CONFIGURATION_CONTROLLED_3D_DESIGN_OF_MECHANICAL_PARTS_AND_ASSEMBLIES_MIM_LF`),
    AP214 (`AUTOMOTIVE_DESIGN`) and AP242
    (`AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF`), geometry and
    topology only; any other `FILE_SCHEMA` is `StepError::Schema`. The
    parser itself is schema-free: Part 21 edition 2 syntax (header, one or
    more `DATA` sections, simple and complex instances, every parameter
    kind, comments, `''` and `\X\`, `\X2\`, `\X4\` string encodings kept
    raw), with a nesting limit and entity numbers checked. *Amended in
    implementation:* AP214's committee drafts (`AUTOMOTIVE_DESIGN_CC1`,
    `_CC2` and the like, named by 61 of the local dataset's files) have
    its geometry and topology and are read; inside a string `\\` and `\S\`
    with the character after it are kept raw (`\S\'` does not end the
    string).
  * **Entities.** Bodies: `MANIFOLD_SOLID_BREP` (a solid), `BREP_WITH_VOIDS`
    (a solid whose `ORIENTED_CLOSED_SHELL`s with orientation `.F.` are
    cavities, as OCCT writes and reads them) and `SHELL_BASED_SURFACE_MODEL`
    (each `OPEN_SHELL` or `CLOSED_SHELL` a free shell: a sheet or a closed
    shell, S6). Topology: `CLOSED_SHELL`, `OPEN_SHELL`,
    `ORIENTED_CLOSED_SHELL`, `ADVANCED_FACE` (and `FACE_SURFACE`),
    `FACE_BOUND`, `FACE_OUTER_BOUND`, `EDGE_LOOP`, `ORIENTED_EDGE`,
    `EDGE_CURVE`, `VERTEX_POINT`, and `VERTEX_LOOP` as the only bound of a
    sphere or torus face (the whole surface, as OCCT reads it). Geometry:
    `CARTESIAN_POINT`, `DIRECTION`, `VECTOR`, `AXIS2_PLACEMENT_3D`, `LINE`,
    `CIRCLE`, `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`,
    `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`; `TRIMMED_CURVE`,
    `SURFACE_CURVE` and `SEAM_CURVE` as their basis or 3D curve (an edge is
    bounded by its vertices; OCCT also takes the 3D curve); later
    sub-steps add `ELLIPSE`, the B-spline curves and surfaces (rational as
    complex instances) and `PCURVE`. Orientations follow
    `StepToTopoDS`: an edge runs along its curve, its vertices swapped when
    `same_sense` is false; a wire is reversed when its bound's orientation
    differs from its face's `same_sense`, and the face is reversed when
    `same_sense` is false. Anything else a body reaches is reported by
    entity name and counted, never approximated.
  * **Bodies found.** Every solid and surface-model item, in entity-number
    order, imported once, in the coordinates of the representation listing
    it. Placements between representations (assemblies, `MAPPED_ITEM`) are
    not applied: such a file reports `AssemblyPlacement` and its bodies are
    still imported in their own coordinates.
  * **Units and uncertainty.** The length unit of the body's
    representation context (`SI_UNIT` with its prefix, or a
    `CONVERSION_BASED_UNIT` through its `LENGTH_MEASURE_WITH_UNIT`) scales
    every coordinate and radius to millimetres by one binary64
    multiplication (none in millimetres), as OCCT's default
    `xstep.cascade.unit`; the plane angle unit (radian, or a conversion
    such as degrees) scales cone semi-angles. The context's
    `UNCERTAINTY_MEASURE_WITH_UNIT` for length, scaled alike, is each
    vertex's, edge's and face's imported tolerance, floored at `1e-7` mm and
    capped at `1` mm (OCCT's `read.precision.mode` "File" and
    `read.maxprecision.val`); without one, `1e-7` mm. A body's resolution
    is then that tolerance, as for a `.brep` body. *Amended in
    implementation:* OCCT uses the uncertainty only as its healing
    precision and builds every entity at `Precision::Confusion` (its
    observed tolerances stay `1e-7` on files declaring `0.1`); the kernel,
    which does not heal, keeps the file's claim, verified by the validator.
    A measure with unit may be a complex instance whose value and unit sit
    in its `MEASURE_WITH_UNIT` record.
  * **Pcurves.** STEP-a derives every pcurve and ignores the file's: on a
    plane, the edge's own data in the plane's coordinates (the `.brep`
    converter's `CurveOnPlane`); on a cylinder, cone, sphere or torus, the
    edge's image under the surface's inverse map when that image is a
    straight segment in `(u, v)`: rulings of cylinders and cones, parallels
    (circles about the axis), and meridians of spheres and tori (circles in
    a plane through the axis). Each use's segment continues the last in the
    universal cover; a seam's second use lies one period from its first;
    at a pole (a sphere's, a cone's apex) the next use's `u` is the nearest
    branch in the direction that keeps the face on the loop's left, and the
    gap becomes a degenerated edge, as OCCT's reader adds one. The
    validator certifies every derived pcurve against its edge within the
    imported enclosure, so a wrong derivation is an import failure, never a
    wrong body. Any other edge on a curved surface is `PCurveNotDerived`,
    and a periodic face whose loops wind without a seam is
    `PeriodicFaceWithoutSeam`, for now. Later sub-steps take a file's
    `PCURVE` where nothing exact can be derived (B-spline surfaces).
  * **Errors and determinism.** Malformed text is `StepError::Syntax`
    (line and expectation), a reference to a missing instance
    `StepError::Reference`, a wrong schema `StepError::Schema`; a body's
    unsupported constructs are its `Rejected::Unsupported` names and are
    also counted per file, as `.brep` import does. Nothing panics on any
    input (the `step` fuzz target). The result is a function of the text:
    bodies by entity number, every table in a fixed order, no hashing.
  * **Entity ids.** Each imported body records its STEP entity number;
    the topology's ids are the `External` derivations
    (`OperationKind::External`, role `External`) that `.brep` import gives,
    one per slot. STEP names, colours, layers and product structure are not
    carried.
  * **Out of scope.** Assemblies and placements, colours, names, layers,
    PMI and validation properties, tessellated and polyloop geometry,
    offset and swept surfaces, Part 21 edition 3 anchors and references,
    and writing STEP.
  * **Sub-steps.** STEP-a: the parser, and solids and sheets on planes,
    cylinders, cones, spheres and tori bounded by lines and circles, with
    units, uncertainty and derived pcurves. STEP-b: ellipses (arcs on
    planes, sinusoid pcurves on cylinders) and B-spline curves and surfaces
    with the file's `PCURVE`s on spline surfaces. STEP-c: a local survey of
    the dataset's 322 `.stp` files (U1: local only, nothing committed) by
    unsupported construct, which orders what follows. Each keeps the
    evidence order: an independent reference (`step_reference.py`) of
    fixtures authored here, then a native `STEPControl_Reader` capture
    (counts, validity, volume, area, centre) before the kernel's importer
    for it exists, then the importer, a probe, `compare_step.py` with
    fingerprinted reviews, the `step` fuzz target and the docs.

  STEP-b decisions, recorded before its code (2026-09-28, by the
  implementing agent of the parallel track; none is a user decision). The
  decisions above hold; these settle what they leave open for ellipses,
  B-splines and the file's pcurves:

  * **Ellipses.** `ELLIPSE(position, semi_axis_1, semi_axis_2)` is
    `c + a1 cos t x + a2 sin t y` in its placement, in either order of the
    semi-axes (OCCT's `StepToGeom` turns the axes a quarter turn when the
    second is longer; the kernel's `Curve3::EllipseArc` takes either, so
    the file's parameter is kept). An edge's range is the vertices' angles
    `atan2(y / a2, x / a1)`, a closed ellipse one turn, as for a circle.
    Pcurves are derived, as in STEP-a: on a plane, the ellipse in the
    plane's coordinates when its axes are the plane's (`Curve2::EllipseArc`,
    S8a.2), else its exact projection onto the plane (`Curve2::Projection`,
    D13); on a cylinder, when the ellipse is the section of the cylinder by
    its own plane (its centre on the axis, its axes projecting onto the
    axis's normal plane as two perpendicular radii), the sinusoid
    `v = a0 + a1 cos u + a2 sin u` of S8a.2 over `u = +-t + phi`, placed on
    the universal cover like the straight segments. An ellipse on a cone,
    sphere or torus, or one that is not such a section, is
    `PCurveNotDerived`.
  * **B-spline curves and surfaces.** `B_SPLINE_CURVE_WITH_KNOTS` and
    `B_SPLINE_SURFACE_WITH_KNOTS` as simple instances, and as complex
    instances whose `B_SPLINE_CURVE` (`B_SPLINE_SURFACE`) record holds the
    degree(s), control points and flags, `B_SPLINE_CURVE_WITH_KNOTS`
    (`..._SURFACE_WITH_KNOTS`) the multiplicities and knots, and
    `RATIONAL_B_SPLINE_CURVE` (`..._SURFACE`) the weights; the other records
    of a complex instance (`BOUNDED_CURVE`, `GEOMETRIC_REPRESENTATION_ITEM`,
    `REPRESENTATION_ITEM` and the like) carry nothing geometric. Control
    points scale with the length unit; knots and weights do not. Surface
    control points are rows in `u` (the kernel's `u`-major order). Curve
    form, closure flags, self-intersection and knot type are not used. A
    B-spline is non-periodic when its multiplicities sum to poles plus
    degree plus one (ISO 10303-42's clamped or unclamped form), periodic
    when they follow OCCT's periodic pattern (`StepToGeom`'s test), and
    otherwise `InvalidBSplineCurve` / `InvalidBSplineSurface`; OCCT also
    makes a closed curve flagged closed periodic (`SetPeriodic`), which
    moves no point. The kernel's constructors decide validity and limits
    (`BSplineControlDataLimit`), as for `.brep` records. `BEZIER_*`,
    `UNIFORM_*` and `QUASI_UNIFORM_*` curves and surfaces are refused by
    name for now (OCCT converts them to knotted B-splines); the survey
    counts them.
  * **Edge ranges on B-splines.** An edge runs along its curve from the
    parameter where its start vertex lies to where its end vertex lies:
    the domain's end when the vertex is within the entity tolerance of it,
    else the nearest point of the curve found by sampling each span and a
    bracketed refinement (OCCT's `ShapeAnalysis_Curve::Project`); a closed
    curve used by one vertex runs its whole domain. A start at or after the
    end on a non-periodic curve is `EdgeRangeInverted`. The validator then
    certifies the vertices on the edge, so a wrong location is an import
    failure with `vertex_off_curve`, never a moved edge. (Amended at the
    fix of CI's scheduled `step` timeout, 2026-10-10: the ends are tested
    on correctly rounded points, the samples and the refinement evaluate
    in binary64 de Boor, the location being only the claim the validator
    certifies.)
  * **The file's pcurves.** On a B-spline surface the edge's `SURFACE_CURVE`
    (or `SEAM_CURVE`) must carry a `PCURVE` whose `basis_surface` is the
    face's surface; its `DEFINITIONAL_REPRESENTATION` holds a 2D `LINE`
    (its point and direction; the vector's magnitude only scales the
    parameter) or a 2D B-spline curve (simple or rational complex, as
    above). Its range is found as the edge's, locating the edge's vertices
    on the pcurve's image `S(c(t))`; the pcurve and the edge must then
    share their parameter up to an affine map, which the validator
    certifies at matching fractions (S4b): a file whose pcurve is
    parameterised otherwise is an import failure with its issues, where
    OCCT reparameterises (`BRepLib::SameParameter`). A pcurve running
    against its edge is `PCurveAgainstEdge`. Without a `PCURVE` on the
    face's surface, or with another 2D curve, the use is
    `PCurveNotDerived` (OCCT projects the edge; the kernel does not
    approximate). The first matching `PCURVE` is used; an edge used twice
    by one spline face is the converter's `SeamOnBSplineSurface`. On
    planes and on the analytic curved surfaces the file's pcurves are still
    ignored and derived as in STEP-a (OCCT, too, recomputes them on
    planes): an edge on a cylinder, cone, sphere or torus that STEP-a and
    the ellipse rule above cannot derive stays `PCurveNotDerived`, even
    when the file carries a pcurve. On a plane a B-spline edge's pcurve is
    its poles in the plane's coordinates (the converter's `CurveOnPlane`).
  * **OCCT's shape structure.** The `occt_brep::Document` gains
    `Curve3::Ellipse` (the file's semi-axes in their order) and
    `Curve2::Sinusoid`, which only the STEP translation produces: the
    `.brep` reader keeps naming ellipse records (`Ellipse`, `Ellipse2d`)
    unsupported, since reading them is a `.brep` change with its own
    evidence (it would move the recorded `.brep` dataset survey), and OCCT
    itself would hold an approximated B-spline pcurve where the kernel
    holds the exact sinusoid. The converter maps them to
    `Curve3::EllipseArc` (one turn snapped exactly, a closed ellipse losing
    its seam vertex becoming a ring like a circle) and `Curve2::Sinusoid`,
    and its plane pcurve of an ellipse whose axes are not the plane's
    becomes the exact projection. Because the reader does not read ellipse
    records, the `step` fuzz target's `.brep` round trip is skipped for a
    body whose written text carries one.
  * **Evidence.** New fixtures of `generate_step_fixtures.py`, each with
    closed forms or mpmath quadrature at 30 or more digits: a half-ellipse
    sheet on the plane of its axes; a cylinder cut by an oblique plane (an
    ellipse on a plane not sharing its axes, and a sinusoid on the
    cylinder); a sheet bounded by a line and a two-span cubic B-spline; a
    prism over that profile whose side is a B-spline surface with the
    file's line pcurves; a two-span bicubic patch sheet; a single-span
    patch trimmed by a curve whose pcurve is a quadratic B-spline (its 3D
    curve the exact composition, degree 9); and a cylinder of two rational
    half-circle surfaces and rational half-circle edges as complex
    instances. `step_reference.py` gains its own B-spline evaluation and
    checks that every edge's curve meets its vertices and lies on each
    face's surface (through the file's pcurve on a spline surface, at
    matching parameters) within `1e-12` of the case's size. The native
    capture of these files (`fixtures/occt-step-b-preimplementation`) is
    taken while `step/spline.rs`, the translation of the new entities,
    does not exist; STEP-a's capture is unchanged.
  * *Amended in implementation:* `PCurveAgainstEdge` applies to B-spline
    pcurves only (a line pcurve's range may run either way), and a line
    pcurve that misses the spline surface's domain is
    `PCurveOutsideSurface`. The rational cylinder's circles, rational
    quarter arcs with knots of multiplicity equal to the degree (the usual
    NURBS circle), are C1 as rational curves but not in homogeneous form,
    which R4 certifies: the kernel translates every entity and its
    validator refuses the body (`edge_not_c1`, `pcurve_not_c1`,
    `face_not_c1`). The fixture and its capture stay as taken;
    `compare_step.py` reviews a validation refusal fingerprinted by its
    issue kinds (`rust_invalid:`), and rational complex instances are shown
    to import by the plate and prism rewritten with equal weights and a
    plate bounded by a rational third of a circle (`tests/step.rs`). A sheet
    of one face keeps the kernel's synthesized counts of a free face (S6:
    no shell) where OCCT's reader keeps the shell; each such fixture has a
    review.

### Decisions pending from the user

* **U6.** Approve the CI budget policy above? **Answered 2026-09-27: yes.**
  Per-push fuzz runs replay every checked-in regression, every seed added
  since the last schedule run and a bounded random sample of the rest, then
  mutate; the daily schedule replays everything; `surface_knots`,
  `degree_elevation` and `surface_editing` run on schedule only, their
  regressions still replayed per push; Windows runs a smoke subset per push
  and the full suite nightly. Acceptance still requires the full replay green
  on the schedule run of the accepted revision, plus the clean local
  600-second campaign. Record the sampling rule and its seed in `FUZZING.md`
  so a per-push run is reproducible.
* **U7.** Intersection curves that are not conics or splines: procedural
  cells that evaluate through their two surfaces with certified enclosures
  (as CGM's edge curves and Parasolid's SP-curves do), approximated by
  splines only for interchange with the approximation's bound recorded; or
  spline approximations as the kernel's own representation, as OCCT does.
  **Answered 2026-09-27: procedural.** Record as D13 in `TOPOLOGY_MODEL.md`:
  an intersection edge stores its two surfaces and a certified
  parameterisation, evaluates by a certified iteration with an enclosure,
  lies on both faces by definition so its pcurves are exact projections,
  and is approximated by a spline only for tessellation and interchange,
  with the approximation's bound recorded beside it. Conics and splines
  stay explicit where the intersection is one.
* **U8.** S5 before S6, application value first; or S6 before S5, upstream
  evidence and Boolean tools first. **Answered 2026-09-27: S5 then S6.**
* **U9 (recorded 2026-09-27 by the implementing agent, T-a of the
  tessellation track).** Tessellation certifies every triangle's deflection
  over the whole triangle, where BRepMesh samples a triangle's centre and
  link middles. On the fixtures that costs 1.0–1.2 times BRepMesh's
  triangles on prisms, 0.2–0.3 on apex cones, 0.9–2.6 on spheres, frusta
  and tori, while BRepMesh exceeds the request on seven of 56 rows (up to 4.4
  times). Keep the certified contract as the only mode, or add a display
  mode with sampled control and no bound (fewer triangles, OCCT's
  behaviour)? Until answered, only the certified mode exists.
* **U10 (recorded 2026-09-28 by the implementing agent, from STEP-b).**
  R4's continuity criterion is homogeneous: a knot of multiplicity equal to
  the degree passes only if the homogeneous curve is C1 there. The usual
  NURBS circle (rational quarter arcs joined at double knots, every STEP
  and IGES writer's circle) is C1 as a rational curve but not in
  homogeneous form, so the kernel refuses it (`edge_not_c1`,
  `pcurve_not_c1`, `face_not_c1`): the STEP-b fixture `rational_cylinder`
  and 21 bodies of the local dataset. Options: (a) keep the homogeneous
  criterion (such files stay refused); (b) test C1 of the rational curve
  exactly at such knots (`(P' w - P w') / w^2` equal from both sides, in
  rationals), keeping OCCT's cells and counts; (c) split edges and faces at
  such knots on import (each piece C-infinity; more cells than OCCT's).
  The implementing agent recommends (b). Until answered, (a) holds.

## Status

* S1 — **done**, accepted at `7e463cb2` (record in
  `IDENTITY_AND_HISTORY.md`). R1: the CI input takes
  8.0–8.2 s locally under AddressSanitizer (4.1 s without), above
  20 s ÷ 2.6, so `surface_editing` has the 60-second budget; both inputs are
  regressions with their times in `fuzz/regressions/README.md`, and
  `FUZZING.md` states the triage rule. R5: `History::steps`, the
  `steps_invalid` check and the two-level test in `history_contracts.rs`.
  R8: the table in `VALIDATION.md`.
* S2 — **done**, accepted at `1a76d29e` (record in `UPSTREAM_TESTS.md`).
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
* S3 — **done**: cone accepted at `166fc905`, sphere at `059616f9`, torus at
  `ac92ec89` (records in `TOPOLOGY_MODEL.md`).
  * Native `MakeCone` capture before any kernel cone code (`dde086c1`);
    the validator's cone and pole rules against the independent reference
    (`1a76d29e`, 88 reports).
  * Builder, ids and history (`Revolve`, role `apex`, the meridian as
    boundary 0), general certified mass properties (U2: exact integration
    in the two tiers rather than quadrature, so the error bound is the
    enclosure itself), `.brep` import and export of cones with OCCT's
    degenerated apex edge, `pcone` and body properties in the DRAW
    adapter, the derived case `pcone_counts`, and cones in four fuzz
    targets.
  * The `MakeRevol` history capture followed the builder: R11, accepted
    by the user. 21 matches, one reviewed difference
    (`nearly_cylinder`: MakeRevol's surface of revolution).
  * Upstream: one more `data/occ` solid imports (30 of 77); six restore-only
    cases no longer report cones, none evaluates yet (`UPSTREAM_TESTS.md`).
  * Clean local 600-second campaigns at `166fc905` (the cone with the fast
    mass engine), AddressSanitizer, standard limits: `brep_validation`
    9,644 executions after 140 s of replay, `identity` 2,013 after 24 s,
    `history` 7,493 after 8 s, `brep_io` 10,091 after 6 s; no artifact.
    At `0c94aa53` CI failed only the revolve bridge's Linux fingerprint (now
    reviewed) and the `identity` and `surface_knots` replays (R12).
  * Sphere implemented after its native captures (`3641a2f7`, before any
    kernel sphere code): the whole sphere as the first face without loops,
    poles on the band's side, loops through a pole, the builder, ids,
    history (22 of 22 matches with `MakeRevol`), mass, interop (fourteen
    more `data/occ` solids import, 44 of 77), `psphere`, `psphere_counts`
    and fuzz. `bugs/modalg_6/bug27264_2` passes on both backends: the
    ledger has 2 mapped-and-verified assertions. Accepted at `059616f9`.
  * Torus implemented after its native captures (`059616f9`, before any
    kernel torus code; OCCT's `inner_half` is inside out, reviewed):
    windings in `v`, the whole torus without loops, v-segments and wedges,
    the builder, ids, history (22 of 22 matches with `MakeRevol`), mass
    (the u↔v exchanged routine for loops wound in `v`), interop (ten more
    `data/occ` solids import, 54 of 77; horn and spindle tori are
    `NonRingToroidalSurface`), `ptorus`, `ptorus_counts` and fuzz
    (`brep_validation` mutation 30, tori in `identity`, `history` and
    `brep_io`). `bug485` now restores and needs only `bfuse`. Accepted at
    `ac92ec89`.
* S4 — **done**: R4 and S4 accepted at `ac92ec89` (records in
  `TOPOLOGY_MODEL.md`).
  * Spline variants in the topology (`Curve3::BSpline`, `Curve2::BSpline`
    with the new `BSplineCurve2`, `Surface::BSpline`); `edge_not_c1`,
    `pcurve_not_c1` and `face_not_c1` by exact removal
    (`topology/validate/continuity.rs`), independently checked by one-sided
    derivatives (`spline_cell_reference.py`) on nineteen fixture cases;
    `brep_validation` mutation 31. Every other check of spline geometry is
    uncertified, as decided above; the reader and writer still refuse
    splines.
  * S4a: native observations of ten spline models (spline sides with
    ruled spline walls, spline geometry on a stadium's cylinder, and
    mutations) before any kernel code certifies spline geometry
    (`fixtures/occt-spline-preimplementation`).
  * S4b-d, first part: exact composition for rational uses
    (`validate/spline_deviation.rs`), enclosed spline areas and plane
    fluxes, and the flux of nonrational spline surfaces
    (`validate/spline_flux.rs`), over shared Bernstein arithmetic in both
    tiers (`validate/bernstein.rs`); the reference samples exactly and
    integrates by quadrature. Six spline cases are now valid, two of them
    spline models that native OCCT also finds valid with the same counts and
    enclosures within its measurements; `brep_validation` mutation 32 moves
    them by exact similarities.
  * S4b-d, second part: Taylor enclosures for every other spline use
    (`validate/spline_taylor.rs`), containment by exact parity, periodic
    areas, a cylinder's flux with constant-u spline pcurves, mass
    properties with spline pcurves on planes, and `SplineSpan` ranges (the
    corpus trims 451 spline edges). 60 fixture cases are valid. Still to
    do: rational spline surfaces' flux, mass properties on spline surfaces
    and with spline pcurves on curved surfaces, periodic spline surfaces,
    and S4e.
  * S4b-e, third part: the flux of every nonperiodic spline surface
    (strips of the `v` domain from the surface jets when the exact route
    does not apply), fluxes along spline pcurves on cylinders, cones,
    spheres and tori, UV gaps on spline surfaces in 3D, pcurves against
    their spline's parameter, and interop: records 7 and 9 in the reader
    and writer, trimmed bases, ranges snapped at printing precision,
    periodic spline surfaces and seams on spline surfaces unsupported as
    decided above. 62 fixture cases are valid; the spline bridge has 7
    matches and 3 reviewed differences (`spline_c0_bulge`: BRepCheck has no
    continuity status; the two shifted pcurves: OCCT reports the face
    unorientable after the edge's InvalidCurveOnSurface, and Closed2d
    accepts the gap). 60 of 77 `data/occ` solids are representable and 58
    certified; `Ball.brep` 108 (`seam_edge`) and `Motor-c.brep` 378
    (`uncertified_containment`) are pinned. The restore-only survey after
    S4 reports no B-spline construct; it found the importer panicking on
    `bug21246`'s internal edges and faces (fixed, with a test and a fuzz
    mutation).
  * S4d, mass: `BRepGProp` captured first (`fc389e8d`), then mass
    properties with every spline geometry but a periodic spline surface, as
    decided above; the reference's quadrature (`brep-spline-mass.tsv`, 16
    valid spline cases) lies in every kernel enclosure, and so do OCCT's
    properties of the five spline models valid on both sides. Enclosure
    widths on the models: volume `3e-14`–`5e-14` on nonrational walls and
    `0.31` (5%) on the rational corner, area `0.05`–`0.13` on nonrational
    walls. Mutation 32 checks moved models' volumes and centroids. Its
    first campaign found the rational corner, shrunk by `2^-10` and moved
    8 away, with an undecided shell orientation after minutes in the exact
    tier: rational jets enclosed as `X(box)/w(box)` in absolute coordinates
    lose the body's scale. Fluxes are now taken relative to a vertex of the
    shell, the strips lift patches translated to it exactly (`X - r w`),
    and loop areas, plane fluxes and plane mass integrals run in UV
    relative to a point of the loop or face; the case is the fixture
    `spline_rounded_corner_far`, and the rational corner stays out of
    mutation 32 (seconds per input under AddressSanitizer).

### Answered by the user (2026-09-27)

* **U1:** the dataset stays local-only; CI does not fetch it.
* **R11:** the cone's `MakeRevol` history capture is accepted as taken
  (after the builder was written, before it was committed).
* **R12:** retained corpora are minimised periodically: `run_fuzz.py
  --minimize` merges each seeded corpus into a fresh directory with
  libFuzzer's `-merge=1` under the target's own input and memory limits and
  replaces it only when the merge completes; the fuzzing workflow runs it
  instead of a campaign every Sunday (and on dispatch with `minimize`), and
  the cache keeps the result. Locally the `brep_io` corpus went from 1,872
  inputs to 1,108 in 190 s. The per-target startup caps stay as a backstop.
* **S4:** periodic spline surfaces and seams on spline surfaces stay
  unsupported, and spline mass enclosures stay first order where not exact;
  a higher-order certified quadrature is later work.
* U6 (parallel track) — implemented: per-push sampled replay with a
  manifest of the last full replay, schedule-only exact tensor targets,
  Windows smoke per push and a nightly kernel schedule (`FUZZING.md`).
* S5 — implemented at `4592ef6c`; kernel CI green there, and the clean
  local 600-second campaigns of `identity`, `modeling`, `history`,
  `brep_io` and `brep_validation` at that revision ran without an artifact
  (1,360 to 7,681 mutation executions each, full replays). Acceptance waits
  for the full replay of the next schedule run, which will also cover S6.
  * Native `MakePrism` observations of thirteen arc prisms captured before
    any kernel arc code (`88ae24df`), with arc paths in the identity
    reference and the native probe.
  * `Boundary::path` with `Segment::Line`/`Segment::Arc`, certified arc
    predicates (`decide/arcs.rs`), closed-form moments, exact ray parity;
    partial cylinder walls in `Topology::prism`; arc pcurves on caps. The
    history bridge's arc family matches 13 of 13, the interop bridge reads
    577 written prisms (13 arc prisms) natively valid with equal properties,
    `identity.rs` matches the reference's ids, `arc_profiles.rs` covers the
    path contract, the builder reproduces the neutral generator's arc
    prisms, the DRAW adapter runs `profile` (derived case `profile_arcs`,
    both backends) and the `identity` fuzz target builds filleted, notched
    paths. Polygon and circle ids are byte-identical (the full suite).
* S6 — implemented; gate pending CI, the schedule replay and the clean
  campaigns.
  * Native observations of eighteen sheet, shell, wire and acorn models and
    of the corpus's 6,223 free shapes captured before any kernel code
    accepted them (`9d7e5dfd`), with the independent validator and reader
    extended.
  * `Topology::class()`, the validator's two-sided faces, wire edges and
    acorn vertices, OCCT's counts (a free face has no shell, a free edge no
    wire, every ring edge one vertex, which also fixed a disc face's and a
    circle wire's count found by the `identity` fuzz target), certified
    area and length with centre (`measure_enclosure`), `Body` with
    `face_from_profile_with` and `wire_from_boundary_with` (`MakeFace`,
    `MakeWire`; roles `Face`, `Edge`, `Vertex`), free-shape import and
    writing. `compare_brep.py --family sheet` gives 17 matches, 1 reviewed
    difference, 0 failures; `compare_brep_io.py` imports 6,188 free shapes
    and reads every one back natively valid, and the eleven builder bodies;
    the DRAW adapter restores free shapes and runs `plane`, `cylinder`,
    `line`, `circle`, `mkface`, `mkedge` (derived case `faces_and_edges`,
    both backends).
  * Decisions taken in implementation, none of them a user decision: a free
    shell is closed when every edge with a curve is used exactly twice; the
    outer loop of a spline face, as of a plane, is the one enclosing the
    largest UV area (a hammer face stored its hole first); a spline over
    `MAX_POLES` is `BSplineControlDataLimit`, and the independent reader
    models that documented limit (terrain, MAT); the hammer's 31 faces the
    validator rejects though BRepCheck accepts them (C0 spline surfaces,
    certified gaps above OCCT's stored tolerance, undecided uses) are pinned
    in `occt_brep.rs`, and twelve restore-only DRAW cases whose free faces
    the validator rejects are failures (`UPSTREAM_TESTS.md`); spline wires
    are not measured.
* S7 — in progress.
  * D13 recorded in `TOPOLOGY_MODEL.md`; the S7 decisions above.
  * S7a implemented (`ANALYTIC_INTERSECTIONS.md`): the exact reference and
    the native `IntAna_QuadQuadGeo` capture of 67 cases came before any
    kernel intersection code (`943463ea`); the kernel's closed forms contain
    the reference on all 67; the bridge gives 63 matches and 4 reviewed
    differences (IntAna's snapping of near-degenerate configurations); the
    `analytic_intersections` fuzz target. Gate pending CI, the schedule
    replay and the clean campaign.
  * S7b.1 implemented: the exact reference and the native `GeomInt_IntSS`
    capture of 20 cylinder/cylinder and cylinder/sphere pairs came before
    any kernel procedural code (`0893aec3`); `ProceduralCurve` with loops,
    rings and figure-eights, certified points on both branches; the kernel
    inside the reference on all 20; the bridge gives 17 matches and 3
    reviewed differences (missed tangent points, one approximated sample);
    the fuzz target samples every component. Gate pending CI, the schedule
    replay and the clean campaign.
  * S7b.2 implemented: the extended reference and a second native capture
    of 15 cylinder/cone and sphere/cone pairs came before any kernel cone
    code (`9d764085`); the sphere/cone loops in closed form, the
    cylinder/cone roots by certified subdivision (mean-value enclosure,
    budget, per-piece rational fallback, bisection to a few ulps); the
    kernel inside the reference on all 35; the bridge gives 31 matches and 4
    reviewed differences. `analytic_intersections` keeps allocation stacks to
    five frames (AddressSanitizer's stack depot otherwise outgrew the 2 GB
    gate). Two cones and a cone whose rational apex is on a sphere stay
    `NotConic`.
  * S7b.3a implemented: the reference extended with the meridian
    parameterisation and a third native capture of 58 torus cases came
    before any kernel torus code; `intersection/toroidal.rs` decides every
    class of a torus and a plane or a sphere exactly (the signs of a
    quadratic at `+-m`, numbers in `Q(sqrt q)`), gives loop ends in closed
    form and the special and coaxial pairs as circles; the kernel inside the
    reference on all 93, the bridge 86 matches and 7 reviewed differences
    (three more missed tangent points). A near-Villarceau case exposed that
    the platform's `hypot` can store a non-Pythagorean normal one unit in the
    last place away from the reference's emulation: the procedural fixture
    now checks stored normals bit for bit, as S7a's does.
  * S7b.3b.1 implemented: the independent reference and a native capture of
    24 torus/cylinder and torus/cone pairs came before any kernel code;
    `intersection/torus_curves.rs` traces the curve as a certified graph
    (folds by subdivision and Krawczyk, branches by certified windows,
    components with winding numbers) and `tangency.rs` decides the
    tangencies exactly from the pipes' spine and axis; the kernel inside the
    reference on all 24, the bridge 18 matches and 6 reviewed differences.
    The fuzz target found a track's end refused after rounding and a point
    evaluation too slow in rational intervals (interval Newton now); both
    are checked-in regressions. Linux CI's scheduled run later timed out
    (42 s) where a fold box's top or bottom edge was crossed by a twin
    fold's branch: the edge subdivision descended to its floor there in
    binary64 and again in rational intervals; a certain sign change of `G`
    between exact points of the edge now reports the crossing at once and
    skips the rational retry, every certified result unchanged (2.3 s,
    `fuzz/regressions/README.md`).
  * S7b.3b.2 implemented: the reference extended with the other torus's
    quartic and both spines' critical pairs, and a native capture of 14
    pairs of tori, came before the kernel intersected two tori; the field's
    jets carry the quartic and `tangency::torus_torus` reduces the
    criticality and distance conditions modulo the first spine's circle; the
    kernel inside the reference on all 38 torus cases, the bridge 28
    matches and 10 reviewed differences.
  * S7b.4 implemented: the projective ruling reference and a native
    capture of 12 pairs came before `intersection/ruled_curves.rs`; the
    traced machinery takes a chart (a torus's meridians, a cone's rulings,
    the apex and twin factors), components count their crossings of
    infinity, the apex is a node, two cones with one apex give their common
    generatrices; the kernel inside the reference on all 12, the bridge 9
    matches and 3 reviewed differences. S7a's `kk_crossing` (two congruent
    cones crossing symmetrically, two conics meeting at nodes) is
    `ComputationLimit`, which its comparison accepts.
  * S7c.1 implemented: the exact sympy reference and a `GeomAPI_IntCS`
    capture of 36 lines and circles against every analytic surface came
    before `intersection/curve_surface.rs`; lines exactly by their
    polynomials, circles by the resultant modulo their conic, both with
    exact multiplicities (tangencies) and containment, a cone's in certified
    intervals; the kernel inside the reference on all 36, the bridge 29
    matches and 7 reviewed differences (a contained curve never a native
    segment); fuzz target `curve_surface`.
  * S7c.2 implemented: the reference (rational parameterisations, exact
    span polynomials) and a `GeomAPI_IntCS` capture of 35 ellipses,
    hyperbolas and splines came before `intersection/conic_surface.rs` and
    `spline_revolved.rs`; conics by S7c.1's resultant (shared as
    `section`), a hyperbola against a cone as a quartic in `e^t`, splines
    against tori exactly and against cones with exact overlaps and apexes
    (`tan a` transcendental); the kernel inside the reference on all 71,
    the bridge 59 matches and 12 reviewed differences. A tilted conic's
    stored `y` showed that `frame_axes` omits `Vec3::normalized`'s scaling:
    the S7c.2 reference emulates it and the fixtures record each conic's
    stored axes (the shared helper stays, its captured inputs unchanged).
  * S7d.1 implemented: the Groebner-basis reference and an
    `IntTools_EdgeEdge` capture of 41 pairs came before
    `intersection/curve_curve.rs`; lines by linear algebra, a line and a
    conic by its plane and equation, conics by S7c's resultant in one plane
    or the gcd along two planes' line; the kernel inside the reference on
    all 41, the bridge 41 matches; fuzz target `curve_curve`.
  * S7d.2 implemented: the span-polynomial reference and an
    `IntTools_EdgeEdge` capture of 13 splines against conics came before
    `intersection/spline_curve.rs`; overlaps, points and tangencies from the
    gcd of the conic's plane and equation on each span; the kernel inside
    the reference on all 13, the bridge 11 matches and 2 reviewed
    differences; `curve_curve` fuzzes splines too.
  * The `lowalgos` group: the Rust DRAW adapter gained DRAW's analytic
    surfaces (`sphere`, `cone`, `torus`), `intersect` (circles, ellipses,
    lines and hyperbolas' branches as DRAW curves, a procedural loop, ring
    or figure-eight and a traced track as parameterised curves), numeric
    variables (`bounds`, `dval`, `renamevar`, `directory`), `dump`,
    geometry-aware `whatis` and `xdistcs` with OCCT's output. Four
    self-contained `intss` cases are registered: three evaluated on both
    backends (a torus with a plane through its axis, twice; a cone with a
    coaxial torus) and one declared gap (antiparallel axes OCCT snaps by
    tolerance); three more run but stay unregistered, the contract
    admitting no failure: they count OCCT's pieces of walking lines
    against the kernel's closed loops. The other self-contained cases need trimmed and
    extruded surfaces, and the data cases B-spline surfaces, beyond S7's
    analytic scope (UPSTREAM_TESTS.md).
* Tessellation (parallel track) — T-a implemented (`TESSELLATION.md`);
  gate pending CI, the schedule replay and the clean campaign.
  * The decisions above, the independent reference
    (`tessellation_reference.py`) and fixtures of 28 bodies at two settings,
    and the native `BRepMesh_IncrementalMesh` capture of all 56 meshes came
    before any kernel tessellation code (`1c4db5a3`); the pinned SDK now
    builds TKMesh (`build_pinned_occt.py --toolkit TKMesh`, also in the
    B-rep oracle job).
  * `tessellation.rs` (edges, charts, faces, the contract check),
    `tessellation/cdt.rs` (constrained Delaunay triangulation with exact
    orientation predicates) and `tessellation/bounds.rs` (certified bounds
    in the `Fast` tier). All 56 fixture meshes pass every check of the
    reference; the 54 certified `data/occ` solids without spline geometry
    mesh closed and oriented with their volumes in the deflection band of
    their mass enclosures; the bridge gives 42 matches and 14 reviewed
    differences (OCCT beyond the request on seven rows, degenerate triangles
    at every apex and pole).
  * Decisions taken in implementation, none of them a user decision: the
    certified map is chosen per triangle among three (linear, tangentially
    corrected, and the fan map at a pole; `MATHEMATICS.md`), added after the
    first measurements showed apex cones at up to 15 times BRepMesh's
    triangles with the linear map alone; a triangle is split at the middle
    of the edge contributing most to its bounds, not its longest in space
    (the longest drove points into an apex without reducing the normal
    turn); the reference's normal test offsets a centroid by twice the larger
    of the request and the triangle's own deviation; native rows reproduce
    with exact counts on every platform (`VALIDATION.md`'s allowance table),
    so a Linux BRepMesh that triangulates differently needs its own capture
    review.
  * Clean local 600-second `tessellation` campaign at `c9dbf4f6`
    (AddressSanitizer, standard limits): 4,855 mutation executions after
    209 s of replay (350 inputs, 14,000 edges, 961 MB peak), no artifact.
  * Open for the user: U9 (a display mode without a bound). T-c
    (procedural edges) pending.
  * T-b implemented (`TESSELLATION.md`, splines); gate pending CI (the
    spline bridge on macOS and Linux; a Linux BRepMesh that meshes the
    splines differently needs its own reviewed capture record), the
    schedule replay and the acceptance of the clean campaign below.
    * The decisions above (`e065cd45`), then the reference's spline bodies
      and twelve spline fixtures at two settings with `--check`
      (`394f1491`), then the native `BRepMesh_IncrementalMesh` capture of
      all 24 meshes (`130cb8b3`), all before any kernel spline
      tessellation code.
    * `tessellation/spline.rs` (cells, their bounds, cones and rates) and
      `tessellation/spline/poly.rs` (tensor Bernstein polynomials in the
      `Fast` tier). All 24 spline fixture meshes pass every check of the
      reference; the bridge gives 21 matches and 3 reviewed differences
      (OCCT beyond the request on the dome and the trimmed sheet at the
      fine setting, its recorded deflection understating the dome's); the
      17 valid B-rep fixture cases with spline geometry and all 58
      certified `data/occ` solids mesh within the contract; the T-a
      bridge is unchanged (56 of 56, 42 matches, 14 reviewed).
    * Decisions taken in implementation, none of them a user decision:
      the first measurements (the rational quarter-circle wall at 34 times
      BRepMesh's nodes, the sharp rational corner at 7,900) led to three
      sharper certificates kept beside the recorded bounds, the smallest
      winning (`MATHEMATICS.md`, numerators): rational cells of low degree
      bound their derivatives through the Bernstein coefficients of the
      quotient rule's numerators, surface cells carry the cone of their
      normal numerator's coefficients, and normal and tangent turns use the
      rates `|M × M_u| / |M|²` per cell; surface cells are halved three
      times per direction rather than twice (the rational sheet's nodes
      fell from 9,402 to 5,387); a spline solid's volume is not checked
      against its certified mass enclosure in the corpus test (minutes per
      solid); the tests' and the fuzz target's distances to spline surfaces
      are projections with their own binary64 Cox-de Boor evaluation from
      the four nearest points of a knot-aware grid (a projection from one
      point of a 33 × 33 grid found the wrong basin on a degree-8 corpus
      surface, where the kernel's own parametric claim held at every
      centroid, checked with its exact evaluation); the fuzz target's
      weights stay within 0.7 to 1.5, its spline deflections at or above a
      64th of the size and angles from 0.35 rad, its spline prisms on
      axis-aligned frames, and a spline prism's volume is checked against a
      quadrature of its profile: uneven weights inflate the rational bounds
      by the weights' ratio to a high power, and the certified mass
      enclosure, the validator's fluxes on tilted spline walls and the
      kernel's exact evaluation cost seconds per input under the sanitizer
      (the first campaign's two slow units, retained as regressions).
    * Clean local 600-second `tessellation` campaign at `d638bd0b`
      (AddressSanitizer, standard limits): 1,217 mutation executions after
      118 s of replay (251 inputs, 19,511 edges, 1,265 MB peak), no
      artifact. The first campaign, before the fuzz target's changes above,
      saved two slow units (retained as regressions,
      `fuzz/regressions/README.md`).
    * `Motor-c.brep` 378 (certified since S8d.2, pinned at the merge as a
      `ComputationLimit`) meshes: its four spline blends (rational by
      weights `1 - 2e-15`, degrees 2 and 5) have cusp corners, where both
      boundary rows leave `(1, 0)` along `-x` (the controls exactly
      collinear), so `S_u × S_v` vanishes there and the normal's limit
      turns 0.45 rad between the two rows. Every turn bound divided by a
      least `|M|` that is zero on a box at the corner (the blends' degree
      also leaves them without cones and rates), so the first segment of
      the blend's edge along `u = 1` stayed infinite at every count and
      twelve doublings (4,096 segments) ended the edge loop; skipping the
      check only moved the failure to the face's refinement. Not a bound
      that fails to shrink but none at all: the thin-triangle check and the
      triangles now also take the corner turn (`MATHEMATICS.md`, cusp
      corners), a second-order Taylor bound of `M` at the corner on its
      patch and three nested boxes; a first-order one (`G y` alone) reached
      only `v < 0.0057` of the patch's `0.0401` and left the edge at 0.2254
      against 0.225 after 4,096 segments. The solid now meshes in about
      11 s (release): 101,696 triangles, the cusp edges at 2,048 and 1,024
      segments (at 1,024, the segment straddling the end of the corner's
      patch, `V / D` = 0.024 from the corner, turns 0.253 under the
      Lipschitz bound), within twelve doublings. The corpus test asserts all 59 certified solids; its
      spline projection adds a compass search where Gauss-Newton stalls at
      the cusp (a centroid at 4.6e-4 against a bound of 3.3e-5 before it).
* STEP import (parallel track) — STEP-a and STEP-b implemented
  (`VALIDATION.md`, `SOURCE_MAP.md`); STEP-a's gate pending CI, the
  schedule replay and the clean campaign.
  * The decisions above (`a5435624`), then the evidence before any importer
    code (`cab82ecb`): the Part 21 reader with its tests, 22 STEP files
    authored for the track (`generate_step_fixtures.py`: boxes in three
    schemas and three units with every orientation flag flipped, prisms, a
    void, cylinders, cones, frusta in degrees, spheres, tori, two solids, an
    open box and a hand-written syntax file), the independent reference
    (`step_reference.py`: its own Part 21 parser, OCCT's counts, closed
    forms) and the native `STEPControl_Reader` capture
    (`fixtures/occt-step-preimplementation`, the importer absent), which
    matches the reference on all 22 within `7.6e-15`. The pinned SDK now
    builds the STEP reader (`build_pinned_occt.py --toolkit TKDESTEP`, CI job
    `pinned-step-oracle`).
  * `step::import` builds OCCT's shape structure and converts it with the
    `.brep` converter; every fixture body imports, validates and contains
    the reference's measures within `1e-12` with OCCT's counts
    (`tests/step.rs`); `compare_step.py` gives 22 matches, 0 reviewed
    differences, 0 failures; the `step` fuzz target (no panic, typed errors,
    determinism, validated bodies, `.brep` round trips).
  * The first clean campaign (at `5c35dd8a`) found the `.brep` writer
    leaving a torus band's rings a period below its latitude seam; fixed
    with a regression at `517b66a8` (`fuzz/regressions/README.md`), after
    which `compare_brep_io.py` (6,834 matches) and `compare_tessellation.py`
    (42 and 14 reviewed) are unchanged. Clean local 600-second campaign at
    `517b66a8` (AddressSanitizer, standard limits): 175,774 mutation
    executions after 149 s of replay (150 seeds), 751 MB peak, no artifact.
  * Local survey of the dataset's 336 STEP files (U1): no panic, 464 bodies
    import; the constructs that stop the rest, in order, are B-spline curves
    and surfaces (577 bodies), seamless periodic faces (106), extrusion and
    revolution surfaces (165), non-straight pcurves on curved surfaces (35)
    and ellipses (11) (`VALIDATION.md`).
  * STEP-b implemented (`VALIDATION.md`, `SOURCE_MAP.md`); gate pending
    CI, the schedule replay and a clean campaign. Its decisions
    (`651da6bd`), then the evidence before its code (`46668ea1`): seven
    fixtures (a half-ellipse sheet, an obliquely cut cylinder, a B-spline
    plate and prism, a two-span patch, a trimmed patch with a quadratic
    B-spline pcurve, a rational quarter-arc cylinder), the reference's
    closed forms, exact integrals and 40-digit quadratures with its own
    geometry check of every file (gaps within `1.1e-16` of the size), and
    the native capture (`fixtures/occt-step-b-preimplementation`, taken
    while `step/spline.rs` did not exist and the importer rejected all
    seven): OCCT's counts and validity agree; its default `BRepGProp`
    misses three bodies' measures, by up to `2.0e-3`.
  * The importer (`be1cb610`, `11979a6b`): ellipses (plane pcurves as the
    ellipse or its exact projection, the cylinder section's exact
    sinusoid), knotted B-spline curves and surfaces, simple and rational
    complex, and the file's pcurves on spline surfaces over the ranges where
    their images meet the vertices. 29 of the 30 fixture bodies import,
    validate and contain the reference's measures within `1e-12` with
    OCCT's counts (a one-face sheet counted as a free face); the rational
    cylinder is refused by R4's C1 rule. `compare_step.py`: 23 matches,
    6 reviewed differences, 0 failures. The `step` fuzz target covers the
    new fixtures; its first smoke run found the `.brep` writer's line
    pcurves on spline surfaces written at unit speed (fixed, with a
    regression).
  * Local survey (U1, counts only): 550 bodies of 131 files import (464 of
    114 before), 760 rejected (`PCurveNotDerived` 458 now leads), 61 fail
    validation (21 through R4's C1 rule on rational quarter arcs, 21 with a
    file's pcurve off its edge's parameter), no panic.
  * Open from STEP-b: (1) rational B-splines with knots of multiplicity equal
    to the degree (the usual NURBS circle) are refused by R4's homogeneous
    C1 test though C1 as rational curves; a rational C1 criterion, or
    splitting such edges and faces at those knots, needs its own decision.
    (2) Pcurves on cylinders, cones, spheres and tori beyond rulings,
    parallels, meridians and cylinder sections (B-spline and other elliptic
    edges there: most of the 458 `PCurveNotDerived`), whether from the
    file's `PCURVE`s placed on the cover or derived as projections (D13).
    (3) A file's pcurve not parameterised like its edge (21 bodies) would
    need an exact reparameterisation; OCCT approximates with
    `SameParameter`. (4) The knotless B-spline forms (25 bodies), the
    `.brep` reader's ellipse records, and a one-face sheet's shell in the
    synthesized counts.
  * Next: seamless periodic faces (windings from the derived pcurves), the
    swept surfaces, the open STEP-b items above, and placements; STEP-c the
    recorded corpus survey.
* S8 — in progress: decisions recorded (2026-09-28).
  * S8a.1 implemented: the quadrature reference and a
    `BRepAlgoAPI_Splitter` capture of 26 prisms came before
    `solid/split.rs`; planes normal to a prism's axis through M3's height
    split, parallel ones through the profile's exact section by the line and
    prisms renamed by provenance (`PlaneSplit`); the kernel inside the
    reference on the 13 such cases, two reviewed count differences (OCCT's
    tangent split and seam); fuzz target `split`.
  * S8a.2 implemented (`solid/split/oblique.rs`): a prism split by a plane
    oblique to its axis, each piece a footprint (the profile's exact section
    by the plane's trace on one cap's plane) with a creased end (its section
    by the trace on the other), built as a general body and validated; new
    kinds `Curve3::EllipseArc`, `Curve2::EllipseArc` (axis-aligned on a
    plane) and `Curve2::Sinusoid` (a plane's section on a cylinder) through
    validation (the harmonic bound certifies an ellipse edge against both
    its pcurves), mass properties, tessellation, `.brep` I/O (OCCT's
    ellipse records; a sinusoid pcurve is `Unwritable`, OCCT having no
    analytic record for it), history (`Family::Ellipse`) and the
    intersections (an ellipse edge's whole ellipse). The kernel inside the
    reference on all 26 cases, 23 native matches and 3 reviewed (a seam
    through the touch point where a plane meets both caps' circles).
    Decisions taken in implementation, none a user decision: a plane
    touching a cap's circle keeps a vertex where the wall's height vanishes
    (a non-winding loop round the cylinder); one touching a cap's arc edge
    between its ends, or whose traces on the two caps lie within the
    resolution of each other (parallel to the axis to binary64, as planes
    "parallel" to a tilted frame's axis are), is `Degenerate`; pieces are
    `Construction::Clipped` (profile, plane in the frame, index, footprint)
    so a rigid motion rebuilds them and they classify points; the mesher
    accepts a vertex twice on one loop when the two uses lie a period apart
    on the cover.
  * S8a.2 accepted locally: a clean 600-second `split` campaign at
    `c5a96c53` (AddressSanitizer: 7,994 mutation executions after 286 s of
    replay, 1,249 MB peak, no artifact) after three fixes it drove (a hole
    grazed within the resolution, a sliver piece, a winding loop's closing
    chord the binary64 mass tier could not sign) and the speed-ups they
    needed (an arc's crossings ordered by orientation signs instead of
    exact arctangents, a binary64 certified arctangent, a rebuild making
    only its own piece).
  * S8c.1 implemented (`solid/split/revolved.rs`): cones, frusta and zones
    by planes normal to their axis (within a quarter of the resolution),
    whole spheres by any plane, pieces the same primitives named by
    provenance; inside the reference on 16 of the 23 primitive cases, counts
    OCCT's. Clean 600-second `split` campaign at `8f0dc704` (4,508
    mutation executions after 331 s of replay, 1,270 MB peak, no artifact).
  * S8c.2 implemented (`solid/split/meridian.rs`): a cone, frustum or zone
    by a plane containing its axis, each half a general body on the
    input's own surfaces (`Construction::Half`); inside the reference on
    the four meridian cases, three reviewed count differences (OCCT's seam
    at angle 0 inside one half). The conic and circle sections (three
    cases) wait for S8d. Clean 600-second `split` campaign at `b789009b`
    (5,148 executions after a 571 s replay slowed by the parallel tracks'
    builds, 1,009 MB peak, no artifact).
  * S8d: decisions recorded (2026-09-28); the reference and a
    `BRepPrimAPI_MakeTorus` capture of 12 whole tori split by planes came
    before any kernel code (`2815ae12`), OCCT's spiric pieces reviewed
    (B-spline sections, 3e-7 to 9e-7 from the reference).
  * S8d.1 implemented (`solid/split/torus.rs`): whole tori by planes normal
    to their axis (tube bands and annuli as general bodies) or containing it
    (half-turn wedges); inside the reference on the 8 such cases, six
    reviewed count differences (OCCT's seams, and its split along the
    circle where a plane touches the tube's top). S8d.2 (D13's engine) next.
  * S8d.2 implemented. D13's engine (`47c27792`): `HyperbolaArc` and
    `ParabolaArc` edges, `Projection` pcurves with certified Taylor jets
    (`jet.rs`, `topology/validate/projection.rs`) through validation, mass,
    tessellation and interop. The builder (`solid/split/conic.rs`): cones,
    frusta, zones and caps by planes neither normal to nor containing their
    axis, as general bodies (`Construction::Half`). Inside the reference on
    S8c's three conic cases and 12 more (`split-conic-cases.txt`); three
    match OCCT, twelve reviewed count differences (its seam at `u = 0`),
    one also BRepGProp's error on a cap. Amendments to the decisions, made
    while building it:
    - A plane touching one rim while crossing the other, or touching a rim
      with a zone's side circle, pinches a wall mid-loop and is
      `Degenerate` (as S8a.2's tangent arc); a closed section touching a
      rim keeps a vertex there (as S8a.2's touch). A plane through an apex
      or pole off the axis is `OutOfDomain`; one within the resolution of a
      frustum's virtual apex cuts its rulings.
    - A zone's circle not round the axis leaves a hole in the other piece's
      band: the validator now certifies an unwound loop in a wound face by
      a signed `+v` ray over the other loops' lines (`BREP_VALIDATION.md`),
      which also certifies `Motor-c.brep` 378 (59 corpus solids certified).
    - The quadrature's order is 14 (interval jets overestimate their high
      coefficients about fourfold per order through `atan2`'s quotient);
      mass moments share one evaluation per piece and take their width
      relative to their size, sign decisions keep it absolute.
    - The 12 added cases' native capture came after `conic.rs` existed and
      is recorded as such (`occt-split-conic-postimplementation/NOTES.md`);
      the capture before any code is S8c's.
    - The `split` fuzz target's per-input limit is 60 s: a cap cut near its
      pole takes 8.8 s under ASan, the quadrature refining where the
      section's angle about the axis turns fast. A parametrization without
      that refinement is a follow-up.
    Two campaigns found a cap's rim crossed within rounding of its tangent
    (`2f60b19f`) and a prism's trace passing a U's vertices within rounding
    (`a4e1c9df`), both now `Degenerate` with their regressions; the third,
    at `a4e1c9df`, was clean: 2,051 mutation executions after a 995 s
    replay, 20,798 edges, 1,903 MB peak (the target now purges the
    allocator as the other heavy ones do), one slow unit (10 s under
    AddressSanitizer, kept as a regression). The conic capture's Linux
    record: CI run 36416172849 (`platform-linux/`).
  * S8d.3 implemented (`solid/split/spiric.rs`): whole tori by any plane in
    spiric sections, `Curve3::Section` graphs over either angle (bands,
    caps with four analytic edges, C-shaped pieces), the torus less discs in
    validation and mass, `v`-wound projection integrals, projection ray
    crossings. Inside the reference on all 17 spiric cases (13 new, captured
    natively before the code, and S8d's four); 17 reviewed differences
    (OCCT's B-spline sections and seams). A section on its own torus
    evaluates by its angles (twice as fast); the quadrature's order is 12.
    The first campaign replayed a torus whose section's angle sat on
    `atan2`'s branch cut (fixed at `9c21cd2a` with its regression); the
    second, at `9c21cd2a`, was clean: 1,787 mutation executions after a
    948 s replay, 21,482 edges, 1,255 MB peak, two slow units kept as
    regressions. The spiric capture's Linux record: CI run 36556520655
    (volumes within 2.1e-12 relative, `platform-linux/`).
  * S8b: decisions recorded (2026-09-28); its evidence came before any
    kernel code (`cb701edc`): the case protocol's spline segment (`B` in an
    `S` row, reversed with its path), the reference slicing spline profiles
    by exact Bezier pieces (straight splines give the polygon's rows, and
    every profile of lines and splines its Green's-theorem moments in exact
    Fractions, within 1e-40), 17 prisms (`split-spline-cases.txt`) and
    `occt-split-spline-preimplementation`: every piece valid, the sides the
    reference's, 10 within 2e-8 and seven reviewed BRepGProp errors (walls
    of one span cut across their rulings, 3.6e-8 to 2.6e-7; faces bounded by
    a three-span spline, up to 1.8e-3, though Green's theorem over OCCT's
    own edges gives the reference's area). The probe cannot read `B` rows
    yet: all 17 `rust_unsupported`. S8b.1 next.
  * S8b.1 and S8b.2 implemented (`decide/splines.rs`, `profile.rs`,
    `topology.rs`): spline profile segments screened for separation (and
    each against itself in pieces turning less than a half-turn), exact
    Bernstein moments, prisms with exact degree-(p, 1) spline walls.
  * S8b.3 implemented (`solid/split/spline.rs`, `oblique.rs`): splits by
    exact spline/line roots, restrictions by exact knot insertion, creases
    on the plane as affine images with exact wall pcurves. Inside the
    reference on all 17 spline cases (103 in all: 52 matches, 51 reviewed,
    one new review: OCCT splits the bulge's cap edges where a plane touches
    them). Pending: the campaign. The spline capture's Linux record: CI run
    36556520655 (`platform-linux/`).
  * S8e: decisions recorded (2026-09-28); its evidence came before any
    kernel code (`5fe8d268`): face and wire bodies in the case protocol
    (`make`, natively a `make` row in place of the prism vector), the
    reference's `planar_rows` (a sheet's side: area, perimeter and centre
    from the profile's section; a wire's: length and centre from its
    boundary cut where the trace crosses or touches it, a piece along the
    trace on the side of the piece before it; rows `side S area perimeter
    cx cy cz` and `side S length 0 cx cy cz`), checked within 1e-38 against
    closed forms, exact clipping, Green's theorem and the sides' sums, 25
    cases (`split-sheet-cases.txt`) and `occt-split-sheet-preimplementation`:
    every piece valid, the sides and wires' runs the reference's, 23 within
    2e-8 and two reviewed BRepGProp errors (a face bounded by two spans of
    the wave, 1.3e-3; `LinearProperties` on the lens's cubic, 3.5e-6). OCCT
    keeps a split wire one wire (the probe groups its runs), splits edges
    at tangencies, keeps circles' seam vertices and splits the sheet whose
    hole the plane touches (the decisions: `Degenerate`). The probe cannot
    build a body case yet: all 25 `rust_unsupported` (128 cases: 75
    matches, 53 reviewed). `Body::split_by_plane` next.
  * S8e implemented (`body/split.rs`, `Topology::open_wire`): sheets by
    their profiles' sections, closed wires into open wires of their runs,
    renamed by provenance; a spline wire's length by certified quadrature;
    a certified perimeter (`Topology::edge_length_enclosure`). All 128 split
    cases inside the reference (70 matches, 58 reviewed); five new reviews:
    OCCT splits edges at the tangent sheet's arc and the tangent wire's
    spline, keeps circles' seam vertices (a disc and a circle wire), and
    splits the sheet whose hole the plane touches, which the kernel refuses
    (`refused`, a documented `Degenerate`, is now a probe row the
    comparison reviews). The DRAW adapter drives `bsplit` with one plane
    face on prisms, sheets and wires (the derived case `split_plane`, on
    both backends) and the upstream `bsplit` group is registered as
    capability sentinels. Amendments: a wire's certified length is its
    stored arcs' `r |sweep|` exactly, whose new vertices are the exact
    split's rounded, so its comparison allows `2^-40` of the case's size;
    a sheet's or wire's arcs and circles take the face frame's axes bit for
    bit (`Frame3::at`: normalizing a unit vector again is not idempotent,
    and the history check compares circle normals exactly). The first S8b
    campaign, at `66c32112`, found no crash but a 66 s slow unit (a lens
    hole cut obliquely, fixed at `55bb0516`); the campaign at `5c50de15`
    (spline prisms, sheets and wires in the target) was clean: 1,787
    mutation executions after a 1,073 s replay, 30,505 edges, 2,020 MB peak
    (the target's RSS cap is 2,048 MB), three slow units kept as
    regressions, the slowest 18 s under AddressSanitizer. The spline and
    sheet captures' Linux records: CI run 36556520655 (`platform-linux/`;
    the sheet's two mirror pieces of `wire_wave_four` in the other order).
* S9 — in progress:
  * S9a: decisions recorded (2026-09-28); its evidence came before any
    kernel code (`71eb169a`, `BOOLEAN.md`): Boolean cases in the case
    protocol (the object's prism rows, `boolean fuse|cut|common ID`, the
    tool's rows; natively the two constructions joined by a `boolean` row),
    the reference `boolean_reference.py` (the tool moved into the object's
    frame by its exactly solved, binary64 offset; each slab's region a union
    of the atoms inside both profiles, the object's only and the tool's
    only, their areas and moments by exact slicing in closed form, their
    boundaries by pieces classified against the other profile; walls, caps
    as symmetric differences of consecutive slabs, solids by union-find,
    touching regions separate), checked within 1.2e-38 against Green's
    theorem, `fuse = A + B - common` and `cut = A - common`, exact Fraction
    clipping and boundary classes of polygons, closed-form lenses and hand
    results, 45 cases (`boolean-cases.txt`: 15 fuses, 17 cuts, 13 commons;
    every class of the decisions in both frames, 13 in the tilted one; 31
    prisms, 6 S9a.2 stacks, 5 empty, 3 degenerate)
    and `occt-boolean-preimplementation`: every `BRepAlgoAPI_Fuse`/`Cut`/
    `Common` result valid with the reference's solid count, all 45 within
    2e-8 (6.1e-15 at worst), no review. Touching conventions (the
    regularized Boolean, as OCCT): a shared wall fuses into one solid, a
    common or cut leaving only a face, an edge or a point is empty. OCCT
    differs from the decisions where the result touches itself: two valid
    solids for prisms touching along a vertical edge, one for a hole
    tangent to the outer circle (the fixtures expect `Degenerate`); and it
    merges the inputs' collinear edges when unified, which the decisions
    leave open for the traced profile. No probe yet: all 45
    `rust_unsupported`. S9a.1 next.
  * S9a.1 implemented (`profile/boolean.rs`, `solid/boolean.rs`,
    `OperationKind::{Fuse, Cut, Common}` coded 11 to 13): all 45 fixtures as
    the reference and OCCT's unified counts (36 results, 5 of them empty, 3
    refused, 6 stacks `OutOfDomain`; 45 matches, no review); the `boolean`
    fuzz target. Amendments, from its implementation: (a) a result's
    boundary keeps no vertex where it does not turn: consecutive collinear
    lines and arcs of one circle in one sense are joined into one segment
    continuing every input segment it holds (OCCT's unified result; the
    evidence left it open), a cycle left as one arc its whole circle; (b) a
    tangency cuts nothing (the pieces on either side lie on one side of the
    other boundary; a result touching itself there is refused when traced or
    validated), and a circle cut at one point only stays whole; (c) a cut's
    tool faces the other way where it bounds the result, so its walls,
    edges and vertices there are `Generated` from the tool's, not continued
    (a `Split` child shares its parent's orientation); (d) several inputs
    reaching several results have no single relation: each such result is
    `Generated` from its parents and those inputs are `Deleted`; (e) a fuse
    of meeting ranges is one prism when the profiles are identical, the
    input holding the other (in 2D and in height), or both inputs when the
    profiles are apart; (f) `decide::arcs` decides an arc's side of a
    direction per orientation (one certain side settles it) and a line's
    crossing of a circle is off an arc whose side excludes its direction
    wherever it lies along the line (a profile's line through an arc's
    centre from a point just off its circle was refused as touching it);
    (g) inputs sharing an entity id (built by one operation, one solid
    twice, or a result with an input whose entities it keeps) are refused
    (`InvalidLabel`): the history names each input's entities by id; a
    result, renamed, is an input again (its entities read off the prism
    built afresh from its profile). The DRAW adapter runs `bfuse`, `bcut`,
    `bcommon`, `btuc`, `bop` with its operations, `bbop` and `bapibop`
    (`UPSTREAM_TESTS.md`): the derived case `boolean_prisms` and 320
    upstream cases evaluate on both backends, none fails; the other 1,258
    native ones wait for S9a.2 (163), S9b (962) or are refused or not
    adapted. Campaigns at `d7e515d9`: `split` clean (1,951 mutation
    executions after a 1,351 s replay of 3,006 inputs, 30,519 edges,
    1,931 MB peak); `boolean` stopped at the 2 GiB gate on an input that
    peaks at 6 MB alone, so the target joined the allocator-purge targets
    (`FUZZING.md`; `split`, listed there, now calls the purge too). An
    earlier `split` campaign found a zone's split rim one ulp off its
    circle (the latitude's sine fused with its cosine into one `sincos`):
    ring heights and radii now come from out-of-line `scaled_sin` and
    `scaled_cos` (`fuzz/regressions/README.md`). The capture's Linux
    record: CI run 36556520655 (within 1.1e-15 relative, `platform-linux/`).
  * S9a.2's stacks implemented (`solid/boolean/stack.rs`,
    `Construction::Stack`): all 45 fixtures as the reference and OCCT's
    unified counts (42 results, 6 of them stacks, 3 refused; 45 matches,
    no review); `tests/booleans.rs` hand stacks (a tower, a pocket, a plug
    filling a hole over part of its height, a cavity, a tool through a
    round wall); the `boolean` target's heights inside and on top of the
    object, its 661-input corpus replaying clean with the history check.
    Amendments, from the implementation: (a) one arrangement of both
    profiles serves every slab (each piece knows which operands lie on its
    left and right; a slab's region is a set function of the two), where
    the decision named each slab boundary's regions; the caps between slabs
    are traced from it as the upward and downward differences; (b) walls
    on one line or circle facing one way join across slab heights and
    piece ends where nothing else meets them, and edges keep no vertex
    where they run straight on between the same faces (OCCT's unified
    result, as S9a.1's amendment (a)); (c) a closed cavity (a tool inside
    the object in 2D and in height) is a second shell of the solid and a
    bounded void generated from the tool's region; one beside several
    solids is `OutOfDomain`; (d) a result touching itself along an edge
    (four faces) or a face meeting itself is `Degenerate`; (e) a horizontal
    face continues the caps at its height facing its way whose region it
    overlaps, decided by an exact 2D common; (f) the order's spline
    profiles (S9 Order: "S9a.2: spline profiles") are S9a.2's second
    part, after the stacks; (g) a stack's rigid motion rebuilds it from its
    construction in the moved frame (as S8's pieces), keeping its ids and
    moving its mass, rather than moving its stored geometry. The DRAW
    adapter runs stacks (unified as built; a stack given to another
    Boolean unsupported until S9b) and the derived case `boolean_stacks`
    (a step, a pocket, a slab cutting a box in two, a closed cavity, a tool
    through a round wall) on both backends; of the 163 stacked upstream
    cases 84 now evaluate (404 of the `boolean` group in all, none failing),
    37 give a stack to a second Boolean and 38 are `Degenerate`. The survey
    found the fuse's containment test propagating a degenerate cut
    (`bopfuse_simple/Z5`, a box inscribed in a cylinder): a refused
    selection now counts as not empty. Profiling a slow fuzz input moved
    certified bounds to the leading-bits bracket before the binary64
    pattern search and a stack cap's overlaps to its pieces' memberships
    (the input six times faster). Pending: the campaign; S9a.2's spline
    profiles next. Their evidence exists before their code (`3951925a`,
    `BOOLEAN.md`): 46 spline cases (36 prisms, 5 stacks, 4 empty, 1
    degenerate; 18 tilted), the reference checked within 3.7e-39 (chords
    2.1e-4), and `occt-boolean-spline-preimplementation`: all valid with the
    reference's solids, 37 within 2e-8, 9 reviewed (BRepGProp on B-spline
    faces, up to 8.0e-4), the kernel `unsupported` on all 46.
  * S9a.2's spline profiles implemented (`profile/boolean/splines.rs`):
    all 46 spline fixtures as the reference (45 results, the degenerate one
    refused) and OCCT's unified counts but four reviewed tangencies (OCCT
    keeps the touching point as vertices and edges); every fixture's
    history checked; the `boolean` target's spline profiles (900 spline
    variants of its corpus replaying clean). Amendments, from the
    implementation: (a) two splines tangent where they meet are refused
    (`Degenerate`): the decisions left that tangency's side undecided and no
    fixture holds one; (b) one curve in either direction is equal poles,
    knots and range, or reversed poles with knots and range mirrored
    exactly; (c) a spline piece's end poles are set to the arrangement's
    vertex (a crossing of two splines is one point; each restriction's own
    end is within rounding of it); (d) a stack's walls on one spline
    segment lie on that segment's whole degree-`(p, 1)` wall, joined across
    its pieces' ends and slab heights; (e) `translated` moves a spline's
    poles exactly (it cloned them: an offset spline tool failed as an
    invalid curve, found by the evidence track). Pending: the campaign.
  * S9b's evidence (a parallel start, before S9a.2's splines end): the
    decisions above, then `polyhedral_reference.py` and 45 fixtures
    (`generate_polyhedral_fixtures.py --check`: S9a's slicing on its 21
    polygon cases within 2.3e-17, the area identity within 3e-41, boxes'
    closed forms exact) and the `BRepAlgoAPI` capture before the kernel
    module (`compare_polyhedral.py`: 45 matches, no review). Amendment: the
    decisions' fixtures of a stack and a plane piece as inputs wait for the
    case protocol to chain operations; the kernel's tests take them.
  * S9b.1 implemented (`solid/boolean/polyhedra.rs`,
    `Construction::Polyhedron`): prisms of line profiles in any relative
    position (and S9a's pairs whose offset or tool profile would round),
    all 45 fixtures as the reference (41 results, 4 refused as declared) and
    OCCT's unified counts but two reviewed imprints; histories checked for
    every fixture; the `boolean` target's turned, leaning and tilted tools,
    its corpus replaying clean. Amendments, from the implementation: (a) the
    decisions' input set is split: S9b.1 takes prisms, S9b.2 the Boolean's
    stacks, plane pieces and S9b results as inputs (their naming needs keys
    beyond a prism's slots); (b) a cap's exact plane is normal to `x * y`,
    not the stored `n` (the stored axes are not exactly orthogonal, so the
    lifted profile points lie on the plane of `x` and `y`); (c) faces are
    fragments of each input face split by every plane of the other's faces
    (not a 2D arrangement per face), joined back into maximal faces; (d) a
    face with a vertex within the resolution of another vertex or of an
    edge not ending there is `Degenerate` (a neck thinner than the
    resolution; the fixture `ell_tilted` fuse and cut declared so after the
    kernel met it), and two solids sharing a vertex are `Degenerate`; (e) a
    result's rigid motion moves its stored geometry (the decisions' rule for
    general bodies): rebuilding its exact model in moved frames changed
    near-coincident results. Pending: the campaign, S9b.2. DRAW survey
    (2026-09-28, `UPSTREAM_TESTS.md`): of the boolean group's 1,802
    self-contained cases Rust evaluates 704 (404 before), none failing; of
    the 647 in frames with different axes 298 evaluate (two boxes, one
    turned; the adapter turns quarter turns exactly), 233 are arcs (S9c),
    114 `Degenerate`, 2 give a polyhedron to another Boolean (S9b.2).
  * S9b.2 implemented: a Boolean's stack or polyhedral result, a plane's
    piece and any solid of planar faces and straight edges are inputs
    (`tests/polyhedral_booleans.rs`: two pockets cut in turn, a turned box
    against a stack, an S9b result against a prism, a plane's piece
    against a box, each result's history checked and the volume
    identities; the `boolean` target's first results chained against a
    turned box). Amendment, from the implementation: such an input is
    decided on its stored geometry, not a construction model (its faces'
    triangles exact, its side by exact ray parity), and its fragments join
    back by the stored face they come from; a line-profile prism keeps its
    construction model. A result's inputs are kept as solids (its rigid
    motion moves them too, for classification). A stack's rigid motion now
    moves its stored geometry too (S9b.1's amendment (e), replacing S9a.2's
    (g)): rebuilding spline stacks in the moved frame was most of the
    `cfc641c9` campaign's timeout. Replaying the corpus with the chained
    stage (debug assertions) found: (a) a tilted stack's stored model open
    where its cap's trapezoids spanned a run of collinear edges without
    their middle vertices (stored vertices are not collinear exactly):
    stored faces are now triangulated on their own vertices (ears clipped,
    the fattest first, holes bridged; the zipped trapezoids a fallback
    under an exact area check), so neighbouring faces share their stored
    edges exactly; (b) a stack touching itself at a vertex (a lens hole's
    corner on the edge line of a block standing on the box) reported as a
    non-manifold vertex (`InvalidTopology`, at `cfc641c9` too): now
    `Degenerate`, S9a's rule; (c) exact fragments of stored models are
    slow: 131 s in release for one chained input before splitting planes
    were filtered by bounding boxes, fragments classified at short
    interior points in one ray cast off the other's planes, and the ray
    casts and collinearity tests given integer forms over common
    denominators after binary64 filters; about 1 to 2 s an operation on a
    stored model of 50 to 150 triangles remains (coordinates of ~900
    bits: three stored planes' intersections), 25 times that under ASan.
    The `boolean` target chains one result of at most 12 faces (cut and
    common against the box), its per-input limit is 60 s, and Bernstein
    degree elevation no longer brackets `i / n` from rationals (the
    certified mass of spline walls, most of a spline stack's cost).
    Campaigns: at `326ad26c` `boolean` clean (1,206 mutation executions
    after the replay of 999 inputs, 324 MB peak) and `split` clean (2,356,
    3,293 inputs, 1,396 MB); at `cfc641c9` `split` (1,898, 3,396 inputs,
    1,900 MB) and `step` (13,991, 743 inputs) clean, `boolean` a timeout in
    its replay (`fuzz/regressions/README.md`). DRAW survey (2026-09-28,
    `UPSTREAM_TESTS.md`): of the boolean group's 1,802 self-contained cases
    Rust evaluates 741 (704 before), none failing; of S9b.1's 39 stacks and
    polyhedra given to another Boolean 35 evaluate (pockets cut one after
    another) and 4 are stacks with cylindrical walls (S9c);
    `bopfuse_simple/H3` and `H4`, S9b.1's directions of zero length,
    evaluate.
  * S9c: decisions recorded (2026-09-28); S9c.1's evidence came before any
    kernel code (`BOOLEAN.md`): `curved_boolean_reference.py` (slices
    parallel to both axes, breakpoints as exact polynomial roots, faces
    swept in their own parameters) and 44 fixtures
    (`generate_curved_boolean_fixtures.py --check`: 35 solid, 2 empty, 7
    degenerate with reasons; closed forms within 9.2e-41 in exact frames and
    1.9e-16 in turned ones, identities 2e-40, Monte Carlo 2.7 sigma, S9a's
    and S9b's references on their 90 fixtures) in frames the kernel stores
    bit for bit (`ROT` is not: its `x` differs in the last bit on macOS
    arm64), and the capture `occt-boolean-curved-preimplementation`
    (`compare_curved_boolean.py`: all 44 valid with the reference's solids,
    within 8.6e-9, no review; the kernel `unsupported` on all 44).
    `compare_polyhedral.py` had checked S9a's capture since the spline set;
    `compare_boolean.make_set` now lets each wrapper choose its set.
    Amendments, from the evidence: (a) a cylinder's exact model is the
    affine one on its frame's stored axes (in a turned frame slightly
    elliptic), so two cylinders are one surface only when their frames' axes
    are equal bit for bit or exactly orthonormal and their exact models
    agree; others within rounding of each other are `Degenerate`, as the
    reference refuses them; (b) two cylinders whose surfaces meet in a curve
    other than lines, circles or ellipses are S9c.2's (`OutOfDomain` until
    then) even when their bounded faces do not touch; (c) a tangency between
    the inputs is refused (`Degenerate`) whether or not the result involves
    it (an empty common, a cut leaving the object unchanged), S9a's rule;
    OCCT returns valid results there, recorded as declared refusals. S9c.1's
    kernel next.
  * S9c.1 implemented (`solid/boolean/curved/`): prisms of line, arc and
    circle profiles in any position, on exact models, as the decisions
    describe (vertices where an edge meets a face, section edges between
    them, each face's pieces traced and kept by exact membership); all 44
    fixtures as declared (38 as the reference, 6 `OutOfDomain`: S9c.2's
    cylinders), every history checked, results deterministic and moved
    rigidly; `compare_curved_boolean.py` 42 matches and 2 reviewed counts
    (OCCT's seam edges on the Steinmetz fuse's and common's faces); the
    `boolean` corpus (1,080 inputs), spline variants and regressions replay
    clean. Amendments, from the implementation: (a) an arc must end on its
    circle exactly (else `OutOfDomain`); full circles are split at rational
    points `(1 - s^2, 2 s) / (1 + s^2)`, each input at another `s`, tried
    again at others when a meeting falls on a seam; (b) evidence amendment
    (b) is narrowed: two cylinders not circular in a common measure are
    S9c.2's only when their faces' bounds meet and their sections are not
    certainly apart; (c) a plane within rounding of a cylinder's axis
    direction (its section's axis past `10^12` radii: frames built from one
    normal round it differently) is `Degenerate`, as are a vertex of one
    input on the other's face, an edge meeting an edge (off one surface),
    and edges of both overlapping on one surface; (d) faces of both inputs
    on one surface (coplanar planes, one cylinder) are taken, not refused:
    each holds the other's edges within it, their pieces facing one way
    join into one face; (e) a piece is classified at a rational point of an
    edge pushed into the piece and off the face (first-order signs, a push
    along a cylinder's circle keeping to it); loops' orientation and
    nesting from binary64 images of the face's parameters with a `1e-9`
    margin (refused within it); (f) pcurves on cylinders are lines and
    sinusoids where they fit the edge at the same fractions, exact
    projections otherwise (their curve and surface now move with a rigid
    motion); the validator decides a sinusoid's +u ray crossings, and
    `occt_counts` adds OCCT's seam vertex splitting a band loop's edge
    where no vertex lies at `u = 0`; (g) the history checker compares
    planes as oriented planes, a spline wall with its reversal in u and the
    other sense, a spline piece with a whole's reversal, and an ellipse
    edge in a plane by its centre and axes' ends. A result with arcs given
    to another Boolean remains `OutOfDomain` (S9c). The campaign at
    `e95fdfbc` (split, brep_validation, tessellation clean; boolean found
    identical prisms with a spline hole given either way round failing the
    history check, amendment (g); `fuzz/regressions/README.md`). DRAW
    survey (2026-09-28, `UPSTREAM_TESTS.md`): of the boolean group's 1,802
    self-contained cases Rust evaluates 804 (741 before; `bfuse_complex/J5`
    and 62 new, all registered); of the 234 arcs in frames with different
    axes 63 evaluate, 48 are S9c.2's cylinders, 114 `Degenerate`,
    `bfuse_simple/E1` an arc ending off its circle, and 8 fail on the
    kernel: `bopfuse_simple`, `bopcut_simple` and `bopcommon_simple` `T7`
    and `Y2` (a box's corners on the cylinder within rounding:
    `InvalidTopology("a hole outside every piece")` where `boptuc` refuses a
    piece thinner than the resolution) and `bopfuse_simple`,
    `bopcut_simple` `ZC5` (a wall with three holes of ellipse arcs: the
    validator's `uncertified_containment`), both fixed with S9c.2a. The
    campaign at `80ba3d9c` was clean: 363 mutation executions after a
    2,607 s replay of 1,081 inputs, 35,013 edges, 486 MB peak, one new slow
    unit (40 s under AddressSanitizer) kept as a regression.
