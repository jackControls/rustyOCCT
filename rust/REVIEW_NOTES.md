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

### Parallel tracks

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
    are checked-in regressions.
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
* STEP import (parallel track) — STEP-a implemented (`VALIDATION.md`,
  `SOURCE_MAP.md`); gate pending CI, the schedule replay and the clean
  campaign.
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
  * Next: STEP-b (ellipses; B-spline curves and surfaces with rational
    complex instances, the file's `PCURVE`s on spline surfaces), then
    seamless periodic faces (windings from the derived pcurves), the swept
    surfaces, and placements; STEP-c the recorded corpus survey.
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
    AddressSanitizer, kept as a regression). Pending: the Linux record of
    the conic capture after CI.
  * S8d.3 implemented (`solid/split/spiric.rs`): whole tori by any plane in
    spiric sections, `Curve3::Section` graphs over either angle (bands,
    caps with four analytic edges, C-shaped pieces), the torus less discs in
    validation and mass, `v`-wound projection integrals, projection ray
    crossings. Inside the reference on all 17 spiric cases (13 new, captured
    natively before the code, and S8d's four); 17 reviewed differences
    (OCCT's B-spline sections and seams). A section on its own torus
    evaluates by its angles (twice as fast); the quadrature's order is 12.
* S9 — pending
