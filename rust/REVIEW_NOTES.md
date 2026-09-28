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
  production-grade for measurement.
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

### S8 — general planar split and face trimming (SplitBody job)

Split any supported solid by an arbitrary plane: face/plane intersection
curves, loop splitting on the universal cover, region classification, the
`Split`, `Generated` and `Deleted` relations M3 introduced but on real
geometry. Native `BRepAlgoAPI_Splitter` bridge; upstream `bsplit` cases.
The first general topology-changing algorithm, and the rehearsal for S9.

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
* **STEP import**, after S6, reusing the converter architecture, with the
  OCCT STEP reader as the native oracle.

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
  * S7c (curve/surface) and S7d (curve/curve) pending.
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
  * Open for the user: U9 (a display mode without a bound). T-b (spline
    edges and faces) and T-c (procedural edges) pending.
* S8 — pending
* S9 — pending
