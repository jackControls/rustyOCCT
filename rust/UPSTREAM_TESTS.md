# Running original OCCT tests against Rust

The bridge evaluates the original DRAW Tcl tests through either the Rust command
adapter or native OCCT's `DRAWEXE`. It uses a real Tcl interpreter, the original
`CheckCommands.tcl` assertions and `_check_log` result parser, plus group/grid
setup, teardown and parsing rules. No expected values or test bodies are
translated, patched, or generated from Rust results.

```text
Original test + begin/end + original assertions
                     |
             real Tcl interpreter
                     |
          +----------+-----------+
          |                      |
   Rust DRAW adapter       native DRAWEXE
          |                      |
    Rust kernel              OCCT kernel
          |                      |
          +--- original result parser --- JSON / JUnit / logs
```

This is an incremental compatibility host, not the full DRAW application. Both
backends run in batch mode without opening a viewer. Native DRAW may link the
distribution's visualization libraries; they are test dependencies only. The
production Rust crate remains dependency-free and forbids unsafe code.

## Run

Requires Python 3.9+, Tcl 8.5+, Rust and a Unix host. Native comparison also needs
OCCT's installed DRAW executable and modeling plugins. Ordinary Rust tests still
run on Windows without Tcl or OCCT.

```sh
# Fetch and verify OCCT's public test dataset once (never done by the bridge):
python3 rust/tools/fetch_occt_test_data.py

# Rust only; explicit expected unsupported/missing cases remain non-passes:
python3 rust/tools/run_upstream_tests.py

# Same original tests on both implementations:
python3 rust/tools/run_upstream_tests.py --backend both --draw-exe /path/to/DRAWEXE

# Verify the harness catches deliberately broken assertions and crashes:
python3 -m unittest discover -s rust/tools -p test_draw_bridge.py -v

# Also cross-check signed boxes, transforms, properties and tolerance boundaries:
RUSTY_TEST_DRAW_EXE=/path/to/DRAWEXE python3 -m unittest discover -s rust/tools -p test_draw_bridge.py -v
```

The runner discovers `DRAWEXE` or Ubuntu's `occt-draw` on PATH; the macOS
fallback is `/opt/homebrew/opt/opencascade/bin/DRAWEXE`. Use
`--case tests/bugs/modalg_7/bug29311_5` to select a registered test,
`--data-dir /path/to/test-data` for further fixtures, and `--timeout 30` to set
the per-case time limit. No data files are downloaded automatically.

`locate_data_file` behaves as upstream's (`TestCommands.tcl`): the case's own
`data` folder first, then `data/` of this repository, the fetched dataset
(`target/occt-test-data/opencascade-dataset-7.9.0`, when present) and each
`--data-dir`, each with all of its subdirectories breadth-first, skipping
dot-directories. A file that cannot be found is `not_fetched` when the
dataset is absent, `private_data` when no public file has its name (Open
Cascade keeps part of its test data confidential), and `missing_fixture`
only when a public file was not found. The case variables upstream's
`_run_test` sets (`casename`, `groupname`, `gridname`, `dirname`,
`imagedir`, `test_image`) are set the same way.

**The dataset.** `fetch_occt_test_data.py` downloads
`opencascade-dataset-7.9.0.tar.xz` from the `V7_9_0_beta2` release of
Open-Cascade-SAS/OCCT (98,739,184 bytes, SHA-256
`a92ed91c3271c299287c1c404bb9454d463251094bfc848acc44aee60e6a026c`), the
file upstream's own CI downloads. It extracts the archive into the ignored
`target/occt-test-data` (3,388 files, 355 MB) and writes an inventory. It
verifies an existing archive instead of downloading it again. No file from it
is committed or redistributed: the archive has no licence file of its own,
and its basis for use is being an asset of an LGPL-2.1-with-exception release
published for upstream's tests. Manifest cases that need it are marked
`"data": true`. Their expectations are stated with the dataset present, and
without it they report `not_fetched`, which the contract accepts. CI does not
fetch it yet (see `REVIEW_NOTES.md` U1).

**Viewer commands** are recorded, not run, on both backends. The list is
`fixtures/draw-viewer-commands.txt`, taken from native DRAW's own command
groups "AIS Viewer" (without `text2brep`, which builds a shape), "DRAW Graphic
Commands" and "geometric display commands", plus `disp`, `donly`, `erase`,
`clear`, `display`, `checkview` and `checkcolor`. A case that would pass and
recorded any of them reports `viewer_skipped`: every geometric assertion was
evaluated, the image commands were not. It is not a pass. The ledger counts
its assertions as it counts a pass's, and no file is edited. `pload` of
`MODELING`, `VISUALIZATION` or `TOPTEST` is a no-op; any other module is run
natively and is a capability gap on the Rust adapter.

Reports are written to `target/upstream-tests/report.json`, `junit.xml`, and
per-case `output.log` files. They record native version/build details, source
revision, dirty Rust files, source/worker hashes, commands, observed query count,
original OCCT adjudication, unsupported commands and missing fixtures. A fresh
process and shape table isolate every test; timeout kills the whole process
group, and stale success records are removed before each run.

## Current coverage contract

| Original case | Rust | Native OCCT 7.9.3 | What it exercises |
| --- | --- | --- | --- |
| `bugs/modalg_7/bug29311_5` | Pass | Pass | Disjoint AABBs, both operand orders |
| `bugs/modalg_7/bug29311_6` | Pass | Pass | AABB containment, both orders |
| `bugs/modalg_7/bug29311_7` | Pass | Pass | AABB overlap, both orders |
| `bugs/modalg_7/bug29311_8` | Unsupported | Pass | Oriented boxes, a capability sentinel |
| `bugs/modalg_6/bug28189_2` | Unsupported | Pass | Boolean common of wire compounds |
| `bugs/modalg_6/bug28189_3` | Unsupported | Pass | Boolean union of wire compounds |
| `bugs/modalg_7/bug29333_1` | Unsupported | Pass | Face fuse, split by an edge, rebuilt with `bbuild` |
| `bugs/modalg_7/bug29333_2` | Unsupported | Pass | Face split by edges, rebuilt and queried with `modified` |
| `bugs/modalg_1/buc60684` | Unsupported | Viewer skipped | Restores a face and runs `prism` without `Copy`; needs the dataset |
| `bugs/modalg_6/bug27264_1` | Pass | Pass | Restores a box: `checknbshapes`, `checkprops -s`, `checkshape`; needs the dataset |
| `bugs/moddata_1/buc60769` | Viewer skipped | Viewer skipped | Restores a solid: `checkshape`; needs the dataset |
| Derived `prism_history_rectangle` | Pass | Pass | Prism history: `generated`, `modified`, `isdeleted` for edges, vertices and the face |
| Derived `prism_history_reversed_triangle` | Pass | Pass | Prism history against a clockwise profile's normal |
| Derived `pcylinder_counts` | Pass | Pass | Seamless cylinder through the count synthesizer: `checkshape`, `checknbshapes`, volume, area and per-use length |
| Derived `pcone_counts` | Pass | Pass | Cones through the count synthesizer: an apex cone, a frustum and a base apex, each with `checkshape`, `checknbshapes`, volume, area and per-use length (S3) |
| Derived `explode_selector` | Pass | Pass | Native selector: a box's faces and edges and a cylinder's faces and rings picked by OCCT index, checked by area, length and centre of gravity; the seam pick is lost |

There are **four original geometry tests passing on both backends** and one
more evaluated on both with its image commands recorded (`buc60769`).
The five derived cases are counted separately (see below).
The 23 bridge self-tests are separate infrastructure checks; they do not count
as more upstream coverage. The existing 66-solid / 2,292-classification native
oracle corpus supplies much broader prism geometry checks independently.

CI runs the bridge and both backends on Ubuntu 24.04 and uploads the reports.
That distribution's OCCT runtime is recorded in each report separately from
the local 7.9.3 oracle and the source fork point. Pin a source-built oracle to
the application's exact SDK before making release-parity claims.

`fixtures/upstream-draw.json` explicitly declares each expected backend status.
An unexpected failure, timeout, capability loss, new success, fixture availability
change or source-file change fails the coverage contract. **Exit zero means
the declared contract holds**, including its declared gaps. JUnit represents
unsupported/missing/known-failure cases as skipped, never as successful tests.
The JSON records which cases actually passed original assertions on both sides;
a Rust-only run makes no new live parity claim.

## Derived cases

No original self-contained test exercises prism history with commands the
adapter can run (`bugs/modalg_6/bug28261` needs `.brep` data, `revol`,
Booleans and chamfers). As `IDENTITY_AND_HISTORY.md` requires, derived cases
live in `fixtures/draw-derived/`, say on their first line that they are not
original OCCT tests, and are never counted as original upstream passes. The
manifest records each one's SHA-256 under `derived_sources` and runs it in a
pinned upstream context (`tests/bugs/modalg_7`: its begin/end scripts and parse
rules). Reports list them under `derived_cases_passing_on_both_backends`. Both
must pass on native DRAW too.

`pcylinder_counts` checks the seamless cylinder against OCCT's seamed one:
native `checknbshapes` counts 2 vertices, 3 edges and 3 wires (the seam, its
two vertices, the wall's one wire), which the Rust adapter synthesizes from a
body with two ring edges and no vertices. Its length assertion is `16π + 10`:
OCCT's `lprops` sums edges per use, so each circle counts once per face and
the seam twice.

`pcone_counts` does the same for cones (S3 of `REVIEW_NOTES.md`): OCCT's apex
cone has 2 vertices (the apex and the base circle's seam vertex), 3 edges
(the circle, the seam and a degenerated edge at the apex) and 2 wires; the
kernel's has a ring edge and the apex as a pole. Its lengths are the circles
twice and the seam twice; the degenerated edge adds nothing. Volumes and
areas come from the kernel's certified mass properties.

The history cases use `prism ... Copy`. Without `Copy`, OCCT builds the prism's end face
as the start face moved by a location, reusing its `TShape`s. DRAW's
`nbshapes` counts those shared shapes once (4 vertices, 8 edges and 5 faces
for a box), while the kernel, like `Copy`, builds distinct end entities.
Native DRAW with `Copy` reports 8, 12 and 6. The adapter rejects `prism`
without `Copy`.

The two `bug29333` cases are the self-contained history candidates for M3's
split and fuse (`IDENTITY_AND_HISTORY.md`). They split and fuse faces made by
`plane`, `mkface`, `line` and `mkedge`, which the Rust adapter cannot build, so
they are capability sentinels: native DRAW must pass them, and Rust reports
them unsupported. `bug21264` also uses splits and Booleans but needs test-local
procedures and general Booleans; it is not registered. The host forwards
their commands to either backend; the Rust worker rejects them. No derived
split or fuse case is registered: OCCT shares a split's cut face and section
edges between the two solids of its compound, while the kernel's pieces are
separate bodies, so compound counts and `generated` results differ by design.
`compare_split_merge.py` compares every split and fuse relation instead.

## Deliberate limits

- Rust signatures: positional `box` with three or six numbers, `pcylinder name
  radius height` on the default axis, two-name `copy`, single-shape
  `ttranslate`/`trotate`, `checkshape`, unique `nbshapes` with synthesized
  seams and seam vertices (`Topology::occt_counts` for whole solids), `vprops`
  with optional positive integration epsilon, `isdraw`, the solid identification
  needed by `checkprops`, and AABB `isbbinterf`. For history: a closed
  `polyline`, `mkplane` of it, `prism name face dx dy dz Copy` normal to the
  profile, `explode` of a profile wire or face into edges or vertices,
  `savehistory`, `generated`, `modified`, `isdeleted`, and `sprops`/`lprops`
  (mass only) on faces, edges, solids and compounds, summing lengths per edge
  use as OCCT's explorer does. `explode` of a kernel solid into faces or
  edges goes through the native selector (below); `sprops` of a selected
  face and `lprops` of a selected edge also report its centre of gravity. `generated` follows OCCT's prism:
  a profile vertex gives its vertical edge, an edge its wall, the face the
  solid; the kernel's start/end copies are OCCT's `FirstShape`/`LastShape`. Everything else is unsupported,
  including OBBs, `sprops`, `nbshapes -t`, option-form boxes and visualization.
- Test geometry remains boxes and rigidly transformed copies. The adapter
  intentionally does not expose every Rust prism capability yet. More command
  adapters must delegate to real kernel functions, never return canned success.
- Metadata `help` registers no UI help; `cpulimit` is handled by the outer wall
  timeout. An existing topology-check registration avoids `bugs/begin` loading
  a viewer. These administrative differences do not skip geometry assertions.
- A safe child interpreter supplies Tcl control flow and prevents tests from
  launching processes or writing files. This is isolation for trusted pinned
  tests, not a security sandbox for arbitrary untrusted native DRAW input.
- Unsupported commands and missing fixtures are latched outside the child, so
  `catch` cannot turn them into passes. No-observation scripts are unverified.
  Upstream BAD/TODO and IMPROVEMENT results are distinct from passes. Original
  assertions retain their own tolerances, sometimes much looser than the
  separate numerical oracle; passing them is not a universal error guarantee.
- OCCT's C++ GoogleTest files are not run against Rust yet. Their APIs cannot be
  linked directly to a Rust crate. A future narrow, test-only C ABI/C++ adapter
  could compile selected original tests unchanged; otherwise reuse their inputs
  and assertions with explicit attribution and label them translated tests.

## Extend coverage

Read the original command implementation and the complete test context. Add
only modeling/query commands needed by noBS-CAD, plus real Rust functionality.
Choose an original test with public fixtures, record the SHA-256 of it and all
of its setup/parser dependencies in the manifest, and require native success
before promoting Rust coverage. Register missing capability/data cases explicitly.
Do not strip viewer calls from a file and then count it as an unchanged passing
test: classify it outside the headless subset or create a separately attributed
derived case. Do not use a success percentage across all OCCT suites as the
readiness metric; measure evidence against the documented kernel capabilities,
numerical domains and failure contracts independently of any application.

## Structure mapping and the coverage ledger

The kernel's topology model (`TOPOLOGY_MODEL.md`) has no seams, degenerate
edges or per-entity tolerances, so original assertions that depend on OCCT's
structure cannot be evaluated on Rust output directly. None of the mappings
rewrites an original file.

* **Count synthesizer** (T1): the adapter's `nbshapes` reports the counts OCCT
  would give for the same body: a seam per periodic direction of a wound
  face, a seam vertex per ring edge, loops as wires, a wound face's loops as
  one wire. Lengths sum per edge use, the seam twice (`lprops`).
* **Native selector** (T2): for a test registered with `"selector": true`,
  the runner first runs it in native DRAW with `RUSTY_DRAW_SELECTOR_OUT` set.
  The host then records every pick of every `explode`: its kind, measure
  (`sprops`/`lprops`, which DRAW prints to six significant digits) and centre
  of gravity (Draw variables, full precision). In the Rust run the adapter
  numbers `explode` calls alike and selects, for each pick, the one kernel
  entity of that kind whose measure agrees within `1e-5` relative and whose
  centre agrees within `1e-7` of the centre's scale. It never reproduces
  OCCT's exploration order. A pick with no entity (a seam) or with several is
  lost: using it is a capability gap, never a pass. A Rust-only run cannot
  select and reports such cases as skipped. Only faces and edges are
  selected; vertex, wire and shell picks are unsupported until a case
  confirms them.
* **`.brep` interop** (T2, `occt_brep`): the converter and writer, verified by
  native round trips in `compare_brep_io.py` (see `VALIDATION.md`).

### The ledger

`run_upstream_tests.py` recomputes the ledger on every run and writes every
assertion with its status to `target/upstream-tests/ledger.json`; the
aggregate goes into the report and must equal the `ledger` block in
`fixtures/upstream-draw.json`. A change fails the contract until it is
reviewed and recorded with `--write-ledger` (`--ledger` recomputes and checks
only). The scan covers every case file under `tests/` (the survey's
17,879); an assertion is a `check*` procedure at a command position or a
custom `puts "Error..."` report.

* `model-independent`: the value does not depend on how OCCT structures the
  body. Validity verdicts, volumes, areas, centres, points, real values,
  images and meshes: `checkshape`, `checkprops` without `-l`, `checkreal`,
  `checkview`, `checktrinfo`, `checkgravitycenter` and the tests' own
  helpers and error reports.
* `mapped-and-verified`: structure-dependent, carried by a mapping (count
  synthesizer for `checknbshapes`, per-use length for `checkprops -l` and
  `checklength`, the native selector for any assertion on a pick), and the
  case is registered with its original assertions passing on both backends.
* `lost`: every other structure-dependent assertion, with its reason: a
  mapping exists but native DRAW has not confirmed it on that case, or none
  exists (`checkfreebounds`, `checksection`, `checkmaxtol` (OCCT's grown
  tolerances are not the kernel's certified enclosures),
  `checkfaults`, `checkloc`, `checkoverlapedges`, `checkcurveonsurf`,
  `check_fsd`). An assertion that names an exploded sub-shape depends on the
  selector whatever its procedure.

Recomputed survey at the pinned revision (17,879 cases): 11,581 load
external data, 7,428 touch the viewer; property assertions in 8,469 files,
sub-shape counts in 4,735, validity in 2,726, tolerance maxima in 961. The
2026-09-25 estimate (11,599; 7,658; 8,415; 4,689; 2,655; 860) used other
search patterns; the ledger's are the definition from now on.

| Status | Assertions |
| --- | ---: |
| model-independent | 22,920 |
| mapped-and-verified | 1 |
| lost | 12,845 |

| Lost because | Kind | Assertions |
| --- | --- | ---: |
| count synthesizer | mapping unverified | 5,399 |
| per-use length synthesis | mapping unverified | 1,975 |
| native selector | mapping unverified | 1,578 |
| section wire and vertex counts follow OCCT's splitting | no mapping | 1,348 |
| OCCT tolerances are grown requests; enclosures are certified bounds | no mapping | 978 |
| free boundaries count seams and degenerate edges | no mapping | 765 |
| per-subshape fault statuses | no mapping | 700 |
| OCCT persistence formats | no mapping | 75 |
| per-edge overlap follows OCCT's edges | no mapping | 17 |
| locations are OCCT structure | no mapping | 5 |
| per-edge deviations follow OCCT's edges | no mapping | 5 |

The one mapped-and-verified assertion is `bug27264_1`'s `checknbshapes`: a
restored box through the count synthesizer, confirmed by native DRAW on the
same case. The other original cases passing on both backends make only
model-independent assertions. Each mapping is also confirmed natively on a
derived case (`mappings_confirmed_natively` in the report: the count
synthesizer and per-use length on `pcylinder_counts`, the selector on
`explode_selector`), but derived cases never count toward the ledger.

**Data-dependent cases.** `restore` goes through the T2 reader and converter
(`occt_brep`). A file's compounds stay compounds and each solid becomes a
body with its OCCT tolerances as `Imported` enclosures. `checkshape`,
`nbshapes`, `sprops`, `lprops` and `explode` (native selector) work on it.
`vprops` waits for general mass properties (S3). A construct the kernel cannot
represent makes the restore `unsupported` by name, and so does a file of
another DRAW type (a saved curve or surface). A solid the converter builds
but the validator rejects is a failure.

`survey_upstream_tests.py --restore-only` runs the 192 cases that only
restore data and run checks on both backends, without registering them. With
the dataset (2026-09-26):

| Rust / native | Cases |
| --- | ---: |
| private data / private data | 98 |
| unsupported / unsupported (the case needs other commands) | 31 |
| unsupported / viewer skipped | 26 |
| unsupported / known failure | 11 |
| unsupported / unverified | 9 |
| unsupported / private data | 5 |
| unsupported / pass | 5 |
| unsupported / failed | 4 |
| known failure / known failure | 1 |
| pass / pass | 1 |
| viewer skipped / viewer skipped | 1 |

Native DRAW evaluates 33 of them completely; Rust evaluates 2, both now
registered. What the other 31 need, by construct (a case may need several):
free faces 22, B-spline curves on surfaces 18, B-spline curves 16, B-spline
surfaces 13, trimmed surfaces 9, Bézier surfaces 3, tori 3, cones 2,
spheres 2, extrusion surfaces 2, and one each of the rest. Two need Booleans.
After the cone (S3), six of these cases no longer report cones. None
evaluates yet: of the two native DRAW evaluates, `bug485` now needs only
tori and `bug437` B-spline curves and a free face; `bug497_2` restores
completely and then needs `bcut`; the other three also need B-spline or
revolution surfaces.
