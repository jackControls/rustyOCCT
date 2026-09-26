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
`--data-dir /path/to/test-data` for existing fixtures, and `--timeout 30` to set
the per-case time limit. No data files are downloaded automatically.

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
| `bugs/modalg_1/buc60684` | Missing fixture | Missing fixture | External `buc60684a.brep` data |
| Derived `prism_history_rectangle` | Pass | Pass | Prism history: `generated`, `modified`, `isdeleted` for edges, vertices and the face |
| Derived `prism_history_reversed_triangle` | Pass | Pass | Prism history against a clockwise profile's normal |
| Derived `pcylinder_counts` | Pass | Pass | Seamless cylinder through the count synthesizer: `checkshape`, `checknbshapes`, volume, area and per-use length |
| Derived `explode_selector` | Pass | Pass | Native selector: a box's faces and edges and a cylinder's faces and rings picked by OCCT index, checked by area, length and centre of gravity; the seam pick is lost |

There are **three original geometry tests passing on both backends**, not nine.
The four derived cases are counted separately (see below).
The 19 bridge self-tests are separate infrastructure checks; they do not count
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
  exists (`checkfreebounds`, `checksection`, `checkmaxtol` until M5,
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
| mapped-and-verified | 0 |
| lost | 12,846 |

| Lost because | Kind | Assertions |
| --- | --- | ---: |
| count synthesizer | mapping unverified | 5,400 |
| per-use length synthesis | mapping unverified | 1,975 |
| native selector | mapping unverified | 1,578 |
| section wire and vertex counts follow OCCT's splitting | no mapping | 1,348 |
| per-entity tolerances (enclosures are M5) | no mapping | 978 |
| free boundaries count seams and degenerate edges | no mapping | 765 |
| per-subshape fault statuses | no mapping | 700 |
| OCCT persistence formats | no mapping | 75 |
| per-edge overlap follows OCCT's edges | no mapping | 17 |
| locations are OCCT structure | no mapping | 5 |
| per-edge deviations follow OCCT's edges | no mapping | 5 |

No original assertion is mapped-and-verified yet. The three original cases
that pass on both backends only make model-independent assertions, and no
other original case with a structure-dependent assertion runs on the adapter.
Each mapping is confirmed natively on a derived case
(`mappings_confirmed_natively` in the report: the count synthesizer and
per-use length on `pcylinder_counts`, the selector on `explode_selector`), but
derived cases never count toward the ledger.

**Data-dependent cases.** The only public data in the repository is
`data/occ` (37 files). Three original cases name its files: one heal-grid
data file whose grid runs shape healing, and two viewer cases. Every other
data-dependent case needs OCCT's external test-data set, which is not in the
repository and is never downloaded. The adapter therefore has no `restore`
yet: it arrives with the first original data-dependent case a configured
`--data-dir` makes runnable on both backends, so that native DRAW confirms it
in the same change.
