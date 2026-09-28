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
| `bugs/modalg_6/bug27264_2` | Pass | Pass | Restores a whole sphere: `checknbshapes`, `checkprops -s`, `checkshape`; needs the dataset (S3) |
| `bugs/moddata_1/buc60769` | Viewer skipped | Viewer skipped | Restores a solid: `checkshape`; needs the dataset |
| Derived `prism_history_rectangle` | Pass | Pass | Prism history: `generated`, `modified`, `isdeleted` for edges, vertices and the face |
| Derived `prism_history_reversed_triangle` | Pass | Pass | Prism history against a clockwise profile's normal |
| Derived `pcylinder_counts` | Pass | Pass | Seamless cylinder through the count synthesizer: `checkshape`, `checknbshapes`, volume, area and per-use length |
| Derived `psphere_counts` | Pass | Pass | Spheres through the count synthesizer: a whole sphere (no loops; OCCT's seam, two pole vertices and degenerated edges), a hemisphere and a zone, with `checkshape`, `checknbshapes`, volume, area and per-use length (S3) |
| Derived `ptorus_counts` | Pass | Pass | Tori through the count synthesizer: a whole torus (no loops; OCCT's two seams through one vertex), a v-segment and a wedge, with `checkshape`, `checknbshapes`, volume, area and per-use length (S3) |
| Derived `pcone_counts` | Pass | Pass | Cones through the count synthesizer: an apex cone, a frustum and a base apex, each with `checkshape`, `checknbshapes`, volume, area and per-use length (S3) |
| Derived `explode_selector` | Pass | Pass | Native selector: a box's faces and edges and a cylinder's faces and rings picked by OCCT index, checked by area, length and centre of gravity; the seam pick is lost |
| Derived `profile_arcs` | Pass | Pass | Prisms of `profile` sketches with arcs: tangent half circles, fillets and a clockwise notch; counts and volumes (S5) |
| Derived `faces_and_edges` | Pass | Pass | Free faces and edges on DRAW geometry: `mkface` on a plane and a cylinder (a patch and the whole band), `mkedge` on a line and a circle, `mkplane` of a `profile` wire with arcs; `checkshape`, `checknbshapes`, areas and lengths (S6) |
| Derived `split_plane` | Pass | Pass | OCCT's splitter (`bclearobjects`, `bcleartools`, `baddobjects`, `baddtools`, `bfillds`, `bsplit`, `bapisplit`) with one plane face as the tool: a prism parallel to its axis, a box cut obliquely, a stadium sheet across an arc, a closed wire, and a prism, a face and a wire at once; `checkshape`, `checknbshapes` with the shared cut entities, volumes, areas and lengths per use (S8e) |
| `lowalgos/intss/bug23177_1` | Viewer skipped | Viewer skipped | A torus and a plane through its axis: two circles, `dump`, `bounds`, `dval`, `xdistcs` on both surfaces (S7) |
| `lowalgos/intss/bug23177_2` | Viewer skipped | Viewer skipped | The same with another plane through the axis (S7) |
| `lowalgos/intss/bug24648` | Viewer skipped | Viewer skipped | A cone and a coaxial torus: every curve a circle (S7) |
| `lowalgos/intss/bug21750` | Unsupported | Unverified | Cylinders antiparallel up to rounding: OCCT snaps them to two lines, the kernel's exact predicates find crossing axes and one closed curve, and the case then calls the undefined `Error:` (S7, a recorded divergence) |
| `bugs/heal/bug29502` | Unsupported | Pass | A whole cylinder band split by a vertex (`bsplit`), then `unifysamedom`; the `heal` group loads `XSDRAW` (S8e) |
| `bugs/heal/bug33171_1` | Unsupported | Viewer skipped | A prism split by four prisms of open polylines, `bopcheck`, `unifysamedom` (S8e) |
| `bugs/modalg_7/bug21264` | Unsupported | Unverified | Progress reports of the Boolean and splitter commands (`XProgress`); no geometric query (S8e) |
| `bugs/modalg_7/bug30092` | Unsupported | Pass | A face on an offset surface split by an edge with a grown tolerance; the restore needs offset surfaces; needs the dataset (S8e) |
| `bugs/modalg_7/bug32578` | Unsupported | Viewer skipped | A face split by many edges with a fuzzy value (`bfuzzyvalue`); needs the dataset (S8e) |
| `bugs/moddata_3/bug31587_1` to `_6` | Unsupported | Pass | A box split by another box's faces, edges, wires, vertices or open polylines, then `removeinternals` (S8e) |
| `boolean/splitter/A5`, `B5`; `bugs/modalg_7/bug28113_1`, `bug28113_2`, `bug29789`, `bug29955`, `bug31201_1` to `_3`, `bug31462`, `bug32644` | Private data | Private data | Splits of restored shapes (`bsplit`, and `bapisplit` in `B5`) whose files are not in the public dataset (S8e) |

There are **five original geometry tests passing on both backends** and four
more evaluated on both with their image commands recorded (`buc60769`, and
S7's `lowalgos/intss` cases `bug23177_1`, `bug23177_2` and `bug24648`).
S8e registers the upstream `bsplit` group as capability sentinels: `boolean/splitter/A5`
and `B5` and the 22 `bugs` cases that call `bsplit` (with `bug29333_1` and
`bug29333_2`, registered before). They need S9's general builder: tools that
are edges, vertices, several faces or shapes of their own, and restored
curved shapes. Of the thirteen that run on public data, native DRAW passes
ten (the two `bug29333` cases among them), evaluates two more with their
image commands recorded (`bug33171_1`, `bug32578`) and makes no geometric
query in `bug21264`; the other eleven need data Open Cascade keeps private.
The Rust adapter reports each of the thirteen unsupported. The host forwards the group's commands (`bbuild`, `bapibuild`,
`bop` and its operations, `bbop`, `bapibop`, `bcut`, `btuc`, `bsection`,
`bopcheck`, `boptions`, `bfuzzyvalue`, `removeinternals`, `unifysamedom`,
`vertex`, `settolerance`, `XProgress`, `dset`, `protect`) to either backend;
the Rust worker rejects all but `dset` (DRAW's numeric variables) and
`protect` (a no-op), which the `boolean` group's `begin` calls.
Three more `intss` cases run on the Rust adapter but are not registered,
because the contract admits no failing status: `bug23178`, `bug28222_2` and
`bug28222_3` count the pieces IntPatch splits its walking lines into (6, 4
and 2), where the kernel returns closed loops (1, 2 and 1) whose `xdistcs`
samples lie within 1.5e-14 of both cylinders.
The ten derived cases are counted separately (see below).
The 24 bridge self-tests are separate infrastructure checks; they do not count
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

`profile_arcs` (S5) builds a stadium, a rectangle with four fillets and a
rectangle with a concave notch with upstream's `profile` sketch command (tangent
`C` arcs, a negative radius for a clockwise arc) and extrudes them: counts
(a polygon's, no seams) and volumes on both backends. The adapter follows
`BRepTest_CurveCommands.cxx`'s `profile` exactly for `F`, `O`, `P`, `X`,
`Y`, `L`, `T`, `R`, `D`, `I`, `C` and `W`; `S` (another face's surface) and
open wires (`WW`) are unsupported.

`faces_and_edges` (S6) makes free faces and edges on DRAW's geometric
objects: `plane`, `cylinder`, `line` and `circle` (3D, with DRAW's automatic
X direction), `mkface` of a plane rectangle (a face body of the kernel's
`Body::face_from_profile`), of a cylinder patch and of the whole band (OCCT
closes it with a seam and a vertex on each circle; the kernel's band has two
ring edges and wound loops), `mkedge` of a line segment, an arc and a whole
circle (one vertex), and `mkplane` of a `profile` wire with arcs. It checks
counts (a free face has no shell, a free edge no wire), areas and lengths on
both backends. `mkface` with a wire, unbounded surfaces and 2D curves are
unsupported.

`split_plane` (S8e) drives OCCT's General Fuse splitter commands with one
plane face (`plane` then `mkface`) as the tool; the adapter splits each
object by the face's plane through `Solid::split_by_plane` and
`Body::split_by_plane`. Its values are computed by hand and hold on both
backends:

* a 6 x 4 x 3 prism (`profile`, `prism`) by `x + y = 3`: 12 vertices, 20
  edges, 11 wires and faces, 2 shells and solids, 1 compound; volume 72,
  area `108 + 18 sqrt 2`, length `128 + 24 sqrt 2`;
* a 4 x 4 x 2 `box` cut obliquely (S8a.2) at a corner tetrahedron with legs
  1, 2 and 10/7: 11 vertices, 18 edges, 10 wires and faces; volume 32, area
  `64 + sqrt 696 / 7`, length `80 + 4 sqrt 5 + 4 (sqrt 149 + sqrt 296) / 7`;
* a stadium sheet (area `6 + pi`) by `x = 3.5` across its right half
  circle, which splits in three: 6 vertices, 7 edges, 2 wires and faces;
  length `6 + 2 pi + 2 sqrt 3`;
* a closed 5 x 2 rectangle wire by `x - y = 1`: 6 vertices and edges, one
  wire and no compound (OCCT's splitter keeps a wire one wire of its split
  edges, and returns a single image itself); length 14;
* a prism, a plane face and a stadium wire apart from each other, by
  `x = 1.5` through `bapisplit`: 24 vertices, 33 edges, 14 wires, 13 faces;
  area 80, length `142 + 2 pi`.

OCCT shares each cut face, chord and crossing vertex between the pieces it
separates, while the kernel's pieces are separate bodies: for a split's
result the adapter counts an entity of one object's pieces that coincides
with an earlier piece's (a vertex at its point, an edge between the same
vertices through its midpoint, a face bounded by the same edges) once, and
measures a shared face in each solid, as OCCT's explorer visits it; a split
wire's open runs count as one wire. The adapter reports unsupported, rather
than splitting: several tools or a tool other than a convex polygonal plane
face; a tool face not covering an object's bounding box's shadow on its
plane (OCCT cuts only where the face reaches); objects whose bounding boxes
overlap (General Fuse intersects them too); a plane through a vertex or
tangent to an arc (where OCCT's sharing is not confirmed); ring edges (OCCT's
seams); solids other than prisms, restored shapes and split results; and
the kernel's `OutOfDomain`, `Degenerate`, `ComputationLimit`,
`LimitExceeded` and `PrecisionLoss` refusals. `vprops` of a split sums its
solids (a face's signed volume is not supported), and `savehistory` after a
split is unsupported.

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
them unsupported (`bug29333_1` moves a face with `ttranslate`, `bug29333_2`
needs `mkvolume`). `bug21264` also uses splits and Booleans; since S8e it is
registered with the `bsplit` group. The host forwards their commands to
either backend. No derived fuse case is registered, and `split_plane`
checks counts and properties only: the adapter keeps no split history for
`savehistory`, since the kernel's history names each piece's own cut face,
not OCCT's shared one. `compare_split_merge.py` compares every split and
fuse relation instead.

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
  For splits (S8e): `bclearobjects`, `bcleartools`, `baddobjects`,
  `baddtools`, `bfillds` and `bsplit name` (or `bapisplit name`) with one
  plane face as the tool, on prisms, planar faces and closed planar wires
  (see `split_plane` above); `bsplit` before `bfillds` prints OCCT's
  message, and after the arguments changed it is unsupported.
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
| mapped-and-verified | 2 |
| lost | 12,844 |

| Lost because | Kind | Assertions |
| --- | --- | ---: |
| count synthesizer | mapping unverified | 5,398 |
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

The two mapped-and-verified assertions are `bug27264_1`'s and
`bug27264_2`'s `checknbshapes`: a restored box and, since the sphere (S3), a
restored whole sphere, both through the count synthesizer and confirmed by
native DRAW on the same case. The other original cases passing on both backends make only
model-independent assertions. Each mapping is also confirmed natively on a
derived case (`mappings_confirmed_natively` in the report: the count
synthesizer and per-use length on `pcylinder_counts`, the selector on
`explode_selector`), but derived cases never count toward the ledger.

**S2 of `REVIEW_NOTES.md` — accepted at `1a76d29e`.** S2 (the fetch
script, `restore` through the `.brep` converter, the viewer-skipped and
not-fetched statuses, the recorded ledger) is `21f34dad`, whose kernel
workflow passed; its fuzzing run was cancelled by the next push. Both
workflows passed at `1a76d29e`, which adds only the cone's validation
rules (S3): twelve kernel jobs, including this bridge on the Rust adapter
and native DRAW with the recorded ledger, and twenty-four fuzz targets. No
workflow downloads the dataset; data cases report `not_fetched` on CI.

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
the dataset, after S4 (2026-09-27):

| Rust / native | Cases |
| --- | ---: |
| private data / private data | 98 |
| unsupported / unsupported (the case needs other commands) | 31 |
| unsupported / viewer skipped | 26 |
| unsupported / known failure | 11 |
| unsupported / unverified | 9 |
| unsupported / private data | 5 |
| unsupported / failed | 4 |
| unsupported / pass | 4 |
| pass / pass | 2 |
| known failure / known failure | 1 |
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
revolution surfaces. After the sphere, `bug27264_2` (a restored whole
sphere) passes on both backends and is registered; spheres no longer
appear among the unsupported constructs. After the torus, `bug485` restores
completely and now needs only `bfuse`; `bug327_3` and `bug27782` no longer
report tori; `bug503`'s one torus is a horn torus (`NonRingToroidalSurface`,
not representable). No further case evaluates on the Rust backend.
After S4 no case reports B-spline curves, 2D curves or surfaces. What the
restores still need, by construct: free faces 46, rectangular trimmed
surfaces 12, ellipses 5, trimmed curves (of other bases) 4, extrusion and
revolution surfaces 4 each, 2D ellipses 4, Bézier surfaces 3, and one or
two each of the rest. `bug21246` has internal edges and faces: the importer
panicked on them (an orientation composed as internal was taken for
forward or reversed); it now reports them unsupported by name
(`InternalOrExternalEdge`, `InternalOrExternalFace`), as the reference
reader does, and `occt_brep.rs` and the `brep_io` fuzz target make
references internal or external. No further case evaluates on the Rust
backend.

After S6 `restore` imports free shells, faces, wires, edges and vertices;
no restore-only case reports free faces any more. The survey (2026-09-27):

| Rust / native | Cases |
| --- | ---: |
| private data / private data | 98 |
| unsupported / unsupported | 31 |
| unsupported / viewer skipped | 20 |
| unsupported / known failure | 10 |
| failed / viewer skipped | 6 |
| unsupported / private data | 5 |
| unsupported / unverified | 4 |
| unverified / unverified | 3 |
| unsupported / failed | 3 |
| failed / unverified | 2 |
| failed / pass | 2 |
| pass / pass | 2 |
| unsupported / pass | 2 |
| failed / known failure, failed / failed, known failure / known failure, viewer skipped / viewer skipped | 1 each |

Twelve cases now restore their free faces and fail on Rust: the validator
rejects a face OCCT accepts, and a restore the validator rejects is a
failure. Each rejection was inspected: uv gaps and pcurves off their edges
beyond the stored tolerance (`OCC221`, `OCC302a`, `OCC399`, `OCC446d`,
`OCC889`, `OCC25558_faulty`, `bug24035`, `bug28499`; up to 0.049 in UV), a
zero-length edge (`OCC432`), a spline face whose only loop runs clockwise
(`OCC161`), a sphere lune whose loop passes both poles between two uses of
one seam (`OCC35`, the pinned `Ball.brep` limitation) and bug 28385's
deliberately degenerate face (two parallel lines one apart closed by
vertices of tolerance 1). None is an import error; they stay failures. What the restores still need, by
construct: rectangular trimmed surfaces 13, degenerated edges outside a
pole 5, ellipses 5, extrusion and revolution surfaces, trimmed curves and
periodic spline surfaces 4 each.
