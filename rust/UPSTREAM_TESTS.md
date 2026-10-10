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
| Derived `boolean_prisms` | Pass | Pass | Booleans of prisms in one frame (`bfuse`, `bcut`, `bcommon`, `bop` with `bopfuse`, `bopcut`, `boptuc`, `bopcommon`, `bbop`, `bapibop`): overlapping boxes, a box and a cylinder (a half cylinder, a hole), a box cut in two, disjoint boxes, one box inside another, stacked boxes of one profile, a stadium prism and a bar; `checkshape`, volumes, areas, and counts and lengths per use after `unifysamedom` (S9a) |
| Derived `boolean_stacks` | Pass | Pass | Booleans of prisms in one frame whose slabs hold different regions, the kernel's stacks (`bfuse`, `bcut`, `bop` with `bopfuse` and `bopcut`, `bbop`): a step, a pocket, a box cut into two solids by a slab, a closed cavity (one solid of two shells) and a tool through a round wall; `checkshape`, volumes, areas, and counts and lengths per use after `unifysamedom` (S9a.2) |
| Derived `boolean_polyhedra` | Pass | Pass | Booleans of prisms in frames with different axes, S9b.1's polyhedra (`bfuse`, `bcut`, `bcommon`, `bop` with `bopfuse` and `boptuc`, `bbop`; tools moved by `tcopy` and `trotate`): a box and its copy turned a quarter turn (an L whose coplanar walls are one face, a plus sign's fuse, common and cuts), a bar turned 45 degrees through a box (an octagon, four corners, a star), a tilted bar cutting a box in two, and a turned box inside another (a closed cavity); `checkshape`, volumes, areas, and counts and lengths per use after `unifysamedom` (S9b.1) |
| Derived `split_plane` | Pass | Pass | OCCT's splitter (`bclearobjects`, `bcleartools`, `baddobjects`, `baddtools`, `bfillds`, `bsplit`, `bapisplit`) with one plane face as the tool: a prism parallel to its axis, a box cut obliquely, a stadium sheet across an arc, a closed wire, and a prism, a face and a wire at once; `checkshape`, `checknbshapes` with the shared cut entities, volumes, areas and lengths per use (S8e) |
| `lowalgos/intss/bug23177_1` | Viewer skipped | Viewer skipped | A torus and a plane through its axis: two circles, `dump`, `bounds`, `dval`, `xdistcs` on both surfaces (S7) |
| `lowalgos/intss/bug23177_2` | Viewer skipped | Viewer skipped | The same with another plane through the axis (S7) |
| `lowalgos/intss/bug24648` | Viewer skipped | Viewer skipped | A cone and a coaxial torus: every curve a circle (S7) |
| `lowalgos/intss/bug21750` | Unsupported | Unverified | Cylinders antiparallel up to rounding: OCCT snaps them to two lines, the kernel's exact predicates find crossing axes and one closed curve, and the case then calls the undefined `Error:` (S7, a recorded divergence) |
| `bugs/heal/bug29502` | Unsupported | Pass | A whole cylinder band split by a vertex (`bsplit`), then `unifysamedom`; the `heal` group loads `XSDRAW` (S8e) |
| `bugs/heal/bug33171_1` | Unsupported | Viewer skipped (OCCT 7.6.3, Linux CI's package: a known failure, its unified shape invalid) | A prism split by four prisms of open polylines, `bopcheck`, `unifysamedom` (S8e) |
| `bugs/modalg_7/bug21264` | Unsupported | Unverified | Progress reports of the Boolean and splitter commands (`XProgress`); no geometric query (S8e) |
| `bugs/modalg_7/bug30092` | Unsupported | Pass | A face on an offset surface split by an edge with a grown tolerance; the restore needs offset surfaces; needs the dataset (S8e) |
| `bugs/modalg_7/bug32578` | Unsupported | Viewer skipped | A face split by many edges with a fuzzy value (`bfuzzyvalue`); needs the dataset (S8e) |
| `bugs/moddata_3/bug31587_1` to `_6` | Unsupported | Pass | A box split by another box's faces, edges, wires, vertices or open polylines, then `removeinternals` (S8e) |
| `boolean/splitter/A5`, `B5`; `bugs/modalg_7/bug28113_1`, `bug28113_2`, `bug29789`, `bug29955`, `bug31201_1` to `_3`, `bug31462`, `bug32644` | Private data | Private data | Splits of restored shapes (`bsplit`, and `bapisplit` in `B5`) whose files are not in the public dataset (S8e) |
| `boolean/bopcommon_simple` (122 cases), `bopcut_simple` (106), `bopfuse_simple` (89), `boptuc_simple` (83), `bcommon_simple/I5`, `J1`, `bfuse_simple/L2`, `E2`, `bcut_simple/G7`, `L8` | Viewer skipped | Viewer skipped | `bop` and its operations, `bfuse`, `bcut` and `bcommon` of two boxes, a box and a `pcylinder`, or two cylinders in one frame (some moved by `ttranslate`, some sized by `dset`): `checkprops -s` (or `-s empty`) and the group's `checkshape`; each records a `checkview` (S9a; 86 of them, `bopfuse_simple/B2`, `Z8` and `ZB3` among them, are S9a.2's stacks) |
| `boolean/bopcommon_simple` (81 more cases), `bopcut_simple` (68), `bopfuse_simple` (69), `boptuc_simple` (80) | Viewer skipped | Viewer skipped | `bop` and its operations on two boxes in frames with different axes (one turned by `trotate` about z, or about x and moved by `ttranslate`, some sized by `dset`; `K3` and `P6` by quarter turns, which the adapter turns exactly): `checkprops -s` (or `-s empty`) and the group's `checkshape`; each records a `checkview` (S9b.1's polyhedra) |
| `boolean/bopcommon_simple/ZL6` | Unsupported | Viewer skipped | A frustum of radii 8 and 4 and height 8 and one of radii 4 and 2 standing on its top disc, on the first's cone, their virtual apexes one point: a cone's apex on the other input's surface (`Degenerate` since S9d.4b's survey, when the adapter built two `pcone`s with ids of their own; solids other than prisms sharing ids before; an S9a sentinel) |
| `boolean/bopcommon_simple/ZP9` | Viewer skipped | Viewer skipped | `bop` and `bopcommon` of a 100 box and a sphere of radius 7.5 placed on a DRAW `plane` (`psphere name plane R`, its centre 1.053 inside the box's wall `x = 100`), the wall cutting a cap off the sphere: `checkprops -s` and the group's `checkshape`; it records a `checkview` (S9d.1's spheres against polyhedral prisms; the adapter places a sphere on a plane since S9d.1's survey; the wall far from tangency, which `b2765f20` refuses within the resolution) |
| `boolean/bcommon_simple/A1`, `bfuse_simple/A4` | Unsupported | Viewer skipped | A unit sphere and a unit box whose corner is at the centre, its far corners on the sphere (an S9a sentinel), the box quarter-turned so a wall is tangent to the sphere: a vertex of one input on the other's face (S9d.2a's `Degenerate`; S9d.1 reported it as a meeting at every seam tried), a tangency between the inputs (S9d.1's `Degenerate`) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZI8`, `ZI9`, `ZJ1`, `ZJ2`, `ZJ3` (20 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` of radius 4 and height 8 and a sphere of radius 6 centred on its top cap (turned by `trotate` about the cap's centre: whole quarter turns about z in `ZI9` to `ZJ2`, about y in `ZJ3`): the cap inside the sphere, the wall meeting it in a parallel at height `8 - 2 sqrt(5)`: `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9d.2a's coaxial rings) |
| `boolean/bopfuse_simple/ZH5` | Unsupported | Viewer skipped | A cylinder of radius 4 and a sphere of radius 4 centred on its top cap: the rim is the sphere's equator, the wall tangent to the sphere along it; the coaxial pair's discriminant vanishes identically: a tangency between the inputs (`Degenerate`; the survey below found it refused as a computation limit, S9d.1 as S9d.2's `OutOfDomain`) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZI4` to `ZI7` (16 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` of radius 4 and height 8 and a sphere of radius 2 centred on its top cap, turned a quarter turn about x (and then one, two or three quarter turns about y in `ZI5`, `ZI6` and `ZI7`), so the cap's plane holds the sphere's axis: the section a great circle through its poles, a vertex at each: `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9d.1's pole follow-up after S9d.3a; `bopfuse_simple/ZI4` was the sentinel for its `PrecisionLoss`; `ZI5`, whose volumes were wrong until a sphere's closing chord at a pole was enclosed narrowly, registered after it: see the S9d.3a survey; their values unchanged bit for bit by S9e.3b's pole vertices, the survey of S9e.3b, S9f.2a and S9f.2b.1) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZF5` to `ZF9`, `ZH1` to `ZH4` (36 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on a box of side 4 and a `pcone` frustum: of radii 1 and 0.5 on the vertical through the box's centre, standing on its top face, inside it (its top disc on the top face or clear of every face) or from the bottom face or below it to the top face, meeting the faces in circles (`ZF5` to `ZF9`); of radii 5 and 4 (5 and 3.5 in `ZH4`) and the box's height, its axis 2 outside a wall and 2 from the walls across it, meeting them in hyperbolas (`ZH1` to `ZH4`; in `ZH3` and `ZH4` turned 30 degrees about its axis): `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9d.3a's cones against polyhedral prisms) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZJ4`, `ZJ6` to `ZJ9`, `ZK5` to `ZK9`, `ZL1` (44 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` of radius 4 and height 8 and a frustum of radii 2 and 1: coaxial, standing on its top cap, inside it with its discs on both caps or clear of them, or through one or both caps (`ZJ4`, `ZJ6` to `ZJ9`: circles), or of height 10 with its axis crossing the cylinder's at right angles at half its height, turned a quarter turn about y or x (`ZK5`, `ZK6`: two quartic rings over the frustum's angle; `ZK7`, `ZK8`: of radii 6 and 1, its wide end through both caps too), or a frustum of radii 1 and 6 and height 8 turned a quarter turn about y either way, its wide end through both caps and its end disc across the wall (`ZK9`, `ZL1`): `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9d.3b.1's cones against cylinders; `bopfuse_simple/ZJ4` was the sentinel for S9d.3b's `OutOfDomain`, `bopfuse_simple/ZK9` for a cone turned by the rounded rotation until the adapter turned it exactly; `ZK7` and `ZK8`, right since S9d.3b.1's survey but past the contract's 30 seconds, registered in the survey of S9d.2c, S9d.3c and S9d.4c) |
| `boolean/bopfuse_simple/ZK1` | Unsupported | Viewer skipped | A `pcylinder` of radius 4 and height 8 and a frustum of radii 4 and 2 standing on its top cap, its base rim the cap's rim and its base disc the cap: a tangency between the inputs (`ZK2` to `ZK4`, the same turned about its axis, alike since the adapter turns a cone exactly; S9d.3b.1's survey, `Degenerate`) |
| `boolean/bopfuse_simple/ZG2`, `ZG4`, `ZG8`, `boptuc_simple/ZG2` | Unsupported | Viewer skipped | A box of side 4 and a frustum of radii 3 and 2 whose axis lies in the box's wall `y = 0`, its top rim through two of the box's corners: a vertex of one input on the other's face; the same turned 30 degrees about its axis: a plane through a cone's apex; with the frustum first: an edge of one input meeting an edge of the other; a frustum of radii 1 and 0.5 standing on the top face, its base circle tangent to two of the face's edges: a tangency between the inputs (S9d.3a's `Degenerate`) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZL2` to `ZL5` (16 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` of radius 4 and height 8 and a coaxial torus of radii 4 and 1 at half its height (turned about the axis by quarter turns in `ZL3` to `ZL5`): the wall through the tube's centre circle, meeting the tube in two circles at heights 3 and 5: `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9d.4b.2a's whole torus against curved faces; `bopfuse_simple/ZL2` was the sentinel for S9d.4b's `OutOfDomain`, evaluated since S9d.4b.2a, the others registered in S9d.4b's survey; their values unchanged bit for bit by S9e.3b's validator for torus bands and its bounds of torus faces, the survey of S9e.3b, S9f.2a and S9f.2b.1) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZM1`, `ZM2`, `ZM4`, `ZM5`, `ZM6` (20 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on two `pcone`s, a frustum of radii 8 and 4 and height 8 and a coaxial one of radii 2 and 1: standing on its top disc, inside it with its top disc on the top disc, from 1 below the bottom disc to the top disc, inside clear of both discs, or through both: `checkprops -s` (or `-s empty`) and the group's `checkshape`; each records a `checkview` (S9d.3b.1's two cones, which the adapter builds with ids of their own since S9d.4b's survey) |
| `boolean/bfuse_complex/J5` | Viewer skipped | Viewer skipped | Two equal cylinders crossed at right angles: their fuse's crossing ellipses (S9c.1) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `U1`, `V3`, `Y5`, `Z7`, `ZA2`, `ZA4`, `ZA7`, `ZB2`, `ZB4`, `ZB6`, `ZB9`, `ZC4`, `ZO7`, `ZO8` (56 cases); `bopcommon_simple/ZC5`, `ZD8`, `ZE1`, `boptuc_simple/ZC5`, `bopfuse_simple/ZD8`, `ZE1` | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` and a box in frames with different axes (the box turned by `trotate` about z by 30, 60 or -30 degrees, 40 of them sized by `dset`; in `ZC5` 45 degrees about a horizontal axis through its corner; in `ZO7` and `ZO8` the cylinder turned about its own axis, `ZO8` by a quarter turn), or two equal cylinders whose axes cross at right angles (`ZD8`, `ZE1`): `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9c.1's prisms with arcs in any position) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZE3`, `ZE4`, `ZE5`, `ZE6` (16 cases); `bopfuse_simple/ZC5`, `bopcut_simple/ZC5` | Viewer skipped | Viewer skipped | `bop` and its operations on a `pcylinder` of radius 0.5 through one of radius 1, their axes crossing at right angles (a quarter turn about x, which the adapter turns exactly; in `ZE4` to `ZE6` then turned about its own axis by one to three quarter turns), meeting in two quartic rings; the fuse and cut of `ZC5` (above), whose cylinder wall has holes bounded by ellipse arcs the validator now places: `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9c.2a's cylinders in exact frames) |
| `boolean/bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`, `boptuc_simple`: `ZF2`, `ZF3` (8 cases) | Viewer skipped | Viewer skipped | `bop` and its operations on two equal parallel cylinders, the second moved 1 along x, turned -120 or 120 degrees about its own axis and then 60 about z (a turned frame, the axes 1 apart), their walls meeting in generatrices at an ellipse's and a circle's crossings on the exact models: `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9c.2b.2's parallel cylinders in turned frames; `bopfuse_simple/ZF2` was S9c.2b.1's sentinel for them) |
| `boolean/bopfuse_simple/ZD9` | Unsupported | Viewer skipped | Equal cylinders whose axes cross at right angles, the second turned 60 degrees about its own axis: a frame not exactly orthonormal, the models touching where their common extent ends, two branches of the section within the resolution of a node (S9c.2b.1's `Degenerate`; S9c.2a refused it as cylinders in turned frames) |
| `boolean/bopfuse_simple/ZE7`; `bopcommon_simple`, `bopcut_simple`, `bopfuse_simple`: `T7`, `Y2` | Unsupported | Viewer skipped | A cylinder of radius 0.5 through one of radius 1 at right angles with its axis moved 0.5 off, the walls touching at a point (a figure-eight section); a box turned 135 degrees about z with two corners on a cylinder within rounding: a tangency between the inputs, a piece thinner than the resolution (S9c.2a's `Degenerate`; S9c.1 failed on `T7` and `Y2`) |
| `boolean/bopfuse_simple/U2`, `S2`, `S3`, `ZD5`, `bopcut_simple/ZD8` | Unsupported | Viewer skipped | A box turned a quarter turn, its wall tangent to a cylinder; a box turned 45 or -45 degrees, its corner on a cylinder within rounding; a cylinder on an equal one turned a quarter turn about their axis; one of two equal crossed cylinders cut from the other: a tangency between the inputs, two meetings within rounding along an arc, a piece thinner than the resolution, a meeting at every seam tried (their rims one circle), solids touching at a vertex (S9c.1's `Degenerate`) |
| `boolean/bopcommon_simple/C3`, `bopcut_simple/F6`, `G8`, `bopfuse_simple/N6` | Unsupported | Viewer skipped | A box turned by 45, 30 or 115 degrees with a corner on the other box's corner, wall or edge within rounding: a face thinner than the resolution, a face using an edge both ways, a face touching itself at a vertex, two solids touching at a point (S9b.1's `Degenerate`) |
| `boolean/bcut_simple/H4` to `L2` (35 cases) | Viewer skipped | Viewer skipped | Pockets cut from a prism one after another by prisms of profiles in planes facing z or -z (`J3`'s last tool moved by `ttranslate`): the first cut is a stack (S9a.2; in `J4` and `J7` a polyhedron, S9b.1), which the next `bcut` takes as its object, on its stored geometry (S9b.2; in `J2` to `J7` and `K7` a third `bcut` takes the second's result): `checkprops -s` and the group's `checkshape`; each records a `checkview` |
| `boolean/bopfuse_simple/H3`, `H4` | Viewer skipped | Viewer skipped | A box turned 45 degrees with its corner on the other's wall within rounding (`H3`'s inside it exactly): refused as a direction of zero length until S9b.2's face frames took each face's whole vector area |
| `boolean/bcut_simple/L3` to `L6` | Viewer skipped | Viewer skipped | DRAW's rollex: a pocket and then a hole cut from a disc prism; the first cut is a stack with cylindrical walls, which the next `bcut` takes as its object with a cylinder standing on the pocket's floor (on a plane facing down in `L3` and `L4`): `checkprops -s` and the group's `checkshape`; each records a `checkview` (S9e.2's stack given to another Boolean on its construction's curved arrangement; `L3` was the sentinel for its `OutOfDomain` since S9b.2, `L4` to `L6` registered in S9e.2's survey; the volume audit of the survey of S9e.1 and S9e.2 the reference's within 3.2e-16 relative, those of S9e.3a and S9f.1 and of S9e.3b, S9f.2a and S9f.2b.1 the same bit for bit) |
| `boolean/bfuse_simple/E1` | Unsupported | Viewer skipped | Prisms with arcs sized by `SCALE`: the tool's profile does not translate exactly into the object's frame, so the kernel decides them on exact models, where an arc must end on its circle exactly; the `profile`'s half circles end off theirs by rounding (S9c) |
| `boolean/bopcommon_simple/C8` | Unsupported | Viewer skipped | An angle `atan2(1,2)*180/pi` in `dset`, which the adapter does not evaluate |
| `boolean/bopcommon_simple/S5`, `bopcut_simple/ZC7`, `S4`, `bopfuse_simple/U7`, `B3`, `boptuc_simple/R1` | Unsupported | Viewer skipped | Box corners on or tangent to a cylinder, boxes touching along an edge: pieces thinner than the resolution, results touching themselves (at a point or along an edge) or each other, a hole touching its boundary (the kernel's `Degenerate`) |
| `boolean/bcut_simple/G8` | Unsupported | Viewer skipped | A semi-infinite prism of an exploded face (`explode`, `SemiInf`) |

There are **five original geometry tests passing on both backends** and
1,050 more evaluated on both with their image commands recorded
(`buc60769`, S7's `lowalgos/intss` cases `bug23177_1`, `bug23177_2` and
`bug24648`, 987 Boolean cases of S9a, S9a.2, S9b.1, S9b.2, S9c.1, S9c.2a,
S9c.2b.2, S9d.1, S9d.2, S9d.3a, S9d.1's pole follow-up, S9d.3b.1,
S9d.4b.2a and S9e.2, and 59 Booleans of restored solids, 16 of S9e.4a, 7
of S9e.4b.1, 4 of S9e.4b.2, 14 of S9e.4b.3c.1, 13 of S9e.4b.3c.2, 1 of
S9e.4b.3c.3a, 1 of S9e.4b.4a and 3 of S9e.4b.4b.1; S9c.2b.1, S9e.3a,
S9f.1, S9e.3b, S9f.2a, S9f.2b.1, S9f.2b.2, S9f.3a, S9f.3b, S9e.4b.3a,
S9e.4b.3b, S9e.4b.3c.3b, S9e.4b.4b.2a, S9e.4b.4b.2b.1 and S9e.4b.4c.1 add
none).
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
`vertex`, `settolerance`, `XProgress`, `dset`, `protect`) to either backend.
The Rust worker runs `dset` (DRAW's numeric variables), `protect` (a
no-op), which the `boolean` group's `begin` calls, and since S9a the
Booleans of one object and one tool (`bfuse`, `bcut`, `bcommon`, `btuc`,
`bop` and `bopfuse`, `bopcut`, `boptuc`, `bopcommon`, `bbop`, `bapibop`)
and `unifysamedom` of a Boolean result (below); it rejects the rest.

**S9a's Boolean group.** The self-contained cases of `tests/boolean`'s
`bfuse`, `bcut`, `bcommon`, `bopfuse`, `bopcut`, `bopcommon` and `boptuc`
grids (`_simple`, `_2d`, `_complex`: 1,802 of their 3,332 cases load no
data; all but two in the `_simple` grids) were run on both backends
(2026-09-28). Every case that evaluates on Rust is registered: 320, each
evaluated by native DRAW too, and each `viewer_skipped` on both, since
every one records a `checkview` (`bop` and one of its operations on two
boxes, a box and a `pcylinder`, or two cylinders in one frame, and
`bfuse_simple/L2`, `bcommon_simple/I5` and `J1`). No case fails on Rust.
Native DRAW evaluates 1,578; the host does not forward the other 224's
commands (`nurbsconvert` 96, `edge` and `wire` 108, `blend`, `orientation`
and a few more). Of the 1,258 native DRAW evaluates and Rust
does not, the kernel reports: frames with different axes (`trotate`) 647,
solids other than prisms (spheres, cones, tori) 315, S9a.2's stacks (a cut
whose tool does not span the object's heights 105, a fuse of different
heights 58), `Degenerate` 78 (a piece or result thinner than the
resolution 64, results touching themselves or each other 10, a hole
touching its boundary 4) and a profile not translating exactly 1; the
adapter does not evaluate `atan2` in `dset` (48) or `explode` a solid
without the native selector (5), nor a `pcylinder` on a plane (1). Fourteen
of them, one or two per reason, are registered as capability sentinels.
Their assertions are `checkprops -s` and `checkshape` (model-independent),
so the ledger does not change.

**S9a.2's stacks in the Boolean group.** The same 1,802 cases were run
again on both backends after S9a.2 (2026-09-28, `survey_upstream_tests.py
--boolean`, the public dataset read through `--data-dir`; no case loads
it). Native DRAW's statuses are unchanged (1,578 evaluated, 224 not
forwarded). Rust now evaluates 404, each evaluated by native DRAW too and
`viewer_skipped` on both; no case fails on Rust. Every one is registered:
84 new (`bop` with `bopcut` 47, `bopfuse` 26, `boptuc` 8, `bcut` 2,
`bfuse` 1, from `bcut_simple`, `bfuse_simple` and the `bop*_simple`
grids; each a stack: two boxes, a box and a `pcylinder`, or two cylinders
whose slabs hold different regions). Against the previous survey's 163
stacked cases, by the change of each count: 84 more evaluate, 37 cut a pocket and then give the stack to the next `bcut` (a
Boolean of a stack, S9b), 38 more are `Degenerate` and 4
more `explode` a Boolean's result, which now exists, into solids. Their assertions are again `checkprops -s` (or `-s empty`) and
`checkshape`: no upstream case in these grids counts a stack's shapes or
measures its lengths; `boolean_stacks` (below) does. Of the 1,174 native
DRAW evaluates and Rust does not, the reasons are now: frames with
different axes 647; solids other than prisms 315 (255 refused by the
kernel, 60 by the adapter, which extrudes again a tool sharing ids with
the object and finds no profile); a stack as a Boolean argument 37;
`Degenerate` 116 (a result thinner than the resolution 62, a profile piece
thinner than it 24, a result touching itself at a point 15 or along an
edge 4, two results touching 4, a result's hole outside its boundary 8);
a profile not translating exactly 1; `atan2` in `dset` 48; `explode` 9 (a
Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a
box 1); a `pcylinder` on a plane `pl1` 1. The sentinels were checked
again: `bopfuse_simple/B2` now evaluates and is registered, replaced by
`bopfuse_simple/B3` (boxes touching along an edge, `Degenerate`);
`bcut_simple/H6` is still refused, now as a Boolean of its first cut's
stack; the other twelve are refused as before.

`bopfuse_simple/Z5` (`pcylinder b1 1 2`, a box `-r -r 0 2*r 2*r 1` with
`dset r sqrt(2)/2`, `bop b1 b2`, `bopfuse result`, `checkprops result -s
18.8496`) is the cylinder (area `6 pi`). Its heights differ, and the fuse
decides containment by exact cuts: the cylinder minus the box, four
circular segments touching at the box's corners (`2 r^2` rounds just above
1), is refused as degenerate, which says it is not empty; the survey found
the first S9a.2 revision propagating that refusal, fixed with a kernel
test (`tests/booleans.rs`).

**S9b.1's polyhedra in the Boolean group.** The same 1,802 cases were run
again on both backends after S9b.1 and the adapter's polyhedra
(2026-09-28, `survey_upstream_tests.py --boolean`, the public dataset read
through `--data-dir`; no case loads it). The host now forwards `tcopy`, so
native DRAW evaluates 1,579 (`bopfuse_simple/ZP6`, three tori, besides) and
does not forward 223. Rust evaluates 704, each evaluated by native DRAW too
and `viewer_skipped` on both; no case fails on Rust. Every one is
registered: 300 new. 298 are S9b.1's polyhedra, `bop` and one of its
operations on two boxes, one turned by `trotate` (about z, 45 of them sized
by `dset`; `Q9` of each grid about x and moved by `ttranslate`), from the
four `bop*_simple` grids (`bopcommon` 81, `boptuc` 80, `bopfuse` 69,
`bopcut` 68). Eight of them turn by a quarter turn: `K3` of each grid
evaluates with the kernel's rounded turn as well, `P6` only with the
adapter's exact one (a box against the other's wall, which the rounded turn
tilts by 6.1e-17 into a face thinner than the resolution). The other two,
`bopfuse_simple/Z8` and `ZB3`, are stacks of a cylinder and a box whose
corners lie on it within rounding: the containment fix made for `Z5`
(above) reached this branch after S9a.2's survey, which had refused them
as results touching themselves at a point. Before the adapter's changes
the same kernel evaluated 700 (the four `P6` refused). Against the
previous survey's 647 cases in frames with different axes: 298 evaluate,
233 are cylinders turned, arcs in frames with different axes (S9c), 114
are S9b.1's `Degenerate` (a turned box's corner on another box's corner,
wall or edge within rounding) and 2 (`bcut_simple/J4`, `J7`) cut a pocket
by a prism in a plane facing -z, a polyhedron, and give it to the next
`bcut`. The 875 cases native DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9a.2's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a tool sharing ids with the object is extruded again, and has no profile) | 316 | +1 (`ZP6`, through `tcopy`) |
| Arcs in frames with different axes (S9c): cylinders turned by `trotate` 233, `bfuse_simple/E1` (its tool's profile does not translate exactly) 1 | 234 | new (from frames with different axes 233, a profile not translating 1) |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | new |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | new |
| `Degenerate`, S9b.1's: a direction of zero length 2, two solids touching at a point 1, a face touching itself at a vertex 1 | 4 | new |
| `Degenerate`: a result thinner than the resolution | 63 | +1 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | -4 (`Z5`, `Z8`, `ZB3` evaluate, one is a result thinner than the resolution) |
| `Degenerate`: a result touching itself along an edge | 6 | +2 (turned boxes) |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| A stack (37) or a polyhedron (2) as a Boolean argument (S9b.2) | 39 | +2 (`J4`, `J7`) |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |

The previous survey's 647 cases in frames with different axes and one
profile not translating exactly are gone from the reasons. The sentinels
were checked again: none evaluates. `bopcommon_simple/C3` (a box turned 45
degrees) is now refused as a face thinner than the resolution,
`bfuse_simple/E1` as arcs decided on exact models (S9c), `bcut_simple/H6`
as a Boolean of a stack (S9b.2); their purposes say so. Five sentinels are
added for S9b.1's refusals: `bopcut_simple/F6` (a face using an edge both
ways), `G8` (a face touching itself at a vertex), `bopfuse_simple/N6` (two
solids touching at a point), `H3` (a direction of zero length) and
`bcut_simple/J4` (a Boolean of a polyhedron).

**S9b.2's Boolean results as arguments in the Boolean group.** The same
1,802 cases were run again on both backends after S9b.2 and the adapter's
change giving stacks and polyhedra to Booleans (2026-09-28,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`; no case loads it). Native DRAW's statuses are unchanged
(1,579 evaluated, 223 not forwarded). Rust evaluates 741, each evaluated by
native DRAW too and `viewer_skipped` on both; no case fails on Rust. Every
one is registered: 37 more than S9b.1's survey (`bcut_simple/H6`, `J4` and
`bopfuse_simple/H3`, registered with the adapter's change, and 34 new). 35
are `bcut_simple/H4` to `L2`, pockets cut from a prism one after another by
prisms of profiles in planes facing z or -z: the first cut is a stack (33)
or a polyhedron (`J4`, `J7`), which the next `bcut` takes as its object on
its stored geometry; in seven (`J2` to `J7`, `K7`) a third `bcut` takes the
second's result. The other two, `bopfuse_simple/H3` and `H4`, are boxes
turned 45 degrees with a corner on the other's wall within rounding, which
S9b.1 refused as a direction of zero length: S9b.2's face frames take each
face's whole vector area. The S9b.2 kernel with the adapter before its
change evaluated 706 (`H3` and `H4` besides S9b.1's 704), every stack or
polyhedron given to a Boolean refused. Of S9b.1's 39 stacks and polyhedra
given to another Boolean, 35 evaluate and 4 (`bcut_simple/L3` to `L6`, a
pocket and then a hole cut from a disc) are stacks with cylindrical walls,
which S9b.2 refuses as curved. The kernel now gives every curved input one
reason, a solid with curved faces or edges in any position (S9c): 493
cases (the spheres, cones and tori it refused, the arcs in frames with
different axes, `bfuse_simple/E1` and these four). After the adapter's
change, two cones, spheres or tori sharing ids (61 cases) were refused as
two stacks or polyhedra sharing ids; the adapter now names a solid other
than a prism sharing ids. The 838 cases native DRAW evaluates and Rust
does not, by reason:

| Reason | Cases | Since S9b.1's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a solid other than a prism sharing ids with the other argument) | 316 | 0 |
| Arcs in frames with different axes (S9c): cylinders turned by `trotate` 233, `bfuse_simple/E1` (its tool's profile does not translate exactly) 1 | 234 | 0 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | new (from a stack as a Boolean argument) |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | -2 (`H3`, `H4` evaluate) |
| `Degenerate`: a result thinner than the resolution | 63 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |

The previous survey's 39 stacks and polyhedra as Boolean arguments are gone
from the reasons. The sentinels were checked again: `bcut_simple/H6`, `J4`
and `bopfuse_simple/H3` evaluate (registered so with the adapter's change;
their purposes say why), the others are refused as before. One sentinel is
added for S9b.2's refusal: `bcut_simple/L3` (a stack with cylindrical walls
given to a Boolean, S9c).

**S9c.1's prisms with arcs in any position in the Boolean group.** The same
1,802 cases were run again on both backends after S9c.1 (2026-09-28,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`; no case loads it). Native DRAW's statuses are unchanged
(1,579 evaluated, 223 not forwarded). Rust evaluates 804, each evaluated by
native DRAW too and `viewer_skipped` on both, and fails 8 (below). Every
case it evaluates is registered: 63 more than S9b.2's survey
(`bfuse_complex/J5`, registered with S9c.1's kernel, and 62 new), all from
the 234 refused as arcs in frames with different axes. 56 are `bop` and its
four operations on a `pcylinder` and a box turned by `trotate` about z (30,
60 or -30 degrees, 40 sized by `dset`) or a cylinder turned about its own
axis (`ZO7`, `ZO8`); 2 a box turned 45 degrees about a horizontal axis
through its corner (`bopcommon_simple/ZC5`, `boptuc_simple/ZC5`); 4 the
fuse and common of two equal cylinders whose axes cross at right angles
(`ZD8`, `ZE1`: S9c.1's crossing ellipses, as `J5`). Of the other 171: 48
are two cylinders S9c.2 takes (radii 1 and 0.5 crossing, 32; equal
cylinders whose turned frame is not exactly orthonormal, crossing, 8, or
parallel, 8: S9c.1 takes two cylinders only when circular in a common
measure), 114 S9c.1's or S9a's `Degenerate` (a box's wall tangent to the
cylinder or its corner on it within rounding, two stacked cylinders' rims
on one circle, the cuts of crossed cylinders touching at a vertex), 1
(`bfuse_simple/E1`) an arc ending off its circle on exact models, and 8
fail. The four stacks with cylindrical walls given to another Boolean
(`bcut_simple/L3` to `L6`) are still refused (S9c). The 775 cases native
DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9b.2's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a solid other than a prism sharing ids with the other argument) | 316 | 0 |
| Two cylinders meeting in curves other than lines and conics (S9c.2): radii 1 and 0.5 crossing 32, equal cylinders with a frame turned off whole quarter turns crossing 8 or parallel 8 | 48 | new (from arcs in frames with different axes) |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | new (from arcs in frames with different axes) |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 34 | new |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | new |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | new |
| `Degenerate`, S9c.1's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn) | 8 | new |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | new |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | +40 (a box's wall tangent to the cylinder within rounding) |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |
| Rust fails: `InvalidTopology` (a hole outside every piece 6, `uncertified_containment` 2) | 8 | new |

The previous survey's 234 arcs in frames with different axes are gone from
the reasons. The eight failures are not registered (the contract admits no
failing status); both are S9c.1 kernel gaps. `bopfuse_simple/T7`, `Y2` and
their `bopcut` and `bopcommon` twins turn a box 135 degrees about z with
two corners on the cylinder within rounding (`dset r sqrt(2)/2`, `2 r^2`
rounding above 1): `curved/graph.rs` finds a hole in no piece
(`InvalidTopology`) where `boptuc` (the box first) and the box turned by
134.9 or 135.0001 degrees refuse a piece thinner than the resolution
(134.9999 fails too), and the unturned box a profile piece thinner than
it; a sliver's binary64 image taken the wrong way round, whose neighbours
share its edges, would do it. `bopfuse_simple/ZC5` and
`bopcut_simple/ZC5`: the result's cylinder wall outside the box has three
holes bounded by ellipse arcs, and the validator's test of a hole in a band
(`signed_cover_crossings`) counts crossings with line pcurves only
(`uncertified_containment`); with the cylinder's radius 1.05 (one hole)
both backends give 20.7828 and 29.1957. The sentinels were checked again:
`bfuse_simple/E1` is now refused as an arc ending off its circle (its
purpose says so), the others as before. Six sentinels are added for
S9c.1's refusals: `bopfuse_simple/ZE3` (S9c.2's cylinders), `U2` (a
tangency), `S2` (two meetings within rounding along an arc), `S3` (a piece
thinner than the resolution), `ZD5` (a meeting at every seam tried) and
`bopcut_simple/ZD8` (solids touching at a vertex).

**S9c.2a's cylinders in exact frames in the Boolean group.** The same
1,802 cases were run again on both backends after S9c.2a (2026-09-28,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`; no case loads it). Native DRAW's statuses are unchanged
(1,579 evaluated, 223 not forwarded). Rust evaluates 822, each evaluated by
native DRAW too and `viewer_skipped` on both; no case fails on Rust. Every
one is registered: 18 more than S9c.1's survey (`bopfuse_simple/ZE3`, a
sentinel until now, and 17 new). 16 are `ZE3` to `ZE6` of the four
`bop*_simple` grids: a `pcylinder` of radius 0.5 through one of radius 1,
their axes crossing at right angles (a quarter turn about x, which the
adapter turns exactly, then none to three quarter turns about its own
axis), meeting in two quartic rings; Rust's areas agree with native DRAW's
to its six digits (fuse 38.0580, cut 35.6606, common 7.49504, `boptuc`
9.89251). The other two are `bopfuse_simple/ZC5` and `bopcut_simple/ZC5`
(19.7221, 27.1316), S9c.1's failures of `uncertified_containment`: the
validator now certifies holes in a band bounded by sinusoid and projection
pcurves. Of the previous survey's 48 cylinders S9c.2 takes, 16 evaluate,
16 (`ZE7` to `ZF1`, the thin cylinder's axis moved 0.5 off, so the walls
touch at one point where the section crosses itself) are a tangency
between the inputs, and 16 are in turned frames (`ZD9`, `ZE2`: equal
cylinders crossing at right angles, the second's frame turned 60 degrees
about its axis, in `ZE2` then 120 about z; `ZF2`, `ZF3`: equal parallel
cylinders, the second moved by turns of -120 or 120 and 60 degrees about z;
S9c.2b).
The 32 counted as quartics before were `ZE3` to `ZF1`; the second half's
quartics have a node. S9c.1's other six failures (`T7`, `Y2` of
`bopfuse_simple`, `bopcut_simple` and `bopcommon_simple`) are now refused
as a piece thinner than the resolution, as `boptuc_simple`'s twins were:
S9c.2a refuses a hole that no piece holds as a sliver. The 757 cases native
DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9c.1's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a solid other than a prism sharing ids with the other argument) | 316 | 0 |
| Two cylinders in turned frames meeting in curves other than lines and conics (S9c.2b): equal cylinders with a frame turned off whole quarter turns crossing 8 or parallel 8 | 16 | -32 (16 evaluate, 16 a tangency) |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | +6 (`T7`, `Y2`, failures before) |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | +16 (`ZE7` to `ZF1`) |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |

The previous survey's eight failures are gone. The sentinels were checked
again: `bopfuse_simple/ZE3` evaluates and is registered so (its purpose
says why), the others are refused as before. Eight sentinels are added:
`bopfuse_simple/ZD9` (cylinders in turned frames, S9c.2b, replacing `ZE3`
for the reason), `ZE7` (S9c.2a's tangency: the thin cylinder's interval
across the axes ending on the thick one's) and the six `T7` and `Y2` (a
piece thinner than the resolution, where S9c.1 failed).

**S9c.2b.1's cylinders in turned frames in the Boolean group.** The same
1,802 cases were run again on both backends after S9c.2b.1 (2026-09-28,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`; no case loads it). Native DRAW's statuses are unchanged
(1,579 evaluated, 223 not forwarded). Rust evaluates 822, as before, each
evaluated by native DRAW too, `viewer_skipped` on both and registered; no
case fails on Rust and no status changes. Of S9c.2a's 16 cylinders in
turned frames, none evaluates: `ZD9` and `ZE2` of the four `bop*_simple`
grids (equal cylinders crossing at right angles, the second's frame
turned 60 degrees about its axis, in `ZE2` then 120 about z) are refused
as `Degenerate`, by decision: their axes meet, so the models touch where
their common extent ends, and two branches of the section lie within the
resolution of a node (the stored Steinmetz pairs' case). `ZF2` and `ZF3`
(equal parallel cylinders, the second moved by turns of -120 or 120 and
60 degrees about z) are refused as before, parallel cylinders in turned
frames being S9c.2b.2's; the message is still S9c.2's general one (two
cylinders meeting in curves other than lines and conics). No other
reason changes. The 757 cases native DRAW evaluates and Rust does not, by
reason:

| Reason | Cases | Since S9c.2a's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a solid other than a prism sharing ids with the other argument) | 316 | 0 |
| Parallel cylinders in turned frames (S9c.2b.2): equal cylinders, one moved and turned off whole quarter turns (`ZF2`, `ZF3`) | 8 | -8 (`ZD9`, `ZE2` now `Degenerate`) |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | +8 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |

The sentinels were checked again: `bopfuse_simple/ZD9` is refused as a
near node now (its purpose says so), the others as before. One sentinel
is added: `bopfuse_simple/ZF2` (parallel cylinders in turned frames,
S9c.2b.2, the reason `ZD9` held until now). The ledger does not change.

**S9c.2b.2's algebraic vertices in the Boolean group.** The same 1,802
cases were run again on both backends after S9c.2b.2 (2026-09-28,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`; no case loads it). Native DRAW's statuses are unchanged
(1,579 evaluated, 223 not forwarded). Rust evaluates 830, each evaluated
by native DRAW too, `viewer_skipped` on both, and each is registered: 8
more than S9c.2b.1's survey, `ZF2` and `ZF3` of the four `bop*_simple`
grids (equal parallel cylinders, the second moved 1 along x and turned by
-120 or 120 degrees about its own axis and then 60 about z, so the axes
are 1 apart). On their exact models the turned wall's cross-section is an
ellipse, and the walls meet in generatrices through its crossings with the
first cylinder's circle. The areas agree with native DRAW's to its printed
digits (fuse 43.62, cut and `btuc` 28.9592, common 19.2119) and with the
closed forms of two unit discs 1 apart (lens `2 pi / 3 - sqrt(3) / 2`,
arcs of `2 pi / 3` inside the other disc) within 6e-15. No case fails on
Rust, and no other status or reason changes. S9c.2's general message
(two cylinders meeting in curves other than lines and conics) no longer
occurs among the Boolean cases.
The 749 cases native DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9c.2b.1's survey |
| --- | --- | --- |
| Solids other than prisms (spheres, cones, tori): refused by the kernel 255, by the adapter 61 (a solid other than a prism sharing ids with the other argument) | 316 | 0 |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | 0 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |
| A `pcylinder` on a plane `pl1` | 1 | 0 |

Parallel cylinders in turned frames (S9c.2b.2), 8 cases before, are gone
from the reasons. The sentinels were checked again:
`bopfuse_simple/ZF2` evaluates and is registered so (its purpose says
why); no sentinel replaces it, since its reason no longer occurs. The
others are refused as before. The ledger does not change.

**S9d.1's spheres in the Boolean group.** The same 1,802 cases were run
again on both backends after S9d.1 (2026-09-29, `survey_upstream_tests.py
--boolean`, the public dataset read through `--data-dir`; no case loads
it). Native DRAW's statuses are unchanged (1,579 evaluated, 223 not
forwarded). The kernel now takes a sphere against a polyhedral prism, and
none of the group's 101 Booleans of a sphere refused until now as solids
other than prisms evaluates: 33 are a unit box whose corner is at a unit
sphere's centre (`A1` to `A5`, `D3` to `D8` of `bfuse_simple` and
`bcommon_simple`, `A1` to `A5`, `F9` to `G5` of `bcut_simple`; the sphere
or the box turned by `trotate`), whose far corners lie on the sphere where
its far walls touch it, and 68 a sphere against a cylinder. S9d.1's decisions
refuse both: an input's vertex on the sphere as `Degenerate` (reported as
a meeting at every seam tried: the arrangement takes a vertex on a
sphere's face for a seam conflict, as on a cylinder's, and retries it at
each seam), a wall tangent to it as a tangency between the inputs (in the
12 with the box turned a quarter turn about y, its wall `x = -1` tangent
to the sphere is met before its corners), and a sphere against a cylinder as S9d.2's `OutOfDomain`. With
the cylinder first the 60 of `ZH5` to `ZJ3` (a cylinder of radius 4 and
height 8, a sphere of radius 4, 2 or 6 centred on its top cap, turned
about the cap's centre) are S9d.2's; with the sphere first
(`boptuc_simple`), the eight of radius 4 (`ZH5` to `ZI3`) are refused
before the pair is: the cylinder's rim, the sphere's equator, lies on the
sphere at every seam tried (the wall tangent to the sphere along it).
The one sphere against a box whose vertices stay off the sphere,
`bopcommon_simple/ZP9` (a 100 box and a sphere of radius 7.5 placed on a
DRAW `plane`, its centre 1.053 inside the wall `x = 100`), was refused by
the adapter, which placed no sphere on a plane (the reason recorded until
now as a `pcylinder` on a plane `pl1`: it is this `psphere`). The adapter
now takes `psphere name plane R [angle1 angle2]` on the plane's frame
(below); the case evaluates, `viewer_skipped` on both, and is registered.
Its area agrees with native DRAW's to its printed digits (576.293) and
with the closed form (the sphere less the cap beyond the wall, plus the
section's disc) within 5.1e-13 (9e-16 relative). Rust evaluates 831, each
evaluated by native DRAW too, `viewer_skipped` on both, and each
registered. No case fails on Rust; no other status changes. The 748 cases
native DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9c.2b.2's survey |
| --- | --- | --- |
| Solids other than prisms (cones, tori): refused by the kernel 154, by the adapter 61 (a solid other than a prism sharing ids with the other argument: two cones 60, two tori 1) | 215 | -101 (spheres, below) |
| A sphere against a cylinder (S9d.2, `OutOfDomain`): `ZH5` to `ZJ3` of the four `bop*_simple` grids, the cylinder first | 60 | new (from solids other than prisms) |
| `Degenerate`, S9d.1's: a meeting at every seam tried, a vertex or an edge of the other input on the sphere (a unit box's corner at a unit sphere's centre, its far corners on the sphere, 21; a cylinder's rim on a sphere's equator, the sphere first, `boptuc_simple/ZH5` to `ZI3`, 8) | 29 | new (from solids other than prisms) |
| `Degenerate`, S9d.1's: a tangency between the inputs (the unit box quarter-turned, its wall `x = -1` tangent to the unit sphere) | 12 | new (from solids other than prisms) |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | 0 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |

The previous table's `pcylinder` on a plane (1) is gone: it was
`bopcommon_simple/ZP9`'s `psphere`, which evaluates. The sentinels were
checked again: `bcommon_simple/A1` (a sphere and a box) is refused as a
meeting at every seam tried now, its purpose updated to say why; the
others are refused as before. Three sentinels are added:
`bopfuse_simple/ZH5` (a sphere against a cylinder, S9d.2),
`bfuse_simple/A4` (a wall tangent to a sphere) and `boptuc_simple/ZH5`
(the same sphere and cylinder as `bopfuse_simple/ZH5` with the sphere
first: the rim on the sphere, refused before the pair). The ledger does
not change.

**S9d.2's spheres against cylinders in the Boolean group.** The same
1,802 cases were run again on both backends after S9d.2a and S9d.2b
(2026-09-29, `survey_upstream_tests.py --boolean`, the public dataset
read through `--data-dir`; no case loads it). Native DRAW's statuses are
unchanged (1,579 evaluated, 223 not forwarded). Of the 60 cases of a
sphere against a cylinder refused until now as S9d.2's `OutOfDomain`
(`ZH5` to `ZJ3` of the four `bop*_simple` grids: a `pcylinder` of radius
4 and height 8, a sphere of radius 4, 2 or 6 centred on its top cap,
turned about the cap's centre by `trotate`), 20 evaluate: the sphere of
radius 6 (`ZI8` to `ZJ3`), which holds the cap and meets the wall in a
parallel at height `8 - 2 sqrt(5)`, a coaxial sphere's ring (S9d.2a),
unturned, turned by whole quarter turns about z or a quarter turn about
y. Each is `viewer_skipped` on both and registered; the areas agree with
native DRAW's to its printed digits (fuse 533.721, common 220.262, cut
196.529, `btuc` 557.453) and with the closed forms (fuse `pi (152 + 8
sqrt(5))`, common `pi (88 - 8 sqrt(5))`, cut `pi (152 - 40 sqrt(5))`,
`btuc` `pi (88 + 40 sqrt(5))`: the sphere's zones by Archimedes, the
wall's bands and the discs) within 3.4e-13 (9.1e-16 relative). The other
40 are refused for reasons S9d.2 did not name. The 24 of radius 4 with
the cylinder first, and with them the 8 of `boptuc_simple` refused before
as a meeting at every seam tried, are refused as a computation limit: the
cylinder's rim is the sphere's equator, the wall tangent to the sphere
along it, so the coaxial pair's discriminant vanishes identically and no
chart with a negative point exists: the search for one gives up before
the tangency test runs (a tangency along a circle, which S9d.2's
decisions name `Degenerate`; the kernel now tests a discriminant vanishing identically first, and these 32 are refused as a tangency between the inputs, `Degenerate`). The 16 of radius 2 (`ZI4` to `ZI7`) are turned
a quarter turn about x first, so the cap's plane holds the sphere's axis
and the section circle, a great circle inside the cap, runs through the
sphere's poles, the singular points of its parameterization: they are
refused as `PrecisionLoss` (coordinates cannot resolve the requested
linear tolerance). The same sphere unturned evaluates, as does a box
against it turned 30 degrees about x; a box's wall through an upright
sphere's centre (a meridian plane) is refused alike, so this is a sphere
against any plane through its poles, S9d.1's domain too. The 21 cases of
a unit box's corner at a unit sphere's centre refused as a meeting at
every seam tried are refused as a vertex of one input on the other's face
now (S9d.2a's amendment (f)). Rust evaluates 851, each evaluated by native
DRAW too, `viewer_skipped` on both, and each registered. No case fails
on Rust and none times out; no other status changes. The 728 cases native
DRAW evaluates and Rust does not, by reason:

| Reason | Cases | Since S9d.1's survey |
| --- | --- | --- |
| Solids other than prisms (cones, tori): refused by the kernel 154, by the adapter 61 (a solid other than a prism sharing ids with the other argument: two cones 60, two tori 1) | 215 | 0 |
| `Degenerate`, S9d.2's: a tangency between the inputs, a coaxial cylinder and sphere whose discriminant vanishes identically (a cylinder's rim on a sphere's equator, the wall tangent to it along the circle: `ZH5` to `ZI3` of the four `bop*_simple` grids; the survey found them refused as a computation limit, corrected after it) | 32 | new (24 from a sphere against a cylinder, 8 from a meeting at every seam tried) |
| `PrecisionLoss`: a section circle through a sphere's poles (a sphere of radius 2 turned so the cap's plane holds its axis: `ZI4` to `ZI7` of the four grids) | 16 | new (from a sphere against a cylinder) |
| `Degenerate`, S9d.2a's: a vertex of one input on the other's face (a unit box's corner at a unit sphere's centre, its far corners on the sphere) | 21 | new (from S9d.1's meeting at every seam tried) |
| `Degenerate`, S9d.1's: a tangency between the inputs (the unit box quarter-turned, its wall `x = -1` tangent to the unit sphere) | 12 | 0 |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | 0 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcylinder` into faces 4, a face of a box 1) | 9 | 0 |

The previous table's sphere against a cylinder (S9d.2, 60) and S9d.1's
meeting at every seam tried with a vertex or an edge on the sphere (29)
are gone. The sentinels were checked again: `bopfuse_simple/ZH5` is
refused as the computation limit now, its purpose updated to say why;
`bcommon_simple/A1` as a vertex on the other's face, its purpose updated;
`boptuc_simple/ZH5` is dropped, its reason gone and its new one
`bopfuse_simple/ZH5`'s; the others are refused as before. One sentinel
is added: `bopfuse_simple/ZI4` (a section through a sphere's poles,
`PrecisionLoss`). The ledger does not change.

**S9d.3a's cones against polyhedral prisms in the Boolean group.** The
same 1,802 cases were run again on both backends after S9d.3a
(2026-09-29, `survey_upstream_tests.py --boolean`, the public dataset
read through `--data-dir`; no case loads it). Native DRAW's statuses are
unchanged (1,579 evaluated, 223 not forwarded). The kernel now takes a
cone or frustum against a polyhedral prism. Of the 154 cases the kernel
refused until now as solids other than prisms, 138 hold a cone: 72 a box
of side 4 and a `pcone` frustum (`ZF5` to `ZH4` of the four `bop*_simple`
grids) and 66 a cone against a cylinder. 36 of the 72 evaluate: a frustum
of radii 1 and 0.5 on the vertical through the box's centre, standing on
its top face, inside it with its top disc on the top face or clear of
every face, from the bottom face to the top face, or from 1 below the
bottom face to the top face (`ZF5` to `ZF9`: circles on the box's
faces), and a frustum of radii 5 and 4 (5 and 3.5 in `ZH4`) of the box's
height whose axis is 2 outside a wall and 2 from the walls across it,
meeting them in hyperbolas (`ZH1` to `ZH4`; in `ZH3` and `ZH4` turned 30
degrees about its axis, a turned frame). Each is `viewer_skipped` on both
and registered. Their areas agree with native DRAW's to its printed
digits (`ZH1`: fuse 273.296, cut 62.1124, common 68.087, `btuc`
241.753), and so do their volumes (the cases run again with `vprops`
added); both agree with closed forms within 2.3e-13 in volume and 8.5e-14
in area (8.1e-16 relative): the frusta's volumes `pi h (R^2 + R r + r^2) /
3` and areas for `ZF5` to `ZF9`, and for `ZH1` to `ZH4` the section of
the box inside the disc of radius `r` about the axis, `2 sqrt(r^2 - 4) +
r^2 asin(2/r) - 8`, the walls' strips between the hyperbolas, `sqrt(r^2 -
4) - 2` wide, and the wall's patch inside the box, `2 r asin(2/r)` per
unit of slant height, integrated over the radius in closed form. The
other 36 are refused as `Degenerate`, each for what its geometry holds:
the frustum's axis lies in a wall of the box (a plane through the apex,
S9d.3's decisions) and a rim of radius 2 runs through two of the box's
corners, tangent there to two of its edges (`ZG1` to `ZG7`), or a base
circle on the box's top face is tangent to two of the face's edges
(`ZG8`, `ZG9`). The reason reported is the first the arrangement meets: a
tangency between the inputs (24), a vertex of one input on the other's
face (the corners on the rim, `ZG2` and `ZG6` with the box first, 6), a
plane through a cone's apex (`ZG4`, `ZG2`'s frustum turned 30 degrees
about its axis, 4), an edge of one input meeting an edge of the other
(`boptuc_simple/ZG2` and `ZG6`, the frustum first, 2). The 66 cones
against cylinders (`ZJ4` to `ZL1` of the four grids: a `pcylinder` of
radius 4 and height 8 and a coaxial cone or frustum, or one turned a
quarter turn about x or y; `bcut_simple/G9` and `H3`, a frustum fused onto
a cylinder) are refused as S9d.3b's `OutOfDomain` now. The 32 cases of a
cylinder's rim on a sphere's equator are refused as a tangency between
the inputs, as the correction after S9d.2's survey says. Rust evaluates
887, each evaluated by native DRAW too, `viewer_skipped` on both, and
each registered. No case fails on Rust and none times out; no other
status changes. The 692 cases native DRAW evaluates and Rust does not, by
reason:

| Reason | Cases | Since S9d.2's survey |
| --- | --- | --- |
| Solids other than prisms (tori; two cones or two tori sharing ids): refused by the kernel 16 (a `pcylinder` and a coaxial torus, `ZL2` to `ZL5` of the four grids), by the adapter 61 (a solid other than a prism sharing ids with the other argument: two cones 60, two tori 1) | 77 | -138 (cones: 36 evaluate, 36 `Degenerate` and 66 S9d.3b's, below) |
| A cone against a prism with arcs (S9d.3b, `OutOfDomain`): a `pcylinder` and a cone or frustum (`ZJ4` to `ZL1` of the four grids, 64; a frustum fused onto a cylinder, `bcut_simple/G9`, `H3`, 2) | 66 | new (from solids other than prisms) |
| `Degenerate`, S9d.3a's: a tangency between the inputs (a box and a frustum whose rim is tangent to two of the box's edges at its corners, its axis in a wall: `ZG1`, `ZG3`, `ZG5`, `ZG7`; a base circle on the box's top face tangent to two of the face's edges: `ZG8`, `ZG9`) | 24 | new (from solids other than prisms) |
| `Degenerate`, S9d.3a's: a vertex of one input on the other's face (the box's corners on the frustum's rim, its axis in the wall `y = 0`: `ZG2`, `ZG6` of the three grids with the box first) | 6 | new (from solids other than prisms) |
| `Degenerate`, S9d.3a's: a plane through a cone's apex (`ZG4`: `ZG2`'s frustum turned 30 degrees about its axis) | 4 | new (from solids other than prisms) |
| `Degenerate`, S9d.3a's: an edge of one input meeting an edge of the other (`boptuc_simple/ZG2`, `ZG6`: the frustum first, the box's edges through its rim) | 2 | new (from solids other than prisms) |
| `Degenerate`, S9d.2's: a tangency between the inputs, a coaxial cylinder and sphere whose discriminant vanishes identically (a cylinder's rim on a sphere's equator, the wall tangent to it along the circle: `ZH5` to `ZI3` of the four `bop*_simple` grids) | 32 | 0 (this survey found them so; S9d.2's found a computation limit, corrected after it) |
| `PrecisionLoss`: a section circle through a sphere's poles (a sphere of radius 2 turned so the cap's plane holds its axis: `ZI4` to `ZI7` of the four grids) | 16 | 0 |
| `Degenerate`, S9d.2a's: a vertex of one input on the other's face (a unit box's corner at a unit sphere's centre, its far corners on the sphere) | 21 | 0 |
| `Degenerate`, S9d.1's: a tangency between the inputs (the unit box quarter-turned, its wall `x = -1` tangent to the unit sphere) | 12 | 0 |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | 0 |
| A stack with cylindrical walls as a Boolean argument (S9c) | 4 | 0 |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcone` into faces 4, a face of a box 1) | 9 | 0 |

The explode row's four faces are of a `pcone` (`explode pc f` in `ZP3`,
`ZP4`, `ZP7`, `ZP8`), which the earlier tables named a `pcylinder`. The
sentinels were checked again: each is refused as before
(`bopfuse_simple/ZH5` as a tangency between the inputs, as its purpose
says since S9d.2's correction). Six sentinels are added:
`bopfuse_simple/ZJ4` (a cone against a cylinder, S9d.3b's
`OutOfDomain`), `bopfuse_simple/ZG8` (a tangency), `bopfuse_simple/ZG2`
(a vertex on the other's face), `bopfuse_simple/ZG4` (a plane through a
cone's apex), `boptuc_simple/ZG2` (an edge meeting an edge) and
`bopfuse_simple/ZL2` (a cylinder and a coaxial torus, the kernel's
remaining solids other than prisms). The ledger does not change.

After S9d.3a's survey two kernel changes landed: two parallel faces
within the resolution are `Degenerate`, and a section through a stored
sphere's pole, or within the resolution of it, gets a vertex there
(S9d.1's follow-up; a section through a pole off the axis's planes is
`OutOfDomain`). The 1,802 cases were run again (2026-09-29). Native
DRAW's statuses are unchanged. Only the 16 `PrecisionLoss` cases change
(`ZI4` to `ZI7` of the four grids, a sphere of radius 2 centred on a
`pcylinder`'s top cap, turned so the cap's plane holds its axis); no case
is refused as parallel faces within the resolution or as a section
through a pole off the axis's planes. 12 evaluate and are registered
(`ZI4`, `ZI6`, `ZI7`; `bopfuse_simple/ZI4`, the sentinel for the
`PrecisionLoss`, re-purposed): areas and volumes native DRAW's to its
printed digits and the closed forms' (fuse and cut area `100 pi`, volume
`128 pi +- 16 pi / 3`; common and `btuc` area `12 pi`, volume `16 pi /
3`) within 1.1e-13 (1.1e-15 relative). `ZI5` (a quarter turn about x and
then one about y, the sphere's stored axis `-y` and reference direction
`-z`; `ZI7` differs only in the reference direction `+z`) is wrong on
Rust and is not registered. Its fuse, cut and common evaluate with areas
native DRAW's (the `checkprops -s` the cases check) and `checkshape`
valid, but volumes off by `32 pi / 9`: the common 27.9253 against
`16 pi / 3` (16.7552 natively), its centre of gravity `(0, -1.333, 8)`
against `(0, 0, 7.25)`; the fuse 407.709 and the cut 374.199 against
418.879 and 385.369. Its `btuc` fails (`invalid topology:
uncertified_shell_orientation`: the kernel's validator cannot certify the
result's shell orientation) where native DRAW evaluates. A box in place
of the cylinder (`box b1 -4 -4 0 8 8 8`, S9d.1's domain) gives the same:
`bcommon` of it and the sphere 27.9253, `bcut` of the sphere and it the
validation failure. Rust evaluates 902 (887 before), 899 registered, and
fails one; 676 cases native DRAW evaluates are refused (692 before), the
`PrecisionLoss` row of the table above gone. The ledger does not change.

`ZI5`'s error was S9d.1's, older than the pole follow-up: a sphere face's
closing chord at a pole between two meridians was enclosed over its `u`
hull, the volume's enclosure 26 wide and its midpoint, the reported
volume, off (`REVIEW_NOTES.md`). After its fix the 1,802 cases were run
again: only the 16 former `PrecisionLoss` cases differ from S9d.3a's
survey, all 16 evaluating; `ZI5`'s four agree with native DRAW's
volumes and centres of gravity to its printed digits and with the closed
forms (volumes within 1.1e-15 relative, centres within 8.9e-16), its
`btuc` evaluates, and all four are registered. Rust evaluates 903, each
registered, and fails none; 676 are refused. A volume audit followed,
since the registrations check areas: every registered case of the group
(903) was run again on both backends with `vprops` before each
`checkprops`. Rust's volumes and centres of gravity agree with native
DRAW's to its printed digits or 1e-6 relative in all but eight, and in
those eight native DRAW is off, not Rust. In seven of `ZE3` to `ZE6` (a
cylinder of radius 0.5 through one of radius 1 at right angles, two
quartic rings) native DRAW's volumes differ from the closed form (the
common `4 int sqrt(1/4 - x^2) sqrt(1 - x^2) dx` over `[-1/2, 1/2]`,
1.5200399881729; the others `pi`, `4 pi` and `5 pi` less it) by up to
7.8e-6 relative, varying with the turn about the axis, and its centres
lie up to 4.0e-6 off the symmetry axes, while Rust's volumes are within
2.0e-15 of it. In `bopcut_simple/A9` (a unit box less a slab across its
middle, two slabs) native DRAW's `vprops` of the result puts its centre
of gravity at `(0.53125, 0.515625, 0.46875)`, though its two solids' own
are `(0.5, 0.125, 0.5)` and `(0.5, 0.875, 0.5)`; Rust's is `(0.5, 0.5,
0.5)`. The ledger does not change.

**S9d.3b.1's cones against cylinders in the Boolean group.** The same
1,802 cases were run again on both backends after S9d.3b.1 (2026-09-29,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`, 120 seconds a case; no case loads it). Native DRAW's
statuses are unchanged (1,579 evaluated, 223 not forwarded). Only the 66
cases refused until now as S9d.3b's `OutOfDomain` change (`ZJ4` to `ZL1`
of the four `bop*_simple` grids: a `pcylinder` of radius 4 and height 8
and a cone or frustum; `bcut_simple/G9`, `H3`); none is refused as
S9d.3b.2's loops. 36 evaluate on Rust: a coaxial frustum of radii 2 and 1
standing on the top cap, inside the cylinder with its discs on both caps
or clear of them, or through one cap or both (`ZJ4`, `ZJ6` to `ZJ9`:
circles), and a frustum of height 10 and radii 2 and 1 (`ZK5`, `ZK6`) or
6 and 1 (`ZK7`, `ZK8`, its wide end through both caps too) whose axis
crosses the cylinder's at right angles at half its height, turned a
quarter turn about y or x (two quartic rings over the frustum's angle).
The cases were run again with `vprops` and `sprops` before each
`checkprops`. Rust's volumes, areas and centres of gravity agree with
closed forms for the coaxial ones (the solids of revolution integrated
piecewise) and with `cones_boolean_reference.py` for the crossing ones
(the cylinder a prism of a circle, the frustum in the turned frame)
within 6.5e-16 relative in volume and area and 2.0e-15 in the centres.
Native DRAW's agree with them to its printed digits on the coaxial ones
and are off by up to 7.8e-6 relative in volume, 7.2e-6 in area and
4.3e-6 in the centres on the crossing ones (BRepGProp on approximated
quartics, as S9d.3b's capture found: `ZK7`'s common 256.47 and 216.599
natively, 256.4710 and 216.6006 exactly); the cases' own expected areas
differ too (`ZK5` 324.106 and `ZK6` 324.207 for congruent fuses,
324.2033 exactly), within `checkprops`' tolerance. 28 are registered
(`ZJ4`, `ZJ6` to `ZJ9`, `ZK5`, `ZK6` of the four grids), each
`viewer_skipped` on both; `bopfuse_simple/ZJ4`, the sentinel for S9d.3b's
`OutOfDomain`, is re-purposed. `ZK7` and `ZK8` are right but slow: run
alone on the debug worker the contract uses, they take 63 to 118 seconds
(`ZK5` and `ZK6` 15 to 25, the coaxial ones under a second), past the
contract's 30 seconds a case; five of the eight timed out in the survey,
eight cases running at once, and none is registered. The other 30 are
refused as `Degenerate` or at a second Boolean. `ZJ5` (a frustum of
radii 1 and 2 inside the cylinder, its top disc on the top cap) has its
virtual apex on the bottom cap's plane: a plane through a cone's apex.
`ZK1` (a frustum of radii 4 and 2 standing on the top cap, its base rim
the cap's rim) is a tangency between the inputs. `ZK2` to `ZK4` (the same
turned about its axis by quarter turns) are a piece thinner than the
resolution, and so is `ZK1` turned a whole turn: the adapter turns a
prism by quarter turns exactly (`BOOLEAN.md`) but a cone by the kernel's
rotation, whose sine and cosine of a quarter turn round, so the rim lies
within rounding of the cap's rim. `ZK9` and `ZL1` (a frustum of radii 1
and 6 and height 8 across the cylinder, turned a quarter turn about y)
are a plane within rounding of a cylinder's direction for the same
reason: the frustum's axis is `(1, 0, 6.1e-17)`, and its disc of radius
6, 3 from the cylinder's axis, cuts the wall. Native DRAW's tolerances
absorb both. `bcut_simple/G9` and `H3` fuse a coaxial frustum onto a
cylinder, which evaluates (volume 881.282 and area 708.632 on both), and
then cut a rod from the fuse: a stack with curved walls given to another
Boolean (S9c). Rust evaluates 934 (903 before), 931 registered, and fails
none; five time out; 640 cases native DRAW evaluates are refused (676
before). The 645 cases native DRAW evaluates and Rust does not, by
reason:

| Reason | Cases | Since S9d.3a's survey |
| --- | --- | --- |
| Solids other than prisms (tori; two cones or two tori sharing ids): refused by the kernel 16 (a `pcylinder` and a coaxial torus, `ZL2` to `ZL5` of the four grids), by the adapter 61 (a solid other than a prism sharing ids with the other argument: two cones 60, two tori 1) | 77 | 0 |
| Timed out at the survey's 120 seconds, eight cases at once (`ZK7`'s cut and common, `ZK8`'s fuse, cut and common: a frustum of radii 6 and 1 across the cylinder, evaluating alone in 63 to 118 seconds, their measures the reference's) | 5 | new (from S9d.3b's) |
| `Degenerate`, S9d.3a's: a tangency between the inputs (a box and a frustum whose rim is tangent to two of the box's edges at its corners, its axis in a wall: `ZG1`, `ZG3`, `ZG5`, `ZG7`; a base circle on the box's top face tangent to two of the face's edges: `ZG8`, `ZG9`) | 24 | 0 |
| `Degenerate`, S9d.3a's: a vertex of one input on the other's face (the box's corners on the frustum's rim, its axis in the wall `y = 0`: `ZG2`, `ZG6` of the three grids with the box first) | 6 | 0 |
| `Degenerate`, S9d.3a's: a plane through a cone's apex (`ZG4`: `ZG2`'s frustum turned 30 degrees about its axis, 4; a frustum of radii 1 and 2 inside the cylinder whose virtual apex lies on the bottom cap's plane, `ZJ5`, 4) | 8 | +4 (from S9d.3b's) |
| `Degenerate`, S9d.3a's: an edge of one input meeting an edge of the other (`boptuc_simple/ZG2`, `ZG6`: the frustum first, the box's edges through its rim) | 2 | 0 |
| `Degenerate`, S9d.3b.1's survey: a tangency between the inputs (`ZK1`: a frustum of radii 4 and 2 standing on the cylinder's top cap, its base rim the cap's rim) | 4 | new (from S9d.3b's) |
| `Degenerate`, S9d.3b.1's survey: a piece thinner than the resolution (`ZK2` to `ZK4`: `ZK1`'s frustum turned about its axis by the kernel's rotation, its rim within rounding of the cap's) | 12 | new (from S9d.3b's) |
| `Degenerate`, S9d.3b.1's survey: a plane within rounding of a cylinder's direction (`ZK9`, `ZL1`: a frustum of radii 1 and 6 across the cylinder, turned a quarter turn about y by the kernel's rotation, its disc of radius 6 cutting the wall) | 8 | new (from S9d.3b's) |
| `Degenerate`, S9d.2's: a tangency between the inputs, a coaxial cylinder and sphere whose discriminant vanishes identically (a cylinder's rim on a sphere's equator, the wall tangent to it along the circle: `ZH5` to `ZI3` of the four `bop*_simple` grids) | 32 | 0 |
| `Degenerate`, S9d.2a's: a vertex of one input on the other's face (a unit box's corner at a unit sphere's centre, its far corners on the sphere) | 21 | 0 |
| `Degenerate`, S9d.1's: a tangency between the inputs (the unit box quarter-turned, its wall `x = -1` tangent to the unit sphere) | 12 | 0 |
| `Degenerate`, S9c.2b.1's: two cylinders' section within the resolution of a node (equal cylinders in a turned frame whose axes meet: `ZD9`, `ZE2`) | 8 | 0 |
| A stack with cylindrical or conical walls as a Boolean argument (S9c; `bcut_simple/G9`, `H3`: a coaxial frustum fused onto a cylinder, then a rod cut from the fuse, 2) | 6 | +2 (from S9d.3b's) |
| An arc ending off its circle on exact models (`bfuse_simple/E1`, S9c) | 1 | 0 |
| `Degenerate`, S9c.1's: a piece thinner than the resolution (a box's corner on the cylinder within rounding) | 40 | 0 |
| `Degenerate`, S9c.1's: two meetings within rounding along an arc (the same) | 16 | 0 |
| `Degenerate`, S9c.1's: a meeting at every seam tried (two stacked cylinders turned about their axis, their rims one circle) | 12 | 0 |
| `Degenerate`, S9c's: a tangency between the inputs (a box's wall tangent to the cylinder after a quarter turn 8; cylinders of radii 1 and 0.5 touching at a point, S9c.2a's, 16) | 24 | 0 |
| `Degenerate`, S9c.1's: solids touching at a vertex (a cut of crossed cylinders) | 4 | 0 |
| `Degenerate`, S9b.1's: a face thinner than the resolution | 99 | 0 |
| `Degenerate`, S9b.1's: a face using an edge both ways | 9 | 0 |
| `Degenerate`, S9b.1's: two solids touching at a point 1, a face touching itself at a vertex 1 | 2 | 0 |
| `Degenerate`: a result thinner than the resolution | 103 | 0 |
| `Degenerate`: a profile piece thinner than the resolution | 24 | 0 |
| `Degenerate`: a result touching itself at a point | 11 | 0 |
| `Degenerate`: a result touching itself along an edge | 6 | 0 |
| `Degenerate`: a result's hole outside its boundary 8, two results touching 4 | 12 | 0 |
| `atan2` in `dset`, which the adapter does not evaluate | 48 | 0 |
| `explode` without the native selector (a Boolean's result into solids 4, a `pcone` into faces 4, a face of a box 1) | 9 | 0 |

The previous table's cone against a prism with arcs (S9d.3b, 66) is gone.
The sentinels were checked again: each is refused as before but
`bopfuse_simple/ZJ4`, which evaluates and is re-purposed. One sentinel is
added: `bopfuse_simple/ZK9` (a plane within rounding of a cylinder's
direction); the plane through a frustum's virtual apex, the rims within
rounding and the stack with a cone's wall keep `bopfuse_simple/ZG4`'s,
`S3`'s and `bcut_simple/L3`'s reasons. The ledger does not change.

After this survey the adapter turns every solid by DRAW's quarter turns
exactly (a signed permutation, as it did prisms), and S9d.3b.2 landed (a
cone and a sphere in loops; carriers whose rulings reach the other's
asymptotic directions). The 1,802 cases were run again (2026-09-29).
Native DRAW's statuses are unchanged. Only cones against cylinders change,
and no case reaches S9d.3b.2's code. `ZK9` and `ZL1` of the four grids
(8) evaluate: the frustum's axis is `(1, 0, 0)` exactly now. Their
volumes, areas and centres of gravity are `cones_boolean_reference.py`'s
within 1.4e-15 relative (native DRAW's up to 4.5e-5 off, `bopcut_simple/ZK9`
128.121 against 128.1268, its mirror image `ZL1` 128.128); they take about
16 seconds each on the debug worker and are registered (`bopfuse_simple/ZK9`,
the sentinel for the rounded turn, re-purposed). `ZK2` to `ZK4` (12) are
refused as `ZK1` is, a tangency between the inputs (the frustum's base rim
the cap's rim), not as a piece thinner than the resolution: the rounding
is gone, the coincident rims are not. Four of the five timeouts evaluate
within the survey's 120 seconds; `bopcommon_simple/ZK8` still times out
there, eight cases running at once. `ZK7` and `ZK8` are faster with
the exact turn but not fast enough: alone on the debug worker they take
23 to 34 seconds at a load average of 5 (`bopfuse_simple/ZK8` 31.5 and
`boptuc_simple/ZK7` 29.9 against 90.0 and 60.2 on the worker before the
turn, side by side), and up to 120 with the machine loaded; on a release
build of the worker, which the runner does not choose, 1.7 to 4.9. They
stay unregistered. Rust evaluates 946 (934 before), 939
registered, and fails none; one times out; 632 cases native DRAW
evaluates are refused (640 before). By reason, against the table above:
the timeouts 1 (-4); `Degenerate`, a tangency between the inputs (`ZK1`
to `ZK4`) 16 (+12); a piece thinner than the resolution (`ZK2` to `ZK4`)
and a plane within rounding of a cylinder's direction (`ZK9`, `ZL1`) gone
(-12, -8); every other row unchanged (633 in all). `bopfuse_simple/ZK1`
(the coincident rims) replaces `bopfuse_simple/ZK9` as a sentinel; the
others are refused as before. The volume audit was run again on every
registered case of the group and the 16 of `ZK7` to `ZL1` (947, `vprops`
and `sprops` before each `checkprops`, both backends). Rust's volumes and
centres of gravity are the earlier audit's within 3.6e-15 on its 903 cases
and the S9d.3b.1 survey's on its 36 cones; the turned spheres (`ZI4` to
`ZJ3`, turned by the rounded rotation before and exactly now) agree with the closed forms within 7.8e-16 relative in volume,
6.8e-16 in area and 6.7e-16 in the centres (native DRAW within 2.5e-6),
the cones with theirs within 1.4e-15. Rust agrees with native DRAW to its
printed digits or 1e-6 relative in all but 35, native DRAW off in each:
`bopcut_simple/A9`'s centre and the crossed cylinders `ZE3` to `ZE6`, as
the earlier audit found (their areas now checked too: native up to 6.1e-6
relative off a quadrature, Rust's within 6.7e-16), and 20 of the crossing
frusta (`ZK5` to `ZL1`). A full contract run of the manifest (both backends, 30 seconds a case)
holds for every case but two, `bopfuse_simple/ZE6` and
`bopcut_simple/ZF2`, which timed out with the machine's load average
near 16 and pass run again (7 and 15 seconds at a load of 5, as on the
workers before S9d.3b.2 and the exact turn). The ledger does not change.

**S9d.4a's whole torus against polyhedral prisms in the Boolean group.**
The same 1,802 cases were run again on both backends after S9d.4a
(2026-09-29, `survey_upstream_tests.py --boolean`, the public dataset
read through `--data-dir`, 120 seconds a case; no case loads it).
Native DRAW's statuses are unchanged. The grids hold no torus against a
box: their tori are the 16 `pcylinder`s of radius 4 and height 8 with a
coaxial torus of radii 4 and 1 at half their height (`ZL2` to `ZL5` of
the four `bop*_simple` grids, turned about the axis by quarter turns),
refused until now as a solid with curved faces in any position (S9c) and
now as S9d.4b's `OutOfDomain`, a torus against a curved face, and
`bopfuse_simple/ZP6`'s three copies of one torus, which the adapter
refuses as solids other than prisms sharing ids. No case reaches a torus
segment or wedge. No case evaluates newly and none fails; the one
timeout of the survey before, `bopcommon_simple/ZK8`, evaluates within
the 120 seconds this time (the eight `ZK7` and `ZK8` cases stay
unregistered, as above). Rust evaluates 947, 939 registered; 632 cases
native DRAW evaluates are refused. Against the table above: solids
other than prisms 61 (-16, the adapter's shared ids only); a torus
against a curved face (S9d.4b, `OutOfDomain`) 16, new; the timeouts 0
(-1). The sentinel `bopfuse_simple/ZL2` is refused as S9d.4b's now, its
purpose updated; the others are refused as before. The `gdml_public`
grid's tori (`A1`, `A2`, `A9`, `B6`, outside the surveyed grids) are
refused by both backends' hosts. The ledger does not change.

**S9d.4b's tori and the adapter's curved primitives in the Boolean
group.** The same 1,802 cases were run again on both backends after
S9d.4b.1 (torus segments and wedges against polyhedral prisms), S9d.4b.2a
(a whole torus against prisms with arcs, spheres and cones) and S9d.4b.2b
(two whole tori), with one change to the adapter (2026-09-30,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`, 120 seconds a case, six at once; no case loads it). The
adapter built every `pcone`, `psphere` and `ptorus` under the unspecified
operation, so two of them in one Boolean shared ids, which it refused as a
solid other than a prism sharing ids. It now builds each under an
operation of its own, and a copy of one (`copy` and `tcopy` keep ids, and
so do `ttranslate` and `trotate`) that shares ids with the other argument
is built again from its constructor's numbers in its own frame, as a prism
is extruded again. Native DRAW's statuses are unchanged. 36 cases
evaluate newly; none fails and none times out:

* `ZL2` to `ZL5` of the four grids (16), S9d.4b.2a's: a `pcylinder` of
  radius 4 and height 8 and a coaxial torus of radii 4 and 1 at half its
  height (turned about the axis by quarter turns), the wall through the
  tube's centre circle, meeting the tube in its top and bottom circles
  (S9d.4a's survey found them S9d.4b's `OutOfDomain`; `bopfuse_simple/ZL2`,
  its sentinel, has evaluated since S9d.4b.2a).
* `ZM1`, `ZM2`, `ZM4`, `ZM5` and `ZM6` of the four grids (20), the
  adapter's: two `pcone`s, a frustum of radii 8 and 4 and height 8 and a
  coaxial one of radii 2 and 1, standing on its top disc, inside it with
  its top disc on the top disc, from 1 below the bottom disc to the top
  disc, inside clear of both discs, or through both (S9d.3b.1's cones;
  the common of `ZM1` and the reversed cuts of `ZM2` and `ZM5` empty).

The 36 were run again with `vprops` and `sprops` before each
`checkprops`. Rust's volumes, areas and centres of gravity agree with
closed forms (the solids of revolution at 40 digits: the torus's halves by
Pappus, the frusta's sections integrated piecewise and their faces as
frusta and annuli) within 2.6e-15 relative in volume, 3.1e-16 in area and
2.7e-15 in the centres; native DRAW's agree with them to its printed
digits (within 3.1e-6 relative), and Rust's with native DRAW's to its
printed digits or 1e-6 relative on every one. In the contract run below
the tori take 8.7 to 8.9 seconds each on the debug worker, the frusta
under 0.4 (the slowest registered cases, `ZK9` and `ZL1`, 12.7). The 35 not
yet in the manifest are registered, each `viewer_skipped` on both, and
`bopfuse_simple/ZL2`'s purpose is theirs now.

The other 41 of the 61 the adapter refused for shared ids reach the
kernel, which refuses them as `Degenerate`: a cone's apex on the other
input's surface 21 (`ZL6` to `ZL9`: a frustum of radii 4 and 2 and height
4 standing on the top disc, on the first's cone, turned about the axis by
quarter turns; `ZM3`: one of radii 2 and 1 and height 8 inside from disc
to disc; both with their virtual apex the first's, `(0, 0, 16)`; and
`boptuc_simple/ZN2`, below), a tangency between the inputs 19 (`ZM7` to
`ZN1` of the four grids: a frustum of radii 4 and 8 and height 4 under the
first, its top rim the first's bottom rim and its top disc the bottom
disc, as `ZK1`'s coincident rims; `ZN2` of the fuse, cut and common grids:
a frustum of radii 8 and 4 and height 4 moved 4 along x, the two bases on
one plane, whose reversed cut finds first that frustum's virtual apex
`(4, 0, 8)` on the first's top rim), and a torus tangent to the other
input's surface 1 (`bopfuse_simple/ZP6`: three copies of a torus of radii
100 and 20, two turned a quarter turn about x and about y, one centre and
perpendicular axes, the tubes of each two touching at four points of the
line their equatorial planes share). Rust evaluates 983 (947 before), 975
registered (940: 939 and `bopfuse_simple/ZL2` since S9d.4b.2a), and fails
none; none times out; 596 cases native DRAW evaluates are refused (632
before). Against the table above: solids other than prisms sharing ids 0
(-61); a torus against a curved face (S9d.4b, `OutOfDomain`) 0 (-16);
`Degenerate`, a cone's apex on the other input's surface 21, new; a
tangency between the inputs, two frusta's rims or bases (`ZM7` to `ZN2`)
19, new; a torus tangent to the other input's surface 1, new; every other
row unchanged. No case reaches S9d.4b.1's segments and wedges, a sphere
against a torus or two tori but `ZP6`. The sentinel
`bopcommon_simple/ZL6` (two cones sharing ids) is refused by the kernel
now, as a cone's apex on the other input's surface, its purpose updated;
the others are refused as before. `ZP6` is no sentinel: native DRAW's
area is a known failure on Linux (its `TODO`). The volume audit was run
again on the 947 cases of the last one (every registered evaluating case
then and `ZK7` to `ZL1`, both backends): Rust's volumes, areas and
centres of gravity are that audit's bit for bit, so neither S9d.4 nor the
adapter's ids changed a registered result (one run lost its output with
the machine loaded and was identical run again), and its disagreements
with native DRAW are the same 35. The `gdml_public` grid's tori (`A1`,
`A2`, `A9`, `B6`) stay refused by both hosts: the Rust adapter at their
first line, `compound result`, and after it at `ptorus name plane R r
angle` (a torus placed on a DRAW `plane`, which the adapter takes for a
`psphere` only); the native host at `add` or `wire`. A full contract run
of the manifest (both backends, 30 seconds a case) holds for every case.
The ledger does not change.

**S9d.2c, S9d.3c and S9d.4c in the Boolean group.** The same 1,802 cases
were run again on both backends after S9d.2c (a turned cap's circles of
unequal axes against a cylinder, sphere-cylinder loops in turned frames),
S9d.3c (cone-cylinder and cone-cone loops in any frames), S9d.4c (a
sphere's cap or zone circle against a torus, torus v-segments and wedges
against prisms with arcs, spheres, cones and tori) and the certified
integrals' speed-up between them (2026-09-30,
`survey_upstream_tests.py --boolean`, the public dataset read through
`--data-dir`, 120 seconds a case, six at once; no case loads it). Native
DRAW's statuses are unchanged, and so are Rust's: no case evaluates newly,
none changes its refusal, none fails and none times out. The group's
`psphere`s and `ptorus`es are whole solids, so no case has a cap, a zone
or a torus part, and its cones meet a cylinder in circles or in rings
about the frustum (S9d.3b.1's), no case in loops. Rust evaluates 983,
596 cases native DRAW evaluates are refused, each for the reason the
survey above found, and the sentinels are refused as before.

What changes is the time. `ZK7` and `ZK8` of the four grids (8), a
frustum of radii 6 and 1 and height 10 across the `pcylinder`, its axis
crossing the cylinder's at right angles at half its height and its wide
end through both caps, were right since S9d.3b.1's survey but took 23 to
34 seconds alone on the debug worker at a load average of 5 (up to 120
with the machine loaded), too near or past the contract's 30. The
certified integrals' speed-up halves them: 8.7 to 12.9 seconds alone at a
load average of 5 to 7 (`bopfuse_simple/ZK8` 8.9 at a load of 1.7,
against 18.9 on the worker before the speed-up, side by side), 8.2 to 9.4
in the contract run below. Their volumes, areas and
centres of gravity are `cones_boolean_reference.py`'s within 1.3e-15
relative in volume, 3.9e-16 in area and 9.7e-16 in the centres (native
DRAW's within 7.8e-6, 7.2e-6 and 4.3e-6, BRepGProp on approximated
quartics, as S9d.3b.1's survey found). They are registered, each
`viewer_skipped` on both: 983 registered (975 before), every case the
group evaluates.

The volume audit was run again on the 983 cases (`vprops` and `sprops`
before each `checkprops`, both backends). Native DRAW's values are the
last audit's to the digit. Rust's are bit for bit on 854 and differ within
rounding on 129, at most 2.2e-15 relative in volume, 7.2e-16 in area and
1.4e-15 in the centres (scaled by the cube root of the volume): 118 of
the cases with curved faces (a cylinder and a turned box, crossed
cylinders, a box and a frustum across its walls, spheres, frusta across
the cylinder, the torus: `bfuse_complex/J5`, 45 of `Y5` to `ZE6`, 11 of
`ZH1` to `ZH4`, `ZI8` to `ZJ3`, `ZK5` to `ZL5`, `bopcommon_simple/ZP9`) by
the certified integrals' speed-up (its binary64 series products summed
under one error bound; the worker built at its parent reproduces the last
audit bit for bit, the worker at it these 118), and the centres of 11 of
`ZI4` to `ZI7` (the sphere turned so the cap's plane holds its axis, its
section through the poles) by up to 2.6e-16 with S9d.4c's kernel (the
worker before it gives the last audit's), their volumes and areas
unchanged. Against the closed forms and the references: the turned
spheres (`ZI4` to `ZJ3`) within 1.4e-15 relative in volume, 1.2e-15 in
area and 9.2e-16 in the centres (7.8e-16, 6.8e-16 and 6.7e-16 before),
the crossed cylinders' areas (`ZE3` to `ZE6`) within 1.0e-15 (6.7e-16),
the crossing frusta (`ZK5` to `ZL1`) within 1.7e-15, 3.9e-16 and 2.8e-15,
the torus and the coaxial frusta (`ZL2` to `ZL5`, `ZM1`, `ZM2`, `ZM4` to
`ZM6`) within
1.3e-15, 5.5e-16 and 2.2e-15. Rust agrees with native DRAW to its printed
digits or 1e-6 relative on all but the same 35, native DRAW off in each.

`bopfuse_simple/ZP6` (three copies of a torus of radii 100 and 20 about
perpendicular axes, the tubes of each two tangent at four points of the
line their equatorial planes share) is refused as before, a torus tangent
to the other input's surface. The `gdml_public` grid's tori (`A1`, `A2`,
`A9`, `B6`) stay refused by both hosts, as in the survey above: the Rust
adapter at `compound result` and, after it, at a `ptorus` on a DRAW
`plane`; the native host at `add` or `wire`. A full contract run of the
manifest (both backends, 30 seconds a case) holds for every case, the
slowest Boolean case 9.4 seconds (`boptuc_simple/ZK8`; the tori `ZL2` to
`ZL5` 2.4 to 2.5, 8.7 to 8.9 before the speed-up). The ledger does not
change.

**S9e.2's given results in the Boolean group.** The boolean group's 25
cases in which a Boolean's result is an argument of another Boolean
command and which were not registered as evaluating (the scan of S9e.1's
survey: `L3` the sentinel among them) were run on both backends after S9e.2
(2026-09-30, `survey_upstream_tests.py --case`, the public dataset read
through `--data-dir`, 120 seconds a case). `bcut_simple/L3` to `L6`, DRAW's
rollex (a disc of radius 60 less a pocket of radius 40 across its rim over
the top 6 of its height 20, a stack, then cut by a cylinder of radius 30
standing on the pocket's floor, its profile on a plane facing down in `L3`
and `L4`), now evaluate on both, each `viewer_skipped`: Rust's `checkprops
-s 30153` reads 30152.9544786437837, `generate_given_boolean_fixtures.py`'s
`rollex_turned_cut` and `rollex_flat_cut` 30152.95447864378141 (the
reference's within 7.6e-17 relative; native DRAW's 30153 to its printed digits).
A volume audit (`vprops` and `sprops` before the `checkprops`, both
backends): Rust's volume 199635.184534622938 in all four, the reference's
199635.1845346228743 (3.2e-16 relative), its centre within 1.5e-16 of the
reference's (scaled by the cube root of the volume); native DRAW's
199635 and 30153 to its printed digits, the capture's
`BRepGProp` values 199635.18453462288 and 30152.954478643784. They are
registered (`L3` no longer a sentinel): 987 registered. The others do not
change: `G9` and `H3` (a cone fused on a cylinder, then cut) are refused
as S9e.3's (a result of solids other than prisms); 6 of the public dataset
give a Boolean a restored shape (not a solid the adapter made), 6 stop at
constructs the adapter does not read (`mkplane`, a `prism` of a restored
face, a restored 2D ellipse and trimmed curve), 5 load private data on
both; `bopcut_simple/ZQ1` (`wire`, unsupported on both) and
`bopfuse_simple/ZP6` (a torus tangent to the other input) as before. None
fails or times out. A contract run of the four (both backends) holds; the
ledger does not change.

**S9e.1 and S9e.2 in the Boolean group.** The same 1,802 cases were run
again on both backends after S9e.1 (a Boolean's result of prisms with
lines, arcs and circles, the only solid of its Boolean, given to another
Boolean on its construction's exact model) and S9e.2 (a stack with arc
walls, an S9b.1 result given with arcs and one solid of a result of
several, on their construction's curved arrangement run again)
(2026-09-30, `survey_upstream_tests.py --boolean`, the public dataset
read through `--data-dir`, 120 seconds a case, four at once; no case
loads it). Native DRAW's statuses are unchanged. Rust's change in the six
cases S9e.2's run of the chained cases (above) found, and in no other:
`bcut_simple/L3` to `L6`, the rollex, evaluate, each `viewer_skipped`, as
registered there; `bcut_simple/G9` and `H3` (a frustum fused onto a
cylinder, then a rod cut from the fuse) are refused as S9e.3's ("a
Boolean's result of solids other than prisms given to another Boolean"),
where the survey before refused them as a solid with curved faces or edges
in any position (S9c). No case evaluates newly besides the rollex, none
fails and none times out. The group's other cases giving a result to a
second Boolean are S9b.2's pockets (registered since S9b.2's survey) or
the 19 others of the 25 S9e.2's run found unchanged (restored arguments,
constructs the adapter does not read, private data, `ZQ1`, `ZP6`). Rust
evaluates 987 (983 before), every one registered; 592 cases native DRAW
evaluates are refused (596), each for the reason the surveys above found,
and 223 are unsupported on both. The previous table's stack with
cylindrical or conical walls as a Boolean argument (S9c, 6) is gone: 4
evaluate, 2 are S9e.3's. The sentinels are refused as before (`L3` no
longer among them).

The volume audit was run again on the 987 cases (`vprops` and `sprops`
before each `checkprops`, both backends). On the 983 of the last audit
both backends' values are that audit's bit for bit (Rust) and to the
digit (native DRAW): S9e.1's and S9e.2's kernels, their amendments
included (a reversed piece joined on one surface, a moved stack's height
range, parallel cylinders met from the second input's arc's side), change
no registered case's values, so no worker had to be built at an earlier
commit to explain a difference. On the rollex Rust's values are S9e.2's
audit's bit for bit: volume 199635.184534622938 and area
30152.9544786437837 in all four, the centre (-1.96306181728239970,
63.5144823959729763, 9.06871924355272263) in `L3` and `L4` and
(-1.96306181728239615, 63.5144823959729834, 9.06871924355272085) in `L5`
and `L6`. Against
`generate_given_boolean_fixtures.py`'s `rollex_turned_cut` and
`rollex_flat_cut` they are within 3.2e-16 relative in volume, 7.5e-17 in
area and 1.5e-16 in the centre (scaled by the cube root of the volume;
8.4e-17 in `L5` and `L6`); native DRAW's printed 199635, 30153 and
(-1.96306, 63.5145, 9.06872) within 9.2e-7, 1.5e-6 and 3.0e-7, its
printed digits. Rust agrees with native DRAW to its printed digits or 1e-6
relative on all but the same 35 as before, native DRAW off in each.

`bopfuse_simple/ZP6` is refused as before, a torus tangent to the other
input's surface. The `gdml_public` grid's tori (`A1`, `A2`, `A9`, `B6`)
stay refused by both hosts: the Rust adapter at `compound result` and,
after it, at a `ptorus` on a DRAW `plane`; the native host at `add` or
`wire`. A full contract run of the manifest (both backends, 30 seconds a
case) holds for every case, the slowest Boolean case 9.8 seconds
(`boptuc_simple/ZK8`), the rollex 1.5 to 1.6. The ledger does not change.

**S9e.3a's given results in the Boolean group.** The boolean group's 20
unregistered cases in which a Boolean's result is an argument of another
Boolean command (S9e.2's 25 less the rollex's four, the scan finding one
fewer of the dataset's restored arguments) were run on both backends after
S9e.3a (2026-10-02, `survey_upstream_tests.py --case`, the public dataset
read through `--data-dir`, 120 seconds a case). `bcut_simple/G9` and `H3`
(a frustum of radii 7 and 6 fused onto a cylinder of radius 9, then a rod
of radius 1 about `(5, 0)` cut from the fuse) are refused now as "a
tangency between the inputs (S9c)": the rod's circle touches the
frustum's top circle of radius 6 at `(6, 0, 4)` (where both seams lie, as
DRAW's comment says), the decisions' refusal; before, the given result was
refused as S9e.3's. Native DRAW's area 727.481 is the reference's
(`generate_given_curved_boolean_fixtures.py`'s `g9_cut`, 727.4813665;
`g9_clear_cut`, the rod moved to `x = 4.5`, evaluates on the kernel with
the same measures). The others do not change: 5 load private data on both,
6 give a Boolean a restored shape (not a solid the adapter made), 5 stop at
constructs the adapter does not read (`mkplane`, a `prism` of a restored
face, a restored 2D ellipse and trimmed curve), `bopcut_simple/ZQ1`
(`wire`, unsupported on both) and `bopfuse_simple/ZP6` (a torus tangent to
the other input). None evaluates newly, fails or times out; none is
registered; the ledger does not change.

**S9e.3a and S9f.1 in the Boolean group.** The same 1,802 cases were run
again on both backends after S9e.3a (a Boolean's result with spheres,
cones and tori given to another Boolean, chains of up to three Booleans,
`curved/chain.rs`), S9f.1 (spline prisms against prisms of lines in any
position, `curved/spline_walls.rs`) and `b2765f20` (a plane crossing a
sphere within the resolution of tangency `Degenerate`, also two spheres'
radical plane) (2026-10-03, `survey_upstream_tests.py --boolean`, the
public dataset read through `--data-dir`, 120 seconds a case, four at
once; no case loads it). No status changes on either backend: Rust
evaluates 987, every one registered, 592 cases native DRAW evaluates are
refused, 223 are unsupported on both. One refusal's reason changes, as
S9e.3a's run of the chained cases (above) found: `bcut_simple/G9` and `H3`
are refused as a tangency between the inputs, where the survey before
refused them as S9e.3's. No other case evaluates newly or changes its
refusal, none fails and none times out, and the sentinels are refused as
before, each for its reason. S9f.1's kernel reaches no case of the group:
its spline solids come from `nurbsconvert` (96 cases), which the adapter
does not run, and its other refusals are unchanged; S9e.3a's reaches only
`G9` and `H3`. `b2765f20` turns no case `Degenerate`: the group's spheres
are cut by planes through their centres (`ZI4` to `ZJ3`) or by
`bopcommon_simple/ZP9`'s wall 1.053 from the centre of a sphere of radius
7.5, the box's other planes missing it by 9.5 or more, and no refused
case's reason becomes the new one.

The volume audit was run again on the 987 cases (`vprops` and `sprops`
before each `checkprops`, both backends): both backends' values are the
last audit's bit for bit (Rust) and to the digit (native DRAW), the
rollex's included, so S9e.3a's and S9f.1's kernels, their amendments (R4's
knot removal before every lift, a plane's section of a cone normal to its
axis given as a circle) and `b2765f20` change no registered case's values,
and no worker had to be built at an earlier commit. The references'
agreement is the last audit's (the rollex within 3.2e-16 relative of
`generate_given_boolean_fixtures.py`'s). Rust agrees with native DRAW to
its printed digits or 1e-6 relative on all but the same 35 as before,
native DRAW off in each. (Three cases first ran with an empty script, the
audit's own two backend threads rewriting one copy of it; run again, they
agree.)

`bopfuse_simple/ZP6` is refused as before, a torus tangent to the other
input's surface. The `gdml_public` grid's tori (`A1`, `A2`, `A9`, `B6`)
stay refused by both hosts as before: the Rust adapter at `compound
result` and, after it, at a `ptorus` on a DRAW `plane`; the native host at
`add` or `wire`. A full contract run of the manifest (both backends, 30
seconds a case) holds for every case, the slowest Boolean case
13.5 seconds (`bopcut_simple/ZK8`, on a machine loaded by other work:
alone 13.0, and 11.3 on a worker built at the last survey's commit at the
same load; 9.8 in that survey), the rollex 1.7 to 1.8. The ledger does not change.

**S9e.3b's given meetings in the Boolean group.** The same 20 unregistered
cases giving a Boolean's result to another Boolean were run on both
backends after S9e.3b's kernel (a given result's meeting of two curved
faces, or a cone's or torus's general section, met by the partner,
`curved/triple.rs`) (2026-10-03, `survey_upstream_tests.py --case`, the
public dataset read through `--data-dir`, 120 seconds a case, four at
once). Every status and reason is S9e.3a's: 5 load private data on both, 6
give a Boolean a restored shape, 5 stop at constructs the adapter does not
read, `bopcut_simple/ZQ1` stops at `wire` on both, `bcut_simple/G9` and
`H3` are refused as a tangency between the inputs and `bopfuse_simple/ZP6`
as a torus tangent to the other input. S9e.3b's kernel reaches none of
them (none of the given results the adapter makes holds a meeting of two
curved faces met by the partner); none evaluates newly, fails or times
out, so none is registered, no volume audit is due and the ledger does not
change.

**S9e.3b, S9f.2a and S9f.2b.1 in the Boolean group.** The same 1,802
cases were run again on both backends after S9e.3b (a given result's
meeting of two curved faces, or a cone's or torus's general section, met
by the partner, `curved/triple.rs`; with it the validator decides a hole
in a torus band wound in v by the +u ray's signed crossings, a result's
bounds take its torus, cylinder and cone faces' bulges, and a vertex at a
sphere face's pole stays), S9f.2a (spline walls against arc, circle and
spline walls of a prism whose axis is exactly parallel,
`curved/spline_parallel.rs`) and S9f.2b.1 (spline walls against cylinder
walls on crossing axes, `curved/spline_crossing.rs`, `Curve3::WallMeet`)
(2026-10-03, `survey_upstream_tests.py --boolean`, the public dataset read
through `--data-dir`, 120 seconds a case, four at once; no case loads it).
Every status and every refusal's reason is the last survey's, on both
backends: Rust evaluates 987, every one registered, 592 cases native DRAW
evaluates are refused, 223 are unsupported on both. No case evaluates
newly, none fails and none times out, and the sentinels are refused as
before, each for its reason. The three kernels reach no case of the group:
its spline solids come from `nurbsconvert` (96 cases), which the adapter
does not run, so no spline wall meets another prism's wall; and of the
results it gives to another Boolean only `bcut_simple/G9` and `H3` reach
the kernel, refused as a tangency between the inputs, and no given result
the adapter makes holds a meeting of two curved faces met by the partner
(S9e.3b's run of the chained cases, above).

The volume audit was run again on the 987 cases (`vprops` and `sprops`
before each `checkprops`, both backends): both backends' values are the
last audit's bit for bit (Rust) and to the digit (native DRAW), among them
the 16 cylinders and coaxial tori (`ZL2` to `ZL5`), the 16 spheres
sectioned through their poles (`ZI4` to `ZI7`) and every case with
cylinder or cone faces. So S9e.3b's validator for torus bands, its bounds
of curved faces and its pole vertices, S9f.2a's mass integrals on
cylinders and cones taking a steep line by its box, and S9f.2a's and
S9f.2b.1's kernels change no registered case's values, and no worker had
to be built at an earlier commit. The references' agreement is the last audit's.
Rust agrees with native DRAW to its printed digits or 1e-6 relative on all
but the same 35 as before, native DRAW off in each.

`bopfuse_simple/ZP6` is refused as before, a torus tangent to the other
input's surface. The `gdml_public` grid's tori (`A1`, `A2`, `A9`, `B6`)
stay refused by both hosts as before: the Rust adapter at `compound
result` and, after it, at a `ptorus` on a DRAW `plane`; the native host at
`add` or `wire`. A full contract run of the manifest (both backends, 30
seconds a case) holds for every case, the slowest Boolean case 12.8
seconds (`boptuc_simple/ZK8`, on a machine at load 5 to 8 with a fuzz
campaign running; 13.5 in the last survey), the rollex 2.1. The ledger
does not change.

**S9e.4a's imported solids: cases restoring a shape for a Boolean.**
Before S9e.4a the adapter refused every restored shape as a Boolean
argument ("an argument other than a solid the adapter made"), and the
Boolean surveys above ran only the boolean group's self-contained cases.
Every case of every group that restores a shape and gives it to a Boolean
command (`bfuse`, `bcut`, `bcommon`, `btuc`, `bop`), 1,814 cases (1,371 of
the boolean group, 420 of `bugs`, 18 of `perf`, 4 of `lowalgos`, 1 of
`blend`), none
registered, was run on the Rust adapter after S9e.4a (2026-10-03, the public
dataset, 120 seconds a case, four at once), its restored solids imported
(`Solid::imported_with`) and the cases they reach run on native DRAW too.
Of the 1,814: 561 load private data; 598 give a Boolean a restored shape
that is not one solid (faces, shells, compounds); 123 stop at constructs
the reader does not represent, 276 at the validator's rejection of a
restored shape (an edge not C1, a pcurve off its edge: failures before
any Boolean, as the restore-only surveys found), 8 time out at their
first restores (before any Boolean), 77 at other commands the adapter does
not read (`bsection`, `halfspace`, `explode`, ...). The import reaches 171:
84 bodies none of the kernel's constructions (S9e.4b), 31 prisms whose arcs
round off their circles in their caps' frames (S9e.4b), 3 with spline
faces (S9f), 37 refused by S9's rules on their constructions (tangencies,
faces within the resolution of one plane, thin faces, a plane through a
cone's apex, a vertex on a face), and 16 evaluate on both backends: boxes,
hexahedra, wedges, prisms with cylindrical walls and cylinders of the
`CTO9xx` and `cts`/`pro` series fused, cut and intersected
(`bcommon_complex/C5`, `bcut_complex/J2`, `J7`, `O9`, `bfuse_complex/F5`,
`H4`, `I4`, `I8`, `J4`, `J7`, `L1`, `L3`, `N3` (two fuses, the first
result of imported solids given to the second), `Q2`,
`bopcommon_complex/K5`) and `bugs/modalg_6/bug21427` (a prism of nine
planes less a box), every check passing on both (`checkprops`,
`checkshape`, `checknbshapes` in `F5` and `Q2`; `viewer_skipped`). A
volume audit (`vprops` and `sprops` before each `checkprops`, both
backends): Rust's volume, area and centre agree with native DRAW's to its
printed digits in all 16 (`L1`'s fused cubes 4,000,000 and 160,000 within
3e-16 relative). They are registered (`data`, `viewer_skipped` on both;
1,089 cases): the contract holds for them on both backends within 30
seconds (Rust 0.2 to 5.6 s on a loaded host), and the ledger records
`F5`'s and `Q2`'s `checknbshapes` confirmed natively: mapped-and-verified 2
to 4, lost 12,844 to 12,842 (`--write-ledger`). One case,
`bcut_complex/I6`, failed when its restored profile's rounded points touched
in the construction: such a body is now S9e.4b's (unsupported).

**S9f.2b.2 and S9e.4b.1: the Boolean group and the restore cases.** Both
sets were run again after S9f.2b.2 (a spline wall's meeting with a
crossing cylinder turning back inside both faces, and a cap circle on a
wall along its axis, `curved/spline_crossing.rs`) and S9e.4b.1 (imported
prisms whose arcs' ends round off their circles in their caps' frames,
each end taken onto its circle, `curved/snapped.rs`) (2026-10-03, at
`93e6fcd0`, the public dataset, 120 seconds a case, four at once). The
1,802 self-contained cases of the Boolean group on both backends: every
status and every refusal's reason is the last survey's, 987 evaluating,
every one registered, 592 refused, 223 unsupported on both; none fails or
times out, the sentinels, `bopfuse_simple/ZP6` and the `gdml_public` tori
are refused as before. The 1,814 restore cases on the Rust adapter, and the
171 the import reaches on native DRAW too: only the 31 prisms whose arcs
round off their circles move. 7 evaluate on both backends with every check
(`bcut_complex/H3`, `K8`, `bfuse_complex/C9`, `E9`, `I6`, `N1`, `N9`:
prisms with one or two cylindrical walls of the `CTO9xx`, `cts` and `pro`
series fused with or cut by cylinders, boxes and prisms of planes); 18 are
refused by S9's rules (12 two faces within the resolution of one plane, 3
an edge of one input meeting an edge of the other, `bcut_complex/J4` and
`bfuse_complex/I9` a tangency, `bugs/modalg_7/bug29807_b1` two cylinders'
section within the resolution of a node) and 6 as S9e.4b.4's (arcs of two
circles meeting at a joint: `bcut_complex/E8`, `P4`, `bfuse_complex/D5`,
`E1`, `bugs/modalg_2/bug4993_1`, `_2`). `bcut_complex/I6` is refused as
S9e.4b's, as the last survey's amendment made it. Of the 171 the import
reaches, 23 now evaluate on both backends, 55 are refused by S9's rules, 84
are bodies none of the kernel's constructions, 6 are S9e.4b.4's and 3 have
spline faces (S9f); native DRAW does not evaluate 24 of them
(`checksection`, `bopargcheck`). Every other status and reason is the last
survey's, the same 8 timing out at their first restores. No case reaches
S9f.2b.2's loops and towers: the group's spline solids come from
`nurbsconvert`, which the adapter does not run, and the three restored
spline bodies (`bcut_complex/L9`, `O1`, `bfuse_complex/N7`) are refused
before any wall meets a cylinder.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends): the 7 new cases' volume, area and centre agree with native
DRAW's to its printed digits; the 1,003 registered cases' values are the
last audits' bit for bit on both backends (the 987 of the Boolean group's
survey at `b0b9adc6`, the 16 of S9e.4a's), so neither kernel moves a
registered value; the same 35 disagreements as before, native off in each.
The 7 are registered (`data`, `viewer_skipped` on both; 1,096 cases): the
contract holds for them on both backends within 30 seconds (Rust 1.8 to 6.5
s on a host at load 25 to 31). The ledger does not change. A full
contract run of the manifest (both backends, the dataset present, 30
seconds a case) holds for every case, the slowest Boolean case 16.2
seconds (`boptuc_simple/ZK8`, on a host at load 7 to 31; 12.8 in the last
survey), the restore cases 0.2 to 5.0, the rollex 4.4.

**S9e.4b.2, S9f.3a and S9f.3b: the Boolean group and the restore cases.**
Both sets were run again after S9e.4b.2 (imported polyhedra other than
prisms decided on their stored vertices, `polyhedra/imported.rs`), S9f.3a
and S9f.3b (spline prisms against spheres, caps and zones, and against
cones and frustums, in any position, `curved/spline_sphere.rs`,
`curved/spline_cone.rs`) and the near-parallel guards (two cylinders, a
spline wall and a cylinder, two tori within rounding of parallel or of one
surface, `Degenerate`) (2026-10-04, at `12b6c176`, the public dataset, 120
seconds a case, four at once). The 1,802 self-contained cases of the
Boolean group on both backends: every status and every refusal's reason
is the last survey's (`93e6fcd0`) field for field, 987 evaluating, every
one registered, 592 refused, 223 unsupported on both; none fails or times
out, the sentinels, `bopfuse_simple/ZP6` and the `gdml_public` tori are
refused as before. The 1,814 restore cases on the Rust adapter, and the
171 the import reaches on native DRAW too: only S9e.4b.2's 7 polyhedra
move, as its trial found. 4 evaluate on both backends with every check
(`bugs/modalg_1/buc60803`, `bug102_1` and `bug102_2`, two frustums of a
pyramid, the second on the first's top, fused in either order;
`bopfuse_complex/K5`, a frustum on a box's top); `bugs/modalg_2/bug578_1`
and `_2` are refused as the next Boolean's argument (frustums whose bases
lie 6.6e-7 to 2.0e-6 apart, above the resolution: their fuse is two
solids) and `bfuse_complex/D9` by S9's rules (a corner the two files share
stored 1e-13 apart: a face thinner than the resolution). Of the 171 the
import reaches, 27 now evaluate on both backends, 56 are refused by S9's
rules, 77 are bodies none of the kernel's constructions, 6 are S9e.4b.4's,
3 have spline faces (S9f) and 2 give the next Boolean several solids;
native DRAW does not evaluate 24 of them (`checksection`, `bopargcheck`).
Every other status and reason is the last survey's, the same 8 timing out
at their first restores. No case reaches S9f.3a's or S9f.3b's spline
walls: the group's spline solids are the 96 cases converting boxes by
`nurbsconvert` (no sphere or cone among them), which neither host runs,
and the three restored spline bodies (`bcut_complex/L9`, `O1`,
`bfuse_complex/N7`) are refused before any wall meets a partner. No
refusal of the near-parallel guards appears in either set.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends): the 4 new cases' volume, area and centre agree with native
DRAW's to its printed digits; the 1,010 registered cases' values are the
last audit's bit for bit on both backends (the 1,003 and the 7 of the
survey at `93e6fcd0`), so neither S9e.4b.2, S9f.3a, S9f.3b nor the
near-parallel guards move a registered value; the same 35 disagreements
as before, native off in each. The 4 are registered (`data`,
`viewer_skipped` on both; 1,100 cases): the contract holds for them on
both backends within 30 seconds (Rust 0.4 to 1.1 s). The ledger does not
change. A full contract run of the manifest (both backends, the dataset
present, 30 seconds a case) holds for every case, the slowest Boolean case
16.9 seconds (`bopcommon_simple/ZK8`, on a host at load 5 to 10;
`boptuc_simple/ZK8` 16.2 in the last survey, 13.3 here), the restore cases
0.1 to 3.9, the rollex 2.1 to 3.0.

**S9e.4b.3a, S9e.4b.3b and the fuzz fixes: the Boolean group and the
restore cases.** Both sets were run again after S9e.4b.3a (imported plane
pieces, a body of one sphere, cylinder or cone face and plane faces as its
primitive common its planes' half-spaces, `curved/pieces.rs`), S9e.4b.3b
(the kernel's own split pieces against curved faces, `curved/splits.rs`)
and the fixes and speedups that landed with them: the loops' certified
integrals, the validator's curve points guarded, the degree-eight
arithmetic (and a sheared projection for points no ruling separates), the
slowest inputs' speedup, prism walls' frames taking their axes bit for bit
(`Frame3::at`), a loop through one pole lifted at it, the three cone-pair
fixes (`ConeSec::lies_on`, coaxial rings nested, a ring on a floor
refused), the narrower `Meet` and `Toric` integrals, a torus cap's round
ends and a plane within rounding of a cone's apex refused (2026-10-05, at
`d665df29`, the public dataset, 120 seconds a case, four at once). The
1,802 self-contained cases of the Boolean group on both backends: every
status, reason and error is the last survey's (`12b6c176`) field for
field, 987 evaluating, every one registered, 592 refused, 223 unsupported
on both; none fails or times out, the sentinels, `bopfuse_simple/ZP6` and
the `gdml_public` tori are refused as before, and the new guards refuse no
case. The 1,814 restore cases on the Rust adapter, and the 171 the import
reaches on native DRAW too: no status moves; 41 refusals' reasons move from
"an imported solid other than a prism, a sphere, a cone or a torus", native
DRAW evaluating all 41. S9e.4b.3a's 39 are refused as its trial found: 18
(every case restoring `so1` or `so4`, whose rims OCCT split at the seam) as
a piece other than its primitive common its planes' half-spaces and 20
(pairs of `so2`, `so3`, `so5`, `so6` and `so7`) as faces of both inputs on
one sphere, both S9e.4b.3c's, and `bugs/modalg_1/buc60926` as a plane
through a cone's apex (`shading_132`'s planes through its virtual apex).
Two CTO parts beyond the trial reach the pieces' recognition now:
`bcut_complex/G4`, whose part is a box with a cylindrical boss (refused as
a piece other than its primitive common its planes, S9e.4b.3c's), and
`bcut_complex/I6`, whose tool is a block less a half cylinder tangent to
the block's own faces `y = 0` and `y = 70`: `Degenerate("a tangency between
the inputs (S9c)")`, raised in the piece's own first arrangement (its
cylinder primitive against its hull), so the reason names the wrong pair
(the tool alone against a far box is refused alike). Of the 171, 27
evaluate on both backends, 58 are refused by S9's rules, 39 as S9e.4b.3c's,
36 are bodies none of the kernel's constructions, 6 are S9e.4b.4's, 3 have
spline faces (S9f) and 2 give the next Boolean several solids; native DRAW
does not evaluate 24 of them. The same 8 time out at their first restores.
S9e.4b.3b's split pieces move no case.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,014 registered Boolean cases): native DRAW's
values are the last audit's bit for bit; the Rust adapter's statuses are
the same and 115 cases' values moved within rounding (49 volumes, 34 areas
and 250 centre coordinates; volumes and areas at most 6.6e-16 relative, a
centre at most 1.6e-15 of the solid's size): 113 at the loops' certified
integrals (`cfeab65d`) and 60 at the narrower `Meet` integrals
(`f87d150d`: two cylinders, or a cylinder and a cone, in `bop*_simple`'s
`ZE3` to `ZE6`, `ZI8` to `ZJ3` and `ZK5` to `ZL1`), 58 at both. None
crosses native DRAW's printed digits: 81 agree to them as before, 34 are
among the same 35 disagreements as before, each flagging the same fields,
native off in each. No case newly evaluates, so none is registered (1,100
cases) and the ledger does not change. A full contract run of the manifest
(both backends, the dataset present, 30 seconds a case) holds for every
case, the slowest Boolean case 6.0 seconds (`bopfuse_simple/ZK8`, on a host
at load 5 to 10; `bopcommon_simple/ZK8` 16.9 in the last survey), the
restore cases 0.1 to 2.6, the rollex 0.4 to 0.7.

**S9e.4b.3c.1 and the switches' speed-up: the Boolean group and the
restore cases.** Both sets were run again after S9e.4b.3c.1 (two inputs on
one sphere in general position, faces on one sphere as faces on one
surface, and sphere pieces whose rims OCCT divided, split at their stored
vertices; a groove refused on import), the speed-up of the fuzz target's
`TORUS_PAIRS` and `TURNED_PARTS` (memos of a circle's resultants with a
torus and of isolated roots, the root refinement, the fields' signs, a
torus meeting's jets shared by its two pcurves, a containment ray
undecided at once beside a fin no boundary test decides) and a prism's
caps keeping its frame bit for bit (`59d0c57b`) (2026-10-05, at
`f4b584f7`, the public dataset, 120 seconds a case, four at once). The
1,802 self-contained cases of the Boolean group on both backends: every
status, reason and error is the last survey's (`d665df29`) field for
field, 987 evaluating, every one registered, 592 refused, 223 unsupported
on both; none fails or times out. The 1,814 restore cases on the Rust
adapter, and the 171 the import reaches on native DRAW too: native DRAW's
statuses and reasons are the last survey's, and on the Rust adapter only
the 38 cases restoring two of `so1` to `so7` (pieces of one sphere of
radius 10) and `bcut_complex/I6` move, as S9e.4b.3c.1's trial found. 14
evaluate on both backends with every check: the hemisphere and the cap
above `z = 5` (`so1` and `so4`, both rims divided at the stored seam)
in `bcommon_complex/B2`, `bcut_complex/C5` and `C6` (the cap less the
hemisphere, empty), `bfuse_complex/B5`, `bugs/modalg_2/bug413_1` and `_2`
and `bugs/moddata_1/bug183_2` and `_3`; the cap and a wedge (`so4` and
`so2`) in `bcommon_complex/B5`, `bcut_complex/D2` and `D3`,
`bfuse_complex/B8` and `bug183_5`; and the common of a wedge and its
turned copy (`so2` and `so6`), `bcommon_complex/B9`. 13 are refused as
S9e.4b.3c.2's (a vertex or a circle of both inputs on one sphere: `so1`
and `so2`, `so2` and `so3`, `so5` and `so2`), 11 with `so6` or `so7` by
S9's rules (two faces within the resolution of one plane, a piece or a
result thinner than the resolution: their corners within rounding of the
partner's at the centre). `bcut_complex/I6`'s notched tool is refused on
import as a groove (S9e.4b.3c's), no longer by its own arrangement's
tangency; `G4` is refused as before. Of the 171, 41 evaluate on both
backends, 68 are refused by S9's rules, 15 are S9e.4b.3c's, 36 are bodies
none of the kernel's constructions, 6 are S9e.4b.4's, 3 have spline faces
(S9f) and 2 give the next Boolean several solids; native DRAW does not
evaluate 24 of them. Of the 8 cases that timed out at their first
restores, 6 still do; `bugs/modalg_1/buc60532_2` and
`bugs/modalg_6/bug23585` ended in 113 and 102 seconds, unsupported (a
`SurfaceOfLinearExtrusion` the reader does not represent; `tolerance p`),
by the host's load: alone, the worker of the last survey's commit takes 77
to 119 seconds on them, this one 78 to 116.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,014 registered Boolean cases and the 14 new ones): native
DRAW's values are the last audit's bit for bit; the Rust adapter's
statuses are the same and 2 cases' values moved within rounding,
`bfuse_complex/F5` and `Q2` (one pair of CTO prisms fused): 4 centre
coordinates, at most 4.7e-16 relative (4.6e-16 of the solid's size), no
volume or area, at the prism caps' frames (`59d0c57b`; the worker of its
parent gives the last audit's bits). Both still agree with native DRAW's
printed digits, and the same 35 disagreements remain, native off in each.
The 14 new cases' volume, area and centre agree with native DRAW's printed
digits. They are registered (`data`, `viewer_skipped` on both; 1,114
cases); the ledger does not change. A full contract run of the manifest
(both backends, the dataset present, 30 seconds a case) holds for every
case, the 14 in 0.6 to 3.3 seconds on the Rust adapter, the slowest
Boolean case 5.3 seconds (`boptuc_simple/ZK8`, on a host at load 7 to 9;
`bopfuse_simple/ZK8` 6.0 in the last survey), the restore cases 0.2 to
4.1, the rollex 0.6 to 0.7.

**S9e.4b.3c.2 and S9e.4b.3c.3a: the Boolean group and the restore
cases.** Both sets were run again after S9e.4b.3c.2 (exact incidences of
two inputs on one sphere: a vertex of both one vertex, a vertex of one
inside the other's edge splitting it, an edge of both one edge of both
inputs' faces) and S9e.4b.3c.3a (an imported body of one curved face and
planes as one Boolean of its primitive and the hull of its other planes: a
groove, a bite or a boss) (2026-10-05, at `b9c7ae1b`, the public dataset,
120 seconds a case, four at once). The 1,802 self-contained cases of the
Boolean group on both backends: every status, reason and error is the
last survey's (`f4b584f7`) field for field, 987 evaluating, every one
registered, 592 refused, 223 unsupported on both; none fails or times
out. The 1,814 restore cases on the Rust adapter, and the 171 the import
reaches on native DRAW too: native DRAW's statuses and reasons are the
last survey's, and on the Rust adapter only 15 cases move, as the two
steps' trials found. 13 that S9e.4b.3c.2 targets evaluate on both
backends with every check: the hemisphere and a wedge (`so1` and `so2`,
the wedge's equator arc on the hemisphere's rim, its pole the
hemisphere's) in `bcommon_complex/B3`, `bcut_complex/C7` and `C8` (the
wedge less the hemisphere, empty), `bfuse_complex/B6` and
`bugs/moddata_1/bug183_4`; two wedges with a corner, a pole and an axis
edge of both and their bases on one plane (`so2` and `so3`) in
`bcommon_complex/B4`, `bcut_complex/C9` and `D1` and `bfuse_complex/B7`;
and two wedges, the higher one's axis edge inside the other's (`so5` and
`so2`), in `bcommon_complex/B6`, `bcut_complex/D4` and `D5` and
`bfuse_complex/B9`. `bcut_complex/G4` (a box with a cylindrical boss less
a parallel cylinder, the boss its hull fused with its primitive) evaluates
on both backends with every check. `bcut_complex/I6` is refused for its
notched tool's own tangency (`Degenerate`: an imported plane piece whose
curved face is tangent to its plane faces), no longer as S9e.4b.3c's;
`bugs/modalg_1/buc60926` stays `Degenerate` (a plane through a cone's
apex) and the 11 `so` cases with `so6` or `so7` by S9's rules. Of the
171, 55 evaluate on both backends, 69 are refused by S9's rules, 36 are
bodies none of the kernel's constructions, 6 are S9e.4b.4's, 3 have
spline faces (S9f) and 2 give the next Boolean several solids; none is
S9e.4b.3c's any more, and native DRAW does not evaluate 24 of them. The
same 6 cases time out at their first restores; `bugs/modalg_1/buc60532_2`
and `bugs/modalg_6/bug23585` end within the 120 s again (104 and 100
seconds), unsupported as before.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,028 audited before and the 14 new cases): both backends'
values are the last audit's bit for bit on all 1,028, the statuses the
same, and the same 35 disagreements remain, native off in each. The 14 new
cases' volume, area and centre agree with native DRAW's printed digits
(`C8`'s cut empty on both; `G4`'s volume 4,812,276.20 and area
193,592.99, native's printed 4.81228e6 and 193593). They are registered
(`data`, `viewer_skipped` on both; 1,128 cases); the ledger does not
change. A full contract run of the manifest (both backends, the dataset
present, 30 seconds a case) holds for every case, the 14 in 0.6 to 1.2
seconds on the Rust adapter (`G4` the fastest, `bfuse_complex/B9` the
slowest), the slowest Boolean case 4.4 seconds (`bopfuse_simple/ZK8`, on a
host at load 5 to 8; `boptuc_simple/ZK8` 5.3 in the last survey), the
restore cases 0.1 to 3.6 (`bcommon_complex/B9`), the rollex 0.5.

**S9e.4b.3c.3b, S9e.4b.4a and the scheduled replay's fixes: the Boolean
group and the restore cases.** Both sets were run again after the fixes
of CI's scheduled full replay (a section within rounding of a sphere's
pole a computation limit; a sheet or wire split by a plane within the
resolution of its own plane `Degenerate`), S9e.4b.3c.3b (an imported body of
one curved face and planes as a Boolean tree of its primitive and several
hulls) and S9e.4b.4a (an imported prism whose arcs of two circles meet at
a joint, each arc taken through its ends; two parallel cylinders within
the resolution of one and not one `Degenerate`) (2026-10-05, at
`16121052`, the public dataset, 120 seconds a case, four at once). The
1,802 self-contained cases of the Boolean group on both backends: every
status, reason and error is the last survey's (`b9c7ae1b`) field for
field, 987 evaluating, every one registered, 592 refused, 223 unsupported
on both; none fails or times out. The 1,814 restore cases on the Rust
adapter, and the 171 the import reaches on native DRAW too: native DRAW's
statuses and reasons are the last survey's, and on the Rust adapter only
S9e.4b.4a's 6 cases of two circles at a joint move, as its trial found.
`bcut_complex/P4` (a box less a prism inside it whose fillets of radii 2
and 0.47 meet tangent from either side) evaluates on both backends with
every check; `bcut_complex/E8`, `bfuse_complex/D5` and `E1` are refused
as two cylinders within the resolution of one cylinder (each tool on a
stored circle of the part's arcs, which the part now takes through their
ends), `bugs/modalg_2/bug4993_1` and `_2` as two faces within the
resolution of one plane. The parallel cylinders' rule moves no registered
case and no self-contained case's status or reason, and refuses only
`E8`, `D5` and `E1`, none evaluating before; with the rule disabled in a
scratch build `E8` evaluates with native DRAW's volume and area to its
printed digits (the rule costs it a result: the open keeping of a circle
the tool holds exactly), `D5` is a computation limit (a face's loop
winding undecided) and `E1` a tangency, so the rule renames their
refusals only. Of the 171, 56 evaluate on both backends, 74 are refused
by S9's rules, 36 are bodies none of the kernel's constructions, 3 have
spline faces (S9f) and 2 give the next Boolean several solids; none is
S9e.4b.4's any more, and native DRAW does not evaluate 24 of them. Six
cases time out at their first restores; the host's load (10 to 64) moved
two of the reader's slowest across the 120 s (`bugs/modalg_1/buc60532_2`
and `bugs/modalg_6/bug23585` timing out, `buc60532` and `buc60532_1`
ending in 114 seconds), and alone the four end in 88 to 117 seconds,
unsupported as before.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,042 audited before and `P4`): both backends' values are
the last audit's bit for bit on all 1,042, the statuses the same, and the
same 35 disagreements remain, native off in each. `P4`'s volume
88,171,481.37, area 1,642,145.80 and centre agree with native DRAW's
printed digits (8.81715e7, 1.64215e6). It is registered (`data`,
`viewer_skipped` on both; 1,129 cases); the ledger does not change. A
full contract run of the manifest (both backends, the dataset present, 30
seconds a case) holds for every case, `P4` in 3.0 seconds on the Rust
adapter, the slowest Boolean case 5.8 seconds (`bopcut_simple/ZK8`, on a
host at load 10 to 45; `bopfuse_simple/ZK8` 4.4 in the last survey), the
restore cases 0.2 to 5.1 (`bcommon_complex/B9`), the rollex 0.6 to 0.7.

**S9e.4b.4b.1: the Boolean group and the restore cases.** Both sets were
run again after S9e.4b.4b.1 (an imported body of several sphere, cylinder
and cone faces whose plane faces are all their primitives' ends, as a
Boolean chain of those primitives) (2026-10-06, at `c62470a3`, the public
dataset, 120 seconds a case, four at once). The 1,802 self-contained
cases of the Boolean group on both backends: every status, reason and
error is the last survey's (`16121052`) field for field, 987 evaluating,
every one registered, 592 refused, 223 unsupported on both; none fails or
times out. The 1,814 restore cases on the Rust adapter, and the 171 the
import reaches on native DRAW too: native DRAW's statuses and reasons are
the last survey's, and on the Rust adapter only the 36 bodies the import
refused as none of its constructions move, as the step's trial found.
`bfuse_complex/E5` (a stepped shaft of two coaxial cylinders fused with a
box on its step's plane), `bcut_complex/G9` and `bugs/modalg_2/bug417` (a
dome, a frustum common a ball, less a pin, a rod common the same ball)
evaluate on both backends with every check; `bugs/modalg_6/bug28773` is
refused as two cylinders' axes within rounding of parallel (its tube a
prism on a disc frame leaning 2.2e-33 off its walls), `bfuse_complex/K1`
as S9e.4b.4b.2's (plane faces other than its primitives' ends); the other
31 keep the general refusal under its new text ("an imported solid other
than a prism, a sphere, a cone, a torus or a Boolean of its primitives"):
`bugs/modalg_2/bug476_1` to `_8` and 23 cases of the `_2d` grids, a solid
of several curved faces (`case_8_solid_repaired`) against wires. Of the
171, 59 evaluate on both backends, 75 are refused by S9's rules, 31 are
bodies none of the kernel's constructions, 1 is S9e.4b.4b.2's, 3 have
spline faces (S9f) and 2 give the next Boolean several solids; native
DRAW does not evaluate 24 of them. Four cases time out at their first
restores, as in the last three surveys (`bugs/modalg_1/buc60531_1`, `_2`,
`bugs/modalg_5/bug23849_1`, `_3`); at this survey's load (3 to 15) the
reader's other slowest end within the 120 seconds (`buc60532_2` and
`bugs/modalg_6/bug23585` in 97 and 94, `buc60532` and `buc60532_1` in
114), unsupported as before.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,043 audited before and the three): both backends' values
are the last audit's bit for bit on all 1,043, the statuses the same, and
the same 35 disagreements remain, native off in each. E5's volume
97,940,759.17, area 1,831,656.94 and centre agree with native DRAW's
printed digits (9.79408e7, 1.83166e6), G9's and `bug417`'s volume
8,941.5371 and area 2,814.0138 with native's (8941.54, 2814.01). The
three are registered (`data`, `viewer_skipped` on both; 1,132 cases); the
ledger does not change. A full contract run of the manifest (both
backends, the dataset present, 30 seconds a case) holds for every case,
E5, G9 and `bug417` in 0.7, 1.4 and 1.3 seconds on the Rust adapter, the
slowest Boolean case 3.9 seconds (`bopfuse_simple/ZK8`, on a host at load
4 to 15; `bopcut_simple/ZK8` 5.8 in the last survey), the restore cases
0.1 to 3.3 (`bcommon_complex/B9`), the rollex 0.4.

**S9e.4b.4b.2a, S9e.4b.4b.2b.1 and the cylinder pairs' fixes: the Boolean
group and the restore cases.** Both sets were run again after S9e.4b.4b.2a
(such bodies led by a prism leaf), S9e.4b.4b.2b.1 (their other plane faces
as a primitive's flat, a prism of one cap or a pocket) and two general
rules of the curved engine, a result touching itself at a vertex refused
and a wider near-node margin for crossing cylinders (2026-10-09, at
`826346b7`, the public dataset, 120 seconds a case, four at once). Neither
rule moves any case. The 1,802 self-contained cases of the Boolean group
on both backends: every status, reason and error is the last survey's
(`c62470a3`) field for field, 987 evaluating, every one registered, 592
refused, 223 unsupported on both; none fails or times out. The 1,814
restore cases on the Rust adapter, and the 171 the import reaches on native
DRAW too: native DRAW's statuses and reasons are the last survey's, and on
the Rust adapter one case moves, `bfuse_complex/K1`, now refused as a
tangency between the inputs (its part imported as its rounded box less its
bore, its tool a rod of the bore's radius whose axis crosses it) where it
was refused as S9e.4b.4b.2's, as S9e.4b.4b.2a's trial found. Replayed at
`1cd7bb3a` (S9e.4b.4b.2a, before either rule), K1 and the 28 cases refused
as a near node or as a result or face touching itself keep the same
refusals. Of the 171, 59 evaluate on both backends, 76 are refused by S9's
rules, 31 are bodies none of the kernel's constructions, 3 have spline
faces (S9f) and 2 give the next Boolean several solids; native DRAW does
not evaluate 24 of them. Eight cases time out at their first restores: the
four of the last four surveys (`bugs/modalg_1/buc60531_1`, `_2`,
`bugs/modalg_5/bug23849_1`, `_3`) and, at this survey's load (8 to 18),
`buc60532`, `buc60532_1`, `buc60532_2` and `bugs/modalg_6/bug23585`, which
end alone in 68 to 73 seconds, unsupported as before.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,046 audited before): both backends' values are the last
audit's bit for bit on all 1,046, the statuses the same, and the same 35
disagreements remain, native off in each. No case newly evaluates, so none
is registered (1,132 cases); the ledger does not change. A full contract
run of the manifest (both backends, the dataset present, 30 seconds a case)
holds for every case, the slowest Boolean case 9.8 seconds
(`bopfuse_simple/ZK8`, on a host at load 7 to 21; 3.9 in the last
survey), the restore cases 0.2 to 7.0 (`bcommon_complex/B9`), the rollex
0.7 to 1.8.

**S9e.4b.4c.1, the near miss of a sphere and the user's crossing
cylinders: the Boolean group and the restore cases.** Both sets were run
again after S9e.4b.4c.1 (imported polyhedra against curved faces, their
cavities, and the assemblies' cavities each a shell of its own), the
user's decision on equal cylinders whose axes cross (a bore against a rod
of its radius a tangency between the inputs in every operation) and a
sphere missing a plane face or another sphere within the resolution
refused as S9d.1's near tangency (2026-10-10, at `44204e18`, the public
dataset, 120 seconds a case, four at once). None of the three changes
moves any case. The 1,802 self-contained cases of the Boolean group on
both backends: every status, reason and error is the last survey's
(`826346b7`) field for field, 987 evaluating, every one registered, 592
refused, 223 unsupported on both; none fails or times out, and none is
refused as a near miss. The 1,814 restore cases on the Rust adapter, and
the 171 the import reaches on native DRAW too: every status and reason is
the last survey's on both backends, but `bugs/modalg_1/buc60532`,
`buc60532_1`, `buc60532_2` and `bugs/modalg_6/bug23585`, which timed out
at the last survey's load and now end in 98 to 104 seconds, unsupported as
in the surveys before it; `bugs/modalg_1/buc60531_1`, `_2`,
`bugs/modalg_5/bug23849_1` and `_3` time out at their first restores as
before. Of the 171, 59 evaluate on both backends,
76 are refused by S9's rules (`bfuse_complex/K1` a tangency between the
inputs, as before), 31 are bodies none of the kernel's constructions, 3
have spline faces (S9f) and 2 give the next Boolean several solids; native
DRAW does not evaluate 24 of them. An instrumented build (not committed)
over the 1,802 and the 171 finds no case reaching the crossing cylinders'
new tangency, a sphere's near miss, an imported polyhedron's triangles or
a cavity, so no case can move with them.

The volume audit (`vprops` and `sprops` before each `checkprops`, both
backends, the 1,046 audited before): both backends' values are the last
audit's bit for bit on all 1,046, the statuses the same, and the same 35
disagreements remain. No case newly evaluates, so none is registered
(1,132 cases); the ledger does not change. A full contract run of the
manifest (both backends, the dataset present, 30 seconds a case) holds for
every case, the slowest Boolean case 5.4 seconds (`bopcut_simple/ZK8`, on
a host at load 6 to 10), the restore cases 0.1 to 3.1
(`bcommon_complex/B9`), the rollex 0.4.

Three more `intss` cases run on the Rust adapter but are not registered,
because the contract admits no failing status: `bug23178`, `bug28222_2` and
`bug28222_3` count the pieces IntPatch splits its walking lines into (6, 4
and 2), where the kernel returns closed loops (1, 2 and 1) whose `xdistcs`
samples lie within 1.5e-14 of both cylinders.
The thirteen derived cases are counted separately (see below).
The 29 bridge self-tests are separate infrastructure checks; they do not count
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

`boolean_prisms` (S9a) runs OCCT's Boolean commands on prisms in one frame;
the adapter maps each to `Solid::fuse`, `cut` or `common` (`btuc`,
`boptuc` and `bbop`/`bapibop` 3 as the tool cut by the object) and returns
OCCT's result, a compound of the result's solids (a compound even of one
solid, and an empty compound for an empty result, as native DRAW prints
them). The result's solids share nothing. OCCT's raw result keeps each
argument's face and edge images (a fused wall in two faces, a stacked
fuse's walls split at the joint), while the kernel merges collinear edges
and coplanar faces as `unifysamedom` does: volumes and areas are checked on
the raw result, counts and lengths per use on `unifysamedom`'s, which the
adapter returns unchanged when each solid is a prism with no two
consecutive profile segments collinear lines or arcs of one circle, or a
stack (S9a.2, below; any other `unifysamedom` is unsupported). Its values
are computed by hand and hold on both backends:

* 4 x 4 x 2 boxes overlapping in a 2 x 2 square: the fuse an octagon prism
  (volume 56, area 104, length 128; 16 vertices, 24 edges, 10 wires and
  faces), the cut a hexagon prism (24, 56, 88; 12, 18, 8, 8), the common a
  2 x 2 x 2 box (8, 24, 48), and again through `bop` with `bopfuse`,
  `bopcut`, `boptuc` (the mirror hexagon) and `bopcommon`, and through
  `bbop 0` and `bapibop 2` on the General Fuse arguments;
* a radius 1 `pcylinder` and a box on its `-x` side: a half cylinder
  (volume `pi`, area `3 pi + 4`, length `4 pi + 16`; 4 vertices, 6 edges, 4
  wires and faces), its arc away from OCCT's seam at `+x`;
* the cylinder cut through a 4 x 4 x 2 box: volume `32 - 2 pi`, area
  `64 + 2 pi`, length `84 + 8 pi` (the seam twice); 10 vertices, 15 edges,
  9 wires, 7 faces, through the count synthesizer;
* a bar spanning a 6 x 2 x 2 box's heights cuts it into two solids (volume
  20, area 56, length 104; 16 vertices, 24 edges, 2 shells and solids, one
  compound);
* disjoint unit boxes: the fuse keeps both, the common is an empty compound
  (`checkprops -s empty -v empty`, one shape);
* a box inside a 10 x 10 x 10 box: the fuse is the outer box, the common
  the inner one; a 2 x 2 bar through its heights leaves a square hole
  (volume 960, area 672; 12 wires);
* stacked 2 x 3 boxes touching at `z = 1`, then overlapping over `[2, 3]`:
  one box each time (volume 18, then 24), a result being an argument again;
* a stadium prism (`profile` with half circles, `Copy`) and a bar across
  it: a cross of 12 profile vertices (volume `20 + 2 pi`, area `48 + 6
  pi`, length `104 + 8 pi`).

`boolean_stacks` (S9a.2) runs the same commands on prisms whose result's
height slabs hold different regions. The kernel builds these results as
stacks (`Construction::Stack`): general bodies, no longer prisms, with
certified mass properties, a closed cavity a second shell of its solid.
`checkshape`, `vprops`, `sprops`, `lprops` and `nbshapes` read the stack's
topology as any Boolean result's (the count synthesizer's seams and seam
vertices for its cylinder walls, `Topology::occt_counts` for its shells).
The kernel's stacks are unified already (walls on one line or circle
facing one way joined across slab heights and piece ends, no vertex where
an edge runs straight on between the same two faces), so `unifysamedom`
returns them unchanged; native DRAW's unified counts agree on every stack
checked, those below and a tower, an L step, a stepped shaft, a plug
filling a hole over part of its height and a cylindrical boss. Its values
are computed by hand and hold on both backends:

* a 2 x 4 x 2 box on half of a 4 x 4 x 2 box (`bfuse`, and `bop` with
  `bopfuse`): an L-shaped prism along y, volume 48, area 88, length 112; 12
  vertices, 18 edges, 8 wires and faces;
* a radius 1 cylinder over `[1, 2]` cut from a 4 x 4 x 2 box (`bcut`, and
  `bopcut`), a pocket: volume `32 - pi`, area `64 + 2 pi`, length `82 + 8
  pi` (the seam of length 1 twice); 10 vertices, 15 edges, 9 wires, 8
  faces;
* a slab over `[2, 3]` across a 3 x 2 x 6 box cuts it into two solids:
  volume 30, area 74, length 120; 16 vertices, 24 edges, 12 wires and
  faces, 2 shells and solids;
* a 2 x 2 x 2 box strictly inside a 4 x 4 x 4 box (`bcut`, and `bbop 2`),
  a closed cavity: volume 56, area 120, length 144; 16 vertices, 24 edges,
  12 wires and faces, 2 shells, one solid;
* a 3 x 1 x 1 box over `[1, 2]` entering a radius 2, height 4 cylinder
  from `-x` to `x = -1` (away from OCCT's seam), a tool through a round
  wall: with `s = sqrt(15) / 2` and `t = asin(1/4)`, volume `16 pi - (s /
  2 + 4 t - 1)`, area `24 pi + 4 t + 3 s - 3`, length `16 pi + 12 + 16 t +
  8 s`; 10 vertices, 15 edges, 9 wires, 8 faces.

`boolean_polyhedra` (S9b.1) runs them on prisms in frames with different
axes, the tools copied by `tcopy` and turned by `trotate`. The kernel
builds these results as polyhedra (`Construction::Polyhedron`): general
bodies decided on the inputs' exact models, unified as built (coplanar
fragments facing one way joined into maximal faces, edges joined where they
run straight on between the same two faces), so `unifysamedom` returns them
unchanged as it does stacks; native DRAW's unified counts agree on each
below. A `trotate` by whole quarter turns about a coordinate axis turns a
prism's frame exactly: the kernel's rotation takes the cosine of the
rounded quarter turn, 6.1e-17, as OCCT's `gp_Trsf` does, and where OCCT's
tolerances make a turned wall coplanar with another, the kernel's exact
decisions would keep a crease (the L below would have 14 vertices, 21
edges and 9 faces unified, where OCCT has 12, 18 and 8). With `r =
sqrt(2)`, its values are computed by hand and hold on both backends (to
`-deps 1e-5` where native DRAW prints six digits of an irrational value):

* a 2 x 1 x 1 box and its copy turned 90 degrees about z, meeting along a
  face: an L-shaped prism (volume 4, area 18, length 52; 12 vertices, 18
  edges, 8 wires and faces), its walls at `y = 0` one face across both
  boxes;
* a 4 x 2 x 1 bar and its copy turned 90 degrees about z: a plus sign
  (volume 12, area 40, length 88; 24 vertices, 36 edges, 14 wires and
  faces), the 2 x 2 middle (4, 16, 40; 8, 12, 6), the two ends (4, 20, 64;
  16, 24, 12, 2 shells and solids), and the other bar's ends (`boptuc`);
* a 2 x 2 x 3 bar over `[-1, 2]` turned 45 degrees about the centre of the
  2 x 2 x 1 box it passes through: the common a regular octagon of
  inradius 1 (volume `8 r - 8`, area `32 r - 32`, length `64 r - 48`; 16
  vertices, 24 edges, 10 wires and faces), the cut four corner triangles
  of legs `2 - r` (volume `12 - 8 r`, area `32 - 16 r`, length 56; 24
  vertices, 36 edges, 20 wires and faces, 4 shells and solids), the fuse
  (volume `24 - 8 r`, area `80 - 32 r`);
* a 1 x 4 x 8 bar turned 45 degrees about the y axis through `(2, 0, 1)`,
  across a 4 x 2 x 2 box: the cut two trapezoid prisms (volume `16 - 4 r`,
  area 40, length `80 + 8 r`; 16 vertices, 24 edges, 12 wires and faces, 2
  shells and solids), the common a parallelogram prism (volume `4 r`, area
  `16 r`, length `24 r + 16`; 8, 12, 6);
* a unit box turned 30 degrees about `(1, 1, 1)` through its centre,
  strictly inside a 4 x 4 x 4 box (`bcut`, `bbop 2`): a closed cavity,
  volume 63, area 102, length 120; 16 vertices, 24 edges, 12 wires and
  faces, 2 shells, one solid; the common is the unit box.

No spline profile is among them: the adapter builds no spline edge
(`bsplinecurve`, and `wire` or `edge` of edges, are unsupported), so
S9a.2's spline profiles are not reachable through DRAW yet.

The adapter reports unsupported, never an answer, whatever the kernel
refuses: two cylinders meeting in quartics (S9c.2), a sphere meeting a
cylinder in a turned frame in a loop (S9d.2b's), a cone against a prism
with arcs, a sphere or a cone (S9d.3b's), tori (still refused as S9c's
curved solids; S9d.4), a plane through a cone's apex, an arc ending off
its circle and a result with arcs as an argument (S9c), results touching
themselves or each other, faces and pieces thinner than the resolution, a
face using an edge both ways or touching itself at a vertex, two solids
sharing a vertex, tangencies and meetings within rounding of arcs in any
position (`Degenerate`), `ComputationLimit`, `LimitExceeded` and `PrecisionLoss`;
also several objects or tools, a section (`bbop`/`bapibop` 4,
`bopsection`), an argument of several solids, two stacks or polyhedra
sharing ids, a cone, a sphere or a torus sharing ids with the other
argument, restored shapes and split results. Two things the kernel's
Boolean requires of its inputs the adapter supplies: their ids apart
(every `box` is a cuboid of the unspecified operation, and a `copy` keeps
its ids), and a construction indexed by profile element (a Boolean
result's entities descend from its inputs); a tool sharing ids with the
object, and a Boolean result that is a prism used as an argument, are
extruded again from their profile under an operation of their own. A stack
or a polyhedron is taken as it is (S9b.2, on its stored geometry), the
other argument extruded again when they share ids. `savehistory` after a
Boolean is unsupported.

The history cases use `prism ... Copy`. Without `Copy`, OCCT builds the prism's end face
as the start face moved by a location, reusing its `TShape`s. DRAW's
`nbshapes` counts those shared shapes once (4 vertices, 8 edges and 5 faces
for a box), while the kernel, like `Copy`, builds distinct end entities.
Native DRAW with `Copy` reports 8, 12 and 6. Since S9a the adapter builds
`prism` without `Copy` (the Boolean group uses it) as an uncopied prism:
volumes, areas and lengths per use agree (OCCT's explorer visits the moved
shapes under their locations), `nbshapes` of it or of a Boolean of it (whose
unsplit faces OCCT reuses, 12 vertices where a copied prism's cut has 16)
and its history are unsupported.

The two `bug29333` cases are the self-contained history candidates for M3's
split and fuse (`IDENTITY_AND_HISTORY.md`). They split and fuse faces made by
`plane`, `mkface`, `line` and `mkedge`, which the Rust adapter cannot build, so
they are capability sentinels: native DRAW must pass them, and Rust reports
them unsupported (`bug29333_1` moves a face with `ttranslate`, `bug29333_2`
needs `mkvolume`). `bug21264` also uses splits and Booleans; since S8e it is
registered with the `bsplit` group. The host forwards their commands to
either backend. `split_plane` and `boolean_prisms` check counts and
properties only: the adapter keeps no split or Boolean history for
`savehistory`, since the kernel's history names each piece's own cut face,
not OCCT's shared one, and merges what OCCT keeps apart.
`compare_split_merge.py` and `compare_boolean.py` compare the relations
instead.

## Deliberate limits

- Rust signatures: positional `box` with three or six numbers, `pcylinder name
  radius height` on the default axis, two-name `copy` and `tcopy`, single-shape
  `ttranslate`/`trotate` (a prism turned by whole quarter turns about a
  coordinate axis in an exactly turned frame), `checkshape`, unique `nbshapes` with synthesized
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
  For Booleans (S9a): `bfuse`, `bcut`, `bcommon`, `btuc name object tool`;
  `bop object tool` then `bopfuse`, `bopcut`, `boptuc` or `bopcommon name`;
  `bbop name 0..3` after `bfillds` and `bapibop name 0..3` with one object
  and one tool, each a solid or a Boolean result of one solid, a stack or a
  polyhedron too (S9b.2); `unifysamedom name result` of a unified Boolean result, a
  stack's or a polyhedron's too (S9a.2, S9b.1; see `boolean_prisms`,
  `boolean_stacks` and `boolean_polyhedra` above); `checkshape`,
  `nbshapes`, `vprops`, `sprops` and `lprops` of any Boolean result, a
  stack's or a polyhedron's too. A Boolean result is not moved
  (`ttranslate` or `trotate` of it is unsupported). `prism` without `Copy` builds an uncopied prism
  (above); `ttranslate` also moves a `profile` sketch; `profile ... C r 360`
  is a whole circle. Numbers of `box`, `pcylinder`, `pcone`, `psphere`,
  `ptorus`, `ttranslate`, `trotate`, `polyline`, `prism` and `profile` may
  be expressions of `dset` variables with `+`, `-`, `*`, `/`, `pi` and
  `sqrt` (Tcl expressions are evaluated by Tcl first); other functions are
  unsupported.
  `psphere name [plane] R [angle1 angle2]` places a sphere, cap or zone
  on the z axis at the origin or, since S9d.1's survey, on a DRAW `plane`'s
  frame (its origin the centre, its normal the axis, as OCCT's `gp_Ax2` of
  the plane); a partial longitude is unsupported.
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
restore data and run checks on both backends, without registering them
(`--boolean` runs the Boolean group's self-contained cases instead; a
dataset fetched in another checkout is read through `--data-dir`, and
`--jobs` runs cases side by side). With
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
